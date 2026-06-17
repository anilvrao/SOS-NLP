      subroutine   xdsle5  ( pvttol, zpcntl, qroot , pnlrow, pnlcol, 
     1                       pospon, offset, invp  , panel , pvtblk, 
     2                       corder, lmaxvl, swap  , ncolf , inrtia, 
     3                       fctops, slvops, ppfmon, xpboxs, psboxs, 
     4                       inpexp, inpsiz, inzsiz, izfail, rtpexp, 
     5                       rtpsiz, rtzsiz, rzfail, ierr )

 
c
c  purpose -- perform a block elimination step, eliminating as many
c             columns from a panel of a front as possible.  We use
c             threshold duff-reid pivoting for stability.  The pivoting
c             is exhaustive only within the panel.  Any columns that are
c             not eliminated do not have stable pivots within the panel,
c             but may have satisfactory pivots within the remainder of
c             the front.
c
c             the low rank outer product modification step is performed
c             only to columns in this panel.
c
c             this is the column oriented panel version
c
c  created            -- 21-feb-97, rgg
c  last modified      -- 27-aug-98, rgg, added zero pivot controls
c                     -- 04-oct-01, jgl, restructured to repair test for
c                                        singularity; error codes
c                                        changed
c                     -- 22-mar-02, jgl, repair previous restructuring
c                                        to properly handle perturbed
c                                        diagonal case
c
c  input variables --
c
c      pvttol -- pivoting tolerance
c      zpcntl -- zero pivot control structures
c      qroot  -- logical flag which if true indicates that this is
c                the root front
c      pnlrow -- number of rows in the panel 
c      pnlcol -- number of columns in the panel
c      pospon -- number of postponed columns in this front
c      offset -- swapping offset
c      invp   -- inverse permutation array.
c      panel  -- panel of the current frontal matrix to be eliminated
c      pvtblk -- integer vector indicates if diagonal element is a
c                1 x 1 or 2 x 2 block pivot.
c      inrtia -- matrix inertia
c      fctops -- factor operation count
c      slvops -- solve operation count
c      ppfmon -- monitor pivot failures or not
c      xpboxs -- number of pivot exponent boxes for pivot monitoring
c      psboxs -- number of panel size boxes for pivot monitoring
c
c  working storage --
c
c      corder -- temporary array of length  pnlcol
c      lmaxvl -- temporary array of length  pnlcol
c
c  output variables --
c
c      invp   -- inverse permutation array.
c      swap   -- swapping record
c      panel  -- panel of the current frontal matrix to be eliminated
c      ncolf  -- the number of columns that were eliminated.  
c      inrtia -- matrix inertia
c      fctops -- factor operation count
c      slvops -- solve operation count
c      inpexp -- base 10 exponents for minimum pivots, interior panels
c      inpsiz -- base 2 exponents for interior panel sizes, pivot fails
c                for nonzero panel block
c      inzsiz -- base 2 exponents for interior panel sizes, pivot fails
c                for identically zero panel block
c      izfail -- count of failures due to zero interior panel block
c      rtpexp -- base 10 exponents for minimum pivots, root panels
c      rtpsiz -- base 2 exponents for root panel sizes, pivot fails
c                for nonzero panel block
c      rtzsiz -- base 2 exponents for root panel sizes, pivot fails
c                for identically zero panel block
c      rzfail -- count of failures due to zero root panel block
c      ierr   -- error return,
c                if ierr = 0, success,
c                        = -1, detected exactly singular matrix
c                              (a column of the reduced matrix is
c                               exactly zero)
c                        = -2, pivot algorithm failure -- no pivots
c                              found in a panel of the root front
c
c     ------------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer            pnlrow, pnlcol, pospon, offset, 
     1                   ncolf , ierr
 
      integer            corder (*),     invp (*),
     1                   pvtblk (*),     inrtia (3),
     2                   swap   (*)

      logical            qroot
 
      double precision   pvttol, fctops, slvops
 
      double precision   panel (pnlrow,pnlcol),  lmaxvl (*),
     1                   zpcntl(*)
 
c     ------------------------------------
c     ... global variables for diagnostics
c     ------------------------------------

      integer            xpboxs, psboxs,  izfail, rzfail,
     1                   inpexp (xpboxs), inpsiz (psboxs),
     2                   inzsiz (psboxs), rtpexp (xpboxs),
     3                   rtpsiz (psboxs), rtzsiz (psboxs)

      logical            ppfmon

c     -------------------
c     ... local variables
c     -------------------
 
      integer            expon , j     , jj    , jpiv  , jpost , k     ,
     1                   kk    , kpiv  , lc    , lg2pcol, ll   , m     ,
     2                   pcol  , psize , triggr

      logical            pmonit

      double precision   diag1 , f     , gamma , lambda, 
     1                   t1    , t2    , t3    , diag2 , ttmax, zero

      double precision   mintol

      parameter        ( zero = 0.0 )

c     --------------------
c     ... subprograms used
c     --------------------
 
      external  xdsle6, xdsle7, xdsle8

c     ------------------------------------------------------------------

      ierr   = 0
      m      = 1
      jpost  = 0
      pmonit = ppfmon
c.debug
c     write(6,'("at start of xdsle5")')
c     call xdslp5('panel', pnlrow*pnlcol, panel, 6 )
c.debug

c     ----------------------------------
c     ... search for an acceptable pivot
c     ----------------------------------
 
   10 continue
      kpiv   = 0
c$$$      if  ( qroot )  then
c$$$         write (*,*) 'NEW SEARCH'
c$$$      endif

c     ----------------------------------------------------------------
c     ... on first iteration, we search postponed columns last.
c         thereafter, we start the search with the column immediately
c         after the last column chosen in the previous step.
c         (could use linked list to avoid next two loops)
c     ----------------------------------------------------------------

      do j = m + jpost, pnlcol
          corder (j-jpost) = j
      enddo

      do j = m, m+jpost-1
          corder (pnlcol-(m+jpost-1)+j) = j
      enddo

      jpost = 0

c     -----------------------------------------------------
c     ... loop through remaining columns searching
c         for a 1x1 or 2x2 pivot
c     -----------------------------------------------------

      triggr = min ( m + 4, pnlcol )

      ttmax  = 1.
      mintol = 0.
 
      do 40 jj = m, pnlcol
c.debug
c     write(6,'("in 40 - jj, m, pnlcol = ", 3i8)')
c    1                    jj, m, pnlcol 
c.debug
  
          j     = corder(jj)
          diag2 = abs ( panel(j,j) )
 
c         --------------------------------------------------------
c         ... find largest off-diagonal element in column j in the
c             supernode+parent (gamma )
c         --------------------------------------------------------
 
          gamma  = -1.

          do k=m,j-1
            gamma = max(gamma,abs(panel(j,k)))
          enddo

          do k=j+1,pnlrow
            gamma = max(gamma,abs(panel(k,j)))
          enddo

          lmaxvl (j) = gamma 

          ttmax = max ( ttmax, gamma )

c$$$          if  ( qroot )  then
c$$$             write (*,*) 'pivot test for column:', j
c$$$             write (*,'(5x, 1p2e12.3)') diag2, gamma
c$$$          end if

          if  ( gamma .gt. zero )  then

             if  ( pmonit )  then
                mintol = max ( mintol, diag2 / gamma )
             end if

c            ---------------------------------------------------------
c            ... general case: first check to see if current column
c                can be used as part of a 2 x 2  pivot with any column
c                preceding it in the interior of the front
c            ---------------------------------------------------------

             do kk = m, jj-1
c.debug
c     write(6,'("in 25 - kk, m, jj     = ", 3i8)')
c    1                    kk, m, jj     
c.debug

                k      = corder(kk)
                diag1  = abs ( panel(k,k) )
                lambda = lmaxvl ( k )

                if ( j .gt. k ) then
                   f = abs ( panel(j,k) )
                else
                   f = abs ( panel(k,j) )
                end if

c               ------------------------------------------------
c               ... use explicit Duff-Reid bound to determine if
c                   2x2 pivot is satisfactory
c               ------------------------------------------------

                t1    = diag2 * lambda + f * gamma
                t2    = diag1 * gamma  + f * lambda
                t3    = abs ( diag1 * diag2 - f ** 2 )

                if ( t3 .gt. zero .and. max (t1,t2) .le. t3/pvttol )
     1          then
 
                   jpiv   = min ( k, j )
                   kpiv   = max ( k, j )

c                  if ( jpiv .eq. m+1 ) then
c                      jpiv = kpiv
c                      kpiv = m+1
c                  end if

                   go to 200

                else

                   if  ( pmonit )  then
                      if  ( max (t1,t2) .gt. zero )  then
                         mintol = max ( mintol, t3 / max (t1, t2) )
                      end if
                   end if

                end if

             enddo

c            ------------------------------------------------------------
c            ... general case:  all tested 2x2 pivots are unsatisfactory,
c                so test whether the diagonal an acceptable 1x1 pivot.
c                to emphasize selection of 2x2 pivots, delay tests for
c                1x1 pivots until after some number of columns have been
c                processed.
c            ------------------------------------------------------------
 
             if  ( jj .gt. triggr )  then

c               ---------------------------------------------------------
c               ... most general case -- test current column as 1x1 pivot
c               ---------------------------------------------------------

                gamma  = lmaxvl (j)

                if ( diag2 .gt. pvttol * gamma ) then

c                  ... accept this column as a 1x1 pivot

                   jpiv = j
                   go to 100

                end if

             else
     1       if  ( jj .eq. triggr )  then

c               ----------------------------------------------------
c               ... test preceding columns as well as current column
c               ----------------------------------------------------

c$$$                do 30 kk = m, triggr-1
                do kk = m, triggr

                   k      = corder(kk)
                   gamma  = lmaxvl (k)

                   if ( abs ( panel(k,k) ) .gt. pvttol * gamma ) then

c                     ... accept this column as a 1x1 pivot

                      jpiv = k
                      go to 100

                   end if

                enddo

             end if

c            -----------------------------------------------
c            ... otherwise, fail through to test next column
c            -----------------------------------------------

          else

c            ---------------------------------------------------------
c            ... process special case with nothing but exact zeros off
c                the main diagonal; also catches special case of
c                singleton root front
c            ---------------------------------------------------------

             jpiv = j

             if  ( diag2 .gt. zero )  then

c               ... non-singular: eliminate as 1x1,
c                   but skip rank 1 modification

                go to 150

             else

c               ... the matrix is exactly singular

                if  ( zpcntl (1) .gt. zero )  then

c                  ... the user wants us to fudge the diagonal

                   panel(j,j) = zpcntl(1)
                   if  ( zpcntl(2) .eq. zero )  then
                      zpcntl(3) = m
                   end if
                   zpcntl(2)  = zpcntl(2) + ( pnlcol - m + 1 )
                   zpcntl(4)  = zpcntl(1)
                   zpcntl(5)  = zpcntl(1)
                   go to 150

                else

c                  ... time to quit with an error return
c                      because this matrix is exactly singular

c                  write (*,*) 'EXACTLY SINGULAR MATRIX!'
                   ierr  = -1
                   ncolf = m - 1
                   go to 500

                end if

             end if

          end if

   40 continue

c     -------------------------------------------------------------
c     ... no acceptable pivot found among the uneliminated columns.
c     -------------------------------------------------------------

c     ... Remaining columns are not identically zero, but have no
c         satisfactory pivot within the current panel.  Postpone their
c         elimination to the next panel or to parent front unless we are
c         in the special case in which we succeeded in eliminating
c         no columns from the root front.

      if  ( m .gt. 1 )  then

c        ... although we were unsuccessful in eliminating every column
c            in the panel, we succeeded in eliminating something.  go
c            on to the next panel

         ncolf = m - 1
         go to 300

      else

c        ... we are in the special case in which we were unable to find
c            any acceptable pivot in the entirety of the current panel
c            In release 4.0 and 4.1, this becomes a fatal error condition
c            if the panel in question is a panel of a root block

         if  ( pmonit )  then

c           ... diagnostic output requested

            lg2pcol = 1
            pcol    = 1
            psize   = pnlcol-m+1

            do j = 2, psboxs
               if  ( psize .ge. 2*pcol )  then
                  lg2pcol = j
                  pcol    = 2*pcol
               end if
            enddo

            if  ( mintol .gt. zero )  then

c              ... could have pivoted with a smaller pivot tolerance

               expon = min ( xpboxs, -int ( log10 ( mintol) ) + 1 )
               if  ( qroot )  then
                  rtpexp ( expon ) = rtpexp ( expon ) + 1
                  rtpsiz (lg2pcol) = rtpsiz (lg2pcol) + 1
               else
                  inpexp ( expon ) = inpexp ( expon ) + 1
                  inpsiz (lg2pcol) = inpsiz (lg2pcol) + 1
               end if

            else

c              ... this panel has an exactly zero leading diagonal block

               if  ( qroot )  then
                  rzfail          = rzfail + 1
                  rtzsiz (lg2pcol) = rtzsiz (lg2pcol) + 1
               else
                  izfail          = izfail + 1
                  inzsiz (lg2pcol) = inzsiz (lg2pcol) + 1
               end if

            endif

         end if

         if  ( .not. qroot  )  then

c           ... move to next front

            ncolf = 0
            go to 500

         else

            if ( zpcntl (1) .gt. zero )  then

c              ------------------------------------------------------
c              ... zero pivot controls have been activated and it is
c                   the root front.  replace the remaining diagonal
c                   entries with zpcntl(1).  then take column m as
c                   a pivot.  THIS IS A HACK.  This does not provide
c                   any guarantees of stability.  It may be that this
c                   should be monitored separately from the monitoring
c                  of identically zero columns in main loop.
c               ------------------------------------------------------

               do j = m, pnlcol
                  t2         = abs ( panel(j,j) ) + ttmax * zpcntl(1)
                  panel(j,j) = sign ( t2, panel(j,j) )
               enddo

               if  ( zpcntl(2) .eq. zero )  then
                  zpcntl(3) = m
               end if
               zpcntl(2) = zpcntl(2) + ( pnlcol - m + 1 )
               zpcntl(4) = zpcntl(1)
               zpcntl(5) = zpcntl(1)

               jpiv = m
               go to 100

            else

c              ... the factorization breaks down because we have
c                  presently no continuation when we fail entirely
c                  on a panel in the root front.

               ierr  = -2
               ncolf = 0
               go to 500

            end if

         end if

c        ... cannot reach this branch

      end if

c     --------------------------------------------------------------
c     ... internal subroutines to perform the actual numerical
c         elimination steps, either  1x1  or  2x2  block elimination
c     --------------------------------------------------------------

c     ----------------------
c     ... 1x1 pivot selected
c     ----------------------

  100 continue
      pmonit = .false.
c.debug
c     write(6,'("after 100 - m, jpiv   = ", 3i8)')
c    1                        m, jpiv   
c     write(6,'("invp(m), invp(jpiv)   = ",4i8)')
c    1            invp(m), invp(jpiv)
c     call xdslp5('panel', pnlrow*pnlcol, panel, 6 )
c.debug

c     ----------------------------------
c     ... interchange columns m and jpiv
c     ----------------------------------
 
      if ( m .ne. jpiv ) call xdsle6 ( m, jpiv, 
     1                                 panel, pnlrow, pnlcol )
c.debug
c     write(6,'("after first xdsle6")')
c.debug

      j          = invp(m)
      invp(m   ) = invp(jpiv)
      invp(jpiv) = j

      swap  (m) = jpiv + offset
 
      pvtblk(m) = 1

c     --------------------------------
c     ... update the remaining columns
c     --------------------------------
 
c.debug
c     write(6,'("before xdsle7")')
c     call xdslp5('panel', pnlrow*pnlcol, panel, 6 )
c.debug
      call xdsle7  ( m, panel, pnlrow, pnlcol, inrtia )
c.debug
c     write(6,'("after  xdsle7")')
c     call xdslp5('panel', pnlrow*pnlcol, panel, 6 )
c.debug

      ll     = pnlrow - m + 1
      lc     = pnlcol - m + 1
      slvops = slvops + 4 * ( ll - 1 ) + 1
      fctops = fctops + ll + ( 2*ll-lc-2 ) * ( lc - 1 )
 
c     --------------------------
c     ... set up for next column
c     --------------------------
 
c.debug
c     write(6,'("before branch back to 10")')
c.debug
 
      if  ( m .lt. pnlcol )  then

c        ... return to check later columns

         m     = m + 1
         jpost = jpiv - m + 1
         go to 10

      else

c        ... all columns have been eliminated

         ncolf = pnlcol
         go to 300

      end if
 
c     ------------------------------------------------------------------

c     --------------------------------------------------------------
c     ... special 1x1 pivot selected (offdiagonals are exactly zero)
c     --------------------------------------------------------------

  150 continue
      pmonit = .false.
c.debug
c     write(6,'("after 100 - m, jpiv   = ", 3i8)')
c    1                        m, jpiv
c     write(6,'("invp(m), invp(jpiv)   = ",4i8)')
c    1            invp(m), invp(jpiv)
c     call xdslp5('panel', pnlrow*pnlcol, panel, 6 )
c.debug

c     ----------------------------------
c     ... interchange columns m and jpiv
c     ----------------------------------
 
      if ( m .ne. jpiv ) call xdsle6 ( m, jpiv, 
     1                                 panel, pnlrow, pnlcol )
c.debug
c     write(6,'("after first xdsle6")')
c.debug

      j          = invp(m)
      invp(m   ) = invp(jpiv)
      invp(jpiv) = j

      swap  (m) = jpiv + offset
 
      pvtblk(m) = 1

c     --------------------------------------------------------
c     ... skip applying low rank modification to the remaining 
c         columns because only the diagonal entry is nonzero
c     --------------------------------------------------------
 
      if ( panel(m,m) .gt. zero ) then
          inrtia(1) = inrtia(1) + 1
      else
          inrtia(2) = inrtia(2) + 1
      end if
 
c     --------------------------
c     ... set up for next column
c     --------------------------

c.debug
c     write(6,'("before branch back to 10")')
c.debug

      if  ( m .lt. pnlcol )  then

c        ... return to check later columns

         m     = m + 1
         jpost = jpiv - m + 1
         go to 10

      else

c        ... all columns have been eliminated

         ncolf = pnlcol
         go to 300

      end if

c     ------------------------------------------------------------------

c     ----------------------
c     ... 2x2 pivot selected
c     ----------------------
c
  200 continue
      pmonit = .false.
c.debug
c     write(6,'("after 200 - m, jpiv, kpiv = ", 3i8)')
c    1                        m, jpiv, kpiv   
c     call xdslp5('panel', pnlrow*pnlcol, panel, 6 )
c.debug

c     ----------------------------------
c     ... interchange columns m and jpiv
c     ----------------------------------
 
      if ( m .ne. jpiv ) call xdsle6 ( m, jpiv, 
     1                                 panel, pnlrow, pnlcol )
c.debug
c     write(6,'("after first xdsle6")')
c     call xdslp5('panel', pnlrow*pnlcol, panel, 6 )
c.debug

      j          = invp(m)
      invp(m   ) = invp(jpiv)
      invp(jpiv) = j
 
c     ------------------------------------
c     ... interchange columns m+1 and kpiv
c     ------------------------------------
 
      if ( m+1 .ne. kpiv ) call xdsle6 ( m+1, kpiv, 
     1                                   panel, pnlrow, pnlcol )
c.debug
c     write(6,'("after second xdsle6")')
c     call xdslp5('panel', pnlrow*pnlcol, panel, 6 )
c.debug

      j          = invp(m+1)
      invp(m+1 ) = invp(kpiv)
      invp(kpiv) = j
 
c     ---------------------------------------------
c     ... compute two columns of the factorization
c     ---------------------------------------------

      swap  (m  ) = jpiv + offset
      swap  (m+1) = kpiv + offset
 
      pvtblk (m)   = 2
      pvtblk (m+1) = 2
 
c     --------------------------------
c     ... update the remaining columns
c     --------------------------------

      call xdsle8  ( m, panel, pnlrow, pnlcol, inrtia )
c.debug
c     write(6,'("after  xdsle8")')
c     call xdslp5('panel', pnlrow*pnlcol, panel, 6 )
c.debug
 
      ll     = pnlrow - m + 1
      lc     = pnlcol - m + 1
      slvops = slvops + 8 * ( ll - 2 ) +  6
      fctops = fctops + 18 + 6 * ( ll - 2 )            
     1                + 2 * ( lc - 2 ) * ( 2*ll - lc - 3 )

c     --------------------------
c     ... set up for next column
c     --------------------------
 
c.debug
c     write(6,'("before branch back to 10")')
c.debug

      if  ( m .le. pnlcol-2 )  then

c        ... return to check other columns

         m     = m + 2
         jpost = max ( jpiv, kpiv ) - m + 1
         go to 10

      else

c        ... fall through if there are no more columns to be eliminated

         ncolf = pnlcol

      end if

c     ------------------------------------------------------------------
 
c     ---------------------------------------------------------------
c     ... set pointers in pvtblk to location of where the offdiagonal
c         entry of the 2x2 block will be in lnz.
c     ---------------------------------------------------------------

  300 continue
      m = 1

  400 continue
      if ( m .le. ncolf ) then

        if ( pvtblk(m) .eq. 2 ) then
            pvtblk(m+1) = - ( 1 + (m-1)*(ncolf-1) - (m-1)*(m-2)/2 )
            m = m + 2
        else
            m = m + 1
        end if

        go to 400
 
      end if
c     ------------------------------------------------------------------
 
  500 continue

      return
      end
