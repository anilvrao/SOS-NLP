      subroutine xdsld5 ( pvttol, zpcntl, qroot , pnlrow, pnlcol,
     1                    pospon, offset, invp  , panell, panelu,
     2                    pvtblk, corder, temp1 , swap  , ncolf ,  
     3                    fctops, slvops, ppfmon, xpboxs, psboxs,
     4                    inpexp, inpsiz, inzsiz, izfail, rtpexp,
     5                    rtpsiz, rtzsiz, rzfail, ierr )

 
c
c  purpose -- perform a block elimination step, eliminating as many
c             columns from a panel of a front as possible.  We use
c             threshold duff-reid pivoting for stability.  The pivoting
c             is exhaustive only within the panel.  Any columns that are
c             not eliminated do not have stable pivots within the panel,
c             but may have satisfactory pivots within the remainder of
c             the front.
c
c             this is the column oriented panel version for
c             unsymmetric matrices
c
c  created            -- 22-may-98, rgg, derived form xdsle5
c  last modified      -- 27-aug-98, rgg, added zero pivot controls
c                     -- 04-oct-01, jgl, restructured to repair test for
c                                        singularity; error codes
c                                        changed
c                     -- 22-mar-02, jgl, repair final low rank 
c                                        modification when we find
c                                        no eliminatable columns
c               
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
c      panell -- panel of the current frontal matrix to be eliminated
c                lower triangle
c      panelu -- panel of the current frontal matrix to be eliminated
c                upper triangle
c      pvtblk -- integer vector indicates if diagonal element is a
c                1 x 1 or 2 x 2 block pivot.
c      fctops -- factor operation count
c      slvops -- solve operation count
c      ppfmon -- monitor pivot failures or not
c      xpboxs -- number of pivot exponent boxes for pivot monitoring
c      psboxs -- number of panel size boxes for pivot monitoring
c
c  working storage --
c
c      corder -- temporary array of length pnlrow
c      temp1  -- temporary array of length pnlrow
c
c  output variable --
c
c      invp   -- inverse permutation array.
c      swap   -- swapping record
c      panel  -- panel of the current frontal matrix to be eliminated
c      ncolf  -- the number of columns that were eliminated.  
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
     1                   pvtblk (*),     
     2                   swap   (*)

      logical            qroot
 
      double precision   pvttol, fctops, slvops
 
      double precision   panell(pnlrow,pnlcol),  
     1                   panelu(pnlrow,pnlcol),  
     2                   temp1 (*), zpcntl(*)
 
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
 
      integer            bstart, bend  , bmax  , expon ,  j    , jpiv  , 
     1                   k     , l     , lc    , lg2pcol, ll   , m     ,
     2                   pcol  , psize

      logical            qnomor, pmonit

      double precision   diag1 , fac   , rdiag1, tmax  , ttmax ,
     1                   t2    , zero  , mintol

      parameter        ( zero = 0.0 )
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      double precision   damax
 
      external           damax, xdsld6, xdsld7
 
c     ------------------------------------------------------------------

      ierr   = 0
      m      = 1
      bstart = 1
      bend   = 0
      bmax   = 2
      qnomor = .false.
      pmonit = ppfmon

c     ----------------------------
c     ... find an acceptable pivot
c     ----------------------------
 
   10 continue
 
c     -----------------------------------------------------
c     ... loop through remaining columns searching
c         for an acceptable pivot
c     -----------------------------------------------------

      ttmax  = 1.
      mintol = 0.
 
      do 40 j = m, pnlcol
c.debug
c     write(6,'("in xdsld5 do 40 loop for j, bstart, bend = ", 3i8)') 
c    1                                     j, bstart, bend
c     call xdslp5('panell at start of do 40', pnlcol*pnlrow, panell, 6 )
c     call xdslp5('panelu at start of do 40', pnlcol*pnlrow, panelu, 6 )
c.debug

c         ----------------------------
c         ... copy column j into temp1
c         ----------------------------

          temp1(m:j-1) = panelu(j,m:j-1)

          temp1(j:pnlrow) = panell(j:pnlrow,j)
c.debug
c     call xdslp5('temp1 after extract', pnlrow-m+1, temp1(m), 6 )
c.debug

c         ---------------------------------------------
c         ... update temp1 with columns not yet applied
c         ---------------------------------------------

          l = pnlrow - m + 1

          do k = bstart, bend
              fac = -panelu(j,k) * panell(k,k) 
              do kk=m,m+l-1
                temp1(kk) = temp1(kk) + fac*panell(kk,k)
              enddo
          enddo
c.debug
c     call xdslp5('temp1 after update', pnlrow-m+1, temp1(m), 6 )
c.debug
 
c         --------------------------------------------
c         ... is the diagonal an acceptable 1x1 pivot?
c         --------------------------------------------

          diag1 = abs ( temp1(j) )
 
          tmax = max ( damax ( j-m     , temp1(m  ), 1 ),
     1                 damax ( pnlrow-j, temp1(j+1), 1 ) )

          ttmax = max ( ttmax, tmax )
 
          if ( diag1 .gt. zero ) then
              if ( diag1 .ge. pvttol * tmax ) then
                  jpiv = j
                  go to 100
              else
                  if  ( pmonit )  then
                     mintol = max ( mintol, diag1 / tmax )
                  end if
              end if
          else
              if ( tmax .eq. zero ) then
                  ierr = -1
                  return
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

      ncolf = m - 1

      if  ( ncolf .gt. 0 )  then

c        ... although we were unsuccessful in eliminating every column
c            in the panel, we succeeded in eliminating something.  go
c            on to the next panel

         qnomor = .true.
         go to 200

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

            qnomor = .true.
            go to 200

         else

            if ( zpcntl (1) .gt. zero )  then

c              -------------------------------------------------------
c              ... zero pivot controls have been activated and it is
c                   the root front.  replace the remaining diagonal
c                   entries with zpcntl(1).  then take column m as
c                   a pivot.  THIS IS A HACK.  This does not provide
c                   any guarantees of stability.  It may be that this
c                   should be monitored separately from the monitoring
c                  of identically zero columns in main loop.
c              -------------------------------------------------------

               do j = m, pnlcol
                   t2          = abs ( panell(j,j) ) + ttmax * zpcntl(1)
                   panell(j,j) = sign ( t2, panell(j,j) ) 
               enddo

               l = pnlrow - m + 1

               temp1(m:m+l-1) = panell(m:m+l-1,m)

               do k = bstart, bend
                   fac = -panelu(m,k) * panell(k,k) 
                   do kk=m,m+l-1
                     temp1(kk) = temp1(kk) + fac*panell(kk,k)
                   enddo
               enddo

               if ( zpcntl(2) .eq. zero ) zpcntl(3) = m
               zpcntl(2) = zpcntl(2) + ( pnlcol - m + 1 )
               zpcntl(4) = zpcntl(1)
               zpcntl(5) = zpcntl(1)

               jpiv = m
               go to 100

            else

c              ... the factorization breaks down because we have
c                  presently no continuation when we fail entirely
c                  on a panel in the root front.

               ierr = -2
               return

            end if

         end if

c        ... cannot reach this branch

      end if
 
c     ------------------------------------------------------------------
 
c     ----------------------
c     ... 1x1 pivot selected
c     ----------------------
 
  100 continue
      pmonit = .false.
 
c     --------------------------------------------------------
c     ... interchange columns m and jpiv.
c         also replace m-th column of a with contents of temp1
c     --------------------------------------------------------
 
      call xdsld6 ( m, jpiv, temp1, panell, panelu, pnlrow, pnlcol )
c.debug
c     write(6,'("interchanging m, jpiv, offset = ", 3i8)') 
c    1      m, jpiv, offset
c     call xdslp5('panell after interchange', pnlcol*pnlrow, panell, 6 )
c     call xdslp5('panelu after interchange', pnlcol*pnlrow, panelu, 6 )
c.debug
 
      j          = invp(m)
      invp(m   ) = invp(jpiv)
      invp(jpiv) = j

      swap  (m) = jpiv + offset
 
      pvtblk(m) = 1
c.debug
c     write(6,'("m,jpiv,invp(m),invp(jpiv),swap(m) = ", 5i8)')
c    1            m,jpiv,invp(m),invp(jpiv),swap(m)
c.debug

c     -------------------------------
c     ... apply updates to row m of u
c     -------------------------------

      do k = bstart, bend

          fac = -panell(m,k) * panell(k,k) 
c.debug
c     write(6,'("in do 110 k, l, fac = ", 2i8,1pd15.5)') k, l, fac 
c.debug
          do kk=m+1,m+l-1
            panelu(kk,m) = panelu(kk,m) + fac*panelu(kk,k)
          enddo
      enddo

c     ---------------------------------
c     ... form column of l and row of u
c     ---------------------------------

      rdiag1 = 1.0 / panell(m,m)

      do k = m+1, pnlrow
          panell(k,m) = rdiag1 * panell(k,m)
          panelu(k,m) = rdiag1 * panelu(k,m)
      enddo
c.debug
c     write(6,'("after forming column and row")') 
c     call xdslp5('panell ', pnlcol*pnlrow, panell, 6 )
c     call xdslp5('panelu ', pnlcol*pnlrow, panelu, 6 )
c.debug

      ll     = pnlrow - m + 1
      lc     = pnlcol - m + 1
      slvops = slvops + 4 * ( ll - 1 ) + 1
      k      = ll + ( 2*ll-lc-2 ) * ( lc - 1 )
      fctops = fctops + 2 * k

      bend  = m

      if ( m .eq. pnlcol ) then
          qnomor = .true.
          ncolf  = pnlcol
      end if
 
c     ------------------------------------------------
c     ... test to see if it is time to do block update
c     ------------------------------------------------

  200 continue
      k = bend - bstart + 1 

      if ( ( qnomor .and. k .gt. 0 )  .or. k .eq. bmax ) then
c.debug
c     write(6,'("before block update-bstart, bend, k, bmax, qnomor= ",
c    1                     4i8,l8)')  bstart, bend, k, bmax, qnomor
c     write(6,'("before block update-m                            = ",
c    1                     4i8,l8)')  m                            
c.debug

          call xdsld7 ( bstart, bend, pnlrow, pnlcol, panell, panelu )

          bstart = m + 1
c.debug
c     write(6,'("after block update")') 
c     call xdslp5('panell ', pnlcol*pnlrow, panell, 6 )
c     call xdslp5('panelu ', pnlcol*pnlrow, panelu, 6 )
c.debug

      end if
 
c     --------------------------
c     ... set up for next column
c     --------------------------
 
      if ( .not. qnomor ) then
          m = m + 1 
          go to 10
      end if

c     ------------------------------------------------------------------
 
      return
      end
