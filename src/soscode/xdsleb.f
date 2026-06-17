      subroutine xdsleb (pvttol, zpcntl, qroot ,
     1                   pnlrow, pnlcol, pospon, offset,
     1                   invp  , panel , pvtblk, corder, temp1 , 
     2                   swap  , ncolf , inrtia, fctops, slvops, 
     3                   ierr  )

c
c  purpose -- to determine the factorization of the column(s) of the
c             node(s) in a supernode and to update the frontal matrix
c             using these factored column(s).
c             threshold bunch-kaufmann pivoting used for stability
c
c             this is the row oriented panel version
c
c  created            -- 21-feb-97, rgg
c                        27-aug-98, rgg, added zero pivot controls
c                                        for root front
c
c  input variables --
c
c      pvttol -- pivoting tolerance
c      zpcntl -- zero pivot control structure
c      qroot  -- logical flag that is true if root front
c      pnlrow -- number of rows in the panel 
c      pnlcol -- number of columns in the panel
c      pospon -- number of postponed columns in this front
c      offset -- swapping recording offset
c      invp   -- inverse permutation array.
c      panel  -- panel of the current frontal matrix to be eliminated
c      pvtblk -- integer vector indicates if diagonal element is a
c                1 x 1 or 2 x 2 block pivot.
c      inrtia -- matrix inertia
c      fctops -- factor operation count
c      slvops -- solve operation count
c
c  working storage --
c
c      corder -- temporary array of length lfront
c      temp1  -- temporary array of length lfront
c
c  output variable --
c
c      invp   -- inverse permutation array.
c      panel  -- panel of the current frontal matrix to be eliminated
c      swap   -- swapping record
c      ncolf  -- the number of columns that were eliminated.  
c      inrtia -- matrix inertia
c      fctops -- factor operation count
c      slvops -- solve operation count
c      ierr   -- error return,
c                if ierr = 0, success,
c                     .ne. 0, zero pivot detected in a 2x2 pivot
c
c     ------------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer            pnlrow, pnlcol, pospon, ncolf , 
     1                   offset, ierr
 
      integer            corder (*),     invp (*),
     1                   pvtblk (*),     inrtia (3),
     2                   swap   (*)
 
      logical            qroot
 
      double precision   pvttol, fctops, slvops
 
      double precision   panel (pnlcol,pnlrow),  temp1 (*),
     1                   zpcntl(*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer            j     , jpiv  , kpiv  , 
     1                   k1    , triggr,
     2                   locf  , m     , 
     3                   jj    , kk    , lc    , ll    , jpost 

      double precision   diag1 , f     , gamma , lambda, 
     1                   t1    , t2    , t3    , diag2 , ttmax
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      external  xdslec, xdsled, xdslee
 
c     ------------------------------------------------------------------

      ierr   = 0
      m      = 1
      ncolf  = pnlcol
      jpost  = 0

c     ----------------------------
c     ... find an acceptable pivot
c     ----------------------------
 
   10 continue
      kpiv   = 0

c     ----------------------------------------------------------------
c     ... search postponed columns last
c         after the first pivot is chosen, start after the last column
c         chosen
c     ----------------------------------------------------------------

      do j = m + jpost, pnlcol
          corder(j-jpost) = j
      enddo

      do j = m, m+jpost-1
          corder(pnlcol-(m+jpost-1)+j) = j
      enddo

      jpost = 0

c     -----------------------------------------------------
c     ... loop through remaining diagonal entries searching
c         for a 1x1 or 2x2 pivot
c     -----------------------------------------------------

      triggr = min ( m + 4, pnlcol )

      ttmax  = 1.
 
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
            gamma = max(gamma,abs(panel(k,j)))
          enddo

          do k=j,pnlrow
            gamma = max(gamma,abs(panel(j,k)))
          enddo

          temp1(j) = gamma 

c         ------------------------------------------------
c         ... for every off-diagonal entry in row j in
c             the interior of the front
c             -- check for 2x2 pivot
c         ------------------------------------------------

          do ll = m, jj-1
c.debug
c     write(6,'("in 25 - ll, m, jj     = ", 3i8)')
c    1                    ll, m, jj     
c.debug

              locf   = corder(ll)

              diag1  = abs ( panel(locf,locf) )

              lambda = temp1 ( locf )

              if ( j .gt. locf ) then 
                  f = abs ( panel(locf,j) )
              else
                  f = abs ( panel(j,locf) )
              end if

c             -----------------------------------------
c             ... choose 2x2 pivot?  use explicit bound
c             -----------------------------------------

              t1    = diag2 * lambda + f * gamma
              t2    = diag1 * gamma  + f * lambda
              t3    = abs ( diag1 * diag2 - f ** 2 )

              if ( t3 .gt. 0. .and. max (t1,t2) .le. t3/pvttol ) then 
 
                  jpiv   = min ( locf, j )
                  kpiv   = max ( locf, j )

c                 if ( jpiv .eq. m+1 ) then
c                     jpiv = kpiv
c                     kpiv = m+1
c                 end if

                  go to 200
 
              end if

          enddo

c         --------------------------------------------
c         ... is the diagonal an acceptable 1x1 pivot?
c             the first of these tests are delayed to 
c             emphasize 2x2 pivots.
c         --------------------------------------------
 
          if ( jj .lt. triggr ) go to 40

          if ( jj .eq. triggr ) then
 
              do kk = m, triggr-1

                  k1     = corder(kk)
                  gamma  = temp1(k1)

                  if ( abs ( panel(k1,k1) ) .gt. pvttol * gamma ) then

                      jpiv = k1

                      go to 100

                  end if

              enddo

          end if

c         -----------------------
c         ... test current column
c         -----------------------

          gamma  = temp1(j)

          if ( diag2 .gt. pvttol * gamma ) then
 
              jpiv = j

              go to 100
 
          end if

   40 continue
 
c     ------------------------------------------------------
c     ... no acceptable pivot found.  remainder of matrix is
c         identically zero.
c     ------------------------------------------------------

      if ( zpcntl(1) .eq. 0. .or. .not. qroot ) then

c         -------------------------------------------------------
c         ... postpone the rest of the columns to the next front.
c         -------------------------------------------------------

          ncolf = m - 1
          go to 300

      else

c         -----------------------------------------------------
c         ... zero pivot controls have been activated and it is
c             the root front.  replace the remaining diagonal
c             entries with zpcntl(1) and restart the search for 
c             a pivot.
c         -----------------------------------------------------

          do j = m, pnlcol
              t2         = abs ( panel(j,j) ) + ttmax * zpcntl(1)
              panel(j,j) = sign ( t2, panel(j,j) )
          enddo

          if ( zpcntl(2) .eq. 0. ) zpcntl(3) = m
          zpcntl(2) = zpcntl(2) + ( pnlcol - m + 1 )
          zpcntl(4) = zpcntl(1)
          zpcntl(5) = zpcntl(1)

          jpiv = m
          go to 100

      end if
 
c     ------------------------------------------------------------------
 
c     ----------------------
c     ... 1x1 pivot selected
c     ----------------------
 
  100 continue
c.debug
c     write(6,'("after 100 - m, jpiv   = ", 3i8)')
c    1                        m, jpiv   
c.debug

c     ----------------------------------
c     ... interchange columns m and jpiv
c     ----------------------------------
 
      if ( m .ne. jpiv ) call xdslec ( m, jpiv, 
     1                                 panel, pnlrow, pnlcol )
c.debug
c     write(6,'("after first xdslec")')
c.debug

      j          = invp(m)
      invp(m   ) = invp(jpiv)
      invp(jpiv) = j

      swap(m)    = jpiv + offset
 
      pvtblk(m) = 1

c     --------------------------------
c     ... update the remaining columns
c     --------------------------------
 
c.debug
c     write(6,'("before xdsled")')
c.debug
      call xdsled  ( m, panel, pnlrow, pnlcol, inrtia )
c.debug
c     write(6,'("after  xdsled")')
c.debug

      ll     = pnlrow - m + 1
      lc     = pnlcol - m + 1
      slvops = slvops + 4 * ( ll - 1 ) + 1
      fctops = fctops + ll + ( 2*ll-lc-2 ) * ( lc - 1 )
 
c     --------------------------
c     ... set up for next column
c     --------------------------
 
      m     = m + 1
      jpost = jpiv - m + 1
c.debug
c     write(6,'("before branch back to 10")')
c.debug
 
      if ( m .le. pnlcol ) go to 10
 
      go to 300
 
c     ------------------------------------------------------------------
 
c     ----------------------
c     ... 2x2 pivot selected
c     ----------------------
c
  200 continue
c.debug
c     write(6,'("after 200 - m, jpiv, kpiv = ", 3i8)')
c    1                        m, jpiv, kpiv   
c     call xdslp5('panel', pnlrow*pnlcol, panel, 6 )
c.debug

c     ----------------------------------
c     ... interchange columns m and jpiv
c     ----------------------------------
 
      if ( m .ne. jpiv ) call xdslec ( m, jpiv, 
     1                                 panel, pnlrow, pnlcol )
c.debug
c     write(6,'("after first xdslec")')
c     call xdslp5('panel', pnlrow*pnlcol, panel, 6 )
c.debug

      j          = invp(m)
      invp(m   ) = invp(jpiv)
      invp(jpiv) = j
 
c     ------------------------------------
c     ... interchange columns m+1 and kpiv
c     ------------------------------------
 
      if ( m+1 .ne. kpiv ) call xdslec ( m+1, kpiv, 
     1                                   panel, pnlrow, pnlcol )
c.debug
c     write(6,'("after second xdslec")')
c     call xdslp5('panel', pnlrow*pnlcol, panel, 6 )
c.debug

      j          = invp(m+1)
      invp(m+1 ) = invp(kpiv)
      invp(kpiv) = j
 
c     ---------------------------------------------
c     ... compute two columns of the factorization
c     ---------------------------------------------

      swap(m  )    = jpiv + offset
      swap(m+1)    = kpiv + offset
 
      pvtblk (m)   = 2
      pvtblk (m+1) = 2
 
c     --------------------------------
c     ... update the remaining columns
c     --------------------------------

      call xdslee  ( m, panel, pnlrow, pnlcol, inrtia )
c.debug
c     write(6,'("after  xdslee")')
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
 
      m     = m + 2
      jpost = max ( jpiv, kpiv ) - m + 1
c.debug
c     write(6,'("before branch back to 10")')
c.debug

      if ( m .le. pnlcol ) go to 10
 
c     ------------------------------------------------------------------
 
  300 continue

c     --------------------------------------------------------------
c     ... set pointer in pvtblk to location of where the offdiagonal
c         entry of the 2x2 block will be in lnz.
c     --------------------------------------------------------------

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
      endif
 
c     ------------------------------------------------------------------

      return
      end
