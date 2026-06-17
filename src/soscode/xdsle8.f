      subroutine xdsle8 ( m, panel, pnlrow, pnlcol, inrtia )
c
c  purpose -- to determine the factorization of the m-th and the
c             m+1 th columns of the current panel of the current 
c             front and to update the remainder of the panel.
c
c             this is the column oriented panel version
c
c  created            -- 21-feb-97, rgg
c  last modifications --
c
c  input variables --
c
c      m      -- column to be eliminated
c      pnlrow -- number of rows    in panel
c      pnlcol -- number of columns in panel
c
c  output variable --
c
c      panel  -- panel of frontal matrix
c      inrtia -- matrix inertia
c
c  =====================================================================
 
      integer           m     , pnlrow, pnlcol
 
      integer           inrtia(3)
 
      double precision  panel (pnlrow,pnlcol)
 
      integer           i     , j     
 
      double precision  diag1 , diag2 , hold  , offdia,
     1                  fac11 , fac21 , fac12 , fac22 , fac13 , fac23
 
      double precision  dinv(2,2)
 
c  =====================================================================
 
c     ----------------------------------------
c     compute the factorization of two columns
c     ----------------------------------------
 
      diag1  = panel(m  ,m  )
      diag2  = panel(m+1,m+1)
      offdia = panel(m+1,m  )
      hold   = diag1 * diag2 - offdia ** 2
 
c     ---------------
c     ... set inertia
c     ---------------
 
      fac11  = ( diag1 + diag2 ) / 2.
 
      if ( fac11 .gt. 0.  .and.  hold .gt. 0. ) then
          inrtia(1) = inrtia(1) + 2
      else if ( fac11 .lt. 0.  .and.  hold .gt. 0. ) then
          inrtia(2) = inrtia(2) + 2
      else
          inrtia(1) = inrtia(1) + 1
          inrtia(2) = inrtia(2) + 1
      end if
 
      dinv (1, 1) = diag2 / hold
      dinv (2, 2) = diag1 / hold
      dinv (1, 2) = -offdia / hold
      dinv (2, 1) = dinv (1, 2)
 
cdir$ ivdep
      do i = m+2, pnlrow
          hold         = panel(i,m)
          panel(i,m  ) = hold * dinv (1, 1) + panel(i,m+1) * dinv (2, 1)
          panel(i,m+1) = hold * dinv (1, 2) + panel(i,m+1) * dinv (2, 2)
      enddo
 
c     ------------------------------------------------
c     update the remaining columns with loops unrolled
c     ------------------------------------------------
 
      do 50 j = m+2, pnlcol, 3

c         -----------------------
c         ... set up first column
c         -----------------------
 
          fac11 = diag1  * panel(j  ,m) + offdia * panel(j  ,m+1)
          fac21 = offdia * panel(j  ,m) + diag2  * panel(j  ,m+1)

          panel(j  ,j  ) = panel(j  ,j  ) - fac11*panel(j  ,m  ) 
     1                                    - fac21*panel(j  ,m+1)

          if ( j .eq. pnlcol ) then

              do i = j+1, pnlrow 
                  panel(i,j  ) = panel(i,j  ) - fac11*panel(i,m  ) 
     1                                        - fac21*panel(i,m+1)
              enddo

              cycle

          end if

c         ------------------------
c         ... set up second column
c         ------------------------

          fac12 = diag1  * panel(j+1,m) + offdia * panel(j+1,m+1)
          fac22 = offdia * panel(j+1,m) + diag2  * panel(j+1,m+1)

          panel(j+1,j  ) = panel(j+1,j  ) - fac11*panel(j+1,m) 
     1                                    - fac21*panel(j+1,m+1)
          panel(j+1,j+1) = panel(j+1,j+1) - fac12*panel(j+1,m) 
     1                                    - fac22*panel(j+1,m+1)

          if ( j+1 .eq. pnlcol ) then

              do i = j+2, pnlrow 
                  panel(i,j  ) = panel(i,j  ) - fac11*panel(i,m  ) 
     1                                        - fac21*panel(i,m+1)
                  panel(i,j+1) = panel(i,j+1) - fac12*panel(i,m  ) 
     1                                        - fac22*panel(i,m+1)
              enddo

              cycle

          end if

c         -----------------------
c         ... set up third column
c         -----------------------

          fac13 = diag1  * panel(j+2,m) + offdia * panel(j+2,m+1)
          fac23 = offdia * panel(j+2,m) + diag2  * panel(j+2,m+1)

c         ---------------------------
c         ... main computational loop
c         ---------------------------
 
cdir$ ivdep
          do i = j+2, pnlrow
              panel(i,j  ) = panel(i,j  ) - fac11*panel(i,m  ) 
     1                                    - fac21*panel(i,m+1)
              panel(i,j+1) = panel(i,j+1) - fac12*panel(i,m  ) 
     1                                    - fac22*panel(i,m+1)
              panel(i,j+2) = panel(i,j+2) - fac13*panel(i,m  ) 
     1                                    - fac23*panel(i,m+1)
          enddo

c         -------------------------------------------
c         ... loop back for next set of three columns
c         -------------------------------------------

   50 continue
 
      return
 
c  =====================================================================
 
      end
