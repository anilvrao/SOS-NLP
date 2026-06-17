      subroutine xdslee ( m, panel, pnlrow, pnlcol, inrtia )
 
c
c  purpose -- to determine the factorization of the m-th and the
c             m+1 th columns of the current panel of the current 
c             front and to update the remainder of the panel.
c
c             this is the row oriented panel version
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
 
      double precision  panel (pnlcol,pnlrow)
 
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
      offdia = panel(m  ,m+1)
c.debug
c     write(6,'("in xdslee - diag1, diag2, offdia = ", 1p3d15.5)')
c    1                        diag1, diag2, offdia 
c.debug
      hold   = diag1 * diag2 - offdia ** 2
c.debug
c     write(6,'("in xdslee - hold                 = ", 1p3d15.5)')
c    1                        hold                 
c.debug
 
c     ---------------
c     ... set inertia
c     ---------------
 
      fac11  = ( diag1 + diag2 ) / 2.
c.debug
c     write(6,'("in xdslee - fac11                = ", 1p3d15.5)')
c    1                        fac11                
c.debug
 
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
c.debug
c     write(6,'("in xdslee - dinv(1,1), dinv(2,2) = ", 1p3d15.5)')
c    1                        dinv(1,1), dinv(2,2) 
c     write(6,'("in xdslee - dinv(2,1)            = ", 1p3d15.5)')
c    1                        dinv(2,1)            
c.debug
 
cdir$ ivdep
      do i = m+2, pnlrow
          hold         = panel(m,i)
          panel(m,i  ) = hold * dinv (1, 1) + panel(m+1,i) * dinv (2, 1)
          panel(m+1,i) = hold * dinv (1, 2) + panel(m+1,i) * dinv (2, 2)
      enddo
c.debug
c     write(6,'("after 10 loop")')
c.debug
 
c     ------------------------------------------------
c     update the remaining columns with loops unrolled
c     ------------------------------------------------
 
      do 50 j = m+2, pnlcol, 3

c         -----------------------
c         ... set up first column
c         -----------------------
 
          fac11 = diag1  * panel(m,j  ) + offdia * panel(m+1,j  )
          fac21 = offdia * panel(m,j  ) + diag2  * panel(m+1,j  )
c.debug
c     write(6,'("j, m, pnlcol                     = ", 3i8)')
c    1            j, m, pnlcol
c     write(6,'("in do 50  - fac11, fac21         = ", 1p3d15.5)')
c    1                        fac11, fac21         
c.debug

          panel(j  ,j  ) = panel(j  ,j  ) - fac11*panel(m  ,j  ) 
     1                                    - fac21*panel(m+1,j  )

          if ( j .eq. pnlcol ) then

              do i = j+1, pnlrow 
                  panel(j  ,i) = panel(j  ,i) - fac11*panel(m  ,i) 
     1                                        - fac21*panel(m+1,i)
              enddo
c.debug
c     write(6,'("after 20 loop")')
c.debug

              go to 50

          end if

c         ------------------------
c         ... set up second column
c         ------------------------

          fac12 = diag1  * panel(m,j+1) + offdia * panel(m+1,j+1)
          fac22 = offdia * panel(m,j+1) + diag2  * panel(m+1,j+1)
c.debug
c     write(6,'("in do 50  - fac12, fac22         = ", 1p3d15.5)')
c    1                        fac12, fac22         
c.debug

          panel(j  ,j+1) = panel(j  ,j+1) - fac11*panel(m  ,j+1) 
     1                                    - fac21*panel(m+1,j+1)
          panel(j+1,j+1) = panel(j+1,j+1) - fac12*panel(m  ,j+1) 
     1                                    - fac22*panel(m+1,j+1)
c.debug
c     write(6,'("before if test for 30 loop")')
c     write(6,'("j, pnlcol = ", 2i8)') j, pnlcol
c.debug

          if ( j+1 .eq. pnlcol ) then
c.debug
c     write(6,'("before 30 loop")')
c.debug

              do i = j+2, pnlrow 
                  panel(j  ,i) = panel(j  ,i) - fac11*panel(m  ,i) 
     1                                        - fac21*panel(m+1,i)
                  panel(j+1,i) = panel(j+1,i) - fac12*panel(m  ,i) 
     1                                        - fac22*panel(m+1,i)
              enddo
c.debug
c     write(6,'("after 30 loop")')
c.debug

              go to 50

          end if

c         -----------------------
c         ... set up third column
c         -----------------------

          fac13 = diag1  * panel(m,j+2) + offdia * panel(m+1,j+2)
          fac23 = offdia * panel(m,j+2) + diag2  * panel(m+1,j+2)
c.debug
c     write(6,'("in do 50  - fac13, fac23         = ", 1p3d15.5)')
c    1                        fac13, fac23         
c.debug

c         ---------------------------
c         ... main computational loop
c         ---------------------------
 
cdir$ ivdep
          do i = j+2, pnlrow
              panel(j  ,i) = panel(j  ,i) - fac11*panel(m  ,i) 
     1                                    - fac21*panel(m+1,i)
              panel(j+1,i) = panel(j+1,i) - fac12*panel(m  ,i) 
     1                                    - fac22*panel(m+1,i)
              panel(j+2,i) = panel(j+2,i) - fac13*panel(m  ,i) 
     1                                    - fac23*panel(m+1,i)
          enddo
c.debug
c     write(6,'("after 40 loop")')
c.debug

c         -------------------------------------------
c         ... loop back for next set of three columns
c         -------------------------------------------

   50 continue
 
      return
 
c  =====================================================================
 
      end
