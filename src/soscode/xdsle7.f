      subroutine xdsle7 ( m, panel, pnlrow, pnlcol, inrtia )
 
c
c  purpose -- to determine the factorization of the m-th column of the
c             current panel of the current front and to update the 
c             remainder of the panel.
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
 
      double precision  diag1 , fac11 , fac12 , rdiag1
 
c  =====================================================================
 
c     ---------------------------
c     scale the column and update
c     ---------------------------
 
      diag1 = panel(m,m)
 
      if ( diag1 .gt. 0.0 ) then
          inrtia(1) = inrtia(1) + 1
      else
          inrtia(2) = inrtia(2) + 1
      end if
 
      rdiag1 = 1. / diag1
 
cdir$ ivdep
      do i = m+1, pnlrow
          panel(i,m) = rdiag1*panel(i,m)
      enddo

c     -----------------------------------------------
c     ... perform the update with some loop unrolling
c     -----------------------------------------------
 
      do j = m+1, pnlcol, 2
c.debug
c     write(6,'("in do 40 in xdsle7")')
c     write(6,'("j, m, pnlcol = ", 3i8)') j, m, pnlcol
c.debug

          if ( j .lt. pnlcol ) then 
 
              fac11 = diag1*panel(j  ,m)
              fac12 = diag1*panel(j+1,m)

              panel(j,j) = panel(j,j) - fac11*panel(j,m)
 
cdir$ ivdep
              do i = j+1, pnlrow 
                  panel(i,j  ) = panel(i,j  ) - fac11*panel(i,m)
                  panel(i,j+1) = panel(i,j+1) - fac12*panel(i,m)
              enddo

          else  

              fac11 = diag1*panel(j  ,m)
 
cdir$ ivdep
              do i = j, pnlrow
                  panel(i,j) = panel(i,j) - fac11*panel(i,m)
              enddo

          end if
 
      enddo
 
      return
 
c  =====================================================================
 
      end
