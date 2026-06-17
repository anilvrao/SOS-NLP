      subroutine xdsled ( m, panel, pnlrow, pnlcol, inrtia )
 
c
c  purpose -- to determine the factorization of the m-th column of the
c             current panel of the current front and to update the 
c             remainder of the panel.
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
          panel(m,i) = rdiag1*panel(m,i)
      enddo

c     -----------------------------------------------
c     ... perform the update with some loop unrolling
c     -----------------------------------------------
 
      do j = m+1, pnlcol, 2

          if ( j .lt. pnlcol ) then 
 
              fac11 = diag1*panel(m,j  )
              fac12 = diag1*panel(m,j+1)

              panel(j,j) = panel(j,j) - fac11*panel(m,j)
 
cdir$ ivdep
              do i = j+1, pnlrow
                  panel(j  ,i) = panel(j  ,i) - fac11*panel(m,i)
                  panel(j+1,i) = panel(j+1,i) - fac12*panel(m,i)
              enddo

          else  

              fac11 = diag1*panel(m,j)
 
cdir$ ivdep
              do i = j, pnlrow
                  panel(j,i) = panel(j,i) - fac11*panel(m,i)
              enddo

          end if
 
      enddo
 
      return
 
c  =====================================================================
 
      end
