      subroutine xdsld7 ( bstart, bend, pnlrow, pnlcol, panell, panelu )

 
c
c  purpose -- to determine the factorization of the first column of the
c             current front and to update the frontal matrix using the
c             factored column.
c
c             this routine updates the reduced matrix two columns at a
c             time.
c
c             this subroutine was derived from updb32.
c
c             unsymmetric version
c
c  created            -- 22-may-98, rgg
c  last modifications --
c
c  input variables --
c
c      bstart -- first column to apply
c      bend   -- last column to apply
c      pnlrow -- number of rows in panel
c      pnlcol -- number of columns in panel
c
c  output variable --
c
c      panell -- panel for lower tri of the factor
c      panelu -- panel for upper tri of the factor
c
c  =====================================================================
 
      integer           bstart, bend, pnlrow, pnlcol
 
      double precision  panell(pnlrow,pnlcol), panelu(pnlrow,pnlcol)
 
      integer           i     , j     , m
 
      double precision  faca11, facu11, faca21, facu21
 
c  =====================================================================

      do 140 m = bstart, bend, 2
c.debug
c     write(6,'("xdsld7-140-m,bstart,bend,pnlrow,pnlcol =", 5i8)')
c    1                       m,bstart,bend,pnlrow,pnlcol        
c.debug

          do j = bend+1, pnlcol

              faca11 = panell(m  ,m  ) * panell(j,m  )
              facu11 = panell(m  ,m  ) * panelu(j,m  )
c.debug
c     write(6,'("xdsld7-130-m,j,pnlrow,pnlcol      =", 4i8)')
c    1                       m,j,pnlrow,pnlcol        
c     write(6,'("xdsld7-130-faca11, facu11         =", 1p2d15.5)')
c    1                       faca11, facu11         
c.debug

              panell(j,j) = panell(j,j) - facu11 * panell(j,m  )

              if ( m .eq. bend ) then

cdir$ ivdep
                  do i = j+1, pnlrow
                      panell(i,j) = panell(i,j) - facu11 * panell(i,m  )
                      panelu(i,j) = panelu(i,j) - faca11 * panelu(i,m  )
                  enddo

              else

                  faca21 = panell(m+1,m+1) * panell(j,m+1)
                  facu21 = panell(m+1,m+1) * panelu(j,m+1)
c.debug
c     write(6,'("xdsld7-130-faca21,facu21          =", 1p2d15.5)')
c    1                       faca21,facu21 
c.debug

                  panell(j,j) = panell(j,j) - facu21 * panell(j,m+1)

cdir$ ivdep
                  do i = j+1, pnlrow
                      panell(i,j) = panell(i,j) - facu11 * panell(i,m  )
     1                                          - facu21 * panell(i,m+1)
                      panelu(i,j) = panelu(i,j) - faca11 * panelu(i,m  )
     1                                          - faca21 * panelu(i,m+1)
                  enddo

              end if
 
          enddo
 
  140 continue
 
c  =====================================================================
 
      return
      end
