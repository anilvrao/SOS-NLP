      subroutine xdslef( ncolf , pnlcol, pnlrow, panel , pvtblk,
     1                   loclfr, front , temp1 , fctops  )
 
c
c  purpose -- to apply the factored columns in panel to the remainder
c             of the front - this version has the columns stored
c             as rows in the panel
c
c  created            -- 07-feb-97, rgg
c  last modifications -- 
c
c  input variables --
c
c      ncolf  -- number of columns that were factored
c      pnlcol -- number of columns in panel
c      pnlrow -- number of rows    in panel
c      panel  -- rectangular array holding factored columns
c      pvtblk -- indicator of 1x1 and 2x2 pivots
c      loclfr -- local size of the front
c      fctops -- factor operation count
c
c  working storage --
c      
c      temp1  -- temporary array of size ncolf by 3
c
c  output variable --
c
c      front  -- the remainder of the front
c      fctops -- factor operation count
c
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           ncolf , pnlcol, pnlrow, loclfr
 
      integer           pvtblk(*)
 
      double precision  fctops
 
      double precision  panel(pnlcol,pnlrow), front(*), temp1(ncolf,3)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer           i     , icol  , ii1   , ii2   , ii3   , irow  , 
     1                  j     , jcol  , jj    , k     ,
     2                  l1    , l2    , l3    , length, 
     3                  lstcol, lstrow, ncol  , nrow

      logical           q1by1

      double precision  dp11  , dp21  , dp31  , dp12  , dp22  ,
     1                  dp32  , dp13  , dp23  , dp33  , 
     2                  fcol  , ffront, ops   ,
     3                  diag1 , diag2 , offdia
 
c     --------------------
c     ... subprograms used
c     --------------------
 
c  =====================================================================
 
      fcol   = ncolf
      ffront = loclfr
 
      ops    = (2*fcol) * ( ffront*(ffront+1) / 2. )
     1       + ffront * fcol

      fctops = fctops + ops

c     --------------------------------
c     ... set 1by1 flag as appropriate
c     --------------------------------

      q1by1  = .true.

      do i = 1, ncolf
          if ( pvtblk(i) .eq. 2 ) then
              q1by1 = .false.
              exit
          end if
      enddo
c.debug
c     write(6,'("q1by1 = ", l5)') q1by1
c.debug
 
c  ----------------------------------------------------------

c     -----------------------------------------
c     ... loop over the remainder of the front.
c     -----------------------------------------

      l3     = 1
      length = loclfr
      lstcol = 0

      do 250 jcol = 1, loclfr, 3

          ncol  = min ( 3, loclfr - jcol + 1 )

c         ----------------------------
c         ... form temp = d * trans(l)
c         ----------------------------

          if ( q1by1 ) then

              do icol = 1, ncol
                  jj = jcol + pnlcol + icol - 1
                  do j = 1, ncolf
                      temp1(j,icol) = panel(j,j) * panel(j,jj)
                  enddo
              enddo

          else

              j = 1

   70         continue
              if ( j .gt. ncolf ) go to 100

              if ( pvtblk(j) .eq. 1 ) then
     
                  do icol = 1, ncol
                      jj = jcol + pnlcol + icol - 1
                      temp1(j,icol) = panel(j,j) * panel(j,jj)
                  enddo

                  j = j + 1

              else

                  diag1  = panel(j  ,j)
                  offdia = panel(j  ,j+1)
                  diag2  = panel(j+1,j+1)
     
                  do icol = 1, ncol
                      jj = jcol + pnlcol + icol - 1
                      temp1(j  ,icol) = diag1  * panel(j  ,jj)
     1                                + offdia * panel(j+1,jj)
                      temp1(j+1,icol) = offdia * panel(j  ,jj)
     1                                + diag2  * panel(j+1,jj)
                  enddo

                  j = j + 2

              end if

              go to 70

  100         continue          

          end if
c.debug
c     write(6,'("ncol, ncolf = ", 2i8)') ncol, ncolf
c     call xdslp5('after 100 - temp1', ncol*ncolf, temp1, 6 )
c.debug

c         -------------------------------------------
c         ... set up to form front = front - l * temp
c         -------------------------------------------

          l1 = l3
          l2 = l1 + length 
          l3 = l2 + length - 1

          ii1 = jcol + pnlcol 
          ii2 = jcol + pnlcol + 1

          lstrow = jcol+1
          lstcol = jcol+2
c.debug
c     write(6,'("in 250 - l1, l2, l3, ii1, ii2 = ", 5i8)') 
c    1                     l1, l2, l3, ii1, ii2 
c     write(6,'("l1  , front(l1  ) = ", i8, 1pd15.5)')
c    1            l1  , front(l1  ) 
c     write(6,'("l1+1, front(l1+1) = ", i8, 1pd15.5)')
c    1            l1+1, front(l1+1) 
c     write(6,'("l2  , front(l2  ) = ", i8, 1pd15.5)')
c    1            l2  , front(l2  ) 
c.debug

          front(l1  ) = front(l1  )
     1                - dot_product(temp1(1:ncolf,1),panel(1:ncolf,ii1))

          if ( jcol   .eq. loclfr ) go to 300

          front(l1+1) = front(l1+1)
     1                - dot_product(temp1(1:ncolf,1),panel(1:ncolf,ii2))

          front(l2  ) = front(l2  )
     1                - dot_product(temp1(1:ncolf,2),panel(1:ncolf,ii2))

          if ( jcol+1 .eq. loclfr ) go to 300

c.debug
c     write(6,'("l1  , front(l1  ) = ", i8, 1pd15.5)')
c    1            l1  , front(l1  ) 
c     write(6,'("l1+1, front(l1+1) = ", i8, 1pd15.5)')
c    1            l1+1, front(l1+1) 
c     write(6,'("l2  , front(l2  ) = ", i8, 1pd15.5)')
c    1            l2  , front(l2  ) 
c.debug

          l1 = l1 + 2
          l2 = l2 + 1

c         ----------------------------------------------------
c         ... start of main computational loop.  this performs
c             a 3x3 array of inner products between the three
c             scaled rows held in temp1 with
c             the next set of three rows held in panel.
c         ----------------------------------------------------

          do 200 irow = jcol+2, loclfr, 3

               nrow = min ( 3, loclfr - irow + 1 )
 
               ii1 = irow + pnlcol
               ii2 = ii1 + 1
               ii3 = ii1 + 2
c.debug
c     write(6,'("in 200 - irow, ii1, ii2, ii3, l1, l2, l3 = ", 7i8)') 
c    1                     irow, ii1, ii2, ii3, l1, l2, l3 
c.debug

               if ( nrow .le. 2 ) then
 
                   front(l1) = front(l1)
     1                - dot_product(temp1(1:ncolf,1),panel(1:ncolf,ii1))
 
                   front(l2) = front(l2)
     1                - dot_product(temp1(1:ncolf,2),panel(1:ncolf,ii1))
 
                   front(l3) = front(l3)
     1                - dot_product(temp1(1:ncolf,3),panel(1:ncolf,ii1))

                   l1 = l1 + 1
                   l2 = l2 + 1
                   l3 = l3 + 1

                   if ( nrow .eq. 1 ) go to 200
 
                   front(l1) = front(l1)
     1                - dot_product(temp1(1:ncolf,1),panel(1:ncolf,ii2))
 
                   front(l2) = front(l2)
     1                - dot_product(temp1(1:ncolf,2),panel(1:ncolf,ii2))
 
                   front(l3) = front(l3)
     1                - dot_product(temp1(1:ncolf,3),panel(1:ncolf,ii2))

                   l1 = l1 + 1
                   l2 = l2 + 1
                   l3 = l3 + 1

                   go to 200

               end if

c              ----------------------------------------------------
c              ... the code for the inner loop could be written as:
c
c              front(l1  ) = front(l1  )
c    1                   - dot ( ncolf, temp1(1,1), 1, panel(1,ii1) 1 )
c
c              front(l1+1) = front(l1+1)
c    1                   - dot ( ncolf, temp1(1,1), 1, panel(1,ii2), 1 )
c
c              front(l1+2) = front(l1+2)
c    1                   - dot ( ncolf, temp1(1,1), 1, panel(1,ii3), 1 )
c
c              front(l2  ) = front(l2  )
c    1                   - dot ( ncolf, temp1(1,2), 1, panel(1,ii1), 1 )
c
c              front(l2+1) = front(l2+1)
c    1                   - dot ( ncolf, temp1(1,2), 1, panel(1,ii2), 1 )
c
c              front(l2+2) = front(l2+2)
c    1                   - dot ( ncolf, temp1(1,2), 1, panel(1,ii3), 1 )
c
c              front(l3  ) = front(l3  )
c    1                   - dot ( ncolf, temp1(1,3), 1, panel(1,ii1), 1 )
c
c              front(l3+1) = front(l3+1)
c    1                   - dot ( ncolf, temp1(1,3), 1, panel(1,ii2), 1 )
c
c              front(l3+2) = front(l3+2)
c    1                   - dot ( ncolf, temp1(1,3), 1, panel(1,ii3), 1 )
c
c              ----------------------------------------------------

               dp11 = 0.
               dp21 = 0.
               dp31 = 0.
               dp12 = 0.
               dp22 = 0.
               dp32 = 0.
               dp13 = 0.
               dp23 = 0.
               dp33 = 0.

               do k = 1, ncolf
                   dp11 = dp11 + temp1(k,1) * panel(k,ii1)
                   dp21 = dp21 + temp1(k,1) * panel(k,ii2)
                   dp31 = dp31 + temp1(k,1) * panel(k,ii3)
                   dp12 = dp12 + temp1(k,2) * panel(k,ii1)
                   dp22 = dp22 + temp1(k,2) * panel(k,ii2)
                   dp32 = dp32 + temp1(k,2) * panel(k,ii3)
                   dp13 = dp13 + temp1(k,3) * panel(k,ii1)
                   dp23 = dp23 + temp1(k,3) * panel(k,ii2)
                   dp33 = dp33 + temp1(k,3) * panel(k,ii3)
               enddo

               front(l1  ) = front(l1  ) - dp11
               front(l1+1) = front(l1+1) - dp21
               front(l1+2) = front(l1+2) - dp31
               front(l2  ) = front(l2  ) - dp12
               front(l2+1) = front(l2+1) - dp22
               front(l2+2) = front(l2+2) - dp32
               front(l3  ) = front(l3  ) - dp13
               front(l3+1) = front(l3+1) - dp23
               front(l3+2) = front(l3+2) - dp33

c              --------------------------------------------------
c              ... adjust the pointers into the front for storing 
c                  the results of inner products
c              --------------------------------------------------

               l1 = l1 + 3
               l2 = l2 + 3
               l3 = l3 + 3

               lstrow = lstrow + 3

  200     continue

c         ------------------------------------
c         ... advance to next set of 3 columns
c         ------------------------------------

          length = length - 3

  250 continue

c.debug
c     write(6,'("after 250 - jcol, lstcol, loclfr = ", 3i8)')
c    1                        jcol, lstcol, loclfr 
c.debug
 
c  =====================================================================
 
  300 continue
 
      return
      end
