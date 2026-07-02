      subroutine xdslej( ncolf , pnlcol, pnlrow, panel , pvtblk,
     1                   loclfr, locbfr, ocbufr, iops5u, wafil5, 
     2                   watrn5, temp1 , fctops, ierr  )
 
c
c  purpose -- to apply the factored columns in panel to the remainder
c             of the front - this version has the columns stored
c             as rows in the panel
c
c  created            -- 07-feb-97, rgg
c  last modifications -- 01-sep-98, rgg, 32 bit integer mods
c
c  input variables --
c
c      ncolf  -- number of columns that were factored
c      pnlcol -- number of columns in panel
c      pnlrow -- number of rows    in panel
c      panel  -- rectangular array holding factored columns
c      pvtblk -- indicator of 1x1 and 2x2 pivots
c      loclfr -- local size of the front
c      locbfr -- length of ocbufr
c      iops5u -- i/o position for the start of the section of the
c                front to be updated
c      wafil5 -- i/o file holding the front
c      fctops -- factor operation count
c
c  working storage --
c      
c      ocbufr -- i/o buffer
c      temp1  -- temporary array of size ncolf by 3
c
c  output variable --
c
c      watrn5 -- i/o transfer count for i/o file wafil5
c      fctops -- factor operation count
c      ierr   -- i/o error return
c
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           ncolf , pnlcol, pnlrow, loclfr, locbfr,
     1                  iops5u, wafil5, ierr
 
      integer           pvtblk(*)
 
      double precision  watrn5, fctops
 
      double precision  panel(pnlcol,pnlrow),   ocbufr(*), 
     1                  temp1(ncolf,3)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer           bfrlen, i     , icol  , ii1   , ii2   , ii3   , 
     1                  iopos , irow  , j     , jcol  , jj    , k     ,
     2                  l1    , l2    , l3    , lend  , length, 
     3                  lstcol, lstrow, ncol  , nrow  , updsiz

      integer           idummy(1)

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
      updsiz = loclfr * ( loclfr + 1 ) / 2
 
      ops    = (2*fcol) * ( ffront*(ffront+1) / 2. )
     1       + ffront * fcol

      fctops = fctops + ops
 
      bfrlen = 0
      iopos  = iops5u

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
      lend   = 0

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

c         ----------------------------
c         ... perform i/o as necessary
c         ----------------------------

          if ( l3 + length - 3 .gt. bfrlen .or. l1 .eq. 1 ) then

c             ----------------------------------------------
c             ... time to read in next section of the front.
c                 first dump current section to i/o file.
c             ----------------------------------------------

              if ( l1 .gt. 1 ) then 
                  
                  call xdslw2 ( wafil5, 2, idummy, idummy, ocbufr,
     1                          iopos, l1-1, ierr )
c.debug
c     write(6,'("dumping lnz out-of-core - iopos, l1-1 = ", 2i8)')
c    1                                      iopos, l1-1
c     call xdslp5 ( 'ocbufr', l1-1, ocbufr, 6 )
c.debug
                  if ( ierr .ne. 0 ) then
                      ierr = -2
                      return
                  end if

                  iopos  = iopos  + l1 - 1
                  watrn5 = watrn5 + l1 - 1
                  updsiz = updsiz - ( l1 - 1 )

              end if

c             --------------------------------------
c             ... read in next section to be updated
c             --------------------------------------

              bfrlen = min ( locbfr, updsiz )
                  
              call xdslw1 ( wafil5, 2, idummy, idummy, ocbufr,
     1                      iopos, bfrlen, ierr )
c.debug
c     write(6,'("reading in next section to update - iopos, bfrlen = ", 
c    1    2i8)') iopos, bfrlen
c     call xdslp5 ( 'ocbufr', bfrlen, ocbufr, 6 )
c.debug
              if ( ierr .ne. 0 ) then
                  ierr = -3
                  return
              end if

              watrn5 = watrn5 + bfrlen

              l1 = 1
              l2 = l1 + length 
              l3 = l2 + length - 1

c.debug
              if ( l3 + length - 3 .gt. bfrlen ) then
                  write(6,'("storage oops in xdslej")')
                  write(6,'("l1, l2, l3, length, bfrlen = ", 5i8)')
     1                        l1, l2, l3, length, bfrlen 
                  RETURN
              end if
c.debug

          end if

c         -------------------------------------------------
c         ... start update of current set of columns of the
c             front.
c         -------------------------------------------------

          ii1 = jcol + pnlcol 
          ii2 = jcol + pnlcol + 1

          lstrow = jcol+1
          lstcol = jcol+2
c.debug
c     write(6,'("in 250 - l1, l2, l3, ii1, ii2 = ", 5i8)') 
c    1                     l1, l2, l3, ii1, ii2 
c     write(6,'("l1  , ocbufr(l1  ) = ", i8, 1pd15.5)')
c    1            l1  , ocbufr(l1  ) 
c     write(6,'("l1+1, ocbufr(l1+1) = ", i8, 1pd15.5)')
c    1            l1+1, ocbufr(l1+1) 
c     write(6,'("l2  , ocbufr(l2  ) = ", i8, 1pd15.5)')
c    1            l2  , ocbufr(l2  ) 
c.debug

          ocbufr(l1  ) = ocbufr(l1  )
     1                - dot_product(temp1(1:ncolf,1),panel(1:ncolf,ii1))

          if ( jcol   .eq. loclfr ) then
              lend = l1
              go to 300
          end if

          ocbufr(l1+1) = ocbufr(l1+1)
     1                - dot_product(temp1(1:ncolf,1),panel(1:ncolf,ii2))

          ocbufr(l2  ) = ocbufr(l2  )
     1                - dot_product(temp1(1:ncolf,2),panel(1:ncolf,ii2))

          if ( jcol+1 .eq. loclfr ) then
              lend = l2
              go to 300
          end if

c.debug
c     write(6,'("l1  , ocbufr(l1  ) = ", i8, 1pd15.5)')
c    1            l1  , ocbufr(l1  ) 
c     write(6,'("l1+1, ocbufr(l1+1) = ", i8, 1pd15.5)')
c    1            l1+1, ocbufr(l1+1) 
c     write(6,'("l2  , ocbufr(l2  ) = ", i8, 1pd15.5)')
c    1            l2  , ocbufr(l2  ) 
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
 
                   ocbufr(l1) = ocbufr(l1)
     1                - dot_product(temp1(1:ncolf,1),panel(1:ncolf,ii1))
 
                   ocbufr(l2) = ocbufr(l2)
     1                - dot_product(temp1(1:ncolf,2),panel(1:ncolf,ii1))
 
                   ocbufr(l3) = ocbufr(l3)
     1                - dot_product(temp1(1:ncolf,3),panel(1:ncolf,ii1))

                   l1 = l1 + 1
                   l2 = l2 + 1
                   l3 = l3 + 1

                   if ( nrow .eq. 1 ) go to 200
 
                   ocbufr(l1) = ocbufr(l1)
     1                - dot_product(temp1(1:ncolf,1),panel(1:ncolf,ii2))
 
                   ocbufr(l2) = ocbufr(l2)
     1                - dot_product(temp1(1:ncolf,2),panel(1:ncolf,ii2))
 
                   ocbufr(l3) = ocbufr(l3)
     1                - dot_product(temp1(1:ncolf,3),panel(1:ncolf,ii2))

                   l1 = l1 + 1
                   l2 = l2 + 1
                   l3 = l3 + 1

                   go to 200

               end if

c              ----------------------------------------------------
c              ... the code for the inner loop could be written as:
c
c              ocbufr(l1  ) = ocbufr(l1  )
c    1                    - dot(ncolf, temp1(1,1), 1, panel(1,ii1) 1 )
c
c              ocbufr(l1+1) = ocbufr(l1+1)
c    1                    - dot(ncolf, temp1(1,1), 1, panel(1,ii2), 1 )
c
c              ocbufr(l1+2) = ocbufr(l1+2)
c    1                    - dot(ncolf, temp1(1,1), 1, panel(1,ii3), 1 )
c
c              ocbufr(l2  ) = ocbufr(l2  )
c    1                    - dot(ncolf, temp1(1,2), 1, panel(1,ii1), 1 )
c
c              ocbufr(l2+1) = ocbufr(l2+1)
c    1                    - dot(ncolf, temp1(1,2), 1, panel(1,ii2), 1 )
c
c              ocbufr(l2+2) = ocbufr(l2+2)
c    1                    - dot(ncolf, temp1(1,2), 1, panel(1,ii3), 1 )
c
c              ocbufr(l3  ) = ocbufr(l3  )
c    1                    - dot(ncolf, temp1(1,3), 1, panel(1,ii1), 1 )
c
c              ocbufr(l3+1) = ocbufr(l3+1)
c    1                    - dot(ncolf, temp1(1,3), 1, panel(1,ii2), 1 )
c
c              ocbufr(l3+2) = ocbufr(l3+2)
c    1                    - dot(ncolf, temp1(1,3), 1, panel(1,ii3), 1 )
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

               ocbufr(l1  ) = ocbufr(l1  ) - dp11
               ocbufr(l1+1) = ocbufr(l1+1) - dp21
               ocbufr(l1+2) = ocbufr(l1+2) - dp31
               ocbufr(l2  ) = ocbufr(l2  ) - dp12
               ocbufr(l2+1) = ocbufr(l2+1) - dp22
               ocbufr(l2+2) = ocbufr(l2+2) - dp32
               ocbufr(l3  ) = ocbufr(l3  ) - dp13
               ocbufr(l3+1) = ocbufr(l3+1) - dp23
               ocbufr(l3+2) = ocbufr(l3+2) - dp33

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

          lend   = l3 - 1

  250 continue

c.debug
c     write(6,'("after 250 - jcol, lstcol, loclfr = ", 3i8)')
c    1                        jcol, lstcol, loclfr 
c.debug
 
c  =====================================================================
 
c     ----------------------------------------------------------
c     ... finish up by writing last section of front to i/o file
c     ----------------------------------------------------------

  300 continue
      if ( lend .gt. 0 ) then

          call xdslw2 ( wafil5, 2, idummy, idummy, ocbufr,
     1                  iopos, lend, ierr )
c.debug
c     write(6,'("dumping lnz out-of-core - iopos, lend = ", 2i8)')
c    1                                      iopos, lend
c     call xdslp5 ( 'ocbufr', lend, ocbufr, 6 )
c.debug
          if ( ierr .ne. 0 ) then
              ierr = -2
              return
          end if

          watrn5 = watrn5 + lend

      end if

c  =====================================================================
 
      return
      end
