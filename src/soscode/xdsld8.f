
      subroutine xdsld8 ( ncolf , pnlcol, pnlrow, panell, panelu,
     1                    fnzlf , fdiag , matsiz, iopos1,
     2                    wafil1, wafil4, ocbufr, lbuffr, ierr    )
 
c
c  purpose -- copy the factored columns from the current panel into
c             permanent storage for factorization.
c             unsymmetric version - front out-of-core version
c
c  created            -- 27-jul-98, rgg, derived from xdsld3
c  last modifications -- 31-jul-98, rgg, added wafil4 
c                        31-aug-98, rgg, modified nzlf and i/o for
c                                        32 bit integer limits
c                        01-oct-01, dkw  set ierr = -4 for wafil4
c
c  input variables --
c
c      ncolf  -- number of columns actually factored for this panel.
c      pnlcol -- number of columns in this panel 
c      pnlrow -- number of rows    in this panel 
c      panell -- panel matrix - lower tri.
c      panelu -- panel matrix - upper tri.
c      matsiz -- size of lower tri. front - offset between lower tri.
c                and upper tri. fronts
c      wafil1 -- i/o file for storing l
c      wafil4 -- i/o file for storing u
c      lbuffr -- length of ocbufr
c
c  working storage --
c
c      ocbufr -- i/o buffer 
c
c  input/output variable --
c
c      fnzlf  -- number of nonzeroes stored in lnz so far
c      iopos1 -- current position on wafil1 for storing l on disk
c
c  output variable --
c
c      fdiag  -- local section of diagonal of the factorization
c      ierr   -- error return,
c                if ierr = 0, success,
c                        =-2, i/o error on wafil1
c                        =-4, i/o error on wafil4
c
c     ------------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer            ncolf , pnlcol, pnlrow, matsiz, 
     1                   wafil1, wafil4, lbuffr, ierr

      integer            iopos1(2)

      double precision   fnzlf
 
      double precision   panell(*),      panelu(*),      fdiag (*),
     1                   ocbufr(*)

c     -------------------
c     ... local variables
c     -------------------

      integer            count , i     , k     , kd    , 
     1                   l     , length, rctbgn, rctstr

      integer            iopos4(2)
 
c     ------------------------------------------------------------------
       
      length = pnlrow
      rctbgn = ncolf + 1
      rctstr = pnlrow

      iopos4(1) = iopos1(1)
      iopos4(2) = iopos1(2)

c     -------------------------------------------------------------
c     ... store lnz on disk.  first copy triangle entries of panell
c         into ocbufr. when ocbufr fills up, empty it to disk and 
c         continue.  also copy diagonal into fdiag.
c     -------------------------------------------------------------

      kd    = 1
      count = 0

      do i = 1, ncolf

          fdiag(i) = panell(kd)
c.debug
c     write(6,'("dumping diag out-of  - i, kd, panell(kd) = ",
c    1       2i8, 1pd15.5)')             i, kd, panell(kd)
c.debug

          l        = ncolf - i

          if ( count + l .gt. lbuffr ) then
              call xdslw6 ( wafil1, ocbufr, iopos1, count, ierr )
c.debug
c     write(6,'("dumping lnz out-of-core - iopos1, count = ", 2i8)') 
c    1                                      iopos1, count
c     call xdslp5 ( 'ocbufr', count, ocbufr, 6 )
c.debug
              if ( ierr .ne. 0 ) then
                  ierr = -2
                  return
              end if 

              call xdslw9 ( iopos1, count )
              fnzlf = fnzlf + count
              count = 0
          end if 
c.debug
c     write(6,'("dumping lnz - 300")')
c     write(6,'("l, kd+1, count+1 = ", 4i8)')
c    1            l, kd+1, count+1
c.debug

          ocbufr(count+1:count+l) = panell(kd+1:kd+l)

          kd    = kd    + length + 1
          count = count + l

      enddo

c     ---------------------------------------------------------
c     ... store the rest of the rectangle of the panel into lnz
c     ---------------------------------------------------------

      k  = rctbgn
      l  = pnlrow - ncolf

      if ( l .gt. 0 ) then 

          do i = 1, ncolf

              if ( count + l .gt. lbuffr ) then
                  call xdslw6 ( wafil1, ocbufr, iopos1, count, ierr )
c.debug
c     write(6,'("dumping lnz out-of-core - iopos1, count = ", 2i8)') 
c    1                                      iopos1, count
c     call xdslp5 ( 'ocbufr', count, ocbufr, 6 )
c.debug
                  if ( ierr .ne. 0 ) then
                      ierr = -2
                      return
                  end if 

                  call xdslw9 ( iopos1, count )
                  fnzlf = fnzlf + count
                  count = 0
              end if 
c.debug
c     write(6,'("dumping lnz - 310")')
c     write(6,'("l, k, count+1 = ", 4i8)')
c    1            l, k, count+1
c.debug

              ocbufr(count+1:count+l) = panell(k:k+l-1)

              k     = k + rctstr
              count = count + l

          enddo

      end if
 
c     -----------------------
c     ... dump rest of buffer
c     -----------------------

      if ( count .gt. 0 ) then
          call xdslw6 ( wafil1, ocbufr, iopos1, count, ierr )
c.debug
c     write(6,'("dumping lnz out-of-core - iopos1, count = ", 2i8)') 
c    1                                      iopos1, count
c     call xdslp5 ( 'ocbufr', count, ocbufr, 6 )
c.debug
          if ( ierr .ne. 0 ) then
              ierr = -2
              return
          end if

          call xdslw9 ( iopos1, count )
          fnzlf = fnzlf + count
          count = 0
      end if 

c     =============================================================

c     -------------------------------------------------------------
c     ... store unz on disk.  first copy triangle entries of panelu
c         into ocbufr. when ocbufr fills up, empty it to disk and 
c         continue.  
c     -------------------------------------------------------------

      kd     = 1
      count  = 0

      do i = 1, ncolf

          l        = ncolf - i

          if ( count + l .gt. lbuffr ) then
              call xdslw6 ( wafil4, ocbufr, iopos4, count, ierr )
c.debug
c     write(6,'("dumping unz out-of-core - iopos4, count = ", 2i8)') 
c    1                                      iopos4, count
c     call xdslp5 ( 'ocbufr', count, ocbufr, 6 )
c.debug
              if ( ierr .ne. 0 ) then
                  ierr = -4
                  return
              end if 

              call xdslw9 ( iopos4, count )
              fnzlf = fnzlf + count
              count = 0
          end if 
c.debug
c     write(6,'("dumping unz - 330")')
c     write(6,'("l, kd+1, count+1 = ", 4i8)')
c    1            l, kd+1, count+1
c.debug

          ocbufr(count+1:count+l) = panelu(kd+1:kd+l)

          kd    = kd    + length + 1
          count = count + l

      enddo

c     ---------------------------------------------------------
c     ... store the rest of the rectangle of the panel into unz
c     ---------------------------------------------------------

      k  = rctbgn
      l  = pnlrow - ncolf

      if ( l .gt. 0 ) then 

          do i = 1, ncolf

              if ( count + l .gt. lbuffr ) then
                  call xdslw6 ( wafil4, ocbufr, iopos4, count, ierr )
c.debug
c     write(6,'("dumping unz out-of-core - iopos4, count = ", 2i8)') 
c    1                                      iopos4, count
c     call xdslp5 ( 'ocbufr', count, ocbufr, 6 )
c.debug
                  if ( ierr .ne. 0 ) then
                      ierr = -4
                      return
                  end if 

                  call xdslw9 ( iopos4, count )
                  fnzlf = fnzlf + count
                  count = 0
              end if 
c.debug
c     write(6,'("dumping unz - 340")')
c     write(6,'("l, k, count+1 = ", 4i8)')
c    1            l, k, count+1
c.debug

              ocbufr(count+1:count+l) = panelu(k:k+l-1)

              k     = k + rctstr
              count = count + l

          enddo

      end if
 
c     -----------------------
c     ... dump rest of buffer
c     -----------------------

      if ( count .gt. 0 ) then
          call xdslw6 ( wafil4, ocbufr, iopos4, count, ierr )
c.debug
c     write(6,'("dumping unz out-of-core - iopos4, count = ", 2i8)') 
c    1                                      iopos4, count
c     call xdslp5 ( 'ocbufr', count, ocbufr, 6 )
c.debug
          if ( ierr .ne. 0 ) then
              ierr = -4
              return
          end if

          call xdslw9 ( iopos4, count )
          fnzlf = fnzlf + count
          count = 0
      end if 
 
c     ------------------------------------------------------------------
  
      return
      end
