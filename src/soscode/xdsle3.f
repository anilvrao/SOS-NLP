      subroutine xdsle3 ( cmajor, ncolf , pnlcol, pnlrow, panel , 
     1                    fnzlf , qincor, fdiag , lnzloc,
     2                    iopos1, wafil1, ocbufr, lbuffr, ierr    )
 
c
c  purpose -- copy the factored columns from the current panel into
c             permanent storage for factorization.
c
c  created            -- 28-feb-97, rgg
c  last modifications -- 01-sep-98, rgg, 32 bit integer mods
c
c  input variables --
c
c      cmajor -- column major flag
c                .eq. 0 row    major form for panel is used
c                .ne. 0 column major form for panel is used
c      ncolf  -- number of columns actually factored for this panel.
c      pnlcol -- number of columns in this panel 
c      pnlrow -- number of rows    in this panel 
c      panel  -- panel matrix
c      qincor -- logical flag on whether l is to be stored in core 
c                or on disk
c      wafil1 -- i/o file for storing l
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
c      lnzloc -- local section of in-core storage for the factorization
c      ierr   -- error return,
c                if ierr = 0, success,
c                        =-2, i/o error on wafil1
c
c     ------------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer            cmajor, ncolf , pnlcol, pnlrow, 
     1                   wafil1, lbuffr, ierr

      integer            iopos1(2)
 
      logical            qincor

      double precision   fnzlf
 
      double precision   panel (*),      fdiag (*),
     1                   lnzloc(*),      ocbufr(*)

c     -------------------
c     ... local variables
c     -------------------

      integer            count , i     , k     , kd    , kl    ,
     1                   l     , length, rctbgn, rctstr, stride

c     ------------------------------------------------------------------
       
      if ( cmajor .eq. 0 ) then
          stride = pnlcol
          length = pnlcol
          rctbgn = pnlcol*ncolf + 1
          rctstr = 1
      else
          stride = 1
          length = pnlrow
          rctbgn = ncolf + 1
          rctstr = pnlrow
      end if

c     ----------------------------------
c     ... branch on in memory or on disk
c     ----------------------------------

      if ( qincor ) then

c         --------------------------------------------------------------
c         ... store lnz in memory.  first copy triangle entries of panel
c             into lnz. also copy diagonal into fdiag.
c         --------------------------------------------------------------

          kd = 1
          kl = 1

          do i = 1, ncolf
              fdiag(i) = panel(kd)
c.debug
c     write(6,'("dumping diag in-core - i, kd, panel(kd) = ", 
c    1       2i8, 1pd15.5)')             i, kd, panel(kd)
c.debug
              l        = ncolf - i
              do j=1,l
                lnzloc(kl+j-1) = panel(kd+j*stride)
              enddo
              kd       = kd + length + 1
              kl       = kl + l
          enddo

c         ---------------------------------------------------------
c         ... store the rest of the rectangle of the panel into lnz
c         ---------------------------------------------------------

          k  = rctbgn
          l  = pnlrow - ncolf

          if ( l .gt. 0 ) then

              do i = 1, ncolf
                  do j=0,l-1
                    lnzloc(kl+j) = panel(k+j*stride)
                  enddo
                  k  = k  + rctstr 
                  kl = kl + l
              enddo

          end if

          fnzlf = fnzlf + kl - 1

      else

c         ------------------------------------------------------------
c         ... store lnz on disk.  first copy triangle entries of panel
c             into ocbufr. when ocbufr fills up, empty it to disk and 
c             continue.  also copy diagonal into fdiag.
c         ------------------------------------------------------------

          kd    = 1
          count = 0

          do i = 1, ncolf

              fdiag(i) = panel(kd)
c.debug
c     write(6,'("dumping diag out-of  - i, kd, panel(kd) = ", 
c    1       2i8, 1pd15.5)')             i, kd, panel(kd)
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
                  fnzlf  = fnzlf  + count
                  count  = 0
              end if 
c.debug
c     write(6,'("dumping lnz - 300")')
c     write(6,'("l, kd+stride, stride, count+1 = ", 4i8)')
c    1            l, kd+stride, stride, count+1
c.debug

              do j=1,l
                ocbufr(count+j) = panel(kd+j*stride)
              enddo

              kd    = kd    + length + 1
              count = count + l

          enddo

c         ---------------------------------------------------------
c         ... store the rest of the rectangle of the panel into lnz
c         ---------------------------------------------------------

          k  = rctbgn
          l  = pnlrow - ncolf

          if ( l .gt. 0 ) then 

              do i = 1, ncolf

                  if ( count + l .gt. lbuffr ) then
                      call xdslw6 ( wafil1, ocbufr, 
     1                              iopos1, count, ierr )
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
                      fnzlf  = fnzlf  + count
                      count  = 0
                  end if 
c.debug
c     write(6,'("dumping lnz - 400")')
c     write(6,'("l, k, stride, count+1 = ", 4i8)')
c    1            l, k, stride, count+1
c.debug

                  do j=0,l-1
                    ocbufr(count+1+j) = panel(k+j*stride)
                  enddo

                  k     = k + rctstr
                  count = count + l

              enddo

          end if
 
c         -----------------------
c         ... dump rest of buffer
c         -----------------------

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
              fnzlf  = fnzlf  + count
          end if 

      end if
 
c     ------------------------------------------------------------------
  
      return
      end
