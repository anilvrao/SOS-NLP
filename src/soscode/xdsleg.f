      subroutine xdsleg ( cmajor, colbgn, colend, loclfr, 
     1                    iops5f, wafil5, watrn5, locbfr, ocbufr,
     2                    pnlcol, pnlrow, panel , error )

 
c
c  purpose -- read in the current panel from the i/o file holding
c             the current front.
c
c  created            -- 24-mar-97, rgg
c  last modifications -- 31-jan-01, jgl -- error handling changed
c
c  input variables --
c
c      cmajor -- column major flag
c                .eq. 0 row    major form for panel is used
c                .ne. 0 column major form for panel is used
c      colbgn -- first column to extract
c      colend -- last  column to extract
c      loclfr -- size of front
c      iops5f -- i/o position for the start of the current panel
c                of the front
c      wafil5 -- i/o file holding the front
c      watrn5 -- amount of i/o transfer on wafil5
c      locbfr -- length of ocbufr
c      pnlcol -- number of columns in panel
c      pnlrow -- number of rows    in panel
c 
c  working storage --
c
c      ocbufr -- i/o buffer
c
c  output variables --
c      panel  -- rectangular the symmetric frontal matrix in packed storage
c      error  -- return code
c                =   0, success
c                =  -1, error reading file wafil5 (single front)
c
c     ------------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer            cmajor, colbgn, colend, loclfr, iops5f,
     1                   wafil5, locbfr, pnlcol, pnlrow, error
 
      double precision   watrn5
 
      double precision   ocbufr (*),     panel (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      parameter  (zero = 0.d0)

      integer            iopos , jcol  , k     , kk    , interr,
     1                   l     , len   , length, ncol  , valend

      integer            idummy(1)
 
c     ------------------------------------------------------------------

      error  = 0
      k      = 1
      kk     = 1
      length = loclfr
      iopos  = iops5f

      if ( cmajor .eq. 0 ) then

c         --------------------
c         ... row major format
c         --------------------

          valend = 0
          l      = pnlrow
          k      = 1
          kk     = 1
          iopos  = iops5f

          do jcol = colbgn, colend

c             -------------------------------------------------
c             ... check if next buffer of data needs to be read
c             -------------------------------------------------

              if ( k + l - 1 .gt. valend ) then

                  ncol   = colend - jcol + 1
                  len    = min ( locbfr, l*ncol - ncol*(ncol-1)/2 )
                  iopos  = iopos + k - 1
c.debug
c     write(6,'("in xdsleg before xdslw1-jcol, len, iopos = ", 3i8)')
c    1                                    jcol, len, iopos 
c.debug

                  call xdslw1 ( wafil5, 2, idummy, idummy, ocbufr, 
     1                          iopos , len, interr )
c.debug
c     write(6,'("in xdsleg after xdslw1 - interr = ", i8)') interr
c.debug

                  if ( interr .ne. 0 ) then
                      error = -1
                      return
                  end if

                  k      = 1
                  valend = len
                  watrn5 = watrn5 + len

              end if
c.debug
c     write(6,'("in xdsleg - jcol, colbgn, colend  = ", 4i8)')
c    1                        jcol, colbgn, colend  
c     write(6,'("            k, kk, l, valend      = ", 4i8)')
c    1                        k, kk, l, valend
c.debug

              do jj=0,l-1
                panel(kk+jj*pnlcol) = ocbufr(k+jj)
              enddo

c.debug
c             j = jcol - colbgn
c             do jj=1,j
c               panel(kk-jj*pnlcol) = zero
c             enddo
c.debug

              kk = kk + pnlcol + 1
              k  = k + l
              l  = l-1

          enddo

      else

c         -----------------------
c         ... column major format
c         -----------------------

          len = pnlcol * pnlrow - ( pnlcol * ( pnlcol-1 ) / 2 )
c.debug
c     write(6,'("in xdsleg before xdslw1 - len  = ", i8)') len 
c.debug

          call xdslw1 ( wafil5, 2, idummy, idummy, panel, iopos,
     1                  len, interr )
c.debug
c     write(6,'("in xdsleg after xdslw1 - interr = ", i8)') interr
c     call xdslp5('panel after read', len, panel, 6 )
c.debug

          if ( interr .ne. 0 ) then
              error = -1
              return
          end if

          watrn5 = watrn5 + len

c         --------------------------------------------
c         ... rearrange triangular packed storage into 
c             rectangular form
c         --------------------------------------------

          k  = len + 1
          kk = pnlcol * pnlrow + 1

          len = pnlrow - pnlcol + 1

          do jcol = colend, colbgn, -1

c.debug
c     write(6,'("in xdsleg - jcol, colbgn, colend  = ", 4i8)')
c    1                        jcol, colbgn, colend  
c     write(6,'("            k, kk, len = ", 4i8)')
c    1                        k, kk, len
c.debug
              k   = k  - len
              kk  = kk - len

              call xdslmv ( len, panel, k, kk ) 


              l   = jcol - colbgn
              kk  = kk   - l
c.debug
              panel(kk:kk+l-1) = zero
c.debug

              len = len + 1

          enddo

c.debug
c         if ( k .ne. 1 .or. kk .ne. 1 ) then
c             write(6,'("oops after 200 in xdsleg - k, kk = ",
c    1                  2i8)') k, kk
c             stop
c         end if
c.debug

      end if

c     ------------------------------------------------------------------
 
      return
      end 
