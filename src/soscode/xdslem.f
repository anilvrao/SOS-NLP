      subroutine xdslem ( cmajor, ncolf , colbgn, colend, loclfr, 
     1                    pnlcol, pnlrow, panel , locbfr, ocbufr,
     2                    iops5u, wafil5, watrn5, ierr )
 
c
c  purpose -- stuff postponed columns back into the front.
c             out-of-memory version
c
c  created            -- 05-aug-98, rgg, derived from xdslek
c  last modifications -- 
c
c  input variables --
c
c      cmajor -- column major flag
c                .eq. 0 row    major form for panel is used
c                .ne. 0 column major form for panel is used
c      ncolf  -- number of columns factored
c      colbgn -- first column of panel  
c      colend -- last  column of panel  
c      loclfr -- size of front
c      pnlcol -- number of columns in panel
c      pnlrow -- number of rows    in panel
c      panel  -- rectangular the symmetric frontal matrix in packed storage
c      locbfr -- length of ocbufr
c      iops5u -- i/o position for the start of the section of the
c                front to be updated
c      wafil5 -- i/o file holding the front
c
c  working storage --
c
c      ocbufr -- i/o buffer
c
c  output variable --
c
c      watrn5 -- i/o transfer count for i/o file wafil5
c      ierr   -- i/o error return
c
c     ------------------------------------------------------------------
 
      integer            cmajor, ncolf , colbgn, colend, loclfr, pnlcol, 
     1                   pnlrow, locbfr, iops5u, wafil5, ierr
 
      double precision   watrn5
 
      double precision   ocbufr (*),     panel (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer            bfruse, iopos , jcol  , 
     1                   kk    , length, stride

      integer            idummy(1)
 
c     --------------------
c     ... subprograms used
c     --------------------
 
c     ------------------------------------------------------------------

      kk     = 1
      length = loclfr
      iopos  = iops5u
      bfruse = 0

      if ( cmajor .eq. 0 ) then
          stride = pnlcol
      else
          stride = 1
      end if

      do 100 jcol = colbgn, colend

c         -------------------------------------------------
c         ... move column jcol from panel into ocbufr which
c             is collecting the current part of the front.
c             note that bfruse points to the last used entry
c             in ocbufr and bfruse+1 will be the location 
c             for the next diagonal entry from the front.
c             kk points to the corresponding entry in the
c             panel.
c         -------------------------------------------------

          if ( jcol .le. ncolf ) go to 50

          if ( bfruse + length .gt. locbfr ) then
                  
              call xdslw2 ( wafil5, 2, idummy, idummy, ocbufr,
     1                      iopos, bfruse, ierr )
c.debug
c     write(6,'("dumping front out-of-core - iopos, bfruse = ", 2i8)')
c    1                                        iopos, bfruse
c     call xdslp5 ( 'ocbufr', bfruse, ocbufr, 6 )
c.debug
              if ( ierr .ne. 0 ) then
                  ierr = -2
                  return
              end if

              iopos  = iopos  + bfruse
              watrn5 = watrn5 + bfruse

              bfruse = 0

          end if

c.debug
c     write(6,'("in xdslem - jcol, colbgn, colend       = ", 4i8)')
c    1                        jcol, colbgn, colend  
c     write(6,'("            length, bfruse, kk, stride = ", 4i8)')
c    1                        length, bfruse, kk, stride
c.debug

          do j=0,length-1
            ocbufr(bfruse+1+j) = panel(kk+j*stride)
          enddo

          bfruse = bfruse + length

c         --------------------------
c         ... adjust for next column
c         --------------------------

   50     continue
          length = length - 1

          if ( cmajor .eq. 0 ) then
              kk = kk + pnlcol + 1
          else
              kk = kk + pnlrow + 1
          end if
          
  100 continue

c     ----------------------------------------------------------
c     ... finish up by writing last section of front to i/o file
c     ----------------------------------------------------------

      if ( bfruse .gt. 0 ) then 
                  
          call xdslw2 ( wafil5, 2, idummy, idummy, ocbufr,
     1                  iopos, bfruse, ierr )
c.debug
c     write(6,'("dumping front out-of-core - iopos, bfruse = ", 2i8)')
c    1                                        iopos, bfruse
c     call xdslp5 ( 'ocbufr', bfruse, ocbufr, 6 )
c.debug
          if ( ierr .ne. 0 ) then
              ierr = -2
              return
          end if

          watrn5 = watrn5 + bfruse

      end if

c     ------------------------------------------------------------------

      return
      end
