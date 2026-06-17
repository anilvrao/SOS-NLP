      subroutine xdsla4( lson  , offset, lindxl, wafil2, iopos2, 
     1                   watrn2, ocbufr, locbfr, colbgn, colend,
     2                   loclfr, panel , fctops, ierr )
c
c  purpose -- assemble section of out-of-core update matrix into the
c             current section of the out-of-core frontal matrix 
c
c  created   -- 19-mar-97, rgg 
c  modified  -- 06-aug-98, rgg, added offset to allow for out-of-core
c                               assembly with postponed columns
c
c  variables
c
c      lson   -- length of update matrix without postponed columns.
c      offset -- offset for indices in lindxl.  = 0 in phase 1.  
c                = number of postponed columns in phase 3.
c      lindxl -- local indices for assembling update matrix into front
c      wafil2 -- i/o file holding update matrix
c                columns
c      locbfr -- length of ocbufr
c      colbgn -- first column of this section of the frontal matrix
c      colend -- last  column of this section of the frontal matrix
c      loclfr -- size of this section of the frontal matrix
c
c  working storge  --
c
c      locbfr -- length of ocbufr array
c      ocbufr -- i/o buffer for reading in parts of the update matrix
c
c  input/output variables
c
c      watrn2 -- amount of i/o transfer for wafil2
c      iopos2 -- i/o position of update matrix after skipping postponed
c      panel  -- array holding numerical values for this section of the
c                front
c      fctops -- factor operation count
c
c  output variables --
c
c      ierr   -- error code
c                =  0 normal return
c                = -1 i/o error on sqfile
c
c  subprograms called
c
c      none
c
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           lson  , offset, wafil2, iopos2, locbfr, colbgn,
     1                  colend, loclfr, ierr

      integer           lindxl(*)
 
      double precision  watrn2, fctops
 
      double precision  ocbufr(*),      panel(*)

c     -------------------
c     ... local variables
c     -------------------
 
      integer           bfrlen, bfrpos, fclson, i     , ic    ,
     1                  j     , jc    , koffst, kpar  , lclson, 
     2                  len   , ncol  , nrow  , parbeg, skip  

      integer           idummy(1)

      double precision  temp  , temp2
c
c  =====================================================================

c     --------------------------------------------------
c     ... scan lindxl for columns to be read in for this
c         section of the front.  reset i/o position.
c         note that offset is 0 for phase 1 but is the
c         number of postponed columns for phase 3.
c     --------------------------------------------------

      fclson = lson + 1
      lclson = 0
c.debug
c     call xislp3 ( 'lindxl at start of xdsla4', lson, lindxl, 6 )
c     write(6,'("colbgn, colend, offset      = ", 4i8)') 
c    1            colbgn, colend, offset
c.debug

      do i = 1, lson

          if ( lindxl(i)+offset .lt. colbgn ) cycle
          if ( lindxl(i)+offset .gt. colend ) cycle

          fclson = min ( fclson, i )
          lclson = max ( lclson, i )

      enddo
c.debug
c     write(6,'("fclson, lclson              = ", 4i8)') 
c    1            fclson, lclson
c.debug

      skip   = ( lson*(fclson-1) ) - ( (fclson-1)*(fclson-2) / 2 )
      iopos2 = iopos2 + skip

c.debug
c     write(6,'("skip, iopos2                = ", 4i8)') 
c    1            skip, iopos2                
c     call xislp3 ( 'lindxl', lson, lindxl, 6 )
c.debug

c     --------------------------------------
c     ... scatter add the entries into front
c     --------------------------------------

      bfrlen = 0
      bfrpos = 1
 
      do 300 jc = fclson, lclson
 
          if ( bfrpos .gt. bfrlen ) then
 
              nrow   = lson - jc + 1
              len    = locbfr
              temp   = max ( ( .5 + nrow ) ** 2 - 2. * len, 0. )
              temp2  = lclson - jc + 1
              ncol   = min ( .5 + nrow - sqrt ( temp ), temp2 )
              bfrlen = nrow * ncol - ncol * ( ncol - 1 ) / 2
c.debug
c     write(6,'("jc, nrow, ncol, bfrlen      = ", 4i8)') 
c    1            jc, nrow, ncol, bfrlen      
c.debug
 
              call xdslw1 ( wafil2, 2, idummy, idummy, ocbufr, 
     1                      iopos2, bfrlen, ierr)
              if ( ierr .ne. 0 ) go to 8000
c.debug
c     write(6,'("in xdsla4 - iopos2, bfrlen = ", 3i8)') 
c    1            iopos2, bfrlen
c     call xdslp5 ( 'ocbufr', bfrlen, ocbufr, 6 )
c.debug
 
              iopos2 = iopos2 + bfrlen
              bfrpos = 1

              fctops = fctops + bfrlen
 
          end if
 
          j      = lindxl(jc) + offset - colbgn + 1
          parbeg = loclfr*(j-1) - (j*(j-1))/2
          koffst = parbeg - colbgn + 1
c.debug
c     write(6,'("j, parbeg                   = ", 4i8)') 
c    1            j, parbeg                   
c.debug

cdir$ ivdep
          do ic = jc, lson
c.debug
c     write(6,'("ic, jc, lson, lindxl(ic)+offset, parbeg, colbgn = ",
c    1  6i8)') ic, jc, lson, lindxl(ic)+offset, parbeg, colbgn 
c.debug
              kpar         = lindxl(ic) + offset + koffst
c.debug
c     write(6,'("ic, kpar, bfrpos, panel(kpar), ocbufr(bfrpos) = ",
c    1  3i8,1p2d15.5)') ic, kpar, bfrpos, panel(kpar), ocbufr(bfrpos)
c.debug

              panel (kpar) = panel (kpar) + ocbufr(bfrpos)
              bfrpos       = bfrpos + 1
          enddo

c.debug
c         if ( bfrpos .gt. bfrlen+1 ) then
c             write(6,'("oops after 200 in xdsla4")')
c             stop
c         end if
c.debug

  300 continue
 
      return
 
c  =====================================================================
 
c     --------------------------------------
c     ... error trap for i/o error on wafil2
c     --------------------------------------
 
 8000 continue
      ierr = -1
      return
 
c  =====================================================================
 
      end
