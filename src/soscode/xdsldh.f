      subroutine xdsldh( ncolf , pnlcol, pnlrow, updsiz, panell,
     1                   panelu, pvtblk, loclfr, locbfr, ocbufr,
     2                   matsiz, iops5u, wafil5, watrn5, 
     3                   temp1 , temp2 , fctops, ierr )
 
c
c  purpose -- to apply the factored columns in panel to the remainder
c             of the front - this version has the columns stored
c             as columns in the panel
c
c             unsymmetric and out-of-core version
c
c  created            -- 01-jun-98, rgg, from xdsld9 and xdsleh
c  last modifications -- 16-apr-03, jgl  restructured to remove buffer
c                                        overrun problem when pivoting
c
c  input variables --
c
c      ncolf  -- number of columns that were factored
c      pnlcol -- number of columns in panel
c      pnlrow -- number of rows    in panel
c      updsiz -- factorization update panel size
c      panell -- rectangular array holding factored columns
c                lower tri.
c      panelu -- rectangular array holding factored columns
c                upper tri.
c      pvtblk -- indicator of 1x1 and 2x2 pivots
c      loclfr -- local size of the front
c      iops5u -- i/o position for the start of the section of the
c                front to be updated
c      matsiz -- size of the lower tri. frontal matrix which is
c                the offset on the i/o file of the lower tri. and
c                upper tri. fronts
c      wafil5 -- i/o file holding the front
c      fctops -- factor operation count
c
c  working storage --
c      
c      ocbufr -- i/o buffer
c      temp1  -- temporary array of size ncolf  by updsiz
c      temp2  -- temporary array of size loclfr by updsiz
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
     1                  iops5u, wafil5, ierr  , updsiz, matsiz
 
      integer           pvtblk(*)
 
      double precision  watrn5, fctops
 
      double precision  panell(pnlrow,pnlcol), 
     1                  panelu(pnlrow,pnlcol), 
     2                  ocbufr(*), 
     4                  temp1(ncolf,updsiz) , 
     5                  temp2(loclfr,updsiz)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer           iopos

      logical           q1by1, lower

      double precision  fcol , ffront, ops
 
c  =====================================================================
 
      fcol   = ncolf
      ffront = loclfr
 
      ops    = (2*fcol) * ( ffront*(ffront+1) / 2. )
     1       + ffront * fcol

      fctops = fctops + 2 * ops

      q1by1 = .true.

c     ------------------------------------------
c     ... apply rank  ncolf  modification to the
c         lower triangle of the reduced matrix
c     ------------------------------------------
 
      iopos = iops5u
      lower = .true.

      call xdsleq  ( ncolf , loclfr, updsiz, pnlrow, pnlcol,
     1               panell, panelu, pvtblk, lower , q1by1 , 
     2               locbfr, ocbufr, iopos , wafil5, watrn5, 
     3               temp1 , temp2 , ierr )

      if ( ierr .ne. 0 ) then
         return
      end if

c     ------------------------------------------
c     ... apply rank  ncolf  modification to the 
c         upper triangle of the reduced matrix
c     ------------------------------------------
 
      iopos = iops5u + matsiz
      lower = .false.
 
      call xdsleq  ( ncolf , loclfr, updsiz, pnlrow, pnlcol,
     1               panell, panelu, pvtblk, lower , q1by1 ,
     2               locbfr, ocbufr, iopos , wafil5, watrn5,
     3               temp1 , temp2 , ierr )

      return

      end
