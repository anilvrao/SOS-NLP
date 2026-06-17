      subroutine xdsleh( ncolf , pnlcol, pnlrow, updsiz, panel ,
     1                   pvtblk,
     1                   Sorder, locbfr, ocbufr, iops5u, wafil5,
     2                   watrn5, temp1 , temp2 , fctops, ierr )
 
c
c  purpose -- to apply the factored columns in panel to the remainder
c             of the front - this version has the columns stored
c             as columns in the panel
c
c  created            -- 07-feb-97, rgg
c  last modifications -- 09-mar-98, rgg, increased update panel size
c                                        from 4 to updsiz
c                     -- 16-apr-03, jgl  restructured to remove buffer
c                                        overrun problem when pivoting
c
c  input variables --
c
c      ncolf  -- number of columns that were factored
c      pnlcol -- number of columns in panel
c      pnlrow -- number of rows    in panel
c      updsiz -- factorization update panel size
c      panel  -- rectangular array holding factored columns
c      pvtblk -- indicator of 1x1 and 2x2 pivots
c      Sorder -- order of the (2,2) portion of the front (the number
c                of colums to receive the outer product modification
c      locbfr -- length of ocbufr
c      iops5u -- i/o position for the start of the section of the
c                front to be updated
c      wafil5 -- i/o file holding the front
c      fctops -- factor operation count
c
c  working storage --
c      
c      ocbufr -- i/o buffer
c      temp1  -- temporary array of size ncolf  by updsiz
c      temp2  -- temporary array of size Sorder by updsiz
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
 
      integer           ncolf , pnlcol, pnlrow, Sorder, locbfr,
     1                  iops5u, wafil5, ierr  , updsiz
 
      integer           pvtblk(*)
 
      double precision  watrn5, fctops
 
      double precision  panel(pnlrow,pnlcol), ocbufr(*),
     1                  temp1(ncolf,updsiz), temp2(Sorder,updsiz)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer           i

      logical           q1by1, lower

      double precision  fcol , ffront, ops
 
c  =====================================================================
 
      fcol   = ncolf
      ffront = Sorder
 
      ops    = (2*fcol) * ( ffront*(ffront+1) / 2. )
     1          + ffront * fcol

      fctops = fctops + ops

c     ----------------------------------------------------------------
c     ... check to see if the eliminated block required any 2x2 pivots
c     ----------------------------------------------------------------

      q1by1 = .true.
      lower = .true.

      do i = 1, ncolf
          if ( pvtblk(i) .eq. 2 ) then
              q1by1 = .false.
              exit
          end if
      enddo
 
      call xdsleq  ( ncolf , Sorder, updsiz, pnlrow, pnlcol,
     1               panel , panel , pvtblk, lower , q1by1 ,
     2               locbfr, ocbufr, iops5u, wafil5, watrn5,
     3               temp1 , temp2 , ierr )

      return

      end
