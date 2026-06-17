      subroutine xdsla3 ( nnode , lfront, nassmb, nstack, actual,
     1                    pospon, lindxl, rstbeg, istack, sup   , 
     2                    locbfr, ocbufr, wafil2, watrn2,
     3                    stfrnt, rstack, fctops, ierr )
 
c
c  purpose -- add postponed columns into expanded frontal matrix.
c             in-core version.
c
c  created   -- 05-feb-97, rgg
c  revisions -- 30-jun-03, dkw -- repaired factorization op count
c                                 calculation for stack out of core
c
c  input variables
c  ---------------
c
c      nnode  -- number of original columns that can be eliminated for
c                this front
c      lfront -- size of current front      
c      nassmb -- number of children to assemble
c      nstack -- total number of children on stack
c      lindxl -- local indices for assembly
c      rstbeg -- array containing update matrix starting positions in 
c                rstack
c      istack -- integer scalar information for each update matrix in
c                the stack
c      locbfr -- length of ocbufr array
c      wafil2 -- i/o file where update matrices are stored
c      stfrnt -- location for the start of the front in rstack
c
c  working storage
c  ---------------
c
c      ocbufr -- buffer for reading in update matrices
c
c  input/output variables
c  ----------------------
c
c      actual -- index to sup vector 
c      sup    -- contains column numbers in each supernode
c      watrn2 -- amount of i/o transfer to and from wafil2
c      rstack -- real working storage for factorization which holds any
c                update matrices in memory + current front.
c      fctops -- factorization operation count
c
c  output variables
c  ----------------
c
c      ierr   -- error code
c                =  0 normal return
c                = -1 i/o error on wafil2
c      
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           nnode , lfront, nassmb, nstack, actual, 
     1                  pospon, locbfr, wafil2, stfrnt, ierr

      integer           lindxl(*),      rstbeg(*),      istack(4,*), 
     1                  sup   (*)      
 
      double precision  watrn2, fctops

      double precision  ocbufr(*),      rstack(*)
 
c     -------------------
c     ... local variables
c     -------------------

      integer           bfrlen, bfrpos, i     , ii    , 
     1                  index , iopos2, irow  , j     , jcol  , 
     2                  jj    , jstack, kpar  , kson  , len   ,
     3                  lndxpt, lson  , lsonp , ncol  , nrow  ,
     4                  sonppb, sonppc, stkpnt

      integer           idummy(1)

      double precision  temp  , temp2
 
c     --------------------
c     ... subprograms used
c     --------------------

      external          icopy
 
c  =====================================================================

      jcol = nnode

      do 100 jstack = 1, nassmb
 
c         ----------------------------------------
c         establish pointers into stack entry for
c         easier access
c         ----------------------------------------
 
          stkpnt = rstbeg (nstack)
 
          lson   = istack (1, nstack)
          lsonp  = istack (3, nstack)
          lndxpt = istack (2, nstack)
c.debug
c     write(6,'("in xdsla3")')
c     write(6,'("nstack, stkpnt, lson, lsonp, lndxpt = ", 5i8)')
c    1            nstack, stkpnt, lson, lsonp, lndxpt 
c     write(6,'("jcol  , nnode , lfront, stfrnt      = ", 5i8)')
c    1            jcol  , nnode , lfront, stfrnt      
c.debug

          if ( lsonp .eq. 0 ) go to 90

          fctops = fctops + lsonp*(lsonp+1)/2 + lsonp*lson

          if ( stkpnt .gt. 0 ) then
 
c             ----------------------------
c             ... update matrix is in-core
c             ----------------------------
 
              sonppc = stkpnt
              sonppb = sonppc + lsonp
 
c             ------------------
c             ... get sup values
c             ------------------

c.debug
c     write(6,'("before sup - actual = ", i8)')
c     call xislp3 ( 'pre - sup ', actual-1, sup, 6 )
c.debug
              call icopy (lsonp, rstack (sonppc), 1, sup (actual), 1 )
              actual = actual + lsonp
c.debug
c     call xislp3 ( 'post - sup ', actual-1, sup, 6 )
c.debug
  
c             --------------------------------------------
c             ... scatter postponed columns into the front
c             --------------------------------------------
 
              kson = sonppb
             
              do 30 j = 1, lsonp

                  jcol = jcol + 1

                  kpar = stfrnt 
     1                 + lfront*(jcol-1) - ( (jcol-1)*(jcol-2) / 2 ) 

                  do i = j, lsonp
c.debug
c     write(6,'("in do 10 - i, j, lsonp, kson, kpar, rstack(kson) = ",
c    1        5i8,1d15.5)')  i, j, lsonp, kson, kpar, rstack(kson) 
c.debug
                      rstack(kpar) = rstack(kson)
                      kson        = kson + 1
                      kpar        = kpar + 1
                  enddo

                  index = lndxpt

                  do i = 1, lson

                      irow = lindxl(index)
                      if ( irow .gt. nnode ) irow = irow + pospon

                      ii   = max ( irow, jcol )
                      jj   = min ( irow, jcol )

                      kpar = stfrnt 
     1                     + lfront*(jj-1) - ( (jj-1)*(jj-2) / 2 ) 
     2                     + ( ii - jj )
c.debug
c     write(6,'("in do 20 - ii,jj,lsonp, kson, kpar, rstack(kson) = ",
c    1        5i8,1d15.5)')  ii,jj,lsonp, kson, kpar, rstack(kson) 
c.debug
 
                      rstack(kpar) = rstack(kson)
 
                      kson  = kson  + 1
                      index = index + 1

                  enddo

   30         continue

          else

c             ----------------------------------
c             ... update matrix is out-of-core.
c                 compute number of columns that
c                 can be read in at this time
c             ----------------------------------
 
              iopos2 = istack(4, nstack) 
 
              nrow   = lsonp + lson
              len    = locbfr - lsonp
              temp   = max ( ( .5 + nrow ) ** 2 - 2. * len, 0. )
              temp2  = lsonp
              ncol   = min ( .5 + nrow - sqrt ( temp ), temp2 )
              bfrlen = lsonp + nrow*ncol - ncol*(ncol-1)/2

              call xdslw1 ( wafil2, 2, idummy, idummy, ocbufr, 
     1                      iopos2, bfrlen, ierr)
              if ( ierr .ne. 0 ) go to 8900

              watrn2 = watrn2 + bfrlen 
 
              iopos2 = iopos2 + bfrlen
              bfrpos = 1
 
c             ------------------
c             ... get sup values
c             ------------------
 
              call icopy (lsonp, ocbufr, 1, sup (actual), 1 )
              actual = actual + lsonp
              bfrpos = bfrpos + lsonp
 
c             ---------------------------------------------
c             ... scatter postponed columns in beginning of
c                 front
c             ---------------------------------------------
 
              do 60 j = 1, lsonp

                  if ( bfrpos .gt. bfrlen ) then

                      nrow   = ( lsonp - j + 1 ) + lson
                      len    = locbfr
                      temp   = max ( ( .5 + nrow ) ** 2 - 2. * len,
     1                               0. )
                      temp2  = lsonp-j+1
                      ncol   = min ( .5 + nrow - sqrt ( temp ), 
     1                               temp2 )
                      bfrlen = nrow*ncol - ncol*(ncol-1)/2
   
                      call xdslw1 ( wafil2, 2, idummy, idummy, 
     1                              ocbufr, iopos2, bfrlen, ierr)
                      if ( ierr .ne. 0 ) go to 8900

                      watrn2 = watrn2 + bfrlen 
   
                      iopos2 = iopos2 + bfrlen
                      bfrpos = 1
  
                  end if

                  jcol = jcol + 1

                  kpar = stfrnt 
     1                 + lfront*(jcol-1) - ( (jcol-1)*(jcol-2) / 2 ) 

                  do i = j, lsonp
                      rstack(kpar) = ocbufr(bfrpos)
                      bfrpos       = bfrpos + 1
                      kpar         = kpar   + 1
                  enddo

                  index = lndxpt

                  do i = 1, lson

                      irow = lindxl(index)
                      if ( irow .gt. nnode ) irow = irow + pospon

                      ii   = max ( irow, jcol )
                      jj   = min ( irow, jcol )

                      kpar = stfrnt 
     1                     + lfront*(jj-1) - ( (jj-1)*(jj-2) / 2 ) 
     2                     + ( ii - jj )
 
                      rstack(kpar) = ocbufr(bfrpos)
 
                      bfrpos = bfrpos + 1
                      index  = index  + 1

                  enddo

   60         continue

          end if
 
   90     continue
          nstack = nstack - 1
 
  100 continue

      return
 
c  =====================================================================
 
c     --------------------------------------
c     ... error trap for i/o error on wafil2
c     --------------------------------------
 
 8900 continue
      ierr = -1
      return
 
c  =====================================================================
 
      end
