      subroutine xdslai( mincor, diag  , xladj , ladjl , lnza  ,
     1                   bassmb, sigma , bmxtyp, bdiag , bxladj, bladjl,
     2                   lnzb  , sup   , nassmb, xlindx, lindxl,
     3                   bgxsup, dspmnt, lfront, loclfr, matsiz,
     4                   mstack, nodbgn, nodnxt, nstack, nwfree, pospon,
     5                   stknod, rstbeg, istack, rstack, stfrnt, sqfil1,
     6                   sqtrn1, sqfil3, sqtrn3, wafil2, watrn2, locbfr,
     7                   ocbufr, mxvlib, mxvlrb, fctops, ierr   )
 
c
c  purpose -- assemble the frontal matrix for the current supernode.
c
c             if sigma .ne. 0. then assemble a - sigma * b.
c
c  created   -- 20-jul-89, rgg
c  revisions -- 12-feb-91, rgg -- added watrn2
c               04-nov-92, rgg -- added use of ocbufr to reduce number
c                                 of i/o access to wafil2.
c               04-dec-92, rgg -- added out-of-memory assembly of a
c                                 and b.
c               04-feb-97, rgg -- modified to move postponed columns 
c                                 to the end of columns to be eliminated.
c
c  variables
c
c      mincor -- minimum core processing switch.  if mincor .ne. 0
c                then either a or b are on i/o file sqfil1 or sqfil3.
c      diag   -- array to hold first the diagonal entries of the
c                original matrix and later those of the factor
c      xladj  -- pointers into indices and entries of the matrix
c      ladjl  -- local indices of the column indices of original matrix
c      lnza   -- column entries of the original matrix
c      sigma  -- shift value
c      bmxtyp -- matrix type of b
c                0 - null matrix
c                1 - general sparse
c                2 - general sparse but same structure as a
c                3 - diagonal
c                4 - identity
c      bdiag  -- diagonal of b
c      bxladj -- xladj for b
c      bladjl -- ladjl for b
c      lnzb   -- column entries of the b matrix
c      nassmb -- number of assemblies for this front
c      xlindx -- pointers to the column indices of the supernodes
c                in the factor
c      lindxl -- local indices of the supernodes with respect to
c                their supernode parent
c      stknod -- integer array to hold the supernode number of
c                each stacked matrix
c      rstbeg -- integer array to hold the starting location of each
c                stacked matrix
c      istack -- 2d integer array to hold the scalar information
c                associated with pivoting for a stacked matrix.
c      rstack -- stack space
c      stfrnt -- pointer into rstack for the start of the current front
c      sqfil1 -- i/o file holding a 
c      sqfil3 -- i/o file holding b
c      wafil2 -- i/o file holding spilled parts of the stack
c
c  working storge  --
c
c      locbfr -- length of ocbufr array
c      ocbufr -- i/o buffer for bringing in update matrices from the
c                spilled stack.
c      mxvlib -- i/o buffer for bringing in scatter indicies from the
c                spilled matrix file, sqfil1 or sqfil3.
c      mxvlrb -- i/o buffer for bringing in matrix values from the
c                spilled matrix file, sqfil1 or sqfil3.
c
c  input/output variables
c
c      nstack -- number of matrices stacked.  decremented by
c                the number of assemblies performed.
c      sqtrn1 -- amount of i/o transfered to and from sqfil1.
c      sqtrn3 -- amount of i/o transfered to and from sqfil3.
c      watrn2 -- amount of i/o transfered to and from wafil2.
c      fctops -- number of floating point operations to compute
c                factorization
c
c  output variables --
c
c      ierr   -- error code
c                =  0 normal return
c                = -1 i/o error on wafil2
c                = -2 i/o error on sqfil1 or sqfil3
c
c  subprograms called
c
c      xdslf7 -- add an update matrix to the frontal matrix
c      xdslf5 -- slide an unpate matrix into the frontal matrix
c      xdslf9 -- scatter/add a vector into another
c
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           bmxtyp, lfront, loclfr, matsiz, mstack,
     1                  nassmb, nodbgn, nodnxt, nstack, nwfree,
     2                  bgxsup, pospon, dspmnt, stfrnt, wafil2,
     3                  locbfr, mincor, sqfil1, sqfil3, ierr  
 
      integer           xladj(*),               ladjl(*),
     1                  bxladj(*),              bladjl(*),
     3                  xlindx(*),              lindxl(*),
     4                  stknod(*),              rstbeg(*),
     5                  istack(4,*),            sup (*),
     6                  mxvlib(*)
 
      logical           bassmb
 
      double precision  sigma , fctops, sqtrn1, sqtrn3, watrn2
 
      double precision  diag(*),                lnza(*),
     1                  bdiag(*),               lnzb(*),
     2                  rstack(*),              ocbufr(*),
     3                  mxvlrb(*)
 
c     -------------------
c     ... local variables
c     -------------------

      integer           actual, bfrlen, bfrpos, endson, iopos2, jstkst,
     1                  k     , jstack, len   , lndxpt, locppc, lson  , 
     2                  lsonp , ncol  , nrow  , paroff, pppont, skip  , 
     3                  son   , sonbeg, sonoff, sonppb, sonppc, stkpnt

      integer           nnode , nstksv

      integer           idummy(1)

      logical           qslide
 
      double precision  temp
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      integer           xdslni
 
      external          xdslni,
     1                  xdslf5, xdslf7, xdslf9, xdslvr, xislvr
c
c  =====================================================================

      nstksv = nstack
c.debug
c     write(6,'("at start of assembling - nassmb = ", i8)') nassmb
c.debug

      qslide = .false.
      if ( nassmb .gt. 0 ) then
        if(rstbeg(nstack) .gt. 0 ) qslide = .true.
      endif
c.debug
c     write(6,'("at start of assembling - qslide = ", l8)') qslide
c.debug
 
c     ------------------------------
c     ... assemble the current front
c     ------------------------------
 
      if ( .not. qslide ) then
 
c         -----------------------------------------------
c         no stacked update matrices in-core to assemble,
c         zero out frontal matrix
c         -----------------------------------------------
c.debug
c     write(6,'("before do 10")')
c.debug
 
          do k = stfrnt, stfrnt + matsiz - 1
              rstack(k) = 0.d0
          enddo

          jstkst = 1
c.debug
c     write(6,'("after  do 10")')
c.debug
 
      else
 
c         --------------------------------------------------------
c         assemble the stacked update matrices, sliding the first
c         into the space for the frontal matrix and scatter/adding
c         the remaining update matrices.
c         --------------------------------------------------------
 
          son    = stknod(nstack)
c.debug
c     write(6,'("sliding first son - nstack, son = ", 2i8)') 
c    1                                nstack, son
c.debug
 
c         ------------------------------------------------------
c         establish pointers into stack entry for easier access
c         ------------------------------------------------------
 
          lson   = istack ( 1, nstack )
          lsonp  = istack ( 3, nstack )
          lndxpt = istack ( 2, nstack )
          sonppc = rstbeg (nstack)
          sonppb = sonppc + lsonp
          sonbeg = sonppb + lsonp*(lsonp+1)/2 + lsonp*lson
          endson = sonbeg + lson *(lson +1)/2
c.debug
c     write(6,'("lson, lsonp, lndxpt, sonppc, sonppb, sonbeg    = ",
c    1    7i8)')  lson, lsonp, lndxpt, sonppc, sonppb, sonbeg    
c     write(6,'("endson                                         = ",
c    1    7i8)')  endson                                         
c.debug
 
c         ---------------------------------
c         zero out area for frontal matrix
c         ---------------------------------
 
          do k = endson, stfrnt+matsiz-1
              rstack(k) = 0.d0
          enddo
 
c         ----------------------------------------
c         place regular update columns into front
c         ----------------------------------------
 
          sonoff = sonbeg - 1
          paroff = stfrnt - 1
c.debug
c     write(6,'("lson,lndxpt,loclfr,lfront,pospon,sonoff,paroff = ",
c    1    7i8)')  lson,lndxpt,loclfr,lfront,pospon,sonoff,paroff 
c     call xislp3 ( 'lindxl', lson, lindxl(lndxpt), 6 )
c     call xdslp5 ( 'son', endson-sonbeg, rstack(sonoff+1), 6 )
c.debug
 
          call xdslf5 ( lson, lindxl (lndxpt), loclfr,
     2                  rstack, sonoff, paroff )

c.debug
c     call xdslp5 ( 'front after sliding first child',
c    1              matsiz, rstack(stfrnt), 6 )
c.debug
 
          fctops = fctops + lson * ( lson + 1 ) / 2
 
          nstack = nstack - 1

          jstkst = 2

      end if

c     ------------------------------------
c     ... add in the other update matrices
c     ------------------------------------
 
      actual = bgxsup
      locppc = 0
      skip   = pospon
      pppont = stfrnt

      do 90 jstack = jstkst, nassmb
 
          son    = stknod(nstack)
c.debug
c     write(6,'("inside do 90 loop - jstack, son = ", 2i8)') 
c    1                                jstack, son
c.debug

c         ----------------------------------------
c         establish pointers into stack entry for
c         easier access
c         ----------------------------------------
 
          stkpnt = rstbeg (nstack)

          lson   = istack (1, nstack)
          lsonp  = istack (3, nstack)
          lndxpt = istack (2, nstack)
c.debug
c     write(6,'("stkpnt, lson, lsonp, lndxpt = ", 4i8)') 
c    1            stkpnt, lson, lsonp, lndxpt 
c.debug
 
          if ( stkpnt .gt. 0 ) then
 
c             ----------------------------
c             ... update matrix is in-core
c             ----------------------------
 
c.debug
c     write(6,'("update is in core")')
c.debug

              sonppc = stkpnt
              sonppb = sonppc + lsonp
              sonbeg = sonppb + lsonp*(lsonp+1)/2 + lsonp*lson
              endson = sonbeg + lson *(lson +1)/2
c.debug
c     write(6,'("stkpnt, lson, lsonp, lndxpt = ", 4i8)') 
c    1            stkpnt, lson, lsonp, lndxpt 
c     write(6,'("sonppc, sonppb, sonbeg, endson = ", 4i8)') 
c    1            sonppc, sonppb, sonbeg, endson 
c     call xdslp5 ( 'update', endson-sonbeg+1, rstack(sonbeg), 6 )
c.debug
 
c             --------------------------------------------
c             add components of regular update into front
c             --------------------------------------------
 
              call xdslf7 ( lson, lindxl (lndxpt),
     1                      rstack (sonbeg), loclfr, 
     2                      rstack (stfrnt))
c.debug
c     write(6,'("after xdslf7")')
c.debug
 
              nwfree = stkpnt
 
          else
 
c             ----------------------------------
c             ... update matrix is out-of-core.
c                 compute number of columns that
c                 can be read in at this time
c             ----------------------------------
 
              skip   = lsonp + lsonp*(lsonp+1)/2 + lsonp*lson
              iopos2 = istack(4, nstack) + skip 

              nrow   = lson
              len    = locbfr
              temp   = max ( ( .5 + nrow ) ** 2 - 2. * len, 0. )
              ncol   = .5 + nrow - sqrt ( temp )
              bfrlen = nrow*ncol - ncol*(ncol-1)/2
c.debug
c     write(6,'("update is out of core")')
c     write(6,'("iopos2, nrow, ncol, bfrlen = ", 4i8)') 
c    1            iopos2, nrow, ncol, bfrlen 
c.debug
 
              call xdslw1 ( wafil2, 2, idummy, idummy, ocbufr, 
     1                      iopos2, bfrlen, ierr)
              if ( ierr .ne. 0 ) go to 8900
c.debug
c     write(6,'("in xdslai - nstack, iopos2, bfrlen = ", 3i8)') 
c    1            nstack, iopos2, bfrlen
c     write(6,'("ocbufr(1), ocbufr(bfrlen) = ", 1p,2d25.15)') 
c    1            ocbufr(1), ocbufr(bfrlen)
c.debug
 
              iopos2 = iopos2 + bfrlen
              bfrpos = 1
 
c             --------------------------------------------
c             add components of regular update into front
c             --------------------------------------------
 
              call xdslf8 ( lson, lindxl (lndxpt),
     1                      wafil2, iopos2, ocbufr, bfrpos,
     2                      bfrlen, locbfr, loclfr,
     3                      rstack (stfrnt), ierr)
c.debug
c     write(6,'("after xdslf8")')
c.debug
              if ( ierr .ne. 0 ) go to 8900
 
              watrn2 = watrn2 + iopos2 - istack (4, nstack)
     1                        - ( sonbeg - sonppc )
 
          end if
 
          fctops = fctops + lson * ( lson + 1 ) / 2

          nstack = nstack - 1
          locppc = locppc + lsonp
c.debug
c     write(6,'("after adding next son")')
c     call xdslp5 ( 'front', matsiz, rstack(stfrnt), 6 )
c.debug
 
   90 continue
c.debug
c     write(6,'("after 90 continue - nstack = ", i8)') nstack 
c.debug

c     --------------------------------------------
c     ... end of adding update matrices into front
c     --------------------------------------------
 
c  =====================================================================

c     ---------------------------------------------------
c     ... add in contributions from the original matrices
c         a and b
c     ---------------------------------------------------
c.debug
c     write(6,'("before adding in org. entries")')
c     call xdslp5 ( 'front', matsiz, rstack(stfrnt), 6 )
c.debug

      call xdsla1 ( mincor, diag  , xladj , ladjl , lnza  ,
     1              bassmb, sigma , bmxtyp, bdiag , bxladj, bladjl,
     2              lnzb  , loclfr, nodbgn, nodnxt, pospon,
     3              rstack(stfrnt), sqfil1, sqtrn1, sqfil3, sqtrn3,
     4              mxvlib, mxvlrb, fctops, ierr )
c.debug
c     write(6,'("after  adding in org. entries")')
c     call xdslp5 ( 'front', matsiz, rstack(stfrnt), 6 )
c.debug

      if ( ierr .ne. 0 ) go to 8800
 
c     ------------------------------------------------------------
c     initialize the sup vector for the original size of the front
c     ------------------------------------------------------------
 
      actual = bgxsup
      do k = nodbgn, nodnxt - 1
          sup (actual) = k
          actual       = actual + 1
      enddo
c.debug
c     write(6,'("after 130 - pospon = ", i8)') pospon
c.debug
c.debug
c     write(6,'("after 130 continue - nstack = ", i8)') nstack 
c.debug

      if ( pospon .eq. 0 ) return
 
c  =====================================================================

c     ----------------------------------------------------------------
c     ... postponed columns are present
c         1.  expand front to make room for postponed columns after
c             the first nodnxt-nodbgn columns.  note:  this also means
c             expanding the first nodnxt-nodbgn columns to make room
c             for pospon extra rows.
c         2.  add in postponed columns from each of the children.
c     ----------------------------------------------------------------

c     ---------------------------------------------------
c     ... expand front to make room for postponed columns
c     ---------------------------------------------------

      nnode = nodnxt - nodbgn
c.debug
c     write(6,'("before expansion for postponed columns")')
c     call xdslp5 ( 'front', matsiz, rstack(stfrnt), 6 )
c.debug

      call xdsla2 ( nnode, pospon, loclfr, rstack(stfrnt) )
c.debug
c     write(6,'("after  expansion for postponed columns")')
c     call xdslp5 ( 'front', matsiz, rstack(stfrnt), 6 )
c.debug

c     ------------------------------------------------------
c     ... add in postponed columns from each of the children
c     ------------------------------------------------------

      call xdsla3 ( nnode , lfront, nassmb, nstksv, actual,
     1              pospon, lindxl, rstbeg, istack, sup   ,
     2              locbfr, ocbufr, wafil2, watrn2,
     3              stfrnt, rstack, fctops, ierr )
c.debug
c     write(6,'("after  inserting postponed columns")')
c     call xdslp5 ( 'front', matsiz, rstack(stfrnt), 6 )
c.debug

      return
 
c  =====================================================================
 
c     ------------------------------------------------
c     ... error trap for i/o error on sqfil1 or sqfil3
c     ------------------------------------------------
 
 8800 continue
      ierr = -2
      return
 
c     --------------------------------------
c     ... error trap for i/o error on wafil2
c     --------------------------------------
 
 8900 continue
      ierr = -1
      return
 
c  =====================================================================
 
      end
