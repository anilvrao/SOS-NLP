      subroutine xdslbo( mincor, diag  , xladj , ladjl , lnza  ,
     1                   bassmb, sigma , bmxtyp, bdiag , bxladj, bladjl,
     2                   lnzb  , sup   , nassmb, xlindx, lindxl,
     3                   bgxsup, dspmnt, n1    , n2    , n3    , matsiz,
     4                   mstack, nodbgn, nodnxt, nstack, 
     5                   stknod, rstbeg, istack, panell, panelu, lpanel,
     6                   sqfil1, sqtrn1, sqfil3, sqtrn3, wafil5,
     6                   iopos2, wafil2, watrn2, 
     7                   locbfr, ocbufr, mxvlib, mxvlrb, fctops, ierr )
 
c
c  purpose -- out-of-core assembly of the frontal matrix for the 
c             current supernode.  at completion the assembled front 
c             is on i/o file wafil5 starting at position 1.
c
c             if sigma .ne. 0. then assemble a - sigma * b.
c
c             this version is for structurally symmetric problems.
c
c  created   -- 22-may-98, rgg, derived from xdslbi
c  revisions -- 17-apr-03, jgl  tracks status of child processing
c                               for multiple postponed columns, to
c                               remove buffer overrun problem
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
c      lpanel -- length of panel array
c      sqfil1 -- i/o file holding a 
c      sqfil3 -- i/o file holding b 
c      wafil5 -- i/o file holding assembled front
c      wafil2 -- i/o file holding spilled parts of the stack
c
c  working storge  --
c
c      panell -- space to hold current section of front - lower tri
c      panelu -- space to hold current section of front - upper tri
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
c      iopos2 -- i/o position on wafil2.
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
c                = -3 i/o error on wafil5
c
c  subprograms called
c
c      xdsla4 -- add an update matrix to the frontal matrix
c      xdsla5 -- add in original matrix entries to the frontal matrix
c      xdsla7 -- insert postponed columns phase 1
c      xdsla8 -- insert postponed columns phase 2
c
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           bmxtyp, n1    , n2    , n3    , matsiz, mstack,
     1                  nassmb, nodbgn, nodnxt, nstack, 
     2                  bgxsup, dspmnt, wafil5, wafil2,
     3                  locbfr, mincor, sqfil1, sqfil3, 
     4                  iopos2, lpanel, ierr
 
      integer           xladj(*),               ladjl(*),
     1                  bxladj(*),              bladjl(*),
     3                  xlindx(*),              lindxl(*),
     4                  stknod(*),              rstbeg(*),
     5                  istack(4,*),            sup (*),
     6                  mxvlib(*)
 
      logical           bassmb, qfirst
 
      double precision  sigma , fctops, sqtrn1, sqtrn3, watrn2
 
      double precision  diag(*),                lnza(*),
     1                  bdiag(*),               lnzb(*),
     2                  panell(*),              panelu(*),
     3                  ocbufr(*),              mxvlrb(*)
 
c     -------------------
c     ... local variables
c     -------------------

      integer           actual, colbgn, colend, i     , iopos5, 
     1                  k     , jstack, lndxpt, lpan  , lson  ,
     2                  lsonp , lstcol, ncol  , nn1   , nrow  ,
     3                  nstksv, phase , skip  , son   , iopsv2

      integer           actsav, lfront, loclfr, lpano ,
     1                  offset, orgsiz, pospon, poffst

      integer           lstchc

      logical           midchl

      integer           idummy(1)
 
      double precision  temp
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      integer           xdslni
 
      external          xdslni,
     1                  xdslf5, xdslf7, xdslf9, xdslvr, xislvr
c
c  =====================================================================

c.debug
c     write(6,'("at start of assembling - nassmb = ", i8)') nassmb
c.debug

      pospon = n2
      loclfr = n1 + n3
      lfront = n1 + n2 + n3
      orgsiz = loclfr * ( loclfr + 1 ) / 2

      nstksv = nstack       
      qfirst = .true.

      phase  = 1
      colend = 0
      poffst = 0
      offset = 0

      if ( n2 .eq. 0 ) then
          lstcol = n1 + n3
      else
          lstcol = n1
      end if

      iopos5 = 1
 
c     ------------------------------------------------------------
c     initialize the sup vector for the original size of the front
c     ------------------------------------------------------------
 
      actual = bgxsup 
      do k = nodbgn, nodnxt - 1
          sup (actual) = k
          actual       = actual + 1
      enddo
c.debug
c     call xislp3 ( 'sup after do 10 in xdslbo', actual-1, sup, 6 )
c.debug

      actsav = actual

c     ---------------------------------------------------------
c     ... compute the end of i/o file wafil2 after the assembly
c         is complete
c     ---------------------------------------------------------

      if ( nassmb .eq. 0 ) then
 
          iopsv2 = iopos2

      else

          if ( nstack .ge. 1 ) then
              iopsv2 = istack(4,nstack)
          else
              iopsv2 = 1
          end if

      end if

      do jstack = 1, nassmb-1
          iopsv2 = min ( iopsv2, istack(4,nstack-jstack) ) 
      enddo
c.debug
c     write(6,'("in xdslbo - iopsv2 = ", i8)') iopsv2
c.debug

c  =====================================================================

c     -------------------------------------------------------
c     ... assembly will be performed in 3 phases
c         phase 1.  assemble columns of current front
c         phase 2.  assemble postponed columns
c         phase 3.  assemble columns of update portion
c         note:  if not postponed columns then there is no
c                phase 2 and phases 1 and 3 will be performed 
c                together
c     -------------------------------------------------------

  100 continue
      colbgn = colend + 1
    
c     -------------------------------------------------------
c     ... compute the number of columns that can be assembled
c     -------------------------------------------------------

      nrow = lfront 

      ncol = lstcol - colend
c.debug
c     ncol = min ( lstcol - colend, 20 )
c.debug

      if ( ncol * nrow - ( ncol * ( ncol-1 ) ) / 2 .le. lpanel ) then

c         ---------------------------------------------------
c         ... remainder of this phase can be done in one shot
c         ---------------------------------------------------

          colend = lstcol
          go to 200

      end if

c     -------------------------------------------------------------
c     ... solve quadratic equation on storage for number of columns
c     -------------------------------------------------------------

      temp = - nrow - 0.5 
      temp = - temp - sqrt ( temp**2 - 2.0 * lpanel )
c.debug
c     write(6,'("after quadratic - temp = ", f15.1)') temp
c.debug

      ncol   = temp

c     ------------------------------------------------------
c     ... now assemble the portion of the front from columns
c         colbgn to colend
c     ... start by zeroing out the panels
c     ------------------------------------------------------

  200 continue
      colend = colbgn + ncol - 1
      lpan   = nrow * ncol - ( ncol * ( ncol-1 ) ) / 2
c.debug
c     write(6,'("at 200 - colbgn, colend, ncol, nrow   = ", 4i8)')
c    1                     colbgn, colend, ncol, nrow 
c     write(6,'("at 200 - phase , lstcol, lpan, lpanel = ", 4i8)')
c    1                     phase , lstcol, lpan, lpanel
c     write(6,'("at 200 - nstack, nassmb               = ", 4i8)')
c    1                     nstack, nassmb              
c.debug

c.debug
c     if ( lpan .gt. lpanel ) then
c         write(6,'("oops - nrow, ncol, lpan, lpanel = ", 4i8)')
c    1                       nrow, ncol, lpan, lpanel 
c         stop
c     end if
c.debug

      panell(1:lpan) = 0.d0
      panelu(1:lpan) = 0.d0

c  =====================================================================
 
c     --------------------------------------------------
c     ... assemble contribution for each of the children
c     --------------------------------------------------
 
      if ( phase .eq. 2 ) go to 400

      nstack = nstksv
 
      do jstack = 1, nassmb
 
          son    = stknod(nstack)
c.debug
c     write(6,'("inside do 300 loop - jstack, nstack, son = ", 3i8)') 
c    1                                 jstack, nstack, son
c.debug

c         ----------------------------------------
c         establish pointers into stack entry for
c         easier access
c         ----------------------------------------
 
 
          lson   = istack (1, nstack)
          lsonp  = istack (3, nstack)
          lndxpt = istack (2, nstack)
 
c.debug
c     write(6,'("lson, lsonp, lndxpt = ", 4i8)') 
c    1            lson, lsonp, lndxpt 
c.debug

c         ------------------------------------------------
c         ... set i/o position just past postponed columns
c         ------------------------------------------------
 
          skip   = lsonp + lsonp*(lsonp+1)/2 + lsonp*lson
          iopos2 = istack(4, nstack) + skip 
c.debug
c     write(6,'("nstack, iopos2      = ", 4i8)') 
c    1            nstack, iopos2                      
c.debug

c         ----------------------------------------------------
c         ... perform out-of-core assembly for this child into
c             this section of the front - lower triangle
c         ----------------------------------------------------

          call xdsla4 ( lson  , offset, lindxl (lndxpt),
     1                  wafil2, iopos2, watrn2, ocbufr, locbfr, 
     2                  colbgn, colend, loclfr, panell, fctops, ierr)

c.debug
c     write(6,'("after xdsla4")')
c     call xdslp5 ( 'panell', lpan, panell, 6 )
c.debug
          if ( ierr .ne. 0 ) go to 8700

c         ------------------------------------------------
c         ... set i/o position just past postponed columns
c         ------------------------------------------------
 
          skip   = lsonp + lsonp*(lsonp+1)/2 + lsonp*lson
     1           + lson*(lson+1)/2 
     2           + lsonp*(lsonp+1)/2 + lsonp*lson
          iopos2 = istack(4, nstack) + skip 
c.debug
c     write(6,'("nstack, iopos2      = ", 4i8)') 
c    1            nstack, iopos2                      
c.debug

c         ----------------------------------------------------
c         ... perform out-of-core assembly for this child into
c             this section of the front - upper triangle
c         ----------------------------------------------------

          call xdsla4 ( lson  , offset, lindxl (lndxpt),
     1                  wafil2, iopos2, watrn2, ocbufr, locbfr, 
     2                  colbgn, colend, loclfr, panelu, fctops, ierr)

c.debug
c     write(6,'("after xdsla4")')
c     call xdslp5 ( 'panelu', lpan, panelu, 6 )
c.debug
          if ( ierr .ne. 0 ) go to 8700

          nstack = nstack - 1
 
      enddo
 
c  =====================================================================

c     ---------------------------------------------------
c     ... add in contributions from the original matrices
c         a and b.  phase 1 only
c     ---------------------------------------------------
c.debug
c     write(6,'("before adding in org. entries")')
c     call xdslp5 ( 'panell', lpan, panell, 6 )
c     call xdslp5 ( 'panelu', lpan, panelu, 6 )
c.debug

      if ( phase .eq. 1 ) then 

          i     = min ( colend, n1 )
          lpano = ( nrow - pospon ) * ncol - ( ncol * ( ncol - 1 ) ) / 2

          call xdslb5 ( mincor, nodbgn, nodnxt, diag  , xladj,  
     1                  ladjl , lnza  ,
     2                  bassmb, sigma , bmxtyp, bdiag , bxladj,
     3                  bladjl, lnzb  , colbgn, i     , poffst,
     4                  lpano , loclfr, orgsiz, panell, panelu, 
     5                  sqfil1, sqtrn1, sqfil3, sqtrn3, 
     6                  qfirst, mxvlib, mxvlrb, fctops, ierr )

          if ( ierr .ne. 0 ) go to 8800

          qfirst = .false.
          poffst = poffst + lpano

      end if
 
c  =====================================================================

c     ---------------------------------------
c     ... deal with postponed columns, if any
c     ---------------------------------------

  400 continue
c411      actual = actsav
c.debug
c     write(6,'("before adjusting for postponed columns")')
c     call xdslp5 ( 'panell', lpan, panell, 6 )
c     call xdslp5 ( 'panelu', lpan, panelu, 6 )
c.debug

      if ( phase .eq. 1 .and. n2 .gt. 0 ) then

c         -------------------------------------------------------
c         ... phase 1 - expand panel and insert postponed columns
c         -------------------------------------------------------

          nstack = nstksv
          nn1    = n1 - colbgn + 1

          call xdslb7 ( nn1   , n2    , n3    , loclfr, lpan  , 
     1                  colbgn, colend, panell, panelu,
     2                  nstack, nassmb, stknod, istack, lindxl,
     3                  wafil2, watrn2, locbfr, ocbufr, ierr  )

      else if ( phase .eq. 2 ) then

c         -------------------------------------------------
c         ... phase 2 - only columns are the postponed ones
c         -------------------------------------------------

          nstack = nstksv

          call xdslb8 ( n1    , n2    , n3    , loclfr, lpan  , 
     1                  colbgn, colend, panell, panelu,
     2                  nstack, nassmb, stknod, istack, lindxl,
     3                  wafil2, watrn2, locbfr, ocbufr, 
     4                  actual, sup   , ierr  , lstchc, midchl )

      end if
      if ( ierr .ne. 0 ) go to 8700
 
c  =====================================================================

c     ---------------------------------------
c     ... write out this section of the front
c     ---------------------------------------
c.debug
c     write(6,'("before writing out the panel-iopos5, lpan, matsiz = "
c    1     3i8)') iopos5, lpan, matsiz
c     call xdslp5 ( 'panell', lpan, panell, 6 )
c     call xdslp5 ( 'panelu', lpan, panelu, 6 )
c.debug

      call xdslw2 ( wafil5, 2, idummy, idummy, panell, 
     1              iopos5, lpan, ierr )

      if ( ierr .ne. 0 ) go to 8900

      call xdslw2 ( wafil5, 2, idummy, idummy, panelu, 
     1              iopos5+matsiz, lpan, ierr )

      if ( ierr .ne. 0 ) go to 8900

      iopos5 = iopos5 + lpan
 
c  =====================================================================

c     -------------------------------------------
c     ... loop back for next section of the front
c     -------------------------------------------

      lfront = lfront - ncol
      if ( phase .ne. 2 ) loclfr = loclfr - ncol

      if ( colend .eq. lstcol ) then
          if ( phase .eq. 1 .and. n2 .ne. 0 ) then
              lstcol = n1 + n2
              phase  = 2
              midchl = .false.
              lstchc = 0
          else if ( phase .eq. 2 ) then
              lstcol = n1 + n2 + n3
              phase  = 3
              offset = n2
          end if
      end if

      if ( phase .eq. 1 .and. colend .gt. n1 ) phase = 3

c.debug
c     write(6,'("before test on branching back for next section")')
c     write(6,'("colend, n1, n2, n3, n1+n2+n3         = ", 5i8)')
c    1            colend, n1, n2, n3, n1+n2+n3
c     write(6,'("phase, lstcol, colend, lfront, loclfr = ", 5i8)')
c    1            phase, lstcol, colend, lfront, loclfr
c.debug

      if ( colend .lt. n1 + n2 + n3 ) go to 100
 
c  =====================================================================

      nstack = nstksv - nassmb
      iopos2 = iopsv2
c.debug
c     write(6,'("leaving xdslbo - iopos2 = ", i8)') iopos2
c.debug

      return
 
c  =====================================================================
 
c     --------------------------------------
c     ... error trap for i/o error on wafil2
c     --------------------------------------
 
 8700 continue
      ierr = -1
      return
 
c     ------------------------------------------------
c     ... error trap for i/o error on sqfil1 or sqfil3
c     ------------------------------------------------
 
 8800 continue
      ierr = -2
      return
 
c     --------------------------------------
c     ... error trap for i/o error on wafil5
c     --------------------------------------
 
 8900 continue
      ierr = -3
      return
 
c  =====================================================================
 
      end
