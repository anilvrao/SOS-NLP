      subroutine   xdslf2   ( pvttol, mincor, fpnlsz, fupdsz, cmajor,
     1                        neqns , diag  , xladj , ladjl , lnza  , 
     2                        sigma , bmxtyp, bdiag , bxladj, bladjl,
     3                        lnzb  , extfil, zpcntl, npcntl, perm  , 
     4                        nsuper, xsup  , sup   , nassmb, npanel,
     5                        xpanel, xlindx, lindxl, lndxg1, nsninx, 
     6                        iwork , rwork , stkbgn, lstack, lnzbgn, 
     7                        llnz,   qincor, mstack, stknod, rstbeg, 
     8                        istack, sqfil1, sqfil3, wafil1, wafil2, 
     9                        wafil5, lbuffr, ocbufr, mxvlrb, mxvlib, 
     A                        fnzlf , nsnind, inrtia, pvtblk, walen1, 
     B                        walen2, walen5, sqtrn1, sqtrn3, watrn1, 
     C                        watrn2, watrn5, mpanel, fctops, slvops, 
     D                        error , ppfmon, output, lrwork, needst )

 
c
c  purpose -- to get a ldlt factorization of a symmetric sparse matrix
c             pivoting using the multifrontal method.
c             a symbolic factorization has been performed and
c             the original row indices are available in local form
c             with respect to the supernodal parent frontal matrix.
c             the supernode indices are available in local form with
c             respect to the parent supernode.
c
c             if sigma .ne. 0. then xdslf2 factors a - sigma * b
c             (used by lanczos) where bmxtyp determines the matrix
c             type of b (0 - null, 1 - general, 4 - general but same
c             structure as a, 3 - diagonal, and 4 - identity).
c
c  created   -- 06-jul-87, cca
c  revisions -- 06-jul-87, cca
c               30-jan-89, rgg, changed xladj to be referenced by node
c                               number to be compatible with matrix
c                               value input requirements.
c               23-feb-89, rgg, changed name to xdslf2 and included
c                               out-of-core features.
c               17-apr-89, djp, incorporated incremental condition
c                               estimator for incore version.
c               17-jun-89, rgg, out-of-core features removed for convex
c                               version
c               03-jul-89, rgg, removed incremental condition
c                               estimator for incore version.
c               11-jul-89, rgg, added -sigma*b assembly
c               12-sep-89, rgg, reinstalled out-of-core features
c               08-feb-91, rgg, added use of xdslil
c               12-feb-91, rgg, allowed disabling of wafil1 and wafil2
c               04-nov-92, rgg, improved i/o performance in xdslf4,
c                               xdslf8, and xdslfh by use of the ocbufr
c                               array
c               04-dec-92, rgg, added minimum core processing (both
c                               a and l out-of-core)
c               19-may-93, rgg, removed generation of lindxg2 when
c                               pvttol = 0. (it is then identical
c                               to lindxg1).
c               08-jun-95, rgg, altered for new pivot selection in
c                               xdslfb
c               23-feb-97, rgg, converted to panel concept
c               18-mar-97, rgg, converted to out-of-core processing for
c                               fronts
c               28-oct-97, rgg, added controls on amount of fill, zero 
c                               and negative pivots
c               30-jul-98, rgg, separated out-of-core front processing
c                               to i/o file wafil5
c               01-sep-98, rgg, 32 bit integer mods
c               08-mar-00, rgg, allowed for possible increase of factor
c                               panel size for factorization of out of
c                               memory fronts (xdsleo) to reduce i/o.
c               01-jun-00, rgg, forced spilling the entire stack
c               25-jul-00, dkw, adjusted storage size for temp1
c               31-jan-01, jgl, error handling modified
c               02-jun-03, jgl, repaired 08-mar-00 fixes to properly
c                               reduce panel size and maximally increase
c                               panel size
c
c  variables
c
c      pvttol -- pivoting tolerance
c      mincor -- minimum core processing switch
c      fpnlsz -- factorization panel size
c      fupdsz -- factorization update panel size
c      cmajor -- column major flag
c      neqns  -- number of equations in the matrix system
c      diag   -- array to hold first the diagonal entries of the
c                original matrix and later those of the factor
c      xladj  -- pointers into indices and entries of the matrix
c      ladjl  -- local indices of the column indices of original matrix
c      lnza   -- column entries of the original matrix
c      sigma  -- shift value
c      bmxtyp -- matrix type of b
c                0 - null matrix
c                1 - general sparse
c                2 - diagonal
c                3 - identity
c                4 - general sparse but same structure as a
c      bdiag  -- diagonal of b
c      bxladj -- xladj for b
c      bladjl -- ladjl for b
c      lnzb   -- column entries of the b matrixc
c      extfil -- amount of extra fill allowed.
c      zpcntl -- array for controlling zero pivots
c      npcntl -- array for controlling negative pivots
c      perm   -- new to old permutation array
c      nsuper -- number of supernodes
c      xsup   -- supernode pointer array. since the matrix is assumed
c                to be ordered by a post-ordered traversal, the nodes
c                for supernode i are xsup(i), ... , xsup(i+1)-1
c      nassmb -- integer work array giving the number of assemblies
c                of stack matrices for each supernode
c      xlindx -- pointers to the column indices of the supernodes
c                in the factor
c      lindxl -- local indices of the supernodes with respect to
c                their supernode parent
c      lndxg1 -- global indices of the supernodes with respect to
c                their supernode parent (pre-factorization)
c      nsninx -- current allowed length of lndxg2 array stored in
c                iwork
c      mstack -- maximum depth of the stack, dimension of
c                stknod and rstbeg
c      sqfil1 -- sequential i/o file for a
c      sqfil3 -- sequential i/o file for b
c      wafil1 -- word addressable i/o file for lnz
c      wafil2 -- word addressable i/o file for stack
c      wafil5 -- word addressable i/o file for out-of-core processing
c                for a front
c      lbuffr -- length of ocbufr
c      ppfmon -- logical flag to monitor panel complete pivot
c                postponement/failure (or not)
c      output -- printer unit number
c
c  working variables --
c
c      stknod -- integer array to hold the supernode number of
c                each stacked matrix
c      rstbeg -- integer array to hold the starting location of each
c                stacked matrix
c      istack -- 2d integer array to hold the scalar information
c                associated with pivoting for a stacked matrix.
c      stkbgn -- location in rwork array where real stack starts.
c      lstack -- length of the stack space
c      lnzbgn -- location in rwork array where real lnz starts.
c      llnz   -- length of lnz space
c      qincor -- logical flag denoting whether lnz is in core or
c                on wafil1.
c
c  working storage --
c
c      iwork, rwork -- two work arrays which are equivalenced by
c                the call to xdslf2.  iwork contains the new global
c                indices (lndxg2) starting at the location 1.  rwork
c                contains the stack starting at location stkbgn.
c                stkbgn can be increased during the factorization to
c                avoid collisions between lndxg2 and stack.  the
c                nonzeros of l (lnz) is stored starting at lnzbgn
c                if qincor = .true..  at the first collision of working
c                space, lnz is written to wafil1 and maintained there
c                for the remainder of the factorization.
c      ocbufr -- buffer for out-of-core processing.  length
c                is lbuffr.
c      mxvlrb -- buffer for out-of-core processing of org. matrix
c                entries.  length is mxanzf. not used if mincor = 0.
c      mxvlib -- buffer for out-of-core processing of org. matrix
c                entries.  length is mxanzf. not used if mincor = 0.
c
c  output variables --
c
c      diag   -- diagonal entries for the factor
c      xsup   -- supernode pointer array. on output, this now points
c                to the first panel number of each supernode.
c      npanel -- number of panels
c      xpanel -- pointer to first node of each panel
c      iwork  -- contains new global indices starting at location 1.
c      fnzlf  -- number of nonzeros in the lower triangle of l.
c      nsnind -- number of supernodal indices (the length of
c                lindxl on entry and the length of lndxg2 on
c                return)
c      inrtia -- matrix inertia
c                inrtia(1) = no. of pos. eigenvalues
c                inrtia(2) = no. of neg. eigenvalues
c                inrtia(3) = no. of zero eigenvalues
c      pvtblk -- integer n vector containing information on the
c                size of the diagonal blocks in the factorization.
c                in the occurrence of a 2x2 block the second
c                corresponding component is equated to the index
c                into a for the off diagonal element in the 2x2
c                pivot.
c      walen1 -- length of i/o file wafil1
c      walen2 -- length of i/o file wafil2
c      walen5 -- length of i/o file wafil5
c      sqtrn1 -- amount of i/o transfered to and from i/o file sqfil1
c      sqtrn3 -- amount of i/o transfered to and from i/o file sqfil3
c      watrn1 -- amount of i/o transfered to and from i/o file wafil1
c      watrn2 -- amount of i/o transfered to and from i/o file wafil2
c      watrn5 -- amount of i/o transfered to and from i/o file wafil5
c      mpanel -- size of largest panel for solve stage
c      needst -- size of stack needed for incore
c      fctops -- number of floating point operations to compute
c                factorization
c      slvops -- number of floating point operations to apply
c                the factorization in a numeric solution.
c      error  -- error flag, error  =   0, success
c                                   =  -1, lstack too small
c                                   =  -2, exactly singular matrix
c                                          (when pivoting)
c                                   =  -4, mstack too small
c                                   =  -5, i/o error on wafil1
c                                   =  -6, i/o error on wafil2
c                                   =  -7, ran out of storage for lnz
c                                          with wafil1 disabled
c                                   =  -8, ran out of storage for stack
c                                          with wafil2 disabled
c                                   =  -9, i/o error on sqfile
c                                   = -10, exceeded amount of allowed
c                                          fill
c                                   = -11, encountered negative pivot
c                                          while monitoring for such
c                                          (when not pivoting)
c                                   = -13, i/o error on wafil5
c                                   = -14, exact zero diagonal
c                                   = -15, failure in panel pivot scheme
c                                          (when pivoting)
c                                          (when not pivoting)
c
c  subprograms called
c
c      xdslfg -- copies factored entries of front to arrays diag and
c                lnz ( starting at lnzbgn in array rwork).
c      xdslfh -- copies factored entries of front to array diag and off-
c                diagonal entries on word addressable i/o file wafil1.
c      xdslmv -- floating point data move
c      icopy  -- vector copy
c      xdslf3 -- spills stack of update matrices to wafil2.
c      xdslai -- in-core frontal matrix assembly routine
c      xdslao -- out-of-core frontal matrix assembly routine
c      xdslei -- in-core elimination routine
c      xdsleo -- out-of-core elimination routine
c      xislfp -- diagnostic output routine
c
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           fpnlsz, fupdsz, cmajor, 
     1                  neqns,  bmxtyp, nsuper, mstack, lstack, nsninx,
     2                  stkbgn, lnzbgn, llnz,   wafil1, wafil2, 
     3                  npanel, nsnind, lbuffr, mincor, sqfil1, sqfil3,
     4                  mpanel, error , wafil5, output, lrwork
 
      integer           xladj (*),     ladjl (*),      perm  (*),
     1                  bxladj(*),     bladjl(*),      pvtblk(*),
     2                  xsup  (*),     xpanel(*),      nassmb(*),
     3                  xlindx(*),     lindxl(*),      inrtia(3),
     4                  stknod(*),     rstbeg(*),
     5                  istack(4,*),   lndxg1(*),      needst   ,
     6                  iwork (*),     sup   (*),      mxvlib(*)
 
      logical           qincor, ppfmon
 
      double precision  sigma , pvttol, fctops, slvops,
     1                  walen1, watrn1, walen2, watrn2, sqtrn1, sqtrn3,
     2                  walen5, watrn5, fnzlf , extfil
 
      double precision  diag(*),            lnza(*),
     1                  bdiag(*),           lnzb(*),
     2                  rwork(lrwork),      ocbufr(*),
     3                  mxvlrb(*),
     4                  zpcntl(*),          npcntl(*)

c     -------------------
c     ... local variables
c     -------------------
 
      integer           actual, bgdxg2, bgxsup, dspmnt, fstore, hldnnd,
     1                  i     , iopos2, isuper, k     , kstack, ierr  ,
     2                  l     , length, lfront, lndxpt, loclfr, lpanel,
     3                  lson  , lsuper, matsiz, mxstck, iopsv2
 
      integer           lextra, n1    , n3    , ncolpp, nnode , nodbgn, 
     1                  nodnxt, nodoff, nstack, nwfree, swap  ,
     2                  panel , pnlsiz, pospon, son   , start , stfrnt,
     3                  stfree, stkoff, stkuse, stneed, temp1 , addfil,
     4                  updsiz, lavail, lnreqd, lenxtr

      integer           iopos1(2)
      integer           lnrwrk
 
      logical           bassmb, qicasm, qnpvt1, qzpvt1 

      double precision  ratio

c.timer
      double precision  t1, t2, t3, t4, t5, w1, w2, w3, w4, w5,
     1                  times(7), f1, f2, f3, fops(4)
c.timer
 
c     --------------------------------------------
c     ... local variables for panel pivot failures
c     --------------------------------------------

      integer           xpboxs, psboxs

      parameter       ( xpboxs = 16, psboxs = 10 )

      integer           izfail, rzfail,
     1                  inpexp (xpboxs), inpsiz (psboxs),
     2                  inzsiz (psboxs), rtpexp (xpboxs),
     3                  rtpsiz (psboxs), rtzsiz (psboxs)

c     --------------------
c     ... subprograms used
c     --------------------
 
      integer           xdslil, xdslni
 
      external          icopy,  xdslil, xdslw6, xdslmv, 
     1                  xdslni, xdslf3, xdslai, xdslao, xdslei, xdsleo
c
c  =====================================================================
 
c  --------------
c  initialization
c  --------------
 
      ierr      = 0
      rstbeg(1) = 1
      stfree    = 1
      stkoff    = stkbgn - 1
      iopos1(1) = 1
      iopos1(2) = 1
      iopos2    = 1
      mxstck    = 0
      nwfree    = 1
      nstack    = 0
      fnzlf     = 0.
      npanel    = 0
      xpanel(1) = 1
      inrtia(1) = 0
      inrtia(2) = 0
      inrtia(3) = 0
      bgxsup    = 1
      bgdxg2    = 1
      actual    = 1
      mpanel    = 0
      fctops    = 0.
      slvops    = 0.
 
      sqtrn1    = 0.
      sqtrn3    = 0.
      walen1    = 0.
      watrn1    = 0.
      walen2    = 0.
      watrn2    = 0.
      walen5    = 0.
      watrn5    = 0.

      bassmb = ( sigma  .ne. 0. ) .and.
     1         ( bmxtyp .ge. 1  ) .and. ( bmxtyp .le. 4  )

      istack(1:4,1:mstack) = 0
      needst = 0

c.timer
      fops(1:4) = 0.d0
      do i = 1, 7
          times(i) = 1.0d-10
      enddo
c.timer

      addfil = 0
      qzpvt1 = .true.
      qnpvt1 = .true.

      do i = 2, 5
          zpcntl(i) = 0.
          npcntl(i) = 0.
      enddo

c  ----------------------------
c  pivot failure initialization
c  ----------------------------

      if ( ppfmon ) then
         izfail = 0
         rzfail = 0
         inpexp(1:xpboxs) = 0
         inpsiz(1:psboxs) = 0
         inzsiz(1:psboxs) = 0
         rtpexp(1:xpboxs) = 0
         rtpsiz(1:psboxs) = 0
         rtzsiz(1:psboxs) = 0
      endif

c  =====================================================================
 
c  ------------------------
c  loop over the supernodes
c  ------------------------
 
      do 1000 isuper = 1,nsuper
c.timer
          call xdslt1 ( t1, w1 )
c.timer
c.debug
c     write(6,'("at start of 1000 loop - isuper, nsuper, nassmb = ", 
c    1        3i8)') isuper, nsuper, nassmb(isuper)
c     call xislp3 ( 'nassmb', nsuper, nassmb, 6 )
c     call xislp2 ( 'istack', 4*nstack, istack, 6 )
c     call xislp3 ( 'sup', bgxsup-1, sup, 6 )
c.debug
 
c  =====================================================================
 
c         --------------
c         assembly stage
c         --------------
 
          pospon = 0
          k      = nstack
          do i = 1, nassmb (isuper)
             pospon = pospon + istack(3,k)
             k      = k - 1
          enddo
          nodbgn = xsup(isuper)
          nodnxt = xsup(isuper+1)
          n1     = nodnxt - nodbgn
          n3     = xlindx(isuper+1) - xlindx(isuper)

          nnode  = n1 + pospon
          lfront = n1 + pospon + n3
          loclfr = n1 + n3
          dspmnt = ((pospon)*(pospon+1))/2
     1           + (lfront-pospon)*pospon
          matsiz = (lfront*(lfront+1))/2
          pnlsiz = min ( fpnlsz, nnode )
          updsiz = max ( min ( fupdsz, loclfr ), 2 )

          if ( extfil .gt. 0. ) then

c             ---------------------------------------------
c             ... if monitoring amount of fill, compute the
c                 amount of extra fill so far
c             ---------------------------------------------

              k = nstack
              do i = 1, nassmb(isuper)
                  lson   = istack(1,k) + istack(3,k)
                  addfil = addfil + istack(3,k) * ( lfront - lson )
                  k = k - 1
              enddo

c             ------------------------------------------------
c             ... if too much fill and not the last supernode,
c                 then abort the factorization
c             ------------------------------------------------

              if ( ( addfil .gt. extfil ) .and.
     1             ( addfil .gt. 1000   ) .and.
     2             ( isuper .lt. nsuper )       ) then
                  error = -10
                  return
              end if

          end if
          
c         ------------------------------------------
c         ... test if front will be processed in- or 
c             out-of-core
c         ------------------------------------------

          fstore = matsiz + ( pnlsiz + updsiz ) * lfront 
     1           + xdslni(lfront)
c.debug
c     write(6,'("isuper, n1    , pospon, n3     = ", 4i8)')
c    1            isuper, n1    , pospon, n3    
c     write(6,'("nnode , lfront, loclfr, matsiz = ", 4i8)')
c    1            nnode , lfront, loclfr, matsiz
c     write(6,'("pnlsiz, fstore, lstack, llnz   = ", 4i8)')
c    1            pnlsiz, fstore, lstack, llnz  
c     write(6,'("updsiz                         = ", 4i8)')
c    1            updsiz                        
c.debug

          if ( fstore .le. ( lstack + llnz ) ) then
              qicasm = .true.
          else
              qicasm = .false.
          end if

c.debug
c     qicasm = .false.
c     write(6,'("qicasm                         = ",  l8)')
c    1            qicasm                        
c.debug
 
c         --------------------------------------------------
c         ... allocate storage for this front from the stack
c         --------------------------------------------------
 
          kstack = nstack

          if ( qicasm .and. nassmb(isuper) .ne. 0 ) then

              if ( rstbeg(nstack) .gt. 0 ) then 

                  l = 0
                  if ( istack(3,nstack) .ne. 0 ) then

c                    -----------------------------------------------
c                    ... update from last son has postponed columns.
c                        this area needs to be preserved.
c                    -----------------------------------------------
 
                     ncolpp = istack(3,nstack)
                     lsuper = istack(1,nstack)
                     l      = ncolpp*lsuper + (ncolpp*(ncolpp+1))/2
                     l      = l + ncolpp
    
                  end if

                  stfree = rstbeg ( nstack ) + l
c.debug
c     write(6,'("at ck pt 1 - stfree = ", i8)') stfree
c.debug
                  kstack = nstack - 1

              else

                  stfree = 1
c.debug
c     write(6,'("at ck pt 1.1 - stfree = ", i8)') stfree
c.debug

              end if

          end if

          stneed = stfree + fstore - 1
c.debug
c     write(6,'("before testing about lnz and stack at start")')
c.debug
 
c         ---------------------------------------
c         ... test if lnz needs to be written out
c         ---------------------------------------

          k = fnzlf + lfront*nnode - nnode*(nnode+1)/2
c.debug
c     write(6,'("qincor, qicasm                 = ", 2l8)')
c    1            qincor, qicasm
c     write(6,'("lstack, stneed, k     , llnz   = ", 4i8)')
c    1            lstack, stneed, k     , llnz   
c     write(6,'("wafil1, wafil2                 = ", 4i8)')
c    1            wafil1, wafil2                 
c.debug
 
          if ( qincor .and. 
     1         ( ( lstack .lt. stneed ) .or. ( k .gt. llnz ) .or.
     2           ( k .lt. 0 ) .or. ( .not. qicasm ) ) ) then 
 
c             -----------------------------
c             ... test if wafil1 is enabled
c             -----------------------------
 
              if ( wafil1 .lt. 0 ) then
                  error = -7
                  needst = stneed + llnz
                  return
              end if
 
c             ---------------------------
c             ... write lnz out to wafil1
c             ---------------------------
 
              k = fnzlf
              call xdslw6 ( wafil1, rwork(lnzbgn), iopos1, k, ierr )
 
              if ( ierr .ne. 0 ) then
                 if  ( ierr .eq. -1 ) then
                    error = -5
                 else
                    error = -99
                 end if
                 return
              endif
 
              lstack = lstack + llnz
              lnzbgn = 0
              llnz   = 0
              qincor = .false.

              call xdslw9 ( iopos1, k )
 
          end if
 
c         -------------------------------------
c         ... test if stack needs to be spilled
c         -------------------------------------
 
  100     continue
          if ( ( .not. qicasm ) .or. ( lstack .lt. stneed ) ) then
 
c             -----------------------------
c             ... test if wafil2 is enabled
c             -----------------------------
 
              if ( wafil2 .lt. 0 ) then
                  error = -8
                  needst = stneed
                  return
              end if
 
c             -------------------
c             ... spill the stack
c             -------------------
 
              kstack = nstack
c.debug
c     write(6,'("spilling the stack - kstack, nstack = ",2i8)') 
c    1            kstack, nstack
c     call xislp2 ( 'stknod',   nstack, stknod, 6 )
c     call xislp2 ( 'rstbeg',   nstack, rstbeg, 6 )
c     call xislp2 ( 'istack', 4*nstack, istack, 6 )
c     call xislp3 ( 'lindx2', bgdxg2-1, iwork, 6 )
c.debug
 
              call xdslf3 ( kstack, stknod, rstbeg, xlindx,
     1                      istack, mstack, iopos2, wafil2,
     2                      walen2, watrn2, rwork(stkbgn), ierr )
 
              if ( ierr .ne. 0 ) then
                 if  ( ierr .eq. -1 )  then
                    error = -6
                 else
                    error = -99
                 end if
                 return
              endif
 
              k      = stfree
              stfree = 1
c.debug
c     write(6,'("at ck pt 2 - stfree = ", i8)') stfree
c.debug
 
              if ( kstack .lt. nstack ) then
 
                  son    = stknod ( nstack )
                  ncolpp = istack ( 3, nstack )
                  lson   = ncolpp + istack ( 1, nstack )
                  l      = ncolpp + ( lson * ( lson + 1 ) ) / 2
                  k      = rstbeg ( nstack )
 
                  call xdslmv ( l, rwork, k+stkoff, stfree+stkoff )
 
                  rstbeg(nstack) = 1

                  stfree = stfree + lson*ncolpp - ( ncolpp*(ncolpp-1))/2
 
              end if
 
          end if
          
c         ------------------------------------------------
c         ... test again if front will be processed in- or 
c             out-of-core
c         ------------------------------------------------

          fstore = matsiz + ( pnlsiz + updsiz ) * lfront 
     1           + xdslni(lfront)
          stneed = stfree + fstore - 1
c.debug
c     write(6,'("isuper, n1    , pospon, n3     = ", 4i8)')
c    1            isuper, n1    , pospon, n3    
c     write(6,'("nnode , lfront, loclfr, matsiz = ", 4i8)')
c    1            nnode , lfront, loclfr, matsiz
c     write(6,'("pnlsiz, fstore, lstack, llnz   = ", 4i8)')
c    1            pnlsiz, fstore, lstack, llnz  
c     write(6,'("updsiz, stneed                 = ", 4i8)')
c    1            updsiz, stneed
c.debug

          if ( stneed .le. lstack ) then
              qicasm = .true.
          else if ( wafil1 .ge. 0 .and. wafil2 .ge. 0 ) then
              qicasm = .false.
          else
              if ( wafil2 .lt. 0 ) error = -8
              if ( wafil1 .lt. 0 .or. wafil5 .lt. 0 ) then
                error = -7
                needst = stneed
              endif
              return
          end if
c.debug
c     qicasm = .false.
c     write(6,'("qicasm                         = ",  l8)')
c    1            qicasm                        
c.debug
 
c  =====================================================================

c         ------------------------------
c         ... assemble the current front
c         ------------------------------
c.timer
          call xdslt2 ( t1, w1, t2, w2 )
          f1 = fctops
c.timer

          if ( qicasm ) then

c             --------------------
c             ... in-core assembly
c             --------------------
 
              swap   = stfree 
              stfrnt = swap   + xdslni ( lfront )
              panel  = stfrnt + matsiz
              temp1  = panel  + pnlsiz * lfront
              stkuse = temp1  + updsiz * lfront
c.debug
c     write(6,'("ica-panel, temp1, stfrnt, stkuse, lstack = ", 5i8)')
c    1                panel, temp1, stfrnt, stkuse, lstack 
c.debug

              mxstck = max ( mxstck, stkuse-1 )
 
c             ----------------------------
c             ... check for stack overflow
c             ----------------------------
c.debug
c     write(6,'("before in-core assembly - lstack, stkuse = ", 
c    1         2i8)') lstack, stkuse
c.debug 
              if ( lstack .lt. stkuse-1 ) then
                  error = -1
                  needst = stkuse
                  return
              endif
c.debug
c     write(6,'("before in-core assembling")')
c.debug
 
              if ( nstack .eq. 0) then 
                  nwfree = stfree
              else
                  if( rstbeg(nstack) .eq. 0 ) then
                    nwfree = stfree
                  else
                    nwfree = rstbeg(nstack)
                  endif
              end if
 
              call xdslai (    mincor, diag  , xladj , ladjl , lnza  ,
     1                 bassmb, sigma , bmxtyp, bdiag , bxladj, bladjl,
     2                 lnzb  , sup   , nassmb(isuper), xlindx, lindxl,
     3                 bgxsup, dspmnt, lfront, loclfr, matsiz,
     4                 mstack, nodbgn, nodnxt, nstack, nwfree, pospon,
     5                 stknod, rstbeg, istack, rwork(stkbgn),  stfrnt,
     6                 sqfil1, sqtrn1, sqfil3, sqtrn3, wafil2, watrn2,
     7                 lbuffr, ocbufr, mxvlib, mxvlrb, fctops, ierr   )
c.debug
c     write(6,'("after  in-core assembling")')
c.debug

              if  ( ierr .ne. 0 )  then
                 if ( ierr .eq. -1 ) then
                    error = -6
                 else if ( ierr .eq. -2 ) then
                    error = -9
                 else
                    error = -99
                 end if
                 return
              end if

          else

c             ------------------------
c             ... out-of-core assembly
c             ------------------------

              panel  = stfree
              lpanel = lstack - panel + 1
              stkuse = lstack
c.debug
c     write(6,'("oca-panel, lpanel, stkuse                = ", 5i8)')
c    1                panel, lpanel, stkuse                
c.debug

              mxstck = max ( mxstck, stkuse )

c.debug
c     write(6,'("before out-of-core assembling")')
c     write(6,'("isuper, nassmb(isuper) = ", 2i8)')
c    1            isuper, nassmb(isuper) 
c     call xislp3 ( 'sup pre assembly', bgxsup-1, sup, 6 )
c     call xislp3 ( 'lindx2 pre assembly', bgdxg2-1, iwork, 6 )
c.debug
 
              nwfree = stfree
c
              call xdslao (    mincor, diag  , xladj , ladjl , lnza  ,
     1                 bassmb, sigma , bmxtyp, bdiag , bxladj, bladjl,
     2                 lnzb  , sup   , nassmb(isuper), xlindx, lindxl,
     3                 bgxsup, dspmnt, n1    , pospon, n3    , matsiz,
     4                 mstack, nodbgn, nodnxt, nstack, stknod, rstbeg,
     5                 istack, rwork(panel+stkoff),    lpanel,
     6                 sqfil1, sqtrn1, sqfil3, sqtrn3, wafil5,
     7                 iopos2, wafil2, watrn2, lbuffr, ocbufr,
     8                 mxvlib, mxvlrb, fctops, ierr )
c.debug
c     call xislp3 ( 'sup post assembly', bgxsup+nnode-1, sup, 6 )
c     call xislp3 ( 'lindx2 post assembly', bgdxg2-1, iwork, 6 )
c.debug

              if  ( ierr .ne. 0 )  then
                 if ( ierr .eq. -1 ) then
                    error = -6
                 else if ( ierr .eq. -2 ) then
                    error = -9
                 else if ( ierr .eq. -3 ) then
                    error = -13
                 else
                    error = -99
                 end if
                 return
              end if

              watrn5 = watrn5 + matsiz
              if ( matsiz .gt. walen5 ) walen5 = matsiz

          end if
 
c  =====================================================================
 
c         -------------------------
c         scaling and updating step
c         -------------------------
 
          hldnnd        = nnode
          xsup (isuper) = npanel + 1
          nodoff = bgxsup-1
c.timer
          call xdslt2 ( t1, w1, t3, w3 )
          f2 = fctops - f1
          f1 = fctops
c         write(6,'("assembly - qicasm, time, ops, rate = ", 
c    1                l5, 3f15.3)') qicasm, t3-t2, f2,
c    2                f2 / ( t3-t2 ) / 1000000.
c.timer
c.debug
c     write(6,'("before scaling and updating")')
c.debug
c.debug
c     write(6,'("isuper, lfront, nnode, qicasm = ", 3i8, l8)') 
c    1            isuper, lfront, nnode, qicasm
c.debug

          if ( qicasm ) then 

c             --------------------------------
c             ... entire front is held in-core
c             --------------------------------

              k = fnzlf 
              if ( .not. qincor ) k = 0

c.debug
c     write(6,'("isuper, lfront, nnode, qicasm = ", 3i8, l8)') 
c    1            isuper, lfront, nnode, qicasm
c     write(6,'("pnlsiz, updsiz                = ", 3i8    )') 
c    1            pnlsiz, updsiz
c     call xdslp5('full front before factor', lfront*(lfront+1)/2, 
c    1            rwork(stfrnt+stkoff), 6 )
c     write(6,'("pre -xdslei - isuper, fnzlf , iopos1 = ", 3i8)')
c    1                          isuper, fnzlf , iopos1        
c.debug
              lnrwrk = max(1,lnzbgn+k)
              lenxtr = lrwork - lnrwrk
              call xdslei ( pvttol, cmajor, lfront, nnode, pnlsiz,
     1                      updsiz, pospon, nodoff, zpcntl, npcntl,
     2                      sup(bgxsup),    rwork (stfrnt+stkoff), 
     3                      pvtblk(bgxsup), rwork(panel+stkoff), 
     4                      rwork(swap+stkoff), rwork(temp1+stkoff),    
     4                      fnzlf,  npanel, 
     5                      xpanel, mpanel, qincor, diag(bgxsup),   
     6                      rwork(lnrwrk),        iopos1, wafil1, 
     7                      ocbufr, lbuffr, inrtia, fctops, 
     8                      slvops, ppfmon, xpboxs, psboxs,
     9                      inpexp, inpsiz, inzsiz, izfail, rtpexp,
     A                      rtpsiz, rtzsiz, rzfail, lenxtr,   ierr )

              if ( ppfmon .and. ierr .lt. 0 ) then
                 call xislfp ( xpboxs, psboxs, inpexp, inpsiz, inzsiz,
     1                         izfail, rtpexp, rtpsiz, rtzsiz, rzfail,
     2                         output )
              endif

c.debug
c     write(6,'("post-xdslei - isuper, fnzlf , iopos1 = ", 4i8)')
c    1                          isuper, fnzlf , iopos1        
c     write(6,'("post-xdslei - ierr                   = ", 3i8)')
c    1                          ierr
c     call xdslp5('full front after factor', lfront*(lfront+1)/2, 
c    1            rwork(stfrnt+stkoff), 6 )
c.debug
              if ( qzpvt1 .and. zpcntl(3) .ne. 0 ) then
                  qzpvt1    = .false.
                  k         = zpcntl(3) + nodoff
                  zpcntl(3) = perm(k)
              end if

              if ( qnpvt1 .and. npcntl(3) .ne. 0 ) then
                  qnpvt1    = .false.
                  k         = npcntl(3) + nodoff
c.debug
c     write(6,'("neg. pivot control - npcntl(3), nodbgn = ",f8.0,i8)')
c    1                                 npcntl(3), nodbgn 
c     write(6,'("neg. pivot control - k, perm(k)        = ", 2i8)')
c    1                                 k, perm(k) 
c.debug
                  npcntl(3) = perm(k)
              end if
 
              if  ( ierr .ne. 0 )  then
                 if  ( ierr .eq. -1 )  then
                    inrtia(3) = 1
                    error = -2
                 else
     1           if  ( ierr .eq. -2 )  then
                    error = -5
                 else
     1           if  ( ierr .eq. -3 )  then
                    error = -15
                 else
     1           if  ( ierr .eq. -4 )  then
                    error = -14
                 else
     1           if  ( ierr .eq. -5 )  then
                    error = -11
                 else
     1           if  ( ierr .eq. -6 )  then
                    needst = 2*lenxtr
                    error = -7 
                 else
                    error = -99
                 end if
                 return
              end if
 
          else

c             --------------------------------------------------------
c             ... front is too big to be eliminated in-core.  
c                 to reduce amount of i/o transferred in out-of-core
c                 elimination, increase the panel size if possible.
c                 if the panel size is already too large for the
c                 space available, reduce both the panel size and
c                 the update panel size together. (Note that
c                 pivoting may have increased the memory requirements
c                 beyond what was predicted, necessitating a smaller
c                 panel size.)
c                 
c                 space available should be the entire stack
c                 minus space for the temporary low rank modification
c                 block and for the pivoting permutation.
c                 (assert stfree == 1 here)
c             -------------------------------------------------------

             lavail = ( lstack - stfree )
     1                 - ( updsiz * lfront + xdslni ( lfront ) )

             if ( lavail .ge. pnlsiz * lfront ) then

c               ------------------------------------------
c               ... enlarge panel size as much as possible
c               ------------------------------------------
                
                pnlsiz = min ( nnode, lavail / lfront )

             else

c               ---------------------------------------------
c               ... shrink both panel size and update size,
c                   ensuring that update size never decreases
c                   below 2.
c               ---------------------------------------------

                lavail = lstack - stfree - xdslni (lfront)
                lnreqd = (pnlsiz + updsiz) * lfront 

                ratio  = dble ( lavail ) / dble ( lnreqd )
                updsiz = ratio * updsiz

c.debug
c      write(6,'("oce-lavail, lnreqd,  ratio = ", 2i8,f12.3)')
c    1                     lavail, lnreqd, ratio
c.debug
                if ( updsiz .ge. 2 ) then

                   pnlsiz = ratio * pnlsiz

                else

                   updsiz = 2
                   lavail =  ( lstack - stfree )
     1                       - ( 2 * lfront + xdslni ( lfront ) )
                   pnlsiz = min ( nnode, lavail / lfront )

                end if

                pnlsiz = max ( 1, pnlsiz )

             end if

c.debug
c            write(6,'("oce-pnlsiz, updsiz = ", 2i8,f12.3)')
c    1                 pnlsiz, updsiz
c.debug


c             -------------------------
c             ... now partition storage
c             -------------------------

              panel  = stfree 
              swap   = panel  + pnlsiz * lfront
              temp1  = swap   + xdslni ( lfront )
              stkuse = temp1  + updsiz * lfront
c.debug
c     write(6,'("oce-panel, temp1 , stkuse, lfront, pnlsiz= ", 5i8)')
c    1                panel, temp1 , stkuse, lfront, pnlsiz
c     write(6,'("oce-stkoff, stkbgn, updsiz, lstack       = ", 5i8)')
c    1                stkoff, stkbgn, updsiz, lstack
c.debug

              mxstck = max ( mxstck, stkuse-1 )
 
c             ----------------------------
c             ... check for stack overflow
c             ----------------------------
 
c.debug
c     write(6,'("before out-of-core elimination - lstack, stkuse = ", 
c    1         2i8)') lstack, stkuse
c.debug 
              if ( lstack .lt. stkuse-1 ) then
                  error = -1
                  return
              endif

c             ------------------------------
c             ... perform actual elimination
c             ------------------------------

c.debug
c     write(6,*) 'pre -xdsleo - isuper, fnzlf , iopos1 = ',
c    1                          isuper, fnzlf , iopos1        
c     write(6,'("pre -xdsleo - iopos2 = ", 3i8)') iopos2                        
c     call xislp3 ( 'sup pre xdsleo', bgxsup+nnode-1, sup, 6 )
c     call xislp3 ( 'lindx2 pre xdsleo', bgdxg2-1, iwork, 6 )
c.debug
              iopsv2 = iopos2
              lextra = lstack - ( stkuse - 1 )
c.debug
c     write(6,'("pre -xdsleo - lbuffr, lextra         = ", 3i8)')
c    1                          lbuffr, lextra                
c.debug

              if ( lbuffr .ge. lextra ) then

              call xdsleo ( pvttol, cmajor, lfront, nnode, pnlsiz,
     1                      updsiz, pospon, nodoff, zpcntl, npcntl,
     2                      sup(bgxsup),            pvtblk(bgxsup),         
     3                      rwork(panel+stkoff),    rwork(swap+stkoff),
     4                      rwork(temp1+stkoff),    fnzlf,  npanel, 
     5                      xpanel, mpanel, diag(bgxsup),   
     5                      wafil5, watrn5,
     6                      iopos1, wafil1, watrn1, iopos2, wafil2,  
     7                      walen2, watrn2, ocbufr, 
     8                      lbuffr, inrtia, fctops, slvops,
     8                      ppfmon, xpboxs, psboxs,
     9                      inpexp, inpsiz, inzsiz, izfail, rtpexp,
     A                      rtpsiz, rtzsiz, rzfail, ierr )

              else

c             ---------------------------------------------------
c             ... the amount of stack left over after the panel
c                 allocation exceeds already allocated storage
c                 for the out-of-core buffer.  use the remainder
c                 of the stack of this buffer for this call only.
c             ---------------------------------------------------

c.debug
c     write(6,'("before call to xdsleo using extra")')
c     write(6,'("lbuffr, lextra         = ", 2i8)')
c    1            lbuffr, lextra
c.debug

              call xdsleo ( pvttol, cmajor, lfront, nnode, pnlsiz,
     1                      updsiz, pospon, nodoff, zpcntl, npcntl,
     2                      sup(bgxsup),            pvtblk(bgxsup),         
     3                      rwork(panel+stkoff),    rwork(swap+stkoff),
     4                      rwork(temp1+stkoff),    fnzlf,  npanel, 
     5                      xpanel, mpanel, diag(bgxsup),   
     5                      wafil5, watrn5,
     6                      iopos1, wafil1, watrn1, iopos2, wafil2,  
     7                      walen2, watrn2, rwork(stkuse+stkoff), 
     8                      lextra, inrtia, fctops, slvops,
     8                      ppfmon, xpboxs, psboxs,
     9                      inpexp, inpsiz, inzsiz, izfail, rtpexp,
     A                      rtpsiz, rtzsiz, rzfail, ierr )

              end if

              if ( ppfmon .and. ierr .lt. 0 ) then
                 call xislfp ( xpboxs, psboxs, inpexp, inpsiz, inzsiz,
     1                         izfail, rtpexp, rtpsiz, rtzsiz, rzfail,
     2                         output )
              endif

c.debug
c     write(6,'("post-xdsleo - isuper, fnzlf , iopos1 = ", 3i8)')
c    1                          isuper, fnzlf , iopos1        
c     call xislp3 ( 'sup post xdsleo', bgxsup+nnode-1, sup, 6 )
c     call xislp3 ( 'lindx2 post xdsleo', bgdxg2-1, iwork, 6 )
c.debug
              if ( qzpvt1 .and. zpcntl(3) .ne. 0 ) then
                  qzpvt1    = .false.
                  k         = zpcntl(3) + nodoff
                  zpcntl(3) = perm(k)
              end if

              if ( qnpvt1 .and. npcntl(3) .ne. 0 ) then
                  qnpvt1    = .false.
                  k         = npcntl(3) + nodoff
                  npcntl(3) = perm(k)
              end if
 
              if  ( ierr .ne. 0 )  then
                 if  ( ierr .eq. -1 )  then
                    inrtia(3) = 1
                    error = -2
                 else
     1           if  ( ierr .eq. -2 )  then
                    error = -5
                 else
     1           if  ( ierr .eq. -3 )  then
                    error = -6
                 else
     1           if  ( ierr .eq. -4 )  then
                    error = -13
                 else
     1           if  ( ierr .eq. -5 )  then
                    error = -15
                 else
     1           if  ( ierr .eq. -6 )  then
                    error = -14
                 else
     1           if  ( ierr .eq. -7 )  then
                    error = -11
                 else
                    error = -99
                 end if
                 return
              end if

          end if

          stkuse = panel - 1
 
c  =====================================================================
 
c         -------------------------------------------------------------
c         record number postponed as well as the global column indices
c         and update the pvtblk vector to indicate the index into lnz
c         where the off diagonal element of the 2x2 block is stored
c         -------------------------------------------------------------
c.debug
c     write(6,'("before recording of postponed columns")')
c     write(6,'("nnode, hldnnd, isuper, bgdxg2, bgxsup = ", 5i8)')
c    1            nnode, hldnnd, isuper, bgdxg2, bgxsup 
c.debug
c.timer
          call xdslt2 ( t1, w1, t4, w4 )
          f3 = fctops - f1
c.timer
 
          ncolpp          = hldnnd - nnode
c.debug
c         if ( ncolpp .ne. 0 ) then
c             error = -2
c             inrtia(3) = 1
c             return
c         endif
c.debug
          lsuper          = xlindx (isuper+1) - xlindx (isuper)
          lndxpt          = xlindx (isuper)
          xlindx (isuper) = bgdxg2
          bgxsup          = bgxsup + nnode

c.debug
c     write(6,'("before recording of postponed columns - ncolpp = ",
c    1        i8)') ncolpp
c.debug
 
          if ( pvttol .eq. 0. ) then
 
              bgdxg2 = bgdxg2 + lsuper
 
          else
 
              if ( bgdxg2 + ncolpp + lsuper - 1 .gt. nsninx ) then
 
c                 --------------------------------------------------
c                 ... collision between new global indices and stack
c                     we need lfront entries tacked onto lindg2 and
c                     we have just freed 2*lfront entries after the
c                     partial cholesky step.  so we just have to
c                     adjust the stack here.  we increase the space
c                     for lindg2 to the min of 1.2*nsninx and what
c                     is available.
c                 --------------------------------------------------
c.debug
c     write(6,'("adjusting stack due to collision")')
c.debug
 
                  temp1  = stkbgn
                  k      = max ( 1.2 * nsninx, bgdxg2+ncolpp+lsuper-1. )
                  l      = nsninx + xdslil ( lstack - stkuse + 1 )
                  nsninx = min ( k, l )
 
                  stkbgn = xdslni ( nsninx ) + 1
                  stkoff = stkbgn - 1
                  lstack = lstack - ( stkbgn - temp1 )
c.debug
c     write(6,'("temp1 , stkbgn, stkuse, lstack = ", 4i8)')
c    1            temp1 , stkbgn, stkuse, lstack 
c.debug
 
                  call xdslmv ( stkuse, rwork, temp1, stkbgn )
 
              end if
 
             iwork(bgdxg2:bgdxg2+ncolpp-1) = sup(bgxsup:bgxsup+ncolpp-1)
              bgdxg2 = bgdxg2 + ncolpp
 
          iwork(bgdxg2:bgdxg2+lsuper-1) = lndxg1(lndxpt:lndxpt+lsuper-1)
              bgdxg2 = bgdxg2 + lsuper
 
          endif
 
c  =====================================================================

c         -------------------------------------------------------
c         ... stack information about the resulting update matrix
c         -------------------------------------------------------

c.debug
c     write(6,'("before stacking matrices")')
c     write(6,'("stacking for isuper = ", i8)') isuper
c     write(6,'("stfree, lsuper, lndxpt, ncolpp = ", 4i8)')
c    1            stfree, lsuper, lndxpt, ncolpp 
c     write(6,'("stkoff                         = ", 4i8)')
c    1            stkoff                         
c.debug

          nstack         = nstack + 1
 
          if ( nstack .gt. mstack ) then
              error = -4
              return
          end if
c.debug
c     write(6,'("nstack                         = ", 4i8)')
c    1            nstack                         
c.debug

          stknod (nstack  ) = isuper
          istack (1,nstack) = lsuper
          istack (2,nstack) = lndxpt
          istack (3,nstack) = ncolpp

          if ( qicasm ) then 
 
c             ------------------------------------------
c             stack the matrix and its indices,
c             moving it up to replace any stack matrices
c             and index lists which were assembled
c             ------------------------------------------
 
              k      = lfront*nnode - (nnode*(nnode-1)/2)
              length = matsiz - k
              start  = stfrnt + matsiz - length
c.debug
c     write(6,'("stacking - k, matsiz, length, start, nstack = ",
c    1            5i8)')     k, matsiz, length, start, nstack 
c.debug
 
              if ( length .eq. 0 ) then
                  if ( nassmb(isuper) .gt. 0 ) stfree = nwfree
                  go to 1000
              end if
 
c.debug
c     write(6,'("stacking - k, matsiz, length, start, nstack = ",
c    1            5i8)')     k, matsiz, length, start, nstack 
c     write(6,'("stacking - ncolpp, stkoff, isuper           = ",
c    1            5i8)')     ncolpp, stkoff, isuper           
c     if ( ncolpp .gt. 0 ) then
c     call xdslp5 ( 'matrix to stack', length, rwork(start+stkoff), 6 )
c     end if
c.debug
 
              if ( nassmb(isuper) .gt. 0 ) then
                  stfree         = nwfree
              endif
c.debug
c     write(6,'("stacking - stfree                           = ",
c    1            5i8)')     stfree                           
c.debug
 
              istack (4,nstack) = 0
              rstbeg (nstack  ) = stfree
 
c             ------------------------------------------
c             stack the postponed global column indices
c             ------------------------------------------
 
              call icopy (ncolpp, iwork (xlindx (isuper)), 1,
     1                            rwork (stfree+stkoff), 1 )

              stfree = stfree + ncolpp
 
c             -----------------
c             stack the update
c             -----------------

c.debug
c     if ( ncolpp .gt. 0 ) then
c     call xdslp5 ( 'matrix before stacking', length, 
c    1              rwork(start +stkoff), 6 )
c     end if
c.debug
              call xdslmv ( length, rwork, start +stkoff, 
     1                                     stfree+stkoff )
c.debug
c     if ( ncolpp .gt. 0 ) then
c     call xdslp5 ( 'matrix after stacking', length, 
c    1              rwork(stfree+stkoff), 6 )
c     end if
c.debug

              stfree = stfree + length

          else

c             ------------------------------------------
c             ... out-of-core processing.  update matrix
c                 has already been moved to wafil2.
c             ------------------------------------------

              rstbeg ( nstack    ) = 0
              istack ( 4, nstack ) = iopsv2
              stfree               = 1
c.debug
c     write(6,'("after out-of-core stacking - iopsv2 = ",i8)') iopsv2
c.debug

          end if

c.debug
c     write(6,'("after  stacking matrices")')
c     write(6,'("nstack                         = ", 4i8)')
c    1            nstack                         
c.debug
 
c  =====================================================================
c.debug
c     write(6,'("before end of 1000 loop")')
c.debug
c.timer
          call xdslt2 ( t1, w1, t5, w5 )
          times(1) = times(1) + t2
          if ( qicasm ) then
              times(2) = times(2) + t3 - t2
              times(4) = times(4) + t4 - t3 
              times(6) = times(6) + t5 - t4
              fops(1)  = fops(1) + f2
              fops(3)  = fops(3) + f3
          else
              times(3) = times(3) + t3 - t2
              times(5) = times(5) + t4 - t3
              times(7) = times(7) + t5 - t4
              fops(2)  = fops(2) + f2
              fops(4)  = fops(4) + f3
          end if
c.timer
 
 1000 continue

      if ( ppfmon ) then
         call xislfp ( xpboxs, psboxs, inpexp, inpsiz, inzsiz,
     1                 izfail, rtpexp, rtpsiz, rtzsiz, rzfail,
     2                 output )
      endif

c.timer
      t1 = times(1) + times(2) + times(3) + times(4)
     1              + times(5) + times(6) + times(7)
      t2 = times(2) + times(3)
      t3 = times(4) + times(5)
      t4 = times(6) + times(7)

c     write(6,68000) times(1), 100. * times(1) / t1
68000 format ( /5x, 'time breakdown for numerical factorization'
     1        /44x, 'time', 5x, '%', 9x, 'time', 5x, '%'
     1       //10x, 'initialization phase      ', 19x, f12.3, 1x, f5.1)

c     write(6,68001) times(2), 100. * times(2) / t2, 
c    1               fops(1)/times(2)/1000000.,
c    2               times(3), 100. * times(3) / t2, 
c    3               fops(2)/times(3)/1000000.,
c    4               t2, 100. * t2 / t1, 
c    5               ( fops(1) + fops(2) ) / t2 / 1000000.
68001 format (/10x, 'assembly       phase i.c. ', f12.3, 1x, f5.1  
     1        /10x, '               rate       ', f12.3
     2        /10x, 'assembly       phase o.c. ', f12.3, 1x, f5.1 
     3        /10x, '               rate       ', f12.3
     4        /10x, 'assembly       phase tot. ', 19x, f12.3, 1x, f5.1
     5        /10x, '               rate       ', 19x, f12.3 )

c     write(6,68002) times(4), 100. * times(4) / t3,
c    1               fops(3)/times(4)/1000000.,
c    2               times(5), 100. * times(5) / t3,
c    3               fops(4)/times(5)/1000000.,
c    4               t3, 100. * t3 / t1, 
c    5               ( fops(3) + fops(4) ) / t3 / 1000000.
68002 format (/10x, 'elimination    phase i.c. ', f12.3, 1x, f5.1  
     1        /10x, '               rate       ', f12.3
     2        /10x, 'elimination    phase o.c. ', f12.3, 1x, f5.1 
     3        /10x, '               rate       ', f12.3
     4        /10x, 'elimination    phase tot. ', 19x, f12.3, 1x, f5.1
     5        /10x, '               rate       ', 19x, f12.3 )

c     write(6,68003) times(6), 100. * times(6) / t4,
c    1               times(7), 100. * times(7) / t4,
c    2               t4, 100. * t4 / t1
68003 format (/10x, 'stacking       phase i.c. ', f12.3, 1x, f5.1  
     1        /10x, 'stacking       phase o.c. ', f12.3, 1x, f5.1 
     2        /10x, 'stacking       phase tot. ', 19x, f12.3, 1x, f5.1)

c     write(6,68004) t1, fctops / t1 / 1000000.
68004 format (/10x, 'total time in xdslf2       ', 18x, f12.3 
     1        /10x, 'overall computational rate ', 18x, f12.3 )
c.timer
 
      error   = 0
      nsnind = bgdxg2 - 1
      lstack = mxstck
 
      xsup  (nsuper+1) = npanel+1
      xpanel(npanel+1) = bgxsup
      xlindx(nsuper+1) = bgdxg2
c.debug
c     write(6,'("leaving xdslf2 - bgxsup = ", i8)') bgxsup
c     call xislp3 ( 'sup', bgxsup-1, sup, 6 )
c.debug
 
      if ( bgxsup - 1 .ne. neqns ) then
          error     = -2
          inrtia(3) = 1
      end if
 
      if ( .not. qincor ) then
          walen1 = fnzlf
          watrn1 = watrn1 + fnzlf
      end if

c  =====================================================================
 
      return
      end
