      subroutine xdslfa ( work,   lwork,  cndnum, inrtia, needs,
     &                    error )
 
c
c     purpose
c     -------
c
c     xdslfa is the top level driver for numeric factorization phase.
c
c     created         30-jan-89   -- rgg --
c     modified        09-jun-89   -- rgg -- modified to allow pivoting
c                                           and condition no. estimation
c     modified        29-oct-90   -- mlc -- calls to xislvr protected
c                                           against nzla .le. 0
c     modified        08-feb-91   -- rgg -- mods to allow no i/o to
c                                           sqfil1
c     modified        12-feb-91   -- rgg -- mods to allow no i/o to
c                                           wafil1 and wafil2
c     modified        03-dec-92   -- rgg -- mods to allow minimum core
c                                           processing
c     modified        19-may-93   -- rgg -- mods to not generate lindxg2
c                                           when pvttol = 0.
c     modified        16-jan-97   -- rgg -- mods to work with release 4
c     modified        19-may-98   -- rgg -- mods to work with release 4
c                                           unsymmetric matrices
c     modified        31-jan-01   -- jgl -- added new error codes to
c                                           reflect underlying code mods
c     modified        02-apr-02   -- dkw -- fixed interr code -4 to
c                                           report -502 instead of -503
c     modified        20-mar-03   -- dkw -- set stage back to 40 instead
c                                           of -1 for error = -503
c
c     input arguments
c     ---------------
c
c     lwork       i   length of work array.
c
c     input/output arguments
c     ----------------------
c
c     work        d   work array.  on input it contains the
c                    .CMNication area and all active arrays.
c
c     output arguments
c     ----------------
c
c     cndnum      d   if user requested condition number estimate to be
c                     computed via xdslsp, then this contains on output 
c                     the 2 norm condition number estimate.
c                     else set to 0. for return
c     inrtia      i   integer array of length 3 holding the number
c                     of positive, negative, and zero eigenvalues in
c                     the three components, respectively.
c     needs       i   amount of workspace required for the next stage.
c     error       i   error flag
c                     =    0  normal return
c                     = -500  incorrect processing path.
c                     = -501  lwork not large enough.
c                     = -502  storage error in numeric factorization.
c                     = -503  numerical failure -- numerically singular
c                             matrix (when using pivoted
c                             factorization)
c                     = -504  i/o error on sqfil1 or sqfil2
c                     = -505  i/o error on wafil1
c                     = -506  i/o error on wafil2
c                     = -507  ran out of memory during factorization to
c                             lnz and unz with wafil1 deactivated
c                     = -508  ran out of memory during factorization to
c                             spill stack with wafil2 deactivated
c                     = -509  fill limit exceeded.
c                     = -510  numerical failure -- negative pivot found
c                             with negative pivot controls activated
c                             (non-pivoting mode)
c                     = -511  i/o error on wafil4
c                     = -512  i/o error on wafil5
c                     = -513  numerical failure -- pivoting failure,
c                             unable to find stable pivot within panel
c                             (symmetric pivoting mode only)
c                     = -514  numerical failure -- pivoting failure,
c                             exact zero pivot in non-pivoting mode
c                     = -599  failure in lower level routine with
c                             undocumented error return code
c
c---------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      integer             lwork,  inrtia(3), needs, error
 
      double precision    work(*),   cndnum
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------
 
      integer             anzp,   arindx, diag,   imstck, inuse,
     2                    invp,   invsup, itemp,  lindxg,
     3                    lindxl, lindx2, llnz,   lnz,    lnzbgn,
     4                    lstack, msglvl, mstack, mfront, mxused,
     5                    ineqns, ineqn1, insnin, insup1, interr,
     6                    istack, front,  mxlfrt, mpanel, k,
     7                    minstk, mxtotf, insup , inzla,  mxanzf,
     8                    l     , lbuffr, ltemp,  fpnlsz, cmajor,
     9                    iofile, fdiag , bmxtyp, fupdsz, needst
 
      integer             nassmb, neqns,  nsnind, nsninx, nsnin2,
     1                    nsuper, nzla,   
     2                    ocbufr, output, perm,   pvtblk, rstbeg,
     3                    sqfil1, stage,  stkbgn, stknod, stkstr,
     4                    slvspc, sqfil2, sup,    temp,   
     5                    wafil1, wafil2, wdynam, wkreqd, wleft,
     6                    xarndx, xlindx, xlind2, xsup,   xsup2,
     7                    xpanel, npanel, amncor, mxvlrb, mxvlib,
     8                    mincor, sqfil3, slvbsz, wafil4, wafil5
 
      integer             idummy(1)
 
      logical             qincor, unsym , ppfmon
 
      double precision    fctops, pvttol, slvops, t1,     t2,
     1                    tnorma, walen1, walen2, watrn1, watrn2, 
     2                    w1    , w2,     xtra,   sqtrn1, sqtrn3, 
     3                    extfil, fnzlf , fnzlf2


      double precision    walen4, walen5, watrn4, watrn5

c.debug
c     double precision    t3, t4, w3, w4  
c.debug
 
      double precision    sigma,  dummy(1)
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      integer             xdslil, xdslni
 
      external            xdslf2, xdslni, xislp3, xdslp5, 
     1                    xdslt1, xdslt2, xdslc2, xdslc3, xislp1,
     2                    xdslw4, xdslil, xdslmv, xislmv
c
c---------------------------------------------------------------------
 
c     ---------------------------------------
c     ... start preparation for factorization
c     ---------------------------------------
 
      error  = 0
      needs  = 0
 
      call xdslt1 ( t1, w1 )
 
c     ----------------------------------------------------
c     ... extract information from the.CMNication area.
c     ----------------------------------------------------
 
      stage  = work ( qstage )
      msglvl = work ( qmsglv )
      output = work ( qoutpu )
      unsym  = work ( qmxtyp ) .eq. 2.
 
      if ( stage .ne. 40 ) go to 8500

      bmxtyp = work ( qbmxty )
 
      neqns  = work ( qneqns )
      nzla   = work ( qnzla  )

      mincor = work ( qmncor )
      amncor = mod ( mincor, 10 )
 
      xarndx = work ( qxadj  )
      perm   = work ( qperm  )
      invp   = work ( qinvp  )
      diag   = work ( qadiag )
 
      xsup   = work ( qxsup  )
      nassmb = work ( qnassm )
      xlindx = work ( qxlind )
      lindxg = work ( qlndxg )
      lindxl = work ( qlndxl )

      if ( amncor .eq. 0 ) then
          anzp   = work ( qofdia )
          arindx = work ( qarndx )
      else
          anzp   = 1
          arindx = 1
      end if

      inuse  = work ( qinuse )
c.debug
c     write(6,'("inuse at start of xdslfa               = ", 5i8)')
c    1            inuse                                 
c.debug
 
c     -----------------------------------------------------
c     ... extract an estimate of the two norm of the matrix
c         computed during value input phase for use by the 
c         condition number estimation
c     -----------------------------------------------------
 
      tnorma = work ( qcndnm )
 
c---------------------------------------------------------------------
 
c     ------------------------------------------------
c     ... start preparation of numerical factorization
c     ------------------------------------------------
 
      if ( msglvl .ge. 2 ) write ( output, 60000 )
 
      inrtia(1) = -1
      inrtia(2) = -1
      inrtia(3) = -1
 
      pvttol = work ( qpvttl )
     
      if ( pvttol .lt. 0. .or. pvttol .gt. 0.5 ) pvttol = .01

      fpnlsz = work ( qfpnls )
      fupdsz = work ( qfupds )
      cmajor = work ( qcmajr )

      fnzlf  = work ( qnzlf  )
      extfil = 0.
      if ( work(qextfl) .gt. 0. ) extfil = work(qextfl) * fnzlf

c     ---------------------
c     ... prepare i/o files
c     ---------------------
 
      sqfil1 = work ( qsqfl1 )
      sqfil2 = work ( qsqfl2 )
      sqfil3 = 0
 
      wafil1 = work ( qwafl1 )
      wafil2 = work ( qwafl2 )
      wafil4 = work ( qwafl4 )
      wafil5 = work ( qwafl5 )
 
      if ( sqfil1 .gt. 0 ) then
          iofile = sqfil1
          call xislrw ( sqfil1, interr )
          if  ( interr .ne. 0 )  then
             if  ( interr .eq. -1  .or.  interr .eq. -2 )  then
                go to 8504
             else
                go to 8599
             end if
          end if
      end if

      nsnind = work ( qnsnin )

      if ( sqfil2 .gt. 0 .and. nsnind .gt. 0 ) then
          iofile = sqfil2
          call xislrw ( sqfil2, interr )
          if  ( interr .ne. 0 ) then
             if  ( interr .eq. -1  .or.  interr .eq. -2 )  then
                go to 8504
             else
                go to 8599
             end if
          end if
      end if
 
c     ----------------------------------------------------
c     ... extract information from the.CMNication area.
c     ----------------------------------------------------
 
      nsuper = work ( qnsupe )
      mstack = work ( qmstac )
      stkstr = work ( qstkst )
      mfront = work ( qmfron )
      mxanzf = work ( qmxanz )
      ppfmon = work ( qppfmn ) .ne. 0.0
 
c     -------------------------------------------------------
c     ... compute maximum value for the order of a front from
c         maximum storage required for a front.
c         then compute stack storage requirements.
c     -------------------------------------------------------
 
      mxlfrt = ( -1. + sqrt ( 1. + 8. * mfront ) ) / 2.
c.debug
c     write(6,'("mfront, mxlfrt                         = ", 5i8)')
c    1            mfront, mxlfrt                        
c.debug

      stkstr = stkstr + ( fpnlsz + fupdsz ) * mxlfrt + xdslni ( mxlfrt )
      minstk =          ( fpnlsz + fupdsz ) * mxlfrt + xdslni ( mxlfrt )
 
      if ( unsym ) then
          stkstr = 2 * stkstr
          minstk = 2 * minstk
      end if

c     --------------------------------------
c     ... bring in local indices from sqfil2
c     --------------------------------------

      insnin = xdslni ( nsnind )

      if ( sqfil2 .gt. 0 ) then

          iofile = sqfil2

          if ( pvttol .eq. 0. ) then

              lindxl = lindxg

              call xislvr ( sqfil2, nsnind, work(lindxl), interr )
              if  ( interr .ne. 0 ) then
                 if  ( interr .eq. -2  .or.  interr .eq. -3 )  then
                    go to 8504
                 else
                    go to 8599
                 end if
              end if

              work(qsqtr2) = work(qsqtr2) + insnin
              
          else

              lindxl = inuse + 1
              inuse  = lindxl + insnin - 1

              call xislvr ( sqfil2, nsnind, work(lindxl), interr )
              if  ( interr .ne. 0 )  then
                 if  ( interr .eq. -2  .or.  interr .eq. -3 )  then
                    go to 8504
                 else
                    go to 8599
                 end if
              end if

              work(qsqtr2) = work(qsqtr2) + insnin

          end if
      
      end if
c.debug
c         call xislp1 ( 'relative indices of l',
c    1                 nsuper, work(xlindx), work(lindxl), output )
c.debug

c     --------------------------------------------------
c     ... allocate storage for factor diagonal if saving
c         matrix a
c     --------------------------------------------------

      if ( work(qsavea) .ne. 0. ) then
          fdiag = inuse + 1
          inuse = inuse + neqns
      else
          fdiag = diag
      end if

      wkreqd = inuse
      if ( wkreqd .gt. lwork ) go to 8501
 
c     ----------------------------------------
c     ... optionally print out the information
c     ----------------------------------------
c.debug
c     write(6,'("xarndx, perm  , invp  , xsup  , nassmb = ", 5i8)')
c    1            xarndx, perm  , invp  , xsup  , nassmb
c     write(6,'("xlindx, lindxg, lindxl, diag  , anzp   = ", 5i8)')
c    1            xlindx, lindxg, lindxl, diag  , anzp  
c     write(6,'("arindx, inuse , lwork , amncor         = ", 5i8)')
c    1            arindx, inuse , lwork , amncor        
c     write(6,'("nsnind, insnin, fdiag , wkreqd         = ", 5i8)')
c    1            nsnind, insnin, fdiag , wkreqd        
c.debug
 
      if ( msglvl .ge. 3 ) then
          call xislp3 ( 'supernodal partition',
     1                  nsuper+1, work(xsup),   output )
          call xislp3 ( 'number of assemblies',
     1                  nsuper,   work(nassmb), output )
          call xdslp5 ( 'values of permuted diag',
     1                  neqns, work(diag), output )
      endif
      if ( msglvl .eq. 3 ) then
          call xislp3 ( 'supernodal indices pointer',
     1                  nsuper+1, work(xlindx), output )
          call xislp3 ( 'relative assembly indices pointer',
     1                  neqns+1, work(xarndx), output )
      endif
      if ( msglvl .ge. 4 ) then
          call xislp1 ( 'relative indices of l',
     1                 nsuper, work(xlindx), work(lindxl), output )
          if ( pvttol .gt. 0. ) then
              call xislp1 ( 'global indices of l',
     1                 nsuper, work(xlindx), work(lindxg), output )
          end if
          if ( amncor .eq. 0 ) then
              call xislp1 ( 'lower adjacency of papt',
     1                     neqns, work(xarndx), work(arindx), output )
              call xdslp5 ( 'values of permuted a',
     1                      nzla, work(anzp), output )
          else
              call xislp3 ( 'relative assembly indices pointer',
     1                      neqns+1, work(xarndx), output )
          end if
      endif
 
c     --------------------------------------------------------------
c     ... compute integer array lengths for storage management logic
c     --------------------------------------------------------------
 
      imstck  = xdslni ( mstack )
      ineqns  = xdslni ( neqns  )
      ineqn1  = xdslni ( neqns + 1 )
      insup   = xdslni ( nsuper )
      insup1  = xdslni ( nsuper + 1 )
      inzla   = xdslni ( nzla )
 
c     -----------------------------------------------
c     ... compute other pointers into the work array.
c     -----------------------------------------------
 
      sup     = inuse + 1
      pvtblk  = sup    + ineqns
      xsup2   = pvtblk + ineqns
      xlind2  = xsup2  + insup1
      xpanel  = xlind2 + insup1
      stknod  = xpanel + ineqn1
      rstbeg  = stknod + imstck
      istack  = rstbeg + imstck
      ocbufr  = istack + 4*imstck
      lbuffr = max ( neqns, fupdsz*mxlfrt )
c.debug
c     write(6,'("lbuffr, neqns , fupdsz, mxlfrt, amncor = ", 5i8)')
c    1            lbuffr, neqns , fupdsz, mxlfrt, amncor
c.debug
 
      if ( amncor .ne. 0 ) then
          mxvlrb = ocbufr + lbuffr
          mxvlib = mxvlrb + mxanzf
          inuse  = mxvlib + xdslni ( mxanzf ) - 1
      else
          mxvlrb = 1
          mxvlib = 1
          inuse  = ocbufr + lbuffr - 1
      end if
c.debug
c     write(6,'("xarndx, perm  , invp  , xsup  , nassmb = ", 5i8)')
c    1            xarndx, perm  , invp  , xsup  , nassmb
c     write(6,'("xlindx, lindxg, lindxl, diag  , anzp   = ", 5i8)')
c    1            xlindx, lindxg, lindxl, diag  , anzp  
c     write(6,'("arindx, sup   , pvtblk, xsup2 , xlind2 = ", 5i8)')
c    1            arindx, sup   , pvtblk, xsup2 , xlind2
c     write(6,'("xpanel, stknod, rstbeg, istack, ocbufr = ", 5i8)')
c    1            xpanel, stknod, rstbeg, istack, ocbufr
c     write(6,'("mxvlrb, mxvlib, inuse , neqns , nzla   = ", 5i8)')
c    1            mxvlrb, mxvlib, inuse , neqns , nzla
c     write(6,'("fnzlf , stkstr, minstk                 = ", 5i8)')
c    1            fnzlf , stkstr, minstk
c     write(6,'("sqfil1, sqfil2                         = ", 5i8)')
c    1            sqfil1, sqfil2          
c     write(6,'("wafil1, wafil2, wafil4, wafil5         = ", 5i8)')
c    1            wafil1, wafil2, wafil4, wafil5
c     write(6,'("in xdslfa - lbuffr, ocbufr, mxvlrb = ", 3i8)')
c    1                        lbuffr, ocbufr, mxvlrb 
c     write(6,'("            mxvlib, neqns , mxanzf = ", 3i8)')
c    1                        mxvlib, neqns , mxanzf 
c     write(6,'("            pvttol                 = ", 1pd15.5)')
c    1                        pvttol                 
c.debug

c     --------------------------------------------------
c     ... attempt to determine if storage is too big and
c         there is a 32 bit integer limitation
c     --------------------------------------------------

      k = 1
      if ( unsym ) k = 2
      l = k * fnzlf

      if ( l .ne. k * fnzlf .and. wafil1 .lt. 0 ) go to 8507
 
c     --------------------------------------------
c     ... determine the amount of in-core required
c     --------------------------------------------
 
      wkreqd = inuse
      wdynam = stkstr + k * fnzlf

      if ( pvttol .gt. 0. ) then
          wkreqd = wkreqd + insnin
          wdynam = wdynam + insnin
      end if
 
c     -------------------------------------------------------------
c     ... determine amount of working storage needed depending on
c         the enabling and disabling of wafil1, wafil2, wafil4, and
c         wafil5.  
c
c         notes:
c         1.  lower triangle factor will be on wafil1
c         2.  multifrontal stack will be on wafil2 - a scratch file
c         3.  upper triangle factor (unsym problems only) will be 
c             on wafil4
c         4.  if a front gets too big to fit into memory, wafil5
c             will be used for the processing of that front.  
c             wafil5 is also a scratch file.
c         5.  out-of-core processing of the stack can be disabled
c             by setting wafil2 to 0.
c         6.  out-of-core storage of the factors can be disabled
c             by setting wafil1 to 0.  If wafil1 is nonzero then
c             wafil4 must be nonzero for unsym problems.
c     -------------------------------------------------------------
 
c         ----------------------
c         ... wafil1 is disabled
c         ----------------------
 
      if ( wafil1 .lt. 0 ) then
 
          if ( unsym ) then
              wkreqd = wkreqd + 2*fnzlf
          else
              wkreqd = wkreqd + fnzlf
          end if
 
      end if
 
c     ----------------------
c     ... wafil2 is disabled
c     ----------------------
 
      if ( wafil2 .lt. 0 ) then
 
          wkreqd = wkreqd + stkstr
 
      else
 
          wkreqd = wkreqd + minstk
 
      end if
c.debug
c     write(6,'("before memory length check no. 2")')
c.debug
 
      if ( wkreqd .gt. lwork .or. wkreqd .lt. 0 ) go to 8501
 
c     --------------------------------------------------
c     ... determine the extra slack allowed for pivoting
c     --------------------------------------------------
 
      wleft   = int(0.95d0*lwork) - ( inuse + wdynam )
      xtra    = dble ( wleft ) / dble ( wdynam )

      if ( wdynam .lt. 0 ) xtra = -1.
 
c     --------------------------------------------------------------
c     ... the remainder of the workspace is divided between lndxg2,
c         the real stack, lnz and unz.  lndxg2 always starts at 1.
c         if lndxg2 grows beyond current allocation the stack is
c         adjusted by xdslf2.
c     --------------------------------------------------------------
c.debug
c     write(6,'("stkstr, minstk, amncor, xtra = ", 3i8, f12.2)')
c    1            stkstr, minstk, amncor, xtra 
c.debug
 
      if ( ( xtra .gt. 2. * pvttol .or. wafil1 .lt. 0 ) .and.
     1     amncor .eq. 0 ) then
 
c         --------------------------------------------------
c         ... there is enough memory to attempt to keep lnz
c             (and unz) in core
c         -------------------------------------------------
 
          nsninx  = ( 1. + xtra ) * nsnind
          if ( pvttol .eq. 0. ) nsninx = 0
 
          insnin  = xdslni ( nsninx )
          lstack  = ( 1. + xtra ) * stkstr
 
c         ------------------------------------------------
c         ... in the case where wleft .lt. 0 which can only
c             happen when wafil2 is enabled then force the
c             stack to be out-of-core.
c         ------------------------------------------------
 
          if ( wleft .lt. 0 .and. wafil2 .ge. 0 ) then
 
              wdynam  = wdynam - stkstr + minstk
              wleft   = lwork - ( inuse + wdynam )
              xtra    = dble ( wleft ) / dble ( wdynam )
 
              nsninx  = ( 1. + xtra ) * nsnind
              if ( pvttol .eq. 0. ) nsninx = 0
 
              insnin  = xdslni ( nsninx )
              lstack  = ( 1. + xtra ) * minstk
 
          end if
 
          lindx2  = inuse + 1
          stkbgn  = insnin + 1
          lnzbgn  = stkbgn + lstack
          llnz    = (lwork - inuse) - (insnin + lstack)
 
          qincor  = .true.
 
      else
 
c         --------------------------------------
c         ... lnz ( and unz ) starts out-of-core
c         --------------------------------------
 
          wdynam = stkstr
          if ( pvttol .gt. 0. ) wdynam = wdynam + insnin
 
          wleft   = lwork - ( inuse + wdynam )
          xtra    = dble ( wleft ) / dble ( wdynam )
 
c         ------------------------------------------------
c         ... in the case where wleft .lt. 0 which can only
c             happen when wafil2 is enabled then force the
c             stack to be out-of-core.
c         ------------------------------------------------
 
          if ( wleft .lt. 0 .and. wafil2 .ge. 0 ) then
 
              wdynam  = minstk
              if ( pvttol .gt. 0. ) wdynam  = insnin + minstk
              wleft   = lwork - ( inuse + wdynam )
              xtra    = dble ( wleft ) / dble ( wdynam )
 
          end if
 
          nsninx  = ( 1. + xtra ) * nsnind
          if ( pvttol .eq. 0. ) nsninx = 0
 
          insnin  = xdslni ( nsninx )
 
          lindx2  = inuse + 1
          stkbgn  = insnin + 1
          lstack  = ( lwork - inuse ) - insnin
 
          lnzbgn  = 0
          llnz    = 0
          qincor  = .false.
 
      end if
c.debug
c     write(6,'("lindx2, lstack, llnz                   = ", 3i8)')
c    1            lindx2, lstack, llnz                 
c     write(6,'("qincor                                 = ",  l8)')
c    1            qincor                               
c.debug
 
      wkreqd  = lwork
 
c     ------------------------------------------
c     ... open word addressable i/o file for lnz
c     ------------------------------------------
 
      if ( wafil1 .ge. 0 ) then
          call xdslw8 ( wafil1, interr )
          if  ( interr .ne. 0 )  then
             if  ( interr .eq. -1 )  then
                go to 8505
             else
                go to 8599
             end if
          end if
      end if
 
c     -----------------------------------------------------
c     ... open word addressable i/o file for spilling stack
c     -----------------------------------------------------
 
      if ( wafil2 .ge. 0 ) then
          call xdslw4 ( wafil2, 2, interr )
          if  ( interr .ne. 0 )  then 
             if  ( interr .eq. -1 )  then
                go to 8506
             else
                go to 8599
             end if
          end if
      end if
 
c     ------------------------------------------
c     ... open word addressable i/o file for unz
c     ------------------------------------------
 
      if ( unsym .and. wafil1 .ge. 0 ) then
          if ( wafil4 .le. 0 ) go to 8505
          call xdslw8 ( wafil4, interr )
          if  ( interr .ne. 0 )   then 
             if  ( interr .eq. -1 )  then
                go to 8511
             else
                go to 8599
             end if
          end if
      end if
 
c     -------------------------------------------------
c     ... open word addressable i/o file for processing 
c         large fronts out-of-memory
c     -------------------------------------------------
 
      if ( wafil5 .ge. 0 ) then
          call xdslw4 ( wafil5, 2, interr )
          if  ( interr .ne. 0 )    then 
             if  ( interr .eq. -1 )  then
                go to 8512
             else
                go to 8599
             end if
          end if
       end if
 
c     ----------------------------------------------
c     ... copy xsup and xlindx into xsup2 and xlind2
c     ----------------------------------------------
 
      call xislmv ( nsuper+1, work, xdslil(xsup  -1)+1, 
     1                              xdslil(xsup2 -1)+1 ) 
      call xislmv ( nsuper+1, work, xdslil(xlindx-1)+1, 
     1                              xdslil(xlind2-1)+1 )

c     ----------------------------------------------------------
c     ... if saving original matrix entries copy diag into fdiag
c     ----------------------------------------------------------
c.debug
c     write(6,'("work(qsavea), diag, fdiag = ", f10.3, 2i8)')
c    1            work(qsavea), diag, fdiag 
c     call xislp1 ( 'relative indices of l - pre-copy',
c    1             nsuper, work(xlindx), work(lindxl), output )
c.debug

      if ( diag .ne. fdiag ) then
          call xdslmv ( neqns, work, diag, fdiag )
      end if
c.debug
c     call xislp1 ( 'relative indices of l - post-copy',
c    1             nsuper, work(xlindx), work(lindxl), output )
c.debug
 
c     ----------------------------------------------
c     ... perform the numerical factorization.
c         input dummy arguments for a second matrix.
c     ----------------------------------------------
 
      interr = 0
      sigma  = 0.

c.debug
c     call xdslt2 ( t1, w1, t3, w3 ) 
c     write(6,'("before actual factorization")')
c     write(6,'("elapsed cpu time = ", f15.6)') t3
c.debug
 
      if ( unsym ) then
 
          call xdslu2 ( pvttol, amncor, fpnlsz, fupdsz, neqns,
     1                  work(fdiag),  work(xarndx), work(arindx), 
     1                  work(anzp),   
     2                  sigma, 0, dummy, idummy, idummy, dummy,
     2                  extfil, work(qzpcnt), work(qnpcnt), 
     2                  work(perm), nsuper, work(xsup2), work(sup), 
     3                  work(nassmb), npanel, work(xpanel),
     4                  work(xlind2), work(lindxl), work(lindxg),
     5                  nsninx, work(lindx2), work(lindx2), stkbgn,
     6                  lstack, lnzbgn, llnz, qincor, mstack,
     7                  work(stknod), work(rstbeg), work(istack),
     8                  sqfil1, sqfil3, wafil1, wafil2, wafil4, wafil5,
     8                  lbuffr, work(ocbufr), 
     9                  work(mxvlrb), work(mxvlib), fnzlf2,
     a                  nsnin2, work(pvtblk), walen1, walen2, 
     b                  walen4, walen5,
     b                  sqtrn1, sqtrn3, watrn1, watrn2, watrn4, watrn5,
     c                  mpanel, fctops, slvops, interr, ppfmon, output )
 
          inrtia(1) = -1
          inrtia(2) = -1
          inrtia(3) = -1

      else

          lenwk = lwork - inuse
          call xdslf2 ( pvttol, amncor, fpnlsz, fupdsz, cmajor, 
     1                  neqns, work(fdiag),
     1                  work(xarndx), work(arindx), work(anzp),
     2                  sigma, 0, dummy, idummy, idummy, dummy,
     2                  extfil, work(qzpcnt), work(qnpcnt), 
     2                  work(perm), nsuper, work(xsup2), work(sup), 
     3                  work(nassmb), npanel, work(xpanel),
     4                  work(xlind2), work(lindxl), work(lindxg),
     5                  nsninx, work(lindx2), work(lindx2), stkbgn,
     6                  lstack, lnzbgn, llnz, qincor, mstack,
     7                  work(stknod), work(rstbeg), work(istack),
     8                  sqfil1, sqfil3, wafil1, wafil2, wafil5, 
     8                  lbuffr, work(ocbufr), 
     9                  work(mxvlrb), work(mxvlib), fnzlf2,
     a                  nsnin2, inrtia, work(pvtblk), walen1, walen2, 
     b                  walen5, sqtrn1, sqtrn3, watrn1, watrn2, watrn5, 
     c                  mpanel, fctops, slvops, interr, ppfmon, output,
     d                  lenwk, needst )
 
          walen4 = 0.
          watrn4 = 0.
 
      end if

c.debug
c     write(6,'("after  actual factorization")')
c     call xdslt2 ( t1, w1, t4, w4 ) 
c     write(6,'("elapsed cpu time = ", f15.6)') t4 - t3
c     write(6,'("fnzlf, fnzlf2, interr = ", 3i8)') fnzlf, fnzlf2, interr
c     write(6,'("watrn1, watrn4   = ", 2f15.6)') watrn1, watrn4
c.debug
 
      iofile = sqfil1

      if  ( interr .ne. 0 )  then
         if       ( interr .eq.  -1 .or. interr .eq. -4 )  then
            go to 8502
         else if  ( interr .eq.  -2 )  then
            go to 8503
         else if  ( interr .eq.  -5 )  then
            go to 8505
         else if  ( interr .eq.  -6 )  then
            go to 8506
         else if  ( interr .eq.  -7 )  then
            wkreqd = lwork + needst
            go to 8507
         else if  ( interr .eq.  -8 )  then
            wkreqd = lwork + needst
            go to 8508
         else if  ( interr .eq.  -9 )  then
            go to 8504
         else if  ( interr .eq. -10 )  then
            go to 8509
         else if  ( interr .eq. -11 )  then
            go to 8510
         else if  ( interr .eq. -12 )  then
            go to 8511
         else if  ( interr .eq. -13 )  then
            go to 8512
         else if  ( interr .eq. -14 )  then
            go to 8514
         else if  ( interr .eq. -15 )  then
            go to 8513
         else
            go to 8599
         end if
      end if

      mxused = work ( qmxuse )
 
      nsninx = nsnin2
      if ( pvttol .eq. 0. ) nsninx = 0
 
      if ( amncor .eq. 0 ) then
          if ( unsym ) then
              l = inuse + xdslni(nsninx) + lstack + 2*fnzlf2 
              l = min ( lwork, l )
          else
              l = inuse + xdslni(nsninx) + lstack + fnzlf2 
              l = min ( lwork, l )
          end if
          mxused = max ( mxused, l )
      else
          l = min ( lwork, inuse + xdslni(nsninx) + lstack )
          mxused = max ( mxused, l )
      end if
c.debug
c     write(6,'("at check point 1")')
c.debug
 
c     ------------------------------------------
c     ... close wafil2 (used to spill the stack)
c     ------------------------------------------
 
      if ( wafil2 .ge. 0 ) then
          call xdslw3 ( wafil2, 2, interr )
          if ( interr .ne. 0 ) go to 8506
      end if
 
c     -----------------------------------------------
c     ... close wafil5 (used to process large fronts)
c     -----------------------------------------------
 
      if ( wafil5 .ge. 0 ) then
          call xdslw3 ( wafil5, 2, interr )
          if ( interr .ne. 0 ) go to 8512
      end if
c.debug
c     write(6,'("at check point 2")')
c.debug

c     ---------------------------------
c     ... find end of storage in xpanel
c     ---------------------------------

      inuse = xpanel + xdslni ( npanel + 1 ) - 1
 
c     ----------------------------------------------------
c     ... put lindx2, if generated, into its final place.
c         note that if pvttol = 0. and sqfil2 .gt. 0 then 
c         sqfil2 is properly positioned for read of global
c         indices.  If sqfil2 .lt. 0 then just change the
c         pointer for lindx2 to the global indicies in 
c         memory.
c     ----------------------------------------------------
 
      itemp  = lindx2
 
      if ( pvttol .eq. 0. ) then
 
          lindx2 = lindxg

          if ( sqfil2 .gt. 0 ) then 

              iofile = sqfil2

              call xislvr ( sqfil2, nsnind, work(lindx2), interr )
              if ( interr .ne. 0 ) go to 8504

              work(qsqtr2) = work(qsqtr2) + insnin

          end if
 
      else
 
          lindx2 = inuse + 1
          call xislmv ( nsnin2, work, xdslil(itemp -1)+1, 
     1                                xdslil(lindx2-1)+1 )
 
          inuse  = lindx2 + xdslni ( nsnin2 ) - 1
 
      end if
c.debug
c     write(6,'("at check point 3")')
c.debug
 
c     ------------------------------------------------------
c     ... determine whether lnz/unz is to be store in or out
c         of core.
c     ------------------------------------------------------
 
      if ( unsym ) then
          l = 2 * fnzlf2
      else
          l = fnzlf2
      end if
 
      wkreqd = inuse + l + ineqns + neqns
 
      if ( wkreqd .gt. lwork  .or.  wkreqd .lt. 0  .or. 
     1     amncor .ne. 0 ) then

          lnz = 0
 
c         -------------------------------------------------
c         ... lnz/unz will be stored on wafil1 and wafil4.
c             if qincor = .false. then it is already there.
c         -------------------------------------------------
 
          if ( qincor ) then
 
              if ( wafil1 .lt. 0 ) go to 8507
 
              itemp = itemp + lnzbgn - 1

              call xdslfb ( unsym , nsuper, npanel, 
     1                      work(xsup2), work(xlind2), work(xpanel),
     2                      work(itemp), wafil1, wafil4, interr )
c.debug
c     write(6,'("after xdslfb - interr = ", i8)')
c.debug

              if ( interr .eq. 1 ) go to 8505
              if ( interr .eq. 2 ) go to 8511
 
              walen1 = fnzlf2
              watrn1 = watrn1 + fnzlf2

              if ( unsym ) then
                  walen4 = fnzlf2
                  watrn4 = watrn4 + fnzlf2
              end if
 
              lnz = 0
              qincor = .false.
 
          end if
 
      else
 
c         --------------------------------------
c         ... lnz/unz are to be stored in memory
c         --------------------------------------
 
          if ( qincor ) then
 
c             ---------------------------------------
c             ... lnz is in core, compress workspace.
c             ---------------------------------------
 
              lnz   = inuse + 1
              inuse = inuse + l
              itemp = itemp + lnzbgn - 1
 
              call xdslmv ( l, work, itemp, lnz )
 
          else
 
c         ----------------------------------------------------
c         ... lnz/unz is out-of-core.
c             restore into memory and close wafil1 and wafil4.
c         ----------------------------------------------------
 
              lnz    = inuse + 1
              inuse  = lnz + l - 1
 
              call xdslfc ( unsym , nsuper, npanel, 
     1                      work(xsup2), work(xlind2), work(xpanel),
     2                      wafil1, wafil4, work(lnz), interr )
c.debug
c     write(6,'("after xdslfc - interr = ", i8)')
c.debug

              if ( interr .eq. 1 ) go to 8505
              if ( interr .eq. 2 ) go to 8511
 
              watrn1 = watrn1 + fnzlf2
              if ( unsym ) watrn4 = watrn4 + fnzlf2
 
              qincor = .true.
 
          end if
 
      end if
 
      mxused = max ( mxused, inuse )
c.debug
c     write(6,'("at check point 4")')
c.debug
 
c     --------------------------------------------------
c     ... if lnz is in core then close wafil1 and wafil4
c     --------------------------------------------------
 
      if ( qincor .and. wafil1 .ge. 0 ) then
 
          call xdslw7 ( wafil1, interr )
          if ( interr .ne. 0 ) go to 8505

          if ( unsym ) then
              call xdslw7 ( wafil4, interr )
              if ( interr .ne. 0 ) go to 8511
          end if
 
      end if
 
c     ----------------------------------------------------------
c     ... build inverse of factorization permutation and correct
c         new global indices
c     ----------------------------------------------------------
 
      invsup = inuse + 1
      inuse  = inuse + ineqns
 
      call xislog ( neqns, work(sup), work(invsup) )
 
      if ( pvttol .gt. 0. )
     1    call xislfi ( nsnin2, work(lindx2), work(invsup) )
c.debug
c     write(6,'("at check point 5 - xlind2, lindx2 = ", 2i8)') 
c    1                               xlind2, lindx2
c     call xdslp5 ( 'diag of the factor', neqns, work(fdiag),
c    1              output )
c     call xislp3 ( 'pivot block information', neqns,
c    1              work(pvtblk), output )
c     call xislp3 ( 'factorization permutation', neqns,
c    1              work(sup), output )
c     call xislp3 ( 'factorization inverse permutation', neqns,
c    1              work(invsup), output )
c     call xislp3 ( 'new supernode partition',
c    1              nsuper+1, work(xsup2), output )
c     call xislp3 ( 'panel partition',
c    1              npanel+1, work(xpanel), output )
c     call xislp3 ( 'new index pointer',
c    1              nsuper+1, work(xlind2), output )
c     call xislp3 ( 'new lindxg', nsnin2, work(lindx2), output )
c     if ( qincor ) then
c         call xdslp5 ( 'lnz', l, work(lnz), output )
c     end if
c.debug
 
c     -------------------------------------
c     ... compute condition number estimate
c     -------------------------------------
 
      temp   = inuse + 1
      mxused = max ( mxused, temp + neqns + 2*mxlfrt -1 )

      if ( work(qcndnc) .eq. 0. ) then

          cndnum = 0.

      else
c.debug
c     write(6,'("before condition number estimation")')
c.debug
 
          if ( qincor ) then
c.debug
c     write(6,'("before xdslc2")')
c.debug
 
              call xdslc2( unsym, neqns, nsuper, work(xsup2),
     1                     work(xpanel), work(xlind2), work(lindx2),
     2                     work(pvtblk), work(fdiag), work(lnz), tnorma,
     3                     work(temp), work(temp+neqns), cndnum   )
c.debug
c     write(6,'("after  xdslc2")')
c.debug
 
          else
 
              front  = temp + neqns + 2*mxlfrt
              wkreqd = front + mpanel - 1
     
c.debug
c     write(6,'("before memory length check no. 3")')
c.debug
              if ( wkreqd .gt. lwork ) go to 8501

              ltemp  = lwork - front + 1
 
              mxused = max ( mxused, wkreqd )
c.debug
c     write(6,'("before xdslc3 - ltemp = ",i8)') ltemp
c.debug
 
              call xdslc3( unsym, neqns, nsuper, work(xsup2),
     1                     work(xpanel), work(xlind2),  work(lindx2),
     2                     work(pvtblk), work(fdiag),   wafil1, 
     2                     wafil4      , ltemp,
     3                     work(front),  tnorma,        work(temp),
     4                     work(temp+neqns), cndnum,  interr )
c.debug
c     write(6,'("after  xdslc3 - ltemp = ",i8)') ltemp
c.debug
 
              if ( interr .eq. -1 ) go to 8505
              if ( interr .eq. -2 ) go to 8511
 
              watrn1 = watrn1 + fnzlf2
              if ( unsym ) watrn4 = watrn4 + fnzlf2
 
          end if

c.debug
c     write(6,'("after  condition number estimation")')
c.debug

      end if
 
c     --------------------
c     ... set up for solve
c     --------------------
 
      call xislfj ( nsuper, work(xsup2), work(xpanel), work(xlind2),
     1              mxlfrt, mfront, mxtotf )
 
      slvbsz = work(qslvbs)
      slvspc = max ( neqns, slvbsz*mxlfrt )
      if ( .not. qincor ) slvspc = max ( neqns, mxtotf + slvbsz*mxlfrt )
 
      work ( qmxnin ) = mxlfrt
      work ( qmfron ) = mfront
      work ( qmxtot ) = mxtotf
 
      work ( qtemp2 ) = inuse + 1
 
      inuse  = inuse + slvspc
      mxused = max ( mxused, inuse )
 
      work ( qneeds ) = inuse
 
c     --------------
c     ... finish up.
c     --------------
 
      if ( msglvl .ge. 3 ) then
          call xdslp5 ( 'diag of the factor', neqns, work(fdiag),
     1                  output )
          call xislp3 ( 'pivot block information', neqns,
     1                  work(pvtblk), output )
          call xislp3 ( 'factorization permutation', neqns,
     1                  work(sup), output )
          call xislp3 ( 'factorization inverse permutation', neqns,
     1                  work(invsup), output )
          call xislp3 ( 'new supernode partition',
     1                  nsuper+1, work(xsup2), output )
          call xislp3 ( 'panel partition',
     1                  npanel+1, work(xpanel), output )
          call xislp3 ( 'new index pointer',
     1                  nsuper+1, work(xlind2), output )
      end if
 
      if ( msglvl .ge. 4 ) then
          call xislp3 ( 'new lindxg', nsnin2, work(lindx2), output )
          if ( qincor ) then
              call xdslp5 ( 'lnz', l, work(lnz), output )
          end if
      end if
 
      call xdslt2 ( t1, w1, t2, w2 )
 
      if ( msglvl .ge. 1 ) write ( output, 71000 ) t2, w2, cndnum
c.debug
c     write(6,'("elapsed time for postprocessing = ", f15.6)') t2-t4
c.debug
 
c     ---------------------------------------------
c     ... store information into.CMNication area
c     ---------------------------------------------
 
      inuse  =  invsup + ineqns - 1
 
      work ( qstage ) = 50
      work ( qinuse ) = inuse
      work ( qmxuse ) = mxused
 
      work ( qnzlf2 ) = fnzlf2
      work ( qnsni2 ) = nsnin2
      work ( qmxstk ) = lstack
 
      work ( qsup   ) = sup
      work ( qinvsu ) = invsup
      work ( qpivot ) = pvtblk
      work ( qxsup2 ) = xsup2
      work ( qxpanl ) = xpanel
      work ( qxlnd2 ) = xlind2
      work ( qlndg2 ) = lindx2
      work ( qdiag  ) = fdiag
      work ( qlnz   ) = lnz
 
      work ( qfctop ) = fctops
      work ( qslvop ) = slvops
 
      work ( qfcttm ) = t2
      work ( qfctwl ) = w2
 
      work ( qcndnm ) = cndnum
 
      work ( qsqtr1 ) = work ( qsqtr1 ) + sqtrn1
      work ( qwaln1 ) = max ( work ( qwaln1 ), walen1 )
      work ( qwatr1 ) = work ( qwatr1 ) + watrn1
      work ( qwaln2 ) = max ( work ( qwaln2 ), walen2 )
      work ( qwatr2 ) = work ( qwatr2 ) + watrn2
      work ( qwaln4 ) = max ( work ( qwaln4 ), walen4 )
      work ( qwatr4 ) = work ( qwatr4 ) + watrn4
      work ( qwaln5 ) = max ( work ( qwaln5 ), walen5 )
      work ( qwatr5 ) = work ( qwatr5 ) + watrn5
 
      if ( qincor ) then
          work(qtemp1) = 0.
      else
          work(qtemp1) = 1.
      end if
 
      go to 9000
 
c---------------------------------------------------------------------
 
c     --------------------------------------------------
c     ... error processing.  last three digits of labels
c         correspond to code returned to user.
c     --------------------------------------------------
 
c     ---------------------------------
c     ... incorrect processing sequence
c     ---------------------------------
 
 8500 continue
      error = -500
      call hherr ( 3, 'xdslfa', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88500 ) error, stage
      go to 9000
 
c     ---------------------------------------
c     ... insufficient storage for this stage
c     ---------------------------------------
 
 8501 continue
      error = -501
      call hherr ( 2, 'xdslfa', error, wkreqd )
      if ( msglvl .gt. 0 ) write ( output, 88501 ) error, wkreqd, lwork
      needs = wkreqd
      go to 9000
 
c     ----------------------------
c     ... storage allocation error
c     ----------------------------
 
 8502 continue
      error = -502
      call hherr ( 3, 'xdslfa', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88502 ) error
      work(qstage) = -1
      go to 9000
 
c     -----------------------------------------------
c     ... exactly zero column found in reduced matrix
c         (pivoting mode only)
c     -----------------------------------------------
 
 8503 continue
      error = -503
      call hherr ( 3, 'xdslfa', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88503 ) error
      work(qstage) = 40
      go to 9000
 
c     ------------------------------------
c     ... i/o error on sequential i/o file
c     ------------------------------------
 
 8504 continue
      error = -504
      call hherr ( 3, 'xdslfa', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88504 ) error, iofile
      go to 9000
 
c     ----------------------------
c     ... i/o error on file wafil1
c     ----------------------------
 
 8505 continue
      error = -505
      call hherr ( 3, 'xdslfa', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 87000 ) error, wafil1
      go to 9000
 
c     ----------------------------
c     ... i/o error on file wafil2
c     ----------------------------
 
 8506 continue
      error = -506
      call hherr ( 3, 'xdslfa', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 87000 ) error, wafil2
      go to 9000
 
c     ----------------------------------------------------
c     ... out of space -- ran out of room to store lnz and 
c         unz with wafil1 disabled.
c     ----------------------------------------------------
 
 8507 continue
      error = -507
      call hherr ( 3, 'xdslfa', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88507 ) error, lwork
      needs = wkreqd
      go to 9000
 
c     ---------------------------------------------
c     ... out of space -- ran out of room to store
c         stack with wafil2 disabled.
c     ---------------------------------------------
 
 8508 continue
      error = -508
      call hherr ( 3, 'xdslfa', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88508 ) error, lwork
      needs = wkreqd
      go to 9000
 
c     --------------------------------------------------------
c     ... amount of allowed fill exceeded when fill monitoring
c         has been activiated.
c     --------------------------------------------------------
 
 8509 continue
      error = -509
      call hherr ( 3, 'xdslfa', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88509 ) error, work(qextfl)
      go to 9000
 
c     --------------------------------------------------------
c     ... negative pivot encountered while monitoring negative
c         pivots (unpivoted factorization only).
c     --------------------------------------------------------
 
 8510 continue
      error = -510
      call hherr ( 3, 'xdslfa', error, 0 )
      itemp = work(qnpcnt+2)
      if ( msglvl .gt. 0 ) write ( output, 88510 ) error, itemp
      go to 9000

c     ----------------------------
c     ... i/o error on file wafil4
c     ----------------------------
 
 8511 continue
      error = -511
      call hherr ( 3, 'xdslfa', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 87000 ) error, wafil4
      go to 9000
 
c     ----------------------------
c     ... i/o error on file wafil5
c     ----------------------------
 
 8512 continue
      error = -512
      call hherr ( 3, 'xdslfa', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 87000 ) error, wafil5
      go to 9000

c     ------------------------------------------------------
c     ... numerical failure -- in symmetric pivoting mode,
c         unable to find stable pivot in panel of root front
c     ------------------------------------------------------

 8513 continue
      error = -513
      call hherr ( 3, 'xdslfa', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88513 ) error
      go to 9000

c     ----------------------------------------------
c     ... numerical failure -- in non-pivoting mode,
c         found exactly zero diagonal entry (pivot) 
c     ----------------------------------------------

 8514 continue
      error = -514
      call hherr ( 3, 'xdslfa', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88514 ) error
      go to 9000

c     -----------------------------------------------
c     ... undocumented error from lower level routine
c     -----------------------------------------------

 8599 continue
      error = -599
      call hherr ( 3, 'xdslfa', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88599 ) error
      go to 9000
 
c---------------------------------------------------------------------
 
c     ------------------------
c     ... end of module xdslfa
c     ------------------------
 
 9000 continue
      return
 
c---------------------------------------------------------------------
 
c     -----------
c     ... formats
c     -----------
 
60000 format ( /1x, '=============================================='
     1         /1x, '= multifrontal numerical factorization phase ='
     2         /1x, '==============================================' )
 
70000 format ( /5x, 'cpu  time for numeric value input       = ', f15.6
     1         /5x, 'wall time for numeric value input       = ', f15.6)
 
71000 format ( /5x, 'cpu  time for numeric factorization    = ', f15.6
     1         /5x, 'wall time for numeric factorization    = ', f15.6
     2         /5x, 'one norm condition number estimate     = ',
     3                                                         1pd15.4 )
 
c     ... common format for errors 505, 506, 511, 512

87000 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslfa encountered i/o error'
     2         /5x, 'on word addressable i/o file no. ', i15 )

c     ... formats for error given by last three digits of label
 
88500 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslfa executed in an'
     2         /5x, 'incorrect sequence.  ',
     3              'current stage = ', i10, 5x, 'should be 40.' )
 
88501 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslfa requires ', i15,
     2         /5x, 'words of workspace and has only',
     3              i15,  ' available.' )
 
88502 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslfa aborted'
     2         /5x, 'with storage allocation error.' )
 
88503 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslfa encountered',
     2         /5x, 'a numerically singular matrix.' )
 
88504 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslfa encountered i/o error'
     2         /5x, 'on i/o file no. ', i15 )
 
88507 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslfa ran out of storage for nonzeroes of ',
     2         /5x, 'factor with i/o file wafil1 disabled.'
     3         /5x, 'words of workspace currenly available is ',
     4              i15,  '.' )
 
88508 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslfa ran out of storage for computational ',
     2         /5x, 'stack with i/o file wafil2 disabled.'
     3         /5x, 'words of workspace currenly available is ',
     4              i15,  '.' )
 
88509 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslfa exceed factorization fill limitation.',
     2         /5x, 'allowed fill growth factor is ', f15.3,  '.' )
 
88510 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslfa negative pivot encounterd with ',
     2         /5x, 'negative pivot monitoring activated.'
     3         /5x, 'first such pivot encountered is ',
     4              i15,  '.' )
 
88513 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslfa unable to find stable pivot in ',
     2         /5x, 'a panel of the root front.'
     3         /5x, 'increase panel size.' )

88514 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslfa found exact zero diagonal in   ',
     2         /5x, 'non-pivoting mode.' )

88599 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslfa received undocumented error    ',
     2         /5x, 'return from lower level routine'
     3         /5x, 'contact customer support.')

c---------------------------------------------------------------------
 
      end
