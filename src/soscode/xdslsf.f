      subroutine xdslsf ( work,   lwork,  needs,  needmn, error )
c
c     purpose
c     -------
c
c     xdslsf is the top level driver for symbolic factorization
c     phase.
c
c     created         26-jan-89   -- rgg --
c     modified        29-oct-90   -- mlc -- calls to xislvw protected
c                                           against nzla .le. 0
c     modified        08-feb-91   -- rgg -- mods to allow no i/o to
c                                           sqfil1
c     modified        08-feb-91   -- rgg -- added use of xdslil
c     modified        17-dec-96   -- rgg -- converted to used compressed
c                                           adjacency structure
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
c     needs       i   amount of workspace required for the
c                     in-core numeric factorization
c     needmn      i   amount of workspace required for the
c                     out-of-core numeric factorization
c
c     error       i   error flag
c                     =    0  normal return
c                     = -300  incorrect processing path.
c                     = -301  lwork not large enough.
c                     = -302  error return from subroutine xisls2
c                     = -303  i/o error on sqfil2.
c
c---------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      integer             lwork,  error,  needs,  needmn
 
      double precision    work(*)
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------
 
      integer             adjncy, alist,  arindx, bmxtyp, collnk,
     1                    ineqn,  ineqn1, insnin,
     2                    inspr1, insupr, inuse,  invp,   inzla,
     3                    length, lindxg, lindxl, marker,
     4                    maxwrk, mrglnk, msglvl, mtxcol, mfront,
     5                    mxused, mxanzf, mxlfrt, l     , mxtotf,
     6                    imstck, mffdyn, mfsdyn, mstack
 
      integer             nadj,   nassmb, neqns,  nsnind, 
     1                    nsuper, nzla,   output, perm,
     2                    relloc, snmson, stage,  stkstr,
     3                    supmap, temp,   wkreqd, xadj,   xlindx,
     4                    xsup  , lwkval, sqfil2, valend, valinp,
     5                    slvbsz, slvspc, fupdsz, kk

      integer             amncor, bmncor, inzlb,  mincor, nzlb

      integer             cmpmap, inzcmp, ncomp,  nzcomp, 
     1                    rowlst, xrowls
 
      double precision    pvttol, t1,     t2,     w1,     w2,
     1                    fnzlf
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      integer             xdslil, xdslni
 
      external            xisls1, xislmv, xdslil, xislog, xislp1,
     1                    xdslni, xislp3, xisls2, xdslt1, xdslt2
 
c---------------------------------------------------------------------
 
      call xdslt1 ( t1, w1 )
      error  = 0
 
c     ----------------------------------------------------
c     ... extract information from the.CMNication area.
c     ----------------------------------------------------
 
      stage  = work ( qstage )
      msglvl = work ( qmsglv )
      output = work ( qoutpu )
 
      if ( msglvl .ge. 2 ) write ( output, 68000 )
 
      if ( stage .ne. 20 ) go to 8000
 
      neqns  = work ( qneqns )
      nzla   = work ( qnzla  )
      ncomp  = work ( qncomp )
      nzcomp = work ( qnzcmp )
      fnzlf  = work ( qnzlf  )
      nsuper = work ( qnsupe )
      nsnind = work ( qnsnin )
      snmson = work ( qsnmso )
      stkstr = work ( qstkst )
      mfront = work ( qmfron )
      mxtotf = work ( qmxtot )

      nadj    = 2 * nzcomp         

      xrowls = work ( qxrwls )
      rowlst = work ( qrwlst )
      cmpmap = work ( qcmpmp )
      xadj   = work ( qxadj  )
      adjncy = work ( qadjnc )
      perm   = work ( qperm  )
      invp   = work ( qinvp  )
      xsup   = work ( qxsup  )
      nassmb = work ( qnassm )
      xlindx = work ( qxlind )

      bmxtyp = work ( qbmxty )
 
      inspr1  = xdslni ( nsuper + 1 )
 
      wkreqd = xlindx + inspr1 - 1
      maxwrk = wkreqd
c.debug
c     write(6,'("at workspace check point no. 1")')
c.debug
      if ( wkreqd .gt. lwork ) go to 8100
 
c     ------------------------------
c     ... optionally print out input
c     ------------------------------
 
      if ( msglvl .ge. 4 ) then
          call xislp1 ( 'compressed row list',
     1                   ncomp, work(xrowls), work(rowlst), output )
          call xislp3 ( 'compression map',
     1                   neqns, work(cmpmap), output )
          call xislp3 ( 'permutation vector',
     1                   neqns, work(perm), output )
          call xislp3 ('supernode partition vector',
     1                   nsuper+1, work(xsup), output )
          call xislp3 ( 'number of assemblies',
     1                   nsuper, work(nassmb), output )
          call xislp3 ('supernodal index pointers',
     1                   nsuper+1, work(xlindx), output )
      endif
 
 
c     -----------------------------------------------
c     ... compute other pointers into the work array.
c     -----------------------------------------------
 
      insupr  = xdslni ( nsuper     )
      inspr1  = xdslni ( nsuper + 1 )
      ineqn   = xdslni ( neqns      )
      ineqn1  = xdslni ( neqns  + 1 )
      inzla   = xdslni ( nzla       )
      inzcmp  = xdslni ( nzcomp     )
      insnin  = xdslni ( nsnind     )
 
      temp    = xlindx + inspr1
      wkreqd  = temp + ineqn - 1
      maxwrk  = max ( maxwrk, wkreqd )
c.debug
c     write(6,'("at workspace check point no. 2")')
c.debug
      if ( wkreqd .gt. lwork ) go to 8100
 
c     ----------------------------------------------------
c     convert full adjacency of original compressed matrix
c     to lower triangle of permuted matrix
c     ----------------------------------------------------
 
      call  xisls1 ( neqns, neqns+1, nadj, nsuper, ncomp, work(xrowls),
     1               work(rowlst), work(cmpmap), work(xsup),
     2               work(perm), work(invp), work(xadj),
     3               work(adjncy), mtxcol, mxanzf, work(temp) )

      if ( work(qmxtyp) .eq. 2. ) mxanzf = 2 * mxanzf
c.debug
c     write(6,'("after xisls1 - mxanzf = ", i8)')
c.debug
 
      if  ( msglvl .ge. 4 )  then
          call xislp1 ( 'lower triangle of papt',
     1                 neqns, work(xadj), work(adjncy), output )
      endif
 
      wkreqd = temp - 1
 
c  =====================================================================
 
c  ---------------------------------------------------------------------
c  compress the data.  note that adjncy is now half as big as previously
c  ---------------------------------------------------------------------
 
      temp   = perm
      perm   = adjncy + inzcmp
      invp   = perm   + ineqn
      xsup   = invp   + ineqn
      nassmb = xsup   + inspr1
      xlindx = nassmb + insupr
 
      wkreqd = xlindx + inspr1 - 1
 
      length = xdslil ( wkreqd - perm + 1 )
 
      call xislmv ( length, work, xdslil(temp-1)+1, xdslil(perm-1)+1 )
 
      if ( msglvl .ge. 4 ) then
          call xislp3 ( 'permutation vector',
     1                  neqns, work(perm), output )
          call xislp3 ( 'inverse permutation vector',
     1                  neqns, work(invp), output )
          call xislp3 ( 'supernode partition',
     1                  nsuper+1, work(xsup), output )
          call xislp3 ( 'number of assemblies',
     1                  nsuper, work(nassmb), output )
          call xislp3 ( 'supernodal index pointer',
     1                  nsuper+1, work(xlindx), output )
          call xislp1 ( 'lower triangle of papt',
     1                 neqns, work(xadj), work(adjncy), output )
      endif
 
c  =====================================================================
 
c     -----------------------------------------------------
c     symbolic factorization
c     ... note that arindx is no longer generated by xisls2
c     -----------------------------------------------------
 
      lindxg = wkreqd + 1
      lindxl = lindxg + insnin
      arindx = lindxl + insnin
      mrglnk = arindx 
      collnk = mrglnk + insupr
      marker = collnk + ineqn1
      alist  = marker + ineqn
      supmap = alist  + ineqn
      relloc = supmap + ineqn
 
      wkreqd = relloc + ineqn - 1
      maxwrk = max ( maxwrk, wkreqd )
 
c.debug
c     write(6,'("at workspace check point no. 3")')
c.debug
      if  ( wkreqd .gt. lwork )  go to 8100
 
      call xisls2 ( neqns, nzcomp, work(xadj), work(adjncy), nsuper,
     .              work(xsup),   nsnind, work(xlindx), work(perm), 
     .              work(cmpmap), work(xrowls), 
     .              work(rowlst), work(lindxg),
     .              work(lindxl), .false., work(arindx), work(mrglnk),
     .              work(collnk), work(marker), mtxcol, work(alist),
     .              work(supmap), work(relloc), error )
 
      if  ( error .ne. 0 )  go to 8200
 
c     ---------------------------------------------
c     optionally print out the symbolic information
c     ---------------------------------------------
 
      if  ( msglvl .ge. 4 )  then
          call xislp1 ( 'global indices of l',
     1                 nsuper, work(xlindx), work(lindxg), output )
          call xislp1 ( 'relative indices of l',
     1                 nsuper, work(xlindx), work(lindxl), output )
      endif

c     -------------------------------------
c     ... write lindxl and lindxg to sqfil2
c     -------------------------------------

      sqfil2 = work ( qsqfl2 )

      if ( sqfil2 .gt. 0 .and. nsnind .gt. 0 ) then

          call xislvo ( sqfil2, error )
          if ( error .ne. 0 ) go to 8300

          call xislvw ( sqfil2, nsnind, work(lindxl), error )
          if ( error .ne. 0 ) go to 8300

          call xislvw ( sqfil2, nsnind, work(lindxg), error )
          if ( error .ne. 0 ) go to 8300

          work(qsqln2) = 2 * insnin
          work(qsqtr2) = 2 * insnin

      end if
 
c     ------------------------------------------------------
c     ... collapse memory to remove compression information
c         and adjacency structure
c     ------------------------------------------------------

      temp   = xadj
      xadj   = lncomm + 1

      call xislmv ( neqns+1, work, xdslil(temp-1)+1, xdslil(xadj-1)+1 )
 
      temp   = perm
      perm   = xadj   + ineqn1
      invp   = perm   + ineqn
      xsup   = invp   + ineqn
      nassmb = xsup   + inspr1
      xlindx = nassmb + insupr

      if ( sqfil2 .gt. 0 ) then
          lindxg = 0
          lindxl = 0
          inuse  = xlindx + inspr1 - 1
      else
          lindxg = xlindx + inspr1
          lindxl = lindxg + insnin
          inuse  = lindxl + insnin - 1
      end if
 
      length = xdslil ( inuse  - perm + 1 )
 
      call xislmv ( length, work, xdslil(temp-1)+1, xdslil(perm-1)+1 )
 
c     -------------------------------------------------------
c     ... set up for matrix value input phase.
c         wkreqd is the minimum amount of storage needed for
c         value input phase.  valinp is how much storage is 
c         required to input the matrices.  valend is how much 
c         storage is required at the end of value input .
c     -------------------------------------------------------
 
      needs  = work ( qlwkif )
      needmn = work ( qlwkof )

      if ( work ( qsavea ) .eq. 0 ) then
          kk = 1
      else
          kk = 2
      end if

      valinp = inuse + 3*neqns + 6*ineqn + 2*xdslni(mxanzf) + mxanzf
      
      if ( bmxtyp .eq. 1 ) then
          valinp = valinp + 2*neqns + 4*ineqn
      else if ( bmxtyp .eq. 3 ) then
          valinp = valinp + neqns
      end if

      mincor = work ( qmncor )
      amncor = mod ( mincor, 10 )
      bmncor = mincor / 10

      valend = inuse + neqns
      if ( bmxtyp .le. 3 ) valend = valend + neqns

      if ( amncor .eq. 0 ) then
          valend = valend + nzla + kk * inzla
      end if

      if ( bmxtyp .eq. 1 ) then
          valend = valend + ineqn1
          if ( bmncor .eq. 0 ) then
              nzlb   = work(qnzlb)
              inzlb  = xdslni(nzlb)
              valend = valend + nzlb + 2*inzlb
          end if
      end if

      if ( sqfil2 .gt. 0 ) valend = valend + insnin

      lwkval = max ( valinp, valend )

      wkreqd = lwkval

      work ( qtemp0 ) = valinp
      work ( qtemp1 ) = valend
c.debug
c     write(6,'("inuse , neqns , ineqn , mxanzf, lwkval = ", 5i8)')
c    1            inuse , neqns , ineqn , mxanzf, lwkval 
c     write(6,'("nzla  , nzlb  , valinp, valend         = ", 5i8)')
c    1            nzla  , nzlb  , valinp, valend         
c     write(6,'("needs , needmn                         = ", 5i8)')
c    1            needs , needmn                         
c.debug

c     -------------------------------------------------
c     ... adjust needs and needmn for solve block size.
c     -------------------------------------------------

      mxlfrt = ( -1. + sqrt ( 1. + 8. * mfront ) ) / 2.
      slvbsz = work(qslvbs)

      l      = slvbsz*mxlfrt - neqns
      if ( l .gt. 0 ) then 
          needs  = needs  + l
          needmn = needmn + l
      end if
 
c     ----------------------------------------------------------
c     ... adjust needs and needmn if saving matrix is activated.
c     ----------------------------------------------------------

      if ( work ( qsavea ) .ne. 0. ) then

          l               = xdslni(nzla)
          if ( needs .ne. 0 ) then
              needs           = needs + l + neqns
          end if
          needmn          = needmn + l + neqns
          work ( qlwkif ) = needs
          work ( qlwkof ) = needmn

      end if
 
c     ---------------------------------------------------------
c     ... adjust needs and needmn if minimum core is activated.
c         note:  a buffer of neqns was already allowed for.
c                length is the additional buffer space needed.
c     ---------------------------------------------------------

      fupdsz = work(qfupds)
 
      if ( work ( qmncor ) .ne. 0 ) then

          length = kk * xdslni(mxanzf) + mxanzf 

          if ( bmxtyp .eq. 1 .and. work(qmncor) .eq. 22 ) 
     1       length = 2 * length

          l = ( ( nzla + kk * xdslni(nzla) ) - length )
     1      + ( neqns - max ( neqns, fupdsz*mxlfrt ) )
c.debug
c     write(6,'("in xdslsf - needs")')
c     write(6,'("length, fnzlf , mfront, mxlfrt         = ", 5i8)')
c    1            length, fnzlf , mfront, mxlfrt
c     write(6,'("needs , l     , wkreqd, nzla           = ", 5i8)')
c    1            needs , l     , wkreqd, nzla
c.debug
 
          if ( needs .ne. 0 ) needs  = max ( needs  - l, wkreqd )

c.debug
c     write(6,'("in xdslsf - needmn")')
c     write(6,'("length, fnzlf , mfront, mxlfrt         = ", 5i8)')
c    1            length, fnzlf , mfront, mxlfrt
c     write(6,'("needmn, l     , wkreqd, nzla           = ", 5i8)')
c    1            needmn, l     , wkreqd, nzla
c.debug

          needmn = max ( needmn - l, wkreqd )
c.debug
c     write(6,'("needmn                                 = ", 5i8)')
c    1            needmn                 
c.debug
 
      end if
 
      call xdslt2 ( t1, w1, t2, w2 )

      if ( msglvl .ge. 1 ) write ( output, 81000 ) mtxcol, mxanzf, 
     1                                             t2, w2, maxwrk,
     2                                             lwkval
          
c     ------------------------------------------------------------
c     ... compute workspace requirements for lanczos.  one is for
c         workspace active during a factorization and the other is
c         for during a recurrence/solve.
c         to prevent overflow of 32 bit integer arithmetics mfsdyn
c         and mffdyn do not include storage for the factorization,
c         the multifrontal stack, or solve temp space.  these
c         are compensated for in xdseex.
c     ------------------------------------------------------------
 
      mstack = work ( qmstac )
      imstck = xdslni ( mstack )
 
      pvttol = work ( qpvttl )

      slvbsz = work(qslvbs)
      slvspc = max ( neqns, slvbsz*mxlfrt )
 
      mfsdyn = 3*ineqn + 2*inspr1 + slvspc
      mffdyn = 6*imstck + 2*neqns + 2*ineqn + 2*inspr1 

      if ( pvttol .gt. 0. ) then 
          mfsdyn = mfsdyn + insnin
          mffdyn = mffdyn + insnin
      end if

      work ( qmffdy ) = mffdyn
      work ( qmfsdy ) = mfsdyn
 
c     ---------------------------------------------
c     ... store information into.CMNication area
c     ---------------------------------------------
 
      stage      = 30
      work ( qstage ) = stage
      work ( qinuse ) = inuse
      work ( qneeds ) = needs
 
      mxused          = work ( qmxuse )
      work ( qmxuse ) = max ( mxused, maxwrk )
 
      work ( qmxanz ) = mxanzf
 
      work ( qxadj  ) = xadj
      work ( qperm  ) = perm
      work ( qinvp  ) = invp
      work ( qxsup  ) = xsup
      work ( qnassm ) = nassmb
      work ( qxlind ) = xlindx
      work ( qlndxl ) = lindxl
      work ( qlndxg ) = lindxg
 
      work ( qsfctm ) = t2
      work ( qsfcwl ) = w2
 
      go to 9000
 
c---------------------------------------------------------------------
 
c     --------------
c     ... error trap
c     --------------
 
c     ---------------------------------
c     ... incorrect processing sequence
c     ---------------------------------
 
 8000 continue
      error = -300
      call hherr ( 3, 'xdslsf', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88000 ) error, stage
      go to 9000
 
c     ---------------------------------------
c     ... insufficient storage for this stage
c     ---------------------------------------
 
 8100 continue
      error = -301
      call hherr ( 2, 'xdslsf', error, wkreqd )
      if ( msglvl .gt. 0 ) write ( output, 88100 ) error, wkreqd, lwork
      go to 9000
 
c     ---------------------
c     ... error from xisls2
c     ---------------------
 
 8200 continue
      error = -302
      call hherr ( 3, 'xdslsf', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88200 ) error
      work(qstage) = -1
      go to 9000
 
c     ----------------------------
c     ... i/o error on file sqfil2
c     ----------------------------
 
 8300 continue
      error = -303
      call hherr ( 3, 'xdslsf', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88300 ) error, sqfil2
      work(qstage) = -1
      go to 9000
 
c---------------------------------------------------------------------
 
c     ------------------------
c     ... end of module xdslsf
c     ------------------------
 
 9000 continue
      return
 
c---------------------------------------------------------------------
 
c     -----------
c     ... formats
c     -----------
 
68000 format ( /1x, '============================================='
     1         /1x, '= multifrontal symbolic factorization phase ='
     2         /1x, '=============================================' )
 
81000 format ( /5x, 'maximum nonzeroes in a column of papt   = ', i15
     1         /5x, 'maximum nonzeroes in a supernode        = ', i15 
     2         /5x, 'cpu  time for symbolic factorization    = ', f15.6
     3         /5x, 'wall time for symbolic factorization    = ', f15.6
     4         /5x, 'max. amount of storage used             = ', i15
     5         /5x, 'amt. of storage needed for value input  = ', i15  )
 
88000 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslsf executed in an'
     2         /5x, 'incorrect sequence.  current stage = ', i10,
     3          5x, 'should be 20.' )
 
88100 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslsf requires ', i15
     2         /5x, 'words of workspace and has only ', i15,
     3              ' available.' )
 
88200 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslsf encountered'
     2         /5x, 'inconsistent matrix structure during symbolic',
     3              'factorization.' )
 
88300 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslsf encountered i/o error'
     2         /5x, 'on i/o file no. ', i15 )
 
c---------------------------------------------------------------------
 
      end
