      subroutine xdslor ( work,   lwork,  needs,  error )
 
c
c     purpose
c     -------
c
c     xdslor is the top level driver for matrix ordering phase.
c
c     created         26-jan-89   -- rgg --
c     modified        08-may-89   -- rgg -- added indexing parameters
c     modified        29-oct-90   -- mlc -- calls to xislvr and xislvw
c                                           protected against zero
c                                           length input
c     modified        08-feb-91   -- rgg -- mods to allow no i/o to sqfi
c     modified        08-feb-91   -- rgg -- added use of xdslil
c     modified        17-dec-96   -- rgg -- modified to used compressed
c                                           (xadj,adjncy) from release 4
c                                           structure input processing
c     modified        11-mar-04   -- dkw -- adjust txsup storage
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
c     needs       i   amount of workspace required for the next stage.
c
c     error       i   error flag
c                     =    0  normal return
c                     = -200  incorrect processing path.
c                     = -201  lwork not large enough.
c                     = -202  i/o error on sqfil1.
c
c---------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      integer             lwork,  error,  needs
 
      double precision    work(*)
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------
 
      integer             adjncy, brosup, clqsiz, delta,  dhead,
     1                    fellow, fleng,  fstsup, inadj,  ineqn,
     2                    ineqn1, insup1, insup2, insp21, isnmsn,
     3                    invp,   isize,  istore, length, list,
     4                    llist,  lnzcol, locxad, locadj, lstfel,
     5                    lstson, lwkifs, lwkofs, lwksym, marker,
     6                    maxint, maxwrk, mfront, mmerge, msglvl,
     7                    msizsn, mstack, mupdmt, mxnind, mxrect,
     8                    mxsons, mxtotf, mxtri , maxzer, mxlfrt
 
      integer             nadj,   nassmb, neqns,  nodls1, nodls2,
     1                    nofsub, nsnind, nsntre, nsons,
     2                    nsup1,  nsup2,  nzla,   output,
     3                    parent, parsup, perm,   qsize,  rstore,
     4                    sndpth, snmson, snroot, snsize, sqfil1,
     5                    stage,  stksti, stkstr, temp,   tperm,
     6                    txsup,  tsldon, uleng,  val,    xadj,
     7                    xcliqu, xlindx, xsup1,  xsup2,  wkreqd,
     8                    zeroes, sqfil2, insp11

      integer             cmpmap, incmp,  incmp1, ncomp,  nzcomp,
     1                    rowlst, xrowls, fpnlsz, fupdsz
 
      logical             incore, unsym
 
      double precision    t1,     t2,     w1,     w2    , pvttol,
     1                    ftemp , fnzlf
 
      double precision    opscnt(8)
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      integer             xdslil, xdslni
 
      external            xislvr, xislvw, xislo1, xislmv, xdslil,
     1                    xislo7, xislp1, xislod, xdslni, xislp3,
     2                    xislo8, xisloe, xdslob, xisloa, xdslt1,
     3                    xdslt2
 
c---------------------------------------------------------------------
 
      call xdslt1 ( t1, w1 )
      error  = 0
 
c     ----------------------------------------------------
c     ... extract information from the.CMNication area.
c     ----------------------------------------------------
 
      msglvl = work ( qmsglv )
      output = work ( qoutpu )
      stage  = work ( qstage )
      if ( stage .ne. 10 .and. stage. ne. 20 ) go to 8000
 
      if ( msglvl .ge. 2 ) write ( output, 68000 )
 
      neqns  = work ( qneqns )
      nzla   = work ( qnzla  )
      ncomp  = work ( qncomp )
      nzcomp = work ( qnzcmp )
 
      xrowls = work ( qxrwls )
      rowlst = work ( qrwlst )
      cmpmap = work ( qcmpmp )
      xadj   = work ( qxadj  )
      adjncy = work ( qadjnc )
 
      if ( msglvl .ge. 4 ) then
          call xislp1 ( 'compressed row list', ncomp,
     1                  work(xrowls), work(rowlst), output )
          call xislp3 ( 'cmpmap', neqns  , work(cmpmap), output )
          call xislp1 ( 'compressed adjacency structure', neqns,
     1                  work(xadj), work(adjncy), output )
      end if
 
c     -----------------------------------------------
c     ... compute other pointers into the work array.
c     -----------------------------------------------

      nadj   = 2 * nzcomp
 
      inadj  = xdslni ( nadj      )
      ineqn  = xdslni ( neqns     )
      ineqn1 = xdslni ( neqns + 1 )
      incmp  = xdslni ( ncomp     )
      incmp1 = xdslni ( ncomp + 1 )
c$$$c.debug
c$$$      write(6,*)' inadj  =',inadj
c$$$      write(6,*)' ineqn  =',ineqn
c$$$      write(6,*)' ineqn1 =',ineqn1
c$$$      write(6,*)' incmp  =',incmp
c$$$      write(6,*)' incmp1 =',incmp1
c$$$c.debug
 
      perm   = adjncy + inadj
      parsup = perm   + ineqn
      lnzcol = parsup + ineqn
      xsup1  = lnzcol + ineqn
      invp   = xsup1  + incmp1
      dhead  = invp   + ineqn
      qsize  = dhead  + ineqn
      llist  = qsize  + incmp
      marker = llist  + ineqn
      xcliqu = marker + ineqn
      clqsiz = xcliqu + ineqn
 
      wkreqd = clqsiz + ineqn - 1
 
c     ---------------------------------------------------------
c     ... determine if additional storage can or should be used
c     ---------------------------------------------------------
 
      sqfil1 = work ( qsqfl1 )
 
      if ( sqfil1 .gt. 0 ) then

          call xislvo ( sqfil1, error )
          if ( error .ne. 0 ) go to 8200
 
      end if
 
      if ( sqfil1 .le. 0 .or. ( wkreqd + ineqn1 + inadj ) .lt. lwork )
     1then
          incore = .true.
          locxad = wkreqd + 1
          locadj = locxad + ineqn1
          wkreqd = locadj + inadj - 1
      else
          incore = .false.
          locxad = xadj
          locadj = adjncy
      endif
 
      maxwrk  = wkreqd
c.debug
c     write(6,'("before work space ck pt 1")') 
c     write(6,'("wkreqd, lwork = ", 2i8)') wkreqd, lwork
c.debug
      if ( wkreqd .gt. lwork ) go to 8100
 
c     -----------------------------------------------------------------
c     ... ordering destroys its copy of (xadj,adjncy).
c         if incore = .true. make a second copy in memory for ordering.
c         otherwise write (xadj,adjncy) to sqfil1 for later retrieval.
c     -----------------------------------------------------------------
 
      if ( incore ) then
 
          call xislmv ( neqns+1, work, xdslil(xadj  -1)+1,   
     1                                 xdslil(locxad-1)+1 )
          call xislmv ( nadj   , work, xdslil(adjncy-1)+1, 
     1                                 xdslil(locadj-1)+1 )
 
      else
 
c         ---------------------------------
c         ... save the adjacency structure.
c         ---------------------------------
 
          call xislrw ( sqfil1, error )
          if ( error .ne. 0 ) go to 8200
 
          if ( neqns .ge. 0 ) then
             call xislvw ( sqfil1, neqns+1, work(xadj  ), error )
             if ( error .ne. 0 ) go to 8200
          endif
          if ( nadj .gt. 0 ) then
             call xislvw ( sqfil1, nadj   , work(adjncy), error )
             if ( error .ne. 0 ) go to 8200
          endif
 
          work ( qsqln1 ) = ineqn1 + inadj
          work ( qsqtr1 ) = ineqn1 + inadj
 
      end if
 
c     --------------------------------------------------
c     ... generate the multiple minimum degree ordering.
c     --------------------------------------------------
 
      maxint = 2**30 - 1
      delta  = 0
 
      call xislo1 ( ncomp, work(xrowls), work(rowlst), work(cmpmap),
     1              neqns, work(locxad), work(locadj), work(invp),
     1              work(perm), nsup1, work(xsup1), work(lnzcol),
     2              work(parsup), delta, work(dhead), work(qsize),
     3              work(llist), work(marker), maxint, nofsub,
     4              nadj, neqns+1, work(xcliqu), work(clqsiz) )
      
      if ( msglvl .ge. 2 ) write ( output, 68100 ) nsup1, nofsub
 
      if ( msglvl .ge. 4 ) then
 
          call xislp3 ( 'old-to-new permutation array', neqns,
     1                  work(invp), output )
          call xislp3 ( 'new-to-old permutation array', neqns,
     1                  work(perm), output )
          call xislp3 ( 'supernodal pointer array', nsup1+1,
     1                  work(xsup1), output )
          call xislp3 ( 'supernode parent array', nsup1,
     1                  work(parsup), output )
          call xislp3 ( 'factor column nonzero counts', neqns,
     1                  work(lnzcol), output )
 
      end if
 
      if ( .not. incore ) then
 
c         ------------------------------------
c         ... restore the adjacency structure.
c         ------------------------------------
 
          call xislrw ( sqfil1, error )
          if ( error .ne. 0 ) go to 8200
 
          if ( neqns .gt. 0 ) then
             call xislvr ( sqfil1, neqns+1, work(xadj  ), error )
             if ( error .ne. 0 ) go to 8200
          endif
          if ( nadj .gt. 0 ) then
             call xislvr ( sqfil1, nadj   , work(adjncy), error )
             if ( error .ne. 0 ) go to 8200
          endif
 
          work ( qsqtr1 ) = work ( qsqtr1 ) + ineqn1 + inadj
 
      end if
c.debug
c     write(6,'("perm at ck pt 01 ")')
c     call xislp3 ( 'new-to-old permutation array', neqns,
c    1              work(perm), output )
c     call xislp1 ( 'compressed row list', ncomp,
c    1              work(xrowls), work(rowlst), output )
c     call xislp3 ( 'cmpmap', neqns  , work(cmpmap), output )
c     call xislp1 ( 'compressed adjacency structure', neqns,
c    1              work(xadj), work(adjncy), output )
c.debug
 
c     ----------------------------
c     ... collapse working storage
c     ----------------------------
 
      insup1 = xdslni ( nsup1 )
      insp11 = xdslni ( nsup1 + 1 )
      wkreqd = xsup1 + insp11 - 1
 
c     ----------------------------------------------------
c     ... allocate storage for initialization of input for
c         the supernode partition amalgamation routines
c     ----------------------------------------------------
 
      fstsup = wkreqd + 1
      brosup = fstsup + insup1
      parent = brosup + insup1
      nodls1 = parent + ineqn
      snsize = nodls1 + ineqn
      uleng  = snsize + insup1
      fleng  = uleng  + insup1
      nsons  = fleng  + insup1
 
c     ----------------------------------
c     ... is there sufficient work space
c     ----------------------------------
 
      wkreqd = nsons + insup1 - 1
      maxwrk = max ( maxwrk, wkreqd )
c.debug
c     write(6,'("before work space ck pt 2")') 
c     write(6,'("wkreqd, lwork = ", 2i8)') wkreqd, lwork
c.debug
      if  ( wkreqd .gt. lwork )  go to 8100
 
c     ----------------------------------------------------------
c     ... initialize input for the supenode amlgamation routines
c     ----------------------------------------------------------
 
      call  xislo7 ( nsup1, neqns, work(parsup), work(xsup1),
     .               work(lnzcol), work(fstsup), work(brosup),
     .               work(parent), work(nodls1), work(snsize),
     .               work(uleng) , work(fleng) , mxsons      ,
     .               work(nsons)   )
 
      if  ( msglvl .ge. 2 ) write ( output, 68150 ) mxsons
 
      if  ( msglvl .ge. 4 )  then
          call xislp3 ( 'first son vector',
     1                  nsup1, work(fstsup), output )
          call xislp3 ( 'brother vector',
     1                  nsup1, work(brosup), output )
          call xislp3 ( 'nodal tree parent vector',
     1                  neqns, work(parent), output )
          call xislp3 ( 'supernode node list vector',
     1                  neqns, work(nodls1), output )
          call xislp3 ( 'supernode sizes',
     1                  nsup1, work(snsize), output )
          call xislp3 ( 'update matrix sizes',
     1                  nsup1, work(uleng),  output )
          call xislp3 ( 'frontal matrix sizes',
     1                  nsup1, work(fleng),  output )
      endif
 
      wkreqd = nsons - 1
c.debug
c     write(6,'("perm at ck pt 02 ")')
c     call xislp3 ( 'new-to-old permutation array', neqns,
c    1              work(perm), output )
c     call xislp1 ( 'compressed row list', ncomp,
c    1              work(xrowls), work(rowlst), output )
c     call xislp3 ( 'cmpmap', neqns  , work(cmpmap), output )
c     call xislp1 ( 'compressed adjacency structure', neqns,
c    1              work(xadj), work(adjncy), output )
c.debug
 
c     -------------------------------------------------------
c     ... allocate storage for supernode amalgamation routine
c     -------------------------------------------------------
 
      fellow = wkreqd + 1
      zeroes = fellow + insup1
      lstfel = zeroes + insup1
      lstson = lstfel + insup1
 
c     ----------------------------
c     ... is there enough storage?
c     ----------------------------
 
      wkreqd = lstson + insup1 - 1
      maxwrk = max ( maxwrk, wkreqd )
c.debug
c     write(6,'("before work space ck pt 3")') 
c     write(6,'("wkreqd, lwork = ", 2i8)') wkreqd, lwork
c.debug
      if  ( wkreqd .gt. lwork )  go to 8100
 
c     -----------------------
c     ... coalesce supernodes
c     -----------------------
 
      snroot = nsup1
      nsup2  = nsup1

      maxzer = work ( qmxzer )
 
      maxzer = max ( maxzer, 0 )
 
      call  xislo8 ( nsup2 , work(parsup), work(fstsup), work(brosup),
     .               snroot, work(snsize), work(uleng) , work(fleng) ,
     .               maxzer, work(zeroes), work(fellow), work(lstfel),
     .               work(lstson)  )
 
      if  ( msglvl .ge. 2 ) write ( output, 68200 ) nsup2, snroot
 
      if  ( msglvl .ge. 4 )  then
          call xislp3 ( 'supernodal parent vector',
     1                  nsup1, work(parsup), output )
          call xislp3 ( 'first son vector',
     1                  nsup1, work(fstsup), output )
          call xislp3 ( 'brother vector',
     1                  nsup1, work(brosup), output )
          call xislp3 ( 'supernode size vector',
     1                  nsup1, work(snsize), output )
          call xislp3 ( 'size of update matrices',
     1                  nsup1, work(uleng),  output )
          call xislp3 ( 'size of frontal matrices',
     1                  nsup1, work(fleng),  output )
          call xislp3 ( 'number of zeroes in the front matrices',
     1                  nsup1, work(zeroes), output )
          call xislp3 ( 'coalescing linked lists',
     1                  nsup1, work(fellow), output )
      endif
 
      wkreqd = zeroes - 1
c.debug
c     write(6,'("perm at ck pt 03")')
c     call xislp3 ( 'new-to-old permutation array', neqns,
c    1              work(perm), output )
c     call xislp1 ( 'compressed row list', ncomp,
c    1              work(xrowls), work(rowlst), output )
c     call xislp3 ( 'cmpmap', neqns  , work(cmpmap), output )
c     call xislp1 ( 'compressed adjacency structure', neqns,
c    1              work(xadj), work(adjncy), output )
c.debug
 
c     ------------------------------------------
c     ... allocate storage for node list routine
c     ------------------------------------------
 
      insup2 = xdslni ( nsup2     )
      insp21 = xdslni ( nsup2 + 1 )
 
      xsup2  = wkreqd + 1
      nodls2 = xsup2  + insp21
 
c     ----------------------------
c     ... is there enough storage?
c     ----------------------------
 
      wkreqd = nodls2 + ineqn - 1
      maxwrk = max ( maxwrk, wkreqd )
c.debug
c     write(6,'("before work space ck pt 4")') 
c     write(6,'("wkreqd, lwork = ", 2i8)') wkreqd, lwork
c.debug
      if  ( wkreqd .gt. lwork )  go to 8100
 
c     ------------------------------------------------------
c     ... generate node list for the new supernode partition
c     ------------------------------------------------------
 
      call xislo9 ( nsup1, work(parsup), work(fstsup), work(brosup),
     .              snroot, work(snsize), work(fellow), nsup2,
     .              work(xsup1), work(nodls1), work(xsup2),
     .              work(nodls2) )
 
      if  ( msglvl .ge. 4 )  then
          call xislp1 ( 'new supernode partition',
     1                 nsup2, work(xsup2), work(nodls2), output )
      endif
c.debug
c     write(6,'("perm at ck pt 04")')
c     call xislp3 ( 'new-to-old permutation array', neqns,
c    1              work(perm), output )
c     call xislp1 ( 'compressed row list', ncomp,
c    1              work(xrowls), work(rowlst), output )
c     call xislp3 ( 'cmpmap', neqns  , work(cmpmap), output )
c     call xislp1 ( 'compressed adjacency structure', neqns,
c    1              work(xadj), work(adjncy), output )
c.debug
 
c     -----------------------------------------------------------
c     ... allocate storage for supernode elimination tree routine
c     -----------------------------------------------------------
 
      tsldon = wkreqd + 1
 
c     ----------------------------
c     ... is there enough storage?
c     ----------------------------
 
      wkreqd = tsldon + ineqn - 1
      maxwrk = max ( maxwrk, wkreqd )
c.debug
c     write(6,'("before work space ck pt 5")') 
c     write(6,'("wkreqd, lwork = ", 2i8)') wkreqd, lwork
c.debug
      if  ( wkreqd .gt. lwork )  go to 8100
 
c     ---------------------------------------
c     ... generate supernode elimination tree
c         for the new supernode partition
c     ---------------------------------------
 
      call xisloa ( neqns, work(parent), nsup2, work(xsup2),
     .              work(nodls2), work(tsldon), work(parsup),
     .              work(fstsup), work(brosup), snroot, nsntre,
     .              sndpth, snmson  )
 
      if  ( msglvl .ge. 2 ) write ( output, 68300 ) nsntre, sndpth,
     1                                             snmson, snroot
 
      if  ( msglvl .ge. 4 )  then
          call xislp3 ( 'new supernode parent vector',
     1                  nsup2, work(parsup), output )
          call xislp3 ( 'new first son vector',
     1                  nsup2, work(fstsup), output )
          call xislp3 ( 'new brother vector',
     1                  nsup2, work(brosup), output )
      endif
 
      wkreqd = tsldon - 1
c.debug
c     write(6,'("perm at ck pt 05")')
c     call xislp3 ( 'new-to-old permutation array', neqns,
c    1              work(perm), output )
c     call xislp1 ( 'compressed row list', ncomp,
c    1              work(xrowls), work(rowlst), output )
c     call xislp3 ( 'cmpmap', neqns  , work(cmpmap), output )
c     call xislp1 ( 'compressed adjacency structure', neqns,
c    1              work(xadj), work(adjncy), output )
c.debug
 
c     -------------------------------------------------------
c     ... generate statistics for the new supernode partition
c     -------------------------------------------------------
 
      isize  = wkreqd + 1
      nassmb = isize  + insup2
 
c     ----------------------------
c     ... is there enough storage?
c     ----------------------------
 
      wkreqd = nassmb + insup2 - 1
      maxwrk = max ( maxwrk, wkreqd )
c.debug
c     write(6,'("before work space ck pt 6")') 
c     write(6,'("wkreqd, lwork = ", 2i8)') wkreqd, lwork
c.debug
      if  ( wkreqd .gt. lwork )  go to 8100
 
      call xdslob ( neqns, work(lnzcol), nsup2, work(xsup2),
     .              work(parsup), work(fstsup), work(brosup),
     .              work(nassmb), work(nodls2), nzla,
     .              msizsn, mxnind, nsnind, fnzlf,
     .              mfront, mupdmt, mxtri , mxrect, mxtotf,
     .              work(fleng), work(uleng), work(isize), opscnt)
 
      if  ( msglvl .ge. 2 ) then
          mxlfrt = ( -1. + sqrt ( 1. + 8. * mfront ) ) / 2.
          write (output,68400) msizsn, mxnind, msizsn, nsnind, fnzlf,
     1                  mxlfrt, mfront, mupdmt, mxtri , mxrect, mxtotf
      end if
 
      if  ( msglvl .eq. 2 ) then
          write (output,68450) opscnt(1)+opscnt(3),opscnt(5)+opscnt(7)
      else if  ( msglvl .ge. 3 )  then
          write (output,68500) opscnt(1),opscnt(2),opscnt(3),opscnt(4),
     1                  opscnt(1)+opscnt(3)
          write (output,68600) opscnt(5),opscnt(6),opscnt(7),opscnt(8),
     1                  opscnt(5)+opscnt(7)
      endif
 
      if  ( msglvl .ge. 4 )  then
          call xislp3 ( 'size of frontal matrices',
     1                  nsup2, work(fleng), output )
          call xislp3 ( 'size of update matrices',
     1                  nsup2, work(uleng), output )
          call xislp3 ( 'integer storage sizes',
     1                  nsup2, work(isize), output )
      endif
 
      work ( qfctop ) = opscnt(1) + opscnt(3)
      work ( qslvop ) = opscnt(5) + opscnt(7)

      ftemp = msizsn
      if ( work(qpvttl) .ne. 0. ) ftemp = neqns
      work ( qfpnls ) = min ( work ( qfpnls ), ftemp )
      work ( qfupds ) = min ( work ( qfupds ), ftemp )
c.debug
c     write(6,'("perm at ck pt 06")')
c     call xislp3 ( 'new-to-old permutation array', neqns,
c    1              work(perm), output )
c     call xislp1 ( 'compressed row list', ncomp,
c    1              work(xrowls), work(rowlst), output )
c     call xislp3 ( 'cmpmap', neqns  , work(cmpmap), output )
c     call xislp1 ( 'compressed adjacency structure', neqns,
c    1              work(xadj), work(adjncy), output )
c.debug
 
c     ----------------------------
c     ... collapse working storage
c     ----------------------------
 
      temp   = fstsup
      fstsup = parsup + insup2
      call xislmv ( nsup2, work, xdslil(temp-1)+1, xdslil(fstsup-1)+1 )
 
      temp   = brosup
      brosup = fstsup + insup2
      call xislmv ( nsup2, work, xdslil(temp-1)+1, xdslil(brosup-1)+1 )
 
      temp   = uleng
      uleng  = brosup + insup2
      call xislmv ( nsup2, work, xdslil(temp-1)+1, xdslil(uleng-1)+1 )
 
      temp   = fleng
      fleng  = uleng  + insup2
      call xislmv ( nsup2, work, xdslil(temp-1)+1, xdslil(fleng-1)+1 )
 
      temp   = xsup2
      xsup2  = fleng  + insup2
      call xislmv ( nsup2+1, work, xdslil(temp-1)+1, xdslil(xsup2-1)+1 )
 
      temp   = nodls2
      nodls2 = xsup2  + insup2 + 1
      call xislmv ( neqns, work, xdslil(temp-1)+1, xdslil(nodls2-1)+1 )
 
      temp   = isize
      isize  = nodls2 + ineqn
      call xislmv ( nsup2, work, xdslil(temp-1)+1, xdslil(isize-1)+1 )
 
      wkreqd = isize + insup2 - 1
c.debug
c     write(6,'("perm at ck pt 07")')
c     call xislp3 ( 'new-to-old permutation array', neqns,
c    1              work(perm), output )
c     call xislp1 ( 'compressed row list', ncomp,
c    1              work(xrowls), work(rowlst), output )
c     call xislp3 ( 'cmpmap', neqns  , work(cmpmap), output )
c     call xislp1 ( 'compressed adjacency structure', neqns,
c    1              work(xadj), work(adjncy), output )
c.debug
 
c     -----------------------------------------
c     ... obtain optimal reordering of children
c         for the multifrontal method
c     -----------------------------------------
 
      isnmsn = xdslni ( snmson )
 
      list   = wkreqd + 1
      val    = list   + isnmsn
      rstore = val    + isnmsn
      istore = rstore + insup2
      wkreqd = istore + insup2 - 1
      maxwrk = max ( maxwrk, wkreqd )
c.debug
c     write(6,'("before work space ck pt 7")') 
c     write(6,'("wkreqd, lwork = ", 2i8)') wkreqd, lwork
c.debug
      if  ( wkreqd .gt. lwork )  go to 8100
 
      call xislod ( nsup2, work(parsup), work(fstsup), work(brosup),
     .              snroot, work(fleng), work(uleng), work(isize),
     .              snmson, work(list), work(val), work(rstore),
     .              work(istore), stkstr, stksti   )
 
      if  ( msglvl .ge. 2 )  write ( output, 68700 ) stkstr, stksti
 
      if  ( msglvl .ge. 4 )  then
          call xislp3 ( 'new supernode parent vector',
     1                  nsup2, work(parsup), output )
          call xislp3 ( 'new first son vector',
     1                  nsup2, work(fstsup), output )
          call xislp3 ( 'new brother vector',
     1                  nsup2, work(brosup), output )
      endif
 
      wkreqd = list - 1
c.debug
c     write(6,'("perm at ck pt 08")')
c     call xislp3 ( 'new-to-old permutation array', neqns,
c    1              work(perm), output )
c     call xislp1 ( 'compressed row list', ncomp,
c    1              work(xrowls), work(rowlst), output )
c     call xislp3 ( 'cmpmap', neqns  , work(cmpmap), output )
c     call xislp1 ( 'compressed adjacency structure', neqns,
c    1              work(xadj), work(adjncy), output )
c.debug
 
c     ------------------------
c     ... obtain post-ordering
c     ------------------------
 
      nassmb = wkreqd + 1
      xlindx = nassmb + insup2
      tperm  = xlindx + insp21
      txsup  = tperm  + ineqn
      wkreqd = txsup  + max ( ineqn, insp21 ) - 1
      maxwrk = max ( maxwrk, wkreqd )
c.debug
c     write(6,'("before work space ck pt 8")') 
c     write(6,'("wkreqd, lwork = ", 2i8)') wkreqd, lwork
c.debug
      if  ( wkreqd .gt. lwork )  go to 8100
 
      call xisloe ( neqns, nsup2, work(xsup2), work(nodls2),
     .              work(perm), work(parsup), work(fstsup),
     .              work(brosup), snroot, work(nassmb), mstack,
     .              mmerge, work(isize), work(xlindx), work(tperm),
     .              work(txsup), ncomp, work(xrowls), work(rowlst) )
 
c     --------------------------------------------------------------
c     ... note that mstack is increased to allow additional storage
c         on the stack for reducible problems where one of the
c         subproblems is singular
c     --------------------------------------------------------------
 
      mstack = mstack + nsntre - 1
 
      if  ( msglvl .ge. 2 ) write ( output, 68800 ) mstack, mmerge
 
      if  ( msglvl .ge. 3 )  then
          call xislp3 ( 'number of assemblies',
     1                  nsup2, work(nassmb), output )
          call xislp3 ( 'new permutation vector',
     1                  neqns, work(perm), output )
          call xislp3 ( 'symbolic factorization pointers',
     1                  nsup2+1, work(xlindx), output )
          call xislp1 ( 'new supernode partition',
     1                 nsup2, work(xsup2), work(nodls2), output )
      endif
 
      wkreqd = tperm - 1
c.debug
c     stop
c.debug
 
c     -----------------------------------------------------
c     ... collapse working storage (leaving space for invp)
c     -----------------------------------------------------
 
      temp   = xsup2
      invp   = perm + ineqn
      xsup2  = invp + ineqn
      call xislmv ( nsup2+1, work, xdslil(temp-1)+1, xdslil(xsup2-1)+1 )
 
      temp   = nassmb
      nassmb = xsup2  + insp21
      xlindx = nassmb + insup2
      length = xdslil ( insup2 + insp21 )
      call xislmv ( length, work, xdslil(temp-1)+1, xdslil(nassmb-1)+1 )
 
      wkreqd = xlindx + insp21 - 1
 
c     ----------------------------
c     ... form inverse permutation
c     ----------------------------
 
      call xislog ( neqns, work(perm), work(invp) )
 
      if  ( msglvl .ge. 3 )  then
          call xislp3 ( 'new inverse permutation vector',
     1                  neqns, work(invp), output )
      endif
 
c  =====================================================================
 
c     --------------------------------------------------------------
c     ... get the storage requirements for the succeeding program
c         segments and store scalars into work array.
 
c         lwksym - storage requirements for symbolic factorization
c         lwkifs - storage requirements for in-core numeric
c                  factor and solve
c         lwkofs - storage requirements for out-of-core numeric
c                  factor and solve
c     --------------------------------------------------------------
 
      work ( qnzlf  ) = fnzlf
      work ( qnsupe ) = nsup2
      work ( qnsnin ) = nsnind
      work ( qmxnin ) = mxnind
      work ( qsnmso ) = snmson
      work ( qmstac ) = mstack
      work ( qstkst ) = stkstr
      work ( qmxtri ) = mxtri
      work ( qmfron ) = mfront
      work ( qmxtot ) = mxtotf
      work ( qmaxze ) = maxzer
 
      unsym = work ( qmxtyp ) .eq. 2.

      fpnlsz = work(qfpnls)
      fupdsz = work(qfupds)
      pvttol = work(qpvttl) 
      sqfil2 = work(qsqfl2)
 
      call xdslof ( unsym,  msglvl, output, neqns,  nzla,   fnzlf,
     1              nsup2,  nsnind, mxnind, snmson, mstack,
     2              stkstr, mxtri,  mfront, mxtotf, maxzer, sqfil1,
     3              sqfil2, fpnlsz, fupdsz, pvttol,
     4              ncomp , nzcomp, lwksym, lwkifs, lwkofs )
 
      work ( qlwksy ) = lwksym
      work ( qlwkif ) = lwkifs
      work ( qlwkof ) = lwkofs
 
      needs  = lwksym
 
      work ( qstage ) = 20
      work ( qinuse ) = wkreqd
      work ( qneeds ) = needs
      work ( qmxuse ) = max ( maxwrk, int ( work ( qmxuse ) ) )
 
c     ----------------------------
c     ... store workspace pointers
c     ----------------------------
 
      work ( qxadj  ) = xadj
      work ( qadjnc ) = adjncy
      work ( qperm  ) = perm
      work ( qinvp  ) = invp
      work ( qxsup  ) = xsup2
      work ( qnassm ) = nassmb
      work ( qxlind ) = xlindx
 
      call xdslt2 ( t1, w1, t2, w2 )
      work ( qordtm ) = t2
      work ( qordwl ) = w2
 
      if ( msglvl .ge. 1 ) write ( output, 81000 ) t2, w2, maxwrk
 
      go to 9000
 
c ======================================================================
 
c     -----------
c     error traps
c     -----------
 
c     ---------------------------------
c     ... incorrect processing sequence
c     ---------------------------------
 
 8000 continue
      error = -200
      call hherr ( 3, 'xdslor', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88000 ) error, stage
      go to 9000
 
c     ---------------------------------------
c     ... insufficient storage for this stage
c     ---------------------------------------
 
 8100 continue
      error = -201
      call hherr ( 2, 'xdslor', error, wkreqd )
      if ( msglvl .gt. 0 ) write ( output, 88100 ) error, wkreqd, lwork
      needs = wkreqd
      go to 9000
 
c     ----------------------------
c     ... i/o error on file sqfil1
c     ----------------------------
 
 8200 continue
      error = -202
      call hherr ( 3, 'xdslor', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88200 ) error, sqfil1
      go to 9000
 
c---------------------------------------------------------------------
 
c     ------------------------
c     ... end of module xdslor
c     ------------------------
 
 9000 continue
      return
 
c---------------------------------------------------------------------
 
c     -----------
c     ... formats
c     -----------
 
68000 format ( /1x, '==============================='
     1         /1x, '= multifrontal ordering phase ='
     2         /1x, '===============================' )
 
68100 format ( /5x, 'number of super nodes                   = ', i15
     1         /5x, 'number of compressed subscripts         = ', i15 )
 
68150 format ( /5x, 'maximum number of sons                  = ', i15 )
 
68200 format ( /5x, 'new number of supernodes                = ', i15
     1         /5x, 'root of supernodal elimination tree     = ', i15 )
 
68300 format ( /5x, 'number of trees in the forest           = ', i15
     1         /5x, 'depth of supernode elimination forest   = ', i15
     2         /5x, 'maximum number of sons                  = ', i15
     3         /5x, 'root of supernode elimination forest    = ', i15 )
 
68400 format ( /5x, 'maximum nodes in a supernode            = ',i15
     1         /5x, 'maximum number of supernodal indices    = ',i15
     1         /5x, 'maximum number of nodes in a supernode  = ',i15
     2         /5x, 'total number of supernodal indices      = ',i15
     3         /5x, 'total number of supernodal nonzeroes    = ',f16.0
     4         /5x, 'maximum order of a frontal matrix       = ',i15
     4         /5x, 'maximum size of a frontal matrix        = ',i15
     5         /5x, 'maximum size of a stacked update matrix = ',i15
     6         /5x, 'maximum entries in the triangle         = ',i15
     7         /5x, 'maximum entries in the rectangle        = ',i15
     8         /5x, 'maximum entries in the trapezoid        = ',i15 )
 
68450 format ( /5x, 'total factor ops                        = ',d15.4
     1         /5x, 'total solve ops                         = ',d15.4 )
 
68500 format ( /5x, 'dense  factor ops                       = ',d15.4
     1         /5x, 'dense  factor vector ops                = ',d15.4
     2         /5x, 'sparse factor ops                       = ',d15.4
     3         /5x, 'sparse factor vector ops                = ',d15.4
     4         /5x, 'total  factor ops                       = ',d15.4 )
 
68600 format ( /5x, 'dense  solve ops                        = ',d15.4
     1         /5x, 'dense  solve vector ops                 = ',d15.4
     2         /5x, 'sparse solve ops                        = ',d15.4
     3         /5x, 'sparse solve vector ops                 = ',d15.4
     4         /5x, 'total  solve ops                        = ',d15.4 )
 
68700 format ( /5x, 'real stack storage                      = ', i15
     1         /5x, 'integer stack storage                   = ', i15 )
 
68800 format ( /5x, 'maximum number of entries in the stack  = ', i15
     1         /5x, 'maximum number of merges                = ', i15 )
 
81000 format ( /5x, 'cpu  time for reordering                = ', f15.6
     1         /5x, 'wall time for reordering                = ', f15.6
     2         /5x, 'max. amount of storage used             = ', i15  )
 
88000 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslor executed in an'
     2         /5x, 'incorrect sequence.  current stage = ', i10,
     3          5x, 'should be 10 or 20.' )
 
88100 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslor requires ', i15
     2         /5x, 'words of workspace and has only ', i15,
     3              ' available.' )
 
88200 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslor encountered i/o error'
     2         /5x, 'on i/o file no. ', i15 )
 
c---------------------------------------------------------------------
 
      end
