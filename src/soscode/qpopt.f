      subroutine   QPOPT     ( n, nclin, ldA, ldH, nHess,
     1                        A, bl, bu, cvec, H,
     2                        qpHess, istate, x,
     3                        prbtyp, cset, minsum, start, 
     4                        inform, obj, 
     5                        msglvl, iPrint, iSumm ,
     6                        Ax, clamda, iw, leniw, w, lenw,
     7                        mxfree, maxact, maxnZ ,
     8                        bigbnd, bigdx , tlcrsh, tolfea, tolOpt,
     9                        tolrnk, tollev,
     A                        fealim, optlim, kdegen, kcheck, 
     B                        fpiter, qpiter, nouter, nlvmod, delta , 
     C                        p1wkmn, p1wkmx, p1wkav,
     D                        p2wkmn, p2wkmx, p2wkav,
     E                        lminzh, lmaxzh, sminaa, smaxaa )
c     ==================================================================
c     ==================================================================
c     ====  QPOPT  /                                                ====
c     ====  qpopt -- top level driver for quadratic prog. solver    ====
c     ==================================================================
c     ==================================================================

      integer            n     , nclin , ldA   , ldH   , nHess , inform,
     1                   msglvl, iPrint, iSumm , leniw , lenw  , mxfree,
     2                   maxact, maxnZ , fealim, optlim, kdegen, kcheck,
     3                   fpiter, qpiter, nouter, nlvmod, p1wkmn, p1wkmx,
     4                   p2wkmn, p2wkmx 

      logical            cset, minsum
      
      double precision   obj   , bigbnd, bigdx , tlcrsh, tolfea, tolOpt,
     1                   tolrnk, tollev, p1wkav, p2wkav, delta , lminzh,
     2                   lmaxzh, sminaa, smaxaa
      
      integer            istate (n+nclin)
      
      integer            iw (leniw)
      
      double precision   A (ldA,*), Ax (*), bl (n+nclin), bu (n+nclin)
      
      double precision   clamda (n+nclin), cvec (*)
      
      double precision   H (ldH,*), x (n)
      
      double precision   w (lenw)

c                             << qpprnt >>
      external           qpHess, QPPRND

c     ==================================================================
c
c     derived from qpopt version 1.09
c     last modification -- 09-May-1996
c         This version of  QPOPT  dated 16-Feb-95.
c         Copyright  1988/1995  Optimates.
c
c     QPOPT  / QPOPT  solves problems of the form
c
c              minimize               f(x)
c                 x
c                                    (  x )
c              subject to    bl  .le.(    ).ge.  bu,
c                                    ( Ax )
c
c     where  '  denotes the transpose of a column vector,  x  denotes
c     the n-vector of parameters and  F(x) is one of the following...
c
c     FP            =              none    (find a feasible point)
c     LP            =    c'x
c     QP1           =          1/2 x'Hx     H n x n symmetric
c     QP2 (default) =    c'x + 1/2 x'Hx     H n x n symmetric
c     QP3           =          1/2 x'H'Hx   H m x n upper trapezoidal
c     QP4           =    c'x + 1/2 x'H'Hx   H m x n upper trapezoidal
c
c     Both  H and  G  are stored in the two-dimensional array  H  of
c     row dimension  ldH.  H  can be entered explicitly as the matrix
c     H,  or implicitly via a user-supplied version of the
c     subroutine qpHess.  If  ldH = 0,  H is not touched.
c
c     The vector  c  is entered in the one-dimensional array  cvec.
c
c     nclin  is the number of general linear constraints (rows of  A).
c     (nclin may be zero.)
c
c     The first  n  components of  bl  and   bu  are lower and upper
c     bounds on the variables.  The next  nclin  components are
c     lower and upper bounds on the general linear constraints.
c
c     The matrix  A  of coefficients in the general linear constraints
c     is entered as the two-dimensional array  A  (of dimension
c     ldA by n).  If nclin = 0, a is not referenced.
c
c     The vector  x  must contain an initial estimate of the solution,
c     and will contain the computed solution on output.
c
c     For more information, see:
c     User's Guide for QPOPT (Version 1.0) by:
c     P. E. Gill, W. Murray and M. A. Saunders.
c
c     Version 1.0-6  Jun 30, 1991. (Nag Mk 16 version). 
c     Version 1.0-7  Mar 21, 1993. Summary file added.
c     Version 1.0-8  Apr 10, 1994. Sum of infeas. added as an option.
c     Version 1.0-9  Jul 15, 1994. Debug output eliminated.
c
c     ==================================================================

      integer            ianrmj, idegen, isdel , iT    , itmax , itnfix, 
     1                   itnsav, j     , jdel  , jadd  , jinf  , jthcol,
     2                   ldQ   , ldR   , ldT   , lines1, lines2, litotl,
     3                   lwtotl, nctotl, miniw , minw  , nact1 , nactiv,
     4                   nartif, ncolT , ndegen, nerror, nfree , ngq   ,
     5                   nmoved, noprnt, nrejtd, nviol , numinf, nZ    ,
     6                   nZr   , phase , totitn  

      integer            xad   , xAnorm, xcq   , xd    , xfeatu, xgq   ,
     1                   xHx   , xkactv, xkx   , xq    , xr    ,
     2                   xrlam , xt    , xwrk  , xwtinf, xxT
      
      integer            nfix (2) 

      logical            inZtHZ, posdef

      double precision   alfa  , Amin  , Asize , condmx, dTmax , dTmin ,
     1                   epsmch, epspt5, epspt9, errmax, feamax, feamin,
     2                   fnlkpa, fnlkph, Hsize , rteps , rthuge, tolinc,
     3                   trulam, xnorm , tdelta, frac

      logical            done  , found , halted, header, 
     1                   named , prnt  , rowerr, Rset  ,
     2                   unitQ , vertex
                        
      character(len=6)   msg
      
      character(len=2)   prbtyp
      
      character(len=4)   start

      character(len=16)  names (1)

      double precision   zero, point3, half, point9, one, hundrd
      
      parameter        ( zero   =    0.0d0 ,
     1                   point3 =   33.0d-2, 
     2                   half   =    5.0d-1,
     3                   point9 =    9.0d-1,
     4                   one    =    1.0d0 ,
     5                   hundrd =  100.0d0 )

      double precision   hdmcon

      external           hdmcon
cjgl
      integer nreset

c     ==================================================================

      call CLKBEG (13)
      
cjgl
      nreset = 0
      phase  = 0
      
      delta  = zero

      lminzh = -one
      lmaxzh = -one

      sminaa = -one
      smaxaa = -one
      
c     Set the machine-dependent constants.

      epsmch = hdmcon (5)
      rteps  = sqrt (epsmch)
      rthuge = sqrt (hdmcon (2))

      epspt5 = rteps
      epspt9 = epsmch**point9

      named  = .false.

      fpiter = 0
      qpiter = 0
      totitn = 0
      nouter = 0  
      nlvmod = 0

      nctotl = n   + nclin

      p1wkmn = nctotl + 1
      p2wkmn = nctotl + 1
      p1wkmx = 0
      p2wkmx = 0
      p1wkav = zero
      p2wkav = zero

      header = .true.
      prnt   = .true.

C-->  condmx = max( one/epspt5, hundrd )
C-->  condmx = max( one/epspt3, hundrd )
      condmx = max( one/epspt5, hundrd )

c     Assign the dimensions of arrays in the parameter list of qpcore
c     (QPCORE).
      
c     Economies of storage are possible if the minimum number of active
c     constraints and the minimum number of fixed variables are known in
c     advance.  The expert user should alter minact and minfxd
c     accordingly.
c     If a linear program is being solved and the matrix of general
c     constraints has fewer rows than columns, i.e.,  nclin .lt. n,  a
c     non-zero value is known for minfxd.  In this case, vertex must
c     be set  .true. .

      vertex =  prbtyp .ne. 'QP'  .and.  nclin  .lt. n

      ldT    = max ( maxact, maxnZ )
      ncolT  = mxfree
      ldR    = ldT
      
      if  ( nclin .eq. 0 )  then
         ldQ = 1
      else
         ldQ = max( 1, mxfree )
      endif

c     ==================================================================
c     Cold start:  Only  x  is provided.
c     Warm start:  Initial working set is specified in  istate.
c     Hot  start:  The work arrays  iw  and  w  are assumed to have been
c                  initialized during a previous run.
c                  The first four components of  iw  contain details
c                  on the dimension of the initial working set.
c     ==================================================================

      
c     Allocate remaining work arrays.

      litotl = 3
      lwtotl = 0
      
c       << qploc >>
      call QPLOC  ( cset, n, nclin, litotl, lwtotl,
     1             ldT, ncolT, ldQ,
     2             miniw, xkactv, xkx,
     3             minw, xfeatu, xAnorm, xad, xHx, xd, xgq,
     4             xcq, xrlam, xR, xT, xQ, xwtinf, xwrk )

c     Check input parameters and storage limits.

c       << qpinit >>
      call QPINIT ( nerror, msglvl, iPrint, start,
     1              leniw, lenw, litotl, lwtotl,
     2              n, nclin,
     3              istate, named, names,
     4              bigbnd, bl, bu )

      if  ( nerror .gt. 0 )  then
         iw (1) = litotl
         iw (2) = lwtotl
         msg = 'errors'
         go to 800
      end if

c     ------------------------------------------------------------------
c     Define the initial feasibility tolerances in clamda.
c     ------------------------------------------------------------------
      
      if  ( tolfea .gt. zero )  then
         w(xfeatu:xfeatu+n+nclin-1) = tolfea
      endif

c       << cmdgen >>
      call CMDGEN ( 'Initialize anti-cycling variables', msglvl,
     1              n, nclin, nmoved, totitn, numinf,
     2              istate, bl, bu, clamda, w(xfeatu), x,
     3              iPrint, iSumm,
     4              tolinc, idegen, kdegen, ndegen,
     5              itnfix, nfix  )

      if  ( start .eq. 'cold'  .or.  start .eq. 'warm' )  then
         
c        ---------------------------------------------------------------
c        Cold or warm start.  Just about everything must be initialized.
c        The only exception is istate during a warm start.
c        ---------------------------------------------------------------

         iAnrmj = xAnorm
         do j = 1, nclin
            rnrm = zero
            do i=1,n
              rnrm = rnrm + A(j,i)**2
            enddo
            w(iAnrmj) = sqrt(rnrm)
            iAnrmj    = iAnrmj + 1
         enddo
         
         if  ( nclin .gt. 0)  then
c             << dcond >>
            call DCOND  ( nclin, w(xAnorm), 1, Asize, Amin )
         endif

c          << dcond >>
         call DCOND  ( nctotl, w(xfeatu), 1, feamax, feamin )
         frac = one/feamin
         w(xwtinf:xwtinf+nctotl-1) = frac*w(xfeatu:xfeatu+nctotl-1)

c        ---------------------------------------------------------------
c        Define the initial working set.
c               nfree ,  nactiv,  kactiv, kx,
c               istate (if start  = 'cold')
c               nartif  ( if vertex = 'true')
c        ---------------------------------------------------------------
         
c          << cmcrsh >>
         call CMCRSH ( start , vertex,
     1                 nclin , nctotl, nactiv, nartif,
     2                 nfree , n     , ldA,
     3                 istate, iw(xkactv), iw(xkx),
     4                 bigbnd, tlcrsh,
     5                 A, Ax, bl, bu, clamda, x, w(xgq), w(xwrk) )

c        ---------------------------------------------------------------
c        Compute the TQ factorization of the working-set matrix.
c        ---------------------------------------------------------------
         
         unitQ = .true.
         nZ    = nfree
         
         if  ( nactiv .gt. 0 )  then
            
            iT     = nactiv + 1
            nact1  = nactiv
            nactiv = 0
            ngq    = 0

c            << rzadds >>
           call RZADDS ( unitQ     , vertex    , 1         , nact1     ,
     1                   iT        , nactiv    , nartif    , nZ        ,
     2                   nfree     , nrejtd    , ngq       , n         ,
     3                   ldQ       , ldA       , ldT       , istate    ,
     4                   iw(xkactv), iw(xkx)   , condmx    , A         ,
     5                   w(xT)     , w(xgq)    , w(xQ)     , w(xwrk)   ,
     6                   w(xd)     , w(xrlam)  , rthuge    , iPrint    ,
     7                   epspt9    , Asize     , dTmax     , dTmin )

        end if

      else
     1if  ( start .eq. 'hot ' )  then
         
c        ---------------------------------------------------------------
c        Arrays  iw  and  w  have been defined in a previous run.
c        The first three elements of  iw  are  unitQ,  nfree and nactiv.
c        ---------------------------------------------------------------

         unitQ  = iw(1) .eq. 1
         nfree  = iw(2)
         nactiv = iw(3)

         nZ     = nfree - nactiv
         
      end if

c     Install the transformed linear term in cq.

      if  ( cset )  then
        w(xcq:xcq+n-1) = cvec(1:n)
c          << cmqmul( >>
         call CMQMUL( 6, n, nZ, nfree, ldQ, unitQ,
     1                iw(xkx), w(xcq), w(xQ), w(xwrk) )
      end if   

      Rset   = .false.
      if  ( prbtyp .eq. 'LP' )  then
         itmax = max ( fealim, optlim )
      else
         itmax = fealim
      end if

      jinf   =  0
      call CLKSUM (13)

C+    When minimizing the sum of infeasibilities,
C+    nZr    =  nZ  implies steepest-descent in the two-norm.
C+    nZr    =  0   implies steepest-descent in the infinity norm.
      nZr    =  0

c     ==================================================================
c     repeat               (until working set residuals are acceptable)
c        -----------------------------------------------
c        Move x onto the constraints in the working set.
c     ==================================================================

 300     continue
         phase  = 1
         nouter = nouter + 1

         call CLKBEG (14)
         
c          << cmsetx >>
         call CMSETX ( rowerr, unitQ, nclin,
     1                 nactiv, nfree, nZ,
     2                 n, ldQ, ldA,
     3                 ldT, istate, iw(xkactv),
     4                 iw(xkx), errmax,
     5                 xnorm, A, Ax,
     6                 bl, bu, clamda,
     7                 w(xT), x, w(xQ),
     8                 w(xd), w(xwrk) )

         if  ( rowerr )  then
            msg    = 'rowerr'
            numinf = 1
            call CLKSUM (14)
            go to 800
         end if

         itnsav = totitn
         
c          << lpcore >>
         call LPCORE ( prbtyp, msg, cset, Rset, unitQ,
     3                 totitn, itmax, jinf, nviol, n, nclin, ldA,
     5                 nactiv, nfree, nZr, nZ, istate, iw(xkactv),
     6                 iw(xkx), QPPRND, obj, numinf, xnorm, A, Ax,
     9                 bl, bu, cvec, clamda, w(xfeatu), x, 
     C                 msglvl, iPrint, iSumm , lines1, lines2,
     D                 header, prnt, ldT   , ldQ, ldR,
     F                 tolOpt, epspt9, Asize , dTmax , dTmin,
     H                 tolinc, idegen, kdegen, ndegen, itnfix, nfix  ,
     J                 alfa  , trulam, isdel , jdel  , jadd,
     K                 kcheck, minsum, bigbnd, bigdx ,
     M                 w(xAnorm), w(xAd), w(xd), w(xgq), w(xcq),
     N                 w(xrlam), w(xR), w(xT), w(xQ), w(xwtinf),
     O                 w(xwrk), p1wkmn, p1wkmx, p1wkav )

         fpiter = fpiter + (totitn - itnsav)

         call CLKSUM (14)
         
         if  ( prbtyp .eq. 'QP'  .and.  msg .eq. 'feasbl')  then
            
            if  ( msglvl .gt. 0 )  then
               if  ( iPrint .gt. 0) write (iPrint, 1000) totitn
               if  ( iSumm  .gt. 0) write (iSumm , 1000) totitn
            end if
            
            phase = 2

            Rset   = .true.
            itmax  = totitn + optlim

c           ------------------------------------------------------------
c           Compute the first QP objective and transformed gradient.
c           ------------------------------------------------------------

            inZtHZ = .true.
            
 400        continue

            call CLKBEG (15)

            obj = zero
            if  ( cset )  then
               do i=1,n
                 obj = obj + cvec(i)*x(i)
               enddo
            end if

            jthcol = 0

            call qpHess ( n, ldH, nHess, jthcol, H, x, delta, w(xHx) )
            
            obj    = obj + half*dot_product(w(xHx:xHx+n-1),x(1:n))
            
            w(xgq:xgq+n-1) = w(xHx:xHx+n-1)
            
c             << cmqmul >>
            call CMQMUL ( 6, n, nZ, nfree, ldQ, unitQ,
     1                    iw(xkx), w(xgq), w(xQ), w(xwrk) )

            if  ( cset)  then
              do j=0,n-1
                w(xgq+j) = w(xgq+j) + w(xcq+j)
              enddo
            endif

c           ----------------------------------------------
c           ... generate the initial reduced Hessian  Z'HZ
c           ----------------------------------------------

            if  ( inZtHZ )  then
c                << qpzthz >>
               call QPZTHZ ( unitQ, qpHess, delta, 
     1                       n, nZ, nfree, nHess,
     2                       ldQ, ldH, ldR,
     3                       iw(xkx), Hsize,
     4                       H, w(xR), w(xQ),
     5                       w(xwrk), w(xrlam) )
            endif
            
c           -------------------------------------------------
c           ... REPEAT until a positive definite QP is solved
c               or until we break down ....
c           -------------------------------------------------
            
            
c           ------------------------------------------------------------
c           Find the Cholesky factor R of an initial reduced Hessian.
c           The magnitudes of the diagonals of  R  are nonincreasing.
c           ------------------------------------------------------------

c              << qpchl0 >>
             call QPCHL0 ( nZ    , nZr   , nZ    , ldR   , Hsize ,
     1                    tolrnk, w(xR) , posdef )
            
            if ( posdef )  then

c              ----------------------------------------------------
c              ... we have a feasible point and a positive definite
c                  initial reduced Hessian, so proceed to find
c                  at least a local stationary point
c              ----------------------------------------------------

               itnsav = totitn
               
c                << qpcore >>
               call QPCORE ( prbtyp, msg, cset, unitQ, totitn, 
     3                 itmax, nviol, n, nclin, nHess, ldA, ldH, 
     5                 nactiv, nfree, nZr, nZ, istate, iw(xkactv),
     6                 iw(xkx), qpHess, delta, QPPRND, obj, 
     8                 xnorm, Hsize, A, Ax, bl, bu,
     A                 cvec, clamda, w(xfeatu), H, x, 
     C                 msglvl, iPrint, iSumm , lines1, lines2,
     D                 header, prnt, ldT , ldQ, ldR, 
     F                 tolOpt, epspt9, Asize , dTmax , dTmin,
     H                 tolinc, idegen, kdegen, ndegen, itnfix, 
     J                 nfix ,alfa  , trulam, isdel , jdel, jadd,
     K                 kcheck, maxnZ, bigbnd, bigdx , tolrnk,
     M                 w(xAnorm), w(xAd), w(xHx), w(xd), w(xgq),
     N                 w(xcq), w(xrlam), w(xR), w(xT), w(xQ),
     O                 w(xwtinf), w(xwrk), p2wkmn, p2wkmx, p2wkav)

               qpiter = qpiter + (totitn - itnsav)
               
               posdef = msg .ne. 'indef '

            end if

            call CLKSUM (15)

            if  ( .not. posdef ) then

c              ----------------------------------------------
c              ... generate the current reduced Hessian  Z'HZ
c                  (we have its factorization, but not the
c                   matrix itself.)
c              ----------------------------------------------

               call CLKBEG (16)
               msg = '      '

c                 << qpzthz >>
                call QPZTHZ ( unitQ, qpHess, delta, 
     1                       n, nZ, nfree, nHess,
     2                       ldQ, ldH, ldR,
     3                       iw(xkx), Hsize,
     4                       H, w(xR), w(xQ),
     5                       w(xHx), w(xwrk) )

c                 << qpmdlv >>
                call QPMDLV ( delta , nZ    , ldR   , epsmch  , tollev, 
     1                       lminzh, lmaxzh, w(xR) , w(xrlam), w(xd) ,
     2                       w(xHx), msglvl, iPrint )

               nlvmod = nlvmod + 1
               inZtHZ = .false.

               call CLKSUM (16)
               go to 400

            endif

         endif

         found  = msg .eq. 'optiml'  .or.
     1            msg .eq. 'feasbl'  .or.
     2            msg .eq. 'deadpt'  .or.
     3            msg .eq. 'weak  '  .or.
     4            msg .eq. 'unbndd'  .or.
     5            msg .eq. 'infeas'
            
         halted = msg .eq. 'itnlim'  .or.
     1            msg .eq. 'Rz2big'  .or.
     2            msg .eq. 'indef '

         if  ( found )  then
c             << cmdgen >>
            call CMDGEN ( 'Optimal', msglvl,
     1                    n, nclin, nmoved, totitn, numinf,
     2                    istate, bl, bu, clamda, w(xfeatu), x,
     3                    iPrint, iSumm,
     4                    tolinc, idegen, kdegen, ndegen,
     5                    itnfix, nfix  )
         end if

         done   = found          .and.
     1            nviol  .eq. 0  .and.
     2            nmoved .eq. 0

c        =============================
c        until      done  .or.  halted
c        =============================
         
         if  ( .not. (done  .or.  halted))  then

            if  ( msg .ne. 'resetx' .and. msg .ne. 'failed' )  then

               write (iPrint, *) 'internal program check: msg =', msg
               msg = 'failed'

            else
cjgl           go to 300
               if(msg.eq.'resetx') nreset = nreset + 1
               if(nreset.lt.2) go to 300
               
            endif

         endif
c
c     --------------------------------------------------
c     Set   clamda.  Print the full solution.
c     Clean up.  Save values for a subsequent hot start.
c     --------------------------------------------------
         
      call CLKBEG (17)

c       << cmwrap >>
      call CMWRAP ( nfree, ldA,
     1              n, nclin, nctotl,
     2              nactiv, istate, iw(xkactv), iw(xkx),
     3              A, bl, bu, x, clamda, w(xfeatu),
     4              w(xwrk), w(xrlam), x )
      
c       << cmprnt >>
      call CMPRNT ( msglvl, iprint, n, nclin, nctotl, bigbnd,
     1             named, names, istate,
     2             bl, bu, clamda, w(xfeatu), w(xwrk) )

      iw(1) = 0
      if  ( unitQ) iw(1) = 1
      iw(2) = nfree
      iw(3) = nactiv

c     --------------------------------------------------
c     ... Set return flag and print concluding messages. 
c     --------------------------------------------------
      
c       << qperpt >>
 800  continue
      call QPERPT ( msglvl, iPrint, msg   , inform, prbtyp,
     1              nZ    , nZr   , maxnZ , nerror, numinf,
     2              obj   , minsum, errmax )

      call CLKSUM (17)
      
c     ------------------------------------------------------
c     ... compute condition number for final reduced Hessian
c     ------------------------------------------------------

      call CLKBEG (18)
         
      if  ( phase .eq. 2 )  then

         if  ( nZr .gt. 0 )  then

            tdelta = delta
            
c             << qpzthz >>
            call QPZTHZ ( unitQ, qpHess, tdelta, 
     1                    n, nZr, nfree, nHess,
     2                    ldQ, ldH, ldR,
     3                    iw(xkx), Hsize,
     4                    H, w(xR), w(xQ),
     5                    w(xHx), w(xwrk) )
         
            noprnt = 0
                     
c             << qpmdlv >>
            call QPMDLV ( tdelta, nZr   , ldR   , epsmch   , tollev,
     1                    lminzh, lmaxzh, w(xR) , w(xwtinf), w(xd) ,
     3                    w(xHx), noprnt, noprnt )

            if  ( lminzh .gt. 0 )  then
               fnlkph = lmaxzh / lminzh
            else
               fnlkph = -one
            endif

         else

            fnlkph = zero
            lminzh = zero
            lmaxzh = zero
            
         endif

      else

         fnlkph = -one
         lminzh = -one
         lmaxzh = -one

      endif

c     ---------------------------------------------------------------
c     ... compute condition number for final active constraint matrix
c     ---------------------------------------------------------------

      if  ( phase .eq. 1  .or.  phase .eq. 2 )  then

         if  ( nactiv .gt. 0 )  then

            if  ( unitQ .or. nactiv .eq. 1 )  then
               
               sminaa = one
               smaxaa = one

            else

               xxT = xT + nZ * ldT
c                << qpcnda >>
               call QPCNDA ( nactiv, ldT   , w(xxT) , sminaa, smaxaa,
     1                       w(xd) , w(xwtinf), 2*(n+nclin)  )

            endif
            
            if  ( sminaa .gt. 0 )  then
               fnlkpa = smaxaa / sminaa
            else
               fnlkpa = -one
            endif

         else

            fnlkpa = zero
            sminaa = zero
            smaxaa = zero
            
         endif

      else

         fnlkpa = -one
         sminaa = -one
         smaxaa = -one

      endif

      call CLKSUM (18)
         
c     ------------------
c     ... closing output
c     ------------------

      if  ( fpiter .gt. 0 )  then
         p1wkav = p1wkav / fpiter
      else
         p1wkmn = 0
      endif

      if  ( qpiter .gt. 0 )  then
         p2wkav = p2wkav / qpiter
      else
         p2wkmn = 0
      endif

      if  ( inform .ne. 6 )  then
         if  ( msglvl .gt. 0 )  then
            
            if  ( iPrint .gt. 0)  then
               write (iPrint, 2000) prbtyp, fpiter, qpiter, inform
            endif
            if  ( iSumm  .gt. 0)  then
               write (iSumm , 2000) prbtyp, fpiter, qpiter, inform
            endif
            
         end if
      end if

      if  ( msglvl .ge. 20 )  then
         
         if  ( lminzh .gt. zero )  then
            fnlkph = lmaxzh / lminzh
         else
            fnlkph = -one
         endif
         
         if  ( sminaa .gt. zero )  then
            fnlkpa = smaxaa / sminaa
         else
            fnlkpa = -one
         endif
         
         write (iPrint, 4000 ) inform, obj   , fpiter, qpiter, 
     1                         nouter, nlvmod, delta , 
     2                         p1wkmn, p1wkmx, p1wkav,
     3                         p2wkmn, p2wkmx, p2wkav,
     4                         fnlkph, fnlkpa

      endif
      
      return
       
 1000 format(  ' Itn', i6, ' -- Feasible point found.' )

 2000 format(/ ' Exit from ', a2, ' problem after ', 
     1         i5, ' Phase 1 iterations,',
     2         i5, ' Phase 2 iterations.',
     3         '  Inform =', i3 )
      
 4000 format(// ' Final results from solution of QP',
     1       /  ' ---------------------------------', 
     2       /  '                           return code:', i8, 
     3       /  '     final value of objective function:', g16.7,
     5       /  '              Phase 1 iterations taken:', i8,
     6       /  '              Phase 2 iterations taken:', i8,
     7       /  '                Outer iterations taken:', i8,
     8       /  '     number of Levenberg modifications:', i8,
     9       /  ' accumulated QP Levenberg modification:', 1pg16.7,
     B       //  ' working set size, Phase I:',
     C       /   '                              minimum:', i8, 
     D       /   '                              maximum:', i8, 
     E       /   '                              average:', g16.7,
     F       //  ' working set size, Phase II:',
     G       /   '                              minimum:', i8, 
     H       /   '                              maximum:', i8, 
     I       /   '                              average:', g16.7,
     K       //  ' condition of final reduced Hessian:', g16.7,
     L       /   ' condition of final active Jacobian:', g16.7 )

c     end of QPOPT  (qpopt)
      
      end
