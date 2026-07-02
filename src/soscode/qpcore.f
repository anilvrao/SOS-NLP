      subroutine   QPCORE   ( prbtyp, msg, cset, unitQ,
     3                        iter, itmax, nviol, n, nclin,
     4                        nHess, ldA, ldH, nactiv, nfree,
     5                        nZr, nZ, istate, kactiv, kx,
     7                        qpHess, delta, itprnt, objqp, xnorm,
     8                        Hsize, A, Ax, bl, bu, cvec, featol,
     A                        featlu, H, x, msglvl, iPrint, iSumm ,
     C                        lines1, lines2, header, prnt, ldT ,
     E                        ldQ, ldR, tolOpt, epspt9, Asize ,
     G                        dTmax , dTmin, tolinc, idegen, kdegen,
     H                        ndegen, itnfix, nfix  ,alfa  , trulam,
     J                        isdel , jdel, jadd, kcheck, maxnZ,
     L                        bigbnd, bigdx , tolrnk, Anorm, Ad, Hx,
     M                        d, gq, cq, rlam, R, T, Q, wtinf, wrk,
     O                        p2wkmn, p2wkmx, p2wkav )

c     ==================================================================
c     ==================================================================
c     ====  QPCORE /                                                ====
c     ====  qpcore -- solve quadratic programming problem           ====
c     ==================================================================
c     ==================================================================

      integer            iter  , itmax , nviol , n  , nclin , nHess,
     1                   ldA, ldH,
     2                   nactiv, nfree , nZr   , nZ    ,
     3                   iPrint, iSumm , lines1, lines2,
     4                   ldT   , ldQ   , ldR   , 
     5                   idegen, kdegen, ndegen, itnfix,
     6                   isdel , jdel  , jadd  ,
     7                   kcheck, maxnZ , msglvl,
     8                   p2wkmn, p2wkmx

      double precision   objqp , xnorm , Hsize, tolOpt, epspt9,
     1                   Asize , dTmax , dTmin,
     2                   tolinc,
     3                   alfa  , trulam,
     4                   bigbnd, bigdx , tolrnk, delta , p2wkav
      
      logical            cset, unitQ,
     1                   header, prnt

      character(len=2)   prbtyp
      
      character(len=6)   msg
      
      integer            istate (n+nclin), kactiv (n), kx (n), nfix (2)
      
      double precision   A (ldA,*), Ax (*), bl (n+nclin), bu (n+nclin),
     1                   cvec (*), featol (n+nclin), featlu (n+nclin),
     2                   H (ldH,*), x (n),
     3                   Anorm (*), Ad (*), Hx (*), d (n), gq (*),
     4                   cq  (*), rlam (*), R (ldR, *), T (ldT, *),
     5                   Q (ldQ, *), wtinf (*), wrk (*)
      
      external           qpHess, itprnt

c     ==================================================================
c
c     derived from qpopt version 1.0
c     last modification -- 24-May-1996
c
c          This version of  qpcore  dated  16-Jan-95.
c          Copyright  1988--1995  Optimates.
c
c     qpcore (QPCORE) is a subroutine for general quadratic programming.
c     On entry, it is assumed that an initial working set of
c     linear constraints and bounds is available.
c     The arrays  istate, kactiv and kx  will have been set accordingly
c     and the arrays  T  and  Q  will contain the TQ factorization of
c     the matrix whose rows are the gradients of the active linear
c     constraints with the columns corresponding to the active bounds
c     removed.  The TQ factorization of the resulting (nactiv by nfree)
c     matrix is  A(free)*Q = (0 T),  where Q is (nfree by nfree) and T
c     is upper-triangular.
c
c     Over a cycle of iterations, the feasibility tolerance featol
c     increases slightly (from tolx0 to tolx1 in steps of tolinc).
c     this ensures that all steps taken will be positive.
c
c     After idegen consecutive iterations, variables within featol of
c     their bounds are set exactly on their bounds and iterative
c     refinement is used to satisfy the constraints in the working set.
c     featol is then reduced to tolx0 for the next cycle of iterations.
c
c     Values of istate(j) for the linear constraints.......
c
c     Istate(j)
c     ---------
c          0    constraint j is not in the working set.
c          1    constraint j is in the working set at its lower bound.
c          2    constraint j is in the working set at its upper bound.
c          3    constraint j is in the working set as an equality.
c
c     Constraint j may be violated by as much as featol(j).
c
c     This version of  qpcore  dated  16-Jan-95.
c     Copyright  1988--1995  Optimates.
c     ==================================================================

      integer            irefn , iadd  , ifix  , inform, issave, 
     1                   it    , jbigst, jdsave, jinf  , jmax  ,
     2                   jsmlst, jthcol, jtiny , kbigst, kdel  ,
     3                   ksmlst, nctotl, ngq   , nmoved, notOpt,
     4                   numinf 
      
      logical            deadpt, delreg, firstv, giveup, hitcon,
     1                   hitlow, minmzr, move  , onbnd , overfl,
     2                   posdef, renewr, Rset  , singlr, statpt,
     3                   unbndd, uncon , unitgZ

      logical            addrun

      integer            runlen

      double precision   alfap , alfhit, bigalf, biggst, condmx,
     1                   condRz, condt , dinky , dnorm , dRzmax,
     2                   dRzmin, Dzz   , errmax, flmax , gfnorm,
     3                   gzdz  , gznorm, gzRnrm, objchg, objsiz,
     4                   smllst, suminf, tinyst, trubig, trusml,
     5                   wssize, zerolm

      double precision   DDIVFN

      external           DDIVFN
      
      character(len=6)   empty
      parameter         (empty = '      ')

      double precision   zero, half, one
      parameter        ( zero = 0.0d0, half = 5.0d-1, one  = 1.0d0 )

      integer            mrefn
      parameter        ( mrefn = 1 )

c     Specify the machine-dependent parameters.

      double precision   hdmcon

      external           hdmcon

c     ==================================================================

      flmax  = hdmcon (2)

      if (cset )  then
         ngq = 2
      else
         ngq = 1
      end if

      iT     = 1


c     Initialize.

      irefn  =   0
      jinf   =   0
      nctotl =   n + nclin
      nviol  =   0
      numinf =   0
      suminf =   zero
      condmx =   flmax

      delreg = .false.
      firstv = .false.
      posdef = .true.
      renewr = .false.
      Rset   = .true.
      singlr = .false.
      uncon  = .false.
      unitgZ = .false.

      notOpt = 0
      Dzz    = one

      addrun = .true.
      runlen = 0

      msg    = empty

c*    ======================Start of main loop==========================
c+    do while (msg .eq. empty)
c*    ======================Start of main loop==========================

 100  continue
      if  (msg .eq. empty )  then


         gznorm = zero
         do i=1,nZ
           gznorm = gznorm + gq(i)**2
         enddo
         gznorm = sqrt(gznorm)

         if  ( nZr .eq. nZ )  then
            gZrnrm = gznorm
         else
            gZrnrm = zero
            do i=1,nZr
              gZrnrm = gZrnrm + gq(i)**2
            enddo
            gZrnrm = sqrt(gZrnrm)
         end if

         gfnorm = gznorm
         if  ( nfree .gt. 0  .and.  nactiv .gt. 0 )  then
            gfnorm = zero
            do i=1,nfree
              gfnorm = gfnorm + gq(i)**2
            enddo
            gfnorm = sqrt(gfnorm)
         endif

         objsiz = one + abs( objqp )
         wssize = zero
         if  ( nactiv .gt. 0)  then
             wssize = dTmax
         endif
     
         dinky  = tolOpt*max( wssize, objsiz, gfnorm )
         if  ( uncon )  then
            unitgZ = gZrnrm .le. dinky
         end if

c        If the reduced gradient (Zr)'g is small and Hz is positive
c        definite,  x is a minimizer on the working set.
c        A maximum number of unconstrained steps is imposed to
c        allow for  dinky  being too large because of bad scaling.

         statpt = gZrnrm .le. dinky
         giveup = irefn  .gt. mrefn

         minmzr = statpt  .and.  posdef
         deadpt = statpt  .and.  singlr

c        ---------------------------------------------------------------
c        Print the details of this iteration.
c        ---------------------------------------------------------------
c        Define small quantities that reflect the size of x, R and
c        the constraints in the working set.

         if  ( prnt )  then
            if  ( nZr .gt. 0 )  then
c                << dcond >>
               call DCOND  ( nZr, R, ldR+1, dRzmax, dRzmin )
               condRz = DDIVFN  ( dRzmax, dRzmin, overfl )
            else
               condRz = one
            end if

            if  ( nactiv .gt. 0 )  then
               condT  = DDIVFN  ( dTmax, dTmin, overfl )
            else
               condT  = one
            end if

c             << itprnt >>
            call itprnt ( prbtyp, header, Rset,
     1                    msglvl, iter,
     2                    isdel, jdel, jadd,
     3                    n, nclin, nactiv,
     4                    nfree, nZ, nZr,
     5                    ldR, ldT, istate,
     6                    iPrint, iSumm , lines1, lines2,
     7                    alfa, condRz, condT,
     8                    Dzz, gZrnrm,
     9                    numinf, suminf, notOpt, objqp, trulam,
     A                    Ax, R, T, x,
     B                    wrk )
         end if

         if  ( minmzr  .or.  giveup )  then

            call CLKBEG (19)
            
c           ============================================================
c           The point  x  is a constrained stationary point.
c           Compute Lagrange multipliers.
c           ============================================================
c           Define what we mean by ``non-optimal'' multipliers.

            notOpt = 0
            jdel   = 0
            zerolm = dinky
            smllst = dinky
            biggst = dinky + one
            tinyst = dinky

c             << cmmul1 >>
            call CMMUL1 ( prbtyp, msglvl,
     1                    n     , ldA   , ldT   ,
     2                    nactiv, nfree , nZ    ,
     3                    istate, kactiv, kx    ,
     4                    zerolm, notOpt, numinf,
     5                    trusml, smllst, jsmlst, ksmlst,
     6                            tinyst, jtiny , jinf  ,
     7                    trubig, biggst, jbigst, kbigst,
     8                    A     , Anorm , gq    ,
     9                    rlam  , T     , wtinf ,
     A                    iPrint ) 

            if  ( nZr .lt. nZ )  then
c                << cmmul2 >>
               call CMMUL2 ( msglvl, n, nZr, nZ,
     1                       zerolm, notOpt, numinf,
     2                       trusml, smllst, jsmlst,
     3                       tinyst, jtiny , gq,
     4                       iPrint )
            end if

            if  ( notOpt .eq. 0  .and.  posdef )  then
               
               msg    = 'optiml'
               call CLKSUM (19)
               go to 100
               
            end if

c           ------------------------------------------------------------
c           Delete one of three types of constraint
c           (1) regular           jsmlst > 0   istate(jsmlst) = 1, 2
c           (2) temporary bound   jsmlst > 0,  istate(jsmlst) = 4 
c           (3) artificial        jsmlst < 0  
c           ------------------------------------------------------------
            trulam = trusml
            jdel   = jsmlst
            delreg = .false.

            if  ( nZr+1 .gt. maxnZ  .and.  jdel .ne. 0 )  then 
               msg    = 'Rz2big'
               call CLKSUM (19)
               go to 100
            end if

            if  ( jdel .gt. 0 )  then

c              Regular constraint or temporary bound.
c              delreg  says that a regular constraint was deleted. 
c              jdsave, issave are only defined if  delreg  is true.

               kdel         = ksmlst
               isdel        = istate(jdel)
               istate(jdel) = 0
               delreg       = isdel .ne. 4
               if  ( delreg )  then
                  jdsave = jdel
                  issave = isdel
               end if
            end if

c           Update the factorizations.

c             << rzdel >>
            call RZDEL  ( unitQ, iT,
     1                   n, nactiv, nfree, ngq, nZ, nZr,
     2                   ldA, ldQ, ldT,
     3                   jdel, kdel, kactiv, kx,
     4                   A, T, gq, Q, d, rlam,
     5                   dTmax, dTmin )

            renewr = .true.
c             << qpcolr >>
            call QPCOLR ( singlr, posdef, renewr, unitQ,
     1                    n, nZr, nfree, nHess, ldQ, ldH, ldR,
     2                    kx, Hsize, Dzz, tolrnk,
     3                    qpHess, delta, H, R, Q,
     4                    wrk, d )
            
            call CLKSUM (19)

            if  ( .not. posdef )  then
               msg = 'indef '
               go to 100
            endif
            
            irefn  =  0
            prnt   = .false.
            uncon  = .false.
            
         else
            
c           ============================================================
c           Compute a search direction.
c           ============================================================

            if  ( .not. posdef )  then
               write (*,*) 'error in QPCORE:'
               write (*,*) 'posdef should be true and isn''t'
            endif
            
            if  ( iter .ge. itmax )  then 
               msg = 'itnlim'
               go to 100
            end if

            call CLKBEG (20)
            
            prnt  = .true.
            iter  = iter  + 1

            p2wkmn = min ( nZr, p2wkmn )
            p2wkmx = max ( nZr, p2wkmx )
            p2wkav = p2wkav + nZr

c           << qpgetd >>
            call CLKBEG (21)
            call QPGETD ( delreg, posdef, statpt, unitgZ, unitQ,
     1                    n, nclin, nfree,
     2                    ldA, ldQ, ldR, nZr,
     3                    issave, jdsave,
     4                    kx, dnorm, gzdz,
     5                    A     , Ad    , d     ,
     6                    gq    , R     , Q     , wrk    )
            call CLKSUM (21)

c           ------------------------------------------------------------
c           Find the constraint we bump into along  d.
c           Update  x  and  Ax  if the step  alfa  is nonZero.
c           ------------------------------------------------------------
c           qpchzr initializes  alfhit  to bigalf. If it is still
c           that value on exit,  it is regarded as infinite.
c                 << ddiv >>
            bigalf = DDIVFN  ( bigdx, dnorm, overfl )

c           << qpchzr >>
            call CLKBEG (22)
            call QPCHZR ( firstv, n, nclin,
     1                    istate, bigalf, bigbnd, dnorm,
     2                    hitlow, move, onbnd, unbndd,
     3                    alfhit, alfap, jadd,
     4                    Anorm , Ad    , Ax,
     5                    bl, bu, featol, featlu, d, x,
     6                    epspt9, tolinc, ndegen )
            call CLKSUM (22)

c           ------------------------------------------------------------
c           If Hz is positive definite,  alfa = 1.0  will be the step
c           to the minimizer of the quadratic on the current working
c           set.  If the unit step does not violate the nearest
c           constraint by more than featol,  the constraint is not
c           added to the working set.
c           ------------------------------------------------------------

            uncon  = alfap .gt. one  .and.  posdef
            hitcon = .not. uncon

            if  ( hitcon )  then
               alfa  = alfhit
               irefn = 0
            else
               irefn  = irefn + 1
               jadd   = 0
               alfa   = one
            end if

           if  ( hitcon  .and.  unbndd )  then 
               msg = 'unbndd'
               call CLKSUM (20)
               go to 100
            end if

c           Predict the change in the QP objective function.

            if  ( posdef )  then
                objchg = alfa*gzdz*(one - half*alfa)
            else
                objchg = alfa*gzdz + half*alfa**2*Dzz
            end if

c           Check for a dead point or unbounded solution.

            if  ( objchg .ge. - epspt9*objsiz  .and.  deadpt )  then
               msg = 'deadpt'
               call CLKSUM (20)
               go to 100
            end if

            if  ( objchg .ge.   epspt9*objsiz )  then
               msg = 'resetx'
               call CLKSUM (20)
               go to 100
            end if

            xnorm = zero
            do j=1,n
              x(j) = x(j) + alfa*d(j)
              xnorm = xnorm + x(j)**2
            enddo
            xnorm = sqrt(xnorm)
            do j=1,nclin
              Ax(j) = Ax(j) + alfa*Ad(j)
            enddo
            
            if  ( hitcon )  then
               
c              -----------------------------------------
c              Add a constraint to the working set.
c              Update the TQ factors of the working set.
c              Use  d  as temporary work space.
c              -----------------------------------------

               if  ( bl(jadd) .eq. bu(jadd) )  then
                  istate(jadd) = 3
               else
     1         if  ( hitlow )  then
                  istate(jadd) = 1
               else
                  istate(jadd) = 2
               end if

               if  ( jadd .gt. n )  then
                  iadd = jadd - n
               else
                  if  ( hitlow )  then
                     x(jadd) = bl(jadd)
                  else
                     x(jadd) = bu(jadd)
                  end if

                  ifixloop: do ifix = 1, nfree
                     if  ( kx(ifix) .eq. jadd) exit ifixloop
                  enddo ifixloop
               end if

c                << rzadd >>

               call CLKBEG (23)
               call RZADD  ( unitQ, Rset,
     1                      inform, ifix, iadd, jadd, it,
     2                      nactiv, nZ, nfree, nZr, ngq,
     3                      n, ldA, ldQ, ldR, ldT,
     4                      kx, condmx, Dzz,
     5                      A, R, T, gq, Q,
     6                      wrk, rlam, d,
     7                      iPrint, epspt9, Asize, dTmax, dTmin )
               call CLKSUM (23)
              
               nZr    = nZr - 1
               nZ     = nZ  - 1

               if  ( jadd .le. n )  then

c                 A simple bound has been added.

                  nfree  = nfree  - 1
               else

c                 A general constraint has been added.

                  nactiv          = nactiv + 1
                  kactiv (nactiv) = iadd
                  
               end if

c              ---------------------------------------------------------
c              Check if  Hz  has become positive definite.
c              Recompute the last column of Rz if unacceptable
c              growth has occurred.
c              --------------------------------------------------------

               if  ( .not. posdef )  then
                  
c                   << qpcolr >>

                  write (*,*) 'this should be dead code around qpcolr'
                  call CLKBEG (26)
                  call QPCOLR ( singlr, posdef, renewr, unitQ,
     1                          n, nZr, nfree, nHess, ldQ, ldH, ldR,
     2                          kx, Hsize, Dzz, tolrnk,
     3                          qpHess, delta, H, R, Q,
     4                          wrk, d )
                  call CLKSUM (26)
                  if  ( .not. posdef )  then
                     msg = 'indef '
                     call CLKSUM (20)
                     go to 100
                  endif

               endif
            
            endif
               
c           Increment featol.

            do j=1,nctotl
              featol(j) = featol(j) + tolinc*featlu(j)
            enddo

            if  ( mod ( iter, kcheck ) .eq. 0 )  then
               
c              ---------------------------------------------------------
c              Check the feasibility of constraints with non-negative 
c              istate  values.  If violations have occurred,  force 
c              iterative refinement and a switch to phase 1.
c              ---------------------------------------------------------

c                << cmfeas >>

               call CLKBEG (24)
               call CMFEAS ( n, nclin, istate,
     1                       bigbnd, nviol, jmax, errmax,
     2                       Ax, bl, bu, featol, x )
               call CLKSUM (24)


               if  ( nviol .gt. 0 )  then
                  if  ( msglvl .gt. 0  .and.  iPrint .gt. 0)  then
                     write (iPrint, 2100) errmax, jmax
                  endif
               end if
            end if

            if  ( mod( iter, idegen ) .eq. 0 )  then

c              Every  idegen  iterations, reset  featol  and
c              move  x  on to the working set if it is close.

c                << cmdgen >>
              
               call CLKBEG (24)
               call CMDGEN ( 'End of cycle', msglvl,
     1                       n, nclin, nmoved, iter, numinf,
     2                       istate, bl, bu, featol, featlu, x,
     3                       iPrint, iSumm,
     4                       tolinc, idegen, kdegen, ndegen,
     5                       itnfix, nfix  )
               call CLKSUM (24)
               
               nviol = nviol + nmoved
               
            end if

            if  ( nviol .gt. 0 )  then 
               msg    = 'resetx'
               call CLKSUM (20)
               go to 100
            end if

c           ------------------------------------------------------------
c           Compute the QP objective and transformed gradient.
c           ------------------------------------------------------------

            if  ( cset )  then
               objqp = dot_product(cvec(1:n),x(1:n))
            else
               objqp = zero
            end if

            jthcol = 0

            call CLKBEG (25)
            call qpHess ( n, ldH, nHess, jthcol, H, x, delta, Hx )
            
            objqp  = objqp + half*dot_product(Hx(1:n),x(1:n))
            
            gq(1:n) = Hx(1:n)
            
c             << cmqmul >>
            call CMQMUL ( 6, n, nZ, nfree, ldQ, unitQ,
     1                   kx, gq, Q, wrk )

            if  ( cset)  then
              do j=1,n
                gq(j) = gq(j) + cq(j)
              enddo
            endif
            call CLKSUM (25)
            call CLKSUM (20)
               
         end if
         

         go to 100
         
c+    end while
      end if
      
c     ======================end of main loop============================

      return

 2100 format(  ' XXX  Iterative refinement.  The maximum violation is ',
     $           1p, e14.2, ' in constraint', i5 )

c     end of QPCORE ( qpcore )
      
      end
