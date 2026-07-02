      subroutine   LPCORE   ( prbtyp, msg,
     1                        cset,Rset, unitQ, 
     3                        iter, itmax, jinf, nviol,
     4                        n, nclin, ldA,nactiv, nfree, nZr, nZ,
     6                        istate, kactiv, kx,itprnt,
     8                        obj, numinf, xnorm,A, Ax, bl, bu,
     A                        cvec, featol, featlu, x, 
     C                        msglvl, iPrint, iSumm , lines1, lines2,
     D                        header, prnt, ldT   , ldQ, ldR,
     F                        tolOpt, epspt9, Asize , dTmax , dTmin,
     H                        tolinc, idegen, kdegen, ndegen, itnfix, 
     J                        nfix  ,alfa  , trulam, isdel , jdel  , 
     K                        jadd, kcheck, minsum, bigbnd, bigdx ,
     M                        Anorm, Ad, d, gq, cq, rlam, R, T, Q,
     N                        wtinf, wrk, p1wkmn, p1wkmx, p1wkav ) 
c     ==================================================================
c     ==================================================================
c     ====  LPCORE /                                                ====
c     ====  lpcore -- solve linear programming problem              ====
c     ==================================================================
c     ==================================================================

      integer            iter, itmax, jinf, nviol,
     1                   n, nclin, ldA,
     2                   nactiv, nfree, nZr, nZ,
     3                   numinf, iPrint, iSumm , lines1, lines2,
     4                   ldT   , ldQ  , ldR , idegen, kdegen, ndegen,
     5                   itnfix, msglvl, kcheck, isdel , jdel  ,
     6                   jadd,
     7                   p1wkmn, p1wkmx

      logical            header, prnt, minsum

      double precision   obj   , xnorm , tolOpt, epspt9,  Asize ,
     1                   dTmax , dTmin , tolinc, bigbnd, bigdx , alfa  ,
     2                   trulam, p1wkav
      
      character(len=2)   prbtyp
      
      character(len=6)   msg
      
      logical            cset, unitQ, Rset

      integer            istate (n+nclin), kactiv (n), kx (n), nfix (2)

      double precision   A (ldA,*), Ax (*), bl (n+nclin), bu (n+nclin),
     1                   cvec (*), featol (n+nclin), featlu (n+nclin),
     2                   x (n)

      double precision   Anorm (*), Ad (*), d (*), gq (*), cq (*),
     1                   rlam (*), R (ldR, *), T (ldT, *), Q(ldQ, *),
     2                   wtinf (*), wrk (*)

      external           itprnt

c     ==================================================================
c
c     derived from qpopt version 1.0
c     last modification -- 26-July-1996
c
c          This version of  lpcore  dated  10-Apr-94.
c
c     Copyright  1988/1994  Optimates.
c     LPCORE / lpcore  is a subroutine for linear programming.
c     On entry, it is assumed that an initial working set of
c     linear constraints and bounds is available.  The arrays  istate,
c     kactiv  and  kx  will have been set accordingly
c     and the arrays  T  and  Q  will contain the TQ factorization of
c     the matrix whose rows are the gradients of the active linear
c     constraints with the columns corresponding to the active bounds
c     removed.  The TQ factorization of the resulting (nactiv by nfree)
c     matrix is  A(free)*Q = (0 T),  where Q is (nfree by nfree) and T
c     is upper-triangular.
c
c     kactiv holds the general constraint indices in the order in which 
c     they were added.  The reverse ordering is used for T since new
c     rows are added at the front of T.  
c
c     Over a cycle of iterations, the feasibility tolerance featol
c     increases slightly (from tolx0 to tolx1 in steps of tolinc).
c     this ensures that all steps taken will be positive.
c
c     After idegen consecutive iterations, variables within featol of
c     their bounds are set exactly on their bounds and iterative
c     refinement is used to satisfy the constraints in the working set.
c     Featol is then reduced to tolx0 for the next cycle of iterations.
c
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
c     This version of  lpcore  dated  10-Apr-94.
c
c     Copyright  1988/1994  Optimates.
c     ==================================================================

      integer            iadd  , ifix  , inform, is    , iT    , j     ,
     1                   jbigst, jmax  , jsmlst, jtiny , kbigst, kdel  ,
     2                   ksmlst, nfixed, ngq   , nmoved, nctotl, notOpt, 
     3                   ntfixd
      
      logical            mnsmvl
      save               mnsmvl
      
      logical            hitlow, lp    , fp    , move  , onbnd,
     1                   overfl, unbndd

      double precision   alfap , alfhit, bigalf, biggst, condmx, condRz, 
     1                   condT , dinky , dnorm , Dzz   , errmax, flmax ,
     2                   gfnorm, gznorm, gZrnrm, objsiz, smllst, suminf,
     3                   tinyst, trubig, trusml, wssize, zerolm
      
      character(len=2)   plabel
      
      character(len=6)   empty
      parameter        ( empty  = '      ')

      double precision   zero, one
      parameter        ( zero = 0.0d0, one = 1.0d0 )

      double precision   DDIVFN, hdmcon

      external           DDIVFN, hdmcon

c     ==================================================================

c     Specify the machine-dependent parameters.

      flmax  = hdmcon (2)

      if  ( cset )  then
         ngq = 2
      else 
         ngq = 1
      end if

      if  ( prbtyp .eq. 'lp'  .or.  prbtyp .eq. 'LP' )  then
         lp     = .true.
         fp     = .false.
         plabel = 'LP'
      else
         lp     = .false.
         fp     = .true.
         plabel = 'FP'
      endif

      iT     = 1

c     -------------------------
c     First entry.  Initialize.
c     -------------------------
      
      if  ( iter .eq. 0 )  then
         jadd    =  0
         jdel    =  0
         isdel   =  0
         mnsmvl  = .false.
         alfa    =  zero
         Dzz     =  one
      end if

      nctotl = n  + nclin
      nviol  = 0

      condmx  =  flmax

      call CLKBEG (30)
c       << lpismi >>  << was cmsinf >>
      call LPSINF ( n, nclin, ldA,
     1              istate, bigbnd, numinf, suminf,
     2              bl, bu, A, featol,
     3              gq, x, wtinf, Ax, wrk )
      call CLKSUM (30)

      if  ( numinf .gt. 0 )  then
c          << cmqmul >>
         call CMQMUL ( 6, n, nZ, nfree, ldQ, unitQ,
     1                 kx, gq, Q, wrk )
      else
     1if  ( lp )  then
         gq(1:n) = cq(1:n)
      end if

      if  ( numinf .eq. 0  .and. lp )  then
         obj = zero
         do i=1,n
           obj = obj + cvec(i)*x(i)
         enddo
      else
         obj = suminf
      end if

      msg    = empty

c*    ======================Start of main loop==========================
c     +    do while (msg .eq. empty)
c*    ==================================================================      
  100 continue
      if  ( msg .eq. empty )  then

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
         if  ( nfree .gt. 0  .and.  nactiv .gt. 0)  then
            gfnorm = zero
            do i=1,nfree
              gfnorm = gfnorm + gq(i)**2
            enddo
            gfnorm = sqrt(gfnorm)
         endif

c        ---------------------------------------------------------------
c        Print the details of this iteration.
c        ---------------------------------------------------------------
c        Define small quantities that reflect the size of x, R and
c        the constraints in the working set.

         if  ( prnt )  then
            condT  = one
            if  ( nactiv .gt. 0)  then
c                     << ddiv >>
               condT  = DDIVFN  ( dTmax, dTmin, overfl )
            endif

c             << itprnt >>
            call itprnt ( plabel, header, Rset, 
     1                    msglvl, iter,
     2                    isdel, jdel, jadd,
     3                    n, nclin, nactiv,
     4                    nfree, nZ, nZr,
     5                    ldR, ldT, istate,
     6                    iPrint, iSumm , lines1, lines2,
     7                    alfa, condRz, condT,
     8                    Dzz, gznorm,
     9                    numinf, suminf, notOpt, obj, trulam,
     A                    Ax, R, T, x,
     B                    wrk )
            jdel  = 0
            jadd  = 0
            alfa  = zero
         end if

         if  ( numinf .gt. 0 )  then
            dinky  = tolOpt*abs( suminf )
         else
            objsiz = one  + abs( obj    )
            wssize = zero
            if  ( nactiv .gt. 0)  then
               wssize = dTmax
            endif
            dinky  = tolOpt*max( wssize, objsiz, gfnorm )
         end if

c        If the reduced gradient Z'g is small enough,
c        Lagrange multipliers will be computed.

         if  ( numinf .eq. 0  .and.  fp )  then 
            msg    = 'feasbl'
            nfixed = n - nfree
            rlam(1:nactiv+nfixed) = zero
            go to 100
         end if

         if  ( gZrnrm .le. dinky )  then
            
c           ============================================================
c           The point  x  is a constrained stationary point.
c           Compute Lagrange multipliers.
c           ============================================================
c           Define what we mean by 'tiny' and non-optimal multipliers.

            call CLKBEG (31)
            notOpt =   0
            jdel   =   0
            zerolm = - dinky
            smllst = - dinky
            biggst =   dinky + one
            tinyst =   dinky

c             << cmmul1 >>
            call CMMUL1 ( plabel, msglvl,
     1                    n, ldA, ldT,
     2                    nactiv, nfree, nZ,
     3                    istate, kactiv, kx,
     4                    zerolm, notOpt, numinf, 
     5                    trusml, smllst, jsmlst, ksmlst,
     6                    tinyst, jtiny, jinf,
     7                    trubig, biggst, jbigst, kbigst,
     8                    A, Anorm, gq,
     9                    rlam, T, wtinf,
     A                    iPrint )

            if  ( nZr .lt. nZ )  then
c                << cmmul2 >>
               call CMMUL2 ( msglvl, n, nZr, nZ,
     1                       zerolm, notOpt, numinf, 
     2                       trusml, smllst, jsmlst,
     3                       tinyst, jtiny, gq,
     4                       iPrint )
            end if

            if  ( abs(jsmlst) .gt. 0 )  then
               
c              ---------------------------------------------------------
c              Delete a constraint.
c              ---------------------------------------------------------

c              cmmul1 (CMMUL1) or  cmmul2 (CMMUL2)  found a non-optimal
c              multiplier.
              

               trulam = trusml
               jdel   = jsmlst

               if  ( jsmlst .gt. 0 )  then

c                 Regular constraint.

                  kdel   = ksmlst
                  isdel  = istate(jdel)
                  istate(jdel) = 0
                  
               end if
               
            else
     1      if  ( minsum )  then
               
               if  ( numinf .gt. 0  .and.  jbigst .gt. 0 )  then

c                 No feasible point exists for the constraints but the
c                 sum of the constraint violations can be reduced by
c                 moving off constraints with multipliers greater than 1.

                  jdel  = jbigst
                  kdel  = kbigst
                  isdel = istate(jdel)
                  
                  if  ( trubig .le. zero)  then
                     is = - 1
                  else
                     is = - 2
                  endif
                  
                  istate (jdel) = is
                  trulam        = trubig
                  mnsmvl        = .true.
                  numinf        = numinf + 1
                  
               end if
               
            end if

            if  ( jdel .eq. 0 )  then
               
               if  ( numinf .gt. 0 )  then 
                  msg = 'infeas'
               else
                  msg = 'optiml'
               end if
               call CLKSUM (31)
               go to 100
               
            end if

c           Constraint  jdel  has been deleted.
c           Update the  TQ  factorization.

c             << rzdel >>

            call RZDEL  ( unitQ, iT,
     1                   n, nactiv, nfree, ngq, nZ, nZr,
     2                   ldA, ldQ, ldT,
     3                   jdel, kdel, kactiv, kx,
     4                   A, T, gq, Q, d, rlam,
     5                   dTmax, dTmin )

            if  ( Rset)  then
c                << lpcolr >>
               call LPCOLR ( nZr, ldR, R, one )
            endif

            prnt    = .false.
            call CLKSUM (31)
            
        else
            
c           ============================================================
c           Compute a search direction.
c           ============================================================

            if  ( iter .ge. itmax )  then 
               msg    = 'itnlim'
               go to 100
            end if

            prnt  = .true.
            iter  = iter  + 1

            p1wkmn = min ( nZr, p1wkmn )
            p1wkmx = max ( nZr, p1wkmx )
            p1wkav = p1wkav + nZr

            dnorm = zero
            do i=1,nZr
              d(i) = -gq(i)
              dnorm = dnorm + d(i)**2
            enddo
            dnorm  = sqrt(dnorm)

c             << cmqmul >>

            call CMQMUL ( 1, n, nZr, nfree, ldQ, unitQ,
     1                    kx, d, Q, wrk )
            call CLKBEG(12)
            call dgemv  ( 'No transpose', nclin, n, one, A, ldA,
     1                                    d, 1, zero, Ad, 1 )
            call CLKSUM(12)

c           ------------------------------------------------------------
c           Find the constraint we bump into along d.
c           Update  x  and  Ax  if the step alfa is nonzero.
c           ------------------------------------------------------------
c           alfhit is initialized to bigalf. If it remains that value
c           after the call to  cmchzr, it is regarded as infinite.

c                  << ddiv >>
            bigalf = DDIVFN  ( bigdx, dnorm, overfl )

c             << lpchzr >>
            call CLKBEG (27)
            call LPCHZR ( mnsmvl, n, nclin,
     1                    istate, bigalf, bigbnd, dnorm,
     2                    hitlow, move, onbnd, unbndd,
     3                    alfhit, alfap, jadd, 
     4                    Anorm, Ad, Ax,
     5                    bl, bu, featol, featlu, d, x,
     6                    epspt9, tolinc, ndegen )
            call CLKSUM (27)

            if  ( unbndd )  then 
               msg    = 'unbndd'
               go to 100
            end if

            alfa   = alfhit
            xnorm  = zero
            do j=1,n
              x(j) = x(j) + alfa*d(j)
              xnorm = xnorm + x(j)**2
            enddo
            xnorm = sqrt(xnorm)

c           ... update  A times x
            
            do j=1,nclin
              Ax(j) = Ax(j) + alfa*Ad(j)
            enddo
               
c           ------------------------------------------------------------
c           Add a constraint to the working set.
c           Update the  TQ  factors of the working set.
c           Use  d  as temporary work space.
c           ------------------------------------------------------------

            if  ( bl(jadd) .eq. bu(jadd) )  then
               istate(jadd) = 3
            else
     1      if  ( hitlow )  then
               istate(jadd) = 1
            else
               istate(jadd) = 2
            end if

            if  ( jadd .gt. n )  then
               
               iadd = jadd - n
               
            else
               
               if  ( alfa .ge. zero )  then
                  if  ( hitlow )  then
                     x(jadd) = bl(jadd)
                  else
                     x(jadd) = bu(jadd)
                  end if
               end if
               
               ifixloop: do ifix = 1, nfree
                  if  ( kx(ifix) .eq. jadd) exit ifixloop
               enddo ifixloop
               
            end if

c             << rzadd >>
            call CLKBEG (28)
            call RZADD  ( unitQ, Rset,
     1                   inform, ifix, iadd, jadd, iT,
     2                   nactiv, nZ, nfree, nZr, ngq,
     3                   n, ldA, ldQ, ldR, ldT,
     4                   kx, condmx, Dzz,
     5                   A, R, T, gq, Q,
     6                   wrk, rlam, d,
     7                   iPrint, epspt9, Asize, dTmax, dTmin )
            call CLKSUM (28)

c     ******************************************************************
c     ******************************************************************
c     *** note -- inform return appears not to be used!!!!          ****
c     ***         do we ever reuse the gradient gq before           ****
c     ***         recomputing?                                      ****
c     ******************************************************************
c     ******************************************************************

            nZ    = nZ  - 1
            nZr   = nZr - 1

            if  ( jadd .le. n )  then

c              A simple bound has been added.

               nfree  = nfree  - 1
               
            else

c              A general constraint has been added.

               nactiv = nactiv + 1
               kactiv(nactiv) = iadd
               
            end if

c           Increment featol.

            do j=1,nctotl
              featol(j) = featol(j) + tolinc*featlu(j)
            enddo

            if  ( mod ( iter, kcheck ) .eq. 0 )  then
               
c              ---------------------------------------------------------
c              Check the feasibility of constraints with non-negative 
c              istate values.  If some violations have occurred, force 
c              iterative refinement and switch to phase 1.
c              ---------------------------------------------------------
               
c                << cmfeas >>
               call CLKBEG (29)
               call CMFEAS ( n, nclin, istate,
     1                       bigbnd, nviol, jmax, errmax,
     2                       Ax, bl, bu, featol, x )
               call CLKBEG (29)

               if  ( nviol .gt. 0 )  then
                  if  ( msglvl .gt. 0  .and.  iPrint .gt. 0 )  then
                     write (iPrint, 2100) errmax, jmax
                  end if
               end if
            end if

            if  ( mod( iter, idegen ) .eq. 0 )  then

c              Every  idegen  iterations, reset  featol  and
c              move  x  on to the working set if it is close.

c                << cmdgen >>
               call CLKBEG (29)
               call CMDGEN ( 'End of cycle', msglvl,
     1                       n, nclin, nmoved, iter, numinf,
     2                       istate, bl, bu, featol, featlu, x,
     3                       iPrint, iSumm,
     4                       tolinc, idegen, kdegen, ndegen,
     5                       itnfix, nfix  )
               call CLKSUM (29)

               nviol = nviol + nmoved
               
            end if

            if  ( nviol .gt. 0 )  then 
               msg    = 'resetx'
               go to 100
            end if

            if  ( numinf .ne. 0 )  then
               
c              << lpismi >>  << was cmsinf >>
               call CLKBEG (30)
               call LPSINF ( n, nclin, ldA,
     1                       istate, bigbnd, numinf, suminf,
     2                       bl, bu, A, featol,
     3                       gq, x, wtinf, Ax, wrk )

               if  ( numinf .gt. 0 )  then
c                   << cmqmul >>
                  call CMQMUL ( 6, n, nZ, nfree, ldQ, unitQ,
     1                          kx, gq, Q, wrk )
               else
     1         if  ( lp )  then
                  gq(1:n) = cq(1:n)
               end if
               call CLKSUM (30)
               
            end if

            if  ( numinf .eq. 0  .and.  lp )  then
               obj = zero
               do i=1,n
                 obj = obj + cvec(i)*x(i)
               enddo
            else
               obj = suminf
            end if
            
         end if
         go to 100
         
      end if
      
c     ======================end of main loop============================
c+    end while
c     ======================end of main loop============================
c
      if  ( msg .eq. 'optiml' )  then
         
         if  ( lp )  then
            
            if  ( nZr .lt. nZ )  then
               msg = 'weak  '
            else 
               ntfixd = 0
               do j = 1, n
                  if  ( istate(j) .eq. 4 )  ntfixd = ntfixd + 1
               enddo
               if  ( ntfixd .gt. 0 )  msg = 'weak  '
            end if
            if  ( abs(jtiny) .gt. 0 )  msg = 'weak  '
            
         end if
         
      else
     1if  ( msg .eq. 'unbndd'  .and. numinf .gt. 0 )  then
         msg = 'infeas'
      end if

      return

 2100 format(  ' XXX  Iterative refinement.  The maximum violation is ',
     $           1p, e14.2, ' in constraint', i5 )

c     end of LPCORE (lpcore)
      
      end
