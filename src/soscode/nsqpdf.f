      subroutine   NSQPDF   ( ndim  , ncon  , prbtyp,
     1                        start , lniwrk, lnrwrk, bigbnd, bigstp,
     2                        crstol, featol, opttol, rnktol, levtol,
     3                        nhess , fealim, optlim, expand,
     5                        lbigbn, lbigst, lcrstl, lfeatl, lopttl,
     6                        lrnktl, llevtl, lnhess, lfealm, loptlm,
     7                        lexpnd, lkchck, cmpcod )
c     ==================================================================
c     ==================================================================
c     ====  NSQPDF /                                                ====
c     ====  nsqpdf -- assign default values to alg. parameters      ====
c     ==================================================================
c     ==================================================================

      integer            ndim  , ncon  ,
     1                   prbtyp, start , lniwrk, lnrwrk,
     2                   nhess , fealim, optlim, expand, 
     3                   lnhess, lfealm, loptlm, lexpnd, lkchck,
     4                   cmpcod

      double precision   bigbnd, bigstp, crstol, featol, opttol, rnktol,
     1                   levtol, lbigbn, lbigst, lcrstl, lfeatl, lopttl,
     2                   lrnktl, llevtl 

c     ... local variables and constants

      double precision   epslon, gigant, hundrd, one,point8, 
     1                   rteps , wrktol, zero

      parameter        ( zero   =   0.0d0,
     1                   one    =   1.0d0,
     2                   hundrd = 100.0d0,
     4                   point8 =   8.0d-1,
     3                   wrktol = one / hundrd )

      double precision   hdmcon

      external           hdmcon

c     ==================================================================

c     ... last modified 24-July-1996

c     ==================================================================
      
c     -------------------------------------------
c     ... assign values to defaultable parameters
c     -------------------------------------------

      lbigbn = cmpcod
      lbigbn = lniwrk
      lbigbn = lnrwrk
      lbigbn = start

      epslon = hdmcon (5)
      
      rteps  = sqrt ( epslon )
      gigant = 1.0 / (hundrd*epslon)
      
c     ... setting for "infinity" in context of bounds
      
      if  ( bigbnd .le. zero )  then
         lbigbn = gigant
      else
         lbigbn = bigbnd
      endif

c     ... setting for "infinity" in context of step length
      
      if  ( bigstp .le. zero )  then
         lbigst = lbigbn
      else
         lbigst = bigstp
      endif

c     ... tolerance for feasibility acceptance in "crash starts"
      
      if  ( crstol .le. zero  .or.  crstol .gt. one )  then
c       "loose"         lcrstl = wrktol
         lcrstl = epslon
      else
         lcrstl = crstol
      endif

c     ... tolerance for absolute feasibility tolerance
      
      if  ( featol .lt. epslon )  then
         lfeatl = rteps
      else
         lfeatl = featol
      endif

c     ... tolerance for determining whether LaGrange multipliers
c         have the correct sign in optimality tests
      
      if ( opttol .lt. epslon )  then
         lopttl = epslon ** point8
      else
         lopttl = opttol
      endif

c     ... tolerance for "positive definite" reduced Hessian
      
      if ( rnktol .le. zero )  then
         lrnktl = epslon * hundrd
      else
         lrnktl = rnktol
      endif

c     ... tolerance for "well-conditioned" reduced Hessian
      
      if ( levtol .lt. one )  then
         llevtl = one / lrnktl ** (0.66)
      else
         llevtl = levtol
      endif

c     ... dimension of non-linear part of problem (number of rows
c         in leading principal minor of  H  that really represent
c         the non-linear problem).  set to zero locally for
c         linear problems

      if  ( prbtyp .lt. 2 ) then
         lnhess = 0
      else
         if  ( nhess .lt. zero )  then
            lnhess = ndim
         else
            lnhess = nhess
         endif
      endif

c     ... feasibility phase iteration limit

      if  ( fealim .lt. 0 )  then
         lfealm = 100 * (ndim+ncon)
      else
         lfealm = fealim
      endif

c     ... optimalility phase iteration limit

      if  ( optlim .lt. 0 )  then
         loptlm = 100 * (ndim+ncon)
      else
         loptlm = optlim
      endif

c     ... setting for anti-cycling frequency expansion parameter

      if  ( expand .le. 0 )  then
         lexpnd = 5
      else
     1     if  ( expand .gt. 9 999 999 )  then
         lexpnd = 9 999 999
      else
         lexpnd = expand
      endif

c     ... setting for checking feasibility of constraints
c         (every  kcheck  iterations).  not passed out to user
c         control currently

      lkchck    = 50

      return

c     end of NSQPDF / nsqpdf

      end
