      subroutine   NLSPQX   ( ndim  , ncon  , ldH   , hmat  , hesstv, 
     1                        cvct  , ldA   , amat  , bupr  , blwr  ,
     2                        xupr  , xlwr  , output, msglvl, prbtyp,
     3                        lintrm, start , minsum, lniwrk, lnrwrk,
     4                        xvec  , iwork , rwork , conmlt, varmlt,
     5                        istatc, istatv, quad  , bigbnd, bigstp,
     6                        crstol, featol, opttol, rnktol, cndtol,
     7                        nhess , fealim, optlim, expand, fpiter,
     8                        qpiter, nouter, nlvmod, delta , p1wkmn,
     9                        p1wkmx, p1wkav, p2wkmn, p2wkmx, p2wkav,
     A                        lminzh, lmaxzh, sminaa, smaxaa, cmpcod  )
c     ==================================================================
c     ==================================================================
c     ====  NLSPQX /                                                ====
c     ====  nlspqx -- null space quadratic programming module       ====
c     ====            with Levenberg-Marquardt modifications        ====
c     ====            for dense Hessians and dense constraints      ====
c     ====            ("expert version")                            ====
c     ==================================================================
c     ==================================================================
c     ===  last modified  03-Sept-1996                              ====
c     ==================================================================

c
c         PURPOSE:
c
c             minimize the quadratic objective 
c
c               quad = .5*(xvec**t)*hmat*xvec + (cvct**t)*xvec
c
c             subject to the ncon general linear constraints
c
c               blwr .le. amat*xvec .le. bupr
c
c             and the ndim bounds
c
c               xlwr .le. xvec .le. xupr
c
c             equality constraints are imposed by setting blwr = bupr.
c             both the hessian HMAT and the jacobian amat are treated
c             as dense matrices.
c
c             the algorithm implements the null space quadratic
c             programming method, proposed by gill, murray and saunders.
c
c     REF.  "User's Guide for QPOPT 1.0:  A Fortran Package
c            for Quadratic Programming"  Philip E. Gill, Walter Murray
c                   and Michael A. Saunders.
c
c            In addition, the implementation includes automatic
c            Levenberg-Marquardt modification to solve the related
c            problem:
c               minimize the quadratic objective 
c
c                  quad = .5*(xvec**t)*(hmat + delta*I)*xvec +
c                         (cvct**t)*xvec
c
c               subject to the same constraints
c            in cases in which  HMAT  is indefinite.
c
c     
c        +--------------------------------+
c        | PARAMETERS:                    |
c        |    are grouped by usage:       |
c        |                                |
c        |    required or standard:       |
c        |        input parameters        |
c        |        input/output parameters |
c        |        output parameters       |
c        |                                |
c        |    optional or diagnostic:     |
c        |        input parameters        |
c        |        output parameters       |
c        +--------------------------------+
c
c     
c        ------------------------------
c        | REQUIRED INPUT PARAMETERS: |
c        ------------------------------
c
c             ----    Problem Data
c
c             NDIM    number of variables
c             NCON    number of general linear constraints
c
c             ----    Objective Function Data
c
c             LDH     leading declared dimension for hessian
c             HMAT    hessian matrix (ldh*ndim)
c             HESSTV  subroutine to perform hessian x vector mult.
c             CVCT    objective function linear term (ndim)
c
c             ----    Constraint Data (Dummy arguments when NCON = 0)
c
c             LDA     leading declared dimension for jacobian
c             AMAT    jacobian matrix (lda*ndim)
c             BUPR    constraint upper bound vector (ncon)
c             BLWR    constraint lower bound vector (ncon)
c
c             ----    Simple Bounds
c
c             XUPR    upper bound for independent variables (ndim)
c             XLWR    lower bound for independent variables (ndim)
c
c             ----    Algorithm Related Data
c
c             OUTPUT  output unit no.
c             MSGLVL  output control flag
c                     = 0 --- no print
c                     msglvl .gt. 100 --- debug print
c             PRBTYP  problem type
c                     = 0 -- find a feasible point
c                     = 1 -- solve linear program
c                     = 2 -- solve quadratic program
c             LINTRM  is linear term present? (logical)
c                     (must be false for  prbtyp .eq. 0  and
c                      must be true  for  prbtyp .eq. 1 )
c             START   integer start option flag
c                     = -1 --- hot start
c                     =  0 --- warm start
c                     =  1 --- cold start
c                     =  2 --- cold start
c             MINSUM  obtain minimum of sums of infeasibilities
c                     in case no feasible point exists (logical)
c
c             LNIWRK  length of integer work array  iwork
c             LNRWRK  length of floating point work array  rwork
c
c     
c         +-----------------------------------+
c         | STANDARD INPUT/OUTPUT PARAMETERS: |
c         +-----------------------------------+
c
c             XVEC    on input, must hold an initial guess for
c                     independent variables; on output, holds the final
c                     values for independent variables (ndim)
c
c                     the following parameters are used as input
c                     parameters only for warm or hot starts.
c                     they are always output parameters.
c     
c             IWORK   work array (lniwrk)
c             RWORK   work array (lnrwrk)
c             ISTATC  integer constraint status (ncon)
c                     = -2 --- infeasible wrt upper bound
c                     = -1 --- infeasible wrt lower bound
c                     = 0  --- free (inactive) inequality
c                     = 1  --- fixed on lower bound
c                     = 2  --- fixed on upper bound
c                     = 3  --- equality constraint
c                     = 4  --- ignore this constraint
c             ISTATV  integer variable status (ndim)
c                   ( = -2 --- infeasible wrt upper bound )
c                   ( = -1 --- infeasible wrt lower bound )
c                     = 0  --- free variable 
c                     = 1  --- fixed on lower bound
c                     = 2  --- fixed on upper bound
c                     = 3  --- fixed permanently (equality)
c                   ( = 4  --- ignore this simple bound )
c                     the first two values are probably not possible
c                     as return values for ISTATV and should not 
c                     generally be used as input states for ISTATV.
c                     the capability to ignore simple bounds is
c                     not yet implemented, so a value of 4 for ISTATV
c                     will be rejected as an error.
c
c     
c         ------------------------------+
c         | STANDARD OUTPUT PARAMETERS: |
c         ------------------------------+
c     
c             CONMLT  lagrange multipliers for constraints (ncon)
c             VARMLT  lagrange multipliers for variables (ndim)
c             QUAD    optimal objective function value
c
c
c         +----------------------------+
c         | OPTIONAL INPUT PARAMETERS: |
c         +----------------------------+
c
c          for all of the optional input parameters, a non-positive
c          input value will cause the default value of the parameter
c          to be used.  note -- the optional parameters are interpretted
c          only on cold or warm starts; in a hot start situation,
c          the optional parameters used in the previous run are used.
c     
c             BIGBND  big bound value (def: 1 / 100*espilon)
c             BIGSTP  big step value (def: BIGBND)
c             CRSTOL  crash tolerance (def:  0.01 )
c             FEATOL  feasibility tolerance (def: square root of eps,
c                     mimimum allowed: eps)
c             OPTTOL  optimality tolerance (def: square root of eps)
c             RNKTOL  rank tolerance (def: 100 epsilon)
c             CNDTOL  well-conditionedness tolerance
c                     (def: 1.0 / ( rnktol ** 0.3) )
c
c             NHESS   number of nonzero rows in hessian (<= nvar)
c                     (def: nvars)
c             FEALIM  feasibility phase iteration limit (def: 50)
c             OPTLIM  optimality phase iteration limit (def: 50)
c             EXPAND  initial value for anti-cycling procedure
c                     (def: 5;  expand > 9 999 999 disables anti-
c                               cycling procedure)
c
c          note -- there are several storage allocation parameters
c                  ( mxfree -- number of non-fixed variables,
c                    mxactv -- maximum number of constraints that
c                              will ever be active,
c                    maxnZ  -- maximum dimension of any nullspace,
c                    minact -- minimum number of constraints ever
c                              active,
c                  that are interrelated and mildly affect the storage
c                  required.  Some reasonable set of these could be
c                  made available for user control.)
c
c         +-------------------------------+
c         | DIAGNOSTIC OUTPUT PARAMETERS: |
c         +-------------------------------+
c
c             FPITER  number of Phase 1 (feasible point) iterations
c             QPITER  number of Phase 2 iterations
c             NOUTER  number of outer iterations (each consisting
c                     of a Phase I and Phase II solution)
c             NLVMOD  number of times the Levenberg parameter
c                     was increased to make a reduced Hessian
c                     positive definite.
c                     (note:  the number of complete Cholesky
c                      factorizations required is the sum of the
c                      preceeding two variables.)
c             DELTA   accumulated modifications to levenberg parameter
c             P1WKMN  smallest working set during Phase I
c             P1WKMX  largest working set during Phase I
c             P1WKMN  average size of working set during Phase I
c             P2WKMN  smallest working set during Phase II
c             P2WKMX  largest working set during Phase II
c             P2WKMN  average size of working set during Phase II
c             LMINZH  smallest eigenvalue of final reduced
c                     hessian (final subproblem)
c                     > 0  --- positive definite
c                              final reduced Hessian
c                     0    --- no variables free at solution
c                              (reduced Hessian is empty matrix)
c                     -1   --- no eigenvalue estimate available
c                              because Phase 1 aborted
c             LMAXZH  largest eigenvalue of final reduced
c                     hessian (final subproblem)
c                     > 0  --- positive definite
c                              final reduced Hessian
c                     0    --- no variables free at solution
c                              (reduced Hessian is empty matrix)
c                     -1   --- no eigenvalue estimate available
c                              because Phase 1 aborted
c             SMINAA  smallest singular value of final active
c                     constraint matrix (Jacobian of final subproblem)
c                     > 0  --- positive definite
c                              final active constraint matrix
c                     0    --- no variables free at solution
c                              (active constraint matrix is empty matrix)
c                     -1   --- no singular value estimate available
c             SMAXAA  largest singular value of final active
c                     constraint matrix (Jacobian of final subproblem)
c                     > 0  --- positive definite
c                              final active constraint matrix
c                     0    --- no variables free at solution
c                              (active constraint matrix is empty matrix)
c                     -1   --- no singular value estimate available
c
c
c          +-------------------------+
c          | COMPLETION RETURN FLAG: |
c          +-------------------------+
c      
c             CMPCOD     integer error return flag
c                     = 0    --- success
c             ------------------------
c            .GT. 0 --- warning errors
c                     = 1     --- weak solution (singular hessian or
c                                 very small lagrange multipliers)
c                     = 2     --- unbounded solution 
c                     = 3     --- infeasible constraints
c                     = 4     --- one of iteration limits reached
c                     = 5     --- insufficient space for reduced hessian
c                                 (solution has fewer active constraints
c                                  than minimum allowed)
c            ---------------------------------------
c            .LE. 0 --- fatal (usually input) errors
c                     = -1001 --- ndim .le. 0
c                     = -1002 --- ncon .lt. 0
c                     = -1003 --- ldh .lt. ndim
c                     = -1004 --- lda .lt. ncon
c                     = -1005 --- invalid value for start
c                     = -1006 --- invalid value for prbtyp
c                     = -1007 --- linear term not present for prbtyp 1
c                                 or linear term present for prbtyp 0
c                     = -1008 --- msglvl  .lt.  0 and output .gt. 0
c                     = -1010 --- bupr(i)  .lt.  blwr(i) or
c                                     bupr(i) .ge. bigbnd and
c                                     blwr(i) .le. -bigbnd
c                     = -1011 --- xupr(i)  .lt.  xlwr(i) 
c                     = -1012 --- mequal + nfixvr .gt. ndim
c                     = -1014 --- invalid entry in input state variable
c     
c                     = -1020 --- qpopt rejected input parameter
c                     = -1021 --- qpopt did not recognize problem type
c     
c                     = -1100 --- insufficient real work storage
c                     = -1101 --- insufficient integer work storage
c
c                     = -1200 --- internal program check
c                     = -1201 --- internal program check
c
c         +---------------+
c         | START OPTIONS |
c         +---------------+
c
c             cold start -- determine everything from scratch; only
c                           an initial guess for the solution  x  is
c                           required
c             warm start -- assume that the ISTAT* vectors represent
c                           a reasonable guess for an initial working
c                           set
c             hot start  -- assume that the ISTAT*, RWORK and IWORK
c                           vectors are unchanged from a previous call
c                           to NLSPQX.  In particular, assume that
c                           the working vectors include the TQ factors
c                           of the active set given by the ISTAT*
c                           vectors.  No error checking is done in a
c                           hot start.  Optional input parameters are
c                           not referenced in a hot start; the values
c                           used in creating the TQ factorization
c                           stored in the working vector are used again.
c                            
c     ===================================================================

c     --------------
c     ... parameters
c     --------------
      
      integer            ndim  , ncon  , ldH   , ldA   , lniwrk, lnrwrk,
     1                   output, msglvl, prbtyp, start , nhess , fealim,
     2                   optlim, expand, fpiter, qpiter, nouter, nlvmod,
     3                   p1wkmn, p1wkmx, p2wkmn, p2wkmx, cmpcod

      logical            lintrm, minsum
      
      double precision   quad  , bigbnd, bigstp, crstol, featol, opttol,
     1                   rnktol, cndtol, delta , p1wkav, p2wkav, lminzh,
     2                   lmaxzh, sminaa, smaxaa

      integer            istatc (*), istatv (ndim), iwork (lniwrk)

      external           hesstv

      double precision   hmat (ldH, ndim), amat (ldA, ndim)

      double precision   cvct  (ndim)  , bupr   (*), blwr   (*),
     1                   xvec  (ndim)  , xupr   (ndim), xlwr   (ndim), 
     2                   rwork (lnrwrk), conmlt (*), varmlt (ndim)

c     -------------------
c     ... local variables
c     -------------------

      integer            inform, lexpnd, lfealm, lkchck, lnhess, lniwk2, 
     1                   lnrwk2, loptlm, maxact, maxnZ , mequal, mode  ,
     2                   mxfree, needed, nclin , nfixvr, sumfil

      integer            xax   , xcexch, xistat, xiwork, xlamda, 
     1                   xlower, xrwork, xupper

      logical            cset, error

      double precision   lbigbn, lbigst, lcrstl, lfeatl, lopttl, lrnktl,
     1                   lcndtl

      character(len=2)   cprbtp
      
      character(len=4)   qstart

      double precision   zero, one, two

      parameter        ( zero = 0.0d0, one = 1.0d0, two = 2.0d0 )

c     ==================================================================

      cmpcod = 0

c     ... start nullspace qp timing clocks

      call CLKBEG (4)
      call CLKBEG (5)

c     --------------------------
c     ... check input parameters
c     --------------------------

      if  ( start .lt. -1  .or.  start .gt. 2 )  then
         cmpcod = -1005
         go to 10000
      else
     1if  ( start .eq. -1 )  then      
         qstart = 'hot '
      else
     1if  ( start .eq. 0 )  then
         qstart = 'warm'
      else
         qstart = 'cold'
      endif

c     -------------------------------------------------------
c     ... check scalar parameters only on cold or warm starts
c     -------------------------------------------------------

      if  ( start .ge. 0 )  then

c          << nsqpcp >>
         call NSQPCP   ( output, msglvl, ndim  , ncon  , nHess, 
     1                   ldH   , ldA   , prbtyp, lintrm, cprbtp,
     2                   cset  , cmpcod )
         
         if  ( cmpcod .ne. 0 )  then
            go to 10000
         endif

      endif
      
c     -------------------------------------------------
c     ... assign values to user-defaulted parameters or
c         restore earlier settings
c     -------------------------------------------------

c       << nsqpdf >>
      call NSQPDF ( ndim  , ncon  , prbtyp,
     1              start , lniwrk, lnrwrk, bigbnd, bigstp,
     2              crstol, featol, opttol, rnktol, cndtol,
     3              nhess , fealim, optlim, expand, 
     4              lbigbn, lbigst, lcrstl, lfeatl, lopttl,
     5              lrnktl, lcndtl, lnhess, lfealm, loptlm,
     6              lexpnd, lkchck, cmpcod )

      if  ( cmpcod .ne. 0 )  then
         go to 10000
      endif

c     --------------------------------------------
c     ... check bounds only on cold or warm starts
c     --------------------------------------------

c       << nsqpcb >>
      call NSQPCB  ( ndim  , ncon  , lbigbn, bupr  , blwr  ,
     1               istatc, istatv, xupr  , xlwr  , nfixvr, 
     2               mequal, nclin , cmpcod )
         
      if  ( cmpcod .ne. 0 )  then
         go to 10000
      endif

c     -----------------------------------
c     ... storage allocation and checking
c     -----------------------------------

c       << nsqpst >>
      call NSQPST ( ndim  , ncon  , nclin , prbtyp, start , lniwrk,
     1              lnrwrk, nfixvr, mequal, xcexch, xistat, xiwork,
     2              xlower, xupper, xlamda, xax   , xrwork, maxact,
     3              mxfree, maxnZ , lniwk2, lnrwk2, iwork , cmpcod,
     4              needed )
      
      if  ( cmpcod .ne. 0 )  then
         go to 10000
      endif


c     ---------------------------------------------------------
c     ... convert problem to internal format by moving ignored
c         constraints to end of list (physically permute bounds
c         vectors and rows of a).  
c     ---------------------------------------------------------
         
      if  ( ncon .gt. nclin )  then
c          << cprmut >>
         call CPRMUT ( ndim, ncon, nclin, amat, lda, istatc, 
     1                 bupr, blwr, iwork(xcexch), error )
         if  ( error )  then
            cmpcod = -1201
            go to 10000
         endif
      endif

c     ----------------------------------------------------------
c     ... concatenate constraint and bound information for qpopt
c         and convert state information.
c     ----------------------------------------------------------

c       << ccncat >>
      call CCNCAT ( ndim, nclin, istatv, istatc, iwork (xistat),
     1              xlwr, blwr, rwork (xlower),
     2              xupr, bupr, rwork (xupper), error )

      if  ( error )  then
         cmpcod = -1014
         go to 10000
      endif
      
c     ----------------
c     ... invoke qpopt
c     ----------------

      sumfil = 0

c       << qpopt >>
      call QPOPT  ( ndim  , nclin , lda   , ldh   , lnHess,
     1              amat  , rwork (xlower), rwork (xupper), cvct, hmat,
     2              hesstv, iwork (xistat), xvec  ,
     3              cprbtp, cset  , minsum, qstart, 
     4              inform, quad  , 
     5              msglvl, output, sumfil, 
     6              rwork (xax   ), rwork (xlamda),
     7              iwork (xiwork), lniwk2, rwork (xrwork), lnrwk2,
     8              mxfree, maxact, maxnZ ,
     9              lbigbn, lbigst, lcrstl, lfeatl, lopttl, lrnktl,
     A              lcndtl,
     B              lfealm, loptlm, lexpnd, lkchck,
     B              fpiter, qpiter, nouter, nlvmod, delta , 
     C              p1wkmn, p1wkmx, p1wkav,
     D              p2wkmn, p2wkmx, p2wkav,
     C              lminzh, lmaxzh, sminaa, smaxaa )

c     ------------------------------
c     ... return solution error flag
c     ------------------------------

      if  ( inform .le. 5 )  then
         cmpcod = inform
      else
         cmpcod = -1200
      endif

c     --------------------------------------------------------
c     ... unpack output arrays and unpermute constraint matrix
c     --------------------------------------------------------
      
c       << cuncat >>
      call CUNCAT ( ndim, nclin, iwork (xistat), istatv, istatc, 
     1              rwork (xlamda), conmlt, varmlt, error )

      if  ( ncon .gt. nclin )  then
         
c          << cunprm >>
         call CUNPRM ( ndim, nclin, amat, lda, istatc, bupr, blwr,
     1                 iwork(xcexch) )

      endif

c     -----------------------------------
c     ... call error handler if necessary
c     -----------------------------------

10000 continue
      if  ( cmpcod.ne.0 )  then
         if  ( cmpcod .gt. 0 )  then
            mode = 0
         else
     1   if  ( cmpcod .gt. -1100 )  then
            mode = 1
         else
     1   if  ( cmpcod .gt. -1200 )  then
            mode = 2
         else
            mode = 3
         endif

         call hherr ( mode, 'NLSPQX  ', cmpcod, needed )

      endif

c     stop schur-qp timing clocks

      call CLKSUM(4)
      call CLKMAX(5)

c     ===================================================================

c     ----------
c     ... return
c     ----------

      return
      end
