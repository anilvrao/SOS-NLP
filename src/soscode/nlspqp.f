      subroutine   NLSPQP   ( ndim  , ncon  , ldH   , hmat  , hesstv, 
     1                        cvct  , ldA   , amat  , bupr  , blwr  ,
     2                        xupr  , xlwr  , bigbnd, output, msglvl,
     3                        lintrm, start , lniwrk, lnrwrk, xvec  ,
     4                        iwork , rwork , conmlt, varmlt, istatc,
     5                        istatv, quad  , fpiter, qpiter, nouter,
     6                        nlvmod, delta , p1wkmn, p1wkmx, p1wkav,
     7                        p2wkmn, p2wkmx, p2wkav, lminzh, lmaxzh,
     8                        sminaa, smaxaa, cmpcod  )
c     ==================================================================
c     ==================================================================
c     ====  NLSPQP /                                                ====
c     ====  nlspqp -- null space quadratic programming module       ====
c     ====            with Levenberg-Marquardt modifications        ====
c     ====            for dense Hessians and dense constraints      ====
c     ==================================================================
c     ==================================================================
c     ===  last modified 03-Sept-1996                               ====
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
c             BIGBND  big bound value; any bounds as large as this
c                     value will be treated as infinite.
c                     (inputting a non-positive value for bigbnd
c                      results in a default value of 1 / 100*espilon)
c
c             ----    Algorithm Related Data
c
c             OUTPUT  output unit no.
c             MSGLVL  output control flag
c                     = 0 --- no print
c                     msglvl .gt. 100 --- debug print
c             LINTRM  is linear term present? (logical)
c             START   integer start option flag
c                     = -1 --- hot start
c                     =  0 --- warm start
c                     =  1 --- cold start
c                     =  2 --- cold start
c             LNIWRK  length of floating point work array  iwork
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
c             P1WKAV  average size of working set during Phase I
c             P2WKMN  smallest working set during Phase II
c             P2WKMX  largest working set during Phase II
c             P2WKAV  average size of working set during Phase II
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
c                     0    --- no constraints active at solution
c                              (active constraint matrix is empty matrix)
c                     -1   --- no singular value estimate available
c             SMAXAA  largest singular value of final active
c                     constraint matrix (Jacobian of final subproblem)
c                     > 0  --- positive definite
c                              final active constraint matrix
c                     0    --- no constraints active at solution
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
c                     = -1006 --- internal program check
c                     = -1007 --- internal program check
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
c                           to NLSPQP.  In particular, assume that
c                           the working vectors include the TQ factors
c                           of the active set given by the ISTAT*
c                           vectors.  No error checking is done in a
c                           hot start.  The TQ factorization
c                           stored in the working vector is used again.
c                            
c     ===================================================================

c     --------------
c     ... parameters
c     --------------
      
      integer            ndim  , ncon  , ldH   , ldA   , lniwrk, lnrwrk,
     1                   output, msglvl, start , fpiter, qpiter, nouter,
     2                   nlvmod, p1wkmn, p1wkmx, p2wkmn, p2wkmx, cmpcod

      logical            lintrm
      
      double precision   bigbnd, quad  , delta , p1wkav, p2wkav, lminzh, 
     1                   lmaxzh, sminaa, smaxaa

      integer            istatc (*), istatv (ndim), iwork (lniwrk)

      external           hesstv

      double precision   hmat (ldH, ndim), amat (ldA, ndim)

      double precision   cvct  (ndim)  , bupr   (*), blwr   (*),
     1                   xvec  (ndim)  , xupr   (ndim), xlwr   (ndim), 
     2                   rwork (lnrwrk), conmlt (*), varmlt (ndim)

c     ===================================================================

c     ... last modified 26-March-1996

c     ... this routine serves as a shell for the expert version 'nlspqx'.
c         additional parameters are all given values to cause their
c         defaults to be used.
c     
c     ===================================================================

c     -------------------
c     ... local variables
c     -------------------

      integer           prbtyp, nhess , fealim, optlim, expand

      logical           minsum

      double precision  bigstp, crstol, featol, opttol, rnktol,
     1                  cndtol

      parameter       ( prbtyp = 2,
     1                  nhess  = -1,
     2                  fealim = -1,
     3                  optlim = -1,
     4                  expand = -1,
     5                  minsum = .true.,
     6                  bigstp = -1.0d0,
     7                  crstol = -1.0d0,
     8                  featol = -1.0d0,
     9                  opttol = -1.0d0,
     A                  rnktol = -1.0d0,
     B                  cndtol = -1.0d0 )
      
c     ===================================================================

c       << nlspqx >>
      call NLSPQX ( ndim  , ncon  , ldH   , hmat  , hesstv, 
     1              cvct  , ldA   , amat  , bupr  , blwr  ,
     2              xupr  , xlwr  , output, msglvl, prbtyp,
     3              lintrm, start , minsum, lniwrk, lnrwrk,
     4              xvec  , iwork , rwork , conmlt, varmlt,
     5              istatc, istatv, quad  , bigbnd, bigstp,
     6              crstol, featol, opttol, rnktol, cndtol,
     7              nhess , fealim, optlim, expand, fpiter,
     8              qpiter, nouter, nlvmod, delta , p1wkmn,
     9              p1wkmx, p1wkav, p2wkmn, p2wkmx, p2wkav,
     A              lminzh, lmaxzh, sminaa, smaxaa, cmpcod  )
      
c     ----------
c     ... return
c     ----------

      return

c     end of NLSPQP / nlspqp
      
      end
