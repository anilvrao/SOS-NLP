      subroutine   LPCHZR  ( mnsmvl, n, nclin,
     1                       istate, bigalf, bigbnd, pnorm,
     2                       hitlow, move, onbnd, unbndd,
     3                       alfa, alfap, jhit,
     4                       Anorm, Ap, Ax,
     5                       bl, bu, featol, featlu, p, x,
     6                       epspt9, tolinc, ndegen )
c     ==================================================================
c     ==================================================================
c     ====  lpchzr / LPCHZR -- choose a constraint to add ...       ====
c     ====                     general case (infeasible point)      ====
c     ==================================================================
c     ==================================================================

      integer            n, nclin, jhit, ndegen
      
      logical            mnsmvl, hitlow, move, onbnd, unbndd

      double precision   bigalf, bigbnd, pnorm, alfa, alfap, epspt9,
     1                   tolinc
      
      integer            istate (n+nclin)
      
      double precision   bl (n+nclin), bu (n+nclin)
      
      double precision   featol (n+nclin), featlu (n+nclin)
      
      double precision   Anorm (*), Ap (*), Ax (*)
      
      double precision   p (n), x (n)

c     ==================================================================
c
c     derived from qpopt version 1.0 cmchzr
c
c     special vectorized version of cmchzr
c     
c     last modification -- 26-July-1996
c
c         Original version written by PEG,  19-April 1988.
c         This version of  cmchzr  dated   6-Jul-1988.
c
c     LPCHZR / lpchzr  finds a step alfa such that the point
c         x + alfa*p reaches one of the linear constraints
c         (including bounds).
c
c     In this version of cmchzr, when x is infeasible, the number of
c     infeasibilities will never increase.  If the number stays the
c     same, the sum of infeasibilities will decrease.  If the number
c     decreases by one or more,  the sum of infeasibilities will usually
c     decrease also, but occasionally it will increase after the step
c     alfa  is taken.  (Convergence is still assured because the number
c     has decreased.)
c
c     Three possible steps are computed as follows:
c
c     alfaf = the maximum step that can be taken without violating
c              one of the constraints that are currently satisfied.
c
c     alfai = reaches a linear constraint that is currently violated.
c              Usually this will be the furthest such constraint along
c              p, subject to the angle between the constraint normal and
c              p being reasonably close to the maximum value among
c              infeasible constraints,  but if mnsmvl = .true. it will
c              be the first one along p.  The latter case applies only
c              when the problem has been determined to be infeasible,
c              and the sum of infeasibilities are being minimized.
c              (Alfai is not defined when x is feasible.)
c
c     Alfai is needed occasionally when infeasible, to prevent
c     going unnecessarily far when alfaf is quite large.  It will
c     always come into effect when x is about to become feasible.
c     (The sum of infeasibilities will decrease initially as alfa
c     increases from zero, but may start increasing for larger steps.
c     Choosing a large alfai allows several elements of  x  to
c     become feasible at the same time.
c
c     In the end, we take  alfa = alfaf  if x is feasible, or if
c     alfai > alfap (where  alfap  is the perturbed step from pass 1).
c     Otherwise,  we take  alfa = alfai.
c
c     Input parameters
c     ----------------
c     mnsmvl normally false, but true when we have discovered that
c            no feasible point exists, and we want to find a feasible 
c            point that minimizes the sum of violations
c     n      total number of variables
c     nclin  total number of linear constraints
c     istate state vector for constraints
c     bigalf defines what should be treated as an unbounded step.
c     bigbnd provides insurance for detecting unboundedness.
c            If alfa reaches a bound as large as bigbnd, it is
c            classed as an unbounded step.
c     pnorm  euclidean norm of  p
c     Anorm  norms of constraint (row-wise norms of A)
c     Ap     product of constraint matrix  A  and search vector   p
c     Ax     product of constraint matrix  A  and current point   x
c     bl     array of lower bounds
c     bu     array of upper bounds
c     featol is the array of current feasibility tolerances used by
c            cmsinf.  Typically in the range 0.5*tolx to 0.99*tolx,
c            where tolx is the featol specified by the user.
c     featlu the user's original feasibility tolerances
c     p      search direction vector
c     x      current iterate
c     epspt9 machine precision to the 0.9th power
c     tolinc (in common) is used to determine stepmn (see below),
c            the minimum positive step.
c     istate is set as follows:
c            istate(j) = -2  if a'x .lt. bl - featol
c                      = -1  if a'x .gt. bu + featol
c                      =  0  if a'x is not in the working set
c                      =  1  if a'x is in the working set at bl
c                      =  2  if a'x is in the working set at bu
c                      =  3  if a'x is in the working set (an equality)
c                      =  4  if x(j) is temporarily fixed.
c            values -2 and -1 do not occur once feasible.
c
c     Output Parameters
c     -----------------
c     hitlow  = true  if a lower bound restricted alfa.
c             = false otherwise.
c     move    = true  if  exact ge stepmn  (defined at end of code).
c     onbnd   = true  if  alfa = exact.  This means that the step  alfa
c                     moves x  exactly onto one of its constraints,
c                     namely  bound.
c             = false if the exact step would be too small
c                     ( exact .lt. stepmn ).
c               (with these definitions,  move = onbnd).
c     unbndd  = true  if alfa = bigalf.  Jhit may possibly be zero.
c               The parameters hitlow, move, onbnd, bound and exact
c               should not be used.
c     jhit    = the index (if any) such that constraint jhit reaches
c               a bound.
c     bound   = the bound value bl(jhit) or bu(jhit) corresponding
c               to hitlow.
c     exact   = the step that would take constraint jhit exactly onto
c               bound.
c     alfa    = an allowable, positive step.
c               if unbndd is true,  alfa = stepmx.
c               otherwise,          alfa = max( stepmn, exact ).
c     ndegen
c
c
c     Cmchzr is based on MINOS 5.2 routine m5chzr, which implements the
c     expand procedure to deal with degeneracy. The step alfaf is
c     chosen as in the two-pass approach of Paula Harris (1973), except
c     that this version insists on returning a positive step, alfa.
c     Two features make this possible:
c
c        1. Featol increases slightly each iteration.
c
c        2. The blocking constraint, when added to the working set,
c           retains the value   Ax(jhit) + alfa * Ap(jhit),
c           even if this is not exactly on the blocking bound.
c
c     For infeasible variables moving towards their bound, we require
c     the rate of change of the chosen constraint to be at least gamma
c     times as large as the biggest available.  This still gives us
c     freedom in pass 2.
c     gamma = 0.1 and 0.01 seemed to inhibit phase 1 somewhat.
c     gamma = 0.001 seems to be safe.
c
c     ==================================================================

      integer            jhitf, jhiti 
      
      logical            blockf, blocki

      double precision   alfai , atp   , atpmxi, atx   , biglow, bigupp, 
     1                   bound , exact , mxreal, stepmn, tolpiv
      double precision   zero  , one   , gamma

      parameter        ( zero  = 0.0d0, one = 1.0d0, gamma = 1.0d-3 )

      double precision   hdmcon

      external           hdmcon

c     ==================================================================
      
c     ... vectorization --
c     
c         the two passes are separated into separate subroutines
c         for convenience.  In addition, the two versions of the
c         second pass (differing on whether mnsmvl is true or false)
c         are represented by separate subroutines, since the compiler
c         was unable to vectorize the code when combined.
c         Within each loop, the special code for each of the relevant
c         states ( -2 = violated lower bound,
c                  -1 = violated upper bound
c                   0 = inactive )
c         is handled as a separate complete branch of an if-then-else.
c         this otherwise awkward and unnecessarily complex approach
c         enables the CRAY compiler to vectorize the loops.  
c         so much for compilers not needing to be massaged.

c     tolpiv is a tolerance to exclude negligible elements of a'p.

      biglow = - bigbnd
      bigupp =   bigbnd

      tolpiv = epspt9*pnorm
      mxreal = hdmcon (2)

c     ... it is important that mxreal is strictly bigger than bigalf

c     ------------------------------------------------------------------
c     First pass -- find steps to perturbed constraints, so that
c     alfap will be slightly larger than the true step.
c     In degenerate cases, this strategy gives us some freedom in the
c     second pass.  The general idea follows that described by P.M.J.
c     Harris, p.21 of Mathematical Programming 5, 1 (1973), 1--28.
c     ------------------------------------------------------------------

c       << lpchz1 >>
      call LPCHZ1 ( n, nclin,
     1              istate, bigalf, bigbnd, 
     3              alfap, 
     4              Anorm, Ap, Ax,
     5              bl, bu, featol, p, x,
     6              tolpiv, mxreal, 
     7              atpmxi )

c     --------------------------------------------------------------
c     ... at end of Pass one, we have set  alfap  to the shortest
c         step to the far side of the feasible region, and just
c         a hair beyond.
c         atpmxi  is the largest angle between the search direction
c         and a currently infeasible constraint for which the search
c         direction is in the direction of the feasible region.
c     --------------------------------------------------------------

c     --------------------------------------------------------------
c     ... Second pass.
c     
c         For feasible variables:
c             recompute steps without perturbation.
c             amongst constraints that are closer than alfap,
c             choose the one that makes the largest angle with
c             the search direction.
c     
c         For infeasible variables:
c             find the largest step 
c             subject to a'p being no smaller than gamma * max(a'p).
c     --------------------------------------------------------------

      atpmxi = gamma*atpmxi
      
      if  ( mnsmvl )  then

c          << lpchz3 >> 
         call LPCHZ3 ( n, nclin,
     1                 istate, bigalf, bigbnd, 
     2                 alfap, alfai, jhitf, jhiti,
     3                 Anorm, Ap, Ax,
     4                 bl, bu, featol, p, x,
     5                 tolpiv, mxreal,
     6                 atpmxi )
         
      else
         
c          << lpchz2 >> 
         call LPCHZ2 ( n, nclin,
     1                 istate, bigalf, bigbnd, 
     2                 alfap, alfai, jhitf, jhiti,
     3                 Anorm, Ap, Ax,
     4                 bl, bu, featol, p, x,
     5                 tolpiv, mxreal,
     6                 atpmxi )
         
      end if

c     ------------------------------------------------------
c     See if a feasible and/or infeasible constraint blocks.
c     ------------------------------------------------------

      blockf = jhitf .gt. 0
      blocki = jhiti .gt. 0
      unbndd = .not. ( blockf  .or.  blocki )

      if  ( unbndd) go to 900

      if  ( blockf )  then
         
c        ---------------------------------------------------------------
c        A constraint is hit which is currently feasible.
c        The corresponding step alfaf is not used, so no need to get it,
c        but we know that alfaf .le. alfap, the step from pass 1.
c        ---------------------------------------------------------------
         
         jhit = jhitf
         if  ( jhit .le. n )  then
            atp = p (jhit)
         else
            atp = Ap (jhit-n)
         end if
         hitlow = atp .lt. zero
         
      end if

c     If there is a choice between alfaf and alfai, it is probably best
c     to take alfai.  However, we can't if alfai is bigger than alfap.

      if  ( blocki  .and.  alfai .le. alfap )  then
         
c        --------------------------------------------------
c        An infeasible variable reaches its violated bound.
c        --------------------------------------------------
         
         jhit = jhiti
         if  ( jhit .le. n )  then
            atp = p (jhit)
         else
            atp = Ap (jhit-n)
         end if
         hitlow = atp .gt. zero
         
      end if

      if  ( jhit .le. n )  then
         atx = x (jhit)
      else
         atx = Ax (jhit-n)
      end if

c     ---------------------------------------------------------------
c     Try to step exactly onto bound, but make sure the exact step
c     is sufficiently positive.  (Exact will be alfaf or alfai.)
c     Since featol increases by  tolinc  each iteration, we know that
c     a step as large as  stepmn  (below) will not cause any feasible
c     variables to become infeasible (where feasibility is measured
c     by the current featol).
c     ---------------------------------------------------------------
      
      if  (  hitlow  )  then
         bound = bl (jhit)
      else
         bound = bu (jhit)
      end if

      unbndd = abs ( bound ) .ge. bigbnd
      if  ( unbndd) go to 900

      stepmn = tolinc * featlu (jhit)  /  abs ( atp )
      exact  = (bound - atx)           /      atp
      alfa   = max ( stepmn, exact )
      onbnd  = alfa  .eq. exact
      move   = exact .ge. stepmn
      if  ( .not. move)  then
         ndegen = ndegen + 1
      endif

      return
      
c     ----------
c     Unbounded.
c     ----------
      
  900 continue
      alfa   = bigalf
      move   = .true.
      onbnd  = .false.

      return

c     end of  LPCHZR (lpchzr).
      end
