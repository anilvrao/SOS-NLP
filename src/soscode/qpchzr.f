      subroutine   QPCHZR  ( firstv, n, nclin,
     1                       istate, bigalf, bigbnd, pnorm,
     2                       hitlow, move, onbnd, unbndd,
     3                       alfa, alfap, jhit,
     4                       anorm, Ap, Ax,
     5                       bl, bu, featol, featlu, p, x,
     6                       epspt9, tolinc, ndegen )
c     ==================================================================
c     ==================================================================
c     ====  qpchzr / QPCHZR -- choose a constraint to add ...       ====
c     ====                     when known to be feasible            ====
c     ==================================================================
c     ==================================================================

      integer            n, nclin, jhit, ndegen
      
      logical            firstv, hitlow, move, onbnd, unbndd

      double precision   bigalf, bigbnd, pnorm, alfa, alfap, epspt9,
     1                   tolinc
      
      integer            istate(n+nclin)
      
      double precision   bl(n+nclin), bu(n+nclin)
      
      double precision   featol(n+nclin), featlu(n+nclin)
      
      double precision   anorm(*), Ap(*), Ax(*)
      
      double precision   p(n), x(n)

c     ==================================================================
c
c     derived from qpopt version 1.0 cmchzr
c
c     special version for phase II of QP, during which all points
c     are known to be feasible
c     
c     last modification --25-June-1996
c
c         Original version of cmchzr written by PEG,  19-April 1988.
c         This version of  cmchzr  dated   6-Jul-1988.
c
c     QPCHZR / qpchzr  finds a step alfa such that the point
c         x + alfa*p reaches one of the linear constraints
c         (including bounds).
c
c     In this version  x is always feasible.
c
c     alfaf = the maximum step that can be taken without violating
c              one of the constraints that are currently satisfied.
c
c     Input parameters
c     ----------------
c     bigalf defines what should be treated as an unbounded step.
c     bigbnd provides insurance for detecting unboundedness.
c            If alfa reaches a bound as large as bigbnd, it is
c            classed as an unbounded step.
c     featol is the array of current feasibility tolerances used by
c            cmsinf.  Typically in the range 0.5*tolx to 0.99*tolx,
c            where tolx is the featol specified by the user.
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
c     bl     the lower bounds on the variables.
c     bu     the upper bounds on ditto.
c     x      the values of       ditto.
c     p      the search direction.
c
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
c
c
c     Qpchzr is based on MINOS 5.2 routine m5chzr, which implements the
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
c     ==================================================================

      integer            j     
      
      double precision   atp   , atpabs, atpmxf, atpscd, atx   , 
     1                   biglow, bigupp, bound , mxreal, exact ,
     2                   res   , stepmn, talfap, tolpiv

      double precision   zero  , one

      parameter        ( zero  = 0.0d0, one = 1.0d0 )

      double precision   hdmcon

      external           hdmcon

c     ==================================================================

      hitlow = firstv
      
c     tolpiv is a tolerance to exclude negligible elements of a'p.

      biglow = - bigbnd
      bigupp =   bigbnd

      tolpiv = epspt9*pnorm

c     ------------------------------------------------------------------
c     First pass -- find steps to perturbed constraints, so that
c     alfap will be slightly larger than the true step.
c     In degenerate cases, this strategy gives us some freedom in the
c     second pass.  The general idea follows that described by P.M.J.
c     Harris, p.21 of Mathematical Programming 5, 1 (1973), 1--28.
c     ------------------------------------------------------------------
c
c     ... vectorization -- reduction test for shortest step
c         (largest angle in pass 2) is moved out of the inner
c         if-then-else  construction and is done for all bounds,
c         not just the inactive ones.  (The test for active bounds
c         is set artificially to fail.)  This otherwise awkward and
c         unnecessarily complex approach enables the CRAY compiler to
c         vectorize the loops
     
      alfap  = bigalf
      mxreal = hdmcon (2)

c     ... it is important that mxreal is strictly bigger than bigalf

      do j = 1, n+nclin

         if  ( istate (j) .eq. 0 )  then

            if  ( j .le. n )  then

c              ... simple bounds
               
               atx    = x (j)
               atp    = p (j)
               atpscd = atp
               
            else
               
c              ... general linear constraints
               
               atx    = Ax (j-n)
               atp    = Ap (j-n)
               atpscd = atp / ( one  +  anorm (j-n) )
               
            endif
            
            if  ( atpscd .lt. -tolpiv  .and.  bl (j) .gt. biglow )  then
               
c              -------------------------------------------------
c              ... a'x  is decreasing and there is a lower bound
c              -------------------------------------------------
               
               talfap = ( atx - bl (j) + featol (j) ) / (-atp)

            else
     1      if  ( atpscd .gt. tolpiv  .and.  bu (j) .lt. bigupp )  then
               
c              --------------------------------------------------
c              ... a'x  is increasing and there is an upper bound
c              --------------------------------------------------
               
               talfap = ( bu (j) - atx + featol (j) )  / atp

            else

c              -----------------------------------------------------
c              either this constraint appears to be constant along p
c              or there is no bound on one side.  In either case it
c              is not used to compute the step.
c              -----------------------------------------------------
               
               talfap = mxreal
               
            end if

c           ---------------------------------------------------------
c           ... test for shortest step to the feasible tolerance just
c               beyond the bound.  
c           ---------------------------------------------------------

            alfap = min (alfap, talfap)
            
         end if
         
      enddo

c     ------------------------------------------------------------------
c     Second pass.
c     Recompute steps without perturbation.
c     amongst constraints that are closer than alfap, choose the one
c     That makes the largest angle with the search direction.
c     ------------------------------------------------------------------

      atpmxf = zero
      jhit   = 0

      do j = 1, n+nclin

         if  ( istate (j) .eq. 0 )  then

            if  ( j .le. n )  then

c              ... simple bounds
               
               atx    = x (j)
               atp    = p (j)
               atpabs = abs (atp)
               atpscd = atp

            else

c              ... general linear constraints

               atx    = Ax (j-n)
               atp    = Ap (j-n)
               atpabs = abs (atp)
               atpscd = atp / ( one  +  anorm (j-n) )

            endif

            if  ( atpscd .lt. -tolpiv  .and.  bl (j) .gt. biglow )  then
               
c              -------------------------------------------------
c              ... a'x  is decreasing and there is a lower bound
c              -------------------------------------------------
               
               res = atx - bl (j)

            else
     1      if  ( atpscd .gt. tolpiv  .and.  bu (j) .lt. bigupp )  then
               
c              --------------------------------------------------
c              ... a'x  is increasing and there is an upper bound
c              --------------------------------------------------
               
               res = bu (j) - atx

            else
               
c              ------------------------------------------------------
c              ... either this is an descent direction for this
c                  constraint, this constraint appears to be constant 
c                  along p or there is no bound on one side.  In any
c                  case it is not used to compute the step.
c              -----------------------------------------------------
               
              res = mxreal

            endif

c           -----------------------------------------------------
c           ... test for bound that makes the largest angle with
c               with the search direction, while having a step no
c               larger than the minimum "almost feasible" step
c               from Pass 1.
c           -----------------------------------------------------

            if  ( res          .le. alfap*atpabs .and.
     1            abs (atpscd) .gt. atpmxf             )  then

               atpmxf = abs (atpscd)
               jhit   = j
               
            endif
            
         endif
         
      enddo

c     ------------------------------------------------------
c     See if a feasible and/or infeasible constraint blocks.
c     ------------------------------------------------------

      unbndd = jhit .eq. 0

      if  ( unbndd) go to 900

c     ---------------------------------------------------------------
c     A constraint is hit which is currently feasible.
c     The corresponding step alfaf is not used, so no need to get it,
c     but we know that alfaf .le. alfap, the step from pass 1.
c     ---------------------------------------------------------------
      
      if  ( jhit .le. n )  then
         atp = p (jhit)
         atx = x (jhit)
      else
         atp = Ap (jhit-n)
         atx = Ax (jhit-n)
      end if
      hitlow = atp .lt. zero
         
c     ---------------------------------------------------------------
c     Try to step exactly onto bound, but make sure the exact step
c     is sufficiently positive.
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

      stepmn = tolinc * featlu (jhit)  /  abs (atp)
      exact  =     (bound - atx)       /     atp
      alfa   = max ( stepmn, exact )
      onbnd  = alfa  .eq. exact
      move   = exact .ge. stepmn
      if  ( .not. move )  then
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

c     end of  QPCHZR (qpchzr).
      end
