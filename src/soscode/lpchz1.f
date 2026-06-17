      subroutine   LPCHZ1  ( n, nclin,
     1                       istate, bigalf, bigbnd, 
     3                       alfap, 
     4                       Anorm, Ap, Ax,
     5                       bl, bu, featol, p, x,
     6                       tolpiv, mxreal, 
     7                       atpmxi )
c     ==================================================================
c     ==================================================================
c     ====  lpchz1 / LPCHZ1 -- choose a constraint to add ...       ====
c     ====                     general case (infeasible point)      ====
c     ====                     Pass 1 of 2                          ====
c     ==================================================================
c     ==================================================================

      integer            n, nclin
      
      double precision   bigalf, bigbnd, alfap, tolpiv, mxreal, atpmxi
      
      integer            istate (n+nclin)
      
      double precision   bl (n+nclin), bu (n+nclin)
      
      double precision   featol (n+nclin) 
      
      double precision   Anorm (*), Ap (*), Ax (*)
      
      double precision   p (n), x (n)

c     ==================================================================
c
c     derived from qpopt version 1.0 cmchzr
c
c     special vectorized version of cmchzr
c     
c     last modification -- 25-June-1996
c
c         Original version written by PEG,  19-April 1988.
c         This version of  cmchzr  dated   6-Jul-1988.
c
c     LPCHZ1 / lpchz1  finds a step alfap such that the point
c         x + alfap*p reaches one of the linear constraints
c         (including bounds).
c
c     LPCHZ1 / lpchz1 is a vectorized implementation of only
c     the first of two passes from  cmchzr
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
c
c     alfap  step to the nearest bound
c     atpmxi largest angle for a step to the closest bound for
c            a currently violated constraint
c     
c     ==================================================================

      integer            j     
      
      double precision   atp   , atpabs, atpscd, 
     1                   atx   , biglow, bigupp, talfap

      double precision   zero  , one

      parameter        ( zero  = 0.0d0, one = 1.0d0 )

c     ==================================================================
      
c     ... vectorization -- loops of original cmchzr are split into
c         separate loops for simple bounds and linear constraints.
c         this otherwise awkward and unnecessarily complex approach
c         enables the CRAY compiler to vectorize the loops.  Further,
c         within each loop, the special code for each of the relevant
c         states ( -2 = violated lower bound,
c                  -1 = violated upper bound
c                   0 = inactive )
c         is handled as a separate complete branch of an if-then-else.
c         so much for compilers not needing to be massaged.

c     tolpiv is a tolerance to exclude negligible elements of a'p.

      biglow = - bigbnd
      bigupp =   bigbnd

c     ------------------------------------------------------------------
c     First pass -- find steps to perturbed constraints, so that
c     alfap will be slightly larger than the true step.
c     In degenerate cases, this strategy gives us some freedom in the
c     second pass.  The general idea follows that described by P.M.J.
c     Harris, p.21 of Mathematical Programming 5, 1 (1973), 1--28.
c     ------------------------------------------------------------------

      atpmxi = zero
      alfap  = bigalf

c     ----------
c     ... Pass 1
c     ----------
      
      do j = 1, n+nclin
         
         talfap = mxreal

         if  ( j .le. n )  then

c           -------------------------
c           ... simple bounds
c           -------------------------
            
            atx    = x (j)
            atp    = p (j)
            atpabs = abs ( atp )
            atpscd = atpabs
            
         else
         
c           -------------------------
c           ... general linear bounds
c           -------------------------
            
            atx    = Ax (j-n)
            atp    = Ap (j-n)
            atpabs = abs ( atp )
            atpscd = atpabs / ( one  +  Anorm (j-n) )

         endif
         
         if  ( istate (j) .eq. 0 )  then

c           -----------------------
c           ... inactive constraint
c           -----------------------
            
            if  (  atp .lt. -tolpiv  .and.  bl (j) .gt. biglow )  then
               
c              ---------------------------------------
c              ... a'x  is decreasing sufficiently
c                  rapidly and there is a lower bound.
c              ---------------------------------------

               talfap = ( atx - bl (j) + featol (j) ) / (-atp)

            else
     1      if  ( atp .gt. tolpiv  .and.  bu (j) .lt. bigupp )  then
               
c              ---------------------------------------
c              ... a'x  is increasing sufficiently
c                  rapidly and there is an upper bound.
c              ----------------------------------------
               
               talfap = ( bu (j) - atx + featol (j) )  / atp

            end if
               
         else
     1   if  ( istate (j) .eq. -1 )  then

c           -------------------------------------
c           ... upper bound is currently violated
c           -------------------------------------
            
            if  ( atp .lt. -tolpiv  .and.  bl (j) .gt. biglow )  then
               
c              -------------------------------------------
c              ... a'x  is decreasing sufficiently rapidly
c                  and there is a lower bound.
c              -------------------------------------------

               talfap = ( atx - bl (j) + featol (j) ) / (-atp)
               atpmxi = max ( atpmxi, atpscd )

            end if
               
         else
     1   if  ( istate (j) .eq. -2 )  then

c           -------------------------------------
c           ... lower bound is currently violated
c           -------------------------------------
            
            if  ( atp  .gt. tolpiv  .and.  bu (j) .lt. bigupp )  then
               
c              ---------------------------------------
c              ... a'x  is increasing sufficiently rapidly
c                  and there is an upper bound.
c              ---------------------------------------
               
               talfap = ( bu (j) - atx + featol (j) )  / atp
               atpmxi = max ( atpmxi, atpscd ) 

            end if

         endif
         
c        -------------------------------------------------------
c        ... test for new shorter step (not done in if-then-else
c            above in order to please the cray compiler)
c        -------------------------------------------------------

         alfap = min (alfap, talfap)
            
      enddo

      return

c     end of  LPCHZ1 (lpchz1).
      end
