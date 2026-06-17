      subroutine   LPCHZ2  ( n, nclin,
     1                       istate, bigalf, bigbnd, 
     2                       alfap, alfai, jhitf, jhiti,
     3                       Anorm, Ap, Ax,
     4                       bl, bu, featol, p, x,
     5                       tolpiv, mxreal, 
     6                       atpmxi )
c     ==================================================================
c     ==================================================================
c     ====  lpchz2 / LPCHZ2 -- choose a constraint to add ...       ====
c     ====                     general case (infeasible point)      ====
c     ====                     Pass Two of Two                      ====
c     ====                     (normal version of Pass Two)         ====
c     ==================================================================
c     ==================================================================

      integer            n, nclin, jhitf, jhiti
      
      double precision   bigalf, bigbnd, alfap, alfai, tolpiv, mxreal,
     1                   atpmxi
      
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
c     last modification -- 26-July-1996
c
c         Original version written by PEG,  19-April 1988.
c         This version of  cmchzr  dated   6-Jul-1988.
c
c     LPCHZ2 / lpchz2  finds a step alfa such that the point
c         x + alfa*p reaches one of the linear constraints
c         (including bounds).
c
c     LPCHZ1 / lpchz1 is a vectorized implementation of only
c     the second of two passes from  cmchzr.  In particular,
c     it implements the second pass for the general case in
c     which  firstv/mnsmvl is false.  This corresponds in QPOPT
c     to the first phase search for a feasible point.
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
c     alfap  step to the nearest bound from Pass 1
c     atpmxi largest angle for a step to the closest bound for
c            a currently violated constraint
c
c     Output Parameters
c     -----------------
c     
c     jhitf   = the index (if any) such that constraint jhitf reaches
c               the farther bound in the direction given by p for any
c               constraint.
c     jhiti   = the index (if any) such that constraint jhitf reaches
c               the nearer bound in the direction given by p for
c               currently violated constraints.
c
c     ==================================================================

      integer            j     
      
      double precision   atp   , atpabs, atpmxf, atpscd, atx   , 
     1                   biglow, bigupp, res  , talfai

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

      biglow = featol(1)
      biglow = bigalf
      biglow = - bigbnd
      bigupp =   bigbnd

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

      alfai = zero

      atpmxf = zero
      jhitf  = 0
      jhiti  = 0

c     ----------
c     ... Pass 2
c     ----------
      
      do j = 1, n+nclin

         if  ( j .le. n )  then

c           -----------------
c           ... simple bounds
c           -----------------
            
            atx    = x (j)
            atp    = p (j)
            atpabs = abs ( atp )
            atpscd = atp

         else

c           -------------------------
c           ... general linear bounds
c           -------------------------
            
            atx    = Ax (j-n)
            atp    = Ap (j-n)
            atpabs = abs ( atp )
            atpscd = atp  / ( one  +  Anorm (j-n) )

         endif

         res    = mxreal
         talfai = - one

         if  ( istate (j) .eq. 0 )  then

c           -----------------------
c           ... inactive constraint
c           -----------------------
   
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

            end if

         else
     1   if  ( istate (j) .eq. -1 )  then

c           -------------------------------------
c           ... upper bound is currently violated
c           -------------------------------------
            
            if  ( atpscd .lt. -tolpiv )  then

c              ----------------------
c              ... a'x  is decreasing
c              ----------------------
               
               talfai  = ( atx - bu (j) ) / atpabs

               if  (  bl (j) .gt. biglow )  then

c                 --------------------------
c                 ... there is a lower bound
c                 --------------------------
               
                  res = atx - bl (j)

               endif

            end if

         else
     1   if  ( istate (j) .eq. -2 )  then

c           -------------------------------------
c           ... lower bound is currently violated
c           -------------------------------------
            
            if  ( atpscd  .gt. tolpiv  )  then
               
c              ----------------------
c              ... a'x  is increasing
c              ----------------------

               talfai = ( bl (j) - atx ) / atpabs

               if  ( bu (j) .lt. bigupp )  then
                  
c                 ---------------------------
c                 ... there is an upper bound
c                 ---------------------------
                  
                  res = bu (j) - atx
                  
               endif
               
            end if

         endif
         
c        ----------------------------------------------------------
c        ... check for the largest angle with  p  for all steps
c            that satisfy the perturbed feasible bound from pass 1.
c        ----------------------------------------------------------
            
         if  ( res          .le. alfap*atpabs .and.
     1         abs (atpscd) .gt. atpmxf             )  then
            atpmxf = abs (atpscd)
            jhitf  = j
         endif

c        -------------------------
c        ... test for bigger alfai
c        -------------------------

         if  ( abs (atpscd) .ge. atpmxi .and.
     1         talfai       .gt. alfai )  then
            
            alfai  = talfai
            jhiti  = j
            
         end if

      enddo

      return
      
c     end of  LPCHZ2 (lpchz2).
      end
