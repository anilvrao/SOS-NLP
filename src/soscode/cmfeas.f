      subroutine   CMFEAS   ( n, nclin, istate,
     1                        bigbnd, nviol, jmax, errmax,
     2                        Ax, bl, bu, featol, x )
c     =================================================================
c     =================================================================
c     ==== CMFEAS /                                                ====
c     ==== cmfeas -- check feasibility against active constraints  ====
c     =================================================================
c     =================================================================

      integer            n, nclin, nviol, jmax

      double precision   bigbnd, errmax
      
      integer            istate (n+nclin)
      
      double precision   Ax (*), bl (n+nclin), bu (n+nclin)
      
      double precision   featol (n+nclin), x (n)

c     ==================================================================
c
c     derived from qpopt version 1.0
c     last modification -- 26-March-1996
c
c          Original version written by PEG,   April    1984.
c          This version of  cmfeas  dated  30-Jun-1988.
c
c     CMFEAS /cmfeas  checks the residuals of the constraints that are
c     believed to be feasible.  The number of constraints violated by 
c     more than featol is computed, along with the maximum constraint
c     violation.
c
c     ==================================================================

      integer            is, j
      
      double precision   biglow, bigupp, con, feasj, res

      double precision   zero 

      parameter        ( zero = 0.0d0 )

c     ==================================================================

      biglow = - bigbnd
      bigupp =   bigbnd

c     ==================================================================
c     Compute the number of constraints (nviol) violated by more than
c     featol and  the maximum constraint violation (errmax).
c     (The residual of a constraint in the working set is treated as if
c     it were an equality constraint fixed at that bound.)
c     ==================================================================

      nviol  = 0
      jmax   = 0
      errmax = zero

      do j = 1, n+nclin
         
         is     = istate(j)

         if  ( is .ge. 0 )  then
            
            feasj  = featol(j)

            if  ( j .le. n )  then
               con =  x(j)
            else
               con = Ax(j-n)
            end if

c           Check for constraint violations.

            if  ( bl(j) .gt. biglow )  then
               res    = bl(j) - con
               if  ( res .gt.   feasj  )  then
                  nviol  = nviol  + 1
                  go to 190
               end if
            end if

            if  ( bu(j) .lt. bigupp )  then
               res    = bu(j) - con
               if  ( res .lt. (-feasj) )  then
                  nviol  =   nviol + 1
                  res    = - res
                  go to 190
               end if
            end if

c           This constraint is satisfied,  but count a large residual
c           as a violation if the constraint is in the working set.

            res   = zero

            if  ( is .eq. 1 )  then
               res = abs( bl(j) - con )
            else
     1      if  ( is .eq. 2 )  then
               res = abs( bu(j) - con )
            else
     1      if  ( is .eq. 3 )  then
               res = abs( bu(j) - con )
            end if

            if  ( res .gt. feasj )  then
               nviol  = nviol  + 1
            endif

  190       continue
            if  ( res .gt. errmax )  then
               jmax   = j
               errmax = res
            end if
            
         end if
         
      enddo

c     end of CMFEAS (cmfeas)
      end         
