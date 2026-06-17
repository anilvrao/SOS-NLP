      subroutine  DROTGC   ( x, y, cs, sn )
c     ==================================================================
c     ==================================================================
c     ====  DROTGC /                                                ====
c     ====  drotgc  -- generate givens rotation (cautiously)        ====
c     ==================================================================
c     ==================================================================

      double precision   x, y, cs, sn

c     ==================================================================

c     derived from qpopt version 1.0
c     last modification -- 26-March-1996
c
c          Systems Optimization Laboratory, Stanford University.
c          Original version dated January 1982.
c          F77 version dated 28-June-1986.
c          This version of DROT3G dated 28-June-1986.

c
c     Note: drotgc/drot3g is different from the Nag routine f06baf.
c     
c     drotgc (DROTGC)  generates a plane rotation that reduces the
c     vector (X, Y) to the vector (A, 0),
c     where A is defined as follows...
c     
c          If both X and Y are negligibly small, or
c          if Y is negligible relative to Y,
c          then  A = X,  and the identity rotation is returned.
c     
c          If X is negligible relative to Y,
c          then  A = Y,  and the swap rotation is returned.
c     
c          Otherwise,  A = sign(X) * sqrt( X**2 + Y**2 ).
c     
c     In all cases,  X and Y are overwritten by A and 0,  and CS will
c     lie in the closed interval (0, 1).  Also,  the absolute value of
c     CS and SN (if nonzero) will be no less than the machine precision,
c     EPS.
c     
c     drotgc (DROTGC) DROT3G  guards against overflow and underflow.
c     It is assumed that  FLMIN .lt. EPS**2  (i.e.  RTMIN .lt. EPS).

c     ==================================================================

      double precision   a, b, eps, rtmin
      
      logical            first
      
      intrinsic          abs, max, sqrt

      double precision   zero, one
      
      parameter        ( zero = 0.0d0, one = 1.0d0 )

      save               first , eps   , rtmin
      
      double precision   hdmcon

      external           hdmcon

      data               first / .true. /

c     ==================================================================

      
      if  ( first )  then
         first = .false.
         eps    = hdmcon (5)
         rtmin  = sqrt ( hdmcon (4) )
      end if

      if  ( y .eq. zero )  then
         
         cs = one
         sn = zero
         
      else
     1if  ( x .eq. zero )  then

         cs = zero
         sn = one
         x  = y
         
      else

         a      = abs(x)
         b      = abs(y)
         
         if  ( max(a,b) .le. rtmin )  then
            
            cs = one
            sn = zero
            
         else
            
            if  ( a .ge. b )  then
               
               if  ( b .le. eps*a )  then
                  cs = one
                  sn = zero
                  go to 900
               else
                  a  = a * sqrt( one + (b/a)**2 )
               end if
               
            else
               
               if  ( a .le. eps*b )  then
                  cs = zero
                  sn = one
                  x  = y
                  go to 900
               else
                  a  = b * sqrt( one + (a/b)**2 )
               end if
               
            end if
            
            if  ( x .lt. zero)  then
               a = - a
            endif
            
            cs = x/a
            sn = y/a
            x  = a
            
         end if
         
      end if

  900 continue
      y  = zero

      
c     end of  DROTGC (drotgc).
      end
