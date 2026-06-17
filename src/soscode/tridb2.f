      subroutine  TRIDB2   ( n, relprc, d, e, e2, leftvl, rghtvl ) 
c
c     ==================================================================
c     ==================================================================
c     ====  TRIDB2 /                                                ====
c     ====  tridb2 -- find the two extreme eigenvalues of a         ====
c     ====            symmetric tridiagonal matrix                  ====
c     ==================================================================
c     ==================================================================
c     ==================================================================
      
      integer           n

      double precision  relprc, leftvl, rghtvl

      double precision  d (n), e (n), e2 (n)
      
      integer           count , i

      double precision  epslon, ei    , eip1 , lowerb, lwrglm, maxd  , 
     1                  mind  , maxe  , range, releps, shift , test1, 
     2                  test2 , upperb, uprglm

      double precision  zero, one, half

      parameter       ( zero = 0.0d0, one = 1.0d0, half = 5.0d-1 )

      integer           hdqpe3
      double precision  hdmcon
      external          hdmcon
      
c     ==================================================================

c     last modification  16-August-1999
      
c     this subroutine is derived from EISPACK tridib, which is a
c     translation of the algol procedure bisect, 
c     num. math. 9, 386-393(1967) by barth, martin, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 249-256(1971) .
c
c     this subroutine finds the largest and smallest eigenvalues of
c     a tridiagonal symmetric matrix using bisection.
c
c     on input
c
c        n is the order of the matrix.
c
c        relprc is an relative error tolerance for the computed
c          eigenvalues.  relprc should be at least as large
c          as machine relative precision.  
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is arbitrary.
c
c     on output
c
c        relprc is unaltered.
c
c        d and e are unaltered.
c
c        elements of e2, corresponding to elements of e regarded
c          as negligible, have been replaced by zero causing the
c          matrix to split into a direct sum of submatrices.
c          e2(1) is also set to zero.
c
c        leftvl and rghtvl are estimates of the smallest and largest
c          eigenvalues, accurate to within 'relerr' relative error
c
c     ------------------------------------------------------------------


      if  ( n .eq. 1  )  then
         leftvl = d (1) 
         rghtvl = d (1) 
         return
      endif

      epslon = hdmcon (5)
      
      if  ( relprc .lt. epslon )  then
         releps = 2*epslon
      else
         releps = relprc
      endif

c     ----------------------------------------------------
c     ... find gerschgorin bounds for matrix and also find
c         find smallest and largest diagonals, to use as
c         initial bounds on intervals.  Also reset  e2
c         for sturm sequence scheme.
c     ----------------------------------------------------
      
      lwrglm = d (1) 
      uprglm = d (1)
      maxd   = d (1)
      mind   = d (1)
      maxe   = zero
      eip1   = zero

      do i = 1, n

         maxd = max ( maxd, d (i) )
         mind = min ( mind, d (i) )
         maxe = max ( maxe, abs (e (i)) )
         
         ei = eip1
         if  ( i .lt. n )  then
            eip1 = abs (e (i+1) )
         else
            eip1 = zero
         endif
         
         lwrglm = min (d (i) - (ei+eip1), lwrglm) 
         uprglm = max (d (i) + (ei+eip1), uprglm)
         
         if  ( i .ne. 1 )  then
            test1 = abs (d (i) ) + abs (d (i-1) ) 
            test2 = test1 + abs (e (i) ) 
            if  ( test2 .le. test1 ) then
               e2 (i) = zero
            endif
         else
            e2 (1) = zero
         endif

         range = max ( abs (lwrglm), abs (uprglm) )
         
      enddo

      
      if  ( maxe .eq. zero )  then

c        ------------------------------------
c        ... special case for diagonal matrix
c        ------------------------------------

         leftvl = mind
         rghtvl = maxd
         return

      else

c        ... find smallest eigenvalue

         lowerb = lwrglm
         upperb = mind


 200     continue
         if  ( upperb - lowerb .gt. releps * range )  then

            shift = lowerb + (upperb - lowerb) * half

c           ... find number of eigenvalues beneath shift

            count = hdqpe3 (  n, d, e, e2, shift, epslon )


            if  ( count .ge. 1 )  then
               upperb = shift
            else
               lowerb = shift
            endif
            go to 200

         else

c           ... cautiously -- return lower bound, not shift
            
            leftvl = lowerb


         endif

c        ... find largest eigenvalue

         lowerb = maxd
         upperb = uprglm


  300     continue
          if  ( upperb - lowerb .gt. releps * range )  then

            shift = lowerb + (upperb - lowerb) * half

            count = hdqpe3 (  n, d, e, e2, shift, epslon )


            if  ( count .eq. n )  then
               upperb = shift
            else
               lowerb = shift
            endif
            go to 300

         else

c           ... cautiously -- return upper bound, not shift
            
            rghtvl = upperb


         endif

      endif

      return

c     end of TRIDB2 (tridb2)
      end
