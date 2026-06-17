      integer function  hdqpe3   ( n, d, e, e2, shift, epslon )
c     ==================================================================
c     ==================================================================
c     ====  hdqpe3 --- sturm sequence count for tridiagonal         ====
c     ==================================================================
c     ==================================================================

      integer           n

      double precision  shift, epslon

      double precision  d (n), e (n), e2 (n)

      integer           i

      double precision  zero, one, qsubi, v

      parameter       ( zero = 0.0d0, one = 1.0d0 )

c     ==================================================================

c     ... last modification  16-August-1999
      
c     ----------------------------------------------------------------
c     ... see Wilkinson, AEP, pg. 300 ff for details of the algorithm,
c         which is based on computing determinants of leading principal
c         minors of the tridiagonal matrix  T - shift*I.
c         
c         the term qsubi is the ratio of the determinant of the ith order
c         leading principal minor to the determinant of the (i-1)st
c         order leading principal minor.
c
c         the number of eigenvalues below "shift" is equal to the
c         number of sign changes in the determinant, hence the number
c         of times the ratio  qsubi  is negative
c
c         "v" is a temporary variable as in EISPACK, to handle the 
c         special case of one of the determinants being exactly zero.
c     ----------------------------------------------------------------
      
      hdqpe3 = 0
     
      do i = 1, n

c        ... check for beginning of new submatrix; note that
c            e2 (1) was set to zero by TRIDB2/tridb2, so qsubi
c            is always defined
         
         if  ( e2 (i) .eq. zero )  then
            qsubi = one
         endif

c        ... special case for zero determinant
         
         if  ( qsubi .eq. zero )   then
            
            if ( e2 (i) .ne. zero )  then
               v = abs ( e (i) ) / epslon
            else
               v = zero
            endif
            
         else
            
            v = e2 (i) / qsubi

         endif
         
         qsubi = d (i) - shift - v
         
         if  ( qsubi .lt. zero )  then
            hdqpe3 = hdqpe3 + 1
         endif
         
      enddo

      return

      end
