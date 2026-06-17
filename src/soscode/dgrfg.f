      subroutine   DGRFG    ( n, alpha, x, incx, tol, zeta )
c     ==================================================================
c     ==================================================================
c     ====  DGRFG  / f06frf /                                       ====
c     ====  dgrfg -- generate generalized Householder reflection    ====
c     ==================================================================
c     ==================================================================

c     .. Scalar Arguments ..
      
      double precision   alpha, tol, zeta
      integer            incx, n
      
c     .. Array Arguments ..
      
      double precision   x ( * )
      
c     ==================================================================
c
c     derived from qpopt version 1.0 f06frf
c     last modification -- 17-April-1996
c        Modified by PEG from NAG F06FRF 9/25/88.
c
c     dgrfg (DGRFG ) generates a generalized Householder reflection
c     such that
c
c     P*( alpha ) = ( beta ),   P'*P = I.
c       (   x   )   (   0  )
c
c  P is given in the form
c
c     P = I - ( zeta )*( zeta  z' ),
c             (   z  )
c
c  where z is an n element vector and zeta is a scalar that satisfies
c
c     1.0 .le. zeta .le. sqrt( 2.0 ).
c
c  zeta is returned in ZETA unless x is such that
c
c     max( abs( x( i ) ) ) .le. max( eps*abs( alpha ), tol )
c
c  where eps is the relative machine precision and tol is the user
c  supplied value TOL, in which case ZETA is returned as 0.0 and P can
c  be taken to be the unit matrix.
c
c  beta is overwritten on alpha and z is overwritten on x.
c  the routine may be called with  n = 0  and advantage is taken of the
c  case where  n = 1.
c
c
c  Nag Fortran 77 O( n ) basic linear algebra routine.
c
c  -- Written on 30-August-1984.
c     Sven Hammarling, Nag Central Office.
c     This version dated 28-September-1984.
      
c     ==================================================================

      double precision   one        , zero

      parameter        ( one = 1.0d0, zero = 0.0d0 )

      double precision   beta, eps, scale, ssq, frac

      logical            first

      double precision   hdmcon

      external           DSUMSQ, hdmcon

      intrinsic          abs, max, sign, sqrt

      save               eps, first

      data               first/ .true. /

c     ==================================================================
      
      if  ( n .lt. 1 )  then
         
         zeta = zero
         
      else
     1if  ( ( n .eq. 1 )  .and.  ( x ( 1 ) .eq. zero ) )  then

         zeta = zero
         
      else

         if  ( first )  then
            first = .false.
            eps   =  hdmcon (5)
         end if
c        
c        Treat case where P is a 2 by 2 matrix specially.
c        
         if  ( n .eq. 1 )  then
c           
c           Deal with cases where  ALPHA = zero  and
c           abs( X( 1 ) ) .le. max( EPS*abs( ALPHA ), TOL )  first.
c           
            if  ( alpha .eq. zero )  then
               
               zeta   =  one
               alpha  =  abs ( x( 1 ) )
               x( 1 ) = -sign( one, x( 1 ) )
               
            else
     1      if  ( abs( x( 1 ) ) .le. max( eps*abs(alpha), tol ) )  then

               zeta   =  zero
               
            else
               
               if  ( abs( alpha ) .ge. abs( x( 1 ) ) )  then
                  beta = abs( alpha ) *sqrt( 1 + ( x( 1 )/alpha )**2 )
               else
                  beta = abs( x( 1 ) )*sqrt( 1 + ( alpha/x( 1 ) )**2 )
               end if
               
               zeta = sqrt( ( abs( alpha ) + beta )/beta )
               if  ( alpha .ge. zero )  then
                  beta = -beta
               endif
               x( 1 ) = -x( 1 )/( zeta*beta )
               alpha  = beta
               
            end if
            
         else
c           
c           Now P is larger than 2 by 2.
c           
            ssq   = one
            scale = zero
            
c             << dsumsq >>
            call DSUMSQ ( n, x, incx, scale, ssq )
c           
c           Treat cases where  SCALE = zero,
c           SCALE .le. max( EPS*abs( ALPHA ), TOL )  and
c           ALPHA = zero  specially.
c           Note that  SCALE = max( abs( X( i ) ) ).
c           
            if  ( ( scale .eq. zero ) .or. 
     1            ( scale .le. max ( eps*abs(alpha), tol ) ) )  then

               zeta  = zero
               
            else
     1      if  ( alpha .eq. zero )  then
               
               zeta  = one
               alpha = scale*sqrt( ssq )
               frac = -1/alpha
               do i=0,n-1
                 x(1+i*incx) = frac*x(1+i*incx)
               enddo
               
            else
               
               if  ( scale .lt. abs( alpha ) )  then
                  beta = abs( alpha )*sqrt( 1 + ssq*( scale/alpha )**2 )
               else
                  beta = scale       *sqrt( ssq +   ( alpha/scale )**2 )
               end if
               
               zeta = sqrt( ( beta + abs( alpha ) )/beta )
               if  ( alpha .gt. zero )  then
                  beta = -beta
               endif
               frac = -1/( zeta*beta )
               do i=0,n-1
                 x(1+i*incx) = frac*x(1+i*incx)
               enddo
               alpha = beta
               
            end if
            
         end if
         
      end if

c     end of DGRFG   / f06frf / dgrfg
      
      end
