      subroutine   DCOND    ( n, x, incx, xmax, xmin )
c     Clean implementation for SOS: return the largest and smallest
c     absolute entries in a strided vector.

      integer            incx, n
      double precision   x( * ), xmax, xmin
      integer            i, ix
      double precision   ax, zero
      parameter        ( zero = 0.0d0 )
      intrinsic          abs

      if  ( n .lt. 1 )  then
         xmax = zero
         xmin = zero
         return
      endif

      ix   = 1
      xmax = abs( x(ix) )
      xmin = xmax

      do i = 2, n
         ix = ix + incx
         ax = abs( x(ix) )
         if  ( ax .gt. xmax )  xmax = ax
         if  ( ax .lt. xmin )  xmin = ax
      enddo

      end
