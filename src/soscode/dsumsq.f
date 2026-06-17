      subroutine   DSUMSQ   ( n, x, incx, scale, sumsq )
c     Clean scaled sum-of-squares update:
c       scale**2*sumsq := scale**2*sumsq + sum_i x_i**2

      integer            n, incx
      double precision   x( * ), scale, sumsq
      integer            i, ix
      double precision   ax, ratio, zero, one
      parameter        ( zero = 0.0d0, one = 1.0d0 )
      intrinsic          abs

      ix = 1
      do i = 1, n
         if  ( x(ix) .ne. zero )  then
            ax = abs( x(ix) )
            if  ( scale .lt. ax )  then
               if  ( scale .eq. zero )  then
                  sumsq = one
               else
                  ratio = scale / ax
                  sumsq = one + sumsq*ratio*ratio
               endif
               scale = ax
            else
               ratio = ax / scale
               sumsq = sumsq + ratio*ratio
            endif
         endif
         ix = ix + incx
      enddo

      end
