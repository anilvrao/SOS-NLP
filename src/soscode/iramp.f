      subroutine  iramp(n,sx,incx)
 
c
c     set the i-th entries of a vector, x, to i.
c
c     ------------------------------------------------------------------
c
      integer              sx(*)
      integer              i,incx,ix,n
 
c     ------------------------------------------------------------------
 
      if(n.le.0)return
 
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      do i = 1,n
        sx(ix) = i
        ix = ix + incx
      enddo
 
c     ------------------------------------------------------------------
 
      return
      end
