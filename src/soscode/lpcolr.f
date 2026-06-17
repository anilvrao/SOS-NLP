      subroutine   LPCOLR   ( nZr, ldR, R, Rzz )
c     Clean implementation for SOS: set the last column of R(1:nZr,1:nZr)
c     to Rzz times the last coordinate vector.

      integer            nZr, ldR, i
      double precision   Rzz, R( ldR, * ), zero
      parameter        ( zero = 0.0d0 )

      if  ( nZr .le. 0 )  return

      do i = 1, nZr
         R( i, nZr ) = zero
      enddo
      R( nZr, nZr ) = Rzz

      end
