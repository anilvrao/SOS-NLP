      subroutine   DRUTHS   ( side0, n, k1, k2, c, s, a, lda )
c     ==================================================================
c     ==================================================================
c     ====  DRUTHS / f06qvf /                                       ====
c     ====  druths -- rotate upper triangular into hessenberg       ====
c     ==================================================================
c     ==================================================================
      
      integer            k1, k2, lda, n
      
      character(len=1)   side, side0

      double precision   a( lda, * ), c( * ), s( * )

c     ==================================================================
c     derived from qpopt version 1.0 f06qvf
c     last modification -- 26-March-1996
c          originally (NAG f06qvf) mark 13 release. nag copyright 1988.
      
C     DRUTHS / F06QVF / DRUTHS applies a  given sequence  of  plane
c     rotations to either the left,  or the right,  of the  n by n
c     upper triangular matrix  U,  to transform U to an  upper
c     Hessenberg matrix. The rotations are applied in planes k1 to k2.
C
C  The upper Hessenberg matrix, H, is formed as
C
C     H = P*U,    when   SIDE = 'L' or 'l',  (  Left-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 )
C
C  and is formed as
C
C     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
C
C  where P is an orthogonal matrix of the form
C
C     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
C
C  P( k ) being a plane rotation matrix for the  ( k, k + 1 ) plane. The
C  cosine and sine that define P( k ), k = k1, k1 + 1, ..., k2 - 1, must
C  be  supplied  in  c( k )  and  s( k )  respectively.  The  two by two
C  rotation part of P( k ), R( k ), is assumed to have the form
C
C     R( k ) = (  c( k )  s( k ) ).
C              ( -s( k )  c( k ) )
C
C  The matrix  U must be supplied in the n by n leading upper triangular
C  part of the array  A, and this is overwritten by the upper triangular
C  part of  H.
C
C  The  sub-diagonal elements of  H, h( k + 1, k ),  are returned in the
C  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
C
C  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
C  greater than n then an immediate return is effected.
C
C
C  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
C
C  -- Written on 13-January-1986.
C     Sven Hammarling, Nag Central Office.
C
C
c     ==================================================================
      
      double precision   one, zero
      
      parameter        ( one = 1.0d0, zero = 0.0d0 )

      double precision   aij, ctemp, stemp, temp

      integer            i, j

      intrinsic          min

c     ==================================================================
      side = side0(1:1)
      
      if  ( ( min( n, k1 ) .lt. 1 )  .or.
     1      ( k2 .le. k1          )  .or.
     2      ( k2.gt.n )                   )  then
         return
      endif
      
      if  ( ( side .eq. 'L' )  .or.  ( side .eq. 'l' ) )  then
c
c        apply the plane rotations to columns n back to k1.
c
         do j = n, k1, -1
            
            if  ( j.ge.k2 )  then
               
               aij = a( k2, j )
               
            else

c              form  the  additional sub-diagonal element  h( j + 1, j )
c              and store it in s( j ).

               aij = c( j )*a( j, j )
               s( j ) = -s( j )*a( j, j )
            end if
            
            do i = min( k2, j ) - 1, k1, -1
               temp = a( i, j )
               a( i + 1, j ) = c( i )*aij - s( i )*temp
               aij = s( i )*aij + c( i )*temp
            enddo
            a( k1, j ) = aij
            
         enddo
         
      else
     1if  ( ( side .eq. 'R' )  .or.  ( side .eq. 'r' ) )  then

c        apply  the  plane rotations  to  columns  k1  up to  ( k2 - 1 )
c        and  form   the   additional  sub-diagonal  elements,   storing
c        h( j + 1, j ) in s( j ).

         do j = k1, k2 - 1
            
            if  ( ( c( j ).ne.one )  .or.  ( s( j ).ne.zero ) )  then
               stemp = s( j )
               ctemp = c( j )
               do i = 1, j
                  temp = a( i, j + 1 )
                  a( i, j + 1 ) = ctemp*temp - stemp*a( i, j )
                  a( i, j ) = stemp*temp + ctemp*a( i, j )
               enddo
               s( j ) = stemp*a( j + 1, j + 1 )
               a( j + 1, j + 1 ) = ctemp*a( j + 1, j + 1 )
            end if
         enddo
         
      end if

      return

c     end of DRUTHS / f06qvf/ druths. 

      end
