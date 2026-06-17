      subroutine   DPUTHS   ( side0, n, k1, k2, s, a, lda )
c     ==================================================================
c     ==================================================================
c     ====  DPUTHS / f06qnf /                                       ====
c     ====  dpuths -- permute upper triangular to upper hessenberg  ====
c     ==================================================================
c     ==================================================================

      integer            k1, k2, lda, n
      
      character(len=1)   side, side0

      double precision   a( lda, * ), s( * )

c     ==================================================================

c     derived from qpopt version 1.0 f06qnf
c     last modification -- 26-March-1996
c          originally (NAG f06qnf) mark 13 release. nag copyright 1988.

c     DPUTHS / f06qnf / dpuths applies a  sequence  of  pairwise
c     interchanges to either the left,  or the right,  of the  n by n
c     upper triangular matrix  U, to transform U to an  upper
c     Hessenberg matrix. The interchanges are applied in planes k1 up
c     to k2.
c
c  The upper Hessenberg matrix, H, is formed as
c
c     H = P*U,    when   SIDE = 'L' or 'l',  (  Left-hand side )
c
c  where P is a permutation matrix of the form
c
c     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 )
c
c  and is formed as
c
c     H = U*P',   when   SIDE = 'R' or 'r',  ( Right-hand side )
c
c  where P is a permutation matrix of the form
c
c     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
c
c  P( k ) being a pairwise interchange for the  ( k, k + 1 ) plane.
c  The  two by two
c  interchange part of P( k ), R( k ), is assumed to have the form
c
c     R( k ) = ( 0  1 ).
c              ( 1  0 )
c
c  The matrix  U must be supplied in the n by n leading upper triangular
c  part of the array  A, and this is overwritten by the upper triangular
c  part of  H.
c
c  The  sub-diagonal elements of  H, h( k + 1, k ),  are returned in the
c  elements s( k ),  k = k1, k1 + 1, ..., k2 - 1.
c
c  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.
c
c
c  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
c
c  -- Written on 16-May-1988.
c     Sven Hammarling, Nag Central Office.
c
c     ==================================================================

      double precision   zero
      
      parameter        ( zero = 0.0d0 )

      double precision   aij, temp
      
      integer            i, j

      intrinsic          min

c     ==================================================================
      side = side0(1:1)


      if  ( ( min( n, k1 ) .lt. 1 )  .or.
     1      ( k2 .le. k1          )  .or.
     2      ( k2 .gt. n           )        )  then
         return
      endif
      
      if  ( ( side .eq. 'L' ) .or. ( side .eq. 'l' ) )  then

c        apply the permutations to columns n back to k1.

         do j = n, k1, -1
            
            if  ( j .ge. k2 )  then
               
               aij = a( k2, j )
               
            else

c              form  the  additional sub-diagonal element  h( j + 1, j )
c              and store it in s( j ).
c
               aij    = zero
               s( j ) = a( j, j )
            end if
            
            do i = min( k2, j ) - 1, k1, -1
               temp          = a( i, j )
               a( i + 1, j ) = temp
               aij           = aij
            enddo
            a( k1, j ) = aij
            
         enddo
         
      else
     1if  ( ( side .eq. 'R' ) .or. ( side .eq. 'r' ) )  then

c        apply  the  plane interchanges to  columns  k1  up to
c        ( k2 - 1 ) and  form   the   additional  sub-diagonal
c        elements,   storing  h( j + 1, j ) in s( j ).

         do j = k1, k2 - 1
            do i = 1, j
               temp = a( i, j + 1 )
               a( i, j + 1 ) = a( i, j )
               a( i, j )     = temp
            enddo
            s( j )            = a( j + 1, j + 1 )
            a( j + 1, j + 1 ) = zero
         enddo
         
      end if

      return

c     end of DPUTHS / f06qnf / dpuths. 

      end
