      subroutine  DRHSTU   ( side0, n, k1, k2, c, s, a, lda )
c     ==================================================================
c     ==================================================================
c     ====  DRHSTU / f06qrf /                                       ====
c     ====  drhstu -- rotate hessenberg to upper triangular         ====
c     ==================================================================
c     ==================================================================
      
      integer            k1, k2, lda, n
      
      character(len=1)   side, side0

      double precision   a( lda, * ), c( * ), s( * )

c     ==================================================================

c     derived from qpopt version 1.0  f06qrf
c     last modification -- 26-March-1996
c
c          originally  NAG  f06qrf mark 13 release. nag copyright 1988.

      
c     DRHSTU /drhstu restores an upper Hessenberg matrix H to upper 
c     triangular form by  applying a sequence of  plane rotations  
c     from either the left, or the right.  The matrix  H  is assumed to
c     have  non-zero  sub-diagonal elements  in  positions
c         h( k + 1, k ),  k = k1, k1 + 1, ..., k2 - 1, only  and
c         h( k + 1, k )  must  be  supplied  in  s( k ).
c
c  H is restored to the upper triangular matrix R either as
c
c     R = P*H,   when   SIDE = 'L' or 'l'  (  Left-hand side )
c
c  where P is an orthogonal matrix of the form
c
c     P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
c
c  or as
c
c     R = H*P',  when   SIDE = 'R' or 'r'  ( Right-hand side )
c
c  where P is an orthogonal matrix of the form
c
c     P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
c
c  in both cases  P( k )  being a  plane rotation  for the  ( k, k + 1 )
c  plane.  The cosine and sine that define P( k ) are returned in c( k )
c  and  s( k )  respectively.  The two by two  rotation part of  P( k ),
c  Q( k ), is of the form
c
c     Q( k ) = (  c( k )  s( k ) ).
c              ( -s( k )  c( k ) )
c
c  The upper triangular part of the matrix  H  must be supplied in the n
c  by n  leading upper triangular part of  A, and this is overwritten by
c  the upper triangular matrix R.
c
c  If n or k1 are less than unity,  or k1 is not less than k2,  or k2 is
c  greater than n then an immediate return is effected.
c
c
c  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
c
c  -- Written on 13-January-1986.
c     Sven Hammarling, Nag Central Office.
c
c     ==================================================================

      double precision   one, zero
      parameter        ( one = 1.0d0, zero = 0.0d0 )

      double precision   aij, ctemp, stemp, subh, temp
      
      integer            i, j

      external           DROTGC

      intrinsic          min

c     ==================================================================
      side = side0(1:1)

      if  ( ( min( n, k1 ).lt.1 ) .or.
     1      ( k2.le.k1          ) .or.
     2      ( k2 .gt. n         )      )  then
         return
      endif
      
      if  ( ( side.eq.'l' ) .or. ( side.eq.'l' ) )  then

c        ---------------------------------------------------------------
c        restore   h  to  upper  triangular  form  by  annihilating  the
c        sub-diagonal elements of h.  the jth rotation is chosen so that
c
c           ( h( j, j ) ) := (  c  s )*( h( j, j )     ).
c           (     0     )    ( -s  c ) ( h( j + 1, j ) )
c
c        apply the rotations in columns k1 up to n.
c        ---------------------------------------------------------------

         do j = k1, n
            
            aij = a( k1, j )
            do i = k1, min( j, k2 ) - 1
               temp      = a( i + 1, j )
               a( i, j ) = s( i )*temp + c( i )*aij
               aij       = c( i )*temp - s( i )*aij
            enddo
            
            if  ( j .lt. k2 )  then

c              set up the rotation.

               subh = s( j )
c                << drotgc >>
               call DROTGC ( aij, subh, c( j ), s( j ) )
               a( j, j ) = aij
               
            else
               
               a( k2, j ) = aij
               
            end if
            
         enddo
         
      else
     1if  ( ( side.eq.'r' ) .or. ( side.eq.'r' ) )  then

c        ---------------------------------------------------------------
c        restore   h  to  upper  triangular  form  by  annihilating  the
c        sub-diagonal elements of h.  the jth rotation is chosen so that
c
c           ( h( j + 1, j + 1 ) ) := (  c  s )*( h( j + 1, j + 1 ) ),
c           (         0         )    ( -s  c ) ( h( j + 1, j )     )
c
c        which can be expressed as
c
c           ( 0  h( j + 1, j + 1 ) ) :=
c
c               ( h( j + 1, j )  h( j + 1, j + 1 ) )*(  c  s ).
c                                                    ( -s  c )
c
c        thus we return  c( j ) = c  and  s( j ) = -s  to make the plane
c        rotation matrix look like
c
c           q( j ) = (  c( j )  s( j ) ).
c                    ( -s( j )  c( j ) )
c        ---------------------------------------------------------------

         do j = k2 - 1, k1, -1
            subh = s( j )
c             << drotgc >>
            call DROTGC ( a( j + 1, j + 1 ), subh, ctemp, stemp )
            stemp = -stemp
            s( j ) = stemp
            c( j ) = ctemp
            
            if  ( ( ctemp.ne.one ) .or. ( stemp.ne.zero ) )  then
               do i = j, 1, -1
                  temp = a( i, j + 1 )
                  a( i, j + 1 ) = ctemp*temp - stemp*a( i, j )
                  a( i, j ) = stemp*temp + ctemp*a( i, j )
               enddo
            end if
            
         enddo
         
      end if

      return

c     end of drhstu / f06qrf / DRHSTU

      end
