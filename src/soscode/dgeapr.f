      subroutine   DGEAPR   ( side0, trans0, n, perm, k, b, ldb )
c     ==================================================================
c     ==================================================================
c     ====  DGEAPR / f06qkf /                                       ====
c     ====  dgeapr -- permute matrix                                ====
c     ==================================================================
c     ==================================================================
      
c     mark 13 release. nag copyright 1988.

      integer            k, ldb, n
      
      character(len=1)   side, trans, side0, trans0

      double precision   perm ( * ), b ( ldb, * )

c     ==================================================================
c
c     derived from qpopt version 1.0 (NAG) f06qkf
c     last modification -- 25-March-1996
c     
c  Purpose
c  =======
c
c  DGEAPR / f06qkf / dgeapr performs one of the transformations
c
c     B := P'*B   or   B := P*B,   where B is an m by k matrix,
c
c  or
c
c     B := B*P'   or   B := B*P,   where B is a k by m matrix,
c
c  P being an m by m permutation matrix of the form
c
c     P = P( 1, index( 1 ) )*P( 2, index( 2 ) )*...*P( n, index( n ) ),
c
c  where  P( i, index( i ) ) is the permutation matrix that interchanges
c  items i and index( i ). That is P( i, index( i ) ) is the unit matrix
c  with rows and columns  i and  index( i )  interchanged. Of course, if
c  index( i ) = i  then  P( i, index( i ) ) = I.
c
c  This  routine is intended  for use in conjunction with  Nag auxiliary
c  routines  that  perform  interchange  operations,  such  as  sorting.
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER(LEN=1).
c  TRANS
c           On entry,  SIDE  ( Left-hand side, or Right-hand side )  and
c           TRANS  ( Transpose, or No transpose )  specify the operation
c           to be performed as follows.
c
c           SIDE = 'L' or 'l'   and   TRANS = 'T' or 't'
c
c              Perform the operation   B := P'*B.
c
c           SIDE = 'L' or 'l'   and   TRANS = 'N' or 'n'
c
c              Perform the operation   B := P*B.
c
c           SIDE = 'R' or 'r'   and   TRANS = 'T' or 't'
c
c              Perform the operation   B := B*P'.
c
c           SIDE = 'R' or 'r'   and   TRANS = 'N' or 'n'
c
c              Perform the operation   B := B*P.
c
c           Unchanged on exit.
c
c  N      - INTEGER.
c
c           On entry, N must specify the value of n.  N must be at least
c           zero.  When  N = 0  then an  immediate  return  is effected.
c
c           Unchanged on exit.
c
c  PERM   - REAL             array of DIMENSION at least ( n ).
c
c           Before  entry,  PERM  must  contain  the  n indices  for the
c           permutation matrices. index( i ) must satisfy
c
c              1 .le. index( i ) .le. m.
c
c           It is usual for index( i ) to be at least i, but this is not
c           necessary for this routine. It is assumed that the statement
c           INDEX = PERM( I )  returns the correct integer in  INDEX, so
c           that,  if necessary,  PERM( I )  should contain a real value
c           slightly larger than  INDEX.
c
c           Unchanged on exit.
c
c  K      - INTEGER.
c
c           On entry with  SIDE = 'L' or 'l',  K must specify the number
c           of columns of B and on entry with  SIDE = 'R' or 'r', K must
c           specify the number of rows of  B.  K must be at least  zero.
c           When  K = 0  then an immediate return is effected.
c
c           Unchanged on exit.
c
c  B      - REAL  array  of  DIMENSION ( LDB, ncolb ),  where  ncolb = k
c           when  SIDE = 'L' or 'l'  and  ncolb = m  when  SIDE = 'R' or
c           'r'.
c
c           Before entry  with  SIDE = 'L' or 'l',  the  leading  m by K
c           part  of  the  array   B  must  contain  the  matrix  to  be
c           transformed  and before  entry with  SIDE = 'R' or 'r',  the
c           leading  K by m part of the array  B must contain the matrix
c           to  be  transformed.  On exit,   B  is  overwritten  by  the
c           transformed matrix.
c
c  LDB    - INTEGER.
c
c           On entry,  LDB  must specify  the  leading dimension  of the
c           array  B  as declared  in the  calling  (sub) program.  When
c           SIDE = 'L' or 'l'   then  LDB  must  be  at  least  m,  when
c           SIDE = 'R' or 'r'   then  LDB  must  be  at  least  k.
c           Unchanged on exit.
c
c
c  Nag Fortran 77 O( n**2 ) basic linear algebra routine.
c
c  -- Written on 11-August-1987.
c     Sven Hammarling, Nag Central Office.
c
c     ==================================================================
      
      logical            left, null, right, trnsp
      
      integer            i, j, l
      
      double precision   temp
      
      intrinsic          min

c     ==================================================================
      side  = side0(1:1)
      trans = trans0(1:1)
      
      if  ( min( n, k ) .eq. 0 )  then
         return
      endif
      
      left  = (  side .eq. 'l' ) .or. (  side .eq. 'L' )
      right = (  side .eq. 'r' ) .or. (  side .eq. 'R' )
      null  = ( trans .eq. 'n' ) .or. ( trans .eq. 'N' )
      trnsp = ( trans .eq. 't' ) .or. ( trans .eq. 'T' )
      
      if  ( left )  then
         
         if  ( trnsp )  then
            
            do i = 1, n
               l = perm( i )
               if  ( l .ne. i )  then
                  do j = 1, k
                     temp = b( i, j )
                     b( i, j ) = b( l, j )
                     b( l, j ) = temp
                  enddo
               end if
            enddo
            
         else
     1   if  ( null )  then
            
            do i = n, 1, -1
               l = perm( i )
               if  ( l .ne. i )  then
                  do j = 1, k
                     temp = b( l, j )
                     b( l, j ) = b( i, j )
                     b( i, j ) = temp
                  enddo
               end if
            enddo
         end if
         
      else
     1if  ( right )  then

         if  ( trnsp )  then
            
            do j = n, 1, -1
               l = perm( j )
               if  ( l .ne. j )  then
                  do i = 1, k
                     temp = b( i, j )
                     b( i, j ) = b( i, l )
                     b( i, l ) = temp
                  enddo
               end if
            enddo
            
         else
     1   if  ( null )  then
            
            do j = 1, n
               l = perm( j )
               if  ( l .ne. j )  then
                  do i = 1, k
                     temp = b( i, l )
                     b( i, l ) = b( i, j )
                     b( i, j ) = temp
                  enddo
               end if
            enddo
            
         end if
         
      end if
c
      return
c
c     end of DGEAPR / f06qkf / dgeapr.
c
      end
