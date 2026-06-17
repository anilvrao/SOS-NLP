      subroutine   DGESRC   ( side0, pivot0, direc0, m, n, k1, k2,
     1                        c, s, a, lda )
c     ==================================================================
c     ==================================================================
c     ====  DGESRC / f06qxf /                                       ====
c     ====  dgesrc -- apply collected rotations                     ====
c     ==================================================================
c     ==================================================================
      
      integer            k1, k2, lda, m, n
      
      character(len=1)   direct, pivot, side, direc0, pivot0, side0

      double precision   a( lda, * ), c( * ), s( * )

c     ==================================================================
      
c     derived from qpopt version 1.0 (NAG) f06qxf
c     last modification -- 26-March-1996
c
c     mark 13 release. nag copyright 1988.
c     Nag Fortran 77 O( n**2 ) basic linear algebra routine.

c     DGESRC  performs the transformation
c
c     A := P*A,   when   SIDE = 'L' or 'l'  (  Left-hand side )
c
c     A := A*P',  when   SIDE = 'R' or 'r'  ( Right-hand side )
c
c  where A is an m by n matrix and P is an orthogonal matrix, consisting
c  of a  sequence  of  plane  rotations,  applied  in  planes  k1 to k2,
c  determined by the parameters PIVOT and DIRECT as follows:
c
c     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
c     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
c
c        P is given as a sequence of plane rotation matrices
c
c           P = P( k2 - 1 )*...*P( k1 + 1 )*P( k1 ),
c
c        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
c        plane.
c
c     When  PIVOT  = 'V' or 'v'  ( Variable pivot )
c     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
c
c        P is given as a sequence of plane rotation matrices
c
c           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
c
c        where  P( k )  is a plane rotation matrix for the  ( k, k + 1 )
c        plane.
c
c     When  PIVOT  = 'T' or 't'  ( Top pivot )
c     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
c
c        P is given as a sequence of plane rotation matrices
c
c           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
c
c        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
c        plane.
c
c     When  PIVOT  = 'T' or 't'  ( Top pivot )
c     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
c
c        P is given as a sequence of plane rotation matrices
c
c           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
c
c        where  P( k )  is a plane rotation matrix for the ( k1, k + 1 )
c        plane.
c
c     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
c     and   DIRECT = 'F' or 'f'  ( Forward sequence ) then
c
c        P is given as a sequence of plane rotation matrices
c
c           P = P( k2 - 1 )*P( k2 - 2 )*...*P( k1 ),
c
c        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
c        plane.
c
c     When  PIVOT  = 'B' or 'b'  ( Bottom pivot )
c     and   DIRECT = 'B' or 'b'  ( Backward sequence ) then
c
c        P is given as a sequence of plane rotation matrices
c
c           P = P( k1 )*P( k1 + 1 )*...*P( k2 - 1 ),
c
c        where  P( k )  is a  plane rotation  matrix  for the  ( k, k2 )
c        plane.
c
c  c( k ) and s( k )  must contain the  cosine and sine  that define the
c  matrix  P( k ).  The  two by two  plane rotation  part of the  matrix
c  P( k ), R( k ), is assumed to be of the form
c
c     R( k ) = (  c( k )  s( k ) ).
c              ( -s( k )  c( k ) )
c
c  If m, n or k1 are less than unity,  or k2 is not greater than k1,  or
c  SIDE = 'L' or 'l'  and  k2  is greater than  m, or  SIDE = 'R' or 'r'
c  and  k2  is greater than  n,  then an  immediate return  is effected.

c     ==================================================================
      
      double precision   one, zero
      
      parameter          ( one = 1.0d0, zero = 0.0d0 )

      double precision   aij, ctemp, stemp, temp
      
      integer            i, j
      
      logical            left, right
      
      intrinsic          min

c     ==================================================================
      side   = side0(1:1)
      pivot  = pivot0(1:1)
      direct = direc0(1:1)
           
      left = ( side .eq. 'L' ) .or. ( side .eq. 'l' )
      right = ( side .eq. 'R' ) .or. ( side .eq. 'r' )
      
      if  ( ( min( m, n, k1 ) .lt. 1  ) .or.
     1       ( k2 .le. k1             ) .or.
     2       (  left .and. k2 .gt. m  ) .or.
     3       (  right .and. k2 .gt. n )      )  then
         return
      endif
      
      if  ( left )  then
         
         if  ( ( pivot .eq. 'V' ) .or. ( pivot .eq. 'v' ) )  then
            
            if  ( ( direct .eq. 'F' ) .or. ( direct .eq. 'f' ) )  then
               
               do j = 1, n
                  aij = a( k1, j )
                  do i = k1, k2 - 1
                     temp = a( i + 1, j )
                     a( i, j ) = s( i )*temp + c( i )*aij
                     aij = c( i )*temp - s( i )*aij
                  enddo
                  a( k2, j ) = aij
               enddo
               
            else
     1      if  ( ( direct .eq. 'B' ) .or. ( direct .eq. 'b' ) )  then

               do j = 1, n
                  aij = a( k2, j )
                  do i = k2 - 1, k1, -1
                     temp = a( i, j )
                     a( i + 1, j ) = c( i )*aij - s( i )*temp
                     aij = s( i )*aij + c( i )*temp
                  enddo
                  a( k1, j ) = aij
               enddo
               
            end if
            
         else
     1   if  ( ( pivot .eq. 'T' ) .or. ( pivot .eq. 't' ) )  then
            
            if  ( ( direct .eq. 'F' ) .or. ( direct .eq. 'f' ) )  then

               do j = 1, n
                  temp = a( k1, j )
                  do i = k1, k2 - 1
                     aij = a( i + 1, j )
                     a( i + 1, j ) = c( i )*aij - s( i )*temp
                     temp = s( i )*aij + c( i )*temp
                  enddo
                  a( k1, j ) = temp
               enddo
               
            else
     1      if  ( ( direct .eq. 'B' ) .or. ( direct .eq. 'b' ) )  then

               do j = 1, n
                  temp = a( k1, j )
                  do i = k2 - 1, k1, -1
                     aij = a( i + 1, j )
                     a( i + 1, j ) = c( i )*aij - s( i )*temp
                     temp = s( i )*aij + c( i )*temp
                  enddo
                  a( k1, j ) = temp
               enddo
               
            end if
            
         else
     1   if  ( ( pivot .eq. 'B' ) .or. ( pivot .eq. 'b' ) )  then

            if  ( ( direct .eq. 'F' ) .or. ( direct .eq. 'f' ) )  then

               do j = 1, n
                  temp = a( k2, j )
                  do i = k1, k2 - 1
                     aij = a( i, j )
                     a( i, j ) = s( i )*temp + c( i )*aij
                     temp = c( i )*temp - s( i )*aij
                  enddo
                  a( k2, j ) = temp
               enddo
               
            else
     1      if  ( ( direct .eq. 'B' ) .or. ( direct .eq. 'b' ) )  then

               do j = 1, n
                  temp = a( k2, j )
                  do i = k2 - 1, k1, -1
                     aij = a( i, j )
                     a( i, j ) = s( i )*temp + c( i )*aij
                     temp = c( i )*temp - s( i )*aij
                  enddo
                  a( k2, j ) = temp
               enddo
               
            end if
            
         end if
         
      else
     1if  ( right )  then

         if  ( ( pivot .eq. 'V' ) .or. ( pivot .eq. 'v' ) )  then

            if  ( ( direct .eq. 'F' ) .or. ( direct .eq. 'f' ) )  then
               
               do j = k1, k2 - 1
                  if  ( ( c( j ) .ne. one  ) .or.
     1                  ( s( j ) .ne. zero )      )  then
                     ctemp = c( j )
                     stemp = s( j )
                     do i = 1, m
                        temp = a( i, j + 1 )
                        a( i, j + 1 ) = ctemp*temp - stemp*a( i, j )
                        a( i, j ) = stemp*temp + ctemp*a( i, j )
                     enddo
                  end if
               enddo
               
            else
     1      if  ( ( direct .eq. 'B' ) .or. ( direct .eq. 'b' ) )  then

               do j = k2 - 1, k1, -1
                  if  ( ( c( j ) .ne. one  ) .or.
     1                  ( s( j ) .ne. zero )      )  then
                     ctemp = c( j )
                     stemp = s( j )
                     do i = m, 1, -1
                        temp = a( i, j + 1 )
                        a( i, j + 1 ) = ctemp*temp - stemp*a( i, j )
                        a( i, j ) = stemp*temp + ctemp*a( i, j )
                     enddo
                  end if
               enddo
               
            end if
            
         else
     1   if  ( ( pivot .eq. 'T' ) .or. ( pivot .eq. 't' ) )  then

            if  ( ( direct .eq. 'F' ) .or. ( direct .eq. 'f' ) )  then

               do j = k1 + 1, k2
                  ctemp = c( j - 1 )
                  stemp = s( j - 1 )
                  if  ( ( ctemp .ne. one ) .or. ( stemp .ne. zero ) )
     1            then
                     do i = 1, m
                        temp = a( i, j )
                        a( i, j ) = ctemp*temp - stemp*a( i, k1 )
                        a( i, k1 ) = stemp*temp + ctemp*a( i, k1 )
                     enddo
                  end if
               enddo
               
            else
     1      if  ( ( direct .eq. 'B' ) .or. ( direct .eq. 'b' ) )  then

               do j = k2, k1 + 1, -1
                  ctemp = c( j - 1 )
                  stemp = s( j - 1 )
                  if  ( ( ctemp .ne. one ) .or. ( stemp .ne. zero ) )
     1            then
                     do i = m, 1, -1
                        temp = a( i, j )
                        a( i, j ) = ctemp*temp - stemp*a( i, k1 )
                        a( i, k1 ) = stemp*temp + ctemp*a( i, k1 )
                     enddo
                  end if
               enddo
               
            end if
            
         else
     1   if  ( ( pivot .eq. 'B' ) .or. ( pivot .eq. 'b' ) )  then

            if  ( ( direct .eq. 'F' ) .or. ( direct .eq. 'f' ) )  then

               do j = k1, k2 - 1
                  if  ( ( c( j ) .ne. one ) .or. ( s( j ) .ne. zero ) )
     1            then
                     ctemp = c( j )
                     stemp = s( j )
                     do i = 1, m
                        temp = a( i, j )
                        a( i, j ) = stemp*a( i, k2 ) + ctemp*temp
                        a( i, k2 ) = ctemp*a( i, k2 ) - stemp*temp
                     enddo
                  end if
               enddo
               
            else
     1      if  ( ( direct .eq. 'B' ) .or. ( direct .eq. 'b' ) )  then

               do j = k2 - 1, k1, -1
                  if  ( ( c( j ) .ne. one ) .or. ( s( j ) .ne. zero ) )
     1            then
                     ctemp = c( j )
                     stemp = s( j )
                     do i = m, 1, -1
                        temp = a( i, j )
                        a( i, j ) = stemp*a( i, k2 ) + ctemp*temp
                        a( i, k2 ) = ctemp*a( i, k2 ) - stemp*temp
                     enddo
                  end if
               enddo
               
            end if
            
         end if
         
      end if
      
c
      return

c     end of DGESRC  / f06qxf / dgesrc

      end
