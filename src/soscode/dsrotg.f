      subroutine   DSROTG  ( pivot0, direc0, n, alpha, x, incx, c, s )
c     =================================================================
c     =================================================================
c     ====  hdqpg /  f06fqf /                                      ====
c     ====  dsrotg  -- rotate vector to unit vector (fancy)        ====
c     =================================================================
c     =================================================================
      
      double precision   alpha

      integer            incx, n
      
      character(len=1)   direct, pivot, direc0, pivot0

      double precision   c (*), s (*), x (*)

c     =================================================================
      
c     derived from qpopt version 1.0 f06fqf
c     last modification -- 26-March-1996
c          originally (NAG f06fqf) mark 12 release. nag copyright 1986.
c     
C  DSROTG generates the parameters of an orthogonal matrix P such that
C
C     when   PIVOT = 'F' or 'f'   and   DIRECT = 'F' or 'f'
C     or     PIVOT = 'V' or 'v'   and   DIRECT = 'B' or 'b'
C
C        P* (alpha) =  (beta),
C           (  x  )    (  0 )
C
C     when   PIVOT = 'F' or 'f'   and   DIRECT = 'B' or 'b'
C     or     PIVOT = 'V' or 'v'   and   DIRECT = 'F' or 'f'
C
C        P* (  x   ) =  (  0 ),
C           (alpha ) =  (beta)
C
C  where alpha is a scalar and x is an n element vector.
C
C  When  PIVOT = 'F' or 'f'  ( fixed pivot )
C  and  DIRECT = 'F' or 'f'  ( forward sequence ) then
C        
C     P is given as the sequence of plane rotation matrices
C
C        P = P (n)*P (n - 1)*...*P (1)
C
C     where P (k) is a plane rotation matrix for the  (1, k + 1) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'V' or 'v'  ( variable pivot )
C  and  DIRECT = 'B' or 'b'  ( backward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P (1)*P (2)*...*P (n)
C
C     where P (k) is a plane rotation matrix for the  (k, k + 1) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'F' or 'f'  ( fixed pivot )
C  and  DIRECT = 'B' or 'b'  ( backward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P (1)*P (2)*...*P (n)
C
C     where P (k) is a plane rotation matrix for the  (k, n + 1) plane
C     designed to annihilate the kth element of x.
C
C  When  PIVOT = 'V' or 'v'  ( variable pivot )
C  and  DIRECT = 'F' or 'f'  ( forward sequence ) then
C
C     P is given as the sequence of plane rotation matrices
C
C        P = P (n)*P (n - 1)*...*P (1)
C
C     where P (k) is a plane rotation matrix for the  (k, k + 1) plane
C     designed to annihilate the kth element of x.
C
C  The routine returns the cosine, c (k), and sine, s (k) that define
C  the matrix P (k), such that the two by two rotation part of P (k),
C  R (k), has the form
C
C     R (k) = (  c (k)  s (k) ).
C             ( -s (k)  c (k) )
C
C  On entry, ALPHA must contain  the scalar alpha and on exit, ALPHA is
C  overwritten by beta. The cosines and sines are returned in the arrays
C  C and S and the vector x is overwritten by the tangents of the plane
C  rotations ( t (k) = s (k) / c (k) ).
C
c  nag fortran 77 o (n) basic linear algebra routine.
c
c  -- written on 19-april-1985.
c     sven hammarling, nag central office.
c
c     =================================================================

      integer            i, ix

      external           DROTGC

c     =================================================================
      pivot  = pivot0(1:1)
      direct = direc0(1:1)

      if  ( n .gt. 0 )  then
         
         if  ( ( direct .eq. 'B' ).or.( direct .eq. 'b' ) )  then
            
            ix = 1 +  (n - 1)*incx
            
            if  ( ( pivot .eq. 'V' ).or.( pivot .eq. 'v' ) )  then
               
               do i = n, 2, -1
c                   << drotgc >>
                  call DROTGC ( x (ix - incx), x (ix), c (i), s (i) )
                  ix = ix - incx
               enddo
               
c                << drotgc >>
               call DROTGC ( alpha, x (ix), c (1), s (1) )
               
            else
     1      if  ( ( pivot .eq. 'F' ).or.( pivot .eq. 'f' ) )  then

c              -------------------------------------------------------
c              here we choose c and s so that
c
c                 ( alpha ) := (  c  s )*( alpha )
c                 (   0   )    ( -s  c ) ( x (i) )
c
c              which is equivalent to
c
c                 (   0   ) := ( c  -s )*( x (i) )
c                 ( alpha )    ( s   c ) ( alpha )
c
c              and so we need to return  s (i) = -s  in order to make
c              r (i) look like
c
c                 r (i) = (  c (i)  s (i) ).
c                         ( -s (i)  c (i) )
c              -------------------------------------------------------

               do i = n, 1, -1
c                   << drotgc >>
                  call DROTGC ( alpha, x (ix), c (i), s (i) )
                  s (i)  = -s (i)
                  x (ix) = -x (ix)
                  ix      =  ix      - incx
               enddo
               
            end if
            
         else
     1   if  ( ( direct .eq. 'F' ).or.( direct .eq. 'f' ) )  then

            ix = 1
            if  ( ( pivot .eq. 'V' ).or.( pivot .eq. 'v' ) )  then

c              -------------------------------------------------------
c              here we choose c and s so that
c
c                 ( x (i + 1) ) := (  c  s )*( x (i + 1) )
c                 (    0      )    ( -s  c ) ( x (i)     )
c
c              which is equivalent to
c
c                 (    0      ) := ( c  -s )*( x (i)     )
c                 ( x (i + 1) )    ( s   c ) ( x (i + 1) )
c
c              and so we need to return  s (i) = -s  in order to make
c              r (i) look like
c
c                 r (i) =  ( c (i)  s (i) ).
c                          (-s (i)  c (i) )
c              -------------------------------------------------------

               do i = 1, n - 1
c                   << drotgc >>
                  call DROTGC ( x (ix + incx), x (ix), c (i), s (i) )
                  s (i)  = -s (i)
                  x (ix) = -x (ix)
                  ix      =  ix      + incx
               enddo
               
c                << drotgc >>
               call DROTGC ( alpha, x (ix), c (n), s (n) )
               s (n)  = -s (n)
               x (ix) = -x (ix)
               
            else
     1      if  ( ( pivot .eq. 'F' ).or.( pivot .eq. 'f' ) )  then

               do i = 1, n
c                   << drotgc >>
                  call DROTGC ( alpha, x (ix), c (i), s (i) )
                  ix = ix + incx
               enddo
               
            end if
            
         end if
         
      end if
c
      return

c     end of DSROTG /  f06fqf / dsrotg

      end
