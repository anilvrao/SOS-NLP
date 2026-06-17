      subroutine   TRED1    ( nm, n, a, d, e, e2 )
      
c     ==================================================================
c     ==================================================================
c     ====  TRED1  /                                                ====
c     ====  tred1 -- eispack tridiagonalize symmetric matrix        ====
c     ==================================================================
c     ==================================================================
      
c     last modification 25-March-1996
      
c     this subroutine is a translation of the algol procedure tred1,
c     num. math. 11, 181-195 (1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226 (1971).
c
c     this subroutine reduces a real symmetric matrix
c     to a symmetric tridiagonal matrix using
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction in its strict lower
c          triangle.  the full upper triangle of a is unaltered.
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e (1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer           n, nm
      
      double precision  a (nm,n), d (n), e (n), e2 (n)
      
      integer           i, ii, j, jp1, k, l
      
      double precision  f, g, h, scale

      double precision  zero

      parameter       ( zero = 0.0d0 )
      
c     ------------------------------------------------------------------
c
      
       do i = 1, n
         d (i)   = a (n,i)
         a (n,i) = a (i,i)
      enddo
      
c     .......... for i=n step -1 until 1 do -- ..........
      
      do 300 ii = 1, n
         
         i     = n + 1 - ii
         l     = i - 1
         h     = zero
         scale = zero
         
         if  ( l .lt.  1) go to 130
         
c        .......... scale row (algol tol then not needed) ..........
         
         do k = 1, l
         scale = scale + abs (d (k))   
         enddo
c
         if  ( scale .ne. zero )  go to 140
c
         do j = 1, l
            d (j) = a (l,j)
            a (l,j) = a (i,j)
            a (i,j) = zero
         enddo
c
  130    continue
         e (i)  = zero
         e2 (i) = zero
         go to 300
c
  140    continue
         do k = 1, l
            d (k) = d (k) / scale
            h     = h + d (k) * d (k)
         enddo
c
         e2 (i) = scale * scale * h
         f      = d (l)
         g      = -sign (sqrt (h),f)
         e (i)  = scale * g
         h      = h - f * g
         d (l)  = f - g
         
         if  ( l .gt. 1 )  then
         
c           .......... form a*u ..........
         
            do j = 1, l
               e (j) = zero
            enddo

            do j = 1, l
               
               f   = d (j)
               g   = e (j) + a (j,j) * f
               jp1 = j + 1
               
               if  ( l .ge. jp1 ) then
                  
                  do k = jp1, l
                     g     = g + a (k,j) * d (k)
                     e (k) = e (k) + a (k,j) * f
                  enddo
                  
               endif
              
               e (j) = g
               
            enddo
         
c           .......... form p ..........
            
            f = zero
           
            do j = 1, l
               e (j) = e (j) / h
               f     = f + e (j) * d (j)
            enddo
           
            h = f /  (h + h)
            
c           .......... form q ..........
            
            do j = 1, l
               e (j) = e (j) - h * d (j)
            enddo
            
c           .......... form reduced a ..........
            
            do j = 1, l
               
               f = d (j)
               g = e (j)
              
               do k = j, l
                  a (k,j) = a (k,j) - f * e (k) - g * d (k)
               enddo
c              
            enddo

         endif

         do j = 1, l
            f       = d (j)
            d (j)   = a (l,j)
            a (l,j) = a (i,j)
            a (i,j) = f * scale
         enddo

  300 continue

      
      return
      
c     end of TRED1  (tred1)
      end
