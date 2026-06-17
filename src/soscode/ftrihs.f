      subroutine   FTRIHS  ( n     , ldH   , nHess, jthcol, H     , 
     1                       x     , delta , Hx    )
c     ==================================================================
c     ==================================================================
c     ====  FTRIHS /                                                ====
c     ====  ftrihs -- dense triangular hessian * vector mult.       ====
c     ==================================================================
c     ==================================================================

c     ... hessian matrix-vector multiplication for standard FORTRAN
c         dense matrix H, stored in triangular (packed) form with no
c         omitted columns or rows (an n x n matrix).  By convention
c         the lower triangle is stored by columns (equivalently the
c         upper triangle is represented by rows)

c     ... last modified 26-March-1996
      
c	================================================================

      integer            n, ldH, nHess, jthcol
      
      double precision   delta
      
      double precision   H (*), Hx (n), x (n)

c     ==================================================================
      
c     FTRIHS / ftrihs  is used to compute the product
c                        (H + delta*I) x, where 
c     H is the QP Hessian matrix stored in H and x is an n-vector.
c
c     parameters:
c        n      -- dimension of Hessian (number of variables)
c        ldH    -- unused in this version
c        nHess  -- number of non-zero rows and columns in H
c                  if  nHess is less than n, only the leading
c                  principal minor of the upper triangle of H
c                  is accessed
c        jthcol -- multiway switch
c                  <= 0  --- compute Hx
c                  >  0  --- compute H times the j-th unit vector
c                           (that is, return the j-th column of H)
c                           (returns zero for j out of range)
c        H      -- Hessian matrix
c        x      -- input vector (not used if jthcol > 0)
c        Hx     -- output vector
c        delta  -- Levenberg parameter
c     
c     This version derived from version of qpHess dated 16-Jan-1995.
      
c     ==================================================================

      integer            i, jstart, jcheck
      
      double precision   one, zero
      
      parameter        ( one = 1.0d0, zero = 0.0d0 )

c     ==================================================================

      i = ldh
      if  ( nHess .lt. n )  then
        Hx(nHess+1:n) = zero
      endif

      if  ( jthcol .gt. 0 )  then
         
c        --------------------------------------------
c        ... Special case -- extract one column of H.
c        --------------------------------------------
         
         if  ( jthcol .gt. nHess )  then
            
           Hx(1:nHess) = zero
            
         else

c           --------------------------------------------------------
c           ... generate upper triangular part of column by marching
c               across a row of lower triangle
c           --------------------------------------------------------

            jstart = jthcol
            do i = 1, jthcol - 1
               Hx (i) = H (jstart)
               jstart = jstart + (nHess - i)
            enddo
            
c           ... lower triangular part of column is accessible
c               directly

            jcheck = 1 + (nHess + 1) * nHess / 2  -
     1               (nHess - jthcol + 2) * (nHess - jthcol + 1) / 2

            if  ( jcheck .ne. jstart )  then
               write (*,*) 'ftrihs didn''t arrive at the right place'
               write (*,*) 'accumulated jstart:', jstart
               write (*,*) 'computed jstart   :', jcheck
            endif
            
            Hx(jthcol:nHess) = H(jstart:jstart+nHess-jthcol)

         end if
         
         Hx (jthcol) = Hx (jthcol) + delta 
            
      else
         
c        ----------------
c        ... Normal case.
c        ----------------
         
         call dspmv ( 'Lower', nHess, one, H, x, 1, zero, Hx, 1 )

         if  ( delta .ne. zero )  then
           do k=1,n
             Hx(k) = Hx(k) + delta*x(k)
           enddo
         endif

      end if

c     end of FTRIHS / ftrihs
      
      end
