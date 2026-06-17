      subroutine   IDENHS  ( n     , ldH   , nHess , jthcol, H     , 
     1                       x     , delta , Hx    )
c     ==================================================================
c     ==================================================================
c     ====  IDENHS /                                                ====
c     ====  idenhs -- identity hessian * vector multiplication      ====
c     ==================================================================
c     ==================================================================

c     ... hessian matrix-vector multiplication for special case of
c         least distance programming, where the hessian matrix is
c         the identity

      integer            n, ldH, jthcol, nHess

      double precision   delta
      
      double precision   H (*), Hx (n), x (n)

c     ==================================================================
c     
c     IDENHS / idenhs  is used to compute the product (H + delta*I) x.
c
c     parameters:
c        n      -- dimension of Hessian (number of variables)
c        ldH    -- ignored
c        jthcol -- multiway switch
c                  <= 0  --- compute Hx
c                  >  0  --- compute I times the j-th unit vector
c                           (that is, return the j-th column of I)
c                           (returns zero for j out of range)
c        nHess  -- number of non-zero rows and columns in H
c                  if  nHess is less than n, only the leading
c                  nHess rows and columns of H are taken to be
c                  the identity, with the rest taken to be zero.
c        H (=I) -- Hessian matrix (ignored)
c        x      -- input vector (not used if jthcol > 0)
c        Hx     -- output vector
c        delta  -- Levenberg parameter
c     
c     This version derived from version of qpHess dated 16-Jan-1995.
      
c     ... last modified 26-March-1996

c     ==================================================================
      
      integer            j
      
      double precision   one, a, zero
      
      parameter        ( one = 1.0d0, zero = 0.0d0 )

      j = ldh
      a = h(1)
      if  ( jthcol .gt. 0 )  then
         
c        ----------------------------------------
c        Special case -- extract one column of H.
c        ----------------------------------------
         
         Hx(1:n) = zero

         if  ( jthcol .le. nHess )  then

            hx (jthcol) = one + delta

         else

            hx (jthcol) = delta

         endif
         
      else
         
c        ------------
c        Normal case.
c        ------------
         
         do j = 1, nHess
            hx (j) = (one + delta) * x (j)
         enddo

         do j = nHess+1, n
            hx (j) = delta * x (j)
         enddo

      end if

c     end of IDENHS / idenhs
      
      end
