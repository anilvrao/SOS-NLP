      subroutine   QPCOLR   ( singlr, posdef, renewr, unitQ,
     1                        n, nZr, nfree, nHess, ldQ, ldH, ldR,
     2                        kx, Hsize, Dzz, tolrnk,
     3                        qpHess, delta, H, R, Q,
     4                        Hz, wrk )

      integer            n, nZr, nfree, nHess, ldQ, ldH, ldR
      
      logical            singlr, posdef, renewr, unitQ

      double precision   Hsize, Dzz, tolrnk, delta
      
      integer            kx (n)

      double precision   H (ldH,*), R (ldR,*), Hz (n), Q (ldQ,*),
     1                   wrk (n)
      
      external           qpHess
      
c     ==================================================================
c     ==================================================================
c     ====  QPCOLR /                                                ====
c     ====  qpcolr -- compute last column of reduced Hessian factor ====
c     ==================================================================
c     ==================================================================
c
c     derived from qpopt version 1.0
c     last modification -- 17-April-1996
c
c          Original f66 version written by PEG,   March-1982.
c          This version of  qpcolr  dated 16-Jan-1995.
c     
c     QPCOLR / QPCOLR  is used to compute the last column of the
c     (nZr x nZr)  triangular factor Rz such that
c                    Hz  =  (Rz)'D(Rz),
c     where Hz is the reduced Hessian Z'HZ, and D is a diagonal
c     matrix.  If  Hz  is positive definite, Rz is the Cholesky factor
c     of  Hz  and  D  is the identity matrix;  otherwise, D(nZr) is
c     negative or small and the last diagonal of Rz is one.
c
c     The element D(nZr) is stored in Dzz.  Dzz is equal to one if
c     posdef is true.
c
c     ==================================================================

      integer            j, jthcol, k, nZr1

      double precision   dRzmax, dRzmin, Dzznew, rznorm, rdsmin,
     1                   Rzz, ztHz
      
      double precision   zero, one, ten
      
      parameter        ( zero = 0.0d0, one = 1.0d0, ten = 10.0d0 )

c     ==================================================================

      if  ( nZr .eq. 0 )  then
         
         posdef = .true.
         renewr = .false.
         singlr = .false.
         Dzz    =  one
         
      else

C        ... TEST OF RENEWR WILL COME OUT--jgl
         
         if  ( renewr)  then
            
c           ------------------------------------------------------------
c           Compute the first nZr-1 elements of the last column of Rz
c           and Dzznew, the square of the last diagonal element.
c           ------------------------------------------------------------

            wrk(1:n) = zero

c           ... probably could avoid rezeroing  wrk  if we tried ---jgl
            
            if  ( unitQ) then

c              Only bounds are in the working set.  The nZr-th column of
c              Z is just a column of the identity matrix.

               jthcol      = kx(nZr)
               wrk(jthcol) = one
               
            else

c              Expand the new column of  Z  into an n-vector.

               do k = 1, nfree
                  j      = kx(k)
                  wrk(j) = Q(k,nZr)
               enddo
               jthcol = 0
            end if

c           Compute the nZr-th column of Z'HZ.

            call qpHess ( n, ldH, nHess, jthcol, H, wrk, delta, Hz )
c             << cmqmul >>
            call CMQMUL ( 4, n, nZr, nfree, ldQ, unitQ,
     1                    kx, Hz, Q, wrk )
            R(1:nZr,nZr) = Hz(1:nZr)

            nZr1   = nZr - 1
            zthz   = R(nZr,nZr)
            Dzznew = zthz
            rznorm = zero
            
            if  ( nZr1 .gt. 0)  then
               call dtrsv ( 'Upper', 'Transpose', 'No transpose',
     1                      nZr1, R, ldR, R(1,nZr), 1 )
               rznorm = zero
               do i=1,nZr1
                 rznorm = rznorm + R(i,nZr)**2
               enddo
               rznorm  = sqrt(rznorm)
               Dzznew  = zthz - rznorm*rznorm
            end if

            R(nZr,nZr) = one
            Dzz        = Dzznew

c           Update the estimate of the norm of the Hessian.

corig       Hsize  = max( Hsize, abs( zthz ) )
            Hsize  = max( Hsize, sqrt (rznorm**2 + zthz**2) )

         end if

         Dzznew = Dzz*R(nZr,nZr)**2

c        ---------------------------------------------------------------
c        Attempt to compute Rzz, the square root of  Dzznew.  The last
c        diagonal of Rz.  The variables posdef and singlr are set here.
c        They are used to indicate if the new Z'HZ is positive definite
c        or singular.  If the required diagonal modification is large
c        the last row and column of Rz are marked for recomputation next
c        iteration.
c        ---------------------------------------------------------------
c        Rdsmin is the square of the smallest allowable diagonal element
c        for a positive-definite Cholesky factor.  Note that the test
c        for positive definiteness is unavoidably scale dependent.

         if  ( nZr .eq. 1)  then
            rdsmin =  tolrnk*Hsize
         else
c             << dcond >>
            call DCOND  ( nZr, R, ldR+1, dRzmax, dRzmin )
            rdsmin = (tolrnk*dRzmax)*dRzmax
         end if

         posdef =  Dzznew .gt. rdsmin

         if  ( posdef) then
            Dzz    = one
            Rzz    = sqrt( Dzznew )
            renewr = .false.
            singlr = .false.
         else
            Dzz    = Dzznew
            Rzz    = one
            singlr = Dzznew .ge. - rdsmin
            renewr = Dzznew .lt. - ten*Hsize
         end if

         R(nZr,nZr) = Rzz
      end if

c     end of QPCOLR (qpcolr)
      
      end
