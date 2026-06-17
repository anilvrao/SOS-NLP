      subroutine   CMSETX   ( rowerr, unitQ,
     1                        nclin, nactiv, nfree, nz,
     2                        n, ldQ, ldA, ldT,
     3                        istate, kactiv, kx,
     4                        errmax, xnorm,
     5                        A, Ax, bl, bu, featol,
     6                        T, x, Q, p, work )
c     ==================================================================
c     ==================================================================
c     ====  CMSETX /                                                ====
c     ====  cmsetx -- find nearest point on working set             ====
c     ==================================================================
c     ==================================================================

      integer            nclin, nactiv, nfree, nz,
     1                   n, ldQ, ldA, ldT
      
      logical            rowerr, unitQ

      double precision   errmax, xnorm
      
      integer            istate (n+nclin), kactiv (n), kx (n)
      
      double precision   A (ldA,*), Ax (*), bl (n+nclin), bu (n+nclin)
      
      double precision   featol (n+nclin), p (n),
     1                   T (ldT,*), x (n), Q (ldQ,*)

      double precision   work (n)

c     ==================================================================
c     derived from qpopt version 1.0 cminit
c     last modification -- 05-June-1996
c         Original version derived from lssetx January-1987.
c         This  version of  cmsetx  dated   5-Jul-1989.
c
c     CMSETX / cmsetx  computes the point on a working set that is 
c     closest in the least-squares sense to the input vector X.
c
c     If the computed point gives a row error of more than the
c     feasibility tolerance, an extra step of iterative refinement is
c     used.  If  X  is still infeasible,  the logical variable ROWERR
c     is set.
c
c     ==================================================================

      integer            i, is, j, jmax, k, ktry
      
      double precision   bnd
      
      integer            ntry

      double precision   zero, one
      
      parameter        ( ntry  = 5 )

      parameter        ( zero  = 0.0d0, one = 1.0d0 )

c     ==================================================================

c     ------------------------------------------------------------------
c     Move  x  onto the simple bounds in the working set.
c     ------------------------------------------------------------------

      do k = nfree+1, n
          j   = kx(k)
          is  = istate(j)
          bnd = bl(j)
          if  ( is .ge. 2 )  bnd  = bu(j)
          if  ( is .ne. 4 )  x(j) = bnd
      enddo

c     ------------------------------------------------------------------
c     Move  x  onto the general constraints in the working set.
c     ntry  attempts are made to get acceptable row errors.
c     ------------------------------------------------------------------

      ktry   = 1
      jmax   = 1
      errmax = zero

c     repeat
  200    continue
         if  ( nactiv .gt. 0 )  then

c           Set work = residuals for constraints in the working set.
c           Solve for P, the smallest correction to x that gives a point
c           on the constraints in the working set.  Define  P = Y*(py),
c           where  py  solves the triangular system  T*(py) = residuals.

            do i = 1, nactiv
               
               k   = kactiv(i)
               j   = n + k
               bnd = bl(j)
               if  ( istate(j) .eq. 2)  then
                  bnd = bu(j)
               endif
               
               work(nactiv-i+1) = bnd - dot_product(A(k,1:n),x(1:n))
               
            enddo

            call dtrsv ( 'Upper', 'No transpose', 'No transpose',
     1                   nactiv, T(1,nz+1), ldT, work, 1 )
            p(1:n) = zero
            p(nz+1:nz+nactiv) = work(1:nactiv)

c             << cmqmul >>
            call CMQMUL ( 2, n, nz, nfree, ldQ, unitQ, kx, p, Q, work )
            do i=1,n
              x(i) = x(i) + p(i)
            enddo

         end if

c        ---------------------------------------------------------------
c        Compute the 2-norm of  x.
c        Initialize  Ax  for all the general constraints.
c        ---------------------------------------------------------------

         xnorm = zero
         do i=1,n
           xnorm = xnorm + x(i)**2
         enddo
         xnorm  = sqrt(xnorm)
         if  ( nclin .gt. 0 )  then
            call dgemv ( 'No transpose', nclin, n, one, A, ldA,
     1                    x, 1, zero, Ax, 1 )
         endif
         
c        ---------------------------------------------------------------
c        Check the row residuals.
c        ---------------------------------------------------------------

         if  ( nactiv .gt. 0 )  then
            
            jmax = 1
            errmax = zero
            do k = 1, nactiv
               i   = kactiv(k)
               j   = n + i
               if  ( istate (j) .eq. 1 )  then
                  work (k) = bl (j) - Ax (i)
               else
     1         if  ( istate (j) .ge. 2 )  then
                  work (k) = bu (j) - Ax (i)
               endif
               if (abs(work(k)).gt.errmax) then
                 jmax = k
                 errmax = abs(work(k))
               endif
            enddo
            
         end if

         ktry = ktry + 1
         
c        until    (errmax .le. featol (jmax) .or. ktry .gt. ntry

c        NOTE -- PHILIP DOES NOT CHECK EACH INFEASIBILITY AGAINST
C        ITS INDIVIDUAL TOLERANCE
         
         if  ( .not. (errmax .le. featol (jmax) .or.
     1               ktry .gt. ntry) )  then
            go to 200
         endif

      rowerr = errmax .gt. featol(jmax)

c     end of  CMSETX (cmsetx)
      
      end
