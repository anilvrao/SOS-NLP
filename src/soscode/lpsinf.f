       subroutine   LPSINF   ( n, nclin, ldA,
     1                        istate, bigbnd, numinf, suminf,
     2                        bl, bu, A, featol,
     3                        grad, x, wtinf,
     4                        Ax, w )
c     ==================================================================
c     ==================================================================
c     ====  LPSINF /                                                ====
c     ====  lpsinf -- compute sum of infeasibilities                ====
c     ==================================================================
c     ==================================================================
C     MODIFICATIONS:   13-APR-99 DSN, ADJUSTIBLE ARRAY DIMENSION FIX
C                                     FOR VAX-VMS

      integer            n, nclin, ldA, numinf

      double precision   bigbnd, suminf

      integer            istate (nclin+n)

      double precision   bl (nclin+n), bu (nclin+n), A (ldA,n),
     1                   featol (nclin+n)

      double precision   grad (n), x (n), wtinf (nclin+n), w (*)

      double precision   Ax (ldA)

c     ==================================================================
c
c     derived from qpopt version 1.0 cmsinf
c     last modification -- 26-July-1996
c
c          Original version written 31-October-1984.
c          This version of cmsinf dated  1-January-1987.
c
c     LPSINF / lpsinf  finds the number and weighted sum of 
c     infeasibilities for the bounds and linear constraints.   
c     An appropriate gradient is returned in grad.
c
c     Positive values of  istate(j)  will not be altered.  These mean
c     the following...
c
c               1             2           3
c           a'x = bl      a'x = bu     bl = bu
c
c     Other values of  istate(j)  will be reset as follows...
c           a'x lt bl     a'x gt bu     a'x free
c              - 2           - 1           0
c
c     This version differs from the original source, cmsinf,
c     in that the values of Ax are input to the routine, not
c     recomputed, and in that two options are allowed for the
c     vectorized form of the gradient
c     
c     ==================================================================

      integer            j

      double precision   biglow, bigupp, ctx, feasj, s, weight

      integer            nmglin

      double precision   one, zero, gamma

      parameter        ( one = 1.0d0, zero = 0.0d0, gamma = 0.3d0 )

c     ==================================================================

      bigupp =   bigbnd
      biglow = - bigbnd

      numinf = 0
      nmglin = 0
      suminf = zero

      grad(1:n) = zero
      w(1:nclin) = zero

c     ---------------------------------------------------------------
c     ... decide which inactive constraints are violated.  (on entry,
c         all inactive constraints have istate = 0.)  Compute sum
c         of infeasibilities and also the portion of the gradient
c         due to violated simple bounds
c     ---------------------------------------------------------------

      do j = 1, n+nclin

         if  ( istate(j) .le. 0 )  then

            if  ( j .le. n )  then
               ctx = x (j)
            else
               ctx = Ax (j-n)
            end if

            istate (j) = 0
            weight     = zero
            feasj      = featol(j)

c           ... See if the lower bound is violated.

            if  ( bl(j) .gt. biglow )  then
               s = bl (j) - ctx
               if  ( s .gt. feasj  )  then
                  istate(j) = - 2
                  weight    = - wtinf(j)
                  suminf    = suminf + wtinf(j) * s
               end if
            end if

c           ... If the lower bound was not violated,
c               see if the upper bound is violated.

            if  ( bu(j) .lt. bigupp  .and.  istate (j) .eq. 0 )  then
               s = ctx - bu (j)
               if  ( s .gt. feasj )  then
                  istate(j) = - 1
                  weight    =   wtinf(j)
                  suminf    = suminf + wtinf(j) * s
               endif
            endif

c           ... Add the infeasibility.

            if  ( istate (j) .lt. 0 )  then
               numinf    = numinf + 1

               if  ( j .le. n )  then
                  grad (j) = weight
               else
                  w (j-n)  = weight
                  nmglin   = nmglin + 1
               end if

            end if

         endif

      enddo

c     --------------------------------------------------
c     ... compute gradient of sum of infeasibilities for
c         general linear constraints
c     --------------------------------------------------

      if  ( nmglin .gt. nclin * gamma )  then

c        ... in cases where the fraction of violated general linear
c            constraints is greater than gamma, use a matrix-vector
c            product form for higher speed.

         call dgemv ( 'Trans', nclin, n, one, A, ldA,
     1                         w, 1, one, grad, 1 )

      else

c        ... in more standard cases, use "axpy", but note that
c            stride is not unity

         do j = 1, nclin
            if  ( istate (j+n) .lt. 0 )  then
              do k=1,n
                grad(k) = grad(k) + w(j)*A(j,k)
              enddo
            endif
         enddo

      endif

      return

c     end of LPSINF (lpsinf)

      end
