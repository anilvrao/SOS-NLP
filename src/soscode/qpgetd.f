      subroutine   QPGETD   ( delreg, posdef, statpt, unitgZ, unitQ,
     1                        n, nclin, nfree,
     2                        ldA, ldQ, ldR, nZr,
     3                        issave, jdsave,
     4                        kx, dnorm, gzdz,
     5                        A, Ad, d,
     6                        gq, R, Q, v )

c     ==================================================================
c     ==================================================================
c     ====  QPGETD /                                                ====
c     ====  qpgetd -- compute search direction  d  and  Ad          ====
c     ==================================================================
c     ==================================================================

      integer            n, nclin, nfree,ldA, ldQ, ldR, nZr,
     1                   issave, jdsave
      
      logical            delreg, posdef, statpt, unitgZ, unitQ

      double precision   dnorm, gzdz

      integer            kx (n)

      double precision   A (ldA,*), Ad (*), d (n),
     1                   gq (n), R (ldR,*), Q (ldQ,*),
     2                   v (n)

c     ==================================================================
c
c     derived from qpopt version 1.0
c     last modification -- 17-April-1996
c
c          Original version written 31-December-1986.
c          This version of  qpgetd  dated 21-Dec-1990.
c
c     QPGETD /qpgetd  computes the following quantities for  qpcore (QPCORE).
c     (1) The search direction d (and its 2-norm).  The vector d is
c         defined as  Z*(dz), where  (dz)  is defined as follows.
c         If Hz is positive definite, (dz) is the solution of the
c         (nZr x nZr)  triangular system  (Rz)'(Rz)*(dz) = - (gz).
c         Otherwise  (dz) is the solution of the triangular system
c         (Rz)(dz) =  gamma ez,  where ez is the nZr-th unit vector and
c         gamma = -sgn(gq(nZr)).
c     (2) The vector Ad,  where A is the matrix of linear constraints.
c
c     ==================================================================

      logical            dellow, revers
      
      double precision   atd
      
      double precision   zero, one
      
      parameter        ( zero = 0.0d0, one  = 1.0d0 )

c     ==================================================================

      
      if  ( posdef )  then
         
         if  ( unitgZ )  then
            if  ( nZr .gt. 1 )  then
              v(1:nZr-1) = zero
            endif
            v(nZr) = - gq(nZr)/R(nZr,nZr)
         else
            v(1:nZr) = -gq(1:nZr)
            call dtrsv ( 'Upper', 'Transpose', 'No transpose',
     1                   nZr, R, ldR, v, 1 )
         end if

         
      else
         
         if  ( nZr .gt. 1)  then
           v(1:nZr-1) = zero
         endif
         if  ( gq(nZr) .gt. zero )  then
            v(nZr) = - one
         else
            v(nZr) =   one
         end if
         
      end if

c     Solve  (Rz)*(dz) =  v.

      d(1:nZr) = v(1:nZr)
      call dtrsv ( 'Upper', 'No transpose', 'No transpose',
     1             nZr, R, ldR, d, 1 )

         
c     Compute  d = Zr*(dz)  and its norm.  Find  gz'dz

      dnorm = zero
      gzdz  = zero
      do i=1,nZr
        dnorm = dnorm + d(i)**2
        gzdz = gzdz + d(i)*gq(i)
      enddo
      dnorm = sqrt(dnorm)

c       << cmqmul >>
      call CMQMUL ( 1, n, nZr, nfree, ldQ, unitQ, kx, d, Q, v )

         
c     Compute  Ad.

      if  ( nclin .gt. 0)  then
         call CLKBEG(11)
         call dgemv ( 'No transpose', nclin, n, one, A, ldA,
     1                 d, 1, zero, Ad, 1 )
         call CLKSUM(11)
      endif
      
         
      if  ( delreg  .and.  (gzdz .gt. zero  .or.  statpt) )  then
         
c        ---------------------------------------------------------------
c        The reduced-gradient norm is small enough that we need to worry
c        about the sign of d.  Make  d  point away from the last deleted
c        constraint.
c        ---------------------------------------------------------------
c        Jdsave  is the index of the last deleted regular constraint.

         
         if  ( jdsave .le. n )  then
            atd =  d(jdsave)
         else
            atd = Ad(jdsave-n)
         end if

         dellow = issave .eq. 1
         if  ( dellow )  then
            revers = atd .lt. zero
         else
            revers = atd .gt. zero
         end if

         if  ( revers )  then
            d(1:n) = -d(1:n)
            if  ( nclin .gt. 0)  then
               Ad(1:nclin) = -Ad(1:nclin)
            endif
            gzdz = - gzdz
         end if
         
      end if

c     end of QPGETD (qpgetd)
      
      end
