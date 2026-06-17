      subroutine   QPMDLV   ( delta , nZr   , ldZtHZ, epslon, levtol, 
     1                        leftvl, rghtvl, ZtHZ  , tdiag , toffd ,
     2                        toffd2, msglvl, output )

c     ==================================================================
c     ==================================================================
c     ====  qpmdlv / QPMDLV -- revise Levenberg parameter           ====
c     ==================================================================
c     ==================================================================

      integer            nZr, ldZtHZ, msglvl, output
      
      double precision   delta, epslon, levtol, leftvl, rghtvl
      
      double precision   ZtHZ (ldZtHZ, *)

      double precision   tdiag (nZr), toffd (nZr), toffd2 (nZr)
      
c     ==================================================================
c     
c     last modification -- 17-April-1996
c     
c     qpmdlv (QPMDLV)  revises the Levenberg parameter so that the current
c     reduced Hessian is adequately positive definite.  It is called
c     whenever the current Levenberg parameter leaves the reduced
c     Hessian either indefinite or numerically singular.
c
c     qpmdlv (QPMDLV)  computes the two extreme eigenvalues of the 
c     reduced Hessian  Z'HZ,  given the (nfree x nZr) matrix Z.
c     Delta is then adjusted by  mu  so that  Z'HZ + mu*I  is
c     adequately positive definite.
c
c     ==================================================================

      integer            i, j 

      double precision   relprc, mu

      double precision   zero, one

      parameter        ( zero = 0.0d0, one = 1.0d0 )

c     ==================================================================

c     ... the problem is to find the leftmost and the rightmost
c         eigenvalues of a dense symmetric matrix.  As of November, 1995
c         we foresee three different approaches to this problem:
c         1.  use LAPACK to reduce the dense matrix to tridiagonal form
c             and then use bisection for the eigenvalues
c         2.  use EISPACK for the same task
c         3.  develop a Lanczos code, probably based on the 1989 version
c             of Beresford Parlett's code.  This problem does NOT require
c             any sort of reorthogonalization, but the "analyze T"
c             procedure would probably give an efficiency boost.
c
c         this initial version does not use Lanczos because it's not
c         clear that there are efficiency reasons for doing the extra
c         programming -- the reduced Hessians should be small.
c     
c         the comments below contain UNTESTED code for calling LAPACK
c         to find the two extreme eigenvalues.
c         LAPACK requires more storage than EISPACK.  LAPACK may
c         have a speed advantage if the order of the reduced hessians
c         is large -- but then Lanczos may be still better.
c        
clapackc     ------------------------------
clapackc     ... reduce to tridiagonal form
clapackc     ------------------------------
clapack
clapack      call dsytrd ( 'Upper', nZ, ZtHZ, ldZtHZ, tdiag, toffd, 
clapack     1              tau, work, lwork, info )
clapack
clapackc     --------------------------------------------------
clapackc     ... find the two extreme eigenvalues by bisection.
clapackc         (tridiagonal should be scaled for lapack.)
clapackc        --------------------------------------------------
clapack
clapack      call dstebz ( 'Index', nZ, rdummy, rdummy, 1, 1, zero,
clapack     1              tdiag, toffd, ifound, nsplit, tau, iblock,
clapack     2              isplit, work, iwork, info )
clapack
clapack      leftvl = tau (1)
clapack
clapack      call dstebz ( 'Index', nZ, rdummy, rdummy, nZ, nZ, zero,
clapack     1              tdiag, toffd, ifound, nsplit, tau, iblock,
clapack     2              isplit, work, iwork, info )
clapack
clapack      rghtvl = tau (nZ)

c     ------------------------------
c     ... reduce to tridiagonal form
c     ------------------------------

c       << tred1 >>
      call TRED1  ( ldZtHZ, nZr, ZtHZ, tdiag, toffd, toffd2 )
      
c     --------------------------------------------------
c     ... find the two extreme eigenvalues by bisection.
c     --------------------------------------------------
      
      if  ( nZr .gt. 1 )  then
         
         relprc = sqrt ( epslon ) 
c          << tridb2 >>
         call TRIDB2 ( nZr   , relprc, tdiag, toffd, toffd2,
     1                 leftvl, rghtvl )

c        ------------------------------------------------------------
c        ... adjust delta -- choose a  mu  so that  Z' H Z + delta*I
c            has condition number  levtol, where levtol  should be
c            chosen in the range  (1, 1/tolrnk).
c        ------------------------------------------------------------

         mu    = max ( ( rghtvl - levtol*leftvl ) / ( levtol - one ),
     1                 one / levtol - leftvl,
     2                 zero )

         if  ( msglvl .ge. 10 .and. output .gt. 0 )  then
            write (output, 60000) delta, mu, delta + mu,
     1                            rghtvl, leftvl, levtol
         endif
         
      else

         leftvl = tdiag (1)
         rghtvl = leftvl
         
         if  ( leftvl .le. zero )  then

c           ---------------------------------------------------------
c           ... special case -- negative definite one by one matrix.
c               adjust delta so  Z' H Z + delta*I  has eigenvalue one
c           ---------------------------------------------------------
            
            mu    = one - leftvl
            
            if  ( msglvl .ge. 10 .and. output .gt. 0 )  then
               write (output, 61000) delta, mu, delta + mu, leftvl
            endif

         else

            mu = zero
           
         endif
         
      endif
      
      delta = delta + mu

c     -------------------------------------------------------
c     ... adjust the diagonals of the current reduced Hessian
c         to reflect the new value of the Levenberg parameter
c         and restore the lower triangle of Z' H Z
c     -------------------------------------------------------
        
      do j = 1, nZr

         ZtHZ (j, j) = ZtHZ (j, j) + mu
         
         do i = j+1, nZr
            ZtHZ (i, j) = ZtHZ (j, i)
         enddo
         
      enddo

      return

60000 format (/ 'Indefinite Quadratic Program -- Modifying ',
     1          'Levenberg Parameter'
     2        // '    Previous Levenberg Parameter:', 1pd12.3
     3        /  ' Increment to Levenberg Parameter:', d12.3
     4        /  '          New Levenberg Parameter:', d12.3
     5        // 'Eigenanalysis for Old Reduced Hessian'
     6        / ' Rightmost Eigenvalue:', d12.3
     7        / '  Leftmost Eigenvalue:', d12.3
     8        / ' New Condition Number:', d12.3/ )

61000 format (/ 'Negative definite 1 by 1 Quadratic Program --', 
     1          ' Modifying Levenberg Parameter'
     2        // '    Previous Levenberg Parameter:', 1pd12.3
     3        /  ' Increment to Levenberg Parameter:', d12.3
     4        /  '          New Levenberg Parameter:', d12.3
     5        // 'Eigenvalue of Old Reduced Hessian:', d12.3 / )

c     end of QPMDLV
      
      end
