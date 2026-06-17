      subroutine   QPCHL0   ( mnZ   , nZr   , nZ    , ldR   , Hsize , 
     1                        tolrnk, R     , posdef )
c     ==================================================================
c     ==================================================================
c     ====  QPCHL0 /                                                ====
c     ====  qpchl0 -- compute Cholesky factor of reduced Hessian    ====
c     ==================================================================
c     ==================================================================

      integer            mnZ, nZr, nZ, ldR
      
      logical            posdef

      double precision   Hsize, tolrnk
      
      double precision   R (ldR,*)
      
c     ==================================================================
c     derived from qpopt version 1.0 cminit, dated 16-Jan-1995.
c     last modification -- 17-April-1996
c     
c     QPCHL0 / qpchl0  computes the Cholesky factor Rz of the reduced 
c     Hessian  Z'HZ,  given the (nfree x nZ) matrix Z.  If the reduced 
c     Hessian is indefinite, the Cholesky factor of the (nZr x nZr) 
c     matrix  Hz  is returned, where Hz is formed from  H  and  nZr
c     columns of Z.
c
c     ==================================================================

      integer            j

      double precision   d, dmin

      double precision   zero, one, frac

      parameter        ( zero = 0.0d0, one = 1.0d0 )

c     ==================================================================

      nZr    = 0
      posdef = .true.
      
      if  ( nZ .gt. 0)  then

c        ----------------------------------------------
c        ... Form the Cholesky factorization R'R = Z'HZ
c        ----------------------------------------------

         dmin = tolrnk*Hsize

         do j = 1, mnZ

c           See if the diagonal is big enough.

            if  ( R (j,j) .gt. dmin )  then

c              Set the diagonal element of R.

               d      = sqrt ( R (j,j) )
               R(j,j) = d
               nZr    = nZr   + 1

               if  ( j .lt. mnZ )  then

c                 Set the super-diagonal elements of this row of R and 
c                 update the elements of the Schur complement.

                  frac = one/d
                  R(j,j+1:mnZ) = frac*R(j,j+1:mnZ)
                  call dsyr  ( 'Upper', mnZ-j, (-one),
     1                          R (j  , j+1), ldR,
     2                          R (j+1, j+1), ldR )

               end if

            else

               posdef = .false.
               go to 200

            endif
            
         enddo

      endif

 200  continue
      return

c     end of  QPCHL0 (qpchl0)
      
      end
