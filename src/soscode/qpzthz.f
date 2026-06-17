      subroutine   QPZTHZ   ( unitQ, qpHess, delta, 
     1                        n, nZ, nfree, nHess,
     2                        ldQ, ldH, ldR,
     3                        kx, Hsize, 
     4                        H, R, Q,
     5                        Hz, wrk )

c     ==================================================================
c     ==================================================================
c     ====  QPZTHZ /                                                ====
c     ====  qpzthz -- compute reduced Hessian                       ====
c     ==================================================================
c     ==================================================================

      integer            n, nZ, nfree, nHess, ldQ, ldH, ldR
      
      logical            unitQ

      double precision   delta, Hsize
      
      integer            kx (n)
      
      double precision   H (*), R (ldR,*),
     1                   Q (ldQ,*), Hz (n), wrk (n)
      
      external           qpHess

c     ==================================================================
c     derived from qpopt version 1.0 qpcrsh, dated 16-Jan-1995.
c     last modification -- 25-March-1996
c     
c     QPZTHZ / qpzthz  computes the reduced Hessian Z'HZ,  given the
c     (nfree x nZ) matrix Z.  The reduced Hessian is stored in the
c     leading columns of  R.
c
c     ==================================================================

      integer            i, j, jthcol, k

      double precision   zero, one

      parameter        ( zero = 0.0d0, one = 1.0d0 )


      integer           diagns

      parameter       ( diagns = 0 )

      double precision  maxdff, mxrldf

c     ==================================================================

      Hsize = zero

      if  ( diagns .gt. 0 )  then

c        SPECIAL DIAGNOSTICS FOR USE WITH FTRIHS

         write (*,*)
         write (*,*) 'FTRIHS version of QPZTHZ'
         write (*,*) 'Full Hessian, by rows of the upper triangle'
         j = 1
         do i = 1, n
            write (*,*) 'row ', i, ':',
     1                  (H(k), k = j, j + n - i)
            j = j + (n - i + 1)
         enddo

      endif

      if  ( nZ .gt. 0 )  then

c        -------------------------------------------------------
c        Compute  Z'HZ  and store the upper-triangular symmetric 
c        part in the first  nZ  columns of R.
c        -------------------------------------------------------

         do k = 1, nZ

            wrk(1:n) = zero

c           -----------------------------------
c           ... first find the k-th column of Z
c           -----------------------------------
            
            if  ( unitQ )  then

c              Only bounds are in the working set.  The k-th column of Z
c              is just a column of the identity matrix.

               jthcol       = kx (k)
               wrk (jthcol) = one
               
            else

c              Expand the column of Z into an n-vector.

               do i = 1, nfree
                  j       = kx (i)
                  wrk (j) = Q (i,k)
               enddo
               jthcol = 0
               
            end if

c           -----------------------------------------
c           Set  R(*,k)  =  top of  Q*H*(column of Z)
c           -----------------------------------------

            call qpHess ( n, ldH, nHess, jthcol, H, wrk, delta, Hz )

c             << cmqmul >>
            call CMQMUL ( 4, n, nZ, nfree, ldQ, unitQ, kx, Hz, Q, wrk )

c           ... (save the entire column of Z'HZ)
            
            cnrm = zero
            do i=1,nZ
              R(i,k) = Hz(i)
              cnrm = cnrm + Hz(i)**2
            enddo
            cnrm = sqrt(cnrm)
    
c           Update an estimate of the size of the reduced Hessian.

corig       Hsize  = max ( Hsize, abs( R(k,k) ) )
            Hsize  = max ( Hsize, cnrm )

         enddo

         if  ( diagns .gt. 0 ) then

            maxdff    = 0.0
            mxrldf = 0.0
            
            do j = 1, nZ
               do i = j+1, nz
                  
                  maxdff = max ( maxdff,
     1                            abs ( R(i,j) - R(j,i) ) )

                  if  ( abs (R(i,j)) .ge. one )  then

                     mxrldf = max ( mxrldf,
     1                                  abs ( R(i,j) - R(j,i) ) /
     2                                  max ( abs (R(i,j)),
     3                                        abs (R(j,i)) ) )

                  else

                     mxrldf = max ( mxrldf,
     1                                  abs ( R(i,j) - R(j,i) ) )

                  endif

               enddo
            enddo

            write (*,*)
            write (*,*) 'diagnostics on symmetry of computed reduced', 
     1                  ' hessian'
            write (*,*)
            write (*,*) ' largest norm-wise relative error',
     1                  maxdff / Hsize
            write (*,*) ' largest entry-wise relative error',
     1                  mxrldf

            if ( diagns .gt. 1 )  then

               write (*,*)
               write (*,*) 'Computed Reduced Hessian'
               do i = 1, nZ
                  write (*,*) 'row ', i, ':', (R(i,j), j = 1, nZ)
               enddo

            endif

         endif

c        ... ensure symmetry by copying the upper triangle into the
c            lower triangle

         do i = 1, nZ
            do j = i+1, nZ
               R(i,j) = R(j,i)
            enddo
         enddo
         
      endif
      
      return

c     end of QPZTHZ (qpzthz)
      
      end
