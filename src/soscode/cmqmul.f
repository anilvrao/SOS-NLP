      subroutine   CMQMUL   ( mode, n, nz, nfree, ldQ, unitQ,
     1                        kx, v, Q, w )
c     ==================================================================
c     ==================================================================
c     ====  CMQMUL /                                                ====
c     ====  cmqmul -- (partial) Q times vector multiply             ====
c     ==================================================================
c     ==================================================================

      integer            mode, n, nz, nfree, ldQ
      
      logical            unitQ
      
      integer            kx (n)
      
      double precision   v (n), Q (ldQ,*), w (n)

c     ==================================================================
c
c     derived from qpopt version 1.0
c     last modification -- 17-April-1996
c         Original F66 version  April 1983.
c         Fortran 77 version written  9-February-1985.
c         Level 2 BLAS added 10-June-1986.
c         This version of cmqmul dated 20-May-1992.
c     
c     CMQMUL / cmqmul  transforms the vector  v  in various ways using
c     the matrix  Q = ( Z  Y )  defined by the input parameters.
c
c        MODE               result
c        ----               ------
c
c          1                v = Z v
c          2                v = Y v
c          3                v = Q v
c
c     On input,  v  is assumed to be ordered as  ( v(free)  v(fixed) ).
c     On output, v  is a full n-vector.
c
c
c          4                v = Z'v
c          5                v = Y'v
c          6                v = Q'v
c
c     On input,  v  is a full n-vector.
c     On output, v  is ordered as  ( v(free)  v(fixed) ).
c
c          7                v = Y'v
c          8                v = Q'v
c
c     On input,  v  is a full n-vector.
c     On output, v  is as in modes 5 and 6 except that v(fixed) is not
c     set.
c
c     Modes  1, 4, 7 and 8  do not involve  v(fixed).
c     ==================================================================

      integer            j, j1, j2, k, l, lenv, nfixed
      
      double precision   zero, one
      
      parameter        ( zero = 0.0d0, one = 1.0d0 )

c     ==================================================================

      nfixed = n - nfree
      j1     = 1
      j2     = nfree
      
      if  ( mode .eq. 1  .or.  mode .eq. 4)  then
         j2 = nz
      endif
      
      if  ( mode .eq. 2  .or.  mode .eq. 5  .or.  mode .eq. 7)  then
         j1 = nz + 1
      endif
      
      lenv   = j2 - j1 + 1
      
      if  ( mode .le. 3 )  then
         
c        ===============================================================
c        Mode = 1, 2  or  3.
c        ===============================================================

         if  ( nfree .gt. 0 )  then
           w(1:nfree) = zero
         endif

c        Copy  v(fixed)  into the end of  wrk.

         if  ( mode .ge. 2  .and.  nfixed .gt. 0)  then
            w(nfree+1:nfree+nfixed) = v(nfree+1:nfree+nfixed)
         endif

c        Set  W  =  relevant part of  Q * V.

         if  ( lenv .gt. 0)  then
            if  ( unitQ )  then
               w(j1:j1+lenv-1) = v(j1:j1+lenv-1)
            else
               call dgemv ( 'No transpose', nfree, j2-j1+1, one, 
     1                      Q(1,j1), ldQ, v(j1), 1, one, w, 1 )
            end if
         end if

c        Expand  w  into  v  as a full n-vector.

         v(1:n) = zero
         do k = 1, nfree
            j      = kx(k)
            v(j)   = w(k)
         enddo

c        Copy  w(fixed)  into the appropriate parts of  v.

         if  ( mode .gt. 1)  then
            do l = 1, nfixed
               j       = kx(nfree+l)
               v(j)    = w(nfree+l)
            enddo
         end if

      else
         
c        ===============================================================
c        Mode = 4, 5, 6, 7  or  8.
c        ===============================================================
c        Put the fixed components of  v  into the end of  w.

         if  ( mode .eq. 5  .or.  mode .eq. 6)  then
            do l = 1, nfixed
               j          = kx(nfree+l)
               w(nfree+l) = v(j)
            enddo
         end if

c        Put the free  components of  v  into the beginning of  w.

         if  ( nfree .gt. 0 )  then
            
            do k = 1, nfree
               j      = kx(k)     
               w(k) = v(j)
            enddo

c           Set  v  =  relevant part of  Q' * w.

            if  ( lenv .gt. 0)  then
               if  ( unitQ )  then
                  v(j1:j1+lenv-1) = w(j1:j1+lenv-1)
               else
                  call dgemv ( 'Transpose', nfree, j2-j1+1, one, 
     1                         Q(1,j1), ldQ, w, 1, zero, v(j1), 1 )
               end if
               
            end if
            
         end if

c        Copy the fixed components of  w  into the end of  v.

         if  ( nfixed .gt. 0  .and.  (mode .eq. 5  .or.  mode .eq. 6))
     1   then
            v(nfree+1:nfree+nfixed) = w(nfree+1:nfree+nfixed)
         endif
         
      end if

c     end of  CMQMUL (cmqmul)
      
      end
