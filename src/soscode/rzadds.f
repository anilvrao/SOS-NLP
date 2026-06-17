      subroutine   RZADDS   ( unitQ, vertex,
     1                        k1, k2, iT, nactiv, nartif, nZ, nfree,
     2                        nrejtd, ngq, n, ldQ, ldA, ldT,
     3                        istate, kactiv, kx,
     4                        condmx,
     5                        A, T, gqm, Q, w, c, s,
     6                        rthuge, iPrint,
     7                        epspt9, Asize, dTmax, dTmin )
c     ==================================================================
c     ==================================================================
c     ====  rzadds / RZADDS  -- generate TQ factorization           ====
c     ==================================================================
c     ==================================================================

      integer            k1, k2, it, nactiv, nartif, nZ, nfree,
     2                   nrejtd, ngq, n, ldQ, ldA, ldT, iPrint
      
      logical            unitQ, vertex

      double precision   condmx, epspt9, Asize, dTmax, dTmin, rthuge
      
      integer            istate (*), kactiv (n), kx (n)
      
      double precision   A (ldA,*), T (ldT,*), gqm (n,*), Q (ldQ,*)
      
      double precision   w (n), c (n), s (n)
                                 
c     ==================================================================
c
c     derived from qpopt version 1.0
c     last modification -- 17-April-1996
c
c          Original version of  rzadds  written by PEG,  October-31-1984.
c          This version of  rzadds  dated  5-Jul-1991.
c
c     rzadds  includes general constraints  k1  thru  k2  as new rows of
c     the  TQ  factorization:
c              A(free) * Q(free)  = (  0 T )
c                        Q(free)  = (  Z Y )
c
c     a) The  nactiv x nactiv  upper-triangular matrix  T  is stored 
c        with its (1,1) element in position  (iT,jT)  of the array  T.
c                                        
c     ==================================================================

      integer            i     , iadd  , iartif, ifix  , inform, iswap ,
     1                   j     , jadd  , jt    , k    , l    , nzadd
      
      logical            overfl, Rset

      double precision   cndmax, cond, delta, dtnew, Dzz, rowmax,
     1                   rnorm , tdTmax, tdTmin

      double precision   DDIVFN

      external           DDIVFN

      double precision   zero, one
      
      parameter        ( zero = 0.0d0, one = 1.0d0 )
                                                 
c     ==================================================================

      jT     = nZ + 1
                                       
c     Estimate the condition number of the constraints already
c     factorized.

      if  ( nactiv .eq. 0 )  then
         
         dTmax = zero         
         dTmin = one
         
         if  ( unitQ )  then

c           First general constraint added.  Set  Q = I.

            call HDQPII ( ldQ, nfree, Q )
            unitQ  = .false.
            
         end if
         
      else
         
c          << dcond >>
         call DCOND  ( nactiv, T(iT,jT), ldT+1, dTmax, dTmin )

      end if

      do 600, k = k1, k2
         
         iadd = kactiv(k)
         jadd = n + iadd
         
         if  ( nactiv .lt. nfree )  then
            
            overfl = .false.  

c           Transform the incoming row of  A  by  Q'.

            w(1:n) = A(iadd,1:n)
c             << cmqmul >>
            call CMQMUL ( 8, n, nZ, nfree, ldQ, unitQ, kx, w, Q, s )

c           Check that the incoming row is not dependent upon those
c           already in the working set.

            dTnew = zero
            do i=1,nZ
              dTnew = dTnew + w(i)**2
            enddo
            dTnew  = sqrt(dTnew)
            if  ( nactiv .eq. 0 )  then

c              This is the first general constraint in the working set.

               cond   = DDIVFN  ( Asize, dTnew, overfl )
               tdTmax = dTnew
               tdTmin = dTnew

            else

c              There are already some general constraints in the working
c              set. Update the estimate of the condition number.

               tdTmax = max( dTnew, dTmax )
               tdTmin = min( dTnew, dTmin )
c                     << ddiv >>
               cond   = DDIVFN  ( tdTmax, tdTmin, overfl )

            end if

            if  ( cond .ge. condmx  .or.  overfl )  then
               
c              ---------------------------------------------------------
c              This constraint appears to be dependent on those already
c              in the working set.  Skip it.
c              ---------------------------------------------------------
               istate(jadd) =   0
               kactiv(k)    = - kactiv(k)
               
            else
               
               if  ( nZ .gt. 1 )  then
                  
c                 ------------------------------------------------------
c                 Use a single column transformation to reduce the first
c                 nZ-1  elements of  w  to zero.
c                 ------------------------------------------------------
c                 Apply the Householder reflection  I  -  w w'.
c                 The reflection is applied to  Z  and gqm so that
c                    y  =    Z  * w,   Z    =  Z    -  y w'  and
c                    y  =  gqm' * w,   gqm  =  gqm  -  w y',
c                 where  w = wrk1 (from Householder),
c                 and    y = wrk2 (workspace).
c
c                 Note that delta  has to be stored after the reflection
c                 is used.

                  delta = w(nZ)
c                   << dgrfg >>
                  call DGRFG  ( nZ-1, delta, w, 1, zero, w(nZ) )
                  
                  if  ( w(nZ) .gt. zero )  then

                     call dgemv ( 'No transpose', nfree, nZ, one,
     1                            Q, ldQ, w, 1, zero, s, 1 )
                     call dger  ( nfree, nZ, (-one), 
     1                            s, 1, w, 1, Q, ldQ )

                     if  ( ngq .gt. 0 )  then
                        
                        call dgemv ( 'Transpose', nZ, ngq, one, gqm, n, 
     1                               w, 1, zero, s, 1 )
                        call dger  ( nZ, ngq, (-one),
     1                               w, 1, s, 1, gqm, n )
                     end if
                     
                  end if
                                   
                  w(nZ) = delta
                  
               end if
               
               iT     = iT     - 1
               jT     = jT     - 1
               nactiv = nactiv + 1
               nZ     = nZ     - 1
               T(iT,jT:jT+nactiv-1) = w(jT:jT+nactiv-1)
               dTmax  = tdTmax
               dTmin  = tdTmin
               
            end if
            
         end if
         
  600 continue
                        
      if  ( nactiv .lt. k2 )  then

c        Some of the constraints were classed as dependent and not
c        included in the factorization.  Re-order the part of  kactiv
c        that holds the indices of the general constraints in the
c        working set.  Move accepted indices to the front and shift
c        rejected indices (with negative values) to the end.

         l = k1 - 1
         do k = k1, k2
            i         = kactiv(k)
            if  ( i .ge. 0 )  then
               l      = l + 1
               if  ( l .ne. k )  then
                  iswap     = kactiv(l)
                  kactiv(l) = i
                  kactiv(k) = iswap
               end if
            end if
         enddo

c        If a vertex is required,  add some temporary bounds.
c        We must accept the resulting condition number of the working
c        set.

         if  ( vertex )  then
            
            Rset  = .false.
            cndmax = rthuge
            Dzz    = one
            nZadd  = nZ
            
            do iartif = 1, nZadd
               
               if  ( unitQ )  then
                  
                  ifix = nfree
                  jadd = kx(ifix)
                  
               else
                  
                  rowmax = zero
                  do i = 1, nfree
                     rnorm = zero
                     do j=1,nZ
                       rnorm = rnorm + Q(i,j)**2
                     enddo
                     rnorm = sqrt(rnorm)
                     if  ( rowmax .lt. rnorm )  then
                        rowmax = rnorm
                        ifix   = i
                     end if
                  enddo
                  jadd = kx(ifix)

c                   << rzadd >>
                  call RZADD  ( unitQ, Rset,
     1                         inform, ifix, iadd, jadd, iT, 
     2                         nactiv, nZ, nfree, nZ, ngq,
     3                         n, ldA, ldQ, ldT, ldT,
     4                         kx, cndmax, Dzz,
     5                         A, T, T, gqm, Q, w, c, s,
     6                         iPrint, epspt9, Asize, dTmax, dTmin )

c     ******************************************************************
c     ******************************************************************
c     *** note -- inform return appears not to be used!!!!          ****
c     ******************************************************************
c     ******************************************************************

               end if
               
               nfree  = nfree  - 1
               nZ     = nZ     - 1
               nartif = nartif + 1
               istate(jadd) = 4
               
            enddo
            
         end if
         
         if  ( iT .gt. 1 )  then
            
c           -----------------------------------------------
c           If some dependent constraints were rejected,  
c           move the  matrix T  to the top of the array  T.
c           -----------------------------------------------
            
            do k = 1, nactiv 
               j = nZ + k
               do i = 1, k 
                  T(i,j) = T(iT+i-1,j)
               enddo
            enddo
            
         end if
         
      end if

      nrejtd = k2 - nactiv

c     end of RZADDS (rzadds)
      end
