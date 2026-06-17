      subroutine   RZDEL    ( unitQ, iT, 
     1                       n, nactiv, nfree, ngq, nZ, nZr,
     2                       ldA, ldQ, ldT,
     3                       jdel, kdel, kactiv, kx,
     4                       A, T, gqm, Q, c, s,
     5                       dTmax, dTmin )

c     ==================================================================
c     ==================================================================
c     ====  rzdel / RZDEL   -- Delete Constraint,                   ====
c     ====                     update TQ factorization              ====
c     ==================================================================
c     ==================================================================

      integer            iT, n, nactiv, nfree, ngq, nZ, nZr,
     1                   ldA, ldQ, ldT, jdel, kdel 
      
      logical            unitQ

      double precision   dTmax, dTmin
      
      integer            kactiv (n), kx (n)
      
      double precision   A (ldA,*), T (ldT,*), gqm (n,*), Q (ldQ,*),
     1                   c (n), s (n)

c     ==================================================================
c
c     derived from qpopt version 1.0
c     last modification -- 25-March-1996
c
c          Original version of  rzdel  written by PEG,  31-October-1984.
c          This version of  rzdel  dated 14-Sep-1992.
c           
c     RZDEL / RZDEL   updates the matrices  Z, Y  and  T  associated 
c             with the factorizations 
c     
c              A(free) * Q(free)  = (  0 T )
c                        Q(free)  = (  Z Y )
c
c     when a regular, temporary or artificial constraint is deleted 
c     from the working set.
c
c     The  nactiv x nactiv  upper-triangular matrix  T  is stored 
c     with its (1,1) element in position  (iT,jT)  of the array  T.  
c           
c     Original version of  rzdel  written by PEG,  31-October-1984.
c     This version of  rzdel  dated 14-Sep-1992.
c           
c     ==================================================================

      integer            i    , ir   , itdel, j    , jart , jt   ,
     1                   k    , l    , nsup , npiv , nzR1

      double precision   cs, sn, gmax
      
      double precision   zero, one
      parameter        ( zero = 0.0d0, one = 1.0d0 )

c     ==================================================================

      jT = nZ + 1

      if  ( jdel .gt. 0 )  then
                                       
c        Regular constraint or temporary bound deleted.

         if  ( jdel .le. n )  then

c           Case 1.  A simple bound has been deleted.
c           =======  Columns  nfree+1  and  ir  of gqm' must be swapped.
             
            ir     = nZ    + kdel
            iTdel  = nactiv + 1
            nfree  = nfree  + 1
            if  ( nfree .lt. ir )  then
               kx(ir)    = kx(nfree)
               kx(nfree) = jdel
               call dswap ( ngq, gqm(nfree,1), n, gqm(ir,1), n )
            end if

            if  ( .not. unitQ )  then

c              Copy the incoming column of  A(free)  into the end of T.

               do k = 1, nactiv
                  i = kactiv(k)
                  T(nactiv-k+1,nfree) = A(i,jdel)
               enddo

c              Expand  Q  by adding a unit row and column.

               if  ( nfree .gt. 1 )  then
                 Q(nfree,1:nfree-1) = zero
                 Q(1:nfree-1,nfree) = zero
               end if
               Q(nfree,nfree) = one
            end if
            
         else

c           Case 2.  A general constraint has been deleted.
c           =======
                                  
c           Delete row  iTdel  of  T  and move up the ones below it.
c           T  becomes lower Hessenberg.

            iTdel = kdel 
            do k = iTdel, nactiv
               j  = jT + k - 1
               do l = iTdel, k-1
                  i      = iT + l - 1
                  T(i,j) = T(i+1,j)
               enddo
            enddo

            do i = nactiv-iTdel+1, nactiv-1
               kactiv(i) = kactiv(i+1)
            enddo
            nactiv = nactiv - 1
         end if

         nZ    = nZ + 1

         if  ( nactiv .eq. 0 )  then
            
            dTmax = one
            dTmin = one
            
         else
            
c           ------------------------------------------------------------
c           Restore the nactiv x (nactiv+1) upper-Hessenberg matrix  T
c           to upper-triangular form.  The  nsup  super-diagonal 
c           elements are removed by a backward sweep of rotations.  
c           The rotation for the  (1,1)-th  element of  T  is generated
c           separately.
c           ------------------------------------------------------------

            nsup   = iTdel - 1
            
            if  ( nsup .gt. 0 )  then
               npiv   = jT + iTdel - 1
               if  ( nsup .gt. 1 )  then
                  
                  call dcopy ( nsup-1, T(iT+1,jT+1), ldT+1, s(jT+1), 1 )
c                   << drhstu >>
                  call DRHSTU ( 'right', nactiv, 1, nsup,
     1                          c(jT+1), s(jT+1), T(iT,jT+1), ldT )

               end if

c                << drotgc >>
               call DROTGC ( T(iT,jT+1), T(iT,jT), cs, sn )
               
               T(iT,jT) = zero
               s(jT)    = - sn
               c(jT)    =   cs
               
c                << dgesrc >>
               call DGESRC ( 'right', 'variable', 'backwards', 
     1                      nfree, nfree, nZ, npiv, c, s, Q, ldQ )
               call DGESRC ( 'left ', 'variable', 'backwards', 
     1                      npiv , ngq  , nZ, npiv, c, s, gqm, n   )
            end if

            jT = jT + 1            
c             << dcond >>
            call DCOND  ( nactiv, T(iT,jT), ldT+1, dTmax, dTmin )
            
         end if
         
      end if

      nZr1 = nZr + 1

      if  ( nZ .gt. nZr )  then
         
         if  ( jdel .gt. 0 )  then
            jArt = nZr1
            gmax = abs(gqm(nZr1,1))
            do i=nZr1+1,nZ
              if (abs(gqm(i,1)).gt.gmax) then
                jArt = i
                gmax = abs(gqm(i,1))
              endif
            enddo
         else
            jArt = - jdel
         end if
          
         if  ( jArt .gt. nZr1 )  then

c           Swap columns  nZr1  and  jArt  of  Q  and  gqm.

            if  ( unitQ )  then
               k        = kx(nZr1)
               kx(nZr1) = kx(jArt)
               kx(jArt) = k
            else
               call dswap ( nfree, Q(1,nZr1), 1, Q(1,jArt), 1 )
            end if

            call dswap ( ngq, gqm(nZr1,1), n, gqm(jArt,1), n )
         end if
         
      end if

      nZr = nZr1

c     end of RZDEL  (rzdel)
      end
