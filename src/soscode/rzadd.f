
      subroutine   RZADD    ( unitQ, Rset,
     1                       inform, ifix, iadd, jadd, iT, 
     2                       nactiv, nZ, nfree, nZr, ngq,
     3                       n, ldA, ldQ, ldR, ldT,
     4                       kx, condmx, Dzz,
     5                       A, R, T, gqm, Q,
     6                       w, c, s,
     7                       iPrint, epspt9, Asize, dTmax, dTmin )
c     ==================================================================
c     ==================================================================
c     ====  rzadd / RZADD  -- add constraint to TQ factorization    ====
c     ==================================================================
c     ==================================================================

      integer            inform, ifix, iadd, jadd, iT, nactiv, nZ,
     1                   nfree, nZr, ngq, n, ldA, ldQ, ldR, ldT, iPrint
      
      logical            unitQ, Rset

      double precision   condmx, Dzz, epspt9, Asize, dTmax, dTmin
      
      integer            kx (n)
      
      double precision   A (ldA,*), R (ldR,*), T (ldT,*),
     1                   gqm (n,*), Q (ldQ,*)

      double precision   w(n), c(n), s(n)

c     ==================================================================
c
c     derived from qpopt version 1.0
c     last modification -- 26-March-1996
c     
c         Original version of  rzadd  written by PEG,  31-October-1984.
c         This version of  rzadd  dated  28-Aug-1991.                  
c
c     rzadd/RZADD   updates the matrices  Z, Y, T, R  and  D  
c                   associated with the factorizations 
c
c              A(free) * Q(free)  = (  0 T )
c                        Q(free)  = (  Z Y )
c                      R' *D * R  =   Hz
c
c     a) The matrices  R  and  T  are upper triangular.
c     b) The arrays  T  and  R  may be the same array.
c     c) The  nactiv x nactiv  upper-triangular matrix  T  is stored 
c        with its (1,1) element in position  (iT,jT) of the 
c        array  T.   The integer  jT  is always  nZ+1.  During regular 
c        changes to the working set,  iT = 1;  when several constraints 
c        are added simultaneously,  it  points to the first row of the
c        existing  T. 
c     d) The matrix  R  is stored in the first  nZr x nZr  rows 
c        and columns of the  nfree x nfree  leading principal submatrix
c        of the array  R.                                    
c     e) If  Rset  is  false,   R  is not touched.
c                  
c     There are three separate cases to consider (although each case
c     shares code with another)...
c
c     (1) A free variable becomes fixed on one of its bounds when there
c         are already some general constraints in the working set.
c
c     (2) A free variable becomes fixed on one of its bounds when there
c         are only bound constraints in the working set.
c                                                                
c     (3) A general constraint (corresponding to row  iadd  of  A) is
c         added to the working set.
c
c     In cases (1) and (2), we assume that  kx(ifix) = jadd.
c     In all cases,  jadd  is the index of the constraint being added.
c
c     If there are no general constraints in the working set,  the
c     matrix  Q = (Z Y)  is the identity and will not be touched.
c
c     If  ngq .gt. 0,  the column transformations are applied to the
c     columns of the  (ngq x n)  matrix  gqm'.
c
c     ==================================================================

      integer            i, j, jt, k, nanew, npiv, nsup
      
      logical            bound , overfl

      double precision   cond, condbd, dtnew, tdtmax, tdtmin

      double precision   zero, one
      
      parameter         (zero = 0.0d0, one = 1.0d0)

      double precision   DDIVFN

      external           DDIVFN
                            
c     ==================================================================

c     If the condition estimator of the updated T is greater than
c     condbd,  a warning message is printed.


      condbd = one / epspt9

      overfl = .false.
      bound  = jadd .le. n 
      jT     = nZ + 1

      
      if  ( bound )  then
         
c        ===============================================================
c        A simple bound has entered the working set.  iadd is not used.
c        ===============================================================
         
         nanew = nactiv

         if  ( unitQ )  then

c           Q is not stored, but  kx  defines an ordering of the columns
c           of the identity matrix that implicitly define Q.
c           Define the sequence of pairwise interchanges P that moves
c           the newly-fixed variable to position  nfree.
c           Reorder  kx  accordingly.

            do i = 1, nfree-1
               if  ( i .ge. ifix )  then
                  w (i) = i + 1
                  kx(i) = kx(i+1)
               else
                  w(i) = i
               end if        
            enddo
            
         else
            
c           ------------------------
c           Q  is stored explicitly.
c           ------------------------
            
c           Set  w = the  (ifix)-th  row of  Q.
c           Move the  (nfree)-th  row of  Q  to position ifix.

            w(1:nfree) = Q(ifix,1:nfree)
            if  ( ifix .lt. nfree )  then
               Q(ifix,1:nfree) = Q(nfree,1:nfree)
               kx(ifix) = kx(nfree)
            end if
            
         end if
         
         kx(nfree) = jadd

         
      else
         
c        =================================================
c        A general constraint has entered the working set.
c        ifix is not used.
c        =================================================

         nanew  = nactiv + 1

c        Transform the incoming row of A by Q'.

         
         w(1:n) = A(iadd,1:n)
c          << cmqmul( >>
         call CMQMUL( 8, n, nZ, nfree, ldQ, unitQ, kx, w, Q, c )

         
c        Check that the incoming row is not dependent upon those
c        already in the working set.
                             
         dTnew = zero
         do i=1,nZ
           dTnew = dTnew + w(i)**2
         enddo
         dTnew  = sqrt(dTnew)
         if  ( nactiv .eq. 0 )  then

c           This is the only general constraint in the working set.

            cond   = DDIVFN  ( Asize, dTnew, overfl )
            tdTmax = dTnew
            tdTmin = dTnew
            
         else

c           There are already some general constraints in the working
c           set.  Update the estimate of the condition number.

            tdTmax = max( dTnew, dTmax )
            tdTmin = min( dTnew, dTmin )
            cond   = DDIVFN  ( tdTmax, tdTmin, overfl )
            
         end if

         if  ( cond .gt. condmx  .or.  overfl)  go to 900

         if  ( unitQ )  then

c           First general constraint added.  Set  Q = I.

            call HDQPII ( ldQ, nfree, Q  )
            unitQ  = .false.
            iT     = 0
            
         end if
         
      end if

         
       if  ( bound )  then
         npiv  = nfree
      else
         npiv  = nZ
      end if      

      if  ( unitQ )  then
         
c        ---------------------------------------------------------------
c        The orthogonal matrix  Q  is not stored explicitly.
c        Apply  P, the sequence of pairwise interchanges that moves the
c        newly-fixed variable to position  nfree.
c        ---------------------------------------------------------------

         if  ( ngq .gt. 0)  then
c             << dgeapr >>
            call DGEAPR ( 'left', 'transpose', nfree-1, w, ngq, gqm, n )
         endif
            
         if  ( Rset )  then

c           Apply the pairwise interchanges to  Rz.
c           The subdiagonal elements generated by this process are
c           stored in  s(ifix), s(2), ..., s(nZr-1).

            nsup = nZr - ifix
c             << dpuths >>
            call DPUTHS ( 'right', nZr, ifix, nZr, s, R, ldR )
            
         end if
         
      else
         
c        ---------------------------------------------------------------
c        The matrix  Q  is stored explicitly.  
c        Define a sweep of plane rotations P such that
c                           Pw = beta*e(npiv).
c        The rotations are applied in the planes (1, 2), (2, 3), ...,
c        (npiv-1, npiv).  The rotations must be applied to Q, gqm', R
c        and T. 
c        ---------------------------------------------------------------

c          << dsrotg >>
         call DSROTG ( 'varble', 'forwrds', npiv-1, w(npiv), w, 1,
     1                 c, s )

         
         if  ( ngq .gt. 0)  then
c             << dgesrc >>
            call DGESRC ( 'left ', 'variable', 'forwards', npiv , ngq,
     1                   1, npiv, c, s, gqm, n )
         endif
         

c           << dgesrc >>
          call DGESRC ( 'right', 'variable', 'forwards', nfree, nfree,
     1                1, npiv, c, s, Q, ldQ )

         
         if  ( Rset )  then
                      
c           Apply the rotations to the triangular part of R.
c           The subdiagonal elements generated by this process are
c           stored in  s(1),  s(2), ..., s(nZr-1).

            nsup = nZr - 1
c             << druths >>
            call DRUTHS ( 'right', nZr, 1, nZr, c, s, R, ldR )
            
         end if
         
      end if

      if  ( Rset )  then
         
c        ---------------------------------------------------------------
c        Eliminate the  nsup  subdiagonal elements of  R  stored in 
c        s(nZr-nsup), ..., s(nZr-1)  with a left-hand sweep of rotations 
c        in planes (nZr-nsup, nZr-nsup+1), ..., (nZr-1, nZr).
c        ---------------------------------------------------------------

c          << drhstu >>
         call DRHSTU ( 'left ', nZr, nZr-nsup, nZr, c, s, R, ldR )

         if  ( nsup .gt. 0  .and.  Dzz .ne. one )  then
            Dzz = c(nZr-1)**2 + Dzz*s(nZr-1)**2

         end if
      end if

      if  ( .not. unitQ )  then
         
         if  ( bound )  then
            
c           ------------------------------------------------------------
c           Bound constraint added.   The rotations affect columns 
c           nZ+1  thru  nfree  of  gqm'  and  T.
c           ------------------------------------------------------------
c           The last row and column of  Q  has been transformed to plus
c           or minus the unit vector  e(nfree).  We can reconstitute the
c           column of gqm' corresponding to the new fixed variable.

            if  ( w(nfree) .lt. zero )  then
               if  ( ngq .gt. 0 )  then
                  gqm(nfree,1:ngq) = -gqm(nfree,1:ngq)
               endif
            end if

            if  ( nactiv .gt. 0 )  then
               
               T(iT,jT-1) = s(jT-1)*T(iT,jT)
               T(iT,jT  ) = c(jT-1)*T(iT,jT)

               if  ( nactiv .gt. 1 )  then
c                   << druths >>
                  call DRUTHS ( 'right', nactiv, 1, nactiv, 
     1                          c(jT), s(jT), T(iT,jT), ldT )
                  call dcopy ( nactiv-1, s(jT), 1, T(iT+1,jT), ldT+1 )
               end if

               jT = jT - 1
c                << dcond >>
               call DCOND  ( nactiv, T(iT,jT), ldT+1, tdTmax, tdTmin )
c                   << ddiv >>
               cond = DDIVFN  ( tdTmax, tdTmin, overfl )
               
            end if
            
         else
            
c           ------------------------------------------------------------
c           General constraint added.  Install  w  at the front of  T. 
c           If there is no room,  shift all the rows down one position. 
c           ------------------------------------------------------------

            iT = iT - 1
            if  ( iT .le. 0 )  then
               iT = 1
               do k = 1, nactiv 
                  j = jT + k - 1
                  do i = k, 1, -1
                     T(i+1,j) = T(i,j)
                  enddo
               enddo
            end if
            
            jT = jT - 1
            T(iT,jT:jT+nanew-1) = w(jT:jT+nanew-1)
            
         end if
         
      end if
   
c     ==================================================================
c     Prepare to exit.  Check the magnitude of the condition estimator.
c     ==================================================================
      
 900  continue
      if  ( nanew .gt. 0 )  then
         
         if  ( cond .lt. condmx  .and.  .not. overfl )  then

c           The factorization has been successfully updated.
               
            inform = 0
            dTmax  = tdTmax
            dTmin  = tdTmin
            if  ( cond .ge. condbd  .and.  iPrint .gt. 0)  then
               write (iPrint, 2000) jadd
            endif
            
         else

c           The proposed working set appears to be linearly dependent.

            inform = 1
            
         end if
      end if

      return

 2000 format(/ ' XXX  Serious ill-conditioning in the working set',
     1         ' after adding constraint ',  i5
     2       / ' XXX  Overflow may occur in subsequent iterations.'//)

c     end of RZADD  (rzadd)
      end
