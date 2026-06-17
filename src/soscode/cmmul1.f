      subroutine   CMMUL1   ( prbtyp, msglvl,
     1                        n     , ldA   , ldT   ,
     2                        nactiv, nfree , nZ    ,
     3                        istate, kactiv, kx    ,
     4                        zerolm, notopt, numinf,
     5                        trusml, smllst, jsmlst, ksmlst,
     6                                tinyst, jtiny , jinf  ,
     7                        trubig, biggst, jbigst, kbigst,
     8                        A     , anorms, gq    ,
     9                        rlamda, T     , wtinf ,
     A                        iPrint )
c     ==================================================================
c     ==================================================================
c     ====  CMMUL1 /                                                ====
c     ====  cmmul1 -- compute lagrange multiplier estimates         ====
c     ==================================================================
c     ==================================================================

      integer            msglvl, n     , ldA   , ldT   , nactiv, nfree ,
     1                   nZ    , notopt, numinf, jsmlst, ksmlst, jtiny ,
     2                   jinf  , jbigst, kbigst, iPrint 
      
      character(len=2)   prbtyp

      double precision   zerolm, trusml, smllst, tinyst, trubig, biggst
      
      integer            istate (*), kactiv (n), kx (n)
      
      double precision   A (ldA,*), anorms (*),
     1                   gq (n), rlamda (n), T (ldT,*), wtinf (*)

c     ==================================================================
c
c
c     derived from qpopt version 1.0
c     last modification -- 17-April-1996
c
c          Original version written 31-October-1984.
c          Based on a version of  lsmuls  dated 30-June-1986.
c          This version of  cmmul1  dated 14-Sep-92.
c
c     cmmul1  first computes the Lagrange multiplier estimates for the
c     given working set.  It then determines the values and indices of
c     certain significant multipliers.  In this process, the multipliers
c     for inequalities at their upper bounds are adjusted so that a
c     negative multiplier for an inequality constraint indicates non-
c     optimality.  All adjusted multipliers are scaled by the 2-norm
c     of the associated constraint row.  In the following, the term
c     minimum refers to the ordering of numbers on the real line,  and
c     not to their magnitude.
c
c     jsmlst          is the index of the constraint whose multiplier is
c                     the minimum of the set of adjusted multipliers
c                     with values less than  small.
c     rlamda(ksmlst)  is the associated multiplier.
c
c     jbigst          is the index of the constraint whose multiplier is
c                     the largest of the set of adjusted multipliers with
c                     values greater than (1 + small).
c     rlamda(kbigst)  is the associated multiplier.
c
c     On exit,  elements  1  thru  nactiv  of  rlamda  contain the
c     unadjusted multipliers for the general constraints.  Elements
c     nactiv  onwards of  rlamda  contain the unadjusted multipliers
c     for the bounds.
c
c     ==================================================================

      integer            i, is, j, k, l, nfixed

      double precision   anormj, blam, rlam, scdlam
      
      double precision   one
      
      parameter        ( one = 1.0d0 )

c     ==================================================================

      nfixed =   n - nfree

      jtiny  =   0
      jsmlst =   0
      ksmlst =   0

      jbigst =   0
      kbigst =   0

c     ------------------------------------------------------------------
c     Compute  jsmlst  for regular constraints and temporary bounds.
c     ------------------------------------------------------------------
c     First, compute the Lagrange multipliers for the general
c     constraints in the working set, by solving  T'*lamda = Y'g.

      if  ( n .gt. nZ )  then
         rlamda(1:n-nZ) = gq(nZ+1:n)
      endif
      
      if  ( nactiv .gt. 0 )  then
         call dtrsv ( 'Upper', 'Transpose', 'No transpose',
     1                nactiv, T(1,nZ+1), ldT, rlamda, 1 )
      endif
      
c     -----------------------------------------------------------------
c     Now set elements  nactiv, nactiv+1,... of  rlamda  equal to
c     the multipliers for the bound constraints.
c     -----------------------------------------------------------------

      do l = 1, nfixed
         j     = kx(nfree+l)
         blam  = rlamda(nactiv+l)
         do k = 1, nactiv
            i    = kactiv(k)
            blam = blam - A(i,j)*rlamda(nactiv-k+1)
         enddo
         rlamda(nactiv+l) = blam
      enddo

c     -----------------------------------------------------------------
c     Find  jsmlst  and  ksmlst.
c     -----------------------------------------------------------------

      do k = 1, n-nZ
         
         if  ( k .gt. nactiv )  then
            j  = kx(nZ+k)
         else
            j  = kactiv(nactiv-k+1) + n
         end if

         is   = istate(j)

         i = j - n
         if  ( j .le. n)   then
            anormj = one
         else
            anormj = anorms(i)
         endif

         rlam = rlamda(k)

c        Change the sign of the estimate if the constraint is in
c        the working set at its upper bound.

         if  ( is .eq. 2 )  then
            rlam =      - rlam
         else
     1   if  ( is .eq. 3 )  then
            rlam =   abs( rlam )
         else
     1   if  ( is .eq. 4 )  then
            rlam = - abs( rlam )
         endif

         if  ( is .ne. 3 )  then
            
            scdlam = rlam*anormj

            if  ( scdlam .lt. zerolm )  then
               
               if  ( numinf .eq. 0) notopt = notopt + 1

               if  ( scdlam .lt. smllst )  then
                  smllst = scdlam
                  trusml = rlamda(k)
                  jsmlst = j
                  ksmlst = k
               end if
               
            else
     1      if  ( scdlam .lt. tinyst )  then
               tinyst = scdlam
               jtiny  = j
            end if

           scdlam = rlam/wtinf(j)
           if  ( scdlam .gt. biggst  .and.  j .gt. jinf )  then
              biggst = scdlam
              trubig = rlamda(k)
              jbigst = j
              kbigst = k
           end if
            
         end if

      enddo

c     -----------------------------------------------------------------
c     If required, print the multipliers.
c     -----------------------------------------------------------------

      if  ( msglvl .ge. 20 )  then
         
         if  ( nfixed .gt. 0)  then
            write (iPrint, 1100) prbtyp, (kx(nfree+k),
     1                           rlamda(nactiv+k), k=1,nfixed)
         endif
         if  ( nactiv .gt. 0)  then
            write (iPrint, 1200) prbtyp, (kactiv(k),
     1                           rlamda(nactiv-k+1), k=1,nactiv)
         end if
         
      end if

      return

 1100 format(/ ' Multipliers for the ', a2, ' bound  constraints   '
     $       / 4(i5, 1pe11.2))
 1200 format(/ ' Multipliers for the ', a2, ' linear constraints   '
     $       / 4(i5, 1pe11.2))

c     end of CMMUL1 (cmmul1)
      
      end
