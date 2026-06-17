      subroutine    CMCRSH   ( start, vertex, 
     1                         nclin, nctotl, nactiv, nartif,
     2                         nfree, n, ldA,
     3                         istate, kactiv, kx, 
     4                         bigbnd, tolact,
     5                         A, Ax, bl, bu, featol, x, wx, work )
c     ==================================================================
c     ==================================================================
c     ====  CMCRSH /                                                ====
c     ====  cmcrsh -- crash start procedure                         ====
c     ==================================================================
c     ==================================================================

      integer            nclin , nctotl, nactiv, nartif, nfree , n     ,
     1                   lda

      character(len=4)   start
      
      logical            vertex

      double precision   bigbnd, tolact
      
      integer            istate (nctotl), kactiv (n), kx (n)
      
      double precision   A (ldA,*), Ax (*), bl (nctotl), bu (nctotl)
      
      double precision   featol (nctotl)
      
      double precision   x (n), wx (n), work (n)

c     ==================================================================
c
c     derived from qpopt version 1.0
c     last modification -- 26-March-1996
c     
c         Original version written by  PEG, 31-October-1984.
c         This version of qpopt cmcrsh  dated 11-May-1995.
c
c     CMCRSH / cmcrsh  computes the quantities  istate (optionally), 
c     kactiv, nactiv,  nz  and  nfree  associated with the working
c     set at x.
c
c     The computation depends upon the value of the input parameter
c     start,  as follows...
c
c     Start = 'cold'  An initial working set will be selected. First,
c                     nearly-satisfied or violated bounds are added.
c                     Next,  general linear constraints are added that
c                     have small residuals.
c
c     Start = 'warm'  The quantities kactiv, nactiv and nfree are
c                     initialized from istate,  specified by the user.
c
c     If vertex is true, an artificial vertex is defined by fixing some
c     variables on their bounds.  Infeasible variables selected for the
c     artificial vertex are fixed at their nearest bound.  Otherwise,
c     the variables are unchanged.
c
c     Values of istate(j)....
c
c        - 2         - 1         0           1          2         3
c     a'x lt bl   a'x gt bu   a'x free   a'x = bl   a'x = bu   bl = bu
c
c     ==================================================================

      integer            i, imin, is, j, k, jmin, jfree, jfix

      double precision   biglow, bigupp, b1, b2, colmin, colsiz, flmax,
     1                   residl, resl, resmin, resu, toobig, tol

      double precision   suminf
      
      double precision   hdmcon

      external           hdmcon
      
      double precision   zero, one
      
      parameter        ( zero = 0.0d0, one = 1.0d0 )

c     ==================================================================

      flmax  =   hdmcon (2)
      
      biglow = - bigbnd
      bigupp =   bigbnd

c     ------------------------------------------------------------------
c     Move the variables inside their bounds.
c     ------------------------------------------------------------------

      
      do j = 1, n
         
         b1    = bl(j)
         b2    = bu(j)
         tol   = featol(j)
            
         if  ( b1 .gt. biglow )  then
            if  ( x(j) .lt. b1 - tol )  x(j) = b1
         end if
            
         if  ( b2 .lt. bigupp )  then
            if  ( x(j) .gt. b2 + tol )  x(j) = b2
         end if
         
      enddo

      wx(1:n) = x(1:n)
      nfree  =   n
      nactiv =   0
      nartif =   0

      if      (start .eq. 'cold' )  then
         
         do j = 1, nctotl
            istate(j) = 0
            if  ( bl(j) .eq. bu(j) )  istate(j) = 3
         enddo

      else if  ( start .eq. 'warm' )  then

C        ... check that existing states make sense (jgl)
         
         do j = 1, nctotl
            if  ( istate(j) .gt. 3  .or.  istate(j) .lt. 0 )  then
               istate(j) = 0
            endif
            
            if  ( bl(j) .le. biglow .and. bu(j) .ge. bigupp )  then
               istate(j) = 0
            endif
            
            if  ( bl(j) .le. biglow .and. istate(j) .eq. 1 )  then
               istate(j) = 0
            endif
            
            if  ( bu(j) .ge. bigupp .and. istate(j) .eq. 2 )  then
               istate(j) = 0
            endif
            
            if  ( bl(j) .ne. bu(j)  .and. istate(j) .eq. 3 )  then
               istate(j) = 0
            endif
            
         enddo

      end if

c     Define nfree and kactiv.
c     Ensure that the number of bounds and general constraints in the
c     working set does not exceed n.

      do j = 1, nctotl
         if  ( nactiv .eq. nfree )  istate(j) = 0

         if  ( istate(j) .gt. 0 )  then
            if  ( j .le. n )  then
               nfree = nfree - 1

               if      (istate(j) .eq. 1 )  then
                  wx(j) = bl(j)
               else if  ( istate(j) .ge. 2 )  then
                  wx(j) = bu(j)
               end if
            else
               nactiv = nactiv + 1
               kactiv(nactiv) = j - n
            end if
         end if
      enddo

c     ------------------------------------------------------------------
c     If a cold start is required,  attempt to add as many
c     constraints as possible to the working set.
c     ------------------------------------------------------------------
      if  ( start .eq. 'cold' )  then
                   
c        See if any bounds are violated or nearly satisfied.
c        If so,  add these bounds to the working set and set the
c        variables exactly on their bounds.

         j = n
         
c+       while (j .ge. 1  .and.  nactiv .lt. nfree) do
 300     continue
         if    (j .ge. 1  .and.  nactiv .lt. nfree )  then
            
            if  ( istate(j) .eq. 0 )  then
               b1     = bl(j)
               b2     = bu(j)
               is     = 0
               
               if  ( b1 .gt. biglow )  then
                  if  ( wx(j) - b1 .le. (one + abs( b1 ))*tolact)  then
                     is = 1
                  endif
               end if
               
               if  ( b2 .lt. bigupp )  then
                  if  ( b2 - wx(j) .le. (one + abs( b2 ))*tolact)  then
                     is = 2
                  endif
               end if
               
               if  ( is .gt. 0 )  then
                  istate(j) = is
                  if  ( is .eq. 1) wx(j) = b1
                  if  ( is .eq. 2) wx(j) = b2
                  nfree = nfree - 1
               end if
            end if
            
            j = j - 1
            go to 300
            
c+       end while
         end if

c        ---------------------------------------------------------------
c        The following loop finds the linear constraint (if any) with
c        smallest residual less than or equal to tolact  and adds it
c        to the working set.  This is repeated until the working set
c        is complete or all the remaining residuals are too large.
c        ---------------------------------------------------------------
c        First, compute the residuals for all the constraints not in the
c        working set.

         suminf = zero
         if  ( nclin .gt. 0  .and.  nactiv .lt. nfree )  then
            do i = 1, nclin
               resl = zero
               resu = zero
               if  ( istate(n+i) .le. 0 )  then
                  Ax (i) = dot_product(A(i,1:n),wx(1:n))
                  if ( bl(i) .gt. biglow )  then
                     resl = min (Ax(i) - bl (i), zero )
                  endif
                  if ( bu(i) .lt. bigupp )  then
                     resu = min ( bu (i) - Ax(i), zero )
                  endif
                  suminf = suminf + resu + resl
               else
               endif
               
            enddo

            is     = 1
            toobig = tolact + tolact

c+          while (is .gt. 0  .and.  nactiv .lt. nfree) do
 500        continue
            if    (is .gt. 0  .and.  nactiv .lt. nfree )  then

C              ... complexity of this search is wrong.  could
C                  sort all residuals smaller than resmin,
C                  looking for the lowest (nfree-nactiv) ones (jgl)
C                  don't need sorted, just lowest ones
               
               is     = 0
               resmin = tolact

               do i = 1, nclin
                  
                  j      = n + i
                  if  ( istate(j) .eq. 0 )  then
                     b1     = bl(j)
                     b2     = bu(j)
                     resl   = toobig
                     resu   = toobig
                     if  ( b1 .gt. biglow)  then
                        resl  = abs( Ax(i) - b1 ) / (one + abs( b1 ))
                     endif
                     
                     if  ( b2 .lt. bigupp)  then
                        resu  = abs( Ax(i) - b2 ) / (one + abs( b2 ))
                     endif
                     
                     residl   = min ( resl, resu )
                     if  ( residl .lt. resmin )  then
                        resmin = residl
                        imin   = i
                        is     = 1
                        if  ( resl .gt. resu )  then
                           is = 2
                        endif
                     end if
                     
                  end if
               enddo

               if  ( is .gt. 0 )  then
                  nactiv = nactiv + 1
                  kactiv(nactiv) = imin
                  j         = n + imin
                  istate(j) = is
               end if
               go to 500
               
c+          end while
            end if
            
         end if              
      end if      

      if  ( vertex  .and.  nactiv .lt. nfree )  then
c        ---------------------------------------------------------------
c        Find an initial vertex by temporarily fixing some variables.
c        ---------------------------------------------------------------
c        Compute lengths of columns of selected linear constraints
c        (just the ones corresponding to variables eligible to be
c        temporarily fixed).        
         
         do j = 1, n             
            if  ( istate(j) .eq. 0 )  then
               colsiz = zero
               do k = 1, nclin
                  if  ( istate(n+k) .gt. 0)  then
                     colsiz = colsiz + abs( A(k,j) )
                  endif
               enddo
               work(j) = colsiz
            end if
         enddo
         
c        Find the  nartif  smallest such columns.
c        This is an expensive loop.  Later we can replace it by a
c        4-pass process (say), accepting the first col that is within
c        t  of  colmin, where  t = 0.0, 0.001, 0.01, 0.1 (say).
c        (This comment written in 1980).
                                    
c+       while (nactiv .lt. nfree) do
 640     continue
         if    (nactiv .lt. nfree )  then
            
            colmin = flmax
            do j = 1, n
               if  ( istate(j) .eq. 0 )  then
                  if  ( nclin .eq. 0)  go to 660
                  colsiz = work(j)
                  if  ( colmin .gt. colsiz )  then
                     colmin = colsiz
                     jmin   = j
                  end if
               end if
            enddo
            
            j         = jmin

c           Fix x(j) at its current value.

  660       continue
            istate(j) = 4
            nartif    = nartif + 1
            nfree     = nfree  - 1
            go to 640
            
c+       end while
         end if
         
      end if
      
      jfree = 1
      jfix  = nfree + 1
      
      do j = 1, n
         if  ( istate(j) .le. 0 )  then
            kx(jfree) = j
            jfree     = jfree + 1
         else
            kx(jfix)  = j
            jfix      = jfix  + 1
         end if
      enddo

      
c     end of CMCRSH (cmcrsh)
      
      end
