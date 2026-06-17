      subroutine   CMWRAP   ( nfree, ldA,
     1                        n, nclin, nctotl,
     2                        nactiv, istate, kactiv, kx,
     3                        A, bl, bu, c, clamda, featol,
     4                        r, rlamda, x )
c     ==================================================================
c     ==================================================================
c     ====  CMWRAP /                                                =====
c     ====  cmwrap -- lagrange multiplier data management           ====
c     ==================================================================
c     ==================================================================

      integer            nfree, ldA, n, nclin, nctotl, nactiv
      
      integer            istate (nctotl), kactiv (n), kx (n)
      
      double precision   A (ldA,*), bl (nctotl), bu (nctotl), c (*),
     1                   clamda (nctotl),featol (nctotl), r (nctotl)
      
      double precision   rlamda (n), x (n)

c     ==================================================================
c     derived from qpopt version 1.0 cmwrap
c     last modification -- 26-March-1996
c
c         Original Fortran 77 version written  05-May-93.
c         This version of  cmwrap  dated  05-May-93.
c
c     cmwrap  creates the expanded Lagrange multiplier vector clamda.
c     and resets istate for the printed solution.
c
c     This version of cmwrap is for upper triangular T.
c     For npopt, qpopt and lpopt, kactiv holds the general constraint
c     indices in the order in which they were added to the working set.
c     The reverse ordering is used for T since new rows are added at
c     the front of T.  
c
c     ==================================================================

      integer            i, is, j, k, nfixed, nplin, nZ

      double precision   b1, b2, rj, rlam, slk1, slk2, tol
      
      double precision   zero
      
      parameter         (zero  = 0.0d0)

c     ==================================================================

      nfixed = n     - nfree
      nplin  = n     + nclin
      nZ     = nfree - nactiv

c     Expand multipliers for bounds, linear and nonlinear constraints
c     into the  clamda  array.

      clamda(1:nctotl) = zero
      do k = 1, nactiv+nfixed
         if  ( k .le. nactiv )  then
            j    = kactiv(k) + n
            rlam = rlamda(nactiv-k+1)
         else
            j    = kx(nz+k)
            rlam = rlamda(k)
         end if
         clamda(j) = rlam
      enddo

c     Reset isate if necessary.

      do j = 1, nctotl
         b1     = bl(j)
         b2     = bu(j)

         if  ( j .le. n )  then
            rj  = x(j)
         else 
     1   if  ( j .le. nplin )  then
            i   = j - n
            rj  = dot_product(A(i,1:n),x(1:n))
         else
            i   = j - nplin
            rj  = c(i)
         end if

         is     = istate(j)
         slk1   = rj - b1
         slk2   = b2 - rj
         tol    = featol(j)
         
         if  (                   slk1 .lt. -tol)  is = - 2
         if  (                   slk2 .lt. -tol)  is = - 1
         if  ( is .eq. 1  .and.  slk1 .gt.  tol)  is =   0
         if  ( is .eq. 2  .and.  slk2 .gt.  tol)  is =   0
         
         istate(j) = is
         r(j)      = rj
      enddo

c     end of CMWRAP (cmwrap)
      end
