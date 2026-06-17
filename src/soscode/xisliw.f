      subroutine   xisliw   ( neqns , perm1 , perm2, itemp )
c
c
c     ==================================================================
c     ====  xisliw -- generate two random permutatations            ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xisliw generates two random permutations by generating an array
c     full of uniform random numbers and then sorting them.  The sort
c     index will then be the random permutation.
c
c     created         04-dec-96   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     neqns       i   number of entries in the local list
c
c     working storage
c     ---------------
c
c     itemp       i   array to hold random numbers
c
c     output arguments
c     ----------------
c
c     perm1       i   array holding first  random permutation
c     perm2       i   array holding second random permutation
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------
 
      integer             neqns
 
      integer             perm1  (*), perm2  (*), itemp  (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             ier   , iseed , kmax
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      external            hjrsun, hjsrtn
 
c     ==================================================================

      iseed = 137539753
      kmax  = 2 ** 30

      call hjrsun ( iseed, kmax, neqns, itemp )

      call hjsrtn ( itemp, neqns, 0, 1, perm1, ier )

      call hjrsun ( iseed, kmax, neqns, itemp )

      call hjsrtn ( itemp, neqns, 0, 1, perm2, ier )

c     ==================================================================

      return
      end
