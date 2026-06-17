      subroutine   xdslv7   ( unsym , neqns , nsuper, nnzero, 
     1                        xsup  , xlindx, lindxg, perm  ,
     2                        coord1, coord2, diag  , values, 
     3                        xarndx, temp  , 
     4                        tnorma, badrow, badcol, error )
c
c
c     ==================================================================
c     ====  xdslv7 -- in memory version of the construction of the  ====
c     ====            relative assembly indicies for the matrix     ====
c     ====            values during the factorization               ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xdslv7 constructs the relative assembly indices for the matrix
c     indices.
c
c     created         08-jan-97   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     unsym       l   symmetric/unsymmetric logical flag
c     neqns       i   number of equations
c     nsuper      i   number of super nodes
c     nnzero      i   size of coord2 and coord1 arrays.
c     xsup        i   supernode partition array
c     xlindx      i   pointer array into the global indices
c     lindxg      i   array of global indices
c     perm        i   new to old permutation array
c     diag        d   array of permuted diagonal entries
c     values      d   array of permuted matrix values
c
c     input/output arguments
c     ----------------------
c
c     coord1,
c       coord2    i   lists of coordinates
c                     coord1 is overwritten with the relative assembly
c                     indices
c
c     working storage
c     ---------------
c
c     temp        d   array of length neqns used to accumulate
c                     abs. row sums for matrix 1-norm computation
c
c     output arguments
c     ----------------
c
c     xarndx      i   pointer array into relative assembly indices
c     tnorma      d   1 norm of papt
c     badrow      i   if error = -1, badrow contains the bad row 
c                     number.  otherwise zero.
c     badcol      i   if error = -1, badrow contains the bad column
c                     number.  otherwise zero.
c     error       i   error flag, = -1 if no entry in the front is
c                                      found
c                                 = -2 if internal error in xdslv9
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------
 
      integer             neqns , nsuper, nnzero, 
     1                    badrow, badcol, error

      logical             unsym
 
      integer             xsup   (*), xlindx (*), lindxg (*),
     1                    perm   (*),
     2                    coord1 (*), coord2 (*), xarndx (*)

      double precision    tnorma

      double precision    diag   (*), values (*), temp   (*)

c     -------------------
c     ... local variables
c     -------------------
 
      integer             i     , ibgn  , iend  , isuper, jbgn  , 
     1                    jcol  , jend  , k     , 
     2                    newchv, nnz   ,oldchv

c     ==================================================================

      error  = 0

c     ------------------------------------------
c     ... prepare temp for computing column sums
c     ------------------------------------------

      do jcol = 1, neqns
          temp(jcol) = abs ( diag(jcol) )
      enddo

c     ----------------
c     ... build xarndx
c     ----------------

      oldchv = 0

      do k = 1, nnzero

          if ( coord1(k) .ne. oldchv ) then

              newchv = coord1(k)
              do i = oldchv+1, newchv
                  xarndx(i) = k
              enddo
              oldchv = newchv

          end if

      enddo

      do i = oldchv+1, neqns+1
          xarndx(i) = nnzero + 1
      enddo

c     -----------------------------
c     ... loop over all super nodes
c     -----------------------------

      k = 1

      do isuper = 1, nsuper

          jbgn   = xsup(isuper)
          jend   = xsup(isuper+1) - 1

          ibgn   = xlindx(isuper)
          iend   = xlindx(isuper+1) - 1

          nnz    = xarndx(jend+1) - xarndx(jbgn)

c.debug
c     write(6,'("in xdslv7 at start of 600 loop")')
c     write(6,'("isuper, k, nnz         = ", 3i8)') 
c    1            isuper, k, nnz        
c     write(6,'("jbgn  , jend           = ", 3i8)') 
c    1            jbgn  , jend          
c     write(6,'("ibgn  , iend           = ", 3i8)') 
c    1            ibgn  , iend          
c.debug

c         ---------------------------------------------------
c         ... compute the relative indices for this supernode
c         ---------------------------------------------------

          call xdslv9 ( unsym, nnz, jbgn, jend, ibgn, iend,
     1                  coord1(k), coord2(k), values(k),
     2                  temp, lindxg, perm, badrow, badcol, error )

          if ( error .ne. 0 ) return

          k = k + nnz

      enddo

c     ------------------------------------------
c     ... complete computation of 1-norm of papt
c     ------------------------------------------

      tnorma = 0.0d0
      do i=1,neqns
        tnorma = max(tnorma,abs(temp(i)))
      enddo
c.debug
c     write(6,'("1 norm of papt = ", 1pd15.5)') tnorma
c.debug

c-----------------------------------------------------------------------

c     -----------------                                                  
c     ... normal return
c     -----------------

      return            
      end
