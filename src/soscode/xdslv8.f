      subroutine   xdslv8   ( unsym , savea , neqns , nsuper, lset, 
     1                        xsup  , xlindx, lindxg, perm  ,
     2                        coord1, coord2, diag  , values, 
     3                        xarndx, temp  , tnorma, badrow, badcol, 
     4                        wafil2, watrn2, sqfile, sqlen,  sqtrn , 
     5                        nrecrd, error )
c
c
c     ==================================================================
c     ====  xdslv8 -- out of memory version of the construction of  ====
c     ====            the relative assembly indicies for the matrix ====
c     ====            values during the factorization               ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xdslv8 constructs the relative assembly indices for the matrix
c     indices when the matrix is out of memory.  it prepares the
c     matrix data for the minimum core factorization.
c
c     created         15-jan-97   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     unsym       l   symmetric/unsymmetric logical flag
c     savea       l   save matrix logical flag
c     neqns       i   number of equations
c     nsuper      i   number of super nodes
c     lset        i   size of ccord1, coord2 and values arrays.
c     xsup        i   supernode partition array
c     xlindx      i   pointer array into the global indices
c     lindxg      i   array of global indices
c     perm        i   new to old permutation array
c     diag        d   array of permuted diagonal entries
c     values      d   buffer to hold the permuted matrix values
c     xarndx      i   pointer array into relative assembly indices
c     wafil2      i   word addressable i/o unit containing representation
c                     of the orginial matrix
c     sqfile      i   sequential i/o unit which on output will hold the
c                     matrix representation required for minimum core
c                     factorization
c
c     input/output arguments
c     ----------------------
c
c     coord1,
c       coord2    i   buffers to hold the lists of coordinates
c                     coord1 is overwritten with the relative assembly
c                     indices
c     watrn2      i   amount of i/o transfers on wafil2
c     sqlen       d   length of sqfile
c     sqtrn       d   amount of i/o transfers on sqfile
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
c     tnorma      d   1 norm of papt
c     badrow      i   if error = -1, badrow contains the bad row 
c                     number.  otherwise zero.
c     badcol      i   if error = -1, badrow contains the bad column
c                     number.  otherwise zero.
c     nrecrd      i   number of records actually written to sqfile.
c     error       i   error flag, = -1 if no entry in the front is
c                                      found
c                                 = -2 internal error in xdslv9
c                                 = -3 i/o error on wafil2
c                                 = -4 i/o error on sqfile
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------
 
      integer             neqns , nsuper, lset, 
     1                    badrow, badcol, error , wafil2, watrn2, 
     2                    sqfile, nrecrd

      logical             unsym , savea
 
      integer             xsup   (*), xlindx (*), lindxg (*),
     1                    perm   (*),
     2                    coord1 (*), coord2 (*), xarndx (*)

      double precision    tnorma, sqlen , sqtrn  

      double precision    diag   (*), values (*), temp   (*)

c     -------------------
c     ... local variables
c     -------------------
 
      integer             i     , ibgn  , iend  , innz  , iopos ,
     1                    isuper, jbgn  , jcol  , jend  , k     , 
     2                    len   , nnz   , valbgn, valend, valmax

      integer             xdslni

      external            xdslni
 
c     ==================================================================

      error  = 0
      valend = 0
      valmax = xarndx ( xsup (nsuper+1) ) - 1

c     ------------------------------------------
c     ... prepare temp for computing column sums
c     ------------------------------------------

      do jcol = 1, neqns
          temp(jcol) = abs ( diag(jcol) )
      enddo

c     -----------------------------
c     ... loop over all super nodes
c     -----------------------------

      do 600 isuper = 1, nsuper

          jbgn   = xsup(isuper)
          jend   = xsup(isuper+1) - 1

          ibgn   = xlindx(isuper)
          iend   = xlindx(isuper+1) - 1

          iopos  = xarndx(jbgn)
          nnz    = xarndx(jend+1) - xarndx(jbgn)

          if ( nnz .eq. 0 ) go to 600

          innz   = xdslni(nnz)

c.debug
c     write(6,'("in xdslv8 at start of 600 loop")')
c     write(6,'("isuper, iopos, nnz     = ", 3i8)')
c    1            isuper, iopos, nnz
c     write(6,'("wafil2, sqfile         = ", 3i8)')
c    1            wafil2, sqfile    
c     write(6,'("jbgn  , jend           = ", 3i8)') 
c    1            jbgn  , jend          
c     write(6,'("ibgn  , iend           = ", 3i8)') 
c    1            ibgn  , iend          
c.debug

          if ( iopos + nnz - 1 .gt. valend ) then

c             ----------------------------
c             ... read in data from wafil2
c             ----------------------------

              valbgn = iopos
              len    = min ( lset, valmax - valbgn + 1 )
              valend = valbgn + len - 1
              k      = 1

              call xdslw1 ( wafil2, 3, coord1, coord2,
     1                      values, iopos, len, error )
              if ( error .ne. 0 ) go to 8000 
c.debug
c     write(6,'("after read in xdslv8 - valbgn, valend = ", 2i8)')
c    1                                   valbgn, valend 
c.debug

          end if
 
c         ---------------------------------------------------
c         ... compute the relative indices for this supernode 
c         ---------------------------------------------------
c.debug
c     write(6,'(/, "before xdslv9 - k = ", i8 )') k
c     call xislp3 ( 'coord1', nnz, coord1(k), 6 )
c     call xislp3 ( 'coord2', nnz, coord2(k), 6 )
c     call xdslp5 ( 'values', nnz, values(k), 6 )
c.debug

          call xdslv9 ( unsym, nnz, jbgn, jend, ibgn, iend,
     1                  coord1(k), coord2(k), values(k),
     2                  temp, lindxg, perm, badrow, badcol, error )
          if ( error .ne. 0 ) return
          
          watrn2 = watrn2 + 2*innz + nnz
c.debug
c     write(6,'(/, "data for super node         ", i8)') isuper
c     write(6,'("after xdslv9 - nnz, savea = ", i8, l8 )')
c    1                              nnz, savea
c     call xislp3 ( 'coord1', nnz, coord1(k), 6 )
c     call xislp3 ( 'coord2', nnz, coord2(k), 6 )
c     call xdslp5 ( 'values', nnz, values(k), 6 )
c.debug

c         -------------------------------------------------
c         ... write out data for minimum core factorization
c         -------------------------------------------------

         if ( nnz .gt. 0 ) then
 
              if ( savea ) then 

                  nrecrd = nrecrd + 1
                  call xislvw ( sqfile, nnz, coord2(k), error )
                  if ( error .ne. 0 ) go to 8100

                  sqtrn = sqtrn + innz 
                  sqlen = sqlen + innz 

              end if

              nrecrd = nrecrd + 1
              call xislvw ( sqfile, nnz, coord1(k), error )
              if ( error .ne. 0 ) go to 8100

              nrecrd = nrecrd + 1
              call xdslvw ( sqfile, nnz, values(k), error )
              if ( error .ne. 0 ) go to 8100
c.debug
c             write(6,'("data for super node ", i8)') isuper
c             call xislp3 ( 'relative indicies', nnz, coord1(k), 6 )
c             call xdslp5 ( 'values           ', nnz, values(k), 6 )
c.debug

              sqtrn = sqtrn + innz + nnz
              sqlen = sqlen + innz + nnz

          end if

          k = k + nnz

  600 continue

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

      go to 9000

c-----------------------------------------------------------------------

c     -----------------------
c     ... i/o error on wafil2
c     -----------------------

 8000 continue
      error = -3
      go to 9000

c     -----------------------
c     ... i/o error on sqfile
c     -----------------------

 8100 continue
      error = -4

c-----------------------------------------------------------------------
 9000 continue
      return
      end
