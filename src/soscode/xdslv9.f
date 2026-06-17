      subroutine   xdslv9   ( unsym , nnz   , jbgn  , jend  , ibgn  , 
     1                        iend  , coord1, coord2, values, temp  ,
     2                        lindxg, perm  , badrow, badcol, error )
c
c
c     ==================================================================
c     ====  xdslv9 -- in memory version of the construction of the  ====
c     ====            relative assembly indicies for the matrix     ====
c     ====            values during the factorization               ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xdslv9 constructs the relative assembly indices for the matrix
c     indices.
c
c     created         08-jan-97   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     unsym       l   symmetric/unsymmetric logical flag
c     nnz         i   number of org. nonzeroes in the current supernode
c     jbgn        i   the first column for the supernode
c     jend        i   the last  column for the supernode
c     ibgn        i   the first supernodal index for the supernode
c     iend        i   the last  supernodal index for the supernode
c     values      d   array of permuted matrix values for the supernode
c     lindxg      i   array of global indices
c     perm        i   new to old permutation array
c
c     input/output arguments
c     ----------------------
c
c     coord1,
c       coord2    i   lists of coordinates for the supernode
c                     coord1 is overwritten with the relative assembly
c                     indices
c     temp        d   array of length neqns used to accumulate
c                     abs. row sums for matrix 1-norm computation
c
c     working storage
c     ---------------
c
c
c     output arguments
c     ----------------
c
c     badrow      i   if error = -1, badrow contains the bad row 
c                     number.  otherwise zero.
c     badcol      i   if error = -1, badrow contains the bad column
c                     number.  otherwise zero.
c     error       i   error flag, = -1 if no entry in the front is
c                                      found
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------
 
      integer             nnz   , jbgn  , jend  , ibgn  , iend  ,
     1                    badrow, badcol, error

      logical             unsym
 
      integer             coord1 (*), coord2 (*), lindxg (*),
     1                    perm   (*)

      double precision    values (*), temp   (*)

c     -------------------
c     ... local variables
c     -------------------
 
      integer             index , irow  , jchvrn, jcol  , k     , 
     1                    kirow , kjcol , lenext, lenint, length, 
     2                    lirow , ljcol , matsiz, offset
 
c     ==================================================================

      error  = 0
      if ( nnz .eq. 0 ) return

      k      = 1
      offset = 0

      lenint = jend - jbgn + 1
      lenext = iend - ibgn + 1
      length = lenint + lenext 

      matsiz = length * ( length + 1 ) / 2

      kirow  = 0
      kjcol  = 0
c.debug
c     write(6,'("at start of xdslv9")')
c     write(6,'("jbgn, jend, ibgn, iend, length, nnz = ", 6i8)')
c    1            jbgn, jend, ibgn, iend, length, nnz 
c     write(6,'("matsiz                              = ", 6i8)')
c    1            matsiz                              
c     call xislp3 ( 'coord1', nnz, coord1, 6 )
c     call xislp3 ( 'coord2', nnz, coord2, 6 )
c     call xislp3 ( 'lindxg', iend-ibgn+1, lindxg(ibgn), 6 )
c.debug

c     --------------------------------------
c     ... for each chevron of the super node
c     --------------------------------------

      do 500 jchvrn = jbgn, jend
c.debug
c     write(6,'("at start of 500 loop")')
c     write(6,'("jchvrn  , k, offset, length = ", 4i8)') 
c    1            jchvrn  , k, offset, length
c.debug

c         ------------------------------------------------------
c         ... for each entry in the chevron corresponding
c             to the current column.
c             note that coord1 is min ( jcol, irow )
c                  and  coord2 is jcol - irow (abs if symmetric)
c         ------------------------------------------------------

  100     continue
c.debug
c     write(6,'("at 100 continue")')
c     write(6,'("k, coord1(k), jchvrn        = ", 4i8)') 
c    1            k, coord1(k), jchvrn
c.debug

          if ( coord1(k) .ne. jchvrn ) go to 400

c         -------------------------------------------
c         ... decode coordinates to get irow and jcol
c         -------------------------------------------

          if ( unsym ) then

              if ( coord2(k) .le. 0 ) then
                  jcol = jchvrn
                  irow = jchvrn - coord2(k)
              else
                  irow = jchvrn
                  jcol = coord2(k) + jchvrn
              end if

          else

              jcol = jchvrn
              irow = coord2(k) + jchvrn

          end if

c.debug
c     write(6,'("after 100 continue")')
c     write(6,'("jchvrn, coord2(k), jcol  , irow  , k = ", 5i8)') 
c    1            jchvrn, coord2(k), jcol  , irow  , k             
c.debug

c         ----------------------------------------
c         ... add contribution to abs. column sums
c         ----------------------------------------

          if ( unsym ) then
              temp(irow) = temp(irow) + abs ( values(k) )
          else
              temp(irow) = temp(irow) + abs ( values(k) )
              temp(jcol) = temp(jcol) + abs ( values(k) )
          end if

c         --------------------------------------------------
c         ... find local row and column numbers in the front
c         --------------------------------------------------

          if ( irow .le. jend ) then
              lirow = irow - jbgn + 1
          else
c.debug
c     write(6,'("before call to xislbs for irow = ", i8)') irow
c     call xislp3 ( 'lindxg', lenext, lindxg(ibgn), 6 )
c.debug
              call xislbs ( irow, lenext, lindxg(ibgn), kirow ) 
c.debug
c     write(6,'("after call to xislbs, kirow    = ", i8)') kirow
c.debug
              if ( kirow .eq. 0 ) go to 800
              lirow = kirow + lenint
          end if

          if ( jcol .le. jend ) then
              ljcol = jcol - jbgn + 1
          else
c.debug
c     write(6,'("before call to xislbs for jcol = ", i8)') jcol
c     call xislp3 ( 'lindxg', lenext, lindxg(ibgn), 6 )
c.debug
              call xislbs ( jcol, lenext, lindxg(ibgn), kjcol ) 
c.debug
c     write(6,'("after call to xislbs, kjcol    = ", i8)') kjcol
c.debug
              if ( kjcol .eq. 0 ) go to 800
              ljcol = kjcol + lenint
          end if

c         -----------------------------------
c         ... compute relative assembly index
c         -----------------------------------

          index = offset + abs ( lirow - ljcol ) + 1 

          if ( unsym .and. ( ljcol .gt. lirow ) ) index = index + matsiz

c         ----------------------------------------
c         ... overwrite coord1 with relative index
c         ----------------------------------------

  300     continue
          coord1(k) = index
          if ( .not. unsym ) coord2(k) = irow
c.debug
c     write(6,'("k, irow, jcol, lirow, ljcol, index = ", 6i8)')
c    1            k, irow, jcol, lirow, ljcol, index 
c.debug

          k = k + 1
          if ( k .gt. nnz ) go to 900
          go to 100

c         ---------------------------
c         ... adjust for next chevron
c         ---------------------------

  400     continue
          offset = offset + length
          length = length - 1

  500 continue

c-----------------------------------------------------------------------

c     -----------------                                                  
c     ... normal return
c     -----------------

  900 continue
      if ( k-1 .ne. nnz ) then
         error = -2
c.debug
c        write(6,'("oops at 900 in xdslv9")')
c        write(6,'("k, nnz = ", 2i8)') k, nnz
c.debug
      end if

c     --------------------------------------------------------
c     ... sort entries by index  
c         already sorted in the symmetric case by construction
c     --------------------------------------------------------

c     if ( unsym ) then
c.debug
c         write(6,'("before sort in xdslv9")')
c         write(6,'("nnz = ", i8)') nnz 
c         call xislp3 ( 'coord1', nnz, coord1, 6 )
c         call xislp3 ( 'coord2', nnz, coord2, 6 )
c         call xdslp5 ( 'values', nnz, values, 6 )
c.debug
c         call xdslq3 ( nnz, coord1, coord2, values )
c     end if

c.debug
c     write(6,'("at end of xdslv9")')
c     write(6,'("nnz = ", i8)') nnz 
c     call xislp3 ( 'coord1', nnz, coord1, 6 )
c     call xislp3 ( 'coord2', nnz, coord2, 6 )
c     call xdslp5 ( 'values', nnz, values, 6 )
c.debug

      return            

c-----------------------------------------------------------------------

c     ---------------------------------
c     ... error return, entry not found
c     ---------------------------------

  800 continue
      error  = -1
      badrow = perm(irow)
      badcol = perm(jcol)
c.debug
c     write(6,'("entry not found in xdslv9")')
c     write(6,'("jcol, irow = ", 3i8)')
c    1            jcol, irow 
c     write(6,'("badrow, badcol     = ", 3i8)')
c    1            badrow, badcol     
c     write(6,'("ibgn, iend         = ", 3i8)')
c    1            ibgn, iend 
c     call xislp3 ( 'section of lindxg', 
c    1              iend-ibgn+1, lindxg(ibgn), 6 )
c.debug

c-----------------------------------------------------------------------

      return
      end
