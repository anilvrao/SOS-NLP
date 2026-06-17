      subroutine   xisliu   ( neqns , nnzero, maxnz , wafile, output,
     1                        colptr, rowlst, collst, ladjnc,
     1                        itemp1, itemp2, itemp3, itemp4, itemp5,
     2                        itemp6, itemp7, 
     3                        ncomp , nzcomp, error  )
c
c
c     ==================================================================
c     ====  xisliu -- builds compressed adjacency structure from    ====
c     ====            rowlst and collst structure out of memory     ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xisliu builds a compressed adjacency structure from the
c     out-of-memory representation of rowlst and collst.
c
c     created         12-dec-96   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     neqns       l   order of uncompressed problem.
c     nnzero      i   the number of row and columns entries on wafile
c     maxnz       i   size of rowlst and collst arrays.
c     wafile      i   i/o file holding row and column entries
c     output      i   output i/o file
c     itemp1      i   array of length neqns holding first random
c                     permutation for use in computing second checksum
c     itemp2      i   array of length neqns holding second random
c                     permutation for use in computing first checksum
c
c     working vectors 
c     ---------------
c
c     colptr      i   pointer to the start of each column in 
c                     rowlst/collst on i/o file wafile
c     rowlst,
c       collst    i   arrays to hold sections of row and column entries 
c                     of the strict lower triangle
c     itemp1      i   working array of length neqns + 1
c     itemp2      i   working array of length neqns + 1
c     itemp3      i   working array of length neqns + 1
c     itemp4      i   working array of length neqns + 1
c     itemp5      i   working array of length neqns + 1
c     itemp6      i   working array of length neqns + 1
c     itemp7      i   working array of length neqns + 1
c
c     output arguments
c     ----------------
c
c     itemp2      i   holds the pointer to itemp3 for the lists
c                     of column numbers for each compress node. 
c                     (xrowls, length ncomp+1)
c     itemp3      i   holds the lists of column numbers for each
c                     compressed node. 
c                     (rowlst, length neqns)
c     itemp4      i   holds the map from the column numbers 
c                     to the compressed nodes. 
c                     (cmpmap, length neqns)
c     itemp5      i   holds the compressed column starts for the
c                     compressed matrix representation. 
c                     (xadj, length neqns+1)
c     ladjnc      i   holds the adjacency structure for the 
c                     compressed matrix representation. 
c                     (ladjnc, length nzcomp)
c                     note:  this is equivalenced with collst and rowlst
c                            through the argument list.
c     ncomp       i   compressed order
c     nzcomp      i   compressed number of nonzeroes
c     error       i   error flag, = -1 ran out of storage.
c                                 = -2 i/o error
c
c     storage allocation assumptions
c     ------------------------------
c
c     the code in do loop 320 assumes that the array itemp5 can be
c     at most 2*neqns entries long instead of neqns.  this code
c     assumes that itemp6 follows itemp5 and may temporarily write 
c     into itemp6.  note that itemp6 is not active.  
c     The sort and compress will remove the overwrite.  so this
c     code implicitly assumes that itemp5 can be extended.
c
c     ==================================================================
 
c     --------------
c     ... parameters
c     --------------

      integer             neqns , nnzero, maxnz , wafile, output,
     1                    ncomp , nzcomp, error
 
      integer             colptr (*), rowlst (*), collst (*),
     1                    ladjnc (*),
     1                    itemp1 (*), itemp2 (*), itemp3 (*),
     2                    itemp4 (*), itemp5 (*), itemp6 (*),
     3                    itemp7 (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             ck1   , ck2   , ck3   , frw   , i     ,
     1                    icomp , ier   , iold  , ipos  , irow  , 
     2                    j     , jbgn  , jbgno ,
     3                    jcol  , jcolo , jcomp , 
     4                    jend  , jlen  , jrep  , k     , k1    ,
     5                    k1bgn , k1end , klen  , kold  , len   , 
     6                    lstloc, maxint, next  , offset, irep

      logical             qalter, qstore

c.debug
c     double precision    t1, t2, w1, w2
c.debug

c     --------------------
c     ... subprograms used
c     --------------------

      integer             jhmcon

      external            jhmcon
 
c     ==================================================================

c.debug
c     write(6,'("in xisliu - neqns, nnzero, maxnz = ", 3i8)')
c    1                        neqns, nnzero, maxnz 
c     do 1 i = 1, nnzero, neqns
c         len  = min ( nnzero - i + 1, neqns )
c         ipos = i 
c         call xislw1 ( wafile, rowlst, collst,
c    1                  ipos, len, error )
c         call xislp3 ( 'collst', len, collst, 6 )
c         call xislp3 ( 'rowlst', len, rowlst, 6 )
c   1 continue
c.debug
c.debug
c     call xdslt1 ( t1, w1 )
c.debug

c     --------------------------------------------
c     ... build checksum and first-row arrays
c         itemp1 holds first permutation
c         itemp2 holds second permutation
c         itemp3 will hold first  checksum
c         itemp4 will hold second checksum
c         itemp5 will hold third  checksum
c         itemp6 will hold first-row-index
c         itemp7 is initialized for the next phase
c     --------------------------------------------

      maxint =  jhmcon(2) - neqns

      do i = 1, neqns
          itemp3(i) = i
          itemp4(i) = itemp1(i)
          itemp5(i) = itemp2(i)
          itemp6(i) = i
          itemp7(i) = i
      enddo

      offset = 0
      jcolo  = 0

c     ------------------------------------------------------
c     ... read in next buffer full of row and column indices
c     ------------------------------------------------------

   20 continue
      ipos   = offset + 1
      len    = min ( maxnz, nnzero - offset )

      call xislw1 ( wafile, rowlst, collst, ipos, len, error )
      if ( error .ne. 0 ) go to 8100         
c.debug
c     write(6,'("after 20 - len, ipos = ", 2i8)') len, ipos 
c     i = min ( 50, len )
c     j = len - i + 1
c     call xislp3 ( 'column list - first 50', i, collst, 6 )
c     call xislp3 ( 'row    list - first 50', i, rowlst, 6 )
c     call xislp3 ( 'column list - last  50', i, collst(j), 6 )
c     call xislp3 ( 'row    list - last  50', i, rowlst(j), 6 )
c.debug
 
c     ---------------------------------------------------
c     ... gather up the statistics for ths row and column 
c         indices now held in memory.
c     --------------------------------------------------- 
  
      call xislit ( maxint, neqns , len   , offset, jcolo , 
     1              rowlst, collst, colptr, itemp1, itemp2,
     2              itemp3, itemp4, itemp5, itemp6  )

      offset = offset + len
c.debug
c     write(6,'("offset, len, nnzero = ", 3i8)')
c    1            offset, len, nnzero 
c.debug
      if ( offset .lt. nnzero ) go to 20

c     -----------------------------------------------------------
c     ... statistics are gathered for all row and column indices.
c         finish up setting of colptr
c     -----------------------------------------------------------

      do k = jcolo+1, neqns+1
          colptr(k) = nnzero + 1
      enddo
c.debug
c     call xislp3 ( 'after 30 - colptr', neqns+1, colptr, 6 )
c     call xislp3 ( 'after 30 - itemp3', neqns, itemp3, 6 )
c     call xislp3 ( 'after 30 - itemp4', neqns, itemp4, 6 )
c     call xislp3 ( 'after 30 - itemp5', neqns, itemp5, 6 )
c     call xislp3 ( 'after 30 - itemp6', neqns, itemp6, 6 )
c     call xislp3 ( 'after 30 - itemp7', neqns, itemp7, 6 )
c     stop
c     i = 0
c     do 32 k = 1, neqns
c         i = max ( i, itemp3(k), itemp4(k), itemp5(k) )
c  32 continue
c     write(6,'("max. integer in do 30 = ", i40)') i
c.debug
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliu after 31     - accum. time = ", 
c    1         f15.6)') t2
c.debug

c     --------------------------------------------
c     ... note that itemp1 and itemp2 are now free
c     --------------------------------------------
 
c     ==================================================================

c     ---------------------------------------------------
c     ... build compression map
c         start by sorting checksums and first-row arrays 
c         (itemp3, itemp4, itemp5, and itemp6)
c         hjsrtn/hdsrtn perform stable sorts and
c         hjprmx/hdprmx apply the permutation to the 
c         other arrays.  
c         itemp7 holds the row numbers for the data.
c     ---------------------------------------------------

      call hjsrtn ( itemp6, neqns, 0, 0, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp3, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp4, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp5, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp7, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
c.debug
c     call xislp3 ( 'after sort 1 - itemp3', neqns, itemp3, 6 )
c     call xislp3 ( 'after sort 1 - itemp4', neqns, itemp4, 6 )
c     call xislp3 ( 'after sort 1 - itemp5', neqns, itemp5, 6 )
c     call xislp3 ( 'after sort 1 - itemp6', neqns, itemp6, 6 )
c     call xislp3 ( 'after sort 1 - itemp7', neqns, itemp7, 6 )
c.debug
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliu after sort 1 - accum. time = ", 
c    1         f15.6)') t2
c.debug

      call hjsrtn ( itemp5, neqns, 0, 0, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp3, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp4, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp6, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp7, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
c.debug
c     call xislp3 ( 'after sort 2 - itemp3', neqns, itemp3, 6 )
c     call xislp3 ( 'after sort 2 - itemp4', neqns, itemp4, 6 )
c     call xislp3 ( 'after sort 2 - itemp5', neqns, itemp5, 6 )
c     call xislp3 ( 'after sort 2 - itemp6', neqns, itemp6, 6 )
c     call xislp3 ( 'after sort 2 - itemp7', neqns, itemp7, 6 )
c.debug
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliu after sort 2 - accum. time = ", 
c    1         f15.6)') t2
c.debug

      call hjsrtn ( itemp4, neqns, 0, 0, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp3, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp5, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp6, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp7, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
c.debug
c     call xislp3 ( 'after sort 3 - itemp3', neqns, itemp3, 6 )
c     call xislp3 ( 'after sort 3 - itemp4', neqns, itemp4, 6 )
c     call xislp3 ( 'after sort 3 - itemp5', neqns, itemp5, 6 )
c     call xislp3 ( 'after sort 3 - itemp6', neqns, itemp6, 6 )
c     call xislp3 ( 'after sort 3 - itemp7', neqns, itemp7, 6 )
c.debug
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliu after sort 3 - accum. time = ", 
c    1         f15.6)') t2
c.debug

      call hjsrtn ( itemp3, neqns, 0, 0, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp4, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp5, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp6, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
      call hjprmx ( itemp7, neqns, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
c.debug
c     call xislp3 ( 'after sort - itemp3', neqns, itemp3, 6 )
c     call xislp3 ( 'after sort - itemp4', neqns, itemp4, 6 )
c     call xislp3 ( 'after sort - itemp5', neqns, itemp5, 6 )
c     call xislp3 ( 'after sort - itemp6', neqns, itemp6, 6 )
c     call xislp3 ( 'after sort - itemp7', neqns, itemp7, 6 )
c.debug
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliu after sort 4 - accum. time = ", 
c    1         f15.6)') t2
c.debug

c     ----------------------------------------------------------
c     ... now build compress map by finding contiguous sets of
c         similar entries in  itemp3, itemp4, itemp5, and itemp6
c         itemp1 will hold the pointer into itemp7 and itemp7 will
c         hold the list of nodes that have identical structure.
c     ----------------------------------------------------------
      
      ck1 = itemp3(1)
      ck2 = itemp4(1)
      ck3 = itemp5(1)
      frw = itemp6(1)
 
      ncomp         = 1
      itemp1(ncomp) = 1

      do i = 2, neqns

          if ( itemp3(i) .eq. ck1  .and.  itemp4(i) .eq. ck2  .and.
     1         itemp5(i) .eq. ck3  .and.  itemp6(i) .eq. frw        )
     2         cycle

c         --------------------------------------
c         ... now at start of new contiguous set
c         --------------------------------------

          ncomp         = ncomp + 1
          itemp1(ncomp) = i
      
          ck1 = itemp3(i)
          ck2 = itemp4(i)
          ck3 = itemp5(i)
          frw = itemp6(i)

      enddo

      itemp1(ncomp+1) = neqns + 1

c.debug 
c     write(6,'(/"number of compressed nodes = ", i8)') ncomp
c     call xislp3 ( 'xrowls', ncomp+1, itemp1, 6 )
c     call xislp3 ( 'rowlst', neqns, itemp7, 6 )
c.debug 
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliu after 100    - accum. time = ", 
c    1         f15.6)') t2
c.debug

c     --------------------------------------------------------------
c     ... note that arrays itemp2, itemp3, itemp4, itemp5 and itemp6
c         are now free
c     --------------------------------------------------------------

c.debug 
c     ==================================================================

c     ------------------------------------------------------------
c     ... sweep through the sets of contiguous sets and check to
c         see if they are sorted
c     ------------------------------------------------------------

c     do 120 i = 1, ncomp

c         jbgn = itemp1(i)
c         jend = itemp1(i+1) - 1

c         do 119 j = jbgn+1, jend
c             if ( itemp7(j-1) .lt. itemp7(j) ) go to 119
c             write(6,'("itemp7 not sorted")')
c             stop
c 119     continue

c 120 continue
c.debug 
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliu after 120    - accum. time = ", 
c    1         f15.6)') t2
c.debug
 
c     ==================================================================

c     --------------------------------------------------------------
c     ... reorder the continguous sets such that they are ordered so
c         that the first column number in each set is in ascending
c         order
c     --------------------------------------------------------------

      do icomp = 1, ncomp

          jbgn = itemp1(icomp)

          itemp2(icomp) = itemp7(jbgn)

      enddo
c.debug 
c     call xislp3 ( 'itemp2 after 200', ncomp, itemp2, 6 )
c.debug 

c     --------------------------------------------------
c     ... sort itemp2 and save the permutation in itemp4
c     --------------------------------------------------

      call hjsrtn ( itemp2, ncomp, 0, 0, itemp4, ier )
      if ( ier .ne. 0 ) go to 8200
c.debug 
c     call xislp3 ( 'itemp2 after sort', ncomp, itemp2, 6 )
c     call xislp3 ( 'itemp4 after sort', ncomp, itemp4, 6 )
c.debug 
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliu after 200    - accum. time = ", 
c    1         f15.6)') t2
c.debug

c     ---------------------------------------------------------
c     ... use itemp3 to permute the contiguous sets so that the
c         first row number in each set are in ascending order
c         for all of the sets.
c         note that itemp1, itemp7 and itemp4 are in use and
c         the other arrays are free.
c         the new representation will be built in itemp2 and
c         itemp3
c     ---------------------------------------------------------

      next = 1

      do icomp = 1, ncomp
      
          iold = itemp4(icomp)

          jbgn = itemp1(iold)
          jlen = itemp1(iold+1) - jbgn
c.debug
c     write(6,'("in 210 - icomp, iold, next, jbgn, jlen = ", 5i8)')
c    1                     icomp, iold, next, jbgn, jlen 
c     call xislp3 ( 'itemp7(jbgn)', jlen, itemp7(jbgn), 6 )
c.debug

          itemp2(icomp) = next

          itemp3(next:next+jlen-1) = itemp7(jbgn:jbgn+jlen-1)

          next = next + jlen

      enddo

      itemp2(ncomp+1) = next
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliu after 210    - accum. time = ", 
c    1         f15.6)') t2
c.debug

c     -----------------------------------------------------------
c     ... build a map of the original nodes to the compress nodes
c         in itemp4. 
c         note that arrays itemp4 through itemp7 are now free
c     -----------------------------------------------------------

      do icomp = 1, ncomp

          jbgn = itemp2(icomp)
          jend = itemp2(icomp+1) - 1
c.debug
c     write(6,'("in 210 - icomp, jbgn, jend = ", 3i8)')
c    1                     icomp, jbgn, jend 
c     call xislp3 ( 'itemp3(jbgn)', jlen, itemp3(jbgn), 6 )
c.debug

          do j = jbgn, jend

              k         = itemp3(j)
              itemp4(k) = icomp

          enddo

      enddo
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliu after 230    - accum. time = ", 
c    1         f15.6)') t2
c.debug

c     --------------------------------------------------------------
c     ... note that itemp2 holds the pointer to itemp3 for the lists
c                          of column numbers for each compress node.
c                   itemp3 holds the lists of column numbers for each
c                          compressed node
c                   itemp4 holds the map from the column numbers 
c                          to the compressed nodes.
c                   itemp1 and itemp 5 through itemp7 are free
c     --------------------------------------------------------------

c.debug 
c     write(6,'("after 230")')
c     call xislp3 ( 'itemp2 - xrowls',              ncomp+1, itemp2, 6 )
c     call xislp3 ( 'itemp3 - rowlst',              neqns  , itemp3, 6 )
c     call xislp3 ( 'itemp4 - map from col. nos. ', neqns  , itemp4, 6 )
c.debug 

c     ==================================================================

c     ------------------------------------------------------
c     ... now build the compressed matrix representation in
c         arrays xadj (in itemp5) and ladjnc (over the top 
c         of collst and rowlst)
c     ------------------------------------------------------

      next   = 1
      lstloc = 0
      offset = 0
      qstore = .true.
      jcolo  = 0

      do 330 icomp = 1, ncomp

          jcol = itemp3(itemp2(icomp))

          k1bgn = colptr(jcol)
          k1end = colptr(jcol+1) - 1

          do k = jcolo+1, jcol
              itemp5(k) = next
          enddo

          jcolo = jcol

          if ( k1bgn .lt. offset+1 .or. k1end .gt. lstloc ) then

c             ------------------------------------------------------
c             ... read in the row and column indices for column jcol
c             ------------------------------------------------------

              ipos = k1bgn 

              len  = min ( neqns, nnzero - k1bgn + 1 )
c.debug
c     write(6,'("in 330 - reading - wafile, ipos, len = ", 3i8)')
c    1                               wafile, ipos, len 
c.debug

              call xislw1 ( wafile, itemp1, itemp6, ipos, len, error )
              if ( error .ne. 0 ) go to 8100

              offset = k1bgn - 1
              lstloc = offset + len

          end if

c         --------------------------------------------------
c         ... place the compress row numbers for the entries 
c             from rowlst to itemp6
c         --------------------------------------------------

c.debug 
c     write(6,'("first column for compressed node")')
c     write(6,'("jcol, k1bgn, k1end, offset, next = ", 5i8)')
c    1            jcol, k1bgn, k1end, offset, next     
c     call xislp3 ( 'rowlst', k1end-k1bgn+1, itemp1(k1bgn-offset), 6 )
c.debug 

          klen = 0

          do k = k1bgn, k1end 
              klen         = klen + 1
              jcomp        = itemp4 ( itemp1(k-offset) )
              itemp6(klen) = itemp3 ( itemp2(jcomp) ) 
          enddo

c         ----------------------------
c         ... sort and compress itemp6
c         ----------------------------

c.debug 
c     call xislp3 ( 'itemp6 before sort', klen, itemp6, 6 )
c     if ( klen .lt. 0 ) write('("oops")')
c.debug 
          if ( klen .gt. 0 ) then
              call xislq1 ( klen, itemp6, ier )
              if ( ier .ne. 0 ) go to 8200
          end if
c.debug 
c     call xislp3 ( 'itemp6 after  sort', klen, itemp6, 6 )
c.debug 

          k1end = klen
          klen  = 0
          kold  = jcol

          do k1 = 1, k1end

              if ( itemp6(k1) .gt. kold  ) then

                klen         = klen + 1
                kold         = itemp6(k1)
                itemp6(klen) = kold      

              endif

          enddo
c.debug 
c     call xislp3 ( 'itemp6 after compression', klen, itemp6, 6 )
c.debug 

c         ---------------------------------------------------
c         ... store compressed row indicies in final location
c         ---------------------------------------------------

          if ( next + klen - 1 .gt. 2 * maxnz ) qstore = .false.

          if ( qstore ) then
              ladjnc(next:next+klen-1) = itemp6(1:klen)
          end if

          next = next + klen

c         ---------------------------------------------------
c         ... loop back to next column of the original matrix
c         ---------------------------------------------------

  330 continue

      nzcomp = next - 1
 
      do k = jcolo+1, neqns+1
          itemp5(k) = next
      enddo
 
c.debug 
c     write(6,'(/"number of compressed nonzeroes = ", i10)') nzcomp
c     call xislp3 ( 'compressed colstr', neqns+1, itemp5, 6 )
c.debug 
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliu after 340    - accum. time = ", 
c    1         f15.6)') t2
c.debug

      if ( .not. qstore ) go to 8000

c.debug 
c     call xislp3 ( 'compressed ladjnc', nzcomp , ladjnc, 6 )
c.debug 

c     ==================================================================

c     -----------------------------------------------------------
c     ... verify that complete matrix structure is represented in
c         the compressed adjncy form
c     -----------------------------------------------------------

      qalter = .false.

      offset = 0
      lstloc = 0

      do 440 jcol = 1, neqns

          k1bgn = colptr(jcol)
          k1end = colptr(jcol+1) - 1

          jcomp = itemp4(jcol)
          jrep  = itemp3(itemp2(jcomp))

          if ( k1bgn .le. offset .or. k1end .gt. lstloc ) then

c             ------------------------------------------------------
c             ... read in next buffer full of row and column indices
c             ------------------------------------------------------

              ipos = k1bgn

              len  = min ( neqns, nnzero - k1bgn + 1 )
              call xislw1 ( wafile, itemp6, itemp7, ipos, len, error )
              if ( error .ne. 0 ) go to 8100

              offset = k1bgn - 1
              lstloc = offset + len

          end if

          klen  = 0

          do k = k1bgn, k1end

              irow  = itemp6(k-offset)

              icomp = itemp4(irow)
              irep  = itemp3(itemp2(icomp))
c.debug
c         write(6,'("jcol, jcomp, jrep, k1bgn, k1end = ", 5i8)')
c    1                jcol, jcomp, jrep, k1bgn, k1end 
c         write(6,'("k, irow, icomp, irep            = ", 5i8)')
c    1                k, irow, icomp, irep            
c.debug

              if ( irep .gt. jrep ) then

                klen  = klen + 1
                itemp1(klen) = irep
  
              endif

          enddo

          if ( klen .eq. 0 ) go to 440

c         ----------------------------
c         ... sort and compress itemp1
c         ----------------------------

c.debug 
c     call xislp3 ( 'itemp1 before sort', klen, itemp1, 6 )
c.debug 
          call xislq1 ( klen, itemp1, ier )
          if ( ier .ne. 0 ) go to 8200
c.debug 
c     call xislp3 ( 'itemp1 after  sort', klen, itemp1, 6 )
c.debug 

          k1end = klen
          klen  = 0
          kold  = 0

          do k1 = 1, k1end

              if ( itemp1(k1) .ne. kold ) then

                klen         = klen + 1
                kold         = itemp1(k1)
                itemp1(klen) = kold      

              endif

          enddo
c.debug 
c     call xislp3 ( 'itemp1 after compression', klen, itemp1, 6 )
c.debug 

c         --------------------------------------------------
c         ... verify that (irep, jrep) is in compress adjncy
c         --------------------------------------------------

          jbgn  = itemp5(jrep)
          jend  = itemp5(jrep+1) - 1

          do 430 k1 = 1, klen

              irow = itemp1(k1)
 
              jbgno = jbgn
c.debug 
c     write(6,'("before 420")')
c     write(6,'("jbgn, jend, k1, irow = ",4i8)')
c    1            jbgn, jend, k1, irow
c     call xislp3 ( 'section of collst', jend-jbgn+1, collst(jbgn), 6 )
c.debug 
              do j = jbgno, jend
c.debug
c     write(6,'("in 420 - j, collst(j) = ", 2i8)')
c    1                     j, collst(j)
c.debug
                  jbgn = j
                  if ( irow .eq. collst(j) ) go to 430
                  if ( irow .lt. collst(j) ) go to 425
              enddo

c             ---------------------------------
c             ... oops, not found, so insert it
c             ---------------------------------

  425         continue
              qalter = .true.

c.debug
c     write(6,'("compression oops - having to insert entry")')
c     write(6,'("jrep, irow = ", 2i8)') jrep, irow
c     stop
c.debug

              if ( nzcomp + 1 .gt. maxnz ) qstore = .false.

              if ( qstore ) then
                  call xislix ( jrep, irow, jbgn, neqns, nzcomp, 
     1                          itemp5, collst )
              end if

  430     continue

  440 continue

      if ( .not. qalter ) go to 9000

c.debug 
      write(6,'("altered number of compressed nonzeroes = ", i10)') 
     1            nzcomp
      call xislp3 ( 'compressed xadj  ', neqns+1, itemp5, 6 )
      call xislp3 ( 'compressed adjncy', nzcomp , collst, 6 )
c.debug 
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliu after 440    - accum. time = ", 
c    1         f15.6)') t2
c.debug

c     ==================================================================
   
   
      go to 9000

c     ==================================================================

c     ----------------------------------------------
c     ... ran out of memory building xadj and ladjnc
c     ----------------------------------------------

 8000 continue
      error = -1
      go to 9000

c     -------------
c     ... i/o error
c     -------------

 8100 continue
      error = -2
      go to 9000

c     -----------------
c     ... sorting error
c     -----------------

 8200 continue
      error = -3
      go to 9000

c     ==================================================================

 9000 continue
      return
      end
