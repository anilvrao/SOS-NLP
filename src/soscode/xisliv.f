      subroutine   xisliv   ( neqns , nnzero, maxnz , 
     1                        colptr, rowlst, collst, 
     2                        itemp1, itemp2, itemp3, itemp4, itemp5,
     3                        itemp6, itemp7, ncomp , nzcomp  )
c
c
c     ==================================================================
c     ====  xisliv -- builds compressed adjacency structure from    ====
c     ====            rowlst and collst structure in memory         ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xisliv builds a compressed adjacency structure from the in-memory
c     representation of rowlst and collst.
c
c     created         04-dec-96   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     neqns       l   order of uncompressed problem.
c     nnzero      i   the number of entries in rowlst and collst
c     maxnz       i   size of rowlst and collst arrays.
c     rowlst,
c       collst    i   lists of row and column entries of the org. matrix
c     itemp1      i   array of length neqns holding first random
c                     permutation for use in computing second checksum
c     itemp2      i   array of length neqns holding second random
c                     permutation for use in computing first checksum
c
c     working vectors 
c     ---------------
c
c     colptr      i   pointer to the start of each column in 
c                     rowlst/collst
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
c     itemp5      i   holds the compressed xadj for the compressed
c                     matrix representation. 
c                     (xadj, length neqns+1)
c     collst      i   holds the adjacency structure for the 
c                     compressed matrix representation. 
c                     (adjncy, length nzcomp)
c     ncomp       i   compressed order
c     nzcomp      i   compressed number of nonzeroes
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

      integer             neqns , nnzero, maxnz , ncomp , nzcomp
 
      integer             colptr (*), rowlst (*), collst (*),
     1                    itemp1 (*), itemp2 (*), itemp3 (*),
     2                    itemp4 (*), itemp5 (*), itemp6 (*),
     3                    itemp7 (*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             ck1   , ck2   , ck3   , frw   , i     ,
     1                    ier   , iold  , irow  , 
     2                    j     , jbgn  , jcol  , jcol1 , jcol2 ,
     3                    jcolo , jcomp , jend  , jlen  , k     , 
     4                    k1    , k1bgn , k1end , k2    , k2bgn ,
     5                    k2end , klen  , kold  , maxint, next  ,
     6                    error , icomp , irep  , jbgno , jrep  ,
     7                    offset

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
c     write(6,'("in xisliv - neqns, nnzero, maxnz = ", 3i8)')
c    1                        neqns, nnzero, maxnz 
c     call xislp3 ( 'collst', nnzero, collst, 6 )
c     call xislp3 ( 'rowlst', nnzero, rowlst, 6 )
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

      jcolo  = 0
      offset = 0

c     ---------------------------------------------------
c     ... gather up the statistics for ths row and column
c         indices
c     ---------------------------------------------------
 
      call xislit ( maxint, neqns , nnzero, offset, jcolo ,
     1              rowlst, collst, colptr, itemp1, itemp2,
     2              itemp3, itemp4, itemp5, itemp6  )

c     -----------------------------------------------------------
c     ... statistics are gathered for all row and column indices.
c         finish up setting of colptr
c     -----------------------------------------------------------

      do k = jcolo+1, neqns+1
          colptr(k) = nnzero + 1
      enddo
c.debug
c     call xislp3 ( 'after 30 - colptr', neqns, colptr, 6 )
c     call xislp3 ( 'after 30 - itemp3', neqns, itemp3, 6 )
c     call xislp3 ( 'after 30 - itemp4', neqns, itemp4, 6 )
c     call xislp3 ( 'after 30 - itemp5', neqns, itemp5, 6 )
c     call xislp3 ( 'after 30 - itemp6', neqns, itemp6, 6 )
c     call xislp3 ( 'after 30 - itemp7', neqns, itemp7, 6 )
c     i = 0
c     do 32 k = 1, neqns
c         i = max ( i, itemp3(k), itemp4(k), itemp5(k) )
c  32 continue
c     write(6,'("max. integer in do 30 = ", i40)') i
c.debug
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliv after 31     - accum. time = ", 
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
c     write(6,'(15x,"in xisliv after sort 1 - accum. time = ", 
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
c     write(6,'(15x,"in xisliv after sort 2 - accum. time = ", 
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
c     write(6,'(15x,"in xisliv after sort 3 - accum. time = ", 
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
c     write(6,'(15x,"in xisliv after sort 4 - accum. time = ", 
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
c     write(6,'("number of compressed nodes = ", i8)') ncomp
c     call xislp3 ( 'pointer into itemp7', ncomp+1, itemp1, 6 )
c     call xislp3 ( 'itemp7', neqns, itemp7, 6 )
c.debug 
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliv after 100    - accum. time = ", 
c    1         f15.6)') t2
c.debug

c     --------------------------------------------------------------
c     ... note that arrays itemp2, itemp3, itemp4, itemp5 and itemp6
c         are now free
c     --------------------------------------------------------------

c     ------------------------------------------------------------
c     ... sweep through the sets of contiguous sets and sort them.
c     ------------------------------------------------------------

c     do 110 i = 1, ncomp

c         jbgn = itemp1(i)

c         if ( i .lt. ncomp ) then
c             jlen = itemp1(i+1) - jbgn
c         else
c             jlen = neqns - jbgn + 1
c         end if

c         call xislq1 ( jlen, itemp7(jbgn), ier )
c             if ( ier .ne. 0 ) go to 8200

c 110 continue

c.debug 
c     call xislp3 ( 'itemp7 after 110', neqns, itemp7, 6 )
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
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliv after 120    - accum. time = ", 
c    1         f15.6)') t2
c.debug

c     ----------------------------
c     ... validate the compression
c     ----------------------------

      do 123 i = 1, ncomp

          jbgn = itemp1(i)
          jend = itemp1(i+1) - 1

          jcol1 = itemp7(jbgn)

          k1bgn = colptr(jcol1)
          k1end = colptr(jcol1+1) - 1

          do j = jbgn+1, jend

              jcol2 = itemp7(j)
              k2bgn = colptr(jcol2)
              k2end = colptr(jcol2+1) - 1

c.debug
c     write(6,'("before 121 loop")')
c     write(6,'("i, jbgn, jend           = ", 5i8)')
c    1            i, jbgn, jend               
c     write(6,'("jcol1, k1bgn, k1end     = ", 5i8)')
c    1            jcol1, k1bgn, k1end
c     write(6,'("jcol2, k2bgn, k2end     = ", 5i8)')
c    1            jcol2, k2bgn, k2end
c     k1 = k1end - k1bgn + 1
c     write(6,'("k1                      = ", 5i8)')
c    1            k1                          
c     call xislp3 ( 'rowlst for 1st col', k1, rowlst(k1bgn), 6 )
c     call xislp3 ( 'rowlst for 2nd col', k1, rowlst(k2bgn), 6 )
c.debug

              k1 = k1end+1
              do k2 = k2end, k2bgn, -1

                  k1 = k1 - 1

                  if ( rowlst(k1) .eq. rowlst(k2) ) cycle

                  write(6,'("compression failure 12")')
                  write(6,'("i, jbgn, jend               = ", 5i8)')
     1                        i, jbgn, jend               
                  write(6,'("jcol1, k1bgn, k1end, k1     = ", 5i8)')
     1                        jcol1, k1bgn, k1end, k1
                  write(6,'("jcol2, k2bgn, k2end, k2     = ", 5i8)')
     1                        jcol2, k2bgn, k2end, k2
                  write(6,'("rowlst(k1), rowlst(k2)      = ", 5i8)')
     1                        rowlst(k1), rowlst(k2)
                  write(6,'("collst(k1), collst(k2)      = ", 5i8)')
     1                        collst(k1), collst(k2)

              enddo

          enddo

  123 continue

c     write(6,'("end of validation loop 123")')
c.debug 
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliv after 123    - accum. time = ", 
c    1         f15.6)') t2
c.debug
 
c     ==================================================================

c     --------------------------------------------------------------
c     ... reorder the continguous sets such that they are ordered so
c         that the first column number in each set is in ascending
c         order
c     --------------------------------------------------------------

      do i = 1, ncomp

          jbgn = itemp1(i)

          itemp2(i) = itemp7(jbgn)

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
c     write(6,'(15x,"in xisliv after 200    - accum. time = ", 
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

      do i = 1, ncomp
      
          iold = itemp4(i)

          jbgn = itemp1(iold)
          jlen = itemp1(iold+1) - jbgn
c.debug
c     write(6,'("in 210 - i, iold, next, jbgn, jlen = ", 5i8)')
c    1                     i, iold, next, jbgn, jlen 
c     call xislp3 ( 'itemp7(jbgn)', jlen, itemp7(jbgn), 6 )
c.debug

          itemp2(i) = next

          itemp3(next:next+jlen-1) = itemp7(jbgn:jbgn+jlen-1)

          next = next + jlen

      enddo

      itemp2(ncomp+1) = next
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliv after 210    - accum. time = ", 
c    1         f15.6)') t2
c.debug

c     -----------------------------------------------------------
c     ... build a map of the original nodes to the compress nodes
c         in itemp4. 
c         note that arrays itemp4 through itemp7 are now free
c     -----------------------------------------------------------

      do i = 1, ncomp

          jbgn = itemp2(i)
          jend = itemp2(i+1) - 1
c.debug
c     write(6,'("in 210 - i, jbgn, jend = ", 3i8)')
c    1                     i, jbgn, jend 
c     call xislp3 ( 'itemp3(jbgn)', jlen, itemp3(jbgn), 6 )
c.debug

          do j = jbgn, jend

              k         = itemp3(j)
              itemp4(k) = i

          enddo

      enddo
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliv after 230    - accum. time = ", 
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

c.debug
c     ==================================================================

c     ----------------------------
c     ... validate the compression
c     ----------------------------

      do 223 i = 1, ncomp

          jbgn = itemp2(i)
          jend = itemp2(i+1) - 1

          jcol1 = itemp3(jbgn)

          if ( i .eq. 1 ) then
              jcolo = 1
              if ( jcol1 .ne. 1 ) then
                  write(6,'("comp. failure 20 - jcol1, jcolo = ",
     1                    2i8)') jcol1, jcolo
              end if
          else
              if ( jcol1 .le. jcolo ) then
                  write(6,'("comp. failure 20 - jcol1, jcolo = ",
     1                    2i8)') jcol1, jcolo
              end if
              jcolo = jcol1
          end if

          k1bgn = colptr(jcol1)
          k1end = colptr(jcol1+1) - 1

          do j = jbgn+1, jend

              jcol2 = itemp3(j)
              k2bgn = colptr(jcol2)
              k2end = colptr(jcol2+1) - 1

c.debug
c     write(6,'("before 221 loop")')
c     write(6,'("i, jbgn, jend           = ", 5i8)')
c    1            i, jbgn, jend               
c     write(6,'("jcol1, k1bgn, k1end     = ", 5i8)')
c    1            jcol1, k1bgn, k1end
c     write(6,'("jcol2, k2bgn, k2end     = ", 5i8)')
c    1            jcol2, k2bgn, k2end
c     k1 = k1end - k1bgn + 1
c     write(6,'("k1                      = ", 5i8)')
c    1            k1                          
c     call xislp3 ( 'rowlst for 1st col', k1, rowlst(k1bgn), 6 )
c     call xislp3 ( 'rowlst for 2nd col', k1, rowlst(k2bgn), 6 )
c.debug

              k1 = k1end+1
              do k2 = k2end, k2bgn, - 1

                  k1 = k1 - 1

                  if ( rowlst(k1) .eq. rowlst(k2) ) cycle

                  write(6,'("compression failure 22")')
                  write(6,'("i, jbgn, jend               = ", 5i8)')
     1                        i, jbgn, jend               
                  write(6,'("jcol1, k1bgn, k1end, k1     = ", 5i8)')
     1                        jcol1, k1bgn, k1end, k1
                  write(6,'("jcol2, k2bgn, k2end, k2     = ", 5i8)')
     1                        jcol2, k2bgn, k2end, k2
                  write(6,'("rowlst(k1), rowlst(k2)      = ", 5i8)')
     1                        rowlst(k1), rowlst(k2)
                  write(6,'("collst(k1), collst(k2)      = ", 5i8)')
     1                        collst(k1), collst(k2)

              enddo

          enddo

  223 continue

c     write(6,'("end of validation loop 223")')
c.debug 
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliv after 223    - accum. time = ", 
c    1         f15.6)') t2
c.debug
 
c     ==================================================================

c     ------------------------------------------------------
c     ... now build the compressed matrix representation in
c         arrays xadj (in itemp5) and adjncy (over the top 
c         of collst)
c     ------------------------------------------------------

      next   = 1
      jcolo  = 0
      qstore = .true.

      do 330 i = 1, ncomp

          jbgn = itemp2(i)
          jend = itemp2(i+1) - 1

          jcol = itemp3(jbgn)

          do k = jcolo+1, jcol
              itemp5(k) = next
          enddo

          jcolo = jcol

          klen = 0
c.debug 
c     write(6,'("in 330 - i, next, jbgn, jend = ", 4i8)')
c    1                     i, next, jbgn, jend 
c     call xislp3 ( 'itemp3', jend - jbgn + 1, itemp3(jbgn), 6 )
c.debug 

          k1bgn = colptr(jcol)
          k1end = colptr(jcol+1) - 1

c.debug 
c     write(6,'("in 320 - j, jcol, k1bgn, k1end = ", 4i8)')
c    1                     j, jcol, k1bgn, k1end     
c.debug 

c         ------------------------------------------------------
c         ... add the representative row numbers for the entries 
c             from rowlst to itemp6
c         ------------------------------------------------------

          do k = k1bgn, k1end
              klen         = klen + 1
c.debug
c     write(6,'("in 305 - k, klen, rowlst(k), itemp4(rowlst(k)) = ",
c    1            4i8)')   k, klen, rowlst(k), itemp4(rowlst(k)) 
c.debug
              jcomp        = itemp4( rowlst(k) )
              itemp6(klen) = itemp3(itemp2(jcomp))
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

              if ( itemp6(k1) .gt. kold ) then

                klen         = klen + 1
                kold         = itemp6(k1)
                itemp6(klen) = kold      

              endif
          enddo
c.debug 
c     call xislp3 ( 'itemp6 after compression', klen, itemp6, 6 )
c.debug 

  320     continue
c.debug 
c     call xislp3 ( 'itemp6 after 320 loop', klen, itemp6, 6 )
c.debug 

c         ---------------------------------------------------
c         ... itemp6 now contains the list of columns for the 
c             compressed structure.
c         ---------------------------------------------------

          if ( next + klen - 1 .gt. maxnz ) qstore = .false.

          if ( qstore ) then
              collst(next:next+klen-1) = itemp6(1:klen)
          end if

          next = next + klen

  330 continue

      nzcomp = next - 1

      do k = jcolo+1, neqns+1
          itemp5(k) = next
      enddo

c.debug 
c     write(6,'("number of compressed nonzeroes = ", i10)') nzcomp
c     call xislp3 ( 'compressed xadj  ', neqns+1, itemp5, 6 )
c.debug 
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliv after 340    - accum. time = ", 
c    1         f15.6)') t2
c.debug

      if ( .not. qstore ) go to 8000

c.debug 
c     call xislp3 ( 'compressed adjncy', nzcomp , collst, 6 )
c.debug 

c     ==================================================================

c     -----------------------------------------------------------
c     ... verify that complete matrix structure is represented in
c         the compressed adjncy form
c     -----------------------------------------------------------

      qalter = .false.

      do 440 jcol = 1, neqns

          k1bgn = colptr(jcol)
          k1end = colptr(jcol+1) - 1

          jcomp = itemp4(jcol)
          jrep  = itemp3(itemp2(jcomp))

          klen  = 0

          do k = k1bgn, k1end

              irow  = rowlst(k)

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
c.debug
c     call xdslt2 ( t1, w1, t2, w2 )
c     write(6,'(15x,"in xisliv after 440    - accum. time = ", 
c    1         f15.6)') t2
c.debug

      if ( .not. qalter ) go to 9000

c.debug 
      write(6,'("altered number of compressed nonzeroes = ", i10)') 
     1            nzcomp
      call xislp3 ( 'compressed xadj  ', neqns+1, itemp5, 6 )
      call xislp3 ( 'compressed adjncy', nzcomp , collst, 6 )
c.debug 

c     ==================================================================
   
      go to 9000

c     ==================================================================

c     ----------------------------------------------
c     ... ran out of memory building xadj and adjncy
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
