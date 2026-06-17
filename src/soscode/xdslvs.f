      subroutine xdslvs ( neqns , nsuper, xsup  , xarndx, nrecrd,
     1                    wafil1, recps1, recln1, walen1, watrn1, 
     2                    wafil2, recps2, recln2, walen2, watrn2, 
     3                    lset  , coor11, coor12, value1, coor21,
     4                    coor22, value2, coor31, coor32, value3,
     5                    nnzero, maxnnz, error )
c
c     purpose
c     -------
c
c     xdslvs performs a sort-merge of the coordinates and value1 for
c     the out-of-core matrix value input phase.

c     wafil1 on input holds nrecrd sorted records.  record irec starts
c     location recps1(irec) and contains recln1(irec) row and column 
c     indices.  each record is sorted first by coordinate 1 and then 
c     by coordinate 2.  on output the file then labeled wafil2 will hold 
c     all of coordinate and value data in one sorted record without 
c     duplicates.
c
c     created         10-jan-97   -- rgg --
c     modified        
c
c     input arguments
c     ---------------
c
c     neqns       i   number of rows in org. matrix
c     nsuper      i   number of super nodes
c     xsup        i   super node partition array
c     nrecrd      i   number of records originally on wafil1
c     lset        i   number of indicies that can be stored in 
c                     memory at any one time.
c
c     input/output arguments
c     ----------------------
c
c     wafil1      i   unit number for first file.
c     wafil2      i   unit number for second file. 
c                     (Note:  these will toggle back and forth at
c                      each step)
c     recps1      i   array of length nrecrd holding i/o positions
c                     at the start of each record for wafil1
c     recln1      i   array of length nrecrd holding length of
c                     of each record for wafil1
c     recps2      i   array of length nrecrd holding i/o positions
c                     at the start of each record for wafil1
c     recln2      i   array of length nrecrd holding length of
c                     of each record for wafil1
c
c     working storage
c     ---------------
c
c     coor11      i   array of length lset to hold coordinate 1 
c                     indicies from set 1 in memory
c     coor12      i   array of length lset to hold coordinate 2 
c                     indicies from set 1 in memory
c     value1      i   array of length lset to hold values from
c                     set 1 in memory
c     coor21      i   array of length lset to hold coordinate 1 
c                     indicies from set 2 in memory
c     coor22      i   array of length lset to hold coordinate 2 
c                     indicies from set 2 in memory
c     value2      i   array of length lset to hold values from
c                     set 2 in memory
c     coor31      i   array of length 2*lset to hold coordinate 1 
c                     indicies from set 3 in memory
c     coor32      i   array of length 2*lset to hold coordinate 2 
c                     indicies from set 3 in memory
c     value3      i   array of length 2*lset to hold values from
c                     set 3 in memory
c
c     output arguments
c     ----------------
c
c     xardnx      i   pointer array for relative assembly indices
c     walen1      i   length of i/o file wafil1
c     watrn1      i   amount of i/o transfer to and from wafil1
c     walen2      i   length of i/o file wafil2
c     watrn2      i   amount of i/o transfer to and from wafil2
c     nnzero      i   number of nonzeroes in final form
c     maxnnz      i   max. number of nonzeroes for a super node
c     error       i   error flag
c                     =  0  normal return
c                     = -1  i/o error on wafil1 or wafil2
 
c     ------------------------------------------------------------------

c     -------------------------------
c     ... global variable declaration
c     -------------------------------

      integer             wafil1, nrecrd, wafil2, lset  , 
     1                    walen1, watrn1, walen2, watrn2,
     2                    neqns , nsuper, nnzero, maxnnz, error

      integer             recps1 (*),     recln1 (*),
     1                    recps2 (*),     recln2 (*),
     2                    coor11 (*),     coor12 (*),
     2                    coor21 (*),     coor22 (*),
     2                    coor31 (*),     coor32 (*),
     4                    xsup   (*),     xarndx (*)

      double precision    value1 (*),     value2 (*),
     1                    value3 (*)

c     ------------------------------
c     ... local variable declaration
c     ------------------------------

      integer             aclen1, aclen2, itemp , k11len, k11pos, 
     1                    k12len, k12pos, k2pos , len1  , len2  ,
     2                    nrcrd1, nrcrd2, orgwa1, rec1  , rec2  ,
     3                    isuper, nnz
c.debug
c     integer             i     , ipos  , len
c.debug

      logical             qfinal

      integer             xdslni
 
      external            xdslni

c     ------------------------------------------------------------------

c     ------------------
c     ... initialization
c     ------------------

      error  = 0
      nrcrd1 = nrecrd
      nnzero = 0

c     -------------------------------------------------------
c     ... record i/o statistics about current state of wafil1
c     -------------------------------------------------------

      orgwa1 = wafil1
 
      len1   = recps1(nrecrd) + recln1(nrecrd) - 1 
      aclen1 = xdslni ( 2*len1 ) + len1
      walen1 = aclen1
      watrn1 = walen1
 
      len2   = 0
      aclen2 = 0
      walen2 = 0
      watrn2 = 0

c     ------------------------------------------------------------------

  100 continue
      nrcrd2 = 0
      k2pos  = 1

c.debug
c     write(6,'(/,"entering 100 loop in xdslvs")') 
c     write(6,'("wafil1, wafil2 = ", 2i8)') wafil1, wafil2
c
c     call xislp3 ( 'record pos.', nrcrd1, recps1, 6 )
c     call xislp3 ( 'record len.', nrcrd1, recln1, 6 )

c     nnzero = 0
c     do 2 i = 1, nrcrd1
c         nnzero = nnzero + recln1(i)
c   2 continue

c     do 1 i = 1, nnzero, 50
c         len  = min ( nnzero - i + 1, 50 )
c         ipos = i 
c         call xdslw1 ( wafil1, 3, coor11, coor12, value1,
c    1                  ipos, nnzero, error )
c         call xislp3 ( 'coor11', len, coor11, 6 )
c         call xislp3 ( 'coor12', len, coor12, 6 )
c         call xdslp5 ( 'value1', len, value1, 6 )
c   1 continue
c.debug

      if ( nrcrd1 .le. 2 ) then
          qfinal = .true.
      else
          qfinal = .false.
      end if

      do rec1 = 1, nrcrd1, 2

          rec2           = rec1 + 1
          k11pos         = recps1(rec1)
          k11len         = recln1(rec1)

          if ( rec2 .le. nrcrd1 ) then
              k12pos         = recps1(rec2)
              k12len         = recln1(rec2)
          else
              k12len         = 0
          end if

          nrcrd2         = nrcrd2 + 1
          recps2(nrcrd2) = k2pos
c.debug
c     write(6,'("in 200 - rec1  , rec2  , k2pos          = ", 4i8)')
c    1                     rec1  , rec2  , k2pos          
c     write(6,'("         k11pos, k11len, k12pos, k12len = ", 4i8)')
c    1                     k11pos, k11len, k12pos, k12len
c     write(6,'("         qfinal                         = ",  l8)')
c    1                     qfinal                        
c.debug

          call xdslvt ( wafil1, k11pos, k11len, k12pos, k12len,
     1                  wafil2, k2pos , lset  , coor11, coor12, 
     2                  value1, coor21, coor22, value2, coor31,
     3                  coor32, value3,
     4                  qfinal, neqns , nnzero, xarndx, error )
c.debug
c     write(6,'("after xdslvt - nrcrd2, k2pos, error = ", 3i8)')
c    1                           nrcrd2, k2pos, error 
c.debug

          recln2(nrcrd2) = k2pos - recps2(nrcrd2) 

          if ( error .ne. 0 ) return

      enddo

c     -------------------------
c     ... record i/o statistics
c     -------------------------

      if ( wafil1 .eq. orgwa1 ) then

          len2   = k2pos - 1
          aclen2 = xdslni ( 2*len2 ) + len2
      
          watrn1 = watrn1 + aclen1

          walen2 = max ( walen2, aclen2 )
          watrn2 = watrn2 + aclen2

      else

          len1   = k2pos - 1
          aclen1 = xdslni ( 2*len1 ) + len1
      
          walen1 = max ( walen1, aclen1 )
          watrn1 = watrn1 + aclen1 

          watrn2 = watrn2 + aclen2

      end if

c     ------------------------------------------------------------------

c     --------------------------
c     ... prepare for next sweep
c     --------------------------

      if ( nrcrd2 .gt. 1 ) then

c         ---------------------------------------
c         ... copy the info from wafil2 to wafil1
c         ---------------------------------------

          itemp  = wafil1
          wafil1 = wafil2
          wafil2 = itemp

          recps1(1:nrcrd2) = recps2(1:nrcrd2)
          recln1(1:nrcrd2) = recln2(1:nrcrd2)

          nrcrd1 = nrcrd2

c         -------------------------------------------
c         ... loop back to 100 and make another sweep
c         -------------------------------------------

          go to 100

      end if

c     ------------------------------------------------------------------

c     --------------------------------------------------------------
c     ... all the data has been sorted and redundant entries removed
c         and left on file currently labeled as wafil2.
c     --------------------------------------------------------------

      maxnnz = 0

      do isuper = 1, nsuper

          nnz    = xarndx(xsup(isuper+1)) - xarndx(xsup(isuper))

          maxnnz = max ( maxnnz, nnz )

      enddo

c.debug
c     call xislp3 ( 'xarndx at end of xdslvs', neqns+1, xarndx, 6 )
c     write(6,'("i/o file at end of xdslvs")')
c     write(6,'("maxnnz, nnzero, wafil2 = ", 3i8)') 
c    1            maxnnz, nnzero, wafil2
c     do 301 i = 1, nnzero, 50
c         len  = min ( nnzero - i + 1, 50 )
c         ipos = i 
c
c         write(6,'(/, "from ", i8, " to ", i8)') i, i+len-1
c
c         call xdslw1 ( wafil2, 3, coor11, coor12, value1,
c    1                  ipos, len, error )
c
c         write(6,'("after xdslw1 - error = ", i8)') error
c
c         call xislp3 ( 'coor11', len, coor11, 6 )
c         call xislp3 ( 'coor12', len, coor12, 6 )
c         call xdslp5 ( 'value1', len, value1, 6 )
c 301 continue
c.debug

c     ------------------------------------------------------------------

      return
      end
