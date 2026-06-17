      subroutine xislis ( nrecrd, sort  , 
     1                    wafil1, recps1, recln1, walen1, watrn1, 
     2                    wafil2, recps2, recln2, walen2, watrn2, 
     3                    lset  , rowls1, colls1, rowls2, colls2,
     4                    rowls3, colls3, nnzero, error )

c
c     purpose
c     -------
c
c     xislis performs a sort-merge of the row and column indicies for
c     the out-of-core matrix structure input phase.

c     wafil1 on input holds nrecrd sorted records.  record irec starts
c     location recps1(irec) and contains recln1(irec) row and column 
c     indices.  each record is sorted first by columns and then by rows.
c     on output the file then labeled wafil2 will hold all of row
c     and column data in one sorted record without duplicates.
c
c     created         11-dec-96   -- rgg --
c     modified        
c
c     input arguments
c     ---------------
c
c     nrecrd      i   number of records originally on wafil1
c     sort        i   flag on whether to perform sort/merge.
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
c     rowls1      i   array of length lset to hold row    
c                     indicies from list 1
c     colls1      i   array of length lset to hold column 
c                     indicies from list 1
c     rowls2      i   array of length lset to hold row    
c                     indicies from list 2
c     colls2      i   array of length lset to hold column 
c                     indicies from list 2
c     rowls3      i   array of length 2*lset to hold row    
c                     indicies from list 3
c     colls3      i   array of length 2*lset to hold column 
c                     indicies from list 3
c
c     output arguments
c     ----------------
c
c     walen1      i   length of i/o file wafil1
c     watrn1      i   amount of i/o transfer to and from wafil1
c     walen2      i   length of i/o file wafil2
c     watrn2      i   amount of i/o transfer to and from wafil2
c     nnzero      i   number of nonzeroes in final form
c     error       i   error flag
c                     =  0  normal return
c                     = -1  i/o error on wafil1 or wafil2
 
c     ------------------------------------------------------------------

c     -------------------------------
c     ... global variable declaration
c     -------------------------------

      integer             wafil1, nrecrd, sort  , wafil2, lset  , 
     1                    nnzero, error

      integer             recps1 (*),     recln1 (*),
     1                    recps2 (*),     recln2 (*),
     2                    rowls1 (*),     colls1 (*),
     3                    rowls2 (*),     colls2 (*),
     4                    rowls3 (*),     colls3 (*)

      integer             walen1, watrn1, walen2, watrn2

c     ------------------------------
c     ... local variable declaration
c     ------------------------------

      integer             aclen1, aclen2, itemp , k11len, k11pos,
     1                    k12len, k12pos, k2pos , len1  , len2  ,
     2                    nrcrd1, nrcrd2, orgwa1, rec1  , rec2  

c.debug
c     integer             i, ipos, j, len
c.debug

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
c.debug
c     write(6,'(/,"entering xislis")') 
c     call xislp3 ( 'record pos.', nrcrd1, recps1, 6 )
c     call xislp3 ( 'record len.', nrcrd1, recln1, 6 )
c.debug

      len1   = recps1(nrecrd) + recln1(nrecrd) - 1 
      aclen1 = 2*len1 
      walen1 = aclen1
      watrn1 = walen1

      len2   = 0
      aclen2 = 0
      walen2 = 0
      watrn2 = 0
c.debug
c     write(6,'(/,"entering xislis - sort, nrcrd1, len1 = ",3i8)') 
c    1                                sort, nrcrd1, len1
c.debug

c     ----------------------------------
c     ... test if sort/merge is required
c     ----------------------------------

      if ( sort .eq. 0 .or. nrcrd1 .eq. 1 ) then

c         -------------------------------------------------------
c         ... sort/merge not required because i/o file is already 
c             in sorted order.
c         -------------------------------------------------------

          nnzero = len1
          wafil1 = wafil2
          wafil2 = orgwa1
          return

      end if

c     ------------------------------------------------------------------

  100 continue
      nrcrd2 = 0
      k2pos  = 1
c.debug
c     write(6,'(/,"entering 100 loop in xislis")') 
c     write(6,'("wafil1, wafil2 = ", 2i8)') wafil1, wafil2
c
c     call xislp3 ( 'record pos.', nrcrd1, recps1, 6 )
c     call xislp3 ( 'record len.', nrcrd1, recln1, 6 )

c     ipos = 1

c     do 2 j = 1, nrcrd1

c         nnzero = recln1(j)
c         write(6,'("j, nnzero = ", 2i8)') j, nnzero

c         do 1 i = 1, nnzero, lset
c             len  = min ( nnzero - i + 1, lset )
c             call xislw1 ( wafil1, rowls1, colls1, 
c    1                      ipos, len, error )
c             ipos = ipos + len
c             call xislp3 ( 'rowls1', len, rowls1, 6 )
c             call xislp3 ( 'colls1', len, colls1, 6 )
c   1     continue

c   2 continue
c.debug

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
c.debug

          call xislir ( wafil1, k11pos, k11len, k12pos, k12len,
     1                  wafil2, k2pos,  lset  , rowls1, colls1, 
     2                  rowls2, colls2, rowls3, colls3, error )
c.debug
c     write(6,'("after xislir - nrcrd2, k2pos, error = ", 3i8)')
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
          aclen2 = 2*len2

          watrn1 = watrn1 + aclen1          
 
          walen2 = max ( walen2, aclen2 )
          watrn2 = watrn2 + aclen2
 
      else
 
          len1   = k2pos - 1
          aclen1 = 2*len1

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

c.debug
c     call xislp3 ( 'at end of xislis - record pos.', 
c    1              nrcrd2, recps2, 6 )
c     call xislp3 ( 'at end of xislis - record len.', 
c    1              nrcrd2, recln2, 6 )
c.debug

      nnzero = recln2(1)

c     ------------------------------------------------------------------

      return
      end
