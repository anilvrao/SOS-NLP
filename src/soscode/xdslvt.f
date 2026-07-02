      subroutine xdslvt ( wafil1, rc1pos, rc1len, rc2pos, rc2len,
     1                    wafil2, k2pos , lset  , coor11, coor12, 
     2                    value1, coor21, coor22, value2, coor31, 
     3                    coor32, value3,
     3                    qfinal, neqns , nnzero, xarndx, error )
c
c     purpose
c     -------
c
c     xdslvt performs the sort-merge of two records on wafil1 and 
c     places the resulting data on wafil2.
c
c     created         11-dec-96   -- rgg --
c     modified        
c
c     input arguments
c     ---------------
c
c     wafil1      i   unit number for first file.
c     rc1pos      i   i/o position at the start of first record
c                     on wafil1
c     rc1len      i   length of first record on wafil1
c     rc2pos      i   i/o position at the start of second record
c                     on wafil1
c     rc2len      i   length of second record on wafil1
c     wafil2      i   unit number for second file. 
c     lset        i   number of indicies that can be stored in 
c                     memory at any one time.
c     qfinal      l   logical flag indicating final pass or not
c     neqns       i   number of rows in org. matrix
c
c     working storage
c     ---------------
c
c     coor11      i   array of length lset to hold coordinate 1
c                     indicies from list 1
c     coor12      i   array of length lset to hold coordinate 2
c                     indicies from list 1
c     value1      d   array of length lset to hold values
c                     from list 1
c     coor21      i   array of length lset to hold coordinate 1
c                     indicies from list 2
c     coor22      i   array of length lset to hold coordinate 2
c                     indicies from list 2
c     value2      d   array of length lset to hold values
c                     from list 2
c     coor31      i   array of length 2*lset to hold coordinate 1
c                     indicies from list 3
c     coor32      i   array of length 2*lset to hold coordinate 2
c                     indicies from list 3
c     value3      d   array of length lset to hold values
c                     from list 3
c
c     input/output arguments
c     ----------------------
c
c     k2pos       i   on input, the starting position to write
c                     on wafil2.  on output, the next free position
c                     on wafil2.
c                     (Note:  these will toggle back and forth at
c                      each step)
c
c     output arguments
c     ----------------
c
c     nnzero      i   number of nonzero values in org. matrix
c     xarndx      i   array holding pointer for relative assembly
c                     indicies for each chevron 
c     error       i   error flag
c                     =  0  normal return
c                     = -1  i/o error on wafil1 or wafil2
c 
c     ------------------------------------------------------------------

c     -------------------------------
c     ... global variable declaration
c     -------------------------------

      integer             wafil1, rc1pos, rc1len, rc2pos, rc2len,
     1                    wafil2, k2pos , lset  , neqns , nnzero,
     2                    error

      logical             qfinal 

      integer             coor11 (*),     coor12 (*),
     1                    coor21 (*),     coor22 (*),
     2                    coor31 (*),     coor32 (*),
     3                    xarndx (*)

      double precision    value1 (*),     value2 (*),
     1                    value3 (*)

c     ------------------------------
c     ... local variable declaration
c     ------------------------------

      integer             count , i     , irow  , jchv  , k     , 
     1                    k1bgn , k1end , k2bgn , k2end , k3    , 
     2                    k3end , k3old , l     , l1    , l2    , 
     3                    lstchv, m1    , m2    , n1    , n2    ,
     4                    read1 , read2 

      logical             q1    , q1end , q2    , q2end

c     ------------------------------------------------------------------

c     ------------------
c     ... initialization
c     ------------------

      l1 = lset 
      l2 = lset

      n1 = 0
      n2 = 0

      q1end = .false.
      q2end = .false.

      read1 = 0
      read2 = 0

      if ( qfinal ) xarndx(1:neqns) = 0
c.debug
c     write(6,'("at start of xdslvt - qfinal = ",l8)') qfinal
c.debug

c     ------------------------------------------------------------------

c     ----------------------------------
c     ... read in next chunk of record 1
c     ----------------------------------

  100 continue
      m1 = min ( l1 - n1, rc1len - read1 )

      if ( m1 .eq. rc1len - read1 ) q1end = .true.
c.debug
c     write(6,'("after 100 - n1, m1, l1, rc1len, read1 = ", 5i8)')
c    1                        n1, m1, l1, rc1len, read1 
c.debug

      if ( m1 .eq. 0 ) go to 120

      call xdslw1 ( wafil1, 3, coor11(n1+1), coor12(n1+1), 
     1              value1(n1+1), rc1pos, m1, error )
      if ( error .ne. 0 ) go to 8000

      rc1pos = rc1pos + m1
      read1  = read1  + m1
      n1     = n1 + m1

c     ----------------------------------
c     ... read in next chunk of record 2
c     ----------------------------------

  120 continue
      m2 = min ( l2 - n2, rc2len - read2 )

      if ( m2 .eq. rc2len - read2 ) q2end = .true.
c.debug
c     write(6,'("after 120 - n2, m2, l2, rc2len, read2 = ", 5i8)')
c    1                        n2, m2, l2, rc2len, read2 
c.debug

      if ( m2 .eq. 0 ) go to 140

      call xdslw1 ( wafil1, 3, coor21(n2+1), coor22(n2+1), 
     1              value2(n2+1), rc2pos, m2, error )
      if ( error .ne. 0 ) go to 8000

      rc2pos = rc2pos + m2
      read2  = read2  + m2
      n2     = n2 + m2

c     ------------------------------------------------------------------

c     -----------------------------
c     ... current data is in memory
c     -----------------------------

  140 continue     
c.debug 
c     write(6,'("at 140 - n1, n2 = ", 2i8)') n1, n2
c     call xislp3 ( 'coor11 list', n1, coor11, 6 )
c     call xislp3 ( 'coor12 list', n1, coor12, 6 )
c     call xdslp5 ( 'value1 list', n1, value1, 6 )
c     call xislp3 ( 'coor21 list', n2, coor21, 6 )
c     call xislp3 ( 'coor22 list', n2, coor22, 6 )
c     call xdslp5 ( 'value2 list', n2, value2, 6 )
c.debug 

      k1bgn = 1
      k2bgn = 1
      k3    = 0

      if ( n1 .eq. 0  .and.  n2 .eq. 0 ) go to 700

      if ( n1 .eq. 0 ) go to 500

      if ( n2 .eq. 0 ) go to 400

      if ( q1end ) then
          lstchv = coor11(n1)
      else
          lstchv = coor11(n1) - 1
      end if

      if ( q2end ) then
          lstchv = min ( lstchv, coor21(n2) )
      else
          lstchv = min ( lstchv, coor21(n2) - 1 )
      end if

      if ( q1end .and. q2end ) then
          lstchv = max ( coor11(n1), coor21(n2) )
      end if

c     ----------------------------------------------------------
c     ... test of bulk moves.
c         1.  list 2 is null, move list 1 to list 3
c         2.  list 1 is null, move list 2 to list 3
c         3.  all chevrons in list 1 are before those in list 2,
c             move list 1 to list 3
c         4.  all chevrons in list 2 are before those in list 1,
c             move list 2 to list 3
c     ----------------------------------------------------------

  150 continue
      if ( k2bgn .gt. n2 .and. q2end ) go to 400

      if ( k1bgn .gt. n1 .and. q1end ) go to 500

      if ( coor11(n1) .lt. coor21(k2bgn) ) go to 400

      if ( coor21(n2) .lt. coor11(k1bgn) ) go to 500

c     ---------------------------------------------
c     ... chevrons from both lists are interleaved.
c         process one chevron at a time.
c     ---------------------------------------------

      jchv   = min ( coor11(k1bgn), coor21(k2bgn) ) 

      if ( jchv .gt. lstchv ) go to 600

c     ------------------------------------------------------------
c     ... copy all entries from first list with chevron index jchv
c         to list 3.
c     ------------------------------------------------------------

  200 continue
      k3old = k3
      q1    = .false.
      q2    = .false.
c.debug 
c     write(6,'(/)')
c     write(6,'("after 200 - jchv, lstchv, k3     = ", 4i8)') 
c    1                        jchv, lstchv, k3
c     write(6,'("            k1bgn, n1, k2bgn, n2 = ", 4i8)') 
c    1                        k1bgn, n1, k2bgn, n2
c.debug 

      k1end = k1bgn - 1

      do k = k1bgn, n1
          if ( coor11(k) .ne. jchv ) exit
              k3         = k3 +1
              coor31(k3) = jchv
              coor32(k3) = coor12(k)
              value3(k3) = value1(k)
              k1end      = k
              q1         = .true.
      enddo

c.debug 
c     write(6,'("after 220 - jchv, k3, k2bgn, n2 = ", 4i8)') 
c    1                        jchv, k3, k2bgn, n2
c     call xislp3 ( 'coor31 list', k3, coor31, 6 )
c     call xislp3 ( 'coor32 list', k3, coor32, 6 )
c     call xdslp5 ( 'value3 list', k3, value3, 6 )
c.debug 

      k2end = k2bgn - 1

      do k = k2bgn, n2
          if ( coor21(k) .ne. jchv ) exit
              k3         = k3 +1
              coor31(k3) = jchv
              coor32(k3) = coor22(k)
              value3(k3) = value2(k)
              k2end      = k
              q2         = .true.
      enddo

c.debug 
c     write(6,'("after 240 - k3, jchv = ", 2i8)') k3, jchv
c     write(6,'("           k2bgn, n2 = ", 2i8)') k2bgn, n2
c     write(6,'("           q1,    q2 = ", 2l8)') q1, q2        
c     call xislp3 ( 'coor31 list', k3, coor31, 6 )
c     call xislp3 ( 'coor32 list', k3, coor32, 6 )
c     call xdslp5 ( 'value3 list', k3, value3, 6 )
c.debug 

c     ------------------------------------------------------
c     ... if entries came from both lists, sort and compress
c     ------------------------------------------------------

      if ( q1 .and. q2 ) then

          k3end = k3
          l     = k3 - k3old
c.debug 
c     write(6,'("before sorting")')
c     call xislp3 ( 'coor31 list', l, coor31(k3old+1), 6 )
c     call xislp3 ( 'coor32 list', l, coor32(k3old+1), 6 )
c     call xdslp5 ( 'value3 list', l, value3(k3old+1), 6 )
c.debug 
     
          call xdslq2 ( l, coor32(k3old+1), value3(k3old+1) )

          irow = 0
          k3   = k3old

          do k = k3old+1, k3end

              if ( coor32(k) .ne. irow ) then

                  k3         = k3 + 1
                  coor32(k3) = coor32(k)
                  irow       = coor32(k)
                  value3(k3) = value3(k)

              else
     
                  value3(k3) = value3(k3) + value3(k)

              end if

          enddo
c.debug 
c     write(6,'("after sort and compress - k3, jchv = ", 2i8)') 
c    1                                      k3, jchv
c     l = k3 - k3old
c     call xislp3 ( 'coor31 list', l, coor31(k3old+1), 6 )
c     call xislp3 ( 'coor32 list', l, coor32(k3old+1), 6 )
c     call xdslp5 ( 'value3 list', l, value3(k3old+1), 6 )
c.debug 

      end if

      k1bgn = k1end + 1
      k2bgn = k2end + 1

c.debug 
c     write(6,'("after sort and compress - k3, jchv = ", 2i8)') 
c    1                                      k3, jchv
c     call xislp3 ( 'coor31 list', k3, coor31, 6 )
c     call xislp3 ( 'coor32 list', k3, coor32, 6 )
c     call xdslp5 ( 'value3 list', k3, value3, 6 )
c.debug 

c     --------------------------------------
c     ... test on loop back for next chevron
c     --------------------------------------

c.debug
c     write(6,'("before loop back test")')
c     write(6,'("k1bgn, n1, k2bgn, n2 = ", 4i8)') 
c    1            k1bgn, n1, k2bgn, n2
c.debug

      if ( k1bgn .le. n1  .and.  k2bgn .le. n2 ) go to 150

c     -----------------------
c     ... test for bulk moves
c     -----------------------

c.debug 
c     write(6,'("k1bgn, k2bgn, q1end, q2end         = ", 2i8,2l8)')
c    1            k1bgn, k2bgn, q1end, q2end
c.debug 

      if ( k1bgn .gt. n1  .and.  q1end ) go to 500

      if ( k2bgn .gt. n2  .and.  q2end ) go to 400

      go to 600

c     -----------------------------------------------------------------

c     --------------------------------------------------------
c     ... if list 2 is null copy remainder of list 1 to list 3
c     --------------------------------------------------------


  400 continue
      l = n1 - k1bgn + 1
c.debug
c     write(6,'("bulk move of list 1 to list 3")')
c     write(6,'("n1, k1bgn, l = ", 3i8)') n1, k1bgn, l
c.debug

      coor31(k3+1:k3+l) = coor11(k1bgn:k1bgn+l-1)
      coor32(k3+1:k3+l) = coor12(k1bgn:k1bgn+l-1)
      value3(k3+1:k3+l) = value1(k1bgn:k1bgn+l-1)

      k1bgn = k1bgn + l
      k3    = k3    + l

      go to 600

c     -----------------------------------------------------------------

c     --------------------------------------------------------
c     ... if list 1 is null copy remainder of list 2 to list 3
c     --------------------------------------------------------

  500 continue

      l = n2 - k2bgn + 1
c.debug
c     write(6,'("bulk move of list 2 to list 3")')
c     write(6,'("n2, k2bgn, l = ", 3i8)') n2, k2bgn, l
c.debug

      coor31(k3+1:k3+l) = coor21(k2bgn:k2bgn+l-1)
      coor32(k3+1:k3+l) = coor22(k2bgn:k2bgn+l-1)
      value3(k3+1:k3+l) = value2(k2bgn:k2bgn+l-1)

      k2bgn = k2bgn + l
      k3    = k3    + l

      go to 600

c     -----------------------------------------------------------------

c     ---------------------------------------------------------
c     ... write out list 3 and slide remainder of lists 1 and 2
c         to beginning.
c     ---------------------------------------------------------

  600 continue

c.debug
c     if ( qfinal ) then
c     write(6,'(// "********************")')
c     write(6,'("before write - k3 = ", i8)') k3
c     write(6,'("********************")')
c     call xislp3 ( 'coor31 list', k3, coor31, 6 )
c     call xislp3 ( 'coor32 list', k3, coor32, 6 )
c     call xdslp5 ( 'value3 list', k3, value3, 6 )
c     end if
c.debug 

      call xdslw2 ( wafil2, 3, coor31, coor32, value3, 
     1              k2pos, k3, error )
      if ( error .ne. 0 ) go to 8000

      k2pos = k2pos + k3

c     -------------------------------------------------------------
c     ... count number of entries that were written in each chevron
c     -------------------------------------------------------------

      if ( qfinal ) then

          jchv  = coor31(1)
          count = 1

          do i = 2, k3

              if ( coor31(i) .eq. jchv ) then
                  count = count + 1
              else
                  xarndx(jchv) = xarndx(jchv) + count
                  count        = 1
                  jchv         = coor31(i)
              end if

          enddo

          xarndx(jchv) = xarndx(jchv) + count

      end if

c     ----------------------------------------------------------
c     ... slide remainder to front of coor11, coor12, value1,
c         coor21, coor22, and value2.
c     ----------------------------------------------------------

      n1 = n1 - k1bgn + 1

      call xislmv ( n1, coor11, k1bgn, 1 )
      call xislmv ( n1, coor12, k1bgn, 1 )
      call xdslmv ( n1, value1, k1bgn, 1 )

      n2 = n2 - k2bgn + 1

      call xislmv ( n2, coor21, k2bgn, 1 )
      call xislmv ( n2, coor22, k2bgn, 1 )
      call xdslmv ( n2, value2, k2bgn, 1 )

      go to 100

c     -----------------------------------------------------------------
 
c     -------------------------------------------------------------
c     ... build the point version of xarndx from the chevron counts
c     -------------------------------------------------------------
 
  700 continue
      k = 1
 
      do i = 1, neqns
          count     = xarndx(i)
          xarndx(i) = k
          k         = k + count
      enddo
 
      xarndx(neqns+1) = k
      nnzero          = k - 1
 
c.debug
c     call xislp3 ( 'xarndx at end of xdslvt', neqns+1, xarndx, 6 )
c.debug
 
      return
 
c     ==================================================================

c     -------------
c     ... i/o error
c     -------------

 8000 continue
      error = -1
      go to 9000

c     ---------------------------
c     ... unforseen sorting error
c     ---------------------------

 8100 continue
      RETURN

c     ==================================================================

 9000 continue
      return
      end
