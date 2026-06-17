      subroutine xislir ( wafil1, rc1pos, rc1len, rc2pos, rc2len,
     1                    wafil2, k2pos , lset  , 
     2                    rowls1, colls1, rowls2, colls2,
     3                    rowls3, colls3, error )

c
c     purpose
c     -------
c
c     xislir performs the sort-merge of two records on wafil1 and 
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
c     error       i   error flag
c                     =  0  normal return
c                     = -1  i/o error on wafil1 or wafil2
c 
c     ------------------------------------------------------------------

c     -------------------------------
c     ... global variable declaration
c     -------------------------------

      integer             wafil1, rc1pos, rc1len, rc2pos, rc2len,
     1                    wafil2, k2pos , lset  , error

      integer             rowls1 (*),     colls1 (*),
     1                    rowls2 (*),     colls2 (*),
     1                    rowls3 (*),     colls3 (*)

c     ------------------------------
c     ... local variable declaration
c     ------------------------------

      integer             irow  , jcol  , k     , k1bgn ,
     1                    k1end , k2bgn , k2end , k3    , k3end ,
     2                    k3old , l     , l1    , l2    , 
     3                    locerr, lstcol, m1    , m2    , 
     3                    n1    , n2    , read1 , read2

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

      if ( m1 .ne. 0 ) then

        call xislw1 ( wafil1, rowls1(n1+1), colls1(n1+1), 
     1                rc1pos, m1, error )
        if ( error .ne. 0 ) go to 8000

        rc1pos = rc1pos + m1
        read1  = read1  + m1
        n1     = n1 + m1

      endif
c     ----------------------------------
c     ... read in next chunk of record 2
c     ----------------------------------

      m2 = min ( l2 - n2, rc2len - read2 )

      if ( m2 .eq. rc2len - read2 ) q2end = .true.
c.debug
c     write(6,'("after 120 - n2, m2, l2, rc2len, read2 = ", 5i8)')
c    1                        n2, m2, l2, rc2len, read2 
c.debug

      if ( m2 .ne. 0 ) then

        call xislw1 ( wafil1, rowls2(n2+1), colls2(n2+1), 
     1                rc2pos, m2, error )
        if ( error .ne. 0 ) go to 8000

        rc2pos = rc2pos + m2
        read2  = read2  + m2
        n2     = n2 + m2

      endif
c     ------------------------------------------------------------------

c     -----------------------------
c     ... current data is in memory
c     -----------------------------

c.debug 
c     write(6,'("at 140 - n1, n2 = ", 2i8)') n1, n2
c     call xislp3 ( 'row    list 1', n1, rowls1, 6 )
c     call xislp3 ( 'column list 1', n1, colls1, 6 )
c     call xislp3 ( 'row    list 2', n2, rowls2, 6 )
c     call xislp3 ( 'column list 2', n2, colls2, 6 )
c.debug 

      k1bgn = 1
      k2bgn = 1
      k3    = 0

      if ( n1 .eq. 0  .and.  n2 .eq. 0 ) return

      if ( n1 .eq. 0 ) go to 500

      if ( n2 .eq. 0 ) go to 400

      if ( q1end ) then
          lstcol = colls1(n1)
      else
          lstcol = colls1(n1) - 1
      end if

      if ( q2end ) then
          lstcol = min ( lstcol, colls2(n2) )
      else
          lstcol = min ( lstcol, colls2(n2) - 1 )
      end if

      if ( q1end .and. q2end ) then
          lstcol = max ( colls1(n1), colls2(n2) )
      end if

c     ---------------------------------------------------------
c     ... test of bulk moves.
c         1.  list 2 is null, move list 1 to list 3
c         2.  list 1 is null, move list 2 to list 3
c         3.  all columns in list 1 are before those in list 2,
c             move list 1 to list 3
c         4.  all columns in list 2 are before those in list 1,
c             move list 2 to list 3
c     ---------------------------------------------------------

  150 continue
      if ( k2bgn .gt. n2 .and. q2end ) go to 400

      if ( k1bgn .gt. n1 .and. q1end ) go to 500

      if ( colls1(n1) .lt. colls2(k2bgn) ) go to 400

      if ( colls2(n2) .lt. colls1(k1bgn) ) go to 500

c     --------------------------------------------
c     ... columns from both lists are interleaved.
c         process one column at a time.
c     --------------------------------------------

      jcol   = min ( colls1(k1bgn), colls2(k2bgn) ) 

      if ( jcol .gt. lstcol ) go to 600

c     -----------------------------------------------------------
c     ... copy all entries from first list with column index jcol
c         to list 3.
c     -----------------------------------------------------------

  200 continue
      k3old = k3
      q1    = .false.
      q2    = .false.
c.debug 
c     write(6,'("after 200 - jcol, lstcol, k3     = ", 4i8)') 
c    1                        jcol, lstcol, k3
c     write(6,'("            k1bgn, n1, k2bgn, n2 = ", 4i8)') 
c    1                        k1bgn, n1, k2bgn, n2
c.debug 

      k1end = k1bgn - 1

      do k = k1bgn, n1
          if ( colls1(k) .ne. jcol ) exit
              k3         = k3 +1
              colls3(k3) = jcol
              rowls3(k3) = rowls1(k)
              k1end      = k
              q1         = .true.
      enddo

c.debug 
c     write(6,'("after 220 - jcol, k3, k2bgn, n2 = ", 4i8)') 
c    1                        jcol, k3, k2bgn, n2
c     call xislp3 ( 'row    list 3', k3, rowls3, 6 )
c     call xislp3 ( 'column list 3', k3, colls3, 6 )
c.debug 

      k2end = k2bgn - 1

      do k = k2bgn, n2
          if ( colls2(k) .ne. jcol ) exit
              k3         = k3 +1
              colls3(k3) = jcol
              rowls3(k3) = rowls2(k)
              k2end      = k
              q2         = .true.
      enddo

c.debug 
c     write(6,'("after 240 - k3, jcol = ", 2i8)') k3, jcol
c     write(6,'("           k2bgn, n2 = ", 2i8)') k2bgn, n2
c     call xislp3 ( 'row    list 3', k3, rowls3, 6 )
c     call xislp3 ( 'column list 3', k3, colls3, 6 )
c.debug 

c     ------------------------------------------------------
c     ... if entries came from both lists, sort and compress
c     ------------------------------------------------------

      if ( q1 .and. q1 ) then

          k3end = k3
          l     = k3 - k3old
     
          call xislq1 ( l, rowls3(k3old+1), locerr )
          if ( locerr .ne. 0 ) go to 8100

          irow = 0
          k3   = k3old

          do k = k3old+1, k3end

              if ( rowls3(k) .eq. irow ) cycle
     
                  k3         = k3 + 1
                  rowls3(k3) = rowls3(k)
                  irow       = rowls3(k)

          enddo

      end if

      k1bgn = k1end + 1
      k2bgn = k2end + 1

c.debug 
c     write(6,'("after sort and compress - k3, jcol = ", 2i8)') 
c    1                                      k3, jcol
c     call xislp3 ( 'row    list 3', k3, rowls3, 6 )
c     call xislp3 ( 'column list 3', k3, colls3, 6 )
c.debug 

c     -------------------------------------
c     ... test on loop back for next column
c     -------------------------------------

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

      colls3(k3+1:k3+l) = colls1(k1bgn:k1bgn+l-1)
      rowls3(k3+1:k3+l) = rowls1(k1bgn:k1bgn+l-1)

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

      colls3(k3+1:k3+l) = colls2(k2bgn:k2bgn+l-1)
      rowls3(k3+1:k3+l) = rowls2(k2bgn:k2bgn+l-1)

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
c     write(6,'(// "********************")')
c     write(6,'("before write - k3 = ", i8)') k3
c     write(6,'("********************")')
c     call xislp3 ( 'row    list 3', k3, rowls3, 6 )
c     call xislp3 ( 'column list 3', k3, colls3, 6 )
c.debug 

      call xislw2 ( wafil2, rowls3, colls3, 
     1              k2pos, k3, error )
      if ( error .ne. 0 ) go to 8000

      k2pos = k2pos + k3

c     ----------------------------------------------------------
c     ... slide remainder to front of rowlst, collst, and marker
c     ----------------------------------------------------------

      n1 = n1 - k1bgn + 1

      call xislmv ( n1, rowls1, k1bgn, 1 )
      call xislmv ( n1, colls1, k1bgn, 1 )

      n2 = n2 - k2bgn + 1

      call xislmv ( n2, rowls2, k2bgn, 1 )
      call xislmv ( n2, colls2, k2bgn, 1 )

      go to 100

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
      stop

c     ==================================================================

 9000 continue
      return
      end
