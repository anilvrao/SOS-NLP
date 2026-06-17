      subroutine xislw2 ( lunit, ix, jx, ipos, n, ier )
 

c     ----------------------------------------------------------------
c     ... xislw4 writes data for the word addressable i/o package
c         comprised of subroutines xislw1, xislw2, xislw3, and xislw4.
c         see xislw4 for conversion notes.
c     ----------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer             lunit,  ipos,   n,      ier

      integer             ix(*),  jx(*)
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     -------------------
c     ... local variables
c     -------------------

      integer             ixbfr(bfrsiz), jxbfr(bfrsiz)
 
      integer             ibgn, irec, jrec, k, l, nrec

c.debug
      logical             qdebug
c.debug
 
c---------------------------------------------------------------------
 
      ier = 0
 
      if (n .le. 0) return
 
c.debug
      qdebug = .false.
c     if ( ipos .le. 825 .and. 825 .le. ipos+n-1 ) qdebug = .true.

      if ( qdebug ) then

          write(6,'("in xislw2 - lunit, ipos, n, bfrsiz = ",
     1               5i8)') lunit, ipos, n, bfrsiz 

          call xislp3 ( 'ix', n, ix, 6 )
          call xislp3 ( 'jx', n, jx, 6 )

      end if

c.debug
 
c     -----------------------------------------------
c     ... compute first segment of bfrsiz and read in
c     -----------------------------------------------
 
      irec = 1 + ( ipos - 1 ) / bfrsiz
 
c     -------------------------------------------------------
c     ... copy data from ix and jx into buffers and write out
c     -------------------------------------------------------
 
      ibgn = mod ( ipos-1, bfrsiz ) + 1
      k    = min ( bfrsiz-ibgn+1, n )
c.debug
      if ( qdebug ) then
          write(6,'("irec, ibgn, k = ", 3i8)') irec, ibgn, k
      end if
c.debug

c     --------------------------------------------------
c     ... first read in first buffer to capture old data
c     --------------------------------------------------
 
      read ( lunit, rec=irec, err=10 ) ixbfr, jxbfr
c.debug
      if ( qdebug ) then
          call xislp3 ( '1st read ixbfr', bfrsiz, ixbfr, 6 )
          call xislp3 ( '1st read jxbfr', bfrsiz, jxbfr, 6 )
      end if
c.debug

c     -----------------------------
c     ... fill buffer with new data
c     -----------------------------

   10 continue
      ixbfr(ibgn:ibgn+k-1) = ix(1:k)
      jxbfr(ibgn:ibgn+k-1) = jx(1:k)
 
c     ------------------------
c     ... write the buffer out
c     ------------------------
 
      write ( lunit, rec=irec, err=800 ) ixbfr, jxbfr
c.debug
      if ( qdebug ) then
          call xislp3 ( 'write ixbfr', bfrsiz, ixbfr, 6 )
          call xislp3 ( 'write jxbfr', bfrsiz, jxbfr, 6 )
      end if
c.debug
 
      if ( k .eq. n ) return
 
c     -----------------------------------------
c     ... write any additional blocks of bfrsiz
c     -----------------------------------------
 
      nrec = ( n - k + bfrsiz - 1 ) / bfrsiz
 
      do jrec = 1, nrec
 
          irec = irec + 1
 
          l = min ( bfrsiz, n - k )
 
c         ---------------------------
c         ... try to read last record
c         ---------------------------
 
          if ( l .ne. bfrsiz ) then
 
              read ( lunit, rec=irec, err=50 ) ixbfr, jxbfr
c.debug
      if ( qdebug ) then
          call xislp3 ( 'last read ixbfr', bfrsiz, ixbfr, 6 )
          call xislp3 ( 'last read jxbfr', bfrsiz, jxbfr, 6 )
      end if
c.debug
 
          end if

c         -------------------------
c         ... fill buffer with data
c         -------------------------

   50     continue
          ixbfr(1:l) = ix(k+1:k+l)
          jxbfr(1:l) = jx(k+1:k+l)
 
          k = k + l

c         --------------------
c         ... write buffer out
c         --------------------
 
          write ( lunit, rec=irec, err=800 ) ixbfr, jxbfr
c.debug
      if ( qdebug ) then
          call xislp3 ( 'last write ixbfr', bfrsiz, ixbfr, 6 )
          call xislp3 ( 'last write jxbfr', bfrsiz, jxbfr, 6 )
      end if
c.debug
 
      enddo
 
      return
 
c---------------------------------------------------------------------
 
  800 continue
      ier = -1
      return
      end
