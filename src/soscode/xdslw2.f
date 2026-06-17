      subroutine xdslw2 ( lunit, option, ix, jx, x, ipos, n, ier )
 

c     ----------------------------------------------------------------
c     ... xdslw2 writes data for the word addressable i/o package
c         comprised of subroutines xdslw1, xdslw2, xdslw3, and xdslw4.
c         see xdslw4 for conversion notes.
c     ----------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer             lunit,  option, ipos,   n,      ier

      integer             ix(*),  jx(*)

      double precision    x(*)
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     -------------------
c     ... local variables
c     -------------------

      integer             ixbfr(bfrsiz), jxbfr(bfrsiz)
 
      double precision    xbfr(bfrsiz)
 
      integer             ibgn, irec, jrec, k, l, nrec

c.debug
      logical             qdebug
c.debug
 
c---------------------------------------------------------------------
 
      ier = 0
 
      if (n .le. 0) return
 
      if ( option .lt. 1 .or. option .gt. 3 ) go to 800 

c.debug
      qdebug = .false.
c     if ( ipos .le. 825 .and. 825 .le. ipos+n-1 
c    1     .and. option .eq. 3 ) qdebug = .true.
c     if ( option .eq. 3 ) qdebug = .true.

      if ( qdebug ) then

          write(6,'("in xdslw2 - lunit, option, ipos, n, bfrsiz = ",
     1               5i8)') lunit, option, ipos, n, bfrsiz 

c         if ( option .eq. 1 .or. option .eq. 3 ) then
c             call xislp3 ( 'ix', n, ix, 6 )
c             call xislp3 ( 'jx', n, jx, 6 )
c         end if

c         if ( option .ge. 2 ) then
c             call xdslp5 ( 'x', n, x, 6 )
c         end if

      end if

c.debug
 
c     -----------------------------------------------
c     ... compute first segment of bfrsiz and read in
c     -----------------------------------------------
 
      irec = 1 + ( ipos - 1 ) / bfrsiz
 
c     ---------------------------------------
c     ... copy data from x to z and write out
c     ---------------------------------------
 
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
 
      if ( option .eq. 1 ) then
 
          read ( lunit, rec=irec, err=10 ) ixbfr, jxbfr
c.debug
c     if ( qdebug ) then
c         call xislp3 ( '1st read ixbfr', bfrsiz, ixbfr, 6 )
c         call xislp3 ( '1st read jxbfr', bfrsiz, jxbfr, 6 )
c     end if
c.debug
 
      else if ( option .eq. 2 ) then

          read ( lunit, rec=irec, err=10 ) xbfr
c.debug
c     if ( qdebug ) then
c         call xdslp5 ( '1st read xbfr ', bfrsiz, xbfr , 6 )
c     end if
c.debug

      else if ( option .eq. 3 ) then

          read ( lunit, rec=irec, err=10 ) ixbfr, jxbfr, xbfr
c.debug
c     if ( qdebug ) then
c         call xislp3 ( '1st read ixbfr', bfrsiz, ixbfr, 6 )
c         call xislp3 ( '1st read jxbfr', bfrsiz, jxbfr, 6 )
c         call xdslp5 ( '1st read xbfr ', bfrsiz, xbfr , 6 )
c     end if
c.debug

      end if

c     -----------------------------
c     ... fill buffer with new data
c     -----------------------------

   10 continue
      if ( option .eq. 1 .or. option .eq. 3 ) then
 
          ixbfr(ibgn:ibgn+k-1) = ix(1:k)
          jxbfr(ibgn:ibgn+k-1) = jx(1:k)
 
      end if
 
      if ( option .ge. 2 ) then
 
          xbfr(ibgn:ibgn+k-1) = x(1:k)
 
      end if

c     ------------------------
c     ... write the buffer out
c     ------------------------
 
      if ( option .eq. 1 ) then
 
          write ( lunit, rec=irec, err=800 ) ixbfr, jxbfr
c.debug
c     if ( qdebug ) then
c         call xislp3 ( 'write ixbfr', bfrsiz, ixbfr, 6 )
c         call xislp3 ( 'write jxbfr', bfrsiz, jxbfr, 6 )
c     end if
c.debug
 
      else if ( option .eq. 2 ) then

          write ( lunit, rec=irec, err=800 ) xbfr
c.debug
c     if ( qdebug ) then
c         call xdslp5 ( 'write xbfr ', bfrsiz, xbfr , 6 )
c     end if
c.debug

      else if ( option .eq. 3 ) then

          write ( lunit, rec=irec, err=800 ) ixbfr, jxbfr, xbfr
c.debug
c     if ( qdebug ) then
c         call xislp3 ( 'write ixbfr', bfrsiz, ixbfr, 6 )
c         call xislp3 ( 'write jxbfr', bfrsiz, jxbfr, 6 )
c         call xdslp5 ( 'write xbfr ', bfrsiz, xbfr , 6 )
c     end if
c.debug

      end if
 
      if ( k .eq. n ) return
 
c     -----------------------------------------
c     ... write any additional blocks of bfrsiz
c     -----------------------------------------
 
      nrec = ( n - k + bfrsiz - 1 ) / bfrsiz
 
      do 100 jrec = 1, nrec
 
          irec = irec + 1
 
          l = min ( bfrsiz, n - k )
 
c         ---------------------------
c         ... try to read last record
c         ---------------------------
 
          if ( l .ne. bfrsiz ) then
 
              if ( option .eq. 1 ) then
 
                  read ( lunit, rec=irec, err=50 ) ixbfr, jxbfr
c.debug
c     if ( qdebug ) then
c         call xislp3 ( 'last read ixbfr', bfrsiz, ixbfr, 6 )
c         call xislp3 ( 'last read jxbfr', bfrsiz, jxbfr, 6 )
c     end if
c.debug
 
              else if ( option .eq. 2 ) then

                  read ( lunit, rec=irec, err=50 ) xbfr
c.debug
c     if ( qdebug ) then
c         call xdslp5 ( 'last read xbfr ', bfrsiz, xbfr , 6 )
c     end if
c.debug

              else if ( option .eq. 3 ) then

                  read ( lunit, rec=irec, err=50 ) ixbfr, jxbfr, xbfr
c.debug
c     if ( qdebug ) then
c         call xislp3 ( 'last read ixbfr', bfrsiz, ixbfr, 6 )
c         call xislp3 ( 'last read jxbfr', bfrsiz, jxbfr, 6 )
c         call xdslp5 ( 'last read xbfr ', bfrsiz, xbfr , 6 )
c     end if
c.debug

              end if

          end if

c         -------------------------
c         ... fill buffer with data
c         -------------------------

   50     continue
          if ( option .eq. 1 .or. option .eq. 3 ) then
 
              ixbfr(1:l) = ix(k+1:k+l)
              jxbfr(1:l) = jx(k+1:k+l)
 
          end if
 
          if ( option .ge. 2 ) then
 
              xbfr(1:l) = x(k+1:k+l)
 
          end if

          k = k + l

c         --------------------
c         ... write buffer out
c         --------------------
 
          if ( option .eq. 1 ) then
     
              write ( lunit, rec=irec, err=800 ) ixbfr, jxbfr
c.debug
c     if ( qdebug ) then
c         call xislp3 ( 'last write ixbfr', bfrsiz, ixbfr, 6 )
c         call xislp3 ( 'last write jxbfr', bfrsiz, jxbfr, 6 )
c     end if
c.debug
 
          else if ( option .eq. 2 ) then

              write ( lunit, rec=irec, err=800 ) xbfr
c.debug
c     if ( qdebug ) then
c         call xdslp5 ( 'last write xbfr ', bfrsiz, xbfr , 6 )
c     end if
c.debug

          else if ( option .eq. 3 ) then

              write ( lunit, rec=irec, err=800 ) ixbfr, jxbfr, xbfr
c.debug
c     if ( qdebug ) then
c         call xislp3 ( 'last write ixbfr', bfrsiz, ixbfr, 6 )
c         call xislp3 ( 'last write jxbfr', bfrsiz, jxbfr, 6 )
c         call xdslp5 ( 'last write xbfr ', bfrsiz, xbfr , 6 )
c     end if
c.debug

          end if
 
  100 continue
 
      return
 
c---------------------------------------------------------------------
 
  800 continue
      ier = -1
      return
      end
