      subroutine xdslw1 ( lunit, option, ix, jx, x, ipos, n, ier )

c     ----------------------------------------------------------------
c     ... xdslw1 reads data for the word addressable i/o package
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
 
      if ( qdebug ) then
 
          write(6,'("in xdslw1 - lunit, option, ipos, n, bfrsiz = ",
     1               5i8)') lunit, option, ipos, n, bfrsiz
 
      end if
c.debug
 
c     ------------------------------------------------------------
c     ... compute first segment of bfrsiz and read in that segment
c     ------------------------------------------------------------
 
      irec = 1 + ( ipos - 1 ) / bfrsiz

      if ( option .eq. 1 ) then

          read ( lunit, rec=irec, err=800 ) ixbfr, jxbfr
c.debug
c     if ( qdebug ) then
c         call xislp3 ( '1st read ixbfr', bfrsiz, ixbfr, 6 )
c         call xislp3 ( '1st read jxbfr', bfrsiz, jxbfr, 6 )
c     end if
c.debug

      else if ( option .eq. 2 ) then

          read ( lunit, rec=irec, err=800 ) xbfr
c.debug
c     if ( qdebug ) then
c         call xdslp5 ( '1st read xbfr ', bfrsiz, xbfr , 6 )
c     end if
c.debug

      else if ( option .eq. 3 ) then

          read ( lunit, rec=irec, err=800 ) ixbfr, jxbfr, xbfr
c.debug
c     if ( qdebug ) then
c         call xislp3 ( '1st read ixbfr', bfrsiz, ixbfr, 6 )
c         call xislp3 ( '1st read jxbfr', bfrsiz, jxbfr, 6 )
c         call xdslp5 ( '1st read xbfr ', bfrsiz, xbfr , 6 )
c     end if
c.debug

      end if
 
c     -----------------------------
c     ... copy the data from z to x
c     -----------------------------
 
      ibgn = mod ( ipos-1, bfrsiz ) + 1
      k    = min ( bfrsiz-ibgn+1, n )
c.debug
      if ( qdebug ) then
          write(6,'("irec, ibgn, k = ", 3i8)') irec, ibgn, k
      end if
c.debug

      if ( option .eq. 1 .or. option .eq. 3 ) then 

          ix(1:k) = ixbfr(ibgn:ibgn+k-1)
          jx(1:k) = jxbfr(ibgn:ibgn+k-1)

      end if

      if ( option .ge. 2 ) then 

          x(1:k) = xbfr(ibgn:ibgn+k-1)

      end if
 
      if ( k .eq. n ) go to 200
 
c     -------------------------------------------------
c     ... if lots of data then read in blocks of bfrsiz
c     -------------------------------------------------
 
      nrec = ( n - k + bfrsiz - 1 ) / bfrsiz
 
      do jrec = 1, nrec
 
          irec = irec + 1

          if ( option .eq. 1 ) then

              read ( lunit, rec=irec, err=800 ) ixbfr, jxbfr
c.debug
c     if ( qdebug ) then
c         call xislp3 ( 'last read ixbfr', bfrsiz, ixbfr, 6 )
c         call xislp3 ( 'last read jxbfr', bfrsiz, jxbfr, 6 )
c     end if
c.debug

          else if ( option .eq. 2 ) then

              read ( lunit, rec=irec, err=800 ) xbfr
c.debug
c     if ( qdebug ) then
c         call xdslp5 ( 'last read xbfr ', bfrsiz, xbfr , 6 )
c     end if
c.debug

          else if ( option .eq. 3 ) then

              read ( lunit, rec=irec, err=800 ) ixbfr, jxbfr, xbfr
c.debug
c     if ( qdebug ) then
c         call xislp3 ( 'last read ixbfr', bfrsiz, ixbfr, 6 )
c         call xislp3 ( 'last read jxbfr', bfrsiz, jxbfr, 6 )
c         call xdslp5 ( 'last read xbfr ', bfrsiz, xbfr , 6 )
c     end if
c.debug

          end if
 
          l = min ( bfrsiz, n - k )

          if ( option .eq. 1 .or. option .eq. 3 ) then 

              ix(k+1:k+l) = ixbfr(1:l)
              jx(k+1:k+l) = jxbfr(1:l)

          end if

          if ( option .ge. 2 ) then 

              x(k+1:k+l) = xbfr(1:l)

          end if
 
          k = k + l
 
      enddo
 
c---------------------------------------------------------------------

  200 continue
c.debug
      if ( qdebug ) then
 
          write(6,'("leaving xdslw1")')
 
c         if ( option .eq. 1 .or. option .eq. 3 ) then
c             call xislp3 ( 'ix', n, ix, 6 )
c             call xislp3 ( 'jx', n, jx, 6 )
c         end if
 
c         if ( option .eq. 2 ) then
c             call xdslp5 ( 'x', n, x, 6 )
c         end if
 
      end if
c.debug
 
      return
 
c---------------------------------------------------------------------
 
  800 continue
      ier = -1
      return
      end
