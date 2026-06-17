      subroutine xdslw5 ( lunit, x, ipos, n, ier )
 

c     ----------------------------------------------------------------
c     ... xdslw5 reads data for the word addressable i/o package
c         comprised of subroutines xdslw5, xdslw6, xdslw7, xdslw8,
c         and xdslw9.  see xdslw8 for conversion notes.
c
c     modified -- 04-0ct-2001  dkw, fixed to read in data directly
c                                   and avoid unnecessary call to copy
c                                   when there is more that 1 record
c
c     ----------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer             lunit,  n,      ier

      integer             ipos(2)
 
      double precision    x(*)
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------

      include '../commons/bcsext4.CMN'

c     -------------------
c     ... local variables
c     -------------------

      double precision    xbfr(bfrsiz)
 
      integer             ibgn, irec, jrec, k, l, nrec

c.debug
      logical             qdebug
c.debug
 
c---------------------------------------------------------------------
 
      ier = 0

      if (n .le. 0) return

c.debug
      qdebug = .false.
 
      if ( qdebug ) then
 
          write(6,'("in xdslw5 - lunit, ipos, n, bfrsiz = ",
     1               5i8)') lunit, ipos, n, bfrsiz
 
      end if
c.debug
 
c     ------------------------------------------------------------
c     ... compute first segment of bfrsiz and read in that segment
c     ------------------------------------------------------------
 
      irec = ipos(1)

      read ( lunit, rec=irec, err=800 ) xbfr
c.debug
c     if ( qdebug ) then
c         call xdslp5 ( '1st read xbfr ', bfrsiz, xbfr , 6 )
c     end if
c.debug

c     --------------------------------
c     ... copy the data from xbfr to x
c     --------------------------------
 
      ibgn = ipos(2)
      k    = min ( bfrsiz-ibgn+1, n )
c.debug
      if ( qdebug ) then
          write(6,'("irec, ibgn, k = ", 3i8)') irec, ibgn, k
      end if
c.debug

      x(1:k) = xbfr(ibgn:ibgn+k-1)

      if ( k .eq. n ) go to 200
 
c     -------------------------------------------------
c     ... if lots of data then read in blocks of bfrsiz
c     -------------------------------------------------
 
      nrec = ( n - k + bfrsiz - 1 ) / bfrsiz
 
      do jrec = 1, nrec - 1

          irec = irec + 1
          call xdslbr ( lunit, x(k+1), irec, bfrsiz, ier )
c.debug
c     if ( qdebug ) then
c         call xdslp5 ( 'last read xbfr ', bfrsiz, xbfr , 6 )
c     end if
c.debug

          if ( ier .ne. 0 ) then
             ier = -2
             return
          end if


          l = min ( bfrsiz, n - k )
          k = k + l
 
      enddo

c     -----------------------------------------------------------
c     ... compute last segment of bfrsiz and read in that segment
c     -----------------------------------------------------------

      if ( nrec .ge. 1 ) then
         irec = irec + 1
         read ( lunit, rec=irec, err=800 )  xbfr
         l = min ( bfrsiz, n - k )
         x(k+1:k+l) = xbfr(1:l)
      end if

c---------------------------------------------------------------------

  200 continue
c.debug
      if ( qdebug ) then
 
          write(6,'("leaving xdslw5")')
 
c         call xdslp5 ( 'x', n, x, 6 )
 
      end if
c.debug
 
      return
 
c---------------------------------------------------------------------
 
  800 continue
      ier = -1
      return
      end
