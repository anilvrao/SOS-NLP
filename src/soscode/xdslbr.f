      subroutine xdslbr ( lunit, xbfr, irec, bfrsiz, ier )
 

c     ----------------------------------------------------------------
c     ... xdslbr reads one block of data with buffer size bfrsiz
c         for the word addressable i/o package.
c         See xdslw4 and xdslw8 for detail information.
c
c     created -- 06-jul-2001, dkw
c
c     ----------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer             lunit, irec, bfrsiz, ier

      double precision    xbfr(bfrsiz)
 
c---------------------------------------------------------------------
 
      ier = 0

      read ( lunit, rec=irec, err=800 ) xbfr

      return
 
c---------------------------------------------------------------------
 
  800 continue
      ier = -1
      return
      end
