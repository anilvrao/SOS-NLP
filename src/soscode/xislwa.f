      subroutine   xislwa   ( mtxunt, format, l, array )
c
c     ==================================================================
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      character(len=*)    format
 
      integer             mtxunt, l
 
      integer             array(*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             i
 
c     ==================================================================
 
      if ( l .le. 0 ) return
 
      write ( mtxunt, fmt = format ) ( array (i), i = 1, l )
 
c     ==================================================================
 
      return
      end
