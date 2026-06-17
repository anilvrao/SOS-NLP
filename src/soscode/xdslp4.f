      subroutine   xdslp4   ( title, n, x, output )
 
c     ==================================================================
c     ====  xdslp4 -- print real vector in table (long format)      ====
c     ==================================================================
c     ==================================================================
c
c     xdslp4 prints out the real vector  x  of length  n  to  logical
c     unit output in a long format.  the character string in  title  is
c     printed as a title for the table
c
c     last modified   --- oct. 29, 1990   -- mlc --
c                     --- july 10, 1989   -- jgl --
c
c     --------------
c     ... parameters
c     --------------
 
      character(len=*)    title
 
      integer             n, output
 
      double precision    x(*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             i, l
 
      character(len=75)   line
 
c     ==================================================================
 
c     ---------------
c     ... write title
c     ---------------
 
      l = min ( len ( title ), 75 )
 
      do i = 1, l
          line(i:i) = '-'
      enddo
 
      do i = l+1, 75
          line(i:i) = ' '
      enddo
 
      write ( output, 2000 ) title (1:l), line (1:l)
 
c     ----------------------
c     ... write out the data
c     ----------------------
 
      write ( output, 2100 ) (x(i), i=1,n)
 
      return
 
c     -----------
c     ... formats
c     -----------
 
 2000 format ( /5x, a  /5x, a  / )
 
 2100 format ( 5x, 1p3d25.14 )
 
      end
