      subroutine xdslw9 ( ipos, length )
 

c     ----------------------------------------------------------------
c     ... xdslw9 increments or decrements the i/o address array for
c         the i/o package comprised of subroutines xdslw5, xdslw6, 
c         xdslw7, and xdslw8.  and xdslw9.  see xdslw8 for notes.
c     ----------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer             ipos(2), length
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c---------------------------------------------------------------------

      ipos(2) = ipos(2) + length

      if ( length .ge. 0 ) then

   10     continue
          if ( ipos(2) .gt. bfrsiz ) then
              ipos(1) = ipos(1) + 1
              ipos(2) = ipos(2) - bfrsiz
              go to 10
          endif

      else

   20     continue
          if ( ipos(2) .le. 0      ) then
              ipos(1) = ipos(1) - 1
              ipos(2) = ipos(2) + bfrsiz
              go to 20
          endif

      end if
 
c---------------------------------------------------------------------
 
      return
      end
