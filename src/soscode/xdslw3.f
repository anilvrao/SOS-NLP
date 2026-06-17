      subroutine xdslw3 ( lunit, option, ierr )
 

c     ----------------------------------------------------------------
c     ... xdslw3 closes files for the word addressable i/o package
c         comprised of subroutines xdslw1, xdslw2, xdslw3, and xdslw4.
c         see xdslw4 for conversion notes.
c     ----------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer             lunit,  option, ierr
 
c---------------------------------------------------------------------
 
      ierr = 0
 
      close ( unit=lunit, status='delete', err=800 )
 
      return
 
c---------------------------------------------------------------------
 
  800 continue
      ierr = -1
      return
      end
