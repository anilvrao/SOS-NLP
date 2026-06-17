      subroutine xislw3 ( lunit, ierr )
 

c     ----------------------------------------------------------------
c     ... xislw4 closes files for the word addressable i/o package
c         comprised of subroutines xislw1, xislw2, xislw3, and xislw4.
c         see xislw4 for conversion notes.
c     ----------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer             lunit,  ierr
 
c---------------------------------------------------------------------
 
      ierr = 0
 
      close ( unit=lunit, status='delete', err=800 )
 
      return
 
c---------------------------------------------------------------------
 
  800 continue
      ierr = -1
      return
      end
