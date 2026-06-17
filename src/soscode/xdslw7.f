      subroutine xdslw7 ( lunit, ierr )
 

c     ----------------------------------------------------------------
c     ... xdslw7 closes files for the word addressable i/o package
c         comprised of subroutines xdslw5, xdslw6, xdslw7, xdslw8,
c         and xdslw9.  see xdslw8 for conversion notes.
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
