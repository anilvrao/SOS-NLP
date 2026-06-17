      subroutine xdslw8  ( lunit, ierr )
 

c     ----------------------------------------------------------------
c     ... xdslw8 opens files for the word addressable i/o package
c         comprised of subroutines xdslw5, xdslw6, xdslw7, xdslw8,
c         see below for conversion notes.
c     ----------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer             lunit,  ierr
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             recfac, reclen
 
      logical             qopen
 
      integer             xdslil, jhccon
 
      external            xdslil, jhccon
 
c---------------------------------------------------------------------

c     --------------------
c     ... conversion notes
c     --------------------
c
c     the subroutines xdslw5, xdslw6, xdslw7, and xdslw8 form a
c     package of subroutines the simulate random access i/o by
c     using fortran-77 direct access i/o with fixed length 
c     records.
c         xdslw5 performs the reads
c         xdslw6 performs the writes
c         xdslw7 closes the file
c         xdslw8 opens the file
c         xdslw9 increments or decrements the i/o address array
c
c     this package differs from xdslw1-xdslw4 in that the 
c     i/o address is a buffer/address in buffer length 2 array.
c
c---------------------------------------------------------------------
 
      ierr   = 0
 
c     -----------------------------------
c     ... inquire on the status of lunit.
c     -----------------------------------
 
      inquire ( unit=lunit, err=800, opened=qopen )
 
      if ( qopen ) then
 
c         ---------------------------------------
c         ... lunit is already opened.  close it.
c         ---------------------------------------
 
          close ( unit=lunit, status='delete', err=800 )
 
      end if
 
c     ----------------------------------------------------
c     ... compute the record length for the open statement
c     ----------------------------------------------------
 
      recfac = jhccon ( 2 )

      reclen = xdslil ( recfac ) * bfrsiz
 
c     ---------------
c     ... open lunit.
c     ---------------

c.debug
c     write(6,'("in xdslw8 opening file - lunit, reclen = ",
c    1       3i8)')                        lunit, reclen 
c.debug

      open ( unit=lunit, access='direct', status='scratch',
     1       form='unformatted', recl=reclen, err=800 )
 
      return
 
c---------------------------------------------------------------------
 
  800 continue
      ierr = -1
      return
      end
