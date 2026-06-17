      subroutine xislw4  ( lunit, ierr )
 

c     ----------------------------------------------------------------
c     ... xislw4 opens files for the word addressable i/o package
c         comprised of subroutines xislw1, xislw2, xislw3, and xislw4.
c         see xislw4 for conversion notes.
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
 
      integer             jhccon
 
      external            jhccon
 
c---------------------------------------------------------------------

c     --------------------
c     ... conversion notes
c     --------------------
c
c     the subroutines xislw1, xislw2, xislw3, and xislw4 form a
c     package of subroutines the simulate random access i/o by
c     using fortran-77 direct access i/o with fixed length 
c     records.
c         xislw1 performs the reads
c         xislw2 performs the writes
c         xislw3 closes the file
c         xislw4 opens the file
c
c     two integer vectors are written per entry.
c
c     this implementation copies the data into integer 
c     buffers and then writes the buffers out without interleaving
c     the buffers.  A true implementation would interleave the data.
 
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

      reclen = 2 * recfac * bfrsiz
 
c     ---------------
c     ... open lunit.
c     ---------------

c.debug
c     write(6,'("in xislw4 opening file - lunit, reclen = ",
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
