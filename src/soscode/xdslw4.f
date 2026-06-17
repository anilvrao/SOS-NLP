      subroutine xdslw4  ( lunit, option, ierr )
 

c     ----------------------------------------------------------------
c     ... xdslw4 opens files for the word addressable i/o package
c         comprised of subroutines xdslw1, xdslw2, xdslw3, and xdslw4.
c         see below for conversion notes.
c     ----------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer             lunit,  option, ierr
 
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
c     the subroutines xdslw1, xdslw2, xdslw3, and xdslw4 form a
c     package of subroutines the simulate random access i/o by
c     using fortran-77 direct access i/o with fixed length 
c     records.
c         xdslw1 performs the reads
c         xdslw2 performs the writes
c         xdslw3 closes the file
c         xdslw4 opens the file
c
c     three types of files are written.
c
c     option = 1 - implies two integer data are written per
c                  entry
c     option = 2 - implies one double precision data is written per
c                  entry
c     option = 3 - implies two integer and one double precision data 
c                  is written per entry
c
c     this implementation copies the data into integer and floating
c     point buffers and then writes the buffers out without interleaving
c     the buffers.  A true implementation would interleave the data
c     depending on the option.
 
c---------------------------------------------------------------------
 
      ierr   = 0
 
      if ( option .lt. 1 .or. option .gt. 3 ) go to 800
 
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

      if ( option .eq. 1 ) then

c         -----------------------------------------
c         ... two integers vectors of length bfrsiz
c         -----------------------------------------

          reclen = 2 * recfac * bfrsiz

      else if ( option .eq. 2 ) then

c         ----------------------------------------------
c         ... one floating point vector of length bfrsiz
c         ----------------------------------------------

          reclen = xdslil ( recfac ) * bfrsiz

      else if ( option .eq. 3 ) then

c         ----------------------------------------------
c         ... two integer vectors and one floating point
c             vector of length bfrsiz
c         ----------------------------------------------

          reclen = ( 2 * recfac + xdslil ( recfac ) ) * bfrsiz

      end if
 
c     ---------------
c     ... open lunit.
c     ---------------

c.debug
c     write(6,'("in xdslw4 opening file - lunit, option, reclen = ",
c    1       3i8)')                        lunit, option, reclen 
c.debug

      open ( unit=lunit, access='direct', status='scratch',
     1       form='unformatted', recl=reclen, err=800 )
 
      return
 
c---------------------------------------------------------------------
 
  800 continue
      ierr = -1
      return
      end
