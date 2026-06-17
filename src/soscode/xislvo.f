      subroutine xislvo ( sqfile, ierr   )
 
c
c     purpose
c     -------
c
c     open a sequential unformatted fortran file.
c
c     created            -- 16-dec-96, rgg
c     last modifications --
c
c     input arguments
c     ---------------
c
c     sqfile      i   sequential file unit number
c
c     output arguments
c     ----------------
c
c     ierr        i   error flag, normal open , ierr =  0
c                                 invalid unit, ierr = -1
c                                 i/o error   , ierr = -2
c
c-----------------------------------------------------------------------
 
      integer            ierr  , sqfile

      logical            qopen
 
c-----------------------------------------------------------------------
 
      ierr = 0
 
c     ------------------------------------
c     ... inquire on the status of sqfile.
c     ------------------------------------
 
      inquire ( unit=sqfile, err=8000, opened=qopen )

      if ( qopen ) then
 
c         ----------------------------------------
c         ... sqfile is already opened.  close it.
c         ----------------------------------------
 
          close ( unit=sqfile, status='delete', err=8000 )
 
      end if
 
c     ---------------
c     ... open lunit.
c     ---------------
 
      open ( unit=sqfile, err=8000, access='sequential',
     1       status='scratch', form='unformatted' )
 
      return
 
c-----------------------------------------------------------------------
 
 8000 continue
      ierr = -1
      return
 
c-----------------------------------------------------------------------
 
      end
