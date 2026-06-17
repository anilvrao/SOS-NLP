      subroutine xislrw ( sqfile, ierr   )
 
c
c     purpose
c     -------
c
c     rewind a sequential unformatted fortran file.
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

c-----------------------------------------------------------------------
 
      ierr = 0
 
      rewind ( unit=sqfile, err=8000 )
 
      return
 
c-----------------------------------------------------------------------
 
 8000 continue
      ierr = -1
      return
 
c-----------------------------------------------------------------------
 
      end
