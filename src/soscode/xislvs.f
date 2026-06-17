      subroutine xislvs ( sqfile, ierr   )
 
c
c     purpose
c     -------
c
c     skip a record on a sequential unformatted fortran file by reading
c     the first piece of data on the file.  the data is discarded.
c
c     created            --  4-mar-98, rgg
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
c                                 i/o error   , ierr = -1
c
c-----------------------------------------------------------------------
 
      integer            ierr  , sqfile, idummy

c-----------------------------------------------------------------------
 
      ierr = 0
 
      read ( unit=sqfile, err=8000 ) idummy
 
      return
 
c-----------------------------------------------------------------------
 
 8000 continue
      ierr = -1
      return
 
c-----------------------------------------------------------------------
 
      end
