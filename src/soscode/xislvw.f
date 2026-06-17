      subroutine xislvw ( sqfile, n    , v     , ierr  )
 
c
c     purpose
c     -------
c
c     write a vector as a single unformatted record on a
c     sequential fortran file.
c
c     created            -- 12-jun-86, jgl
c     last modifications -- 15-jun-86, cca
c                           20-may-87, cca
c
c     input arguments
c     ---------------
c
c     sqfile      i   sequential file unit number
c     n           i   length of the vector v
c     v           i   integer vector to be written
c
c     output arguments
c     ----------------
c
c     ierr        i   error flag, normal write, ierr =  0
c                                 i/o error   , ierr = -2
c
c----------------------------------------------------------------------
 
      integer   ierr  , n     , sqfile
 
      integer   v(n)
 
c-----------------------------------------------------------------------
 
      ierr = 0
c.debug
c     write(6,'("in xislvw - sqfile, n = ", 2i10)') sqfile, n
c.debug
 
      if ( n .gt. 0 ) write ( sqfile, err=100 ) v
c.debug
c     write(6,'("in xislvr - read successful")') 
c.debug
 
      return
 
c-----------------------------------------------------------------------
 
  100 continue
      ierr = -2
c.debug
c     write(6,'("in xislvr - error encountered")') 
c.debug
      return
 
c-----------------------------------------------------------------------
 
      end
