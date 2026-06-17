      subroutine xislvr( sqfile, n     , v     , ierr   )
 
c
c     purpose
c     -------
c
c     read a vector as a single unformatted record on a
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
c     v           i   integer vector to be read
c
c     output arguments
c     ----------------
c
c     ierr        i   error flag, normal write, ierr =  0
c                                 i/o error   , ierr = -2
c                                 eof reached , ierr = -3
c
c-----------------------------------------------------------------------
 
      integer   ierr  , n     , sqfile
 
      integer   v(n)
 
c-----------------------------------------------------------------------
 
      ierr = 0
c.debug
c     write(6,'("in xislvr - sqfile, n = ", 2i10)') sqfile, n
c.debug
 
      if ( n .gt. 0 ) read ( sqfile, err=100, end=200 ) v
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
 
  200 continue
      ierr = -3
c.debug
c     write(6,'("in xislvr - end encountered")') 
c.debug
      return
 
c-----------------------------------------------------------------------
 
      end
