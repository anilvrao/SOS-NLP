      subroutine xdslvr( sqfile, n     , v     , ierr   )
 
c
c     purpose
c     -------
c
c     read a vector as a single unformatted record on a
c     sequential fortran file.
c
c     created            -- 04-dec-92, rggj
c     last modifications --
c
c     input arguments
c     ---------------
c
c     sqfile      i   sequential file unit number
c     n           i   length of the vector v
c     v           d   integer vector to be read
c
c     output arguments
c     ----------------
c
c     ierr        i   error flag, normal write, ierr =  0
c                                 i/o error   , ierr = -2
c                                 eof reached , ierr = -3
c
c-----------------------------------------------------------------------
 
      integer            ierr  , n     , sqfile
 
      double precision   v(n)
 
c-----------------------------------------------------------------------
 
      ierr = 0
 
      if ( n .gt. 0 ) read ( sqfile, err=100, end=200 ) v
 
      return
 
c-----------------------------------------------------------------------
 
  100 continue
      ierr = -2
      return
 
  200 continue
      ierr = -3
      return
 
c-----------------------------------------------------------------------
 
      end
