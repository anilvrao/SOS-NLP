      subroutine xdslsr ( work, inuse, mxused, time, opcnts )
 
c
c     purpose
c     -------
c
c     xdslsr retrieves statistics from the global.CMNication array.
c
c     created         26-jan-89   -- rgg --
c     last modified
c
c     input arguments
c     ---------------
c
c     work        d   work array.  on input it contains the
c                    .CMNication area and all active arrays.
c
c     output arguments
c     ----------------
c
c     inuse       i   amount of workspace currently in use.
c     mxused      i   maximum amount of workspace used.
c     time        d   array of length 6 containing the time
c                     for the most current execution of each
c                     stage
c     opcnts      d   array of length 2 containing the operation
c                     counts for the factorization and solve.
c
c---------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      integer             inuse , mxused
 
      double precision    work(*), time(6), opcnts(2)
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c---------------------------------------------------------------------
c
 
c     -----------------------
c     ... retrieve statistics
c     -----------------------
 
      inuse     = work (qinuse)
      mxused    = work (qmxuse)
 
      opcnts(1) = work (qfctop)
      opcnts(2) = work (qslvop)
 
      time  (1) = work (qinptm)
      time  (2) = work (qordtm)
      time  (3) = work (qsfctm)
      time  (4) = work (qvaltm)
      time  (5) = work (qfcttm)
      time  (6) = work (qslvtm)
 
c---------------------------------------------------------------------
 
      return
      end
