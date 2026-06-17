      subroutine xdslfk ( sqfile, unsym , savea , nsuper, xsup  ,
     1                    xadj  , rowlst, lclidx, anzp  , 
     2                    nrecrd, error )
 
c
c  purpose -- to spill the original matrix to sqfile to start minimum
c             memory processing.
c
c  created   -- 04-dec-92, rgg
c  revisions --
c
c  variables
c
c      sqfile -- sequential i/o file position at end-of-information
c      unsym  -- unsymmetric matrix switch.  
c      nsuper -- number of supernodes
c      xsup   -- supernode pointer array. since the matrix is assumed
c                to be ordered by a post-ordered traversal, the nodes
c                for supernode i are xsup(i), ... , xsup(i+1)-1
c      xadj   -- pointer array into lclidx, anzp
c      lclidx -- scatter indices for original matrix.
c      anzp   -- lower triangular values for original matrix.
c
c  output variables --
c
c      nrecrd -- number of records written to sqfile
c      error  -- error flag, error =  0, success
c                                  = -1, i/o error on sqfile
c
c  subprograms called
c
c      xislvw -- integer vector write
c      xdslvw -- vector write
c
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           sqfile, nsuper, nrecrd, error
 
      integer           xsup(*),        xadj(*),
     1                  rowlst(*),      lclidx(*)
 
      logical           unsym , savea
 
      double precision  anzp(*)
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer           isuper, jbgn  , jend  , kbgn  , klong
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      external          xislvw, xdslvw
c
c  =====================================================================
 
c     ---------------------------------
c     ... write entire matrix to sqfile
c     ---------------------------------
c.debug
c     write(6,'("in xdslfk - sqfile, unsym, savea  = ", i8, 2l8)')
c    1                        sqfile, unsym, savea      
c.debug
 
      do isuper = 1, nsuper
 
          jbgn  = xsup(isuper)
          jend  = xsup(isuper+1) - 1
c.debug
c     write(6,'("in xdslfk - isuper, jbgn, jend   = ", 3i8)')
c    1                        isuper, jbgn, jend 
c.debug
 
          kbgn  = xadj(jbgn)
          klong = xadj(jend+1) - kbgn
c.debug
c     write(6,'("            kbgn  , klong         = ", 3i8)')
c    1                        kbgn  , klong             
c.debug

 
          if ( klong .le. 0 ) cycle
 
          if ( savea ) then 
c.debug
c     call xislp3 ( 'rowlst', klong, rowlst(kbgn), 6 )
c.debug
              nrecrd = nrecrd + 1
              call xislvw ( sqfile, klong, rowlst(kbgn), error )
              if ( error .ne. 0 ) go to 8000
          end if
 
c.debug
c     call xislp3 ( 'lclidx', klong, lclidx(kbgn), 6 )
c.debug
          nrecrd = nrecrd + 1
          call xislvw ( sqfile, klong, lclidx(kbgn), error )
          if ( error .ne. 0 ) go to 8000
 
c.debug
c     call xdslp5 ( 'anzp  ', klong, anzp  (kbgn), 6 )
c.debug
          nrecrd = nrecrd + 1
          call xdslvw ( sqfile, klong, anzp(kbgn), error )
          if ( error .ne. 0 ) go to 8000
 
      enddo
 
      return
 
c  =====================================================================
 
c     ------------------
c     ... i/o error trap
c     ------------------
 
 8000 continue
      error = -1
      return
 
c  =====================================================================
 
      end
