      subroutine   NLSPST   ( ndim  , ncon  , estiwk, estrwk )
c     ==================================================================
c     ==================================================================
c     ====  NLSPST /                                                ====
c     ====  nlspst -- top level storage allocation estimator        ====
c     ==================================================================
c     ==================================================================

c     input:
c
c             NDIM    number of variables
c             NCON    number of general linear constraints
c
c
c     output:
c             ESTIWK  estimate (upper bound) for integer storage
c             ESTRWK  estimate (upper bound) for real storage

c     ==================================================================

c     ... parameters
      
      integer           ndim  , ncon  , estiwk, estrwk

c     ... local variables

      integer           idrive, iqpopt, ldQ   , ldT   , lenRT , maxact, 
     1                  mxfree, maxnZ , ncolT , rdrive, rqpopt

c     ==================================================================
c     ... last modification 24-July-1996
c     ==================================================================

c     ******************************************************************
c     ******************************************************************
c     ... NB:  ANY CHANGES MADE TO THIS SUBROUTINE SHOULD BE ECHOED
c              BY CHANGES TO  nsqpst.f 
c     ******************************************************************
c     ******************************************************************

c     ==================================================================

c     ... estimate number of free variables
c         (variables not fixed on their bounds)
         
      mxfree = ndim

c     ... largest possible number of active constraints

      maxact = max (1, min (mxfree, ncon) )

c     ... largest possible dimension of nullspace (could incorporate
c         number of general equality constraints into this
c         calculation if we are sure that they are not linearly
c         dependent.  That is,
c               maxnZ = mxfree - minact = mxfree - mequal.  
c         Since we don't know that in general, we use a pessimistic
c         bound.)  The special case is from  qpdflt.

      maxnZ  = mxfree

      
c     -----------------------------------------------------
c     ... integer storage required for the top level driver
c     -----------------------------------------------------

      idrive = 2*ncon + ndim
      iqpopt = 2*ndim + 3
      
      estiwk = idrive + iqpopt
      
c     --------------------------------------------------
c     ... real storage required for the top level driver
c     --------------------------------------------------

      rdrive = 4*ncon + 3*ndim

      if  ( ncon .eq. 0 )  then
         ldQ = 1
      else
         ldQ    = max (1, mxfree)
      endif
      
      ldT    = max ( maxnZ, maxact )
      ncolT  = mxfree
      lenRT  = ldT *ncolT
      rqpopt = 5*ncon + 8*ndim + lenRT + ldQ**2 

      estrwk = rdrive + rqpopt
      
c     end of  NLSPST / nlspst
      
      end
