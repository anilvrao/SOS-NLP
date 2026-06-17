      subroutine   NSQPST   ( ndim  , ncon  , nclin , prbtyp, start , 
     1                        lniwrk, lnrwrk, nfixvr, mequal, xcexch,
     2                        xistat, xiwork, xlower, xupper, xlamda,
     3                        xax   , xrwork, maxact, mxfree, maxnZ ,
     4                        lniwk2, lnrwk2, iwork , cmpcod, needed )
c     ==================================================================
c     ==================================================================
c     ====  NSQPST /                                                ====
c     ====  nsqpst -- top level storage allocation                  ====
c     ==================================================================
c     ==================================================================

c     ... parameters
      
      integer           ndim  , ncon  , nclin , prbtyp, start , lniwrk, 
     1                  lnrwrk, nfixvr, mequal, xcexch, xistat, xiwork,
     2                  xlower, xupper, xlamda, xax   , xrwork, maxact,
     3                  mxfree, maxnZ , lniwk2, lnrwk2, cmpcod, needed

      integer           iwork (lniwrk)

c     ... local variables

      integer           ifree, iqpopt, ldQ, ldT, lenRT, ncolT, rqpopt

c     ==================================================================
c     ... last modification 24-July-1996
c     ==================================================================

c     ******************************************************************
c     ******************************************************************
c     ... NB:  ANY CHANGES MADE TO THIS SUBROUTINE SHOULD BE ECHOED
c              BY CHANGES TO  NLSPST.F
c     ******************************************************************
c     ******************************************************************

c     ==================================================================

      mxfree = iwork(1)
      mxfree = mequal
      mxfree = start

c     ... number of free variables
c         (variables not fixed on their bounds)
         
      mxfree = ndim - nfixvr

c     ... largest possible number of active constraints

      maxact = max (1, min (mxfree, nclin) )

c     ... largest possible dimension of nullspace (could incorporate
c         number of general equality constraints into this
c         calculation if we are sure that they are not linearly
c         dependent.  That is,
c               maxnZ = mxfree - minact = mxfree - mequal.  
c         Since we don't know that in general, we use a pessimistic
c         bound.)  The special case is from  qpdflt.

      maxnZ  = mxfree

      if  ( nclin .lt. ndim  .and.  prbtyp .lt. 2 )  then
         mxfree = nclin + 1
         maxnZ  = mxfree
      endif

c     ----------------------------------------------
c     ... storage allocation for the top level driver
c     ----------------------------------------------
      
      if  ( ncon .gt. nclin )  then
         xcexch = 1
         ifree  = ncon + 1
      else
         ifree = 1
      endif
         
      xlower = 1
      xupper = xlower + ndim + nclin
      xlamda = xupper + ndim + nclin
      xax    = xlamda + ndim + nclin
      xrwork = xax    + nclin


      if  ( nclin .eq. 0 )  then
         ldQ = 1
      else
         ldQ    = max (1, mxfree)
      endif
      
      ldT    = max ( maxnZ, maxact )
      ncolT  = mxfree
      lenRT  = ldT *ncolT
      rqpopt = 5*nclin + 8*ndim + lenRT + ldQ**2 
      
      if  ( xrwork + rqpopt  .gt. lnrwrk + 1 )  then
         cmpcod = -1100
         needed = xrwork + rqpopt  - 1
         go to 10000
      endif
      
      xistat = ifree
      xiwork = xistat + ndim + nclin

      iqpopt = 2*ndim + 3
      
      if  ( xiwork + iqpopt  .gt. lniwrk + 1 )  then
         cmpcod = -1101
         needed = xiwork + iqpopt  - 1

         go to 10000
      endif

      lnrwk2 = lnrwrk - xrwork + 1
      lniwk2 = lniwrk - xiwork + 1

10000 continue
      return

c     end of  NSQPST / nsqpst
      
      end
