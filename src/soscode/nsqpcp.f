      subroutine   NSQPCP   ( output, msglvl, ndim  , ncon  , nHess ,
     1                        ldH   , ldA   , prbtyp, lintrm, cprbtp,
     2                        cset  , cmpcod )

c     ==================================================================
c     ==================================================================
c     ====  NSQPCP /                                                ====
c     ====  nsqpcp -- check scalar input parameters                 ====
c     ==================================================================
c     ==================================================================

c     ... parameters
      
      integer            output, msglvl, ndim  , ncon  , nHess , ldH   ,
     1                   ldA   , prbtyp, cmpcod

      logical            cset, lintrm
      
      character(len=2)   cprbtp
      
c     ... local variables

      integer            ntest   
      
c     ==================================================================

c     ... last modified 26-March-1996

c     ==================================================================

      if  ( ndim .le. 0 )  then
         cmpcod = -1001
         go to 10000
      endif
      
      if  ( ncon .lt. 0 )  then
         cmpcod = -1002
         go to 10000
      endif
      
      if  ( nHess .gt. 0 )  then
         ntest = nHess
      else
         ntest = ndim
      endif
      
      if  ( ldh .lt. ntest )  then
         cmpcod = -1003
         go to 10000
      endif
      
      if  ( lda .lt. ncon )  then
         cmpcod = -1004
         go to 10000
      endif

      if  ( prbtyp .lt. 0  .or.  prbtyp .gt. 2 )  then
         cmpcod = -1006
         go to 10000
      else
         if  ( prbtyp .eq. 0 )  then
            cprbtp = 'FP'
            if  ( lintrm )  then
               cmpcod = -1007
               go to 10000
            else
               cset   = .false.
            endif
         else
     1   if  ( prbtyp .eq. 1 )  then
            cprbtp = 'LP'
            if  ( lintrm )  then
               cset   = .true.
            else
               cmpcod = -1007
               go to 10000
            endif
         else
            cprbtp = 'QP'
            cset   = lintrm
         endif
      endif

      if  ( msglvl .lt. 0  .and.  output .gt. 0 )  then
         cmpcod = -1008
         go to 10000
      endif

10000 continue
      return

c     end of NSQPCP / nsqpcp

      end
