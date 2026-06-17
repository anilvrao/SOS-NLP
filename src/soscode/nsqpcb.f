      subroutine   NSQPCB   ( ndim  , ncon  , bigbnd, bupr  , blwr  ,
     1                        istatc, istatv, xupr  , xlwr  , nfixvr,
     2                        mequal, nclin , cmpcod )
c     ==================================================================
c     ==================================================================
c     ====  NSQPCB /                                                ====
c     ====  nsqpcb -- check bounds specifications                   ====
c     ==================================================================
c     ==================================================================

c     ... parameters
      
      integer            ndim  , ncon  , nfixvr, mequal, nclin ,
     1                   cmpcod

      double precision   bigbnd
      
      integer            istatc (*), istatv (ndim)

      double precision   bupr   (*), blwr   (*),
     1                   xupr   (ndim), xlwr   (ndim)

c     ... local variables

      integer            i, mignor  
      
c     ==================================================================

c     ... last modified 03-Sept-1996

c     ==================================================================

      mequal = 0
      mignor = 0

c     -----------------------------------------------------------------
c     ... look for general linear equality constraints and check bounds
c     -----------------------------------------------------------------
c 
      do i = 1, ncon

         if  ( istatc (i) .lt. -2  .or.  istatc (i) .gt. 4 )  then

            cmpcod = -1014
            go to 10000

         else
     1   if  ( istatc (i) .ne. 4 )  then
              
            if  ( bupr(i) .lt. blwr(i)   .or.   
     1            ( bupr(i) .ge. bigbnd  .and.  blwr(i) .le. -bigbnd ) )
     2      then

c              ... constraint bounds are incorrect

               cmpcod = -1010
               go to 10000
               
            else
     1      if  ( bupr(i) .eq. blwr(i) )  then
c                
c              ... constraint is an equality
c              
               mequal = mequal + 1
c              
            endif
c           
         else
            
c           ... count ignored constraints
            
            mignor = mignor + 1
c           
         endif
         
      enddo

c     -------------------------------------------------------
c     ... check variable bounds and check for fixed variables
c     -------------------------------------------------------
      
      nfixvr = 0
      
      do i = 1, ndim
         
         if  ( istatv (i) .lt. -2  .or.  istatv (i) .gt. 3 )  then

            cmpcod = -1014
            go to 10000

         else
     1   if  ( xupr(i) .lt. xlwr(i) )  then

c           ... variable bounds are incorrect

            cmpcod = -1011
            go to 10000

         else
     1   if  ( xupr(i) .eq. xlwr(i) )  then

c           ... variable bound is fixed 

            nfixvr = nfixvr + 1

         endif

      enddo

      nclin = ncon - mignor

10000 continue
      return

c     end of NSQPCB / nsqpcb
      
      end
