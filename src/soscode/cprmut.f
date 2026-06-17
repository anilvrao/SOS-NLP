      subroutine   CPRMUT  ( nvbl, ncon, nclin, a, lda, state,
     1                       bupr, blwr, cexchv, error )
c     ==================================================================
c     ==================================================================
c     ====  hddqqp /                                                ====
c     ====  cprmut -- permute constraint matrix and vectors         ====
c     ==================================================================
c     ==================================================================

      integer            nvbl, ncon, nclin, lda

      logical            error

      integer            state (*), cexchv (*)

      double precision   a (lda, nvbl)

      double precision   bupr (*), blwr (*)

c     ... local variables

      integer            i, lastus, nextig, nonign

      double precision   t

c     ===================================================================

c     ... last modified 03-Sept-1996

c     ==================================================================
      
c     -------------------------------------------------------------
c     ... of the  ncon  general linear constraints, only  nclin
c         are not ignored.  rearrange the rows of the constraint
c         matrix so that the  nclin  non-ignored constraints are
c         listed first.
c
c         constraints to be ignored are signalled by state (i) = 4
c
c         rearrange is based on pairwise exchanges, swapping the
c         lowest numbered constraint to be ignored with the highest
c         numbered non-ignored constraint, assuming that this pair
c         is in fact out of order.
c     -------------------------------------------------------------

      do i = 1, ncon
         cexchv (i) = i
      enddo
      nonign = 0

      
      nextig = 0
      lastus    = ncon + 1

c     ========================================
c     repeat ... until  lastus < nextig
c     ========================================

c     ... find next constraint to be ignored
      
 200  continue
      nextig = nextig + 1
      if  ( nextig .lt. lastus )  then
         
         if  ( state (nextig) .ne. 4 )  then
            nonign = nonign + 1
            go to 200
         endif
         
      endif

c     ... find last number constraint to be used
 300  continue
      if  ( nextig .lt. lastus )  then

         lastus = lastus - 1
         if  ( state (lastus) .eq. 4 )  then
            lastus = lastus - 1
            go to 300

         endif

      endif

c     ... if we found such an out-of-order pair, exchange them

      if  ( nextig .lt. lastus )  then

         nonign = nonign + 1
         
         i                   = state (nextig)
         state (nextig) = state (lastus)
         state (lastus)    = i

         t                  = bupr (nextig)
         bupr (nextig) = bupr (lastus)
         bupr (lastus)    = t

         t                  = blwr (nextig)
         blwr (nextig) = blwr (lastus)
         blwr (lastus)    = t

         call dswap ( nvbl, a (nextig, 1), lda,
     1                      a (lastus   , 1), lda )

         cexchv (nextig) = lastus

         go to 200

      endif

c     =================================
c     ... until  lastus < nextig
c     =================================

c     ... check that we agreed on the number ...
      
      error = nonign .ne. nclin

      return
      
c     end of hddqqp / cprmut
      
      end
