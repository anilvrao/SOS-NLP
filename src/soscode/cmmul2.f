      subroutine   CMMUL2   ( msglvl, n, nZr, nZ,
     1                        zerolm, notopt, numinf, 
     2                        trusml, smllst, jsmlst,
     3                        tinyst, jtiny , gq,
     4                        iPrint )
c     ==================================================================
c     ==================================================================
c     ====  CMMUL2 /                                                ====
c     ====  cmmul2 -- "multipliers" for artificial constraints      ====
c     ==================================================================
c     ==================================================================

      integer            msglvl, n, nZr, nZ, notopt, numinf, jsmlst,
     1                   jtiny, iPrint

      double precision   zerolm, trusml, smllst, tinyst

      double precision   gq (n)

c     ==================================================================
c      
c     CMMUL2 / cmmul2  updates jsmlst and smllst when there are 
c     artificial constraints.
c
c     On input,  jsmlst  is the index of the minimum of the set of
c     adjusted multipliers.
c     On output, a negative jsmlst defines the index in Q'g of the
c     artificial constraint to be deleted.
c
c     Original version written 17-Jan-1988.
c     This version of cmmul2 dated  23-Jul-1991.
c
c     last modification 25-March-1996
c      
c     ==================================================================

      integer            j, k

      double precision   rlam

c     ==================================================================

      do j = nZr+1, nZ
         rlam = - abs( gq(j) )

         if  ( rlam .lt. zerolm )  then
            
            if  ( numinf .eq. 0)  then
               notopt = notopt + 1
            endif

            if  ( rlam .lt. smllst )  then
               trusml =   gq(j)
               smllst =   rlam
               jsmlst = - j
            end if

         else
     1   if  ( rlam .lt. tinyst )  then
            tinyst =   rlam
            jtiny  = - j
         end if
         
      enddo

      if  ( msglvl .ge. 20)  then
         write (iPrint, 1000) (gq(k), k=nZr+1,nZ)
      endif

      return

 1000 format(/ ' Multipliers for the artificial constraints        '
     $       / 4(5x, 1pe11.2))

c     end of hdqpm (cmmul2)
      
      end
