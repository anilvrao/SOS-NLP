      subroutine   CUNPRM  ( nvbl, ncon, a, lda, state, bupr, blwr,
     1                       cexchv )
c     ==================================================================
c     ==================================================================
c     ====  CUNPRM /                                                ====
c     ====  cunprm -- unpermute constraint matrix and vectors       ====
c     ==================================================================
c     ==================================================================

      integer            nvbl, ncon, lda

      integer            cexchv (*), state (*)

      double precision   a (lda, nvbl)

      double precision   bupr (*), blwr (*)

c     ... local variables

      integer            i, j

      double precision   t

c     ================================================================

c     ... last modified 26-March-1996
      
c     =================================================================

      do i = 1, ncon

         if  ( cexchv (i) .ne. i )  then

            j                  = state (i)
            state (i)          = state (cexchv (i))
            state (cexchv (i)) = j

            t                 = bupr (i)
            bupr (i)          = bupr (cexchv (i))
            bupr (cexchv (i)) = t

            t                 = blwr (i)
            blwr (i)          = blwr (cexchv (i))
            blwr (cexchv (i))  = t

            call dswap ( nvbl, a (i, 1), lda,
     1                         a (cexchv (i), 1), lda )

         endif
            
      enddo

      return

c     end of  CUNPRM / cunprm

      end
