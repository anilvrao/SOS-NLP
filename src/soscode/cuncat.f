      subroutine   CUNCAT   ( ndim, mcon, istate, istatv, istatc, 
     1                        clamda, conmlt, varmlt, error )
c     ==================================================================
c     ==================================================================
c     ====  CUNCAT /                                                ====
c     ====  cuncat -- interface bcs and qpopt data conventions      ====
c     ==================================================================
c     ==================================================================

      integer            ndim, mcon

      logical            error

      integer            istatv (ndim), istatc (*), istate (*)

      double precision   varmlt (ndim), conmlt (*), clamda (*)

c     ================================================================

c	... last modified 03-Sept-1996

c     ================================================================

      integer            i

      integer            sttmap (-2:4)

      data               sttmap / -1, -2, 0, 1, 2, 3, 4 /

c     ==================================================================

c     ... the output data from qpopt, in vectors of length
c          ndim + mcon, is split into pairs of vectors, of lengths
c          ndim and mcon, respectively.
c
c         in addition, the meaning of the state variables are changed
c         according to the following map
c
c                0 (inactive) ->  0
c                1 (at lower) ->  1
c                2 (at upper) ->  2
c                3 (equality) ->  3
c                4 (ignore  ) ->  4      (should never occur)
c               -2 (infeas.l) -> -1
c               -1 (infeas.u) -> -2
c     
c         we suspect that codes -1 and -2 (infeasible) cannot
c         occur for simple bounds.
c
c     ==================================================================
      
      error = .false.
      
      do i = 1, ndim
         
         if (  istate (i) .lt. -2  .or.  istate (i) .gt. 3 )  then
             error = .true. 
         else
            istatv (i) = sttmap ( istate (i) )
            varmlt (i) = clamda (i)
         endif
         
      enddo

      do i = ndim + 1, ndim + mcon

         if (  istate (i) .lt. -2  .or.  istate (i) .gt. 3 )  then
             error = .true. 
         else
             istatc (i - ndim) = sttmap ( istate (i) )
             conmlt (i - ndim) = clamda (i)
         endif

      enddo

      return
c     end of CUNCAT / cuncat

      end
