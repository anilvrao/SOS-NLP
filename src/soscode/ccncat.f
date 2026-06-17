      subroutine   CCNCAT   ( ndim, mcon, istatv, istatc, istate,
     1                        xlwr, blwr, lower, xupr, bupr, upper,
     2                        error )
c     ==================================================================
c     ==================================================================
c     ====  CCNCAT /                                                ====
c     ====  ccncat -- interface bcs and qpopt data conventions      ====
c     ==================================================================
c     ==================================================================

      integer            ndim, mcon

      logical            error

      integer            istatv (ndim), istatc (*),
     1                   istate (ndim+mcon)

      double precision   xlwr (ndim), blwr (*), lower (ndim+mcon),
     1                   xupr (ndim), bupr (*), upper (ndim+mcon)

c     ================================================================
c	... last modified 03-sept-1996
c     ================================================================

      integer            i

      integer            sttmap (-2:4)

      data               sttmap / -1, -2, 0, 1, 2, 3, 4 /

c     ==================================================================

c     ... the input data in pairs of vectors, of lengths ndim and mcon,
c         respectively, are concatenated into output vectors of length
c         ndim + mcon.
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

         if (istatv (i) .lt. -2  .or.  istatv (i) .gt. 3 )  then
            error            = .true.
            istate (i+ ndim) = 0
         else
            istate (i) = sttmap ( istatv (i) )
            lower  (i) = xlwr   (i)
            upper  (i) = xupr   (i)
         endif
         
      enddo

      do i = 1, mcon

         if (istatc (i) .lt. -2  .or.  istatc (i) .gt. 3 )  then
            error            = .true.
            istate (i+ ndim) = 0
         else
             istate (i + ndim) = sttmap ( istatc (i) )
             lower  (i + ndim) = blwr   (i)
             upper  (i + ndim) = bupr   (i)
          endif

      enddo

      return

c     end of CCNCAT / ccncat

      end
