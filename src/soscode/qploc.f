      subroutine   QPLOC    ( cset, n, nclin, litotl, lwtotl,
     1                        ldT, ncolT, ldQ,
     2                        miniw, lkactv, lkx,
     3                        minw, lfeatu, lAnorm, lAd, lHx, ld, lgq,
     4                        lcq, lrlam, lR, lT, lQ, lwtinf, lwrk )

c     ==================================================================
c     ==================================================================
c     ====  QPLOC  /                                               =====
c     ====  qploc -- storage allocation for qpopt (QPOPT )         =====
c     ==================================================================
c     ==================================================================

      integer            n, nclin, litotl, lwtotl,
     1                   ldT, ncolT, ldQ,
     2                   miniw, lkactv, lkx,
     3                   minw, lfeatu, lAnorm, lAd, lHx, ld, lgq,
     4                   lcq, lrlam, lR, lT, lQ, lwtinf, lwrk
      
      logical            cset

c     ==================================================================
c     QPLOC  / qploc  allocates the addresses of the work arrays for
c                     qpcore (QPCORE).
c
c     derived from qpopt version 1.0
c     last modification -- 25-March-1996
c
c          Original version written  2-January-1987.
c         This version of qploc dated  18-Nov-1990.
c
c     Note that the arrays ( gq, cq ) lie in contiguous areas of
c     workspace.
c
c     ==================================================================

      integer            lenRT, lenQ, lencq
      
c     ==================================================================

c     ------------------------------------------------------------------
c     Refer to the first free space in the work arrays.
c     ------------------------------------------------------------------

      miniw     = litotl + 1
      minw      = lwtotl + 1

      
c     -----------------------------
c     Integer workspace allocation.
c     -----------------------------

      lkactv    = miniw
      lkx       = lkactv + n
      miniw     = lkx    + n

c     ------------------------------------------------------------------
c     Real workspace.
c     Assign array lengths that depend upon the problem dimensions.
c     ------------------------------------------------------------------

      lenRT     = ldT *ncolT
      if (nclin .eq. 0) then
         lenQ  = 0
      else
         lenQ  = ldQ*ldQ
      end if

      if  ( cset )  then
         lencq  = n
      else
         lencq  = 0
      end if

c     ------------------------------------------------------------------
c     We start with arrays that can be preloaded by smart users.
c     ------------------------------------------------------------------

      lfeatu    = minw
      minw      = lfeatu + nclin + n

c     Next comes stuff used by  lpcore (LPCORE) and  qpcore (QPCORE).

      lAnorm    = minw
      lAd       = lAnorm + nclin
      lHx       = lAd    + nclin
      ld        = lHx    + n
      lgq       = ld     + n
      lcq       = lgq    + n
      lrlam     = lcq    + lencq
      lR        = lrlam  + n
      lT        = lR
      lQ        = lT     + lenRT
      lwtinf    = lQ     + lenQ
      lwrk      = lwtinf + n  + nclin
      minw      = lwrk   + n  + nclin

      litotl    = miniw - 1
      lwtotl    = minw  - 1

c     end of QPLOC  (qploc)
      
      end
