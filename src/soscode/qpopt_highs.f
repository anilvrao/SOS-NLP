      subroutine   QPOPT     ( n, nclin, ldA, ldH, nHess,
     1                        A, bl, bu, cvec, H,
     2                        qpHess, istate, x,
     3                        prbtyp, cset, minsum, start,
     4                        inform, obj,
     5                        msglvl, iPrint, iSumm,
     6                        Ax, clamda, iw, leniw, w, lenw,
     7                        mxfree, maxact, maxnZ,
     8                        bigbnd, bigdx, tlcrsh, tolfea, tolOpt,
     9                        tolrnk, tollev,
     A                        fealim, optlim, kdegen, kcheck,
     B                        fpiter, qpiter, nouter, nlvmod, delta,
     C                        p1wkmn, p1wkmx, p1wkav,
     D                        p2wkmn, p2wkmx, p2wkav,
     E                        lminzh, lmaxzh, sminaa, smaxaa )
c     Experimental QPOPT-compatible entry point backed by HiGHS.
c     This preserves the SOS Fortran call signature while moving the
c     dense convex QP solve into an MIT-licensed external solver.

      integer            n, nclin, ldA, ldH, nHess, inform
      integer            msglvl, iPrint, iSumm, leniw, lenw
      integer            mxfree, maxact, maxnZ
      integer            fealim, optlim, kdegen, kcheck
      integer            fpiter, qpiter, nouter, nlvmod
      integer            p1wkmn, p1wkmx, p2wkmn, p2wkmx
      integer            istate ( n+nclin ), iw ( leniw )

      logical            cset, minsum

      character*(*)      prbtyp, start

      double precision   A( ldA, * ), H( ldH, * )
      double precision   bl( n+nclin ), bu( n+nclin )
      double precision   cvec( * ), x( n ), Ax( * ), clamda( n+nclin )
      double precision   w( lenw )
      double precision   obj, bigbnd, bigdx, tlcrsh, tolfea, tolOpt
      double precision   tolrnk, tollev, delta, p1wkav, p2wkav
      double precision   lminzh, lmaxzh, sminaa, smaxaa

      external           qpHess

      fpiter = 0
      qpiter = 0
      nouter = 0
      nlvmod = 0
      delta  = 0.0d0
      p1wkmn = 0
      p1wkmx = 0
      p2wkmn = 0
      p2wkmx = 0
      p1wkav = 0.0d0
      p2wkav = 0.0d0
      lminzh = -1.0d0
      lmaxzh = -1.0d0
      sminaa = -1.0d0
      smaxaa = -1.0d0

      call SOS_HIGHS_QP ( n, nclin, ldA, ldH, nHess,
     1                    A, bl, bu, cvec, H, istate, x,
     2                    inform, obj, msglvl, Ax, clamda,
     3                    bigbnd, tolfea, tolOpt )

      return
      end
