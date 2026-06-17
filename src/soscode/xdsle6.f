      subroutine xdsle6 ( icol, jcol, panel, pnlrow, pnlcol )
c
      integer              icol  , jcol  , pnlrow, pnlcol

      double precision     panel(pnlrow,pnlcol)
 
      integer              ic    , jc    

      double precision     t
 
c     ------------------------------------------------------------------

      ic = min ( icol, jcol )
      jc = max ( icol, jcol )

      if ( ic .eq. jc ) return

c     ------------------------------------------------------------------

c     -------------------------------------------- 
c     ... swap rows ic and jc up to ic-th diagonal
c     -------------------------------------------- 

      call dswap ( ic-1, panel(ic,1), pnlrow, panel(jc,1), pnlrow )

c     ------------------------------------------
c     ... swap ic-th diagonal and jc-th diagonal
c     ------------------------------------------

      t            = panel(ic,ic)
      panel(ic,ic) = panel(jc,jc)
      panel(jc,jc) = t

c     -------------------------------------------- 
c     ... swap rows ic and jc up to ic-th diagonal
c     -------------------------------------------- 

      call dswap ( jc-ic-1, panel(ic+1,ic), 1, panel(jc,ic+1), pnlrow )

c     -------------------------------
c     ... swap last of rows ic and jc
c     -------------------------------

      if(pnlrow-jc.gt.0) then
        call dswap ( pnlrow-jc, panel(jc+1,ic), 1, panel(jc+1,jc), 1 )
      endif
 
c     ------------------------------------------------------------------
 
      return
      end
