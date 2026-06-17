      subroutine xdslec ( icol, jcol, panel, pnlrow, pnlcol )
 
c
      integer              icol  , jcol  , pnlrow, pnlcol

      double precision     panel(pnlcol,pnlrow)
 
      integer              ic    , jc    

      double precision     t
 
c     ------------------------------------------------------------------
c.debug
c     write(6,'("in xdslec - icol, jcol = ", 2i8)') icol, jcol
c.debug

      ic = min ( icol, jcol )
      jc = max ( icol, jcol )

      if ( ic .eq. jc ) return
c.debug
c     write(6,'("in xdslec - ic, jc     = ", 2i8)') ic, jc   
c.debug

c     ------------------------------------------------------------------

c     -------------------------------------------- 
c     ... swap rows ic and jc up to ic-th diagonal
c     -------------------------------------------- 

      call dswap ( ic-1, panel(1,ic), 1, panel(1,jc), 1 )
c.debug
c     write(6,'("after first dswap")')
c.debug

c     ------------------------------------------
c     ... swap ic-th diagonal and jc-th diagonal
c     ------------------------------------------

      t            = panel(ic,ic)
      panel(ic,ic) = panel(jc,jc)
      panel(jc,jc) = t
c.debug
c     write(6,'("after diagonal swap")')
c.debug

c     -------------------------------------------- 
c     ... swap rows ic and jc up to ic-th diagonal
c     -------------------------------------------- 

      call dswap ( jc-ic-1, panel(ic,ic+1), pnlcol, panel(ic+1,jc), 1 )
c.debug
c     write(6,'("after second dswap")')
c.debug

c     -------------------------------
c     ... swap last of rows ic and jc
c     -------------------------------

      call dswap ( pnlrow-jc, panel(ic,jc+1), pnlcol, 
     1                        panel(jc,jc+1), pnlcol )
c.debug
c     write(6,'("after third dswap")')
c.debug

c     ------------------------------------------------------------------
 
      return
      end
