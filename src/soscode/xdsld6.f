      subroutine xdsld6 ( icol, jcol, temp, panell, panelu, pnlrow,
     1                    pnlcol )
 
c
      integer              icol  , jcol  , pnlrow, pnlcol

      double precision     temp(*)       , panell(pnlrow,pnlcol),
     1                     panelu(pnlrow,pnlcol)
 
      integer              ic    , jc    , k
 
      double precision     t
 
c     ------------------------------------------------------------------

      ic = min ( icol, jcol )
      jc = max ( icol, jcol )

      if ( ic .eq. jc ) then
          k    = pnlrow - ic + 1
          panell(ic:ic+k-1,ic) = temp(ic:ic+k-1)
          return
      end if

c.debug
c     write(6,'("in xdsld6 - ic, jc = ", 2i8)') ic, jc
c     call xdslp5('temp1                    ', pnlrow       , temp  ,6)
c     call xdslp5('panell before interchange', pnlcol*pnlrow, panell,6)
c     call xdslp5('panelu before interchange', pnlcol*pnlrow, panelu,6)
c.debug

c     ------------------------------------------------------------------

c     ----------------------------------------------------------------
c     ... interchange columns ic and jc from diagonal to end of matrix
c     ----------------------------------------------------------------

c     ----------------------------------------------------------
c     ... first interchange rows 1 to ic-1 of columns ic and jc
c     ----------------------------------------------------------

      do k = 1, ic-1
          t            = panelu(jc,k)
          panelu(jc,k) = panelu(ic,k)
          panelu(ic,k) = t        
      enddo

c     ----------------------------------------------------------
c     ... now interchange rows ic to jc-1 of columns ic and jc
c         note that temp1 holds column jc of l.
c     ----------------------------------------------------------

      do k = ic, jc-1
          panelu(jc,k) = panell(k,ic)
          panell(k,ic) = temp(k)
      enddo

c     -----------------------------------------------------------
c     ... now interchange rows jc to pnlrow of columns ic and jc.
c     -----------------------------------------------------------

      do k = jc, pnlrow
          panell(k,jc) = panell(k,ic)
          panell(k,ic) = temp(k)      
      enddo

c     ------------------------------------------------------------------

c     ----------------------------------------------------------------
c     ... interchange rows ic and jc from diagonal to end of matrix
c     ----------------------------------------------------------------

c     --------------------------------------------------------
c     ... first interchange columns 1 to ic of rows ic and jc.
c     --------------------------------------------------------

      do k = 1, ic
          t            = panell(jc,k)
          panell(jc,k) = panell(ic,k)
          panell(ic,k) = t
      enddo

c     ---------------------------------------------------------
c     ... now interchange columns ic+1 to jc of rows ic and jc.
c     ---------------------------------------------------------

      do k = ic+1, jc
          t            = panell(jc,k) 
          panell(jc,k) = panelu(k,ic)
          panelu(k,ic) = t
      enddo

c     ------------------------------------------------------------
c     ... now interchange columns jc+1 to pnlrow of row ic and jc.
c     ------------------------------------------------------------

      do k = jc+1, pnlrow
          t            = panelu(k,jc) 
          panelu(k,jc) = panelu(k,ic)
          panelu(k,ic) = t
      enddo
c.debug
c     call xdslp5('panell after interchange', pnlcol*pnlrow, panell, 6 )
c     call xdslp5('panelu after interchange', pnlcol*pnlrow, panelu, 6 )
c     stop
c.debug

c     ------------------------------------------------------------------
 
      return
      end
