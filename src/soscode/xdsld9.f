      subroutine xdsld9( ncolf , pnlcol, pnlrow, updsiz, panell,
     1                   panelu, pvtblk, loclfr, frontl, frontu,
     2                   temp1 , temp2 , fctops )
 
c
c  purpose -- to apply the factored columns in panel to the remainder
c             of the front - this version has the columns stored
c             as columns in the panel
c
c             unsymmetric version
c
c  created            -- 01-jun-98, rgg, from xdsle9
c  last modifications -- 
c                                        from 4 to updsiz
c
c  input variables --
c
c      ncolf  -- number of columns that were factored
c      pnlcol -- number of columns in panel
c      pnlrow -- number of rows    in panel
c      updsiz -- factorization update panel size
c      panell -- rectangular array holding factored columns
c                lower tri.
c      panelu -- rectangular array holding factored columns
c                upper tri.
c      pvtblk -- indicator of 1x1 and 2x2 pivots
c      loclfr -- local size of the front
c      fctops -- factor operation count
c
c  working storage --
c      
c      temp1  -- temporary array of size ncolf  by updsiz
c      temp2  -- temporary array of size loclfr by updsiz
c
c  output variable --
c
c      frontl -- the remainder of the front - lower tri.
c      frontu -- the remainder of the front - upper tri.
c      fctops -- factor operation count
c
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           ncolf , pnlcol, pnlrow, updsiz, loclfr
 
      integer           pvtblk(*)
 
      double precision  fctops
 
      double precision  panell(pnlrow,pnlcol), 
     1                  panelu(pnlrow,pnlcol), 
     2                  frontl(*), 
     3                  frontu(*), 
     4                  temp1(ncolf,updsiz) , 
     5                  temp2(loclfr,updsiz)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer           icol  , irow  ,
     1                  j     , jcol  , jj    , 
     2                  l     , ncol  , nrow

      double precision  fcol  , ffront, mone  , 
     1                  one   , ops   , zero
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      external          dgemm
 
c  =====================================================================
c.debug
c     write(6,'("in xdsld9 - ncolf, loclfr = ",2i8)') ncolf, loclfr
c.debug
 
      fcol   = ncolf
      ffront = loclfr
 
      ops    = (2*fcol) * ( ffront*(ffront+1) / 2. )
     1       + ffront * fcol

      fctops = fctops + 2 * ops

      one    =  1.0d0
      mone   = -1.0d0
      zero   =  0.0d0
 
c  ----------------------------------------------------------

c     ----------------------------------------------------
c     ... loop over the remainder of the lower tri. front.
c     ----------------------------------------------------

      l      = 1

      do 200 jcol = 1, loclfr, updsiz

          ncol  = min ( updsiz, loclfr - jcol + 1 )

c         ---------------------
c         ... form temp = d * u
c         ---------------------

          do icol = 1, ncol
              jj = jcol + pnlcol + icol - 1
              do j = 1, ncolf
                  temp1(j,icol) = panell(j,j) * panelu(jj,j)
              enddo
          enddo

c         ------------------------------
c         ... main computational loop.
c             form temp2 = panell * temp1
c         ------------------------------

          nrow = loclfr - jcol + 1
          irow = pnlcol + jcol
c.debug
c     write(6,'("before dgemm - nrow, ncol, ncolf, pnlrow = ", 4i8)')
c    1                           nrow, ncol, ncolf, pnlrow 
c     write(6,'("before dgemm - irow                      = ", 4i8)')
c    1                           irow                      
c     call xdslp5 ( 'panell', pnlrow*ncolf, panell, 6 )
c     call xdslp5 ( 'temp1', ncolf*ncol  , temp1, 6 )
c.debug

          call dgemm ( 'n', 'n', nrow, ncol, ncolf,
     1                 one, panell(irow,1), pnlrow,
     2                 temp1, ncolf, zero, temp2, loclfr )
c.debug
c     call xdslp5('after dgemm - temp2', updsiz*loclfr, temp2, 6 )
c.debug

c         --------------------------------------------------
c         ... subtract update from current columns in front.
c         --------------------------------------------------

          do icol = 1, ncol
c.debug
c     write(6,'("in 150 - icol, nrow, l = ", 3i8)') 
c    1                     icol, nrow, l 
c.debug

              do k=0,nrow-1
                frontl(l+k) = frontl(l+k) - temp2(icol+k,icol)
              enddo

              l    = l + nrow
              nrow = nrow - 1

          enddo

  200 continue

c  ----------------------------------------------------------

c     ----------------------------------------------------
c     ... loop over the remainder of the upper tri. front.
c     ----------------------------------------------------

      l      = 1

      do 400 jcol = 1, loclfr, updsiz

          ncol  = min ( updsiz, loclfr - jcol + 1 )

c         ----------------------------
c         ... form temp = d * trans(l)
c         ----------------------------

          do icol = 1, ncol
              jj = jcol + pnlcol + icol - 1
              do j = 1, ncolf
                  temp1(j,icol) = panell(j,j) * panell(jj,j)
              enddo
          enddo

c         ------------------------------
c         ... main computational loop.
c             form temp2 = temp1 * panelu
c         ------------------------------

          nrow = loclfr - jcol + 1
          irow = pnlcol + jcol
c.debug
c     write(6,'("before dgemm - nrow, ncol, ncolf, pnlrow = ", 4i8)')
c    1                           nrow, ncol, ncolf, pnlrow 
c     write(6,'("before dgemm - irow                      = ", 4i8)')
c    1                           irow                      
c     call xdslp5 ( 'panelu', pnlrow*ncolf, panelu, 6 )
c     call xdslp5 ( 'temp1', ncolf*ncol  , temp1, 6 )
c.debug

          call dgemm ( 'n', 'n', nrow, ncol, ncolf,
     1                 one, panelu(irow,1), pnlrow,
     2                 temp1, ncolf, zero, temp2, loclfr )
c.debug
c     call xdslp5('after dgemm - temp2', updsiz*loclfr, temp2, 6 )
c.debug

c         --------------------------------------------------
c         ... subtract update from current columns in front.
c         --------------------------------------------------

          do icol = 1, ncol
c.debug
c     write(6,'("in 350 - icol, nrow, l = ", 3i8)') 
c    1                     icol, nrow, l 
c.debug

              do k=0,nrow-1
                frontu(l+k) = frontu(l+k) - temp2(icol+k,icol)
              enddo

              l    = l + nrow
              nrow = nrow - 1

          enddo

  400 continue

c  =====================================================================
 
      return
      end
