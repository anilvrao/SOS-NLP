      subroutine xdsle9( ncolf , pnlcol, pnlrow, updsiz, panel ,
     1                   pvtblk, loclfr, front , temp1 , temp2 , 
     2                   fctops )
 
c
c  purpose -- to apply the factored columns in panel to the remainder
c             of the front - this version has the columns stored
c             as columns in the panel
c
c  created            -- 07-feb-97, rgg
c  last modifications -- 09-mar-98, rgg, increased update panel size
c                                        from 4 to updsiz
c
c  input variables --
c
c      ncolf  -- number of columns that were factored
c      pnlcol -- number of columns in panel
c      pnlrow -- number of rows    in panel
c      updsiz -- factorization update panel size
c      panel  -- rectangular array holding factored columns
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
c      front  -- the remainder of the front
c      fctops -- factor operation count
c
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           ncolf , pnlcol, pnlrow, updsiz, loclfr
 
      integer           pvtblk(*)
 
      double precision  fctops
 
      double precision  panel(pnlrow,pnlcol), front(*), 
     1                  temp1(ncolf,updsiz) , temp2(loclfr,updsiz)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer           i     , icol  , irow  ,
     1                  j     , jcol  , jj    , 
     2                  l     , ncol  , nrow

      logical           q1by1

      double precision  diag1 , diag2 , fcol  , ffront, mone  , 
     1                  offdia, one   , ops   , zero
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      external          dgemm
 
c  =====================================================================
c.debug
c     write(6,'("in xdsle9 - ncolf, loclfr = ",2i8)') ncolf, loclfr
c.debug
 
      fcol   = ncolf
      ffront = loclfr
 
      ops    = (2*fcol) * ( ffront*(ffront+1) / 2. )
     1       + ffront * fcol

      fctops = fctops + ops

      one    =  1.0d0
      mone   = -1.0d0
      zero   =  0.0d0

c     --------------------------------
c     ... set 1by1 flag as appropriate
c     --------------------------------

      q1by1  = .true.

      do i = 1, ncolf
          if ( pvtblk(i) .eq. 2 ) then
              q1by1 = .false.
              exit
          end if
      enddo
c.debug
c     write(6,'("q1by1 = ", l5)') q1by1
c.debug
 
c  ----------------------------------------------------------

c     -----------------------------------------
c     ... loop over the remainder of the front.
c     -----------------------------------------

      l      = 1

      do 200 jcol = 1, loclfr, updsiz

          ncol  = min ( updsiz, loclfr - jcol + 1 )

c         ----------------------------
c         ... form temp = d * trans(l)
c         ----------------------------

          if ( q1by1 ) then

              do icol = 1, ncol
                  jj = jcol + pnlcol + icol - 1
                  do j = 1, ncolf
                      temp1(j,icol) = panel(j,j) * panel(jj,j)
                  enddo
              enddo

          else

            j = 1

   70       continue
            if ( j .le. ncolf ) then

              if ( pvtblk(j) .eq. 1 ) then
     
                  do icol = 1, ncol
                      jj = jcol + pnlcol + icol - 1
                      temp1(j,icol) = panel(j,j) * panel(jj,j)
                  enddo

                  j = j + 1

              else

                  diag1  = panel(j,j)
                  offdia = panel(j+1,j)
                  diag2  = panel(j+1,j+1)
     
                  do icol = 1, ncol
                      jj = jcol + pnlcol + icol - 1
                      temp1(j  ,icol) = diag1  * panel(jj,j  )
     1                                + offdia * panel(jj,j+1)
                      temp1(j+1,icol) = offdia * panel(jj,j  )
     1                                + diag2  * panel(jj,j+1)
                  enddo

                  j = j + 2

              end if

              go to 70

            endif

          end if
c.debug
c     write(6,'("ncol, ncolf = ", 2i8)') ncol, ncolf
c     call xdslp5('after 100 - temp1', ncol*ncolf, temp1, 6 )
c.debug

c         ------------------------------
c         ... main computational loop.
c             form temp2 = panel * temp1
c         ------------------------------

          nrow = loclfr - jcol + 1
          irow = pnlcol + jcol
c.debug
c     write(6,'("before dgemm - nrow, ncol, ncolf, pnlrow = ", 4i8)')
c    1                           nrow, ncol, ncolf, pnlrow 
c     write(6,'("before dgemm - irow                      = ", 4i8)')
c    1                           irow                      
c     call xdslp5 ( 'panel', pnlrow*ncolf, panel, 6 )
c     call xdslp5 ( 'temp1', ncolf*ncol  , temp1, 6 )
c.debug

          call dgemm ( 'n', 'n', nrow, ncol, ncolf,
     1                 one, panel(irow,1), pnlrow,
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

              do kk=0,nrow-1
                front(l+kk) = front(l+kk) - temp2(icol+kk,icol)
              enddo

              l    = l + nrow
              nrow = nrow - 1

          enddo

  200 continue

c  =====================================================================

      return
      end
