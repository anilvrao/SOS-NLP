      subroutine  xdsleq    ( nelim , Sorder, updsiz, pnlrow, pnlcol,
     1                        panell, panelu, pvtblk, lower , q1by1 , 
     2                        locbfr, ocbufr, iops5u, wafil5, watrn5, 
     3                        temp1 , temp2 , ierr )

 
c
c  purpose -- to apply the rank  nelim  outer product modification from
c             a block of  nelim  factored columns (one panel) to one 
c             triangle of the remainder of the front - this version has
c             the columns stored as columns in the panel
c
c  created            -- 16-apr-03, jgl  from xdsleh and xdsldh
c  last modifications -- 
c
c  input variables --
c
c      nelim  -- number of columns that were factored
c      Sorder -- order of the (2,2) portion of the front (the number 
c                of colums to receive the outer product modification
c      updsiz -- factorization update panel size
c      pnlrow -- number of rows    in panel
c      pnlcol -- number of columns in panel
c      panell -- rectangular array holding factored columns from L
c      panelu -- rectangular array holding factored columns from U
c      pvtblk -- indicator of 1x1 and 2x2 pivots (used only in 
c                the symmetric case)
c      lower  -- logical specifying whether to modify the lower or
c                the upper triangle of the reduced matrix
c      q1by1  -- logical specifying whether all pivots are 1x1
c                (always true in structurally symmetric case)
c      locbfr -- length of ocbufr
c      iops5u -- i/o position for the start of the section of the
c                front to be updated
c      wafil5 -- i/o file holding the front
c
c  working storage --
c      
c      ocbufr -- i/o buffer (must be large enough to hold the
c                longest column of the reduced matrix)
c      temp1  -- temporary array of size nelim  by updsiz
c      temp2  -- temporary array of size Sorder by updsiz
c
c  output variable --
c
c      watrn5 -- i/o transfer count for i/o file wafil5
c      ierr   -- i/o error return
c
c  =====================================================================
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           nelim , pnlcol, pnlrow, Sorder, locbfr,
     1                  iops5u, wafil5, ierr  , updsiz
 
      logical           lower, q1by1

      integer           pvtblk (nelim)
 
      double precision  watrn5
 
      double precision  panell (pnlrow,pnlcol), 
     1                  panelu (pnlrow,pnlcol),
     2                  ocbufr (locbfr),
     3                  temp1 (nelim,updsiz), 
     4                  temp2 (Sorder,updsiz)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer           bfrlen, Bncols, B1stcl, clsfit, iopos , j     , 
     1                  k     , kk    , lenclk, lnBcl1, P1strw, x_Skk  

      integer           idummy(1)

      double precision  diag1 , diag2 , mone  , offdia, one   , zero

      parameter       ( one = 1.0d0, mone = -1.0d0, zero = 0.0d0 )
 
c  =====================================================================

      bfrlen = 0
      iopos  = iops5u

c     --------------------------------------------------------------
c     ... We have determined a block of  nelim  leading columns 
c         (and rows) of an LDL' or LDU factorization of a dense
c         matrix.  We now perform a rank  nelim  outer product to 
c         compute the remaining columns of the lower triangle of the 
c         reduced matrix (or the remaining rows of the upper triangle
c         of the reduced matrix).  We do so in a series of block 
c         operationsa, in which the triangle of the reduced matrix
c         is modified in blocks of at most  updsiz  columns (rows
c         for the upper triangle).  
c         The outer product has the form L * D * L' or L D U,
c         where  L  denotes the block of columns just
c         eliminated (and U the block of rows).  
c
c         We form this product as
c              temp1 = (D * L')     (symmetric case)
c              temp2 = L * temp1
c         or
c              temp1 = (D * U)      (lower triangle, unsymmetric case)
c              temp2 = L * temp1
c         or
c              temp1 = (D * L')     (upper triangle, unsymmetric case)
c              temp2 = U' * temp1
c         in order to use matrix-matrix multiply as the main kernel.
c
c         In all of these, the matrix temp1 is square, of order nelim.
c         The last of these follows naturally because we store the
c         upper triangle U by rows rather than by columns.  Thus U'
c         is available with stride one.
c
c         The symmetric case follows for the same reason by supplying
c         L as the argument for both L and U.  The diagonal is 
c         stored in L.
c     --------------------------------------------------------------

      x_Skk = 1

      do 200 B1stcl = 1, Sorder, updsiz

          Bncols  = min ( updsiz, Sorder - B1stcl + 1 )
          lnBcl1 = Sorder - B1stcl + 1

c         ---------------------------------
c         ... form temp1 = D * L'  or D * U
c         ---------------------------------

          if ( q1by1 .and. lower ) then

c            --------------------------------
c            ... D  has only 1x1 pivot blocks
c                symmetric or general lower
c                triangle case
c            --------------------------------

              do k = 1, Bncols
                  kk = B1stcl + pnlcol + k - 1
                  do j = 1, nelim
                      temp1(j,k) = panell(j,j) * panelu(kk,j)
                  enddo
              enddo

          elseif ( q1by1 .and. .not. lower ) then

c            --------------------------------
c            ... D  has only 1x1 pivot blocks
c                upper triangle case
c            --------------------------------

              do k = 1, Bncols
                  kk = B1stcl + pnlcol + k - 1
                  do j = 1, nelim
                      temp1(j,k) = panell(j,j) * panell(kk,j)
                  enddo
              enddo

          else

c            -----------------------------------------
c            ... D has some 2x2 pivot blocks
c                (applies only in the symmetric case).
c            -----------------------------------------

            j = 1

   70       continue
            if ( j .le. nelim ) then

              if ( pvtblk(j) .eq. 1 ) then
     
                  do k = 1, Bncols
                      kk         = B1stcl + pnlcol + k - 1
                      temp1(j,k) = panell(j,j) * panell(kk,j)
                  enddo

                  j = j + 1

              else

                  diag1  = panell(j,j)
                  offdia = panell(j+1,j)
                  diag2  = panell(j+1,j+1)
     
                  do k = 1, Bncols
                      kk = B1stcl + pnlcol + k - 1
                      temp1(j  ,k) = diag1  * panell(kk,j  ) +
     1                               offdia * panell(kk,j+1)
                      temp1(j+1,k) = offdia * panell(kk,j  ) +
     1                               diag2  * panell(kk,j+1)
                  enddo

                  j = j + 2

              end if

              go to 70

            endif

          end if
c.debug
c     write(6,'("Bncols, nelim = ", 2i8)') Bncols, nelim
c     call xdslp5('after 100 - temp1', Bncols*nelim, temp1, 6 )
c.debug

c         ------------------------------------------
c         ... main computational loop.  
c             form temp2 = L * temp1  or  U' * temp1
c         ------------------------------------------

          P1strw = pnlcol + B1stcl 

c.debug
c     write(6,'("in 200 - x_Skk, P1strw           = ", 4i8)') 
c    1                     x_Skk, P1strw          
c.debug

c.debug
c     write(6,
c    1  '("before dgemm - lnBcl1, Bncols, nelim, pnlrow = ", 4i8)')
c    1                           lnBcl1, Bncols, nelim, pnlrow 
c     write(6,'("in 200 - x_Skk                 = ", 4i8)') 
c    1                     x_Skk                
c     call xdslp5 ( 'panel', pnlrow*nelim, panel, 6 )
c     call xdslp5 ( 'temp1', nelim*Bncols  , temp1, 6 )
c.debug

          if ( lower ) then

             call dgemm ( 'n', 'n', lnBcl1, Bncols, nelim,
     1                    one, panell(P1strw,1), pnlrow,
     2                    temp1, nelim, zero, temp2, Sorder )

          else

             call dgemm ( 'n', 'n', lnBcl1, Bncols, nelim,
     1                    one, panelu(P1strw,1), pnlrow,
     2                    temp1, nelim, zero, temp2, Sorder )

          endif

c.debug
c     call xdslp5('after dgemm - temp2', Bncols*Sorder, temp2, 6 )
c.debug

c         -----------------------------------------------------------
c         ... subtract the low rank modification as a block from 
c             the next  Bncols  columns in the front.  
c             These columns are on secondary storage, and are read in
c             order into the buffer  ocbufr.
c         ------------------------------------------------------------

          do 150 k = 1, Bncols

c            --------------------------------------------------
c            ... check to see if this column is already in
c                memory.  if not, write the current contents of
c                the secondary storage buffer, and refill with
c                as many columns as will fit.
c            --------------------------------------------------

             lenclk = lnBcl1 - k + 1

             if ( x_Skk + lenclk  .gt.  bfrlen + 1 ) then

c               ----------------------------------------------
c               ... time to read in next section of the front.
c                   first dump current section to i/o file.
c               ----------------------------------------------

                if ( x_Skk .gt. 1 ) then 
                   
                   call xdslw2 ( wafil5, 2, idummy, idummy, ocbufr,
     1                           iopos, x_Skk-1, ierr )
                   if ( ierr .ne. 0 ) then
                      ierr = -2
                      return
                   end if

                   iopos  = iopos  + x_Skk - 1
                   watrn5 = watrn5 + x_Skk - 1

                end if

c               ---------------------------------------
c               ... read in as many columns as will fit
c                   note that we may bring in columns 
c                   from the next block now.
c               ---------------------------------------

                call xislex ( lenclk, locbfr, clsfit, bfrlen )
                if  ( bfrlen .le. 0 )  then
                   ierr = -1
                   return
                end if
                
                call xdslw1 ( wafil5, 2, idummy, idummy, ocbufr,
     1                        iopos, bfrlen, ierr )
                if ( ierr .ne. 0 ) then
                   ierr = -3
                   return
                end if

                watrn5 = watrn5 + bfrlen

                x_Skk  = 1

             end if

             do kk=0,lenclk-1
               ocbufr(x_Skk+kk) = ocbufr(x_Skk+kk) - temp2(k+kk,k)
             enddo

             x_Skk  = x_Skk + lenclk
             lenclk = lenclk - 1

 150      continue

 200   continue


c     ----------------------------------------
c     ... finish up by writing last section of 
c         one triangle of front to i/o file
c     ----------------------------------------

      if ( x_Skk - 1 .ne. bfrlen ) then 

         ierr = -1
         return

      else
                  
         call xdslw2 ( wafil5, 2, idummy, idummy, ocbufr,
     1                 iopos, bfrlen, ierr )
         
         if ( ierr .ne. 0 ) then
            ierr = -2
            return
         end if          

         watrn5 = watrn5 + bfrlen

      end if

      return

      end
