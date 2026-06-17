      subroutine xdslb7 ( nn1   , n2    , n3    , loclfr, lpan  , 
     1                    colbgn, colend, panell, panelu,
     2                    nstack, nassmb, stknod, istack, lindxl,
     3                    wafil2, watrn2, locbfr, ocbufr, ierr  )
 
c  =====================================================================

c     --------------------
c     ... global variables
c     --------------------

      integer                 nn1   , n2    , n3    , loclfr, lpan  ,
     1                        colbgn, colend, nstack, nassmb, wafil2,
     2                        locbfr, ierr

      integer                 stknod(*),      istack(4,*),
     1                        lindxl(*)

      double precision        watrn2

      double precision        panell(*),      panelu(*),
     1                        ocbufr(*)

c     -------------------
c     ... local variables
c     -------------------
 
      integer                  bfrlen, bfrpos, i     , iopos2, 
     1                         j     , jcol  , jstack, k     , l     ,
     2                         lfront, lndxpt, lpano , lson  , lsonp , 
     3                         pnlcol, pnlrow, son   , space , totlen,
     4                         pnlrsv

      integer                  idummy(1)

c     ----------------------
c     ... statement function
c     ----------------------

      integer                  pnlind

      pnlind(i,j) = ( lfront*(j-1) ) - ( j*(j-1)/2 ) + i

c  =====================================================================

c     -----------------------------------------------
c     ... expand panel to allow for postponed columns
c     -----------------------------------------------

      space = ( colend - colbgn + 1 ) * n2
      lpano = lpan - space

      k     = lpano + 1	
      l     = loclfr - ( colend - colbgn ) - n3

      do jcol = colend, colbgn, -1

c         ----------------------------------------------------
c         ... slide n3 locations for column jcol down by space
c             amount
c         ----------------------------------------------------

          k = k - n3
c.debug
c     write(6,'("in xdslb7 loop 100")')
c     write(6,'("colbgn, colend, loclfr, lpan  , space   = ", 5i8)')
c    1            colbgn, colend, loclfr, lpan  , space
c     write(6,'("n2    , n3    , k     , l     , jcol    = ", 5i8)')
c    1            n2    , n3    , k     , l     , jcol
c.debug

          call xdslmv ( n3, panell, k, k+space )
          call xdslmv ( n3, panelu, k, k+space )

c         -------------------------
c         ... zero out n2 locations
c         -------------------------

          space = space - n2
c.debug
c     write(6,'("before zeroing out pad - k+space, n2 = ", 2i8)')
c    1                                     k+space, n2
c.debug
 
          panell(k+space:k+space+n2-1) = 0.d0
          panelu(k+space:k+space+n2-1) = 0.d0

c         ----------------------------------
c         ... slide down remainder of column
c         ----------------------------------

          k = k - l
c.debug
c     write(6,'("before slide           - k, space, l = ", 3i8)')
c    1                                     k, space, l
c.debug
          
          call xdslmv ( l, panell, k, k+space )
          call xdslmv ( l, panelu, k, k+space )

          l = l + 1

      enddo

      if ( k .ne. 1 .or. space .ne. 0 ) then
          write(6,'("in xdslb7 - oops after 100")')
          write(6,'("k, space = ", 2i8)') k, space
          stop
      end if
 
c  =====================================================================

c     ----------------------------------------------------
c     ... now put information from postponed columns into 
c         expanded space
c         note:  xdslb7 is only called when processing the
c                first n1 columns and there are postponed
c                columns.  the 230 loop is only adding the
c                information into the rows associated with
c                the postponed rows/columns.  this allows 
c                the do 210 loop to skip over the diagonal
c                block of the postponed columns.  see
c                xdsla8 for code that handles the last
c                n2 + n3 columns.
c     ----------------------------------------------------

      pnlrow = nn1
      lfront = loclfr + n2

      do 300 jstack = 1, nassmb
 
          son    = stknod(nstack)
c.debug
c     write(6,'("inside do 300 loop - jstack, nstack, son = ", 3i8)') 
c    1                                 jstack, nstack, son
c.debug

c         ---------------------------------------------------------
c         ... get info for postponed columns for this update matrix
c         ---------------------------------------------------------
 
          lson   = istack (1, nstack)
          lndxpt = istack (2, nstack) - 1
          lsonp  = istack (3, nstack)
          iopos2 = istack (4, nstack) 
 
c.debug
c     write(6,'("lson, lsonp, lndxpt = ", 3i8)') 
c    1            lson, lsonp, lndxpt 
c.debug

          if ( lsonp .eq. 0 ) go to 290

c         ---------------------------------------
c         ... set up for processing these columns
c         ---------------------------------------

          totlen = lsonp + lsonp*(lsonp+1)/2 + lsonp*lson
          bfrlen = 0
          bfrpos = lsonp + 1

c         ----------------------------          
c         ... add in postponed columns
c         ----------------------------          

          pnlrsv = pnlrow

          do 210 jcol = 1, lsonp

              pnlrow = pnlrow + 1

              bfrpos = bfrpos + ( lsonp - jcol + 1 )

              do j = 1, lson

                  pnlcol = lindxl(lndxpt+j)
c.debug
c     write(6,'("jcol  , pnlrow, bfrpos, j     , pnlcol   = ", 5i8)') 
c    1            jcol  , pnlrow, bfrpos, j     , pnlcol          
c     write(6,'("colbgn, colend                           = ", 5i8)') 
c    1            colbgn, colend                                  
c.debug

                  if ( pnlcol .ge. colbgn  .and.
     1                 pnlcol .le. colend ) then

                      if ( bfrpos .gt. bfrlen ) then

                          k      = ( bfrpos - 1 ) - bfrlen
                          iopos2 = iopos2 + k 
                          totlen = totlen - k
                          bfrlen = min ( locbfr, totlen )
                          bfrpos = 1
c.debug
c     write(6,'("reading from wafil2")') 
c     write(6,'("iopos2, bfrlen, totlen, locbfr           = ", 5i8)') 
c    1            iopos2, bfrlen, totlen, locbfr                  
c.debug

                          call xdslw1 ( wafil2, 2, idummy, idummy, 
     1                                  ocbufr, iopos2, bfrlen, ierr)
                          if ( ierr .ne. 0 ) go to 8000

                          watrn2 = watrn2 + bfrlen
 
                          iopos2 = iopos2 + bfrlen
                          totlen = totlen - bfrlen
c.debug
c     call xdslp5 ( 'ocbufr', bfrlen, ocbufr, 6 )
c.debug

                      end if

                      k = pnlind(pnlrow,pnlcol-colbgn+1)

                      if ( pnlrow .le. pnlcol ) then
                          panell(k) = ocbufr(bfrpos)
                      else
                          panelu(k) = ocbufr(bfrpos)
                      end if
c.debug
c     write(6,'("pnlrow, pnlcol, k     , bfrpos           = ", 5i8)') 
c    1            pnlrow, pnlcol, k     , bfrpos                  
c     write(6,'("ocbufr(bfrpos), panell(k), panelu(k)  = ", 3d25.15)') 
c    1            ocbufr(bfrpos), panell(k), panelu(k)
c.debug

                  end if

                  bfrpos = bfrpos + 1

              enddo

  210     continue 
c.debug
c     write(6,'("after adding in postponed columns in xdslb7")')
c     l = lpano + ( colend - colbgn + 1 ) * n2
c     call xdslp5 ( 'panell', l, panell, 6 )
c.debug

c         -----------------------------
c         ... now handle upper triangle
c         -----------------------------

          iopos2 = istack (4, nstack) 
     1           + lsonp + lsonp*(lsonp+1)/2 + lsonp*lson
     2           + lson*(lson+1)/2 
          totlen = lsonp*(lsonp+1)/2 + lsonp*lson

          bfrlen = 0
          bfrpos = 1

c         ----------------------------          
c         ... add in postponed columns
c         ----------------------------          

          pnlrow = pnlrsv

          do 230 jcol = 1, lsonp

              pnlrow = pnlrow + 1

              bfrpos = bfrpos + ( lsonp - jcol + 1 )

              do j = 1, lson

                  pnlcol = lindxl(lndxpt+j)
c.debug
c     write(6,'("jcol  , pnlrow, bfrpos, j     , pnlcol   = ", 5i8)') 
c    1            jcol  , pnlrow, bfrpos, j     , pnlcol          
c     write(6,'("colbgn, colend                           = ", 5i8)') 
c    1            colbgn, colend                                  
c.debug

                  if ( pnlcol .ge. colbgn  .and.
     1                 pnlcol .le. colend ) then

                      if ( bfrpos .gt. bfrlen ) then

                          k      = ( bfrpos - 1 ) - bfrlen
                          iopos2 = iopos2 + k 
                          totlen = totlen - k
                          bfrlen = min ( locbfr, totlen )
                          bfrpos = 1
c.debug
c     write(6,'("reading from wafil2")') 
c     write(6,'("iopos2, bfrlen, totlen, locbfr           = ", 5i8)') 
c    1            iopos2, bfrlen, totlen, locbfr                  
c.debug

                          call xdslw1 ( wafil2, 2, idummy, idummy, 
     1                                  ocbufr, iopos2, bfrlen, ierr)
                          if ( ierr .ne. 0 ) go to 8000

                          watrn2 = watrn2 + bfrlen
 
                          iopos2 = iopos2 + bfrlen
                          totlen = totlen - bfrlen
c.debug
c     call xdslp5 ( 'ocbufr', bfrlen, ocbufr, 6 )
c.debug

                      end if

                      k = pnlind(pnlrow,pnlcol-colbgn+1)

                      if ( pnlrow .le. pnlcol ) then
                          panelu(k) = ocbufr(bfrpos)
                      else
                          panell(k) = ocbufr(bfrpos)
                      end if
c.debug
c     write(6,'("pnlrow, pnlcol, k     , bfrpos           = ", 5i8)') 
c    1            pnlrow, pnlcol, k     , bfrpos                  
c     write(6,'("ocbufr(bfrpos), panell(k), panelu(k)  = ", 3d25.15)') 
c    1            ocbufr(bfrpos), panell(k), panelu(k)
c.debug

                  end if

                  bfrpos = bfrpos + 1

              enddo

  230     continue 
c.debug
c     write(6,'("after adding in postponed columns in xdslb7")')
c     l = lpano + ( colend - colbgn + 1 ) * n2
c     call xdslp5 ( 'panelu', l, panelu, 6 )
c.debug

c         ---------------------------------
c         ... set up for next update matrix
c         ---------------------------------

  290     continue
          nstack = nstack - 1

  300 continue

      return
 
c  =====================================================================
 
c     --------------------------------------
c     ... error trap for i/o error on wafil2
c     --------------------------------------
 
 8000 continue
      ierr = -1
      return

c  =====================================================================
 
      end
