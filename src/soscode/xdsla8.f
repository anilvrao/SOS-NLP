      subroutine xdsla8 ( n1    , n2    , n3    , loclfr, lpan  , 
     1                    colbgn, colend, panel ,
     2                    nstack, nassmb, stknod, istack, lindxl,
     3                    wafil2, watrn2, locbfr, ocbufr, 
     4                    actual, sup   , ierr  , lstchc, midchl )

c
c  purpose -- add postponed columns into expanded frontal matrix.
c             symmetric out-of-core version.
c
c  created   -- unknown  , rgg
c  revisions -- 14-apr-03, jgl  retained status of child processing
c                               so that code can properly handle
c                               multiple children in multiple
c                               buffer loads 
c
c  variables
c
c     n1     -- number of expected nodes in supernode (number of
c               fully assembled columns)
c     n2     -- number of postponed columns being added
c     n3     -- dimension of Schur complement / modification block
c     loclfr -- order of some partial front matrix
c     lpan   -- active length of panel array
c     colbgn -- index of first column in front of active panel
c     colend -- index of last  column in front of active panel
c     panel  -- a block of columns of the symmetric packed frontal
c               matrix, corresponding in this case to some or all
c               of the columns postponed to this elimination step
c     nstack -- number of modification matrices on stack.  Also
c               index of top of stack
c     nassmb -- number of children on this supernode
c     stknod -- original supernode number of children
c     istack -- integerscalar information for each update matrix in
c                the stack
c     lindxl -- direct assembly indices for expected front
c     wafil2 -- direct access file for stack
c     watrn2 -- i/o transfer count for wafil2
c     locbfr -- length of Out of core buffer
c     ocbufr -- out of core buffer
c     actual -- pointer to next unfilled slot in old to new
c               stability permutation ("sup")
c     sup    -- old to new permutation vector for interchanges
c               resulting from stability pivoting
c     ierr   -- 
c     lstchc -- last child whose postponed columns have been added
c               completely into the front.  (lstchc+1 is the child
c               either about to be or already in process of adding
c               postponed columns into front.)
c     midchl -- am or am not partly through merging the active child
 
c  =====================================================================

c     --------------------
c     ... global variables
c     --------------------

      integer                 n1    , n2    , n3    , loclfr, lpan  ,
     1                        colbgn, colend, nstack, nassmb, wafil2,
     2                        locbfr, actual, ierr

      integer                 lstchc

      logical                 midchl

      integer                 stknod(*),      istack(4,*),
     1                        lindxl(*),      sup(*)

      double precision        watrn2

      double precision        panel(*),       ocbufr(*)

c     -------------------
c     ... local variables
c     -------------------
 
      integer                  bfrlen, bfrpos, i     , iopos2,
     1                         j     , jcol  , jstack, k     , 
     2                         l     , lfront, lndxpt, lson  , lsonp ,
     3                         nn2   , pnlcol, pnloff, pnlrow,
     4                         son   , totlen

      integer                  clsfit, fstlen

      integer                  idummy(1)

c     ----------------------
c     ... statement function
c     ----------------------

      integer                  pnlind

      pnlind(i,j) = ( lfront*(j-colbgn) ) 
     1            - ( (j-colbgn)*(j-colbgn-1)/2 ) 
     2            + (i-j+1)

c  =====================================================================

c     -----------------------------------------------
c     ... put information from postponed columns into 
c         columns reserved for them.
c     -----------------------------------------------

      pnlcol = n1
      nn2    = n2
      lfront = loclfr + n2 - ( colbgn - ( n1 + 1 ) ) 

      do jstack = 1, lstchc
         nn2    = nn2 - istack (3, nstack)
         pnlcol = pnlcol + istack (3, nstack)
         nstack = nstack - 1
      enddo

      do 230 jstack = lstchc + 1, nassmb
 
          son    = stknod(nstack)
c.debug
c     write(6,'("inside do 230 loop - jstack, nstack, son = ", 3i8)') 
c    1                                 jstack, nstack, son
c     write(6,'("                     lfront              = ", 3i8)') 
c    1                                 lfront
c.debug

c         ---------------------------------------------------------
c         ... get info for postponed columns for this update matrix
c         ---------------------------------------------------------
 
          lson   = istack (1, nstack)
          lndxpt = istack (2, nstack) - 1
          lsonp  = istack (3, nstack)
          iopos2 = istack (4, nstack) 
c.debug
c                write(6,'("lson, lsonp, lndxpt, iopos2 = ", 4i8)') 
c    1                       lson, lsonp, lndxpt, iopos2 
c.debug

          if ( lsonp .eq. 0 ) then
             lstchc = jstack
             go to 220
          end if

          totlen = lsonp + lsonp*(lsonp+1)/2 + lsonp*lson
          if ( .not. midchl ) then

c             --------------------------------------
c             ... read in first block of information
c             --------------------------------------
             
              clsfit = colend - pnlcol
              if ( clsfit .lt. lsonp ) then
                  fstlen = lsonp + clsfit * (lsonp + lson) 
     1                           - (clsfit * (clsfit-1)) / 2
                  bfrlen = min ( locbfr, fstlen )
              else
                  bfrlen = min ( locbfr, totlen )
              end if
c.debug
c     write(6,'("totlen, bfrlen              = ", 4i8)') 
c    1            totlen, bfrlen 
c.debug

              call xdslw1 ( wafil2, 2, idummy, idummy, ocbufr, 
     1                      iopos2, bfrlen, ierr)
              if ( ierr .ne. 0 ) go to 8000
c.debug
c     call xislp3 ( 'ocbufr-integers', lsonp, ocbufr, 6 )
c     call xdslp5 ( 'ocbufr-floats', bfrlen-lsonp, ocbufr(lsonp+1), 6 )
c.debug
 
              iopos2 = iopos2 + bfrlen
              totlen = totlen - bfrlen
              bfrpos = lsonp + 1
              watrn2 = watrn2 + bfrlen

              call icopy ( lsonp, ocbufr, 1, sup(actual), 1 )
c.debug
c     write(6,'("copying to sup - lsonp, actual = ", 2i8)') 
c    1            lsonp, actual  
c                call xislp3 ( 'new perm values', lsonp, 
c    1                          sup (actual), 6 )
c     call xislp3 ( 'sup', actual+lsonp-1, sup, 6 )
c.debug
              actual = actual + lsonp

          else

              bfrpos = lsonp + 1
              bfrlen = 0
              
          end if

c         ----------------------------          
c         ... add in postponed columns
c         ----------------------------          

          midchl = .true.
          nn2    = nn2 - lsonp

          do 210 jcol = 1, lsonp

              l      = lsonp - jcol + 1
              pnlcol = pnlcol + 1
c.debug
c     write(6,'("in 210 - jcol, l, pnlcol, bfrpos = ", 4i8)') 
c    1            jcol, l, pnlcol, bfrpos
c.debug

              if ( pnlcol .lt. colbgn ) then
                  bfrpos = bfrpos + l + lson
                  go to 210
              end if

              if ( pnlcol .gt. colend ) then
                 return
              end if

c             ----------------------------------------------
c             ... read data for current column, if necessary
c             ----------------------------------------------

              if ( bfrpos-1+l+lson .gt. bfrlen ) then

                  k      = (bfrpos-1) - bfrlen
                  iopos2 = iopos2 + k 
                  totlen = totlen - k
                  bfrlen = min ( locbfr, totlen )
c.debug
c     write(6,'("totlen, bfrlen, iopos2      = ", 4i8)') 
c    1            totlen, bfrlen, iopos2
c.debug

                  call xdslw1 ( wafil2, 2, idummy, idummy, 
     1                          ocbufr, iopos2, bfrlen, ierr)
                  if ( ierr .ne. 0 ) go to 8000
c.debug
c     call xdslp5 ( 'ocbufr', bfrlen, ocbufr, 6 )
c.debug

                  watrn2 = watrn2 + bfrlen
                  iopos2 = iopos2 + bfrlen
                  totlen = totlen - bfrlen
                  bfrpos = 1

              end if

c             --------------------------
c             ... copy in diagonal block
c             --------------------------

              k        = pnlind(pnlcol,pnlcol)
c.debug
c     write(6,'("pnlcol, k, l, bfrpos        = ", 4i8)') 
c    1            pnlcol, k, l, bfrpos        
c.debug

              panel(k:k+l-1) = ocbufr(bfrpos:bfrpos+l-1)

              bfrpos = bfrpos + l

              k      = k + l

c             -------------------------------------------
c             ... zero out rows associated with postponed
c                 columns from other children
c             -------------------------------------------
c.debug
c     write(6,'("k, nn2                      = ", 4i8)') 
c    1            k, nn2                      
c.debug
   
              panel(k:k+nn2-1) = 0.d0

c             ------------------------------------------
c             ... add in information for last n3 columns
c             ------------------------------------------

              pnloff = k + nn2 - 1

              do j = 1, lson

                  pnlrow = lindxl(lndxpt+j)
c.debug
c     write(6,'("in 200 - j, pnlrow, n1      = ", 4i8)') 
c    1            j, pnlrow, n1               
c     write(6,'("in 200 - pnloff, bfrpos     = ", 4i8)') 
c    1            pnloff, bfrpos              
c.debug

                  if ( pnlrow .gt. n1 ) then

                      panel(pnloff+pnlrow-n1) = ocbufr(bfrpos)

                  end if

                  bfrpos = bfrpos + 1

              enddo

  210     continue 

          midchl = .false.
          lstchc = jstack

c         ---------------------------------
c         ... set up for next update matrix
c         ---------------------------------

  220     continue
          if  ( pnlcol .eq. colend + 1 )  then
             return
          else
             nstack = nstack - 1
          end if
      
  230 continue
c.debug
c     call xdslp5 ( 'panel at end of xdsla8', lpan, panel, 6 )
c.debug

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
