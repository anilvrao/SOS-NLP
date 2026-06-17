      subroutine xdsla1( mincor, diag  , xladj , ladjl , lnza  ,
     1                   bassmb, sigma , bmxtyp, bdiag , bxladj, bladjl,
     2                   lnzb  , loclfr, nodbgn, nodnxt, pospon,
     3                   front , sqfil1, sqtrn1, sqfil3, sqtrn3,
     4                   mxvlib, mxvlrb, fctops, ierr )
 
c
c  purpose -- add the original matrix entries in the frontal matrix 
c             for the current supernode.  if sigma .ne. 0. then 
c             assemble a - sigma * b.
c
c  created   -- 04-feb-97, rgg -- extracted from old version of xdslf4
c                                 for better modularity of code.
c
c  variables
c
c      mincor -- minimum core processing switch.  if mincor .ne. 0
c                then either a or b are on i/o file sqfil1 or sqfil3.
c      diag   -- array to hold first the diagonal entries of the
c                original matrix and later those of the factor
c      xladj  -- pointers into indices and entries of the matrix
c      ladjl  -- local indices of the column indices of original matrix
c      lnza   -- column entries of the original matrix
c      sigma  -- shift value
c      bmxtyp -- matrix type of b
c                0 - null matrix
c                1 - general sparse
c                2 - general sparse but same structure as a
c                3 - diagonal
c                4 - identity
c      bdiag  -- diagonal of b
c      bxladj -- xladj for b
c      bladjl -- ladjl for b
c      lnzb   -- column entries of the b matrix
c      front  -- array holding representation of current front
c      sqfil1 -- i/o file holding a
c      sqfil3 -- i/o file holding b
c
c  working storge  --
c
c      mxvlib -- i/o buffer for bringing in scatter indicies from the
c                spilled matrix file, sqfil1 or sqfil3.
c      mxvlrb -- i/o buffer for bringing in bringing in the matrix
c                values if mincor .ne. 0.
c
c  input/output variables
c
c      sqtrn1 -- amount of i/o transfered to and from sqfil1.
c      sqtrn3 -- amount of i/o transfered to and from sqfil3.
c
c  output variables --
c
c      fctops -- number of floating point operations to compute
c                factorization
c      ierr   -- error code
c                =  0 normal return
c                = -2 i/o error on sqfil1 or sqfil3
c
c  subprograms called
c
c      xdslf9 -- scatter/add a vector into another
c
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           bmxtyp, loclfr, nodbgn, nodnxt, 
     1                  pospon, mincor, sqfil1, sqfil3, ierr
 
      integer           xladj(*),       ladjl(*),
     1                  bxladj(*),      bladjl(*),
     2                  mxvlib(*)
 
      logical           bassmb
 
      double precision  sigma , fctops, sqtrn1, sqtrn3
 
      double precision  diag(*),        lnza(*),
     1                  bdiag(*),       lnzb(*),
     2                  front(*),       mxvlrb(*)

c     -------------------
c     ... local variables
c     -------------------
 
      integer           amncor, asize,  bmncor, bsize, ii,
     1                  k     , length, node  , 
     2                  starta, startb 
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      integer           xdslni
 
      external          xdslni, xdslf9, xdslvr, xislvr
c
c  =====================================================================
 

c     ---------------------------------
c     ... add in original matrix values
c     ---------------------------------
 
      starta = xladj(nodbgn)
      asize  = xladj(nodnxt) - starta
c.debug
c     write(6,'("in xdsla1 - starta, asize, mincor = ", 3i8)')
c    1                        starta, asize, mincor 
c.debug

      amncor = mod ( mincor, 10 )

      if ( amncor .ne. 0 ) then
 
c         ----------------------------------------------
c         ... a is out-of-core.  read in off-diag values
c             and add into front
c         ----------------------------------------------

          if ( asize .gt. 0 ) then
c.debug
c     write(6,'("before reads for a")')
c     write(6,'("sqfil1, asize, amncor             = ", 3i8)')
c    1            sqfil1, asize, amncor             
c.debug
 
              if ( amncor .eq. 2 ) then
                  call xislvr ( sqfil1, asize, mxvlib, ierr )
                  if ( ierr .ne. 0 ) go to 8800
                  sqtrn1 = sqtrn1 + xdslni ( asize )
              end if

              call xislvr ( sqfil1, asize, mxvlib, ierr )
              if ( ierr .ne. 0 ) go to 8800
              sqtrn1 = sqtrn1 + xdslni ( asize )
c.debug
c     write(6,'("after integer read in xdsla1")')
c     call xislp3 ( 'org. entries', asize, mxvlib, 6 )
c.debug

              call xdslvr ( sqfil1, asize, mxvlrb, ierr )
              if ( ierr .ne. 0 ) go to 8800
              sqtrn1 = sqtrn1 + asize

c.debug
c     write(6,'("before out-of-core call to xdslf9")')
c     call xdslp5 ( 'org. entries', asize, mxvlrb, 6 )
c     call xislp3 ( 'org. entries', asize, mxvlib, 6 )
c.debug
              call xdslf9 ( asize, mxvlrb, mxvlib, front )
c.debug
c     write(6,'("after call to xdslf9")')
c.debug

          end if
 
      else
 
c         --------------------------------------------------
c         ... a is in-core.  add off-diag entries into front
c         --------------------------------------------------
 
c.debug
c     write(6,'("before in-core call to xdslf9")')
c     call xdslp5 ( 'org. entries', asize, lnza (starta), 6 )
c     call xislp3 ( 'org. entries', asize, ladjl(starta), 6 )
c.debug
          call xdslf9 ( asize, lnza (starta), ladjl (starta), front )
c.debug
c     write(6,'("after call to xdslf9")')
c.debug
 
      end if
 
      fctops = fctops + asize

c     -----------------------------------------------
c     ... add in original diagonal entries into front
c     -----------------------------------------------
 
      k      = 1
      length = loclfr 
 
      do node = nodbgn, nodnxt - 1
          front (k) = front (k) + diag (node)
          k         = k + length
          length    = length - 1
      enddo
c.debug
c     write(6,'("after 100 loop")')
c.debug

c     ---------------------------------------
c     ... test if b is to be included as well
c     ---------------------------------------
c.debug
c     write(6,'("bassmb, bmxtyp = ", l8, i8 )') bassmb, bmxtyp
c.debug 
      if ( bassmb ) then

c         ------------------------------------
c         ... handle the diagonal entries of b
c         ------------------------------------
 
          k      = 1
          length = loclfr 

          if ( bmxtyp .le. 3 ) then

c             -----------------------------------------
c             ... b is has non-trivial diagonal entries
c             -----------------------------------------

              do node = nodbgn, nodnxt - 1
                  front(k) = front(k) - sigma * bdiag(node)
                  k        = k + length
                  length   = length - 1
              enddo
c.debug
c     write(6,'("after 110 loop")')
c.debug
 
              fctops = fctops + 2 * ( nodnxt - nodbgn )
 
          else if ( bmxtyp .eq. 4 ) then

c             ----------------------------
c             ... b is the identity matrix
c             ----------------------------
 
              do node = nodbgn, nodnxt - 1
                  front(k) = front(k) - sigma
                  k        = k + length
                  length   = length - 1
              enddo
c.debug
c     write(6,'("after 120 loop")')
c.debug
 
              fctops = fctops + ( nodnxt - nodbgn )
 
          end if
 
          if ( bmxtyp .eq. 1 ) then

c             -------------------------
c             ... b is a general matrix
c             -------------------------
 
              startb = bxladj(nodbgn)
              bsize  = bxladj(nodnxt) - startb
              bmncor = mincor / 10
c.debug
c     write(6,'("startb, bsize, bmncor = ", 3i8)')
c    1            startb, bsize, bmncor 
c.debug
 
              if ( bmncor .ne. 0 ) then
 
c                 --------------------
c                 ... b is out-of-core
c                 --------------------
 
                  if ( bsize .gt. 0 ) then

                      if ( bmncor .eq. 2 ) then
                          call xislvr ( sqfil3, bsize, mxvlib, ierr )
                          if ( ierr .ne. 0 ) go to 8800
                          sqtrn3 = sqtrn3 + xdslni ( bsize )
                      end if

                      call xislvr ( sqfil3, bsize, mxvlib, ierr )
                      if ( ierr .ne. 0 ) go to 8800
                      sqtrn3 = sqtrn3 + xdslni ( bsize )

                      call xdslvr ( sqfil3, bsize, mxvlrb, ierr )
                      if ( ierr .ne. 0 ) go to 8800
                      sqtrn3 = sqtrn3 + bsize
c.debug
c     write(6,'("before out-of-core call")')
c     call xdslp5 ( 'org. entries', bsize, mxvlrb, 6 )
c     call xislp3 ( 'org. entries', bsize, mxvlib, 6 )
c.debug
 
                      do ii=1,bsize
                        k = mxvlib(ii)
                        front(k) = front(k) - sigma*mxvlrb(ii)
                      enddo
c.debug
c     write(6,'("after  out-of-core call")')
c.debug
 
                  end if
 
              else
 
c                 ----------------
c                 ... b is in-core
c                 ----------------
 
c.debug
c     write(6,'("before in-core call")')
c     call xdslp5 ( 'org. entries', bsize, lnzb  (startb), 6 )
c     call xislp3 ( 'org. entries', bsize, bladjl(startb), 6 )
c.debug
                  do ii=startb,bxladj(nodnxt)-1
                    k = bladjl(ii)
                    front(k) = front(k) - sigma*lnzb(ii)
                  enddo
c.debug
c     write(6,'("after  in-core call")')
c.debug
 
              end if
 
              fctops = fctops + 2 * bsize
 
          end if
 
      end if

      return
 
c  =====================================================================
 
c     ------------------------------------------------
c     ... error trap for i/o error on sqfil1 or sqfil3
c     ------------------------------------------------
 
 8800 continue
      ierr = -2
      return
 
c  =====================================================================
 
      end
