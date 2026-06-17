      subroutine xdsla5( mincor, nodbgn, nodnxt, 
     1                   diag  , xladj , ladjl , lnza  ,
     1                   bassmb, sigma , bmxtyp, bdiag , bxladj, 
     2                   bladjl, lnzb  , colbgn, colend, poffst, 
     3                   lpan  , loclfr, panel , sqfil1, sqtrn1, 
     4                   sqfil3, sqtrn3, qfirst, mxvlib, mxvlrb, 
     5                   fctops, ierr )
c
c  purpose -- add the original matrix entries in the frontal matrix 
c             for the current supernode.  if sigma .ne. 0. then 
c             assemble a - sigma * b.
c
c  created   -- 04-feb-97, rgg -- extracted from old version of xdslai
c                                 for better modularity of code.
c
c  variables
c
c      mincor -- minimum core processing switch.  if mincor .ne. 0
c                then either a or b are on i/o file sqfil1 or sqfil3.
c      nodbgn -- first column of org. matrix for this front
c      nodnxt -- first column of org. matrix for next front
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
c      colbgn -- first column of this section of the front
c      colend -- last column of this section of the front
c      poffst -- offset in front for the start of this panel
c      lpan   -- length of this panel
c      loclfr -- size of current section of front
c      sqfil1 -- i/o file holding a 
c      sqfil3 -- i/o file holding b 
c      qagain -- logical variable.  if true this is not the first
c                time through this subroutine for this front.
c
c  working storge  --
c
c      mxvlib -- i/o buffer for bringing in scatter indicies from the
c                spilled matrix file, sqfil1 or sqfil3.
c      mxvlrb -- i/o buffer for bringing in matrix values from the
c                spilled stack.
c
c  input/output variables
c
c      panel  -- array holding numerical values for section of 
c                current front
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
c      xdsla6 -- scatter/add a vector into another
c
c  =====================================================================
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer           nodbgn, nodnxt, bmxtyp, loclfr, poffst,
     1                  lpan  , colbgn, colend, mincor, sqfil1, 
     2                  sqfil3, ierr
 
      integer           xladj(*),       ladjl(*),
     1                  bxladj(*),      bladjl(*),
     2                  mxvlib(*)
 
      logical           bassmb, qfirst
 
      double precision  sigma , fctops, sqtrn1, sqtrn3
 
      double precision  diag(*),        lnza(*),
     1                  bdiag(*),       lnzb(*),
     2                  mxvlrb(*),      panel(*)

c     -------------------
c     ... local variables
c     -------------------
 
      integer           amncor, asize , bmncor, bsize , k     ,
     1                  length, node  , nodoff, starta, startb 
 
      double precision  one                   
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      integer           xdslni
 
      external          xdslni, xdslf9, xdslvr, xislvr
 
c  =====================================================================
 
c     ---------------------------------
c     ... add in original matrix values
c     ---------------------------------

      one    = 1.0d0
      nodoff = nodbgn - 1
 
      starta = xladj(nodbgn)
      asize  = xladj(nodnxt) - starta
      amncor = mod ( mincor, 10 )
c.debug
c     write(6,'("in xdsla5 - starta, asize, mincor = ", 3i8)')
c    1                        starta, asize, mincor 
c     write(6,'("            qfirst                = ",  l8)')
c    1                        qfirst                
c.debug

c     ----------------------------------------------
c     ... a is out-of-core.  read in off-diag values
c         and add into front
c     ----------------------------------------------

      if ( amncor .ne. 0 ) then 

          if ( asize .gt. 0 ) then

              if ( qfirst ) then
 
                  if ( amncor .eq. 2 ) then
                      call xislvr ( sqfil1, asize, mxvlib, ierr )
                      if ( ierr .ne. 0 ) go to 8800
                      sqtrn1 = sqtrn1 + xdslni ( asize )
                  end if
c.debug
c     write(6,'("after second call to xislvr")')
c.debug
    
                  call xislvr ( sqfil1, asize, mxvlib, ierr )
                  if ( ierr .ne. 0 ) go to 8800
                  sqtrn1 = sqtrn1 + xdslni ( asize )
c.debug
c     write(6,'("after integer read in xdsla5")')
c     call xislp3 ( 'org. entries', asize, ocbufi, 6 )
c.debug

                  call xdslvr ( sqfil1, asize, mxvlrb, ierr )
                  if ( ierr .ne. 0 ) go to 8800
                  sqtrn1 = sqtrn1 + asize

              end if

c.debug
c     write(6,'("before out-of-core call to xdsla6")')
c     write(6,'("poffst, lpan = ", 2i8)') poffst, lpan
c     call xdslp5 ( 'org. entries', asize, mxvlrb, 6 )
c     call xislp3 ( 'org. entries', asize, mxvlib, 6 )
c.debug
              call xdsla6 ( one, asize, mxvlrb, mxvlib, 
     1                      poffst, lpan, panel ) 
c.debug
c     write(6,'("after call to xdsla6")')
c     call xdslp5 ( 'panel', lpan, panel, 6 )
c.debug

          end if

      else       

c         --------------------------------------------------
c         ... a is in-core.  add off-diag entries into front
c         --------------------------------------------------

c.debug                                                      
c     write(6,'("before in-core call to xdsla6")')
c     call xdslp5 ( 'org. entries', asize, lnza (starta), 6 )
c     call xislp3 ( 'org. entries', asize, ladjl(starta), 6 )
c.debug
          call xdsla6 ( one, asize, lnza(starta), ladjl(starta),
     1                  poffst, lpan, panel ) 
c.debug
c     write(6,'("after call to xdsla6")')
c     call xdslp5 ( 'panel', lpan, panel, 6 )
c.debug

      end if
 
      fctops = fctops + asize

c     -----------------------------------------------
c     ... add in original diagonal entries into front
c     -----------------------------------------------
 
      k      = 1
      length = loclfr 
c.debug
c     write(6,'("before 100 loop - colbgn, colend = ",2i8)')
c    1                              colbgn, colend
c.debug
 
      do node = colbgn, colend
c.debug
c         write(6,'("node,nodoff,k,diag(node+nodoff),panel(k)=", 
c    1                3i8, 1p2d15.5)')
c    1                node,nodoff,k,diag(node+nodoff),panel(k)
c.debug
          panel (k) = panel (k) + diag (node+nodoff)
          k         = k + length
          length    = length - 1
      enddo

      fctops = fctops + ( colend - colbgn + 1 )

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

              do node = colbgn, colend 
                  panel(k) = panel(k) - sigma * bdiag(node+nodoff)
                  k        = k + length
                  length   = length - 1
              enddo
 
              fctops = fctops + 2 * ( colend - colbgn + 1 )
 
          else

c             ----------------------------
c             ... b is the identity matrix
c             ----------------------------
 
              do node = colbgn, colend 
                  panel(k) = panel(k) - sigma
                  k        = k + length
                  length   = length - 1
              enddo
 
              fctops = fctops + ( colend - colbgn + 1 )
 
          end if
 
          if ( bmxtyp .eq. 1 ) then

c             -------------------------
c             ... b is a general matrix
c             -------------------------
 
              startb = bxladj(nodbgn)
              bsize  = bxladj(nodnxt) - startb
              bmncor = mincor / 10
 
              if ( bmncor .ne. 0 ) then
 
c                 --------------------
c                 ... b is out-of-core
c                 --------------------
 
                  if ( bsize .gt. 0 ) then

                      if ( qfirst ) then

                          if ( bmncor .eq. 2 ) then 
                              call xislvr ( sqfil3, bsize,
     1                                      mxvlib(asize+1), ierr )
                              if ( ierr .ne. 0 ) go to 8800
                              sqtrn3 = sqtrn3 + xdslni ( bsize )
                          end if

                          call xislvr ( sqfil3, bsize, mxvlib(asize+1), 
     1                                  ierr )
                          if ( ierr .ne. 0 ) go to 8800
                          sqtrn3 = sqtrn3 + xdslni ( bsize )

                          call xdslvr ( sqfil3, bsize, mxvlrb(asize+1), 
     1                                  ierr )
                          if ( ierr .ne. 0 ) go to 8800
                          sqtrn3 = sqtrn3 + bsize

                      end if
c.debug
c     write(6,'("before call to xdsal6")')
c     call xdslp5 ( 'org. entries', bsize, mxvlrb(asize+1), 6 )
c     call xislp3 ( 'org. entries', bsize, mxvlib(asize+1), 6 )
c.debug
 
                      call xdsla6 ( -sigma, bsize,
     1                               mxvlrb(asize+1), mxvlib(asize+1), 
     2                               poffst, lpan, panel ) 
c.debug
c     write(6,'("after call to xdsla6")')
c     call xdslp5 ( 'panel', lpan, panel, 6 )
c.debug
 
                  end if
 
              else
 
c                 ----------------
c                 ... b is in-core
c                 ----------------
c.debug
c     write(6,'("before call to xdsal6")')
c     write(6,'("asize, bsize, startb = ", 3i8)') 
c    1            asize, bsize, startb
c     call xdslp5 ( 'org. entries', bsize, mxvlrb(asize+1), 6 )
c     call xislp3 ( 'org. entries', bsize, mxvlib(asize+1), 6 )
c.debug
 
                  call xdsla6 ( -sigma, bsize, lnzb(startb),
     1                          bladjl(startb), poffst, lpan, panel ) 
c.debug
c     write(6,'("after call to xdsla6")')
c     call xdslp5 ( 'panel', lpan, panel, 6 )
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
