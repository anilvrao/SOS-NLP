      subroutine xdslic ( matrix, jcol,   nzcol,  jrowin, work,   lwork,
     1                    error )
 
c
c     purpose
c     -------
c
c     xdslic is the top level driver for structural input of the
c     matrix by columns.
c
c     created         27-jan-89   -- rgg --
c     modified        08-feb-91   -- rgg -- mods to allow no i/o to
c                                           sqfile
c     modified        08-feb-91   -- rgg -- added use of xdslil
c     modified        26-feb-91   -- rgg -- mods to allow unsymmetric
c                                           input
c                     12-dec-96   -- rgg -- rewritten to allow out
c                                           of memory assembly
c                     05-sep-97   -- rgg -- changes to allow combine
c                                           input processing routines
c                                           for solve and eigensolve
c                     01-nov-01   -- dkw -- changed error -106 to -105
c
c     input arguments
c     ---------------
c
c     matrix      c   character string denoting which matrix
c                     is being input.
c     jcol        i   number of column being input
c     nzcol       i   number of nonzeroes in lower triangle of the
c                     current column which may or not include the
c                     diagonal entry
c     jrowin      i   array containing the lists of row indices
c                     for the current column of a.
c     lwork       i   length of work array.
c
c     input/output arguments
c     ----------------------
c
c     work        d   work array.  if jcol =1 its contents have been
c                     initialized by subroutine xdslin.
c                     for jcol .gt. 1 its contents hold the partially
c                     assembled matrix structure.  when jcol .eq. neqns
c                     work contains the fully assembled matrix structure
c                     ordering phase to begin.
c
c     output arguments
c     ----------------
c
c     error       i   error flag
c                     =    0  normal return
c                     = -100  incorrect processing path
c                     = -101  lwork not large enough.
c                     = -103  nzcol .le. 0 or nzcol .gt. neqns-jcol+1
c                     = -104  illegal jcol
c                     = -105  illegal irow
c                     = -109  illegal value for matrix
c                     = -110  i/o error on wafil1
c
c     ------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------

      character(len=*)    matrix
 
      integer             jcol,   nzcol, lwork,  error
 
      integer             jrowin(*)
 
      double precision    work(*)
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------
 
      integer             acolls, anrecd, anzero, arowls, 
     1                    i     , ineqn1, irow,   itemp1, itemp2,
     2                    itemp3, k     , lstcol, maxnz , msglvl, 
     3                    neqns , output, sort  , srtlst, stage , 
     4                    wafil1, wkreqd, maxrec, imxrec
 
      double precision    t1,     w1
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      integer             xdslil, xdslni

      logical             lsame
 
      external            lsame , xdslil, xisli3, xdslni,
     1                    xdslt1, xdslt2, xislp1
 
c---------------------------------------------------------------------
 
      error  = 0
 
      stage  = work ( qstage )
      msglvl = work ( qmsglv )
      output = work ( qoutpu )
 
      if ( stage .ne. 1  .and.  stage .ne. 2 ) go to 8000
 
      neqns  = work ( qneqns )
      ineqn1 = xdslni(neqns+1)

      maxrec = min ( neqns, mxrecd )
      imxrec = xdslni(maxrec)

c     --------------------------
c     ... test validity of entry
c     --------------------------

      if ( .not. ( lsame ( matrix(1:1), 'A' ) .or.
     1             lsame ( matrix(1:1), 'K' ) .or.          
     2             lsame ( matrix(1:1), 'B' ) .or.
     3             lsame ( matrix(1:1), 'M' )      ) ) go to 8850
 
      if ( jcol .lt. 1 .or. jcol .gt. neqns ) go to 8400
 
      if ( nzcol .lt. 0 ) go to 8300 
 
      do i = 1, nzcol
 
          irow = jrowin(i)
 
          if ( irow .lt.    1  ) go to 8500
          if ( irow .gt. neqns ) go to 8500

          if ( lsame ( matrix(1:1), 'B' ) .or.
     1         lsame ( matrix(1:1), 'M' )      ) then
                   if ( irow .ne. jcol ) work(qnzlb) = work(qnzlb) + 1
          end if

      enddo 

      if ( stage .eq. 1 ) then

c         ------------------------------------------
c         ... first call so initialize input package
c         ------------------------------------------

          call xdslt1 ( t1, w1 )
          work(qinptm) = t1
          work(qinpwl) = w1

          k      = lwork - ( lncomm + 2*imxrec + ineqn1 ) 
          maxnz  = xdslil ( k / 3 )
c.debug
c         maxnz  = 150
c         maxnz  = 100
c         maxnz  =  50
c     write(6,'("in xdsliv - maxnz                 = ",i8)')
c    1                        maxnz 
c.debug

          wkreqd = 8*ineqn1 + 2*imxrec + lncomm

          if ( maxnz .lt. neqns ) go to 8100

c         --------------------------
c         ... partition up workspace
c         --------------------------

          acolls = lncomm + 2*imxrec + ineqn1 + 1
          arowls = acolls + xdslni ( maxnz )
          srtlst = arowls + xdslni ( maxnz )
c.debug
c     write(6,'("in xdslic - lncomm, neqns, ineqn1, acolls = ", 4i8)')
c    1                        lncomm, neqns, ineqn1, acolls
c.debug

          anzero = 0
          anrecd = 0
          lstcol = 0
          sort   = 0

          work(qnzla ) = anzero
          work(qnzab ) = maxnz
          work(qtemp1) = anrecd
          work(qtemp2) = lstcol
          work(qtemp3) = sort  
          work(qarowi) = arowls
          work(qacols) = acolls
          work(qtemp4) = srtlst

          work(qinuse) = srtlst + xdslni ( maxnz ) - 1
          work(qmxuse) = work(qinuse)

      else

c         ------------------------------------------------
c         ... not first call so extract needed information
c             from global.CMNication area
c         ------------------------------------------------

          maxnz  = work(qnzab )
          anzero = work(qnzla )
          anrecd = work(qtemp1)
          lstcol = work(qtemp2)
          sort   = work(qtemp3)
          
          arowls = work(qarowi)
          acolls = work(qacols) 
          srtlst = work(qtemp4)

      end if
 
c     ----------------------------------------------------
c     ... insert (irow,jcol) into row and column lists
c         representing the matrix structure.  If lists get
c         too big they are spilled to wafil1.
c     ----------------------------------------------------
 
      wafil1   = work(qwafl1)
 
      itemp1   = lncomm + 1 
      itemp2   = itemp1 + imxrec 
      itemp3   = itemp2 + imxrec
 
      call xisli3 ( jcol, nzcol, jrowin, neqns,
     1              anzero, maxnz, anrecd, maxrec, lstcol, sort,
     2              wafil1, work(arowls), work(acolls), work(srtlst),
     3              work(itemp1), work(itemp2), work(itemp3), error )
c.debug
c     write(6,'("after  xisli3 - error = ",i8)') error
c.debug
 
      if ( error .eq. -1 ) go to 8900
      if ( error .ne.  0 ) go to 8100
 
      work ( qstage ) = 2.
      work ( qnzla  ) = anzero
      work ( qtemp1 ) = anrecd
      work ( qtemp2 ) = lstcol
      work ( qtemp3 ) = sort  
 
      go to 9000
 
c---------------------------------------------------------------------
 
c     --------------
c     ... error trap
c     --------------
 
c     ---------------------------------
c     ... incorrect processing sequence
c     ---------------------------------
 
 8000 continue
      error = -100
      call hherr ( 3, 'xdslic', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88000 ) error, stage
      go to 8999
 
c     ---------------------------------------
c     ... insufficient storage for this stage
c     ---------------------------------------
 
 8100 continue
      error = -101
      call hherr ( 2, 'xdslic', error, wkreqd )
      if ( msglvl .gt. 0 ) write ( output, 88100 ) error, wkreqd, lwork
      go to 8999
 
c     -----------------------------
c     ... nzcol .lt. 0 .or. too big
c     -----------------------------
 
 8300 continue
      error = -103
      call hherr ( 1, 'xdslic', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88300 ) error, nzcol,
     1                                             jcol,  neqns
      go to 8999
 
c     --------------------------
c     ... illegal value for jcol
c     --------------------------
 
 8400 continue
      error = -104
      call hherr ( 1, 'xdslic', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88400 ) error, jcol,
     1                                             neqns
      go to 8999
 
c     --------------------------
c     ... illegal value for irow
c     --------------------------
 
 8500 continue
      error = -105
      call hherr ( 1, 'xdslic', error, 0 )
      if ( msglvl .gt. 0 ) then
          write ( output, 88500 ) error, irow, jcol, neqns
          call xislp3 ( 'row indices', nzcol, jrowin, output )
      end if
      go to 8999
 
c     ----------------------------
c     ... illegal value for matrix
c     ----------------------------
 
 8850 continue
      error = -109
      call hherr ( 1, 'xdslic', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88850 ) error, matrix
      go to 8999

c     ----------------------------
c     ... i/o error on file wafil1
c     ----------------------------

 8900 continue
      error = -110                 
      call hherr ( 3, 'xdslic', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88900 ) error, wafil1
      go to 8999

c     --------------------------------------
c     ... set stage to prevent further calls
c     --------------------------------------
 
 8999 continue
      work(qstage) = -1
 
c---------------------------------------------------------------------
 
c     ------------------------
c     ... end of module xdslic
c     ------------------------
 
 9000 continue
      return
 
c---------------------------------------------------------------------
 
c     -----------
c     ... formats
c     -----------
 
80000 format ( /5x, 'cpu  time for structural input          = ', f15.6
     1         /5x, 'wall time for structural input          = ', f15.6)
 
88000 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslic executed in an'
     2         /5x, 'incorrect sequence.  current stage = ',
     3              i10, 5x, 'should be 1 or 5.')
 
88100 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslie requires ', i15
     2         /5x, 'words of workspace and has only ', i15,
     3              ' available.' )
 
88300 format ( /5x, '*** fatal error no. ', i5, ' *** number of ',
     1              'nonzeroes input to'
     2         /5x, 'subroutine xdslic is incorrect.'
     3        /10x, 'nzcol = ', i10, '  jcol = ', i10,
     4              '  neqns = ', i10 )
 
88400 format ( /5x, '*** fatal error no. ', i5, ' *** invalid ',
     1              'value of jcol as to'
     2         /5x, 'subroutine xdslic.'
     3        /10x, 'jcol = ', i10, '  neqns = ', i10 )
 
88500 format ( /5x, '*** fatal error no. ', i5, ' *** invalid ',
     1              'value of irow input to'
     2         /5x, 'subroutine xdslic.  irow = ', i10,
     3              '  jcol = ', i10, ' neqns = ', i10 )
 
88850 format ( /5x, '*** fatal error no. ', i5, ' *** invalid ',
     1              'value of matrix input to'
     2         /5x, 'subroutine xdslic.'
     3        /10x, 'matrix = ', a )

88900 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslic encountered i/o error'
     2         /5x, 'on i/o file no. ', i15 )
 
c---------------------------------------------------------------------
 
      end
