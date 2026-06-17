      subroutine   xdslwm   ( mode, utitle, mtxunt, work, lwork, error )
c
c     ==================================================================
c     ====  xdslwm -- writes matrices in either ascii or binary     ====
c     ====            format to i/o unit                            ====
c     ==================================================================
c     ==================================================================
c
c     purpose
c     -------
c
c     xdslwm is the user callable subroutine to write the a and b
c     matrices to i/o unit mtxunt
c
c     created         24-feb-93   -- rgg --
c     last modified   12-aug-98   -- rgg -- modified for release 4
c                                           added choice of ascii or
c                                           binary.
c                     10-sep-01   -- dkw -- removed duplicate format
c                                           statement.
c                     01-nov-01   -- dkw -- changed error -400 to
c                                           -451 and -500
c
c     input arguments
c     ---------------
c
c     mode        ch  character string denoting ascii or binary
c     utitle      ch  character string denoting the title for the matrix
c     mtxunt      i   i/o unit to use.
c     lwork       i   length of work array.
c     work        d   work array.  on input it contains the
c                    .CMNication area and all active arrays.
 
c     output arguments
c     ----------------
c
c     error       i   error flag
c                     = -401  lwork not large enough.
c                     = -402  i/o error on mtxunt.
c                     = -410  i/o error on sqfil1 or sqfil3
c                     = -451  original matrices are not saved
c                     = -500  incorrect processing path.
c
c     ==================================================================
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      character(len=*)    mode, utitle
 
      integer             mtxunt, lwork, error
 
      double precision    work   (*)
 
c     -----------------------------
c     ... global.CMNication area
c     -----------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     -------------------
c     ... local variables
c     -------------------
 
      character(len=80)   title

      integer             acolst, adiag , amncor, amxtyp, arowin,
     1                    astage, avalue, bcolst, bdiag , bmncor,
     2                    bmxtyp, browin, bvalue, inuse ,
     3                    mincor, msglvl, mtxtyp, mxanzf, mxvlib,
     3                    mxvlrb, neqns , nsuper, nzla  , nzlb  ,
     4                    output, perm  , sqfile, wkreqd, xsup 

      logical             qascii

c     ---------------------
c     ... external function
c     ---------------------

      integer             xdslni

      logical             lsame

      external            lsame, xdslni
 
c     ==================================================================
 
c     -----------------------------------------------------
c     ... get stage information from the.CMNication area
c     -----------------------------------------------------
 
      error  = 0
      astage = work (qstage)
      msglvl = work (qmsglv)
      output = work (qoutpu)
 
      if  ( astage .ne. 40 )  go to 8450

      if ( work(qsavea) .eq. 0. ) go to 8451

      mtxtyp = work ( qmxtyp )

      mincor = work ( qmncor )
      amncor = mod ( mincor, 10 ) 
      bmncor = mincor / 10
 
      neqns  = work (qneqns)
      bmxtyp = work (qbmxty)
 
      perm   = work (qperm)
 
      nzla   = work (qnzla )
      acolst = work (qxadj )
      arowin = work (qarowi)
      adiag  = work (qadiag)
      avalue = work (qofdia)
 
      if ( bmxtyp .eq. 1 ) then 
          nzlb   = work (qnzlb )
          bcolst = work (qbcols)
          browin = work (qbrowi)
          bdiag  = work (qbdiag)
          bvalue = work (qbvalu)
      else if ( bmxtyp .eq. 3 ) then
          bdiag  = work (qbdiag)
          nzlb   = 0
      else
          nzlb   = 0
      end if

      qascii = lsame ( mode(1:1), 'A' ) 
 
c     ==================================================================
 
c     --------------------------------------------------
c     ... write the matrices out to mtxunt.
c         first write header and new to old permutation.
c     --------------------------------------------------
 
      title = utitle
      nzla  = nzla + neqns
      if ( bmxtyp .eq. 1 ) then
          nzlb = nzlb + neqns
      else if ( bmxtyp .eq. 3 ) then
          nzlb = neqns
      else
          nzlb = 0
      end if

      if ( qascii ) then
 
c         ----------------
c         ... ascii format
c         ----------------

          write ( mtxunt, '(a)'    ) title
          write ( mtxunt, '(5i10)' ) mtxtyp, neqns, nzla, bmxtyp, nzlb
          call xislwa ( mtxunt, '(8i10)', neqns, work(perm) )

      else
 
c         -----------------
c         ... binary format
c         -----------------

          write ( mtxunt ) title
          write ( mtxunt ) mtxtyp, neqns, nzla, bmxtyp, nzlb
          call xislvw ( mtxunt, neqns, work(perm), error )
          if ( error .ne. 0 ) go to 8900

      end if
 
c     ==================================================================

c     ----------------------
c     ... write out matrix a
c     ----------------------

      amxtyp = 1

      if ( amncor .eq. 0 ) then 

          call xdslwb ( mtxunt, qascii, mtxtyp, amxtyp, neqns,
     1                  work(perm), work(adiag), work(acolst), 
     2                  work(arowin), work(avalue) )

      else

          inuse  = work ( qinuse )
          mxanzf = work ( qmxanz )
          mxvlib = inuse + 1
          mxvlrb = mxvlib + xdslni ( mxanzf )
          wkreqd = mxvlrb + mxanzf - 1
       
          if ( lwork .lt. wkreqd ) go to 8100

          nsuper = work ( qnsupe )
          xsup   = work ( qxsup  )
          sqfile = work ( qsqfl1 )
 
          call xdslwc ( mtxunt, qascii, mtxtyp, amxtyp, neqns, 
     1                  nsuper, sqfile,
     2                  work(perm  ), work(xsup  ), 
     3                  work(adiag ), work(acolst), 
     4                  work(mxvlib), work(mxvlrb), error )

          if ( error .ne. 0 ) go to 8900

      end if
 
c     ==================================================================

c     ----------------------
c     ... write out matrix b
c     ----------------------

      if ( bmxtyp .ne. 1 .and. bmxtyp .ne. 3 ) go to 8999

      if ( bmncor .eq. 0 ) then 

          call xdslwb ( mtxunt, qascii, mtxtyp, bmxtyp, neqns,
     1                  work(perm  ), work(bdiag ), work(bcolst),
     2                  work(browin), work(bvalue) )

      else

          inuse  = work ( qinuse )
          mxanzf = work ( qmxanz )
          mxvlib = inuse + 1
          mxvlrb = mxvlib + xdslni ( mxanzf )
          wkreqd = mxvlrb + mxanzf - 1
       
          if ( lwork .lt. wkreqd ) go to 8100

          nsuper = work ( qnsupe )
          xsup   = work ( qxsup  )
          sqfile = work ( qsqfl3 )

          call xdslwc ( mtxunt, qascii, mtxtyp, bmxtyp, neqns, 
     1                  nsuper, sqfile,
     2                  work(perm  ), work(xsup  ), 
     3                  work(bdiag ), work(bcolst), 
     4                  work(mxvlib), work(mxvlrb), error )

          if ( error .ne. 0 ) go to 8900

      end if
 
c     ==================================================================
 
      go to 8999
 
c     ==================================================================
 
c     ------------------
c     ... error handling
c     ------------------
 
c     ---------------------------------------
c     ... insufficient storage for this stage
c     ---------------------------------------
 
 8100 continue
      error = -401
      call hherr ( 2, 'xdslwm', error, wkreqd )
      if ( msglvl .gt. 0 ) write ( output, 88100 ) error, wkreqd, lwork
      go to 8999  
 
c     ---------------------------------
c     ... incorrect processing sequence
c     ---------------------------------
 
 8450 continue
      error = -500
      call hherr ( 3, 'xdslwm', error, 0 )
      if  ( msglvl .gt. 0 )  then
          if ( astage .ne. 40 ) write ( output, 88450 ) error, astage
      endif
      go to 8999

c     ---------------------------------
c     ... matrices have not been saved
c     ---------------------------------

 8451 continue
      error = -451
      call hherr ( 3, 'xdslwm', error, 0 )
      if  ( msglvl .gt. 0 )  then
          if ( work(qsavea) .eq. 0. ) write ( output, 88451 )
      endif
      go to 8999

c     ------------------------------------
c     ... i/o error while reading matrices
c     ------------------------------------
 
 8900 continue
      error = -410
      call hherr ( 3, 'xdslwm', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88900 ) error, sqfile
      go to 8999
 
c     ==================================================================

 8999 continue
      return
 
c     ==================================================================
 
c     -----------
c     ... formats
c     -----------
 
88100 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslwm requires ', i15
     2         /5x, 'words of workspace and has only ', i15,
     3              ' available.' )
 
88450 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslwm  executed out of sequence'
     2        /10x, 'current stage for a: ', i10
     3        /10x, 'should be 40' )      
 
88451 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslwm executed when the original matrices are '
     2              'not saved -- see xdslsp.' )
 
88499 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslwm  -- i/o error on mtxunt = ', i10 )

88900 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslwm encountered i/o error'
     2         /5x, 'on i/o file no. ', i15 )
 
      end
