      subroutine xdslvc ( matrix, jcol  , nzcol , jrowin, value,
     1                    work  , lwork , error )
 
c
c     purpose
c     -------
c
c     xdslvc is the top level driver for numeric input of a column
c     of data at a time.
c
c     created         27-jan-89   -- rgg --
c     modified        29-oct-90   -- mlc -- call to xislvr protected
c                                           against nzla .le. 0
c     modified        08-feb-91   -- rgg -- mods to allow no i/o to sqfi
c     modified        19-dec-96   -- rgg -- mods for out-of-core
c                                           assembly
c     modified        05-sep-97   -- rgg -- mods to combine solver
c                                           and eigensolver input.
c     modified        08-jun-06   -- dkw -- fixed error testing from 1
c                                           to -1 after calling xdslv3
c
c     input arguments
c     ---------------
c
c     matrix      c   character string denoting which matrix is
c                     being referenced.
c     jcol        i   column number.
c     nzcol       i   no. of nonzeroes in the column.
c     jrowin      i   row indicies for the column.
c     value       d   values for the column.
c     lwork       i   length of work array.
c
c     input/output arguments
c     ----------------------
c
c     work        d   work array.  on input it contains the
c                    .CMNication area and all active arrays.
c
c     output arguments
c     ----------------
c
c     error       i   error flag
c                     =    0  normal return
c                     = -400  incorrect processing path.
c                     = -401  lwork not large enough.
c                     = -402  error in matrix input.
c                     = -409  illegal value for matrix
c                     = -410  i/o error on wafile
c
c---------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------

      character(len=*)    matrix

      integer             jcol,   nzcol , lwork,  error

      integer             jrowin(*)
 
      double precision    value (*),      work(*)
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------

      integer             acolls, anrecd, anzero, arowls, anzp  , 
     1                    adiag , arclen, arcpos, 
     2                    bcolls, bnrecd, bnzero, browls, bnzp  , 
     3                    bdiag , brclen, brcpos, bmxtyp,
     4                    colstr, i     , imaxnz, imxrec, ineqns, 
     5                    inuse , invp  , irow  ,
     5                    k     , kmatrx, lbound, maxnz , maxrec,
     6                    msglvl, neqns , output, 
     7                    srtlst, srtval, stage , wafile, wkreqd

      logical             unsym
 
      double precision    kmult , t1,     w1
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      integer             xdslni
 
      logical             lsame
 
      external            lsame , xdslv3, xdslni, xdslt1
 
c---------------------------------------------------------------------
 
      error  = 0
 
c     ----------------------------------------------------
c     ... extract information from the.CMNication area.
c     ----------------------------------------------------
 
      stage  = work ( qstage )
      msglvl = work ( qmsglv )
      output = work ( qoutpu )
      unsym  = work ( qmxtyp ) .eq. 2.
      bmxtyp = work ( qbmxty )
 
      neqns  = work ( qneqns )
      ineqns = xdslni ( neqns )

      maxrec = min ( neqns, mxrecd )
      imxrec = xdslni ( maxrec )
 
      if ( stage .ne. 31 ) then
 
          if ( stage .ne. 30 .and. stage .ne. 40
     1                       .and. stage. ne. 50 ) go to 8000
 
c         ------------------------------------------
c         ... initialize matrix value input by entry
c         ------------------------------------------
 
          if ( msglvl .ge. 2 ) write ( output, 68000 )
 
          if ( stage .eq. 30 ) then
              inuse  = work ( qinuse )
          else
              inuse  = work ( qadiag ) - 1
          end if

          adiag  = inuse + 1
          inuse  = adiag + neqns - 1

          if ( bmxtyp .le. 3 ) then
              bdiag  = inuse + 1
              inuse  = bdiag + neqns - 1
          else
              bdiag  = 1
          end if

          arcpos = inuse  + 1
          arclen = arcpos + imxrec
          inuse  = arclen + imxrec - 1
          kmatrx = 1

          if ( bmxtyp .eq. 1 ) then
              brcpos = inuse  + 1
              brclen = brcpos + imxrec 
              inuse  = brclen + imxrec - 1
              kmatrx = 2
          else
              brcpos = 1
              brclen = 1
          end if
 
          k      = lwork - inuse - ineqns
          kmult  = (2*kmatrx+1.) * xdslni(64) + (kmatrx+1.) * 64 
          kmult  = kmult / 64.
          maxnz  = ( k - kmult + 1 ) / kmult
c.debug
c         write(6,'("value input - k, maxnz, kmult = ", 2i8,f8.2)')
c    1                              k, maxnz, kmult 
c.debug

          if ( maxnz .lt. neqns ) then
              wkreqd = inuse + (2*kmatrx+2)*ineqns
     1                       + (  kmatrx+1)* neqns
              go to 8100
          end if

          imaxnz = xdslni ( maxnz )

          acolls = inuse  + 1
          arowls = acolls + imaxnz
          anzp   = arowls + imaxnz
          inuse  = anzp   + maxnz - 1

          if ( bmxtyp .eq. 1 ) then
              bcolls = inuse  + 1
              browls = bcolls + imaxnz
              bnzp   = browls + imaxnz
              inuse  = bnzp   + maxnz - 1
          else
              bcolls = 1
              browls = 1
              bnzp   = 1
          end if
   
          srtlst = inuse  + 1
          srtval = srtlst + imaxnz
          colstr = srtval + maxnz

          inuse  = colstr + ineqns - 1
c.debug
          if ( inuse .gt. lwork ) then
              write(6,'("storage oops in xdslvc")')
              RETURN
          end if
c.debug
         
          call xdslt1 ( t1, w1 )
          work ( qvaltm ) = t1
          work ( qvalwl ) = w1
 
          work ( qstage ) = 31
          work ( qinuse ) = inuse 
          work ( qmxuse ) = max ( work(qmxuse), work(qinuse) )
 
          work(adiag:adiag+neqns-1) = 0.d0
          if ( bmxtyp .le. 3 ) work(bdiag:bdiag+neqns-1) = 0.d0

          anzero = 0
          anrecd = 0
          bnzero = 0
          bnrecd = 0

          work(qnzab ) = maxnz
          work(qtemp5) = srtlst

          work(qacols) = acolls
          work(qarowi) = arowls
          work(qadiag) = adiag
          work(qofdia) = anzp 

          work(qtemp1) = anzero
          work(qtemp2) = anrecd
          work(qtemp3) = arcpos
          work(qtemp4) = arclen

          work(qbcols) = bcolls
          work(qbrowi) = browls
          work(qbdiag) = bdiag
          work(qbvalu) = bnzp 

          work(qtemp6) = bnzero
          work(qtemp7) = bnrecd
          work(qtemp8) = brcpos
          work(qtemp9) = brclen

      else

          maxnz  = work(qnzab ) 
          imaxnz = xdslni ( maxnz )

          acolls = work(qacols) 
          arowls = work(qarowi)
          adiag  = work(qadiag) 
          anzp   = work(qofdia)

          anzero = work(qtemp1) 
          anrecd = work(qtemp2)
          arcpos = work(qtemp3) 
          arclen = work(qtemp4)

          bcolls = work(qbcols) 
          browls = work(qbrowi)
          bdiag  = work(qbdiag) 
          bnzp   = work(qbvalu)

          bnzero = work(qtemp6) 
          bnrecd = work(qtemp7)
          brcpos = work(qtemp8) 
          brclen = work(qtemp9)
 
          srtlst = work(qtemp5) 
          srtval = srtlst + imaxnz
          colstr = srtval + maxnz
 
      end if

c     -------------------------
c     ... validate the entries.
c     -------------------------
 
      if ( jcol .lt. 1      ) go to 8200
      if ( jcol .gt. neqns  ) go to 8200

      if ( unsym ) then                  
          lbound = 1
      else
          lbound = jcol
      end if

      do i = 1, nzcol

          irow = jrowin(i)

          if ( irow .lt. lbound ) go to 8200
          if ( irow .gt. neqns  ) go to 8200
 
      enddo
 
c     ----------------------
c     ... insert the entries
c     ----------------------
 
      invp     = work(qinvp)

      if ( lsame ( matrix(1:1), 'A' ) .or. 
     1     lsame ( matrix(1:1), 'K' )    ) then

          wafile   = work(qwafl1)

          call xdslv3 ( jcol  , nzcol , jrowin, value , unsym, neqns,
     1                  anzero, maxnz , anrecd, maxrec, wafile, 
     1                  work(invp),
     2                  work(acolls), work(arowls), work(adiag ),
     3                  work(anzp  ), work(arcpos), work(arclen),
     4                  work(srtlst), work(srtval), work(colstr),
     5                  error )

          if ( error .eq. -1 ) go to 8900                   
          if ( error .ne.  0 ) go to 8200

          work ( qtemp1 ) = anzero
          work ( qtemp2 ) = anrecd
 
      else if ( lsame ( matrix(1:1), 'B' ) .or. 
     1          lsame ( matrix(1:1), 'M' )      ) then

          if ( bmxtyp .gt. 3 ) go to 8500
 
          if ( bmxtyp .eq. 2 ) then
              if ( irow .ne. jcol ) go to 8500
          end if

          wafile   = work(qwafl3)

          call xdslv3 ( jcol  , nzcol , jrowin, value , unsym, neqns,
     1                  bnzero, maxnz , bnrecd, maxrec, wafile, 
     1                  work(invp),
     2                  work(bcolls), work(browls), work(bdiag ),
     3                  work(bnzp  ), work(brcpos), work(brclen),
     4                  work(srtlst), work(srtval), work(colstr),
     5                  error )

          if ( error .eq. -1 ) go to 8900                   
          if ( error .ne.  0 ) go to 8200

          work ( qtemp6 ) = bnzero
          work ( qtemp7 ) = bnrecd

      else

          go to 8850

      end if

      go to 9000
 
c---------------------------------------------------------------------
 
c     --------------
c     ... error trap
c     --------------
 
c     ---------------------------------
c     ... incorrect processing sequence
c     ---------------------------------
 
 8000 continue
      error = -400
      call hherr ( 3, 'xdslvc', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88000 ) error, stage
      go to 9000
 
c     ---------------------------------------
c     ... insufficient storage for this stage
c     ---------------------------------------
 
 8100 continue
      error = -401
      call hherr ( 2, 'xdslvc', error, wkreqd )
      if ( msglvl .gt. 0 ) write ( output, 88100 ) error, wkreqd, lwork
      go to 9000
 
c     ---------------------
c     ... error from xdslv3
c     ---------------------
 
 8200 continue
      error = -402
      call hherr ( 1, 'xdslvc', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88200 ) error
      work(qstage) = -1
      go to 9000
 
c     -----------------------------------------------
c     ... conflict between declaration of b and input
c     -----------------------------------------------
 
 8500 continue
      error = -408
      call hherr ( 1, 'xdslvc', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88500 ) error
      work(qstage) = -1
      go to 9000
 
c     ----------------------------
c     ... illegal value for matrix
c     ----------------------------
 
 8850 continue
      error = -409
      call hherr ( 1, 'xdslvc', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88850 ) error, matrix
      work(qstage) = -1
      go to 9000
 
c     ----------------------------
c     ... i/o error on file wafile
c     ----------------------------
 
 8900 continue
      error = -410
      call hherr ( 3, 'xdslvc', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88900 ) error, wafile
      work(qstage) = -1
 
c---------------------------------------------------------------------
 
c     ------------------------
c     ... end of module xdslvc
c     ------------------------
 
 9000 continue
      return
 
c---------------------------------------------------------------------
 
c     -----------
c     ... formats
c     -----------
 
68000 format ( /1x, '=================================='
     1         /1x, '= multifrontal value input phase ='
     2         /1x, '==================================' )
 
88000 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslvc executed in an'
     2         /5x, 'incorrect sequence.  current stage = ', i10,
     3          5x, 'should be 30.' )
 
88100 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslvc requires ', i15
     2         /5x, 'words of workspace and has only ', i15,
     3              ' available.' )
 
88200 format ( /5x, '*** fatal error no. ', i5, ' *** value input ',
     1              ' to subroutine xdslvc for'
     2         /5x, 'an entry not input during structural input.' )
 
88500 format ( /5x, '*** fatal error no. ', i5, ' *** value input ',
     1              ' to subroutine xdslvc for'
     2         /5x, 'entries in conflict with type of matrix b.' )

88850 format ( /5x, '*** fatal error no. ', i5, ' *** invalid ',
     1              'value of matrix input to'
     2         /5x, 'subroutine xdslvc.'
     3        /10x, 'matrix = ', a )
 
88900 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslvc encountered i/o error'
     2         /5x, 'on i/o file no. ', i15 )
 
c---------------------------------------------------------------------
 
      end
