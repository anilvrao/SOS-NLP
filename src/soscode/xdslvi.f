      subroutine xdslvi ( unsym , neqns , nsuper, nnzero, mincor,
     1                    savemx, sqfile, nrecrd, sqflln, sqfltr,
     2                    msglvl, output, xsup  , xlindx, lindxg,
     2                    perm  , diag  , xrelin, temp  , work  , 
     3                    collst, rowlst, mtxval, relind,
     4                    lstusd, tnorm , badrow, badcol, error   )
 
c
c     purpose
c     -------
c
c     xdslvi performs the finalization of matrix value input when the
c     matrix is in memory.
c
c     created         10-sep-97   -- rgg --
c     modified        03-apr-00   -- rgg -- replaced rewind with open
c                                           for sqfile
c
c     ------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      integer             neqns , nsuper, nnzero, mincor, sqfile,  
     1                    nrecrd, msglvl, output, collst, rowlst,
     2                    mtxval, relind, lstusd, 
     3                    badrow, badcol, error

      integer             xsup   (*),     xlindx (*),     lindxg (*),
     1                    perm   (*),     xrelin (*)

      logical             unsym , savemx
 
      double precision    tnorm , sqflln, sqfltr 
 
      double precision    diag   (*),     temp   (*),     work(*)
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------
 
      integer             k1,     k2,     oldloc

c     --------------------
c     ... subprograms used
c     --------------------

      integer             xdslni
 
      external            xdslfk, xdslmv, xdslni, xdslp5, xdslv7,
     1                    xislmv, xislp1, xislp2 
 
c---------------------------------------------------------------------

c     -----------------------------------------------
c     ... construct assembly indicies for a with data
c         held in memory.  relind is built on top of
c         collst.
c     ------------------------------------------------
c.debug
c     write(6,'("entering xdslvi - neqns, nnzero = ",2i8)')
c    1        neqns, nnzero
c     call xislp3 ( 'col list', nnzero, work(collst), output )
c     call xislp3 ( 'row list', nnzero, work(rowlst), output )
c     call xdslp5 ( 'diagp', neqns, diag, output )
c     call xdslp5 ( 'mtxval', nnzero, work(mtxval), output )
c.debug

      call xdslv7 ( unsym , neqns , nsuper, nnzero, 
     1              xsup  , xlindx, lindxg, perm  ,
     2              work(collst), work(rowlst), diag, 
     3              work(mtxval),   xrelin, temp  ,
     4              tnorm , badrow, badcol, error )

      if ( error .ne. 0 ) return

      relind = collst
c.debug
c     write(6,'("in xdslvi - msglvl, output, neqns = ",3i8)')
c    1                        msglvl, output, neqns
c.debug
c.debug
c     write(6,'("in xdslvi after xdslv7 - neqns, nnzero = ",2i8)')
c    1        neqns, nnzero
c     call xislp3 ( 'col list', nnzero, work(collst), output )
c     call xislp3 ( 'row list', nnzero, work(rowlst), output )
c     call xdslp5 ( 'diagp', neqns, diag, output )
c     call xdslp5 ( 'mtxval', nnzero, work(mtxval), output )
c.debug

      if ( msglvl .ge. 3 ) then
          call xdslp5 ( 'diagp', neqns, diag, output )
      end if
    
      if ( msglvl .eq. 3 ) then
          call xislp3 ( 'xrelin', neqns+1, xrelin, output )
      else if ( msglvl .ge. 4 ) then
          call xislp1 ( 'relative indicies of permuted matrix', 
     1                  neqns, xrelin, work(relind), output )
          call xislp1 ( 'row indicies of permuted matrix', 
     1                  neqns, xrelin, work(rowlst), output )
      end if

      if ( msglvl .ge. 4 ) then
          call xdslp5 ( 'mtxval', nnzero, work(mtxval), output )
      end if
         
c     --------------------------------------------------------
c     ... prepare for data of matrix a for numeric phases.
c         depends of mincor settings for a and whether a
c         is to be preserved for matrix multiplication.
c     --------------------------------------------------------

      if ( mincor .eq. 0 ) then

c         ----------------------------------------
c         ... if not preserving a the purge rowlst
c             and move matrix values for a.
c         ----------------------------------------

          if ( .not. savemx ) then 

              oldloc = mtxval
              mtxval = rowlst
    
              call xdslmv ( nnzero, work, oldloc, mtxval )

          end if

          lstusd = mtxval + nnzero - 1 

c         ------------------
c         ... optional print
c         ------------------

          if ( msglvl .ge. 3 ) then
              call xdslp5 ( 'diagp', neqns, diag, output )
          end if
    
          if ( msglvl .eq. 3 ) then
              call xislp3 ( 'xrelin', neqns+1, xrelin, output )
          end if
 
          if ( msglvl .ge. 4 ) then
              if ( savemx ) 
     1        call xislp1 ( 'row indices of permuted matrix',
     2                      neqns, xrelin, work(rowlst), output )
              call xislp1 ( 'relative indices of permuted matrix',
     1                      neqns, xrelin, work(relind), output )
              call xdslp5 ( 'mtxval', nnzero, work(mtxval), output )
          end if

      else

c         ----------------------------------------------------
c         ... minimum core activated for the matrix.  so write
c             original matrix information to sqfile.
c         ----------------------------------------------------
c.debug
c     write(6,'("in xdslvi before call to xdslfk")')
c.debug

          if ( sqfile .le. 0 ) then
              error = -3
              return
          end if

          call xislvo ( sqfile, error )
          if ( error .ne. 0 ) then      
              error = -3
              return
          end if

          call xdslfk ( sqfile, unsym, savemx, nsuper, xsup,
     1                  xrelin, work(rowlst), work(relind),
     2                  work(mtxval), nrecrd, error )
c.debug
c     write(6,'("in xdslvi after  call to xdslfk")')
c.debug
          if ( error .ne. 0 ) then      
              error = -3
              return
          end if

          k1 = 1
          k2 = 1
          if ( savemx ) k1 = 2
          if ( unsym  ) k2 = 2

          sqflln = sqflln + k1*xdslni(nnzero) + k2*nnzero
          sqfltr = sqfltr + k1*xdslni(nnzero) + k2*nnzero
c.debug
c     write(6,'("in xdslvi after  call to xdslfk")')
c     write(6,'("sqflln, sqfltr = ", 2f25.0)') sqflln, sqfltr
c.debug

          lstusd = relind - 1
          relind = 0
          rowlst = 0
          mtxval = 0

c         ------------------
c         ... optional print
c         ------------------

          if ( msglvl .ge. 3 ) then
              call xdslp5 ( 'diagp' , neqns  , diag  , output )
              call xislp3 ( 'xrelin', neqns+1, xrelin, output )
          end if

c         ---------------------------------------------------
c         ... end of processing for in-core representation of 
c             original matrix
c         ---------------------------------------------------

      end if
 
c---------------------------------------------------------------------
 
c     ------------------------
c     ... end of module xdslvi
c     ------------------------
 
      return
      end
