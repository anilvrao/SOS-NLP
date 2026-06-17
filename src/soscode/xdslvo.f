      subroutine xdslvo ( unsym , neqns , nsuper, nsnind, nnzero, 
     1                    savemx, msglvl, output, xsup  , xlindx, 
     2                    perm  , diag  , xrelin, sqfile, sqnrec,
     3                    sqflln, sqfltr, sqfil2, sqtrn2, 
     4                    nrecrd, reclen, recpos, wafile, walen ,
     5                    watrn , wafil2, walen2, watrn2,
     6                    inuse , work  , lwork , wkreqd, lindxg,
     7                    lstusd, tnorm , iofile, badrow, badcol, 
     8                    error   )
 
c
c     purpose
c     -------
c
c     xdslvo performs the finalization of matrix value input when
c     the matrix is held on disk.
c
c     created         10-sep-97   -- rgg --
c     modified        14-apr-99   -- rgg -- corrected error handling
c                                           at end of subroutine
c     modified        03-apr-00   -- rgg -- replaced rewind with open
c                                           for sqfile
c
c     ------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      integer             neqns , nsuper, nsnind, nnzero, 
     1                    msglvl, output, sqfile, sqnrec,
     1                    sqfil2, nrecrd, wafile, wafil2,
     2                    inuse , lwork , wkreqd, lindxg, lstusd, 
     3                    iofile, badrow, badcol, error

      integer             xsup   (*),     xlindx (*),     
     1                    perm   (*),     xrelin (*),
     2                    reclen (*),     recpos (*)

      logical             unsym , savemx
 
      double precision    tnorm , sqflln, sqfltr , sqtrn2, walen ,
     1                    watrn , walen2, watrn2 
 
      double precision    diag   (*),     work   (*)
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------
 
      integer             coor21, coor22, coor31, coor32, 
     1                    coord1, coord2,
     2                    ilset , inrecd, insnin, iwaln1, iwaln2,
     3                    iwatr1, iwatr2, kmult , l     , lset  ,   
     4                    maxnnz, recln2, recps2, temp  ,
     5                    values, value2, value3

      double precision    tlen

c     --------------------
c     ... subprograms used
c     --------------------

      integer             xdslni

      external            xdslni, xdslp5, xdslv8, xdslvs, 
     1                    xislp1, xislp2 
 
c---------------------------------------------------------------------

c     ------------------------------
c     ... allocate temporary storage
c     ------------------------------

      inrecd = xdslni ( nrecrd )

      recps2 = inuse  + 1
      recln2 = recps2 + inrecd
      coord1 = recln2 + inrecd

      l      = lwork - coord1 + 1
      kmult  = 4. * ( 64. + 2. * xdslni(64) ) / 64.

      lset  = ( l - kmult + 1 ) / kmult
c.debug
c     lset  = 50
c.debug
      ilset = xdslni ( lset )         

      coord2 = coord1 + ilset
      values = coord2 + ilset
      coor21 = values +  lset
      coor22 = coor21 + ilset
      value2 = coor22 + ilset
      coor31 = value2 +  lset
      coor32 = coor31 + 2 * ilset
      value3 = coor32 + 2 * ilset
c.debug
c     write(6,'("inuse , nrecrd, inrecd = ", 3i8)')
c    1            inuse , nrecrd, inrecd 
c     write(6,'("lset  , ilset          = ", 3i8)')
c    1            lset  , ilset          
c     write(6,'("recps2, recln2, coord1 = ", 3i8)')
c    1            recps2, recln2, coord1 
c     write(6,'("coord2, values, coor21 = ", 3i8)')
c    1            coord2, values, coor21 
c     write(6,'("coor21, coor22, value2 = ", 3i8)')
c    1            coor21, coor22, value2 
c     write(6,'("coor31, coor32, value3 = ", 3i8)')
c    1            coor31, coor32, value3 
c.debug

c.debug
c     if ( value3 + 2 * lset - 1 .gt. lwork ) then
c         write(6,'("storage oops in xdslvo before xdslvs")')
c         stop
c     end if
c.debug

c     ---------------------------------------
c     ... perform the actual sort/merge for a
c     ---------------------------------------

c.debug
c     write(6,'(15x,"before call to xdslvs  - lset, nrecrd = ",2i8)')
c    1                                         lset, nrecrd
c     write(6,'("before xdslvs - wafile, wafil2        = ",3i8)')
c    1                            wafile, wafil2        
c.debug
 
      call xdslvs ( neqns , nsuper, xsup  , xrelin, nrecrd, 
     1              wafile, recpos, reclen, iwaln1, iwatr1,
     2              wafil2, work(recps2), work(recln2),
     3              iwaln2, iwatr2, lset,
     4              work(coord1), work(coord2), work(values), 
     5              work(coor21), work(coor22), work(value2), 
     6              work(coor31), work(coor32), work(value3), 
     7              nnzero, maxnnz, error )

c.debug
c     write(6,'("after  xdslvs - wafile, wafil2, error = ",3i8)')
c    1                            wafile, wafil2, error 
c.debug
c.debug
c     t1 = work ( qvaltm )
c     w1 = work ( qvalwl )
c     call xdslt2 ( t1, w1, t3, w3 )
c     write(6,'(15x,"after xdslvs           - accum. time = ", 
c    1         f15.6)') t3
c.debug
 
      if ( error .ne. 0 ) then
          iofile = wafile
c.debug
c     write(6,'("sqfile i/o - 1 - iofile, error = ",2i8)')
c    1                             iofile, error 
c.debug
          error  = -3
          return 
      end if
 
      tlen         = iwaln1
      walen  = max ( walen, tlen )
      watrn  = watrn + iwatr1
 
      tlen         = iwaln2
      walen2 = max ( walen2, tlen )
      watrn2 = watrn2 + iwatr2

c     ------------------
c     ... optional print
c     ------------------

      if ( msglvl .ge. 3 ) then
          call xdslp5 ( 'diagp',  neqns  , diag  , output )
          call xislp3 ( 'xrelin', neqns+1, xrelin, output )
      end if

c     ------------------
c     ... read in lindxg
c     ------------------

      insnin = xdslni(nsnind)

      if ( sqfil2 .gt. 0 ) then

          lindxg = inuse + 1
          lstusd = lindxg + insnin - 1

          if ( lstusd .gt. lwork ) then
              wkreqd = lstusd + 2*xdslni(maxnnz) + l
              error = -4
              return
          end if

          call xislrw ( sqfil2, error )
          if ( error .ne. 0 ) then
              iofile = sqfil2
c.debug
c     write(6,'("sqfile i/o - 2 - iofile, error = ",2i8)')
c    1                             iofile, error 
c.debug
              error  = -3
              return 
          end if
          
          call xislvr ( sqfil2, nsnind, work(lindxg), error )
          if ( error .ne. 0 ) then
              iofile = sqfil2
c.debug
c     write(6,'("sqfile i/o - 3 - iofile, error = ",2i8)')
c    1                             iofile, error 
c.debug
              error  = -3
              return 
          end if
          
          call xislvr ( sqfil2, nsnind, work(lindxg), error )
          if ( error .ne. 0 ) then
              iofile = sqfil2
c.debug
c     write(6,'("sqfile i/o - 4 - iofile, error = ",2i8)')
c    1                             iofile, error 
c.debug
              error  = -3
              return 
          end if

          sqtrn2 = sqtrn2 + 2 * insnin

      else

          lstusd = inuse

      end if
c.debug
c     call xislp3 ( 'lindxg after read in xdslvo', nsnind, 
c    1              work(lindxg), output )
c.debug

c     -----------------------------------
c     ... reallocate temporary work space
c     -----------------------------------

      temp   = lstusd + 1
      coord1 = temp  + neqns

      l      = lwork - coord1 + 1
      kmult  = ( 64. + 2. * xdslni(64) ) / 64.

      lset  = ( l - kmult + 1 ) / kmult
      ilset = xdslni ( lset )         
 
      coord2 = coord1 + ilset
      values = coord2 + ilset
c.debug
c     write(6,'("inuse , lstusd, lindxg = ", 3i8)')
c    1            inuse , lstusd, lindxg 
c     write(6,'("temp  , coord1, coord2 = ", 3i8)')
c    1            temp  , coord1, coord2 
c     write(6,'("values, lwork          = ", 3i8)')
c    1            values, lwork          
c     write(6,'("kmult , lset  , ilset  = ", 3i8)')
c    1            kmult , lset  , ilset  
c.debug

c.debug
c     if ( values + lset - 1 .gt. lwork ) then
c         write(6,'("storage oops in xdslvf before xdslvs")')
c         stop
c     end if
c.debug

c     ------------------------------------------------------
c     ... check if enough storage is available for next step
c     ------------------------------------------------------

      if ( lset .lt. maxnnz ) then
          l      = maxnnz - lset
          wkreqd = lwork + 2*xdslni(l) + l
c.debug
c     write(6,'("in xdslvf before storage ck pt 02")')
c.debug
          error = -4
          return
      end if

c     --------------------------------------------------
c     ... compute the relative indices for a and prepare
c         matrix data for a for minimum core processing.
c     --------------------------------------------------

      call xislvo ( sqfile, error )
      if ( error .ne. 0 ) then
          iofile = sqfile
c.debug
c     write(6,'("sqfile i/o - 5 - iofile, error = ",2i8)')
c    1                             iofile, error 
c.debug
          error  = -3
          return 
      end if

c.debug
c     t1 = work ( qvaltm )
c     w1 = work ( qvalwl )
c     call xdslt2 ( t1, w1, t3, w3 )
c     write ( output, '("time before xdslv8 = ", 1pd15.5)') t3
c     call xislp3 ( 'xsup   before xdslv8', nsuper+1, xsup  , output )
c     call xislp3 ( 'xlindx before xdslv8', nsuper+1, xlindx, output )
c.debug

      iwatr2 = 0

      call xdslv8 ( unsym , savemx, neqns , nsuper, lset, 
     1              xsup  , xlindx, work(lindxg), 
     2              perm  , work(coord1), work(coord2),
     3              diag  , work(values), xrelin, work(temp  ),
     5              tnorm , badrow, badcol, wafil2, iwatr2, 
     6              sqfile, sqflln, sqfltr, sqnrec, error )
c.debug
c     t1 = work ( qvaltm )
c     w1 = work ( qvalwl )
c     call xdslt2 ( t1, w1, t3, w3 )
c     write(6,'(15x,"after xdslv8           - accum. time = ", 
c    1         f15.6)') t3
c     write(6,'("in xdslvo after  call to xdslv8")')
c     write(6,'("sqflln, sqfltr = ", 2f25.0)') sqflln, sqfltr
c.debug

      if ( error .eq. -4 ) then
          iofile = sqfile
          error = -3
          return
      end if

      if ( error .eq. -3 ) then
          iofile = wafil2
          error = -3
          return
      end if

      if ( error .ne. 0 ) return

      watrn2 = watrn2 + iwatr2
 
c---------------------------------------------------------------------
 
c     ------------------------
c     ... end of module xdslvo
c     ------------------------
 
      return
      end
