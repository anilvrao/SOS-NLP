      subroutine xdslvf ( work,   lwork,  error )
c
c     purpose
c     -------
c
c     xdslvf is the top level driver for finalization of matrix value
c     input.
c
c     created         08-jan-97   -- rgg --
c     modified        01-nov-01   -- dkw -- changed error -402 to -404
c                                           and test for invalid stage
c                                           error = -450
c
c     input arguments
c     ---------------
c
c     lwork       i   length of work array.
c
c     input/output arguments
c     ----------------------
c
c     work        d   work array.
c
c     output arguments
c     ----------------
c
c     error       i   error flag
c                     =    0  normal return
c                     = -401  lwork not large enough
c                     = -403  internal data structure error
c                     = -404  error in matrix input
c                     = -410  i/o error on wafil1, wafil2, 
c                             sqfil1, or sqfil2
c                     = -450  incorrect processing path
c
c     ------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      integer             lwork,  error
 
      double precision    work(*)
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------
 
      integer             acolls, anrecd, anzero, anzp  , arindx,
     1                    arowls, adiag , ianzer, xarndx,
     1                    arcpos, arclen, amncor,
     1                    bcolls, bnrecd, bnzero, bnzp  , brindx,
     1                    browls, bdiag , ibnzer, xbrndx,
     1                    brcpos, brclen, bmncor,
     1                    ineqns, ineqn1, 
     2                    inuse , lindxg, maxnz , msglvl, mxused, 
     3                    maxwrk, neqns , nsuper, output, 
     4                    stage,  temp  , temp2 , wafil1, wafil2, 
     5                    wkreqd, xlindx, xsup  , nsnind, insnin

      integer             badcol, badrow, ianrec, ibnrec,
     1                    imaxnz, iofile, l,      lstusd,
     3                    mincor, perm  , sqfil1, sqfil2, sqfil3,
     4                    nrecrd, colstr, srtlst, srtval

      integer             bmxtyp, imxrec, maxrec, wafil3

      integer             idummy(1)

      logical             unsym , savea , saveb , qincor

 
      double precision    tnorma, tnormb, t1,     t2,     t3,
     1                    w1,     w2,     w3
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      integer             xdslil, xdslni
 
      external            icopy,  xdslvs, xdslvz, xdslvi, xdslvo,
     1                    xdslil, xdslni, xdslt2, xislp1, xislp2,
     2                    xdslmv, xislmv
 
c---------------------------------------------------------------------

c     -----------------------------------------------
c     ... extract data from global.CMNication area
c     -----------------------------------------------
 
      error  = 0
 
      stage  = work ( qstage )
      msglvl = work ( qmsglv )
      output = work ( qoutpu )
      if ( stage .ne. 31 ) go to 8000

      unsym  = work ( qmxtyp ) .eq. 2.
      bmxtyp = work ( qbmxty )
 
      neqns  = work ( qneqns )
      ineqns = xdslni ( neqns )
      ineqn1 = xdslni ( neqns + 1 )

      maxrec = min ( neqns, mxrecd )
      imxrec = xdslni ( maxrec )

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

      nsnind = work ( qnsnin )
      nsuper = work ( qnsupe )
      mincor = work ( qmncor ) 

      amncor = mod ( mincor, 10 ) 
      bmncor = mincor / 10

      insnin = xdslni ( nsnind )

      savea  = .false.
      if ( work(qsavea) .ne. 0. ) savea = .true.

      saveb  = .true.

      sqfil1 = work ( qsqfl1 )
      sqfil2 = work ( qsqfl2 )
      sqfil3 = work ( qsqfl3 )

c.debug
c         write(6,'("in xdslvf")')
c         write(6,'("adiag,  arcpos, arclen = ", 3i8)')
c    1                adiag,  arcpos, arclen
c         write(6,'("acolls, arowls, anzp   = ", 3i8)')
c    1                acolls, arowls, anzp 
c         write(6,'("bdiag,  brcpos, brclen = ", 3i8)')
c    1                bdiag,  brcpos, brclen
c         write(6,'("bcolls, browls, bnzp   = ", 3i8)')
c    1                bcolls, browls, bnzp 
c         write(6,'("mincor, amncor, bmncor = ", 3i8)')
c    1                mincor, amncor, bmncor
c.debug

      t1 = work ( qvaltm )
      w1 = work ( qvalwl )
      call xdslt2 ( t1, w1, t2, w2 )

c.debug
c     write(6,'(15x,"at start of xdslvf     - accum. time = ", 
c    1         f15.6)') t2
c.debug
 
c     --------------------------------
c     ... sort and compress data for a
c     --------------------------------

      wafil1 = work(qwafl1)

c.debug
c     write(6,'(/15x, "in xdslvf before xdslvz")')
c     write(6,'(15x,"anzero, maxnz  = ", 2i8)') anzero, maxnz
c     write(6,'(15x,"anrecd, wafil1 = ", 2i8)') anrecd, wafil1
c     write(6,'(15x,"arowls, acolls = ", 2i8)') arowls, acolls
c     write(6,'(15x,"anzp           = ", 2i8)') anzp
c     write(6,'(15x,"arcpos, arclen = ", 2i8)') arcpos, arclen
c     write(6,'(15x,"srtlst, srtval = ", 2i8)') srtlst, srtval
c     write(6,'(15x,"colstr         = ", 2i8)') colstr
c
c     call xislp3 ( 'coord1 list', anzero, work(acolls), output )
c     call xislp3 ( 'coord2 list', anzero, work(arowls), output )
c     call xdslp5 ( 'values list', anzero, work(anzp  ), output )
c.debug

      call xdslvz ( .true., neqns, anzero, maxnz, anrecd, maxrec, 
     1              wafil1, work(acolls), work(arowls), work(anzp),
     2              work(arcpos), work(arclen), 
     3              work(srtlst), work(srtval), work(colstr),
     4              error )

      if ( error .ne. 0 ) then
          iofile = wafil1
c.debug
c         write(6,'("in xdslvf post xdslvz-iofile, error = ", 2i8)')
c    1                                      iofile, error
c.debug
          go to 8900
      end if
 
c     --------------------------------
c     ... sort and compress data for b
c     --------------------------------

      if ( bmxtyp .eq. 1 ) then 

          wafil3 = work(qwafl3)

c.debug
c     write(6,'(/15x, "in xdslvf before xdslvz")')
c     write(6,'(15x,"bnzero, maxnz  = ", 2i8)') bnzero, maxnz
c     write(6,'(15x,"bnrecd, wafil3 = ", 2i8)') bnrecd, wafil3
c     write(6,'("browls, bcolls = ", 2i8)') browls, bcolls
c     write(6,'("bnzp           = ", 2i8)') bnzp
c
c     call xislp3 ( 'coord1 list', bnzero, work(bcolls), output )
c     call xislp3 ( 'coord2 list', bnzero, work(browls), output )
c     call xdslp5 ( 'values list', bnzero, work(bnzp  ), output )
c.debug

          call xdslvz ( .true., neqns, bnzero, maxnz, bnrecd, maxrec, 
     1                  wafil3, work(bcolls), work(browls), work(bnzp),
     2                  work(brcpos), work(brclen), 
     3                  work(srtlst), work(srtval), work(colstr),
     4                  error )

          if ( error .ne. 0 ) then
              iofile = wafil3
c.debug
c         write(6,'("in xdslvf post xdslvz-iofile, error = ", 2i8)')
c    1                                      iofile, error
c.debug
              go to 8900
          end if

      end if

c     -------------------------------------------------
c     ... test to see if all data can be kept in memory
c     -------------------------------------------------

      ianzer = xdslni ( anzero )
      ibnzer = xdslni ( bnzero )
      
      l = 2*ianzer + anzero + 2*ibnzer + bnzero + neqns
      if ( sqfil2 .gt. 0 ) l = l + insnin
      if ( bmxtyp .eq. 1 ) l = l + ineqn1

      inuse = adiag + neqns - 1
      if ( bmxtyp .le. 3 ) inuse = bdiag + neqns - 1

      if ( inuse + l .le. lwork .and.
     1     anrecd .eq. 0        .and.
     2     bnrecd .eq. 0              ) then
          qincor = .true.
      else
          qincor = .false.
      end if
c.debug
c     write(6,'("in xdslvf - qincor = ", l8)') qincor
c.debug

c     --------------------------------------------------------
c     ... if out-of-memory, write a and b to disk as necessary
c     --------------------------------------------------------

      if ( .not. qincor .and. anrecd .eq. 0 ) then

c         ---------------
c         ... write out a
c         ---------------

          iofile = wafil1

          call xdslw4 ( wafil1, 3, error )
          if ( error .ne. 0 ) go to 8900

          call xdslw2 ( wafil1, 3, work(acolls), work(arowls), 
     1                  work(anzp), 1, anzero, error )
          if ( error .ne. 0 ) go to 8900

          anrecd = 1

          idummy(1) = 1
          call icopy ( 1, idummy, 1, work(arcpos), 1 )

          idummy(1) = anzero
          call icopy ( 1, idummy, 1, work(arclen), 1 )

          anzero = 0

      end if

      if ( .not. qincor .and. bnrecd .eq. 0 .and. bmxtyp .eq. 1 ) then

c         ---------------
c         ... write out b
c         ---------------

          iofile = wafil3

          call xdslw4 ( wafil3, 3, error )
          if ( error .ne. 0 ) go to 8900

          call xdslw2 ( wafil3, 3, work(bcolls), work(browls), 
     1                  work(bnzp), 1, bnzero, error )
          if ( error .ne. 0 ) go to 8900

          bnrecd = 1

          idummy(1) = 1
          call icopy ( 1, idummy, 1, work(brcpos), 1 )
    
          idummy(1) = bnzero
          call icopy ( 1, idummy, 1, work(brclen), 1 )

          bnzero = 0

       end if

c.debug
c     write(6,'(/ "in xdslvf after  xdslvz")')
c     write(6,'("anzero, maxnz  = ", 2i8)') anzero, maxnz
c     write(6,'("anrecd, wafil1 = ", 2i8)') anrecd, wafil1
c     write(6,'("arowls, acolls = ", 2i8)') arowls, acolls
c     write(6,'("anzp           = ", 2i8)') anzp
c     write(6,'("bnzero, maxnz  = ", 2i8)') bnzero, maxnz
c     write(6,'("bnrecd, wafil3 = ", 2i8)') bnrecd, wafil3
c     write(6,'("browls, bcolls = ", 2i8)') browls, bcolls
c     write(6,'("bnzp           = ", 2i8)') bnzp
c     call xislp3 ( 'arcpos', anrecd, work(arcpos), output )
c     call xislp3 ( 'arclen', anrecd, work(arclen), output )
c     call xislp3 ( 'brcpos', bnrecd, work(brcpos), output )
c     call xislp3 ( 'brclen', bnrecd, work(brclen), output )
c
c     call xislp3 ( 'a-coord1 list', anzero, work(acolls), output )
c     call xislp3 ( 'a-coord2 list', anzero, work(arowls), output )
c     call xdslp5 ( 'a-values list', anzero, work(anzp  ), output )
c     call xislp3 ( 'b-coord1 list', bnzero, work(bcolls), output )
c     call xislp3 ( 'b-coord2 list', bnzero, work(browls), output )
c     call xdslp5 ( 'b-values list', bnzero, work(bnzp  ), output )
c.debug
c.debug
c     t1 = work ( qvaltm )
c     w1 = work ( qvalwl )
c     call xdslt2 ( t1, w1, t3, w3 )
c     write(6,'(15x,"after xdslvz           - accum. time = ", 
c    1         f15.6)') t3
c.debug

c     ------------------------------------------------------------------
 
c     -------------------------------------------------
c     ... remainder of processing depends on whether
c         the coord1/coord2 structure is in memory or 
c         not.
c         in either case the result is the relative 
c         indices required for the matrix assembly 
c         step of the numerical factorization.
c     -------------------------------------------------

      xarndx = work ( qxadj  )
      perm   = work ( qperm  )
      xsup   = work ( qxsup  )
      xlindx = work ( qxlind )

      inuse = adiag + neqns - 1
      if ( bmxtyp .le. 3 ) inuse = bdiag + neqns - 1
c.debug
c     write(6,'("after assignment of adiag and bdiag")')
c     write(6,'("inuse = ", i8)') inuse 
c.debug

      work ( qtemp6 ) = 0
      work ( qtemp7 ) = 0

c     ------------------------------------------------------------------

      if ( qincor ) then

c         ----------------------------------------------
c         ... matrix data for a and b are held in-memory
c         ----------------------------------------------
c.debug
c     write(6,'("in xdslvf - processing matrices in memory", )')
c.debug

c         --------------------------
c         ... compress storage for a
c         --------------------------

          temp   = acolls
          acolls = inuse + 1
          call xislmv ( anzero, work, xdslil(temp-1)+1,
     1                                xdslil(acolls-1)+1 )

          temp   = arowls
          arowls = acolls + ianzer
          call xislmv ( anzero, work, xdslil(temp-1)+1,
     1                                xdslil(arowls-1)+1 )

          temp   = anzp  
          anzp   = arowls + ianzer
          call xdslmv ( anzero, work, temp, anzp )
 
          inuse  = anzp + anzero - 1

c         --------------------------
c         ... compress storage for b
c         --------------------------

          if ( bmxtyp .eq. 1 ) then

              temp   = bcolls
              bcolls = inuse + 1
              call xislmv ( bnzero, work, xdslil(temp-1)+1,
     1                                    xdslil(bcolls-1)+1 )

              temp   = browls
              browls = bcolls + ibnzer
              call xislmv ( bnzero, work, xdslil(temp-1)+1,
     1                                    xdslil(browls-1)+1 )

              temp   = bnzp  
              bnzp   = browls + ibnzer
              call xdslmv ( bnzero, work, temp, bnzp )

              inuse  = bnzp + bnzero - 1

          end if
c.debug
c     write(6,'(/ "in xdslvf after  compression")')
c     write(6,'("anzero         = ", 2i8)') anzero
c     write(6,'("arowls, acolls = ", 2i8)') arowls, acolls
c     write(6,'("anzp           = ", 2i8)') anzp
c     write(6,'("bnzero         = ", 2i8)') bnzero 
c     write(6,'("browls, bcolls = ", 2i8)') browls, bcolls
c     write(6,'("bnzp           = ", 2i8)') bnzp
c     write(6,'("after compression")')
c     write(6,'("inuse = ", i8)') inuse 
c.debug
 
c         ------------------------------
c         ... read in lindxg from sqfil2 
c         ------------------------------

          if ( sqfil2 .gt. 0 ) then

c.debug
c     write(6,'("need to read lindxg in")')
c     write(6,'("inuse , insnin = ", 2i8)') inuse , insnin
c.debug

              lindxg = inuse + 1
              inuse  = lindxg + insnin - 1

              if ( nsnind .gt. 0 ) then

                  call xislrw ( sqfil2, error )
                  if ( error .ne. 0 ) go to 8900
    
                  iofile = sqfil2
              
                  call xislvr ( sqfil2, nsnind, work(lindxg), error )
                  if ( error .ne. 0 ) go to 8900
              
                  call xislvr ( sqfil2, nsnind, work(lindxg), error )
                  if ( error .ne. 0 ) go to 8900

                  work(qsqtr2) = work(qsqtr2) + 2 * insnin

              end if

          else

              lindxg = work ( qlndxg )
c.debug
c     write(6,'("lindxg is in memory")')
c     write(6,'("lindxg         = ", 2i8)') lindxg        
c.debug

          end if

c         ------------------------------------------
c         ... prepare to construct assembly indicies 
c             for a and b
c         ------------------------------------------

          if ( bmxtyp .eq. 1 ) then
              xbrndx = inuse + 1
              inuse  = xbrndx + ineqn1
          end if

          temp2  = inuse + 1
          inuse  = temp2 + neqns - 1

          wkreqd = inuse
          mxused = inuse

          if ( wkreqd .gt. lwork ) go to 8100
c.debug
c     write(6,'("lindxg, xbrndx = ", 2i8)') lindxg, xbrndx
c     write(6,'("temp2 , inuse  = ", 2i8)') temp2 , inuse 
c.debug

c         --------------------------------------------
c         ... perform in-core finalization of data for
c             matrix a
c         --------------------------------------------

          nrecrd = 0
          iofile = sqfil1
c.debug
c     write(6,'("before xdslvi for a")')
c     call xislp3 ( 'b-coord1 list', bnzero, work(bcolls), output )
c.debug

          call xdslvi ( unsym , neqns , nsuper, anzero, amncor,
     1                  savea , sqfil1, nrecrd, work(qsqln1),   
     2                  work(qsqtr1), msglvl, output,
     3                  work(xsup),     work(xlindx),
     4                  work(lindxg),   work(perm),     work(adiag),
     5                  work(xarndx),   work(temp2),    work,
     6                  acolls, arowls, anzp  , arindx,
     7                  lstusd, tnorma, badrow, badcol, error )
c.debug
c     write(6,'("after  xdslvi for a")')
c     write(6,'("error          = ", 2i8)') error 
c     write(6,'("anzero         = ", 2i8)') anzero
c     write(6,'("arowls, acolls = ", 2i8)') arowls, acolls
c     write(6,'("anzp,   arindx = ", 2i8)') anzp,   arindx
c     write(6,'("bnzero         = ", 2i8)') bnzero 
c     write(6,'("browls, bcolls = ", 2i8)') browls, bcolls
c     write(6,'("bnzp           = ", 2i8)') bnzp
c     write(6,'("lstusd         = ", 2i8)') lstusd
c     write(6,'("setting qtemp6 = ", 2i8)') nrecrd
c     call xislp3 ( 'b-coord1 list', bnzero, work(bcolls), output )
c.debug

          if ( error .eq. -1 ) go to 8200
          if ( error .eq. -2 ) go to 8300
          if ( error .eq. -3 ) go to 8900

          work(qtemp6) = nrecrd

c         ------------------------------------------
c         ... compress storage for b as appropriate
c             b is a general matrix, will be in-core
c             and there is gas to be compressed out.
c         ------------------------------------------

          if ( bmxtyp .eq. 1 .and. 
     1         bmncor .eq. 0 .and.
     2         bcolls .ne. lstusd + 1 ) then

              temp   = bcolls
              bcolls = lstusd + 1
              call xislmv ( bnzero, work, xdslil(temp-1)+1,
     1                                    xdslil(bcolls-1)+1 )

              temp   = browls
              browls = bcolls + ibnzer
              call xislmv ( bnzero, work, xdslil(temp-1)+1,
     1                                    xdslil(browls-1)+1 )

              temp   = bnzp  
              bnzp   = browls + ibnzer
              call xdslmv ( bnzero, work, temp, bnzp )

              lstusd = bnzp   + bnzero - 1

          end if
         
          if ( bmxtyp .eq. 1 ) then 

c             --------------------------------------------
c             ... perform in-core finalization of data for
c                 matrix b
c             --------------------------------------------

              nrecrd = 0 
              iofile = sqfil3
 
              if ( bmncor .ne. 0 ) then
                  call xislvo ( sqfil3, error )
                  if ( error .ne. 0 ) go to 8900
              end if
c.debug
c.debug
c     write(6,'("before xdslvi for b")')
c     write(6,'("xbrndx, bdiag , bcolls, browls, bnzp   = ", 5i8)')
c    1            xbrndx, bdiag , bcolls, browls, bnzp  
c     write(6,'("bnzero, bmncor, bmxtyp, xlindx, lindxg = ", 5i8)')
c    1            bnzero, bmncor, bmxtyp, xlindx, lindxg
c     call xislp3 ( 'b-coord1 list', bnzero, work(bcolls), output )
c     call xislp3 ( 'b-coord2 list', bnzero, work(browls), output )
c     call xdslp5 ( 'b-values list', bnzero, work(bnzp  ), output )
c     call xislp3 ( 'xlindx        ', nsuper+1, work(xlindx), output )
c     call xislp3 ( 'global indices', nsnind, work(lindxg), output )
c.debug
 
              call xdslvi ( unsym , neqns , nsuper, bnzero, bmncor,
     1                  saveb , sqfil3, nrecrd, work(qsqln3),
     2                  work(qsqtr3), msglvl, output,
     3                  work(xsup),     work(xlindx),
     4                  work(lindxg),   work(perm),     work(bdiag),
     5                  work(xbrndx),   work(temp2),    work,
     6                  bcolls, browls, bnzp  , brindx,
     7                  lstusd, tnormb, badrow, badcol, error )

              if ( bmncor .ne. 0 .and. amncor .ne. 0 ) then
                  lstusd = acolls - 1
              end if

c.debug
c     write(6,'("after  xdslvi for b")')
c     write(6,'("xbrndx, bdiag , bcolls, browls, bnzp   = ", 5i8)')
c    1            xbrndx, bdiag , bcolls, browls, bnzp  
c     write(6,'("bnzero, bmncor, bmxtyp, brindx, error  = ", 5i8)')
c    1            bnzero, bmncor, bmxtyp, brindx, error
c     write(6,'("setting qtemp7 = ", 2i8)') nrecrd
c.debug

              if ( error .eq. -1 ) go to 8200
              if ( error .eq. -2 ) go to 8300
              if ( error .eq. -3 ) go to 8900

              work(qtemp7) = nrecrd

          end if

c         ------------------
c         ... optional print
c         ------------------

          if ( bmxtyp .eq. 3 .and. msglvl .ge. 3 ) then

              call xdslp5 ( 'diagp for b', neqns, work(bdiag), 
     1                      output )

          end if

c         ------------------------------------------
c         ... compress the remainder of the storage 
c             as appropriate.
c         ------------------------------------------

          if ( sqfil2 .gt. 0 ) then

              temp   = lindxg
              lindxg = lstusd + 1
              lstusd = lindxg + xdslni(nsnind) - 1

              if ( lindxg .ne. temp ) then

                  call xislmv ( nsnind, work, xdslil(temp  -1)+1,
     1                                        xdslil(lindxg-1)+1 )

                  if ( msglvl .ge. 3 ) then
                      call xislp3 ( 'global indices', nsnind, 
     1                              work(lindxg), output )
                  end if

              end if

          end if
c.debug
c     write(6,'("after  setting up for lindxg")')
c     write(6,'("lindxg, nsnind, lstusd                 = ", 5i8)')
c    1            lindxg, nsnind, lstusd  
c.debug

          if ( bmxtyp .eq. 1 ) then

              if ( xbrndx .ne. lstusd + 1 ) then

                  temp   = xbrndx
                  xbrndx = lstusd + 1
                  lstusd = lstusd + ineqn1

                  call xislmv ( neqns+1, work, xdslil(temp-1)+1,
     1                                         xdslil(xbrndx-1)+1 )

                  if ( msglvl .ge. 3 ) then
                      call xislp3 ( 'xbrndx', neqns+1, 
     1                              work(xbrndx), output )
                  end if

              else 

                  lstusd = xbrndx + ineqn1 - 1

              end if

          end if

          inuse  = lstusd
c.debug
c     write(6,'("after  setting up for xbrndx")')
c     write(6,'("xbrndx, neqns , inuse                  = ", 5i8)')
c    1            xbrndx, neqns , inuse   
c.debug

c         ---------------------------------------------------
c         ... end of processing for in-core representation of 
c             original matrix
c         ---------------------------------------------------

      else

c         ---------------------------------------------------
c         ... matrix data for a and b are held out-of-memory
c             note that the last calls to xdslvz flushed
c             all data to wafil1 and wafil3.
c             work(arcpos) holds record position info for a
c             work(arclen) holds record length   info for a
c             work(brcpos) holds record position info for b
c             work(brclen) holds record length   info for b
c         --------------------------------------------------
          lindxg = work ( qlndxg )
          if ( sqfil2 .gt. 0 ) lindxg = 0
c.debug
c     write(6,'("in xdslvf - processing matrices out of memory", )')
c     write(6,'("lindxg                                 = ", 5i8)')
c    1            lindxg
c.debug

          mxused = lwork

          ianrec = xdslni ( anrecd )
          ibnrec = xdslni ( bnrecd )

          temp   = arclen
          arclen = arcpos + ianrec
          call xislmv ( anrecd, work, xdslil(temp  -1)+1,
     1                                xdslil(arclen-1)+1 )

          temp   = brcpos
          brcpos = arclen + ianrec
          call xislmv ( bnrecd, work, xdslil(temp  -1)+1,
     1                                xdslil(brcpos-1)+1 )

          temp   = brclen
          brclen = brcpos + ibnrec
          call xislmv ( bnrecd, work, xdslil(temp  -1)+1,
     1                                xdslil(brclen-1)+1 )

          inuse  = brclen + ibnrec - 1

c         -----------------------------------------------
c         ... prepare of out-of-core processing of matrix
c             data
c         -----------------------------------------------

          wafil2 = work(qwafl2)

          call xdslw4 ( wafil2, 3, error )

          if ( error .ne. 0 ) then
              iofile = wafil2
              go to 8900
          end if

c         --------------------------------------------------------
c         ... complete processing of data for matrix a.
c             finalize sort/merge, compute assembly indices,
c             and place data on sqfil1.
c         --------------------------------------------------------

          nrecrd = 0

          call xdslvo ( unsym , neqns , nsuper, nsnind, anzero,
     1                  savea , msglvl, output, work(xsup),
     2                  work(xlindx), 
     3                  work(perm), work(adiag), work(xarndx),
     4                  sqfil1, nrecrd, work(qsqln1), work(qsqtr1), 
     5                  sqfil2, work(qsqtr2), 
     6                  anrecd, work(arclen), work(arcpos),
     7                  wafil1, work(qwaln1), work(qwatr1),
     8                  wafil2, work(qwaln2), work(qwatr2),
     9                  inuse , work  , lwork , wkreqd, lindxg,
     a                  lstusd, tnorma, 
     b                  iofile, badrow, badcol, error )
c.debug
c     write(6,'("after  xdslvo for a - error = ",i8)') error
c.debug
          if ( error .eq. -1 ) go to 8200
          if ( error .eq. -2 ) go to 8300
          if ( error .eq. -3 ) go to 8900
          if ( error .eq. -4 ) go to 8100

          arindx = 0
          acolls = 0
          arowls = 0
          anzp   = 0

c.debug
c     write(6,'("after  xdslvo for a")')
c     write(6,'("anzero         = ", 2i8)') anzero
c     write(6,'("arowls, acolls = ", 2i8)') arowls, acolls
c     write(6,'("anzp,   arindx = ", 2i8)') anzp,   arindx
c     write(6,'("bnzero         = ", 2i8)') bnzero 
c     write(6,'("setting qtemp6 = ", 2i8)') nrecrd
c.debug

          work(qtemp6) = nrecrd

c         --------------------------------------------------------
c         ... complete processing of data for matrix b.
c             finalize sort/merge, compute assembly indices,
c             and place data on sqfil3.
c         --------------------------------------------------------

          if ( bmxtyp .eq. 1 ) then 

c             --------------------------------------------
c             ... perform in-core finalization of data for
c                 matrix b
c             --------------------------------------------

              nrecrd = 0

              xbrndx = inuse + 1
              inuse  = xbrndx + ineqn1 - 1

              iofile = sqfil3
              call xislvo ( sqfil3, error )
              if ( error .ne. 0 ) go to 8900

              call xdslvo ( unsym , neqns , nsuper, nsnind, bnzero,
     1                  saveb , msglvl, output, work(xsup),
     1                  work(xlindx), 
     2                  work(perm), work(bdiag), work(xbrndx),
     3                  sqfil3, nrecrd, work(qsqln3), work(qsqtr3), 
     4                  sqfil2, work(qsqtr2), 
     5                  bnrecd, work(brclen), work(brcpos),
     6                  wafil3, work(qwaln3), work(qwatr3),
     7                  wafil2, work(qwaln2), work(qwatr2),
     8                  inuse , work  , lwork , wkreqd, lindxg,
     9                  lstusd, tnormb, 
     a                  iofile, badrow, badcol, error )
c.debug
c     write(6,'("after  xdslvo for b - error =",i8)') error
c.debug

              if ( error .eq. -1 ) go to 8200
              if ( error .eq. -2 ) go to 8300
              if ( error .eq. -3 ) go to 8900
              if ( error .eq. -4 ) go to 8100

              brindx = 0
              bcolls = 0
              browls = 0
              bnzp   = 0

c.debug
c     write(6,'("after  xdslvo for b")')
c     write(6,'("xbrndx, bdiag , bcolls, browls, bnzp   = ", 5i8)')
c    1            xbrndx, bdiag , bcolls, browls, bnzp  
c     write(6,'("bnzero, bmncor, bmxtyp                 = ", 5i8)')
c    1            bnzero, bmncor, bmxtyp
c     write(6,'("setting qtemp7 = ", 2i8)') nrecrd
c.debug

              work(qtemp7) = nrecrd

              temp   = xbrndx
              xbrndx = arcpos
              call xislmv ( neqns+1, work, xdslil(temp  -1)+1,
     1                                     xdslil(xbrndx-1)+1 )

              inuse  = xbrndx + ineqn1 - 1

          else
        
              inuse  = arcpos - 1

          end if

c         ----------------------------------------
c         ... optional print for diagonal matrices
c         ----------------------------------------

          if ( bmxtyp .eq. 3 .and. msglvl .ge. 3 ) then

              call xdslp5 ( 'diagp for b', neqns, work(bdiag), 
     1                      output )

          end if

c         ------------------------------
c         ... read in lindxg from sqfil2 
c         ------------------------------

          if ( sqfil2 .gt. 0 ) then

              lindxg = inuse + 1
              inuse  = lindxg + insnin - 1

              call xislrw ( sqfil2, error )
              if ( error .ne. 0 ) go to 8900

              iofile = sqfil2
          
              call xislvr ( sqfil2, nsnind, work(lindxg), error )
              if ( error .ne. 0 ) go to 8900
          
              call xislvr ( sqfil2, nsnind, work(lindxg), error )
              if ( error .ne. 0 ) go to 8900

              work(qsqtr2) = work(qsqtr2) + 2 * insnin

          else

              lindxg = work ( qlndxg )

          end if

c         ----------------------------------------------
c         ... prepare to bring matrices back into memory
c             according to settings of mincor
c         ----------------------------------------------

          wkreqd = inuse
    
          if ( amncor .eq. 0 ) then
              ianzer = xdslni(anzero)
              wkreqd = wkreqd + anzero + ianzer
              if ( savea ) wkreqd = wkreqd + ianzer
          end if

          if ( bmxtyp .eq. 1 .and. bmncor .eq. 0 ) then
              ibnzer = xdslni ( bnzero )
              wkreqd = wkreqd + bnzero + 2*ibnzer
          end if

          if ( wkreqd .gt. lwork ) go to 8100

c         ---------------------------------------------
c         ... bring matrices back into memory according
c             to settings of mincor.
c         ---------------------------------------------

          if ( amncor .eq. 0 ) then

              arindx = inuse + 1

              if ( savea ) then
                  arowls = arindx + ianzer
              else
                  arowls = arindx
              end if

              anzp   = arowls + ianzer

              inuse  = anzp + anzero - 1

              iofile = sqfil1
c.debug
c     write(6,'("before restoring a")')
c     write(6,'("arindx, arowls, anzp  , anzero, inuse  = ", 5i8)')
c    1            arindx, arowls, anzp  , anzero, inuse 
c     write(6,'("savea                                  = ", 5l8)')
c    1            savea
c     call xislp3 ( 'xarndx', neqns+1, work(xarndx), output )
c.debug

              call xdslvg ( unsym , savea , nsuper, anzero, 
     1                      work(xsup  ), work(xarndx), work(arindx),
     2                      work(arowls), work(anzp  ),
     3                      sqfil1, work(qsqtr1), error )
c.debug
c     write(6,'("after restoring a - error =",i8)') error
c.debug

              if ( error .ne. 0 ) go to 8900

              work(qtemp6) = 0

          end if

          if ( bmxtyp .eq. 1 .and. bmncor .eq. 0 ) then

              brindx = inuse + 1
              browls = brindx + ibnzer
              bnzp   = browls + ibnzer

              inuse  = bnzp + bnzero - 1

              iofile = sqfil3
c.debug
c     write(6,'("before restoring b")')
c     write(6,'("brindx, browls, bnzp  , bnzero, inuse  = ", 5i8)')
c    1            brindx, browls, bnzp  , bnzero, inuse 
c     write(6,'("saveb                                  = ", 5l8)')
c    1            saveb
c     call xislp3 ( 'xbrndx', neqns+1, work(xbrndx), output )
c.debug

              call xdslvg ( unsym , saveb , nsuper, bnzero, 
     1                      work(xsup  ), work(xbrndx), work(brindx),
     2                      work(browls), work(bnzp  ),
     3                      sqfil3, work(qsqtr3), error )
c.debug
c     write(6,'("after restoring b - error =",i8)') error
c.debug

              if ( error .ne. 0 ) go to 8900

              work(qtemp7) = 0

          end if
c.debug
c     write(6,'("after restoring a and b")')
c     write(6,'("arindx, arowls, anzp  , anzero, inuse  = ", 5i8)')
c    1            arindx, arowls, anzp  , anzero, inuse 
c     write(6,'("brindx, browls, bnzp  , bnzero         = ", 5i8)')
c    1            brindx, browls, bnzp  , bnzero
c.debug
 
c         ---------------------------------------------
c         ... clean up at end of out-of-core processing
c         ---------------------------------------------

          iofile = wafil1
          call xdslw3 ( wafil1, 3, error )
          if ( error .ne. 0 ) go to 8900  

          iofile = wafil2
          call xdslw3 ( wafil2, 3, error )
          if ( error .ne. 0 ) go to 8900  

          iofile = wafil3
          if ( bmxtyp .eq. 1 ) call xdslw3 ( wafil3, 3, error )
          if ( error .ne. 0 ) go to 8900  
 
      end if

c     ---------------------------------------------------------
c     ... the relative assembly indices have now been built.
c         finish processing
c     ---------------------------------------------------------
 
      t1 = work ( qvaltm )
      w1 = work ( qvalwl )
 
      call xdslt2 ( t1, w1, t3, w3 )
 
c     ----------------------------------
c     ... print optional timing messages
c     ----------------------------------
 
      if ( msglvl .ge. 1 ) then 
          write ( output, 80000 ) t2, w2, t3, w3
      end if
 
c     ---------------------------------------------
c     ... store information into.CMNication area
c     ---------------------------------------------

      work ( qstage ) = 40
c.debug
c     write(6,'("at end of xdslvf - inuse  = ", i8)') inuse
c.debug

      work ( qinuse ) = inuse 

      maxwrk          = work ( qmxuse )
      work ( qmxuse ) = max ( maxwrk, mxused )

      work ( qcndnm ) = tnorma

      work ( qtemp1 ) = 0.
      work ( qtemp2 ) = 0.
      work ( qtemp3 ) = 0.
      work ( qtemp4 ) = 0.
      work ( qtemp5 ) = 0.
      work ( qtemp8 ) = 0.
      work ( qtemp9 ) = 0.

      work ( qvaltm ) = t3
      work ( qvalwl ) = w3

      work ( qlndxg ) = lindxg

      work ( qnzla  ) = anzero

      work ( qadiag ) = adiag
      work ( qarowi ) = arowls
      work ( qofdia ) = anzp

      work ( qxadj  ) = xarndx
      work ( qarndx ) = arindx

      if ( bmxtyp .gt. 1 ) then
          bnzero = 0
          browls = 1
          bnzp   = 1
          xbrndx = 1
          brindx = 1
      end if

      if ( bmxtyp .gt. 3 ) then
          bdiag  = 1
      end if

      work ( qnzlb  ) = bnzero

      work ( qbdiag ) = bdiag
      work ( qbrowi ) = browls
      work ( qbvalu ) = bnzp

      work ( qbcols ) = xbrndx
      work ( qlndxb ) = brindx

c.debug
c     write(6,'("at end of xdslvf")')
c     write(6,'("xarndx, adiag , arindx, arowls, anzp   = ", 5i8)')
c    1            xarndx, adiag , arindx, arowls, anzp  
c     write(6,'("anzero, amncor                         = ", 5i8)')
c    1            anzero, amncor
c     write(6,'("xbrndx, bdiag , brindx, browls, bnzp   = ", 5i8)')
c    1            xbrndx, bdiag , brindx, browls, bnzp  
c     write(6,'("bnzero, bmncor, bmxtyp, lindxg         = ", 5i8)')
c    1            bnzero, bmncor, bmxtyp, lindxg
c.debug
     
      go to 9000
 
c---------------------------------------------------------------------
 
c     --------------
c     ... error trap
c     --------------
 
c     ---------------------------------
c     ... incorrect processing sequence
c     ---------------------------------
 
 8000 continue
      error = -450
      call hherr ( 3, 'xdslvf', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88000 ) error, stage
      go to 8999
 
c     ---------------------------------------
c     ... insufficient storage for this stage
c     ---------------------------------------
 
 8100 continue
      error = -401
      call hherr ( 2, 'xdslvf', error, wkreqd )
      if ( msglvl .gt. 0 ) write ( output, 88100 ) error, wkreqd, lwork
      go to 8999
 
c     -------------------------
c     ... error on matrix input
c     -------------------------
 
 8200 continue
      error = -404
      call hherr ( 1, 'xdslvf', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88200 ) error, badrow, badcol
      work(qstage) = -1
      go to 8999
 
c     ---------------------------------
c     ... internal data structure error
c     ---------------------------------
 
 8300 continue
      error = -403
      call hherr ( 3, 'xdslvf', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88300 ) error
      work(qstage) = -1
      go to 8999

c     ----------------------------
c     ... i/o error on file wafil1
c     ----------------------------
 
 8900 continue
      error = -410
      call hherr ( 3, 'xdslvf', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88900 ) error, iofile
      go to 8999

c     --------------------------------------
c     ... set stage to prevent further calls
c     --------------------------------------
 
 8999 continue
      work(qstage) = -1
 
c---------------------------------------------------------------------
 
c     ------------------------
c     ... end of module xdslvf
c     ------------------------
 
 9000 continue
      return
 
c---------------------------------------------------------------------
 
c     -----------
c     ... formats
c     -----------
 
80000 format ( /5x, 'cpu  time for numeric value input       = ', f15.6
     1         /5x, 'wall time for numeric value input       = ', f15.6 
     2         /5x, 'cpu  time for numeric value processing  = ', f15.6 
     3         /5x, 'wall time for numeric value processing  = ', f15.6)
 
81000 format ( /5x, 'number of rows in matrix                = ', i15
     1         /5x, 'number of nonzeroes in full matrix      = ', i15
     2         /5x, 'number of compressed nodes in matrix    = ', i15
     3         /5x, 'number of compressed nonzeroes          = ', i15 )
 
88000 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslvf executed in an'
     2         /5x, 'incorrect sequence.  current stage = ', i10,
     3          5x, 'should be 31.')
 
88100 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslvf requires ', i15
     2         /5x, 'words of workspace and has only ', i15,
     3              ' available.' )
 
88200 format ( /5x, '*** fatal error no. ', i5, ' *** invalid ',
     1              'matrix input' 
     2         /5x, 'no entry found for row = ', i10,
     3              ' and column = ', i10, '.' )
 
88300 format ( /5x, '*** fatal error no. ', i5, ' *** internal ',
     1              'data structure error in value input.')

88900 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslvf encountered i/o error'
     2         /5x, 'on i/o file no. ', i15 )
 
c---------------------------------------------------------------------
 
      end
