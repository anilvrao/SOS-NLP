      subroutine xdslps ( work )
 
c
c     purpose
c     -------
c
c     xdslps prints statistics
c
c     created         26-jan-89   -- rgg --
c     last modified   16-dec-96   -- rgg -- for compressed nodes
c                     19-oct-99   -- rgg -- added version number 
c                     01-jan-01   -- dkw -- print i/o statistics in
c                                           megabytes format
c
c     input arguments
c     ---------------
c
c     work        d   work array.  on input it contains the
c                    .CMNication area and all active arrays.
c
c---------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      double precision    work(*)
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------

      character(len=55)   libver
 
      integer             inuse,  mxstck, mxtotf, mxused, ncomp ,
     1                    nzcomp, needs,  neqns,  nsnind, nsnin2, 
     2                    nsuper, nzla,   output,
     3                    stkstr, stage,  mfront, mxlfrt, fpnlsz,
     4                    fupdsz, slvbsz
 
      logical             qfirst
 
      double precision    cndnum, rate  ,
     1                    fnzlf , fnzlf2, fctops, slvops,
     2                    ioamnt, ioleng, iotran, iounit
 
c---------------------------------------------------------------------
 
c     ----------------------------
c     ... write scalar information
c     ----------------------------
 
      stage  = work (qstage)
      output = work (qoutpu)
      if ( stage .lt. 10 ) return
 
      inuse  = work (qinuse)
      needs  = work (qneeds)
      mxused = work (qmxuse)
 
      neqns  = work (qneqns)
      nzla   = work (qnzla )
      ncomp  = work ( qncomp ) 
      nzcomp = work ( qnzcmp )

      call xislvn ( libver )
 
      write ( output, 10000 ) libver

      write ( output, 10001 ) neqns, nzla, ncomp, nzcomp, inuse
 
      if ( stage .gt. 10 ) write ( output, 11000 ) mxused
 
      if ( stage .lt. 40 ) write ( output, 12000 ) needs
 
      if ( stage .eq. 10 ) go to 100
 
      fnzlf  = work (qnzlf )
      nsuper = work (qnsupe)
      nsnind = work (qnsnin)
      stkstr = work (qstkst)
      fctops = work (qfctop)
      slvops = work (qslvop)
      cndnum = work (qcndnm)
 
      write ( output, 13000 ) fnzlf, nsuper, nsnind, stkstr
 
      if ( stage .gt. 20 ) then
 
          mxtotf = work (qmxtot)
          mfront = work (qmfron)
          mxlfrt = ( -1. + sqrt ( 1. + 8. * mfront ) ) / 2.
 
          write ( output, 14000 ) mxlfrt, mfront, mxtotf
 
      end if
 
      write ( output, 15000 ) fctops, slvops

      if ( stage .lt. 50 ) go to 100
       
      if ( work(qfcttm) .ne. 0.0 ) then
          fnzlf2 = work (qnzlf2)
          nsnin2 = work (qnsni2)
          mxstck = work (qmxstk)
          write ( output, 16000 ) fnzlf2, nsnin2, mxstck
      end if

      if ( work(qzpcnt) .ne. 0.0 ) then
          write(output,17000) work(qzpcnt+1)
          if ( work(qzpcnt+1) .gt. 0. ) then
              write(output,17010) work(qzpcnt+2), work(qzpcnt+3),
     1                            work(qzpcnt+4)
          end if
      end if

      if ( work(qnpcnt) .ne. 0.0 ) then
          write(output,18000) work(qnpcnt+1)
          if ( work(qnpcnt+1) .gt. 0. ) then
              write(output,18010) work(qnpcnt+2), work(qnpcnt+3),
     1                            work(qnpcnt+4)
          end if
      end if

c     ---------------------------------------------------
c     ... write out parameters controlling performance of
c         factorization and solve
c     ---------------------------------------------------

      fpnlsz = work ( qfpnls )
      fupdsz = work ( qfupds )
      slvbsz = work ( qslvbs )
 
      write ( output, 19000 ) fpnlsz, fupdsz, slvbsz
 
c     ----------------------------
c     ... write timing information
c     ----------------------------
 
  100 continue
      write ( output, 20000 ) work(qinptm)
 
      if ( stage .ge. 20 ) write ( output, 21000 ) work(qordtm)
 
      if ( stage .ge. 30 ) write ( output, 22000 ) work(qsfctm)
 
      if ( stage .ge. 40 ) write ( output, 24000 ) work(qvaltm)
 
      if ( work(qfcttm) .ne. 0.0 ) then
          rate = work(qfctop) / work(qfcttm) / 1000000.0
          write ( output, 25000 ) work(qfcttm), rate, cndnum
      end if
 
      if ( work(qslvtm) .ne. 0.0 ) then
          rate   = work(qslvop) / work(qslvtm) / 1000000.0
          write ( output, 26000 ) work(qslvtm), rate
      end if
 
c     ------------------------
c     ... write i/o statistics
c     ------------------------
 
      qfirst = .true.
 
      iotran = work ( qsqtr1 )
      if ( iotran .gt. 0 ) then
          if ( qfirst ) then
              write ( output, 40000 )
              qfirst = .false.
          end if
          iounit = work ( qsqfl1 )
          ioleng = 8. * work ( qsqln1 ) / 1000000.
          ioamnt = 8. * work ( qsqtr1 ) / 1000000.
          write ( output, 41000 ) iounit, ioleng, ioamnt 
      end if
 
      iotran = work ( qsqtr2 )
      if ( iotran .gt. 0 ) then
          if ( qfirst ) then
              write ( output, 40000 )
              qfirst = .false.
          end if
          iounit = work ( qsqfl2 )
          ioleng = 8. * work ( qsqln2 ) / 1000000.
          ioamnt = 8. * work ( qsqtr2 ) / 1000000.
          write ( output, 41000 ) iounit, ioleng, ioamnt 
      end if
 
      iotran = work ( qsqtr3 )
      if ( iotran .gt. 0 ) then
          if ( qfirst ) then
              write ( output, 40000 )
              qfirst = .false.
          end if
          iounit = work ( qsqfl3 )
          ioleng = 8. * work ( qsqln3 ) / 1000000.
          ioamnt = 8. * work ( qsqtr3 ) / 1000000.
          write ( output, 41000 ) iounit, ioleng, ioamnt 
      end if
 
      iotran = work ( qsqtr4 )
      if ( iotran .gt. 0 ) then
          if ( qfirst ) then
              write ( output, 40000 )
              qfirst = .false.
          end if
          iounit = work ( qsqfl4 )
          ioleng = 8. * work ( qsqln4 ) / 1000000.
          ioamnt = 8. * work ( qsqtr4 ) / 1000000.
          write ( output, 41000 ) iounit, ioleng, ioamnt 
      end if
 
      iotran = work ( qsqtr5 )
      if ( iotran .gt. 0 ) then
          if ( qfirst ) then
              write ( output, 40000 )
              qfirst = .false.
          end if
          iounit = work ( qsqfl5 )
          ioleng = 8. * work ( qsqln5 ) / 1000000.
          ioamnt = 8. * work ( qsqtr5 ) / 1000000.
          write ( output, 41000 ) iounit, ioleng, ioamnt 
      end if
 
      iotran = work ( qwatr1 )
      if ( iotran .gt. 0 ) then
          if ( qfirst ) then
              write ( output, 40000 )
              qfirst = .false.
          end if
          iounit = work ( qwafl1 )
          ioleng = 8. * work ( qwaln1 ) / 1000000.
          ioamnt = 8. * work ( qwatr1 ) / 1000000.
          write ( output, 41000 ) iounit, ioleng, ioamnt 
      end if
 
      iotran = work ( qwatr2 )
      if ( iotran .gt. 0 ) then
          if ( qfirst ) then
              write ( output, 40000 )
              qfirst = .false.
          end if
          iounit = work ( qwafl2 )
          ioleng = 8. * work ( qwaln2 ) / 1000000.
          ioamnt = 8. * work ( qwatr2 ) / 1000000.
          write ( output, 41000 ) iounit, ioleng, ioamnt 
      end if
 
      iotran = work ( qwatr3 )
      if ( iotran .gt. 0 ) then
          if ( qfirst ) then
              write ( output, 40000 )
              qfirst = .false.
          end if
          iounit = work ( qwafl3 )
          ioleng = 8. * work ( qwaln3 ) / 1000000.
          ioamnt = 8. * work ( qwatr3 ) / 1000000.
          write ( output, 41000 ) iounit, ioleng, ioamnt 
      end if
 
      iotran = work ( qwatr4 )
      if ( iotran .gt. 0 ) then
          if ( qfirst ) then
              write ( output, 40000 )
              qfirst = .false.
          end if
          iounit = work ( qwafl4 )
          ioleng = 8. * work ( qwaln4 ) / 1000000.
          ioamnt = 8. * work ( qwatr4 ) / 1000000.
          write ( output, 41000 ) iounit, ioleng, ioamnt 
      end if
 
      iotran = work ( qwatr5 )
      if ( iotran .gt. 0 ) then
          if ( qfirst ) then
              write ( output, 40000 )
              qfirst = .false.
          end if
          iounit = work ( qwafl5 )
          ioleng = 8. * work ( qwaln5 ) / 1000000.
          ioamnt = 8. * work ( qwatr5 ) / 1000000.
          write ( output, 41000 ) iounit, ioleng, ioamnt 
      end if
 
      if ( qfirst ) write ( output, 42000 )
 
c---------------------------------------------------------------------
 
c     ------------------------
c     ... end of module xdslps
c     ------------------------
 
      return
 
c---------------------------------------------------------------------
 
c     -----------
c     ... formats
c     -----------
 
10000 format ( 
     1 /1x,'==========================================================='
     2 /1x,'= ', a55,                                               ' ='
     3 /1x,'= multifrontal statistics                                 ='
     4 /1x,'==========================================================='
     5       )
 
10001 format ( /5x, 'number of equations                     = ', i15
     4         /5x, 'no. of nonzeroes in lower triangle of a = ', i15
     5         /5x, 'number of compressed nodes              = ', i15
     6         /5x, 'no. of compressed nonzeroes in l. tri.  = ', i15
     7         /5x, 'amount of workspace currently in use    = ', i15 )
 
11000 format (  5x, 'max. amt. of workspace used             = ', i15 )
 
12000 format (  5x, 'amt. of workspace needed for next stage = ', i15 )
 
13000 format (  5x, 'no. of nonzeroes in the factor l        = ', f16.0
     1         /5x, 'number of super nodes                   = ', i15
     2         /5x, 'number of compressed subscripts         = ', i15
     3         /5x, 'size of stack storage                   = ', i15 )
 
14000 format (  5x, 'maximum order of a front matrix         = ', i15 
     1         /5x, 'maximum size of a front matrix          = ', i15 
     2         /5x, 'maximum size of a front trapezoid       = ', i15 )
 
15000 format (  5x, 'no. of floating point ops for factor    = ',
     1                                                        1pd15.4
     2         /5x, 'no. of floating point ops for solve     = ',
     3                                                        1pd15.4 )
 
16000 format (  5x, 'actual no. of nonzeroes in the factor l = ', f16.0
     1         /5x, 'actual number of compressed subscripts  = ', i15
     1         /5x, 'actual size of stack storage used       = ', i15 )

17000 format (  5x, 'near zero pivot monitoring activated'
     1         /5x, 'number of pivots adjusted               = ',f15.0)

17010 format (  5x, 'first pivot adjusted                    = ', f15.0
     1         /5x, 'minimum adjustment                      = ', 
     2                                                        1pd15.4
     3         /5x, 'maximum adjustment                      = ', 
     4                                                        1pd15.4 )

18000 format (  5x, 'negative pivot monitoring activated'
     1         /5x, 'number of negative pivots encountered   = ',f15.0)

18010 format (  5x, 'first negative pivot                    = ', f15.0
     1         /5x, 'minimum negative pivot                  = ', 
     2                                                        1pd15.4
     3         /5x, 'maximum negative pivot                  = ', 
     4                                                        1pd15.4 )
 
19000 format (  5x, 'factorization panel size                = ', i15 
     1         /5x, 'factorization update panel size         = ', i15 
     2         /5x, 'solution block size                     = ', i15 )
 
20000 format (  5x, 'time (in seconds) for structure input   = ',
     1                                                          f15.6 )
 
21000 format (  5x, 'time (in seconds) for ordering          = ',
     1                                                          f15.6 )
 
22000 format (  5x, 'time (in seconds) for symbolic factor   = ',
     1                                                          f15.6 )
 
24000 format (  5x, 'time (in seconds) for value input       = ',
     1                                                          f15.6 )
 
25000 format (  5x, 'time (in seconds) for numeric factor    = ',
     1                                                          f15.6,
     2         /5x, 'computational rate (mflops) for factor  = ',
     3                                                          f15.6,
     4         /5x, 'condition number estimate               = ',
     5                                                        1pd15.4 )
 
26000 format (  5x, 'time (in seconds) for numeric solve     = ',
     1                                                          f15.6,
     2         /5x, 'computational rate (mflops) for solve   = ',
     3                                                          f15.6 )
 
40000 format ( /5x, 'i/o statistics:', 3x, 'unit number',
     1          8x, 'length', 8x, 'amount'
     2        /41x,                    '(Mbytes)', 6x, '(Mbytes)'   
     3        /23x, '-----------', 7x, '--------', 6x, '--------'    / )
 
41000 format ( 20x, f14.0, 1x, 2f14.2 )
 
42000 format ( /5x, 'no input or output performed' )
 
c---------------------------------------------------------------------
 
      end
