      subroutine xdslsl ( nrhs,   rhs,    ldrhs,  work,   lwork,
     1                    wkreqd, error )
 
c
c     purpose
c     -------
c
c     xdslsl is the top level driver for numeric solution phase.
c
c     created         30-jan-89   -- rgg --
c     last modified   09-may-89   -- rgg -- added indexing parameters
c                     02-nov-90   -- mlc -- changed usage of nrhs in
c                                           dimensioning to *
c                     09-nov-95   -- rgg -- added temp1 and temp2
c                                           to reflect mods to xdsls5
c                                           and xdsls6
c                     27-feb-97   -- rgg -- converted to panels
c                                           switched to dgemv and dgemm
c                     10-mar-97   -- rgg -- reconstructed i/o for
c                                           out-of-core solves
c                     27-mar-98   -- rgg -- converted to processing 
c                                           slvbsz rhs at one time
c
c     input arguments
c     ---------------
c
c     nrhs        i   the number of right-hand-sides.
c     ldrhs       i   leading dimension of the right-hand-side array.
c     lwork       i   length of work array.
c
c     input/output arguments
c     ----------------------
c
c     rhs         d   right hand side array.
c
c     work        d   work array.  on input it contains the
c                    .CMNication area and all active arrays.
c
c     output arguments
c     ----------------
c
c     error       i   error flag
c                     =    0  normal return
c                     = -600  incorrect processing path.
c                     = -601  lwork not large enough.
c                     = -602  nrhs .le. 0
c                     = -603  ldrhs .lt. neqns
c                     = -604  i/o error on wafil1
c                     = -605  i/o error on wafil4
c
c---------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      integer             nrhs,   ldrhs,  lwork, error
 
      double precision    rhs(ldrhs,*),        work(*)
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------
 
      integer             diag,   i,      invp,   invsup, lindxg,
     1                    lnz,    ltemp3, mxlfrt, mxtotf, msglvl
 
      integer             neqns,  nsuper, 
     1                    output, perm,   pvtblk, stage,  sup,
     2                    slvbsz, temp1,  temp3 , wafil1, wafil4,
     3                    wkreqd, xlindx, xpanel, xsup
 
      double precision    t1,     t2,     w1,     w2,     fnzlf
 
      logical             unsym
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      external            dsctr,  xdsls3, xdsls4,
     1                    xdslp4, xdslt1, xdslt2
 
c
c---------------------------------------------------------------------
 
      call xdslt1 ( t1, w1 )
      error  = 0
 
      msglvl = work ( qmsglv )
      output = work ( qoutpu )
      unsym  = work ( qmxtyp ) .eq. 2.
 
      if ( msglvl .ge. 2 ) write ( output, 68000 )
 
      if ( nrhs .le. 0 ) go to 8200
 
c     ----------------------------------------------------
c     ... extract information from the.CMNication area.
c     ----------------------------------------------------
 
      stage  = work ( qstage )
      if ( stage .ne. 50 ) go to 8000
 
      neqns  = work ( qneqns )
      fnzlf  = work ( qnzlf2 )
      nsuper = work ( qnsupe )
      mxtotf = work ( qmxtot )
      mxlfrt = work ( qmxnin ) 
 
      if ( ldrhs .lt. neqns ) go to 8300
 
      perm   = work ( qperm  )
      invp   = work ( qinvp  )
      xsup   = work ( qxsup2 )
      xpanel = work ( qxpanl )
      xlindx = work ( qxlnd2 )
      lindxg = work ( qlndg2 )
      diag   = work ( qdiag  )
      lnz    = work ( qlnz   )
      pvtblk = work ( qpivot )
      sup    = work ( qsup   )
      invsup = work ( qinvsu )
 
c     -----------------------------------------------
c     ... compute other pointers into the work array.
c     -----------------------------------------------

      slvbsz = work ( qslvbs )
      temp1  = work ( qtemp2 )
 
      if ( work (qtemp1) .eq. 0. ) then
          wkreqd  = temp1  + max ( neqns, slvbsz*mxlfrt ) - 1
      else
          wkreqd  = temp1  + max ( neqns, mxtotf+slvbsz*mxlfrt ) - 1
      end if
c.debug
c     write(6,'("in xdslsl")')
c     write(6,'("perm  , invp  , xsup  , xpanel, xlindx = ",5i8)')
c    1            perm  , invp  , xsup  , xpanel, xlindx 
c     write(6,'("lindxg, diag  , lnz   , pvtblk, sup    = ",5i8)')
c    1            lindxg, diag  , lnz   , pvtblk, sup   
c     write(6,'("invsup, slvbsz, temp1 , work(qtemp1)   = ",5i8)')
c    1            invsup, slvbsz, temp1 , work(qtemp1) 
c     write(6,'("mxlfrt, mxtotf, wkreqd, lwork          = ",5i8)')
c    1            mxlfrt, mxtotf, wkreqd, lwork       
c.debug
 
      if ( wkreqd .gt. lwork ) go to 8100
 
c     --------------------------------------------
c     ... permute the rhs, solve, and permute back
c     --------------------------------------------
 
      do i = 1, nrhs
 
          if ( msglvl .ge. 3 ) then
              call xdslp4 ( 'org. rhs', neqns, rhs(1,i), output )
          end if
 
          call dsctr ( neqns, rhs(1,i), work(invp), work(temp1) )
          call dsctr ( neqns, work(temp1), work(invsup), rhs(1,i) )
 
          if ( msglvl .ge. 4 ) then
              call xdslp4 ( 'perm. rhs', neqns, rhs(1,i), output )
          end if
 
      enddo
 
      if ( work ( qtemp1 ) .eq. 0.d0 ) then
 
c          -------------------------
c          ... perform in-core solve
c          -------------------------
c.debug
c     write(6,'("before call to xdsls3")')
c     write(6,'("slvbsz = ", i8)') slvbsz
c.debug

           call xdsls3 ( unsym, neqns, nsuper, nrhs, ldrhs,
     1                   work(xsup),   work(xpanel), work(xlindx),
     2                   work(lindxg), work(pvtblk), work(diag),
     3                   work(lnz),    rhs,          
     4                   slvbsz,       work(temp1) )
 
      else
 
c         -----------------------------
c         ... perform out-of-core solve
c         -----------------------------
 
          temp3  = temp1 + slvbsz*mxlfrt
          ltemp3 = lwork - temp3 + 1

          wafil1 = work ( qwafl1 )
          wafil4 = work ( qwafl4 )
c.debug
c     write(6,'("before call to xdsls4 - wafil1, mxlfrt, ltemp3 = ", 
c    1        3I8)')                      wafil1, mxlfrt, ltemp3
c     write(6,'("                        slvbsz, wafil4         = ", 
c    1        2i8)')                      slvbsz, wafil4
c.debug
 
          call xdsls4 ( unsym , wafil1, wafil4, neqns, nrhs, ldrhs,  
     1                  rhs   , nsuper,   work(xsup),   work(xpanel),
     2                  work(xlindx), work(lindxg), work(diag),
     3                  work(pvtblk), slvbsz,       work(temp1),  
     4                  ltemp3,       work(temp3),  error    )
 
          if ( error .eq. -1 ) go to 8400
          if ( error .eq. -2 ) go to 8500
 
          if ( unsym ) then
              work ( qwatr1 ) = work ( qwatr1 ) + fnzlf
              work ( qwatr4 ) = work ( qwatr4 ) + fnzlf
          else
              work ( qwatr1 ) = work ( qwatr1 ) + 2 * fnzlf
          end if
 
      end if
 
c     -----------------------------------------------
c     ... permute the solution back to original order
c     -----------------------------------------------
 
      do i = 1, nrhs
 
          if ( msglvl .ge. 4 ) then
              call xdslp4 ( 'perm. sol', neqns, rhs(1,i), output )
          end if
 
          call dsctr ( neqns, rhs(1,i), work(sup), work(temp1) )
          call dsctr ( neqns, work(temp1), work(perm), rhs(1,i) )
 
          if ( msglvl .ge. 3 ) then
              call xdslp4 ( 'comp. sol', neqns, rhs(1,i), output )
          end if
 
      enddo
 
      call xdslt2 ( t1, w1, t2, w2 )
 
      if ( msglvl .ge. 1 ) write ( output, 80000 ) t2, w2
 
c     ---------------------------------------------
c     ... store information into.CMNication area
c     ---------------------------------------------
 
      work ( qslvtm ) = t2 / nrhs
      work ( qslvwl ) = w2 / nrhs
 
      go to 9000
 
c---------------------------------------------------------------------
 
c     --------------
c     ... error trap
c     --------------
 
c     ---------------------------------
c     ... incorrect processing sequence
c     ---------------------------------
 
 8000 continue
      error = -600
      call hherr ( 3, 'xdslsl', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88000 ) error, stage
      go to 9000
 
c     ---------------------------------------
c     ... insufficient storage for this stage
c     ---------------------------------------
 
 8100 continue
      error = -601
      call hherr ( 2, 'xdslsl', error, wkreqd )
      if ( msglvl .gt. 0 ) write ( output, 88100 ) error, wkreqd, lwork
      go to 9000
 
c     ----------------
c     ... nrhs .le. 0
c     ----------------
 
 8200 continue
      error = -602
      call hherr ( 1, 'xdslsl', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88200 ) error, nrhs
      go to 9000
 
c     --------------------
c     ... ldrhs .lt. neqns
c     --------------------
 
 8300 continue
      error = -603
      call hherr ( 1, 'xdslsl', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88300 ) error, ldrhs, neqns
      go to 9000
 
c     -----------------------
c     ... i/o error on wafil1
c     -----------------------
 
 8400 continue
      error = -604
      call hherr ( 3, 'xdslsl', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88400 ) error, wafil1
      go to 9000
 
c     -----------------------
c     ... i/o error on wafil4
c     -----------------------
 
 8500 continue
      error = -605
      call hherr ( 3, 'xdslsl', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88400 ) error, wafil4
      go to 9000
 
c---------------------------------------------------------------------
 
c     ------------------------
c     ... end of module xdslsl
c     ------------------------
 
 9000 continue
      return
 
c---------------------------------------------------------------------
 
c     -----------
c     ... formats
c     -----------
 
68000 format ( /1x, '========================================='
     1         /1x, '= multifrontal numerical solution phase ='
     2         /1x, '=========================================' )
 
80000 format ( /5x, 'cpu  time for numeric solution          = ', f15.6
     1         /5x, 'wall time for numeric solution          = ', f15.6)
 
88000 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslsl executed in an'
     2         /5x, 'incorrect sequence.  current stage = ', i10,
     3          5x, 'should be 50.' )
 
88100 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslsl requires ', i15
     2         /5x, 'words of workspace and has only ', i15,
     3              ' available.' )
 
88200 format ( /5x, '*** fatal error no. ', i5, ' *** number of ',
     1              'right-hand-sides input to'
     2         /5x, 'subroutine xdslsl is nonpositive.  nrhs = ',
     3              i10 )
 
88300 format ( /5x, '*** fatal error no. ', i5, ' *** leading ',
     1              'dimension of rhs array'
     2         /5x, 'input to subroutine xdslsl is less than ',
     3              'the number of equations.  '
     4        /10x, 'ldrhs = ', i10, '  neqns = ', i10 )
 
88400 format ( /5x, '*** fatal error no. ', i5, ' *** i/o error',
     1              'encountered in subroutine'
     2         /5x, 'xdslsl on word addressable i/o file no. ',
     3              i10 )
 
c---------------------------------------------------------------------
 
      end
