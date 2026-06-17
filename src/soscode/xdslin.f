      subroutine xdslin ( neqns,  mtxtyp, msglvl, output, sqfil1,
     1                    sqfil2, wafil1, wafil2, wafil4, wafil5,
     2                    work,   lwork,  error )
 
c
c     purpose
c     -------
c
c     xdslin is the initiallization subroutine for structural input
c     of the matrix.  it must be called before any other calls.
c
c     created         27-jan-89   -- rgg --
c     last modified   03-dec-96   -- rgg -- allowed for 
c                                           out-of-core assembly
c                     30-jul-98   -- rgg -- added wafil4 and wafil5
c                     10-oct-01   -- dkw -- initialize qppfmn
c
c     input arguments
c     ---------------
c
c     neqns       i   number of equations
c     mtxtyp      c   'u' or 'u' - unsymmetric
c                     otherwise symmteric
c     msglvl      i   message level.
c     output      i   output unit.
c     sqfil1      i   unit number for sequential i/o file no. 1
c     sqfil2      i   unit number for sequential i/o file no. 2
c     wafil1      i   unit number for word add. i/o file no. 1
c     wafil2      i   unit number for word add. i/o file no. 2
c     wafil4      i   unit number for word add. i/o file no. 4
c     wafil5      i   unit number for word add. i/o file no. 5
c     lwork       i   length of work array.
c
c     output arguments
c     ----------------
c
c     error       i   error flag
c                     =    0  normal return
c                     = -101  lwork not large enough.
c                     = -102  neqns .le. 0
c
c---------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      character(len=*)    mtxtyp
 
      integer             neqns,  msglvl, output, lwork,
     1                    sqfil1, sqfil2, wafil1, wafil2, 
     2                    wafil4, wafil5, error
 
      double precision    work(*)
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------

      character(len=55)   libver
 
      integer             i,      ineqn1, wkreqd
 
      logical             lsame
 
      external            lsame
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      integer             xdslni
 
      external            xdslni
 
c---------------------------------------------------------------------
 
      error  = 0

      call xislvn ( libver )
 
      if ( msglvl .ge. 2 ) write ( output, 68000 ) libver


      if ( neqns .le. 0 ) go to 8200

      ineqn1 = xdslni(neqns+1)
 
      wkreqd = 10*ineqn1 + lncomm
      if ( lwork .le. wkreqd ) go to 8100
 
c     ---------------------------------------------
c     ... store information into.CMNication area
c     ---------------------------------------------
 
      do i = 1, lncomm
          work(i) = 0.0d0
      enddo
 
      work ( qstage ) = 1.
      work ( qmsglv ) = msglvl
      work ( qoutpu ) = output
 
      work ( qneqns ) = neqns
 
      if ( lsame ( mtxtyp(1:1), 'U' ) ) then
          work ( qmxtyp ) = 2.
      else
          work ( qmxtyp ) = 1.
      end if

      work ( qbmxty ) = 4.
 
      work ( qsqfl1 ) = sqfil1
      work ( qsqfl2 ) = sqfil2
      work ( qwafl1 ) = wafil1
      work ( qwafl2 ) = wafil2
      work ( qwafl4 ) = wafil4
      work ( qwafl5 ) = wafil5

c     -----------------------------------------------------
c     ... set defaults for factorization and solve controls
c     -----------------------------------------------------

      work ( qmxzer ) = max ( 100., .01 * neqns )
      work ( qmncor ) = 0.
      work ( qpvttl ) = .01
      work ( qcmajr ) = 1.
      work ( qfpnls ) = 128.
      work ( qfupds ) =  32.
      work ( qslvbs ) =   8.
      work ( qcndnc ) = 1.
      work ( qsavea ) = 0.
      work ( qextfl ) = 0.
      work ( qzpcnt ) = 0.
      work ( qnpcnt ) = 0.
      work ( qppfmn ) = 0.
 
      go to 9000
 
c---------------------------------------------------------------------
 
c     --------------
c     ... error trap
c     --------------
 
c     ---------------------------------------
c     ... insufficient storage for this stage
c     ---------------------------------------

 8100 continue
      error = -101
      call hherr ( 2, 'xdslin', error, wkreqd )
      if ( msglvl .gt. 0 ) write ( output, 88100 ) error, wkreqd, lwork
      go to 8900
 
c     ----------------
c     ... neqns .le. 0
c     ----------------
 
 8200 continue
      error = -102
      call hherr ( 1, 'xdslin', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88200 ) error, neqns
      go to 8900
 
c     --------------------------------------
c     ... set stage to prevent further calls
c     --------------------------------------
 
 8900 continue
      work(qstage) = -1
 
c---------------------------------------------------------------------
 
c     ------------------------
c     ... end of module xdslin
c     ------------------------

 9000 continue
      return
 
c---------------------------------------------------------------------
 
c     -----------
c     ... formats
c     -----------
 
68000 format ( 
     1 /1x,'==========================================================='
     2 /1x,'= ', a55,                                               ' ='
     3 /1x,'= multifrontal structural input phase                     ='
     4 /1x,'==========================================================='
     5       )
 
88100 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslin requires ', i15
     2         /5x, 'words of workspace and has only ', i15,
     3              ' available.' )
 
88200 format ( /5x, '*** fatal error no. ', i5, ' *** number of ',
     1              'equations input to'
     2        /5x,  'subroutine xdslin is nonpositive.  neqns = ',
     3              i10 )
 
c---------------------------------------------------------------------
 
      end
