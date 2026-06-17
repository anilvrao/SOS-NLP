      subroutine xdslsp ( char  , jinp  , finp  , work  , error )
 
c
c     purpose
c     -------
c
c     xdslsp is the subroutine to change various control parameters for
c     the sparse factor and solve capabilities in bcslib-ext
c
c     created         13-mar-97   -- rgg --
c     last modified   02-may-00   -- rgg -- added full reorthog. control
c                     05-may-00   -- rgg -- case desensitized by using
c                                           xistrc
c                     01-nov-01   -- dkw -- added monitor panel pivots
c                                           failure control
c                     13-nov-01   -- dkw -- changed error code -1,-2,-3
c                                           to -911,-912,-913 and added
c                                           stage test error -100
c                     13-aug-02   -- jgl -- added seminorm check control
c                     29-jan-04   -- wrf -- added one shift only control
c
c     input arguments
c     ---------------
c
c     char        c   character string denoting control parameter to set
c     jinp        i   integer input parameter associated with char
c     finp        d   floating point input parameter associated with char
c
c     output arguments
c     ----------------
c
c     work        i   work array.
c     error       i   error flag
c                     =    0  normal return
c                     = -100  incorrect processing path
c                     = -911  character string not recognized.
c                     = -912  illegal value of jinp for the specified
c                             control parameter 
c                     = -913  illegal value of finp for the specified
c                             control parameter 
c
c---------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      character(len=*)    char            
 
      integer             jinp  , error
 
      double precision    finp
 
      double precision    work(*)
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------
 
      integer             msglvl, output, stage

c     ---------------------
c     ... external function
c     ---------------------

      logical             xistrc

      external            xistrc
 
c---------------------------------------------------------------------
 
      error  = 0
 
      stage  = work ( qstage )
      msglvl = work ( qmsglv )
      output = work ( qoutpu )

      if ( stage .lt. 1 ) go to 1000
 
      if ( msglvl .ge. 2 ) write ( output, 60000 ) char, jinp, finp

c     ------------------------------------------------------------
c     ... test of which control parameter and take the appropriate 
c         action
c     ------------------------------------------------------------

      if ( xistrc ( char, 'supernode relaxation parameter' ) ) then

c         --------------------------------------
c         ... set supernode relaxation parameter
c         --------------------------------------
 
          work ( qmxzer ) = jinp

      else if ( xistrc ( char, 'minimum core eigen a and b' ) ) then

c         ------------------------------------------------------
c         ... activate minimum core processing for eigenanalysis
c             forcing both a and b to disk
c         ------------------------------------------------------
 
          work ( qmncor ) = 22.

      else if ( xistrc ( char, 'minimum core eigen a' ) ) then

c         ------------------------------------------------------
c         ... activate minimum core processing for eigenanalysis
c             forcing a to disk
c         ------------------------------------------------------
 
          work ( qmncor ) = 2.

      else if ( xistrc ( char, 'minimum core eigen b' ) ) then

c         ------------------------------------------------------
c         ... activate minimum core processing for eigenanalysis
c             forcing both b to disk
c         ------------------------------------------------------
 
          work ( qmncor ) = 20.

      else if ( xistrc ( char, 'minimum core' ) ) then

c         ----------------------------------------
c         ... activate minimum core processing for 
c             factorization.
c         ----------------------------------------
 
          if ( work ( qsavea ) .eq. 0 ) then
              work ( qmncor ) = 1.
          else
              work ( qmncor ) = 2.
          end if

      else if ( xistrc ( char, 'save original matrix' ) ) then

c         -----------------------------------------
c         ... set flags to preserve original matrix
c         -----------------------------------------

          work ( qsavea ) = 1.
          if ( work ( qmncor ) .ne. 0 ) work ( qmncor ) = 2.

      else if ( xistrc ( char, 'pivot tolerance' ) ) then

c         -----------------------
c         ... set pivot tolerance
c         -----------------------

          if ( finp .lt. 0.  .or.  finp .gt. 0.5 ) then
              error = -913
              go to 8000
          end if

          work(qpvttl) = finp

      else if ( xistrc ( char, 'column major' ) ) then

c         ----------------------------------------------------
c         ... set factorization storage scheme to column major
c         ----------------------------------------------------

          work(qcmajr) = 1.

      else if ( xistrc ( char, 'row major' ) ) then

c         -------------------------------------------------
c         ... set factorization storage scheme to row major
c         -------------------------------------------------

          work(qcmajr) = 0.

      else if ( xistrc ( char, 'panel size' ) ) then

c         ----------------------------------------
c         ... set maximum factorization panel size
c         ----------------------------------------

          if ( jinp .le. 0. ) then
              error = -912
              go to 8000
          end if

          work(qfpnls) = jinp

      else if ( xistrc ( char, 'update size' ) ) then

c         -----------------------------------------
c         ... set maximum factorization update size
c         -----------------------------------------

          if ( jinp .le. 0. ) then
              error = -912
              go to 8000
          end if

          work(qfupds) = max ( jinp, 4 )

      else if ( xistrc ( char, 'solve block size' ) ) then

c         ------------------------------------------
c         ... set maximum number of right-hand-sides
c             processed at one time
c         ------------------------------------------

          if ( jinp .le. 0. ) then
              error = -912
              go to 8000
          end if

          work(qslvbs) = max ( jinp, 2 )


      else if ( xistrc ( char, 'estimate condition number' ) ) then

c         -----------------------------------------
c         ... set flag to estimate condition number
c         -----------------------------------------

          work(qcndnc) = 1.

      else if ( xistrc ( char, 'no condition number' ) ) then

c         ---------------------------------------------
c         ... set flag not to estimate condition number
c         ---------------------------------------------

          work(qcndnc) = 0.

      else if ( xistrc ( char, 'message level' ) ) then

c         ---------------------
c         ... set message level
c         ---------------------

          work(qmsglv) = jinp

      else if ( xistrc ( char, 'output unit' ) ) then

c         ---------------------
c         ... reset output unit
c         ---------------------

          work(qoutpu) = jinp

      else if ( xistrc ( char, 'limit fill' ) ) then

c         -------------------------------------------------
c         ... set control to limit fill during the numeric
c             factorization due to pivoting
c             note that = 0. means no limit
c         -------------------------------------------------

          if ( finp .lt. 0. ) then
              error = -913
              go to 8000
          end if

          work(qextfl) = finp

      else if ( xistrc ( char, 'no pivot controls' ) ) then

c         -------------------------------------------------
c         ... set control on pivots to default
c         -------------------------------------------------

          work(qzpcnt) = 0.
          work(qnpcnt) = 0.

      else if ( xistrc ( char, 'adjust near zero pivots' ) ) then

c         -------------------------------------------------
c         ... set control to perturb zero pivots during a
c             no pivoting factorization.
c             statistics can be retrieved with xdslpc
c         -------------------------------------------------

          if ( finp .lt. 0. ) then
              error = -913
              go to 8000
          end if

          work(qzpcnt) = finp

      else if ( xistrc ( char, 'monitor negative pivots' ) ) then

c         ---------------------------------------------------
c         ... set control to monitor negative pivots during a
c             no pivoting factorization.
c             statistics can be retrieved with xdslpc
c         ---------------------------------------------------

          work(qnpcnt) = 1.

      else if ( xistrc ( char, 'abort on negative pivots' ) ) then

c         -----------------------------------------------------
c         ... set control to abort on negative a pivot during a
c             no pivoting factorization.
c             statistics can be retrieved with xdslpc
c         -----------------------------------------------------

          work(qnpcnt) = 2.

      else if ( xistrc ( char, 'monitor panel pivots failure' ) ) then

c         ---------------------------------------------------
c         ... set control to monitor panel pivots failurering
c         ---------------------------------------------------

          work ( qppfmn ) = 1.

      else if ( xistrc ( char, 'full reorthogonalization' ) ) then

c         ----------------------------------------------------
c         ... set full reorthogonalization control for lanczos
c         ----------------------------------------------------

          work ( qfullr ) = 1.

      else if ( xistrc ( char, 'partial reorthogonalization' ) ) then

c         -------------------------------------------------------
c         ... set partial reorthogonalization control for lanczos
c         -------------------------------------------------------

          work ( qfullr ) = 0.

      else if ( xistrc ( char, 'enable seminorm check' ) ) then

c         ----------------------------------------------------
c         ... set semi-norm check control for lanczos
c         ----------------------------------------------------

          work ( qsmchk ) = 1.

      else if ( xistrc ( char, 'disable seminorm check' ) ) then

c         -------------------------------------------------------
c         ... set no semi-norm check control for lanczos
c         -------------------------------------------------------

          work ( qsmchk ) = 0.

      else if ( xistrc ( char, 'enable one shift only' ) ) then

c         ----------------------------------------------------
c         ... set only one shift control for lanczos
c         ----------------------------------------------------

          work ( qonesh ) = 1.

      else if ( xistrc ( char, 'disable one shift only' ) ) then

c         -------------------------------------------------------
c         ... set default shift strategy control for lanczos
c         -------------------------------------------------------

          work ( qonesh ) = 0.

      else

          error = -911
          go to 8000

      end if
 
      go to 9000

c---------------------------------------------------------------------
 
c     --------------
c     ... error trap
c     --------------
 
 1000 continue
      error = -100
      call hherr ( 3, 'xdslsp', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 81000 ) error, stage
      go to 9000

 8000 continue
      call hherr ( 3, 'xdslsp', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88000 ) error, char,
     1                                             jinp , finp
      go to 9000
 
c---------------------------------------------------------------------
 
c     ------------------------
c     ... end of module xdslsp
c     ------------------------
 
 9000 continue
      return
 
c---------------------------------------------------------------------
 
c     -----------
c     ... formats
c     -----------
 
60000 format ( / 5x, 'sparse factorization and solve default ',
     1              'control parameters being reset.'
     2         /10x, 'control parameter = ', a
     3         /10x, 'integer input     = ', i15
     4         /10x, 'real input        = ', 1pd15.5 )

81000 format ( / 5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslsp executed in an'
     2         /5x, 'incorrect sequence.  current stage = ',
     3              i10, 5x, 'should be greater or equal to 1.')

88000 format ( / 5x, '*** input to subroutine xdslsp ignored - error ',
     1              'no. = ', i5
     2         /10x, 'control parameter = ', a
     3         /10x, 'integer input     = ', i15
     4         /10x, 'real input        = ', 1pd15.5 )

c---------------------------------------------------------------------
 
      end
