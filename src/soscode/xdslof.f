      subroutine xdslof( unsym , msglvl, msgunt, neqns,  nzla,   fnzlf,
     1                   nsuper, nsnind, mxnind, snmson, mstack,
     2                   stkstr, mxtri,  mfront, mxtotf, maxzer, sqfil1,
     3                   sqfil2, fpnlsz, fupdsz, pvttol, 
     4                   ncomp , nzcomp, lwksym, lwkifs, lwkofs )
 
 
c
c  purpose -- to determine the storage needed for
c             the succeeding program segments
c
c  created            -- 24-jun-87, cca
c  last modifications -- 24-jun-87, cca
c                        08-may-89, rgg  split intdat into scalar
c                                        quantities
c                        08-feb-91, rgg  mods to allow no i/o to sqfil1
c
c  input variables
c
c      unsym  -- unsymmetric flag
c      msglvl -- message level
c      msgunt -- message unit
c      neqns  -- number of equations
c      nzla   -- number of nonzeroes in original matrix lower triangle
c      fnzlf  -- exact number of nonzeroes in factor's lower triangle
c      nsuper -- number of supernodes in this partition
c      nsnind -- number of supernode indices
c      mxnind -- maximum number of indices for merge routine
c      snmson -- maximum number of merges for a supernode
c      mstack -- maximum depth of the stack
c      stkstr -- maximum stack storage
c      mxtri  -- maximum size of triangle of a supernode's entries
c      mfront -- maximum size of total front matrix
c      mxtotf -- maximum size of trapezoid in the front
c      maxzer -- supernodal relaxation factor
c      sqfil1 -- i/o unit for integer sequential i/o.  if .le. 0 then
c                more things are stored in memory
c      sqfil2 -- i/o unit for integer sequential i/o.  if .le. 0 then
c                more things are stored in memory
c      fpnlsz -- factorization panel size
c      fupdsz -- factorization update panel size
c      pvttol -- pivot tolerance
c      ncomp  -- number of compressed nodes
c      nzcomp -- size of compressed adjacency structure
c
c  output variables --
c
c      lwksym -- storage needed for the symbolic factorization
c      lwkifs -- storage needed for in-core numerical factorization
c                and solve
c      lwkofs -- storage needed for out-of-core numerical factorization
c                and solve
c
c  =====================================================================
 
      integer             msglvl, msgunt
 
      integer             lwkifs, lwkofs, lwksym, mxlfrt, twork,
     1                    maxzer, mfront, snmson, mstack, mxnind,
     2                    mxtotf, mxtri , neqns , nsnind, 
     3                    nsuper, nzla  , sqfil1, sqfil2,
     4                    stkstr, fpnlsz, fupdsz, ncomp , nzcomp,
     5                    locbfr
 
      logical             unsym

      double precision    fnzlf, pvttol
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
      integer             k, l

      logical             qincor
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      integer             xdslni
 
      external            xdslni
c  =====================================================================
 
c     -------------------------------------------------------
c     ... compute maximum value for the order of a front from
c         maximum storage required for a front
c     -------------------------------------------------------
 
      mxlfrt = ( -1. + sqrt ( 1. + 8. * mfront ) ) / 2.

c  ---------------------------------------------------------
c  optionally print out the current values of the parameters
c  ---------------------------------------------------------
 
      if ( msglvl .ge. 2 ) then
          write ( msgunt, 68100 )
     1          neqns , nzla  , fnzlf , nsuper, nsnind, mxnind,
     2          snmson, mstack, stkstr, mxtri , mfront, mxtotf, maxzer,
     3          pvttol, ncomp , nzcomp, fpnlsz, fupdsz, mxlfrt
      endif

c     ------------------------------------------------------------
c     ... determine the number of versions of supernodal indices
c         that have to be stored during the numeric factorization.
c     ------------------------------------------------------------

      if ( sqfil2 .gt. 0 ) then
          if ( pvttol .eq. 0. ) then
              k = 1
          else 
              k = 3
          end if
      else
          if ( pvttol .eq. 0. ) then
              k = 2
          else 
              k = 3
          end if
      end if
c.debug
c     write(6,'("in xdslof - k             = ", 2i10)')
c    1                        k             
c.debug

c     -------------------------------------------------
c     ... test for 32 bit integer limitation of in-core
c         factorization
c     -------------------------------------------------

      qincor = .true.

      if ( unsym ) then
          l = 2 * fnzlf
          if ( l .ne. 2 * fnzlf ) qincor = .false.
      else
          l = fnzlf
          if ( l .ne.     fnzlf ) qincor = .false.
      end if

c  ------------------------------------------------------------
c     ... determine the work space needed for later modules.
c         the variable twork represents the amount of temporary
c         workspace needed during factorization and condition
c         number estimation.  the maximum of the two is chosen.
c  ------------------------------------------------------------
 
      l = max ( 2*nzcomp +  6*neqns + 3*nsuper + ncomp,
     1            nzcomp + 10*neqns + 4*nsuper + ncomp + 2*nsnind )

      lwksym = xdslni ( l + 21 ) + lncomm

      locbfr = max ( neqns, fupdsz*mxlfrt )

      if ( unsym ) then
 
          twork  = max ( 2*stkstr + 2*(fpnlsz+fupdsz)*mxlfrt 
     1                   + 2*xdslni(mxlfrt) + 6*mstack,
     1                   xdslni(neqns+mxlfrt) )
 
          lwkifs = 2*nzla + 2*fnzlf + twork + neqns 
     1           + locbfr + lncomm
     1           + xdslni( 2*nzla + 6*neqns + 6*nsuper + k*nsnind + 20 )
c.debug
c     write(6,'("in xdslof - twork, lwkifs = ", 2i10)')
c    1                        twork, lwkifs 
c.debug
 
          if ( sqfil1 .lt. 0 ) lwkifs = lwkifs + xdslni ( nzla )
 
          twork  = max ( 2*mfront + 2*(fpnlsz+fupdsz)*mxlfrt 
     1                   + 2*xdslni(mxlfrt) + 6*mstack,
     1                   xdslni(neqns+mxlfrt) + neqns + mfront )
 
          lwkofs = 2*nzla + twork + neqns + locbfr + lncomm
     1           + xdslni( 2*nzla + 6*neqns + 6*nsuper + k*nsnind + 20 )
c.debug
c     write(6,'("in xdslof - twork, lwkofs = ", 2i10)')
c    1                        twork, lwkofs 
c.debug
 
          if ( sqfil1 .lt. 0 ) lwkofs = lwkofs + xdslni ( nzla )
 
      else
 
          twork  = max ( stkstr + (fpnlsz+fupdsz)*mxlfrt 
     1                   + xdslni(mxlfrt) + 6*mstack, 
     1                   xdslni(neqns+mxlfrt) )

          lwkifs = nzla + fnzlf + twork + neqns + locbfr + lncomm
     1           + xdslni ( nzla + 6*neqns + 6*nsuper + k*nsnind + 20 )
 
          if ( sqfil1 .lt. 0 ) lwkifs = lwkifs + xdslni ( nzla )
c.debug
c     write(6,'("in xdslof - twork, lwkifs = ", 2i10)')
c    1                        twork, lwkifs 
c.debug

          mfront = (fpnlsz+fupdsz)*mxlfrt
 
          twork  = max ( mfront + xdslni(mxlfrt) + 6*mstack,
     1                   xdslni(neqns+mxlfrt) + neqns + mfront )
 
          lwkofs = nzla + twork + neqns + locbfr + lncomm
     1           + xdslni ( nzla + 6*neqns + 6*nsuper + k*nsnind + 20 )
 
          if ( sqfil1 .lt. 0 ) lwkofs = lwkofs + xdslni ( nzla )
c.debug
c     write(6,'("in xdslof - twork, lwkofs = ", 2i10)')
c    1                        twork, lwkofs 
c.debug
 
      end if

c     --------------------------------------------------
c     ... test if 32 bit integer limitation has happened
c     --------------------------------------------------

      if ( .not. qincor .or. lwkifs .lt. lwkofs ) lwkifs = 0
 
c  ---------------------------------------------------------
c  optionally print out the work space for following modules
c  ---------------------------------------------------------
 
      if ( msglvl .ge. 2 ) then
          write ( msgunt, 68200 ) lwksym, lwkifs, lwkofs
      end if
 
c  =====================================================================
 
      return
 
c  =====================================================================
 
c     -----------
c     ... formats
c     -----------
 
68100 format ( /5x, '  neqns   nzla          nzlf nsuper nsnind mxnind'
     1         /5x, 2i7, f14.0, 3i7
     2        //5x, ' snmson mstack stkstr  mxtri mfront mxtotf maxzer'
     3         /5x, 7i7  
     4        //5x, ' pvttol  ncomp nzcomp fpnlsz fupdsz mxlfrt'
     5         /5x, f7.3, 6i7 )
 
68200 format ( /5x, 'work space for symbolic factorization   = ', i15
     1         /5x, 'work space for in-core     num.fact.    = ', i15 
     2         /5x, 'work space for out-of-core num.fact.    = ', i15 )
 
      end
