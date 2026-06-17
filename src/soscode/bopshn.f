
      SUBROUTINE BOPSHN(MODEB)
C
C ======================================================================
C     BOPSHN===>bopshn   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
C
C         PURPOSE:  DISPLAY THE Barrier OPTIONS
C
C-------------------------------------------------------------
C
      PARAMETER (NRSYM=11,NISYM=35,NCSYM=3)
      PARAMETER (NBRSYM=5,NBISYM=3,NBCSYM=2)
C
C======================================================================
C Array equivalences into include file NLPSPR.CMN
C
      INCLUDE '../commons/NLPSPR.CMN'
C
C     COMMON /NPSPRR/ ALFLWR, ...      double precision
C     COMMON /NPSPRI/ INNPER, ...      integer
C     COMMON /NPSPRC/ ALGOPT, ...      character(len=6)
C     COMMON /NPSALG/ ALGNAM           character(len=6) (no equivalence)
C
      DOUBLE PRECISION RNPSPR(NRSYM)
      INTEGER          INPSPR(NISYM+1)
      CHARACTER(LEN=6)  CNPSPR(NCSYM)
C
      EQUIVALENCE (RNPSPR(1),ALFLWR)
      EQUIVALENCE (INPSPR(1),INNPER)
      EQUIVALENCE (CNPSPR(1),ALGOPT)
C
      INCLUDE '../commons/BARNLP.CMN'
C
C     COMMON /BNPSPR/ BIGCON, ...      double precision
C     COMMON /BNPSPI/ IMAXMU, ...      integer
C
      DOUBLE PRECISION RBNPSP(NBRSYM)
      INTEGER          IBNPSP(NBISYM)
C
      EQUIVALENCE (RBNPSP(1),BIGCON)
      EQUIVALENCE (IBNPSP(1),IMAXMU)
C
C End of array equivalences
C======================================================================
C
C-------------------------------------------------------------
C
C         DEFAULT VALUE COMMON
C
      COMMON /NPSDF1/ RDFLT(11),IDFLT(35)
      COMMON /NBSDF1/ RBDFLT(5),IBDFLT(3)
      COMMON /NPSDF2/ CDFLT(3),DFAWLT
      CHARACTER(LEN=6) CDFLT
      CHARACTER(LEN=17) DFAWLT
C
C-------------------------------------------------------------
c       purpose:   display barrier NLP options
c
c       argument:
c
c         modeb    display mode option
c                  = 0   Abbreviated (Terse) format
c                  = 1   Options format
c                  = 2   Full options format
C
C-------------------------------------------------------------
C
c       output character variable format
c
c         outvar = varnam|vty|T|...d...e
c
c       where:
c
c         varnam = 6 character variable name (e.g PGDTOL)
c         vty    = 3 character symbol of the form zxxx 
c                  where first character
c                  z = s ---> sqp variable
c                  z = b ---> barrier variable
c                  and the second two characters
c                  xx ---> number of variable in default list
c         T      = 1 character symbol defining variable type
c                  T = r   real (double precision)
c                  T = i   integer
c                  T = c   character
c         d      = character string containing variable description
c                the string is subdivided with each subdivision 
c                delimited by one of two embedded character strings:
c                %> ---> denotes a line break
c                %< ---> denotes the end of the string
c                
c       nbrvar = number of variables in barrier option display
c       nmlpvr = maximum number of lines per variable in description
c       lnotvr = length of outvar 
C
C-------------------------------------------------------------
C
      parameter (nmlpvr=10,nbrvar=41)
      parameter (lnotvr=10 + nmlpvr*50)
      character(LEN=lnotvr) :: outvar(nbrvar)
c
      dimension ibrief(nbrvar),ishow(nbrvar)
C
      DATA (outvar(I),I=1,10) /
     $'ALFLWRs01rLower Bound on ALFA (Line Search Diagnostic Plot)%<',
     $'ALFUPRs02rUpper Bound on ALFA (Line Search Diagnostic Plot)%<',
     $'BIGCONb01rUpper Bound on Equality Constraint Error%<',
     $'CONTOLs03rConstraint Tolerance%<',
     $'EPSRLFs04rRelative Perturbation Size Parameter%<',
     $'FEATOLb02rInitial Variable Offset for Feasibility%<',
     $'OBJTOLs05rObjective Function Tolerance%<',
     $'PGDTOLs06rProjected Gradient Tolerance%>Convergence Requires:%>(A
     $) Max Absolute Error in Active Constraints and%>    Bounds .LT. CO
     $NTOL%>(B) Max Absolute Error in KT Conditions .LT.%>    PGDTOL*MAX
     $(1,|DELF|)%>(C) |F(X) - FMIN| .LT. OBJTOL%>(D) |Steplength| = ALFA
     $*|P| .LT.%>    SQRT( OBJTOL/(1 + |F(X)|) )*(1 + |X|)%>(E) Correct 
     $Sign For All Lagrange Multipliers%<',
     $'PMULWRb03rLower Bound on Initial Barrier Parameter%<',
     $'PTHTOLb04rCentral Path Convergence Tolerance%<'/
c
      DATA (outvar(I),I=11,20) /
     $'RHOLWRb05rLower Bound on Initial Relaxation Parameter%<',
     $'TOLFILs09rMultifrontal Fill Tolerance%<',
     $'TOLKTCs10rKT Condition Number Tolerance%<',
     $'TOLPVTs11rMultifrontal Pivot Tolerance%<',
     $'IHESHNs01iHessian Matrix Evaluation Option%> 0  = IHESHN      Fin
     $ite Difference or Analytic%>|1| = IHESHN      SR1 (Symmetric Rank 
     $One)%>|2| = IHESHN      BFGS (Symmetric Positive Def.)%>|3| = IHES
     $HN      SSQN (Self-Scaling Quasi-Newton)%>IHESHN < 0 Finite Differ
     $ence Initialization%<',
     $'IMAXMUb01iMaximum Iterations with Fixed Barrier Parameter%<',
     $'IOFLAGs02iOutput Level%>0  = IOFLAG       No Output%>0  < IOFLAG 
     $< 10  Terse Output%>9  < IOFLAG < 20  Standard Output%>19 < IOFLAG
     $ < 30  Interpretive Output%>     IOFLAG = 30  Diagnostic Output%<
     $',
     $'IOFMFRs04iMultifrontal Output Level (0,1,2,3,4)%<',
     $'IOFPATs05iOutput Sparsity Pattern%>.GE. 10 Overrides Default%<',
     $'IPUMF1s10iMultifrontal I/O Unit%<' /
c
      DATA (outvar(I),I=21,30) /
     $'IPUMF2s11iMultifrontal I/O Unit%<',
     $'IPUMF3s12iMultifrontal I/O Unit%<',
     $'IPUMF4s13iMultifrontal I/O Unit%<',
     $'IPUMF5s14iMultifrontal I/O Unit%<',
     $'IPUMF6s15iMultifrontal I/O Unit%<',
     $'IPUMF7s16iMultifrontal I/O Unit%<',
     $'IPUNLPs17iOutput Unit Number%<',
     $'LYNFNCs24iFunction Number (Line Search Diagnostic Plot)%<',
     $'LYNOUTs25iOutput Unit Number (Line Search Diagnostic Plot)%<',
     $'LYNPLTs26iIteration Number For Line Search Diagnostic Plot%<'/
c
      DATA (outvar(I),I=31,40) /
     $'LYNPNTs27iNo. of Plot Points (Line Search Diagnostic Plot)%<',
     $'LYNVARs28iVariable Number (Line Search Diagnostic Plot)%<',
     $'MAXLYNs29iMaximum Line Limit for Array Output%<',
     $'MAXNFEs30iMaximum Number of Function Evaluations%<',
     $'MUCALCb02iBarrier/Multiplier Initialization Option%<',
     $'MXQPITb03iMaximum Number of Barrier QP Iterations%<',
     $'NEWTONs32iNewton Option (0,1,2) = (Default,Newton,Gauss)%<',
     $'NITMAXs33iMaximum Number of Iterations%<',
     $'NITMINs34iMinimum Number of Iterations%<',
     $'ALGOPTs01cAlgorithm Control Option%>=  FM      Find Feasible Poin
     $t Then Minimize%>=  M       Minimize From The Initial Point%>=  F 
     $      Find Feasible Point Only%<'/
c
      DATA (outvar(I),I=41,nbrvar) /
     $'KTOPTNs02cKT Matrix Factorization Option%>=  SMALL   Condensed KK
     $T System%>=  LARGE   Full KKT System%<' /
C
c        nonzero ibrief(kk) to display variable in brief format
c
      DATA (ibrief(I),I=1,nbrvar) /
     $  0, 0, 0, 1, 0, 0, 1, 1, 0, 0,
     $  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
     $  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     $  0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 
     $  0/
C
      logical endvar,endlin
c
      CHARACTER(LEN=110) BLNKLN,LYNE
      DATA BLNKLN(1:2) / ' '/
      DATA BLNKLN(4:105) / ' '/
      DATA BLNKLN(3:3) / '*'/
      DATA BLNKLN(106:106) / '*'/
c
      character(len=1) algrm,vtype
C
C-------------------------------------------------------------
C
      CHARACTER(LEN=10) OUT1,OUT2
      LOGICAL SUMRY,full
C
C-------------------------------------------------------------
C
c         summary format flag
c
c           sumry = true    ----> abbreviated format
c           sumry = false   ----> full format
c
      SUMRY = .FALSE.
C
C         NOTE: IPUNLP = INPSPR(18) = IPUOPS
C
      IPUOPS = INPSPR(18)
c
      IF(modeb.eq.0) THEN
C
        SUMRY = .TRUE.
        WRITE(IPUOPS,1003) DFAWLT
C
      ELSEIF(modeb.EQ.1) THEN
C
C         SHORT OPTIONS TITLE
C
        WRITE(IPUOPS,1001) ALGNAM
C
      ELSEIF(modeb.eq.2) THEN
C
C         FULL OPTIONS TITLE
C
        WRITE(IPUOPS,1002) ALGNAM
C
      else
        print *,'invalid modeb = ',modeb
        stop
      ENDIF
c
      full = .not.sumry
c
c         set variable display flag
c
      ishow(1:nbrvar) = 1
      if(sumry) ishow(1:nbrvar) = ibrief(1:nbrvar)
C
C-------------------------------------------------------------
C-----Loop Over All Variables for Barrier Dispaly-------------
C-------------------------------------------------------------
c
      write(ipuops,'(a)') blnkln
c
      if(full) then
c
        lyne = blnkln
        lyne(9:99) = repeat('-',91)
        write(ipuops,'(a)') lyne
c
        write(ipuops,'(a)') blnkln
c
        lyne = blnkln
        lyne(9:14) = 'SYMBOL'
        lyne(23:29) = 'DEFAULT'
        lyne(40:44) = 'VALUE'
        lyne(50:61) = 'DESCRIPTION'
        write(ipuops,'(a)') lyne
c
        write(ipuops,'(a)') blnkln
c
        lyne = blnkln
        lyne(9:99) = repeat('-',91)
        write(ipuops,'(a)') lyne
c
        write(ipuops,'(a)') blnkln
c
      endif
c
      varloop: do ivar=1,nbrvar
c
c       ----separation between real/integer & integer/character
c
        if((ivar.eq.15.or.ivar.eq.40).and.full) then
c
          write(ipuops,'(a)') blnkln
c
          lyne = blnkln
          lyne(9:99) = repeat('-',91)
          write(ipuops,'(a)') lyne
c
          write(ipuops,'(a)') blnkln
c
        endif
c
c       ----blank line 
c
        lyne = blnkln
c
c       ----symbol name
c
        lyne(9:14) = outvar(ivar)(1:6)
c
c       ----algorithm type information
c           algrm = s  sqp     variable
c           algrm = b  barrier variable
c
        algrm = outvar(ivar)(7:7)
c
c       ----location of default value in default array
c
        read(outvar(ivar)(8:9),'(i2)') lcdflt 
c
c       ----variable type information
c           vtype = r  real (double precision) variable
c           vtype = i  integer variable
c           vtype = c  character variable
c
        vtype = outvar(ivar)(10:10)
c
C-------------------------------------------------------------
C-------------------------------------------------------------
C-------------------------------------------------------------
        select case (algrm)
        case('s')
c
C=============================================================
C=========SQP Algorithm=======================================
C=============================================================
c
          select case (vtype)
          case('r')
c
c           construct the variable and default strings
c
            WRITE(OUT1(1:10),'(1PG10.3)') RDFLT(LCDFLT)
            WRITE(OUT2(1:10),'(1PG10.3)') RNPSPR(LCDFLT)
            out1(1:10) = adjustr(out1(1:10))
            out2(1:10) = adjustr(out2(1:10))
            IF(OUT1.ne.OUT2) lyne(8:8) = '*'
            lyne(20:29) = out1(1:10)
            lyne(35:44) = out2(1:10)
c
          case('i')
c
            WRITE(OUT1(1:10),'(I10)') IDFLT(LCDFLT)
            WRITE(OUT2(1:10),'(I10)') INPSPR(LCDFLT+1)
            out1(1:10) = adjustr(out1(1:10))
            out2(1:10) = adjustr(out2(1:10))
            IF(OUT1.ne.OUT2) lyne(8:8) = '*'
            lyne(20:29) = out1(1:10)
            lyne(35:44) = out2(1:10)
c
          case('c')
c
            OUT1(1:10) = '          '
            OUT2(1:10) = OUT1(1:10)
            OUT1(5:10) = CDFLT(lcdflt)
            OUT2(5:10) = CNPSPR(lcdflt)
            out1(1:10) = adjustr(out1(1:10))
            out2(1:10) = adjustr(out2(1:10))
            IF(OUT1.ne.OUT2) lyne(8:8) = '*'
            lyne(20:29) = out1(1:10)
            lyne(35:44) = out2(1:10)
c
          end select
c
        case('b')
c
C=============================================================
C=========Barrier Algorithm===================================
C=============================================================
c
          select case (vtype)
          case('r')
c
c           construct the variable and default strings
c
            WRITE(OUT1(1:10),'(1PG10.3)') RBDFLT(LCDFLT)
            WRITE(OUT2(1:10),'(1PG10.3)') RBNPSP(LCDFLT)
            out1(1:10) = adjustr(out1(1:10))
            out2(1:10) = adjustr(out2(1:10))
            IF(OUT1.ne.OUT2) lyne(8:8) = '*'
            lyne(20:29) = out1(1:10)
            lyne(35:44) = out2(1:10)
c
          case('i')
c
            WRITE(OUT1(1:10),'(I10)') IBDFLT(LCDFLT)
            WRITE(OUT2(1:10),'(I10)') IBNPSP(LCDFLT)
            out1(1:10) = adjustr(out1(1:10))
            out2(1:10) = adjustr(out2(1:10))
            IF(OUT1.ne.OUT2) lyne(8:8) = '*'
            lyne(20:29) = out1(1:10)
            lyne(35:44) = out2(1:10)
c
          case('c')
c
            OUT1(1:10) = '          '
            OUT2(1:10) = OUT1(1:10)
            OUT1(5:10) = CDFLT(lcdflt)
            OUT2(5:10) = CNPSPR(lcdflt)
            out1(1:10) = adjustr(out1(1:10))
            out2(1:10) = adjustr(out2(1:10))
            IF(OUT1.ne.OUT2) lyne(8:8) = '*'
            lyne(20:29) = out1(1:10)
            lyne(35:44) = out2(1:10)
c
          end select
c
        end select
c
C-------------------------------------------------------------
C-------------------------------------------------------------
C-------------------------------------------------------------
c
c         display description line(s)
c
        jj = 10
        disloop: do kk = 49,106
          jj = jj + 1
          endlin = '%>'.eq.outvar(ivar)(jj+1:jj+2)
          endvar = '%<'.eq.outvar(ivar)(jj+1:jj+2)
          if(sumry.and.endlin) endvar = .true.
          lyne(kk:kk) = outvar(ivar)(jj:jj)
          if(endlin.or.endvar) exit disloop
        enddo disloop
c
c       ----display complete line; omit "brief" variables with default values
c
        if(ishow(ivar).eq.0.and.lyne(8:8).eq.' ') cycle varloop
c
        write(ipuops,'(a)') lyne
c
        if(endvar) cycle varloop
c
        dlinloop: do ll = 1,nmlpvr
c        
          lyne = blnkln
          jj = jj + 2

          lynloop: do kk = 49,106
            jj = jj + 1
            endlin = '%>'.eq.outvar(ivar)(jj+1:jj+2)
            endvar = '%<'.eq.outvar(ivar)(jj+1:jj+2)
            lyne(kk:kk) = outvar(ivar)(jj:jj)
            if(endlin.or.endvar) exit lynloop
          enddo lynloop
c
c       ----display complete line
c
          write(ipuops,'(a)') lyne
c
          if(endvar) cycle varloop
c
        enddo dlinloop
c
      enddo varloop
C
 1001 FORMAT(T3,'*',T106,'*'/T3,'*',T33,A8,  ' OPTIONS  (ALSO SEE ''FULL
     $'' OPTIONS)',T106,'*')
 1002 FORMAT(T3,'*',T106,'*'/T3,'*',T33,A8,   ' FULL OPTIONS  (ALSO SEE 
     $''OPTIONS'')',T106,'*')
 1003 FORMAT(T3,'*',T106,'*'/T3,'*',T11,'ALGORITHM CONTROL PARAMETERS',
     * 2X,A17,T106,'*')
      RETURN
      END
