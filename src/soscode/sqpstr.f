
      SUBROUTINE SQPSTR(MCON,MINEQL,NDIM,NONZH,NONZG,IWRKBD,IRWRKB)
C
C ======================================================================
C     SQPSTR===>sqpstr   J.T. BETTS
C ======================================================================
C
C     ==================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:   COMPUTE AN ESTIMATE FOR THE REAL AND INTEGER
C                    WORK ARRAY STORAGE.  THESE ARE LOWER BOUNDS.
C
C         INPUT:     
C
C           MCON     NUMBER OF CONSTRAINTS
C           MINEQL   NUMBER OF INEQUALITY CONSTRAINTS
C           NDIM     NUMBER OF VARIABLES
C           NONZH    NONZEROS IN THE HESSIAN
C           NONZG    NONZEROS IN THE JACOBIAN
C
C         OUTPUT:
C
C           IWRKBD   LOWER BOUND FOR THE INTEGER WORK ARRAY SIZE
C           IRWRKB   LOWER BOUND FOR THE REAL WORK ARRAY SIZE
C
C
C
C     KTOPTN DETERMINES WHETHER TO USE OLD OR NEW SPARSE QP. THE NEW
C     SPARSE QP HAS ONLY ACTIVE CONSTRS IN KKT MATRIX.
C
C-------------------------------------------------------------
      INCLUDE '../commons/NLPSPR.CMN'
C-------------------------------------------------------------
C
C
      IF ( KTOPTN .EQ. 'SMALL '  ) THEN
C       NEQN <-- MAX SIZE OF KKT SYSTEM.
C
        NEQN =  1 + 2*NDIM
C
C       NEQBIG <-- MAX SIZE KKT SYTEM WOULD HAVE IF FREE SLACKS AND
C                THE CORRESPONDING INACTIVE INEQUALITIES WERE INCLUDED.
C
        NEQBIG = 1 + MINEQL + NDIM + MCON
C
        XNDIM = NDIM
        XNSQR = XNDIM**2
        XNONZG = NONZG
        NKTNZ = 1 + NDIM + NONZH + MIN(XNONZG,XNSQR)
C
C           QPSKT STORAGE REQUIREMENTS:
C
        IQPSKT = 0
        IRQPSK = 2*NKTNZ + 2*NEQN + 200
C
C           QPSALG STORAGE REQUIREMENTS:
C
        KMAX = MIN(100,2*MAX(NDIM,MCON))
        NSTD = MINEQL + 1 + NDIM
        XMCON = MAX(1,MCON)
        XNDMC = XNDIM*XMCON
        XKMAX = KMAX
        XNEQN = NEQN
        XNONZH = NONZH
        NONZU = XKMAX*XNEQN*(XNONZH/XNSQR + 
     $      XNONZG/XNDMC)*1.1D0 + MIN(MCON,NDIM)
        NONZU = MAX(1,NONZU)
C      
        IQPSLG = 4*NSTD + 4*KMAX + 1 + 2*NEQN + 2*MCON + NONZU
        IRQPSL = 6*NSTD + 5*KMAX + 2*KMAX**2 + 6*NEQN + 
     $           NEQBIG + NONZU + MCON
C
C             SHURQP STORAGE REQUIREMENTS:
C
        ISHRQP = 2*MCON + NDIM + MAX(2*MCON,NONZG,NDIM) + NONZG + 1
        IRSHRQ = 3*NSTD + 2*MCON + 1
C 
C             SHURQP STORAGE TOTALS:
C
        IWRKBD = IQPSKT + IQPSLG + ISHRQP
        IRWRKB = IRQPSK + IRQPSL + IRSHRQ
C
      ELSEIF( KTOPTN .EQ. 'LARGE '  ) THEN
C
        NEQN = MINEQL + 1 + NDIM + MCON
C
        NKTNZ = MINEQL + 1 + MCON + NONZH + NONZG
C
C           QPSKT STORAGE REQUIREMENTS:
C
        IQPSKT = 0
        IRQPSK = 2*NKTNZ + 2*NEQN + 200
C
C           QPSALG STORAGE REQUIREMENTS:
C
        KMAX = MIN(100,2*MAX(NDIM,MCON))
        NSTD = MINEQL + 1 + NDIM
        XNDIM = NDIM
        XNSQR = XNDIM**2
        XMCON = MAX(1,MCON)
        XNDMC = XNDIM*XMCON
        XKMAX = KMAX
        XNEQN = NEQN
        XNONZH = NONZH
        XNONZG = NONZG
        NONZU = XKMAX*XNEQN*(XNONZH/XNSQR + 
     $      XNONZG/XNDMC)*1.1D0 + MIN(MCON,NDIM)
        NONZU = MAX(1,NONZU)
C      
        IQPSLG = 4*NSTD + 4*KMAX + 1 + 2*NEQN + NONZU
        IRQPSL = 6*NSTD + 5*KMAX + 2*KMAX**2 + 7*NEQN + NONZU + MCON
C
C           SHURQP STORAGE REQUIREMENTS:
C
        ISHRQP = 2*MCON + NDIM + MAX(2*MCON,NONZG,NDIM) + NONZG + 1
        IRSHRQ = 3*NSTD + 2*MCON + 1
C
C           SHURQP STORAGE TOTALS:
C
        IWRKBD = IQPSKT + IQPSLG + ISHRQP
        IRWRKB = IRQPSK + IRQPSL + IRSHRQ
      ELSE
        PRINT *,'INVALID KTOPTN =',KTOPTN
        STOP
      ENDIF
C      
C     ===================================================================
C
C     ----------
C     ... RETURN
C     ----------
C
      RETURN
      END
