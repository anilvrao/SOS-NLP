

      SUBROUTINE BIGGLM(AMAT,IROWA,JCOLA,NONZA,ISTATC,NROW,
     $    BIGELM,IROWBG,JCOLBG,LENBIG,INDEX,IERBIG)
C
C ======================================================================
C     BIGGLM===>bigelm   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  COMPUTE THE LARGEST ELEMENTS IN A SPARSE MATRIX
C
C         ARGUMENTS:
C
C            AMAT   ARRAY CONTAINING FIRST MATRIX ELEMENTS
C            IROWA  ROW INDICES FOR AMAT
C            JCOLA  COLUMN INDICES FOR AMAT
C            NONZA  LENGTH OF AMAT, IROWA, JCOLA
C            ISTATC INTEGER ARRAY DEFINING IGNORED ROWS (NROW)
C                   = 4 FOR IGNORED CONSTRAINT
C            NROW   NUMBER OF ROWS IN AMAT
C            BIGELM ARRAY CONTAINING LARGEST (MAGNITUDE) ELEMENTS IN AMAT 
C                   (LENBIG)
C            IROWBG ROW NUMBER FOR BIGGEST ELEMENTS (LENBIG)
C            JCOLBG COLUMN NUMBER FOR BIGGEST ELEMENTS (LENBIG)
C            INDEX  INTEGER WORK ARRAY (LENBIG)
C            LENBIG LENGTH OF BIGGEST ELEMENT ARRAY
C            IERBIG NONZERO ERROR RETURN FLAG
C
C
C     *******************************************************
C
      PARAMETER (ZERO=0.0D0,ONE=1.D0)
      DIMENSION AMAT(NONZA),IROWA(NONZA),JCOLA(NONZA),ISTATC(NROW),
     $          BIGELM(LENBIG),IROWBG(LENBIG),JCOLBG(LENBIG),
     $          INDEX(LENBIG)
C
C     *******************************************************
C
      DO I=1,LENBIG
        BIGELM(I) = ZERO
        IROWBG(I) = 0
        JCOLBG(I) = 0
      ENDDO
      LBIGST = MIN(LENBIG,NONZA)
      DO I=1,NONZA
        IF(IROWA(I).LE.0.OR.IROWA(I).GT.NROW) THEN
          IERBIG = -147
          RETURN
        ENDIF
      ENDDO
C
      DO II = 1,LBIGST
C
        ABIG = ZERO
        IBIG = 0
        DO I = 1,NONZA
          IROW = IROWA(I)
          IF(ISTATC(IROW).NE.4) THEN
            ABSA = ABS(AMAT(I))
            IF(ABSA.GT.ABIG) THEN
              ABIG = ABSA
              IBIG = I
            ENDIF
          ENDIF
        ENDDO
C
        IF(IBIG.EQ.0) THEN
          LBIGST = II - 1
          EXIT
        ENDIF
        BIGELM(II) = AMAT(IBIG)
        IROWBG(II) = IROWA(IBIG)
        JCOLBG(II) = JCOLA(IBIG)
        INDEX(II)  = IBIG
        AMAT(IBIG) = ZERO
C
      ENDDO
C
      DO II = 1,LBIGST
C
        IBIG = INDEX(II)
        AMAT(IBIG) = BIGELM(II)
C
      ENDDO
C
      RETURN
      END
