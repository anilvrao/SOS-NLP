

      SUBROUTINE EVSYMR(SYMBOL,VALSYM,NSYM,SYMR,RELSYM,FOUND,IER)
C
C ======================================================================
C     EVSYMR===>evsymr   J.T. BETTS
C ======================================================================
C
C
C     PURPOSE:  Find SYMBOL in the list SYMR and load the value in
C               VALSYM string into RELSYM.  WARNING: The length
C               of the strings SYMBOL and SYMR should be the same.
C
C     INPUT:
C
C        SYMBOL    CHARACTER STRING CONTAINING SYMBOL.
C        VALSYM    CHARACTER STRING CONTAINING VALUE.
C        SYMR      LIST OF VALID NAMES.
C
C     OUTPUT:
C
C        RELSYM    ARRAY OF VALUES TO BE UPDATED.
C        FOUND     LOGICAL SPECIFYING WHETHER SYMBOL WAS FOUND IN SYMI.
C        IER       ERROR RETURN FLAG.
C
C     ------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER(LEN=*) SYMBOL,VALSYM
      CHARACTER(LEN=6) SYMR(NSYM)
      CHARACTER(LEN=1) TYPE
      DIMENSION  RELSYM(NSYM)
      LOGICAL  FOUND
      REAL SVAL
C
C     ------------------------------------------------------------------
C
      IER = 0
      FOUND = .FALSE.
C
C             COMPARE THE INPUT SYMBOL AGAINST THE LIST OF POSSIBLE
C             CANDIDATES.
C
      DO ISYM = 1,NSYM
        IF(SYMBOL.EQ.SYMR(ISYM)) THEN
          FOUND = .TRUE.
C
C             NOW FIND THE NUMERIC VALUE IN VALSYM.
C             AFTER THE CALL THE NUMERIC VALUE IS FOUND
C             IN VALSYM(LCNRL:LCNRR), AND ITS TYPE IS 
C             GIVEN BY THE VARIABLE TYPE.
C
          CALL HHFCNR(VALSYM,' ',LCNRL,LCNRR,TYPE)
C
C             CONVERT THE CHARACTER VARIABLE TO ITS NUMERIC
C             REPRESENTATION AS INDICATED BY THE TYPE FLAG.
C
          IF (TYPE.EQ.'F' .OR. TYPE.EQ.'E')
     &      THEN
C
              CALL HHCNSP(VALSYM(LCNRL:LCNRR),SVAL,IERCNV)
              DVAL = SVAL
              IF (IERCNV.EQ.0)
     &          THEN
                  RELSYM(ISYM) = DVAL
                ELSE
                  IER = -1
              ENDIF
C
            ELSEIF(TYPE.EQ.'D') THEN
C
              CALL HHCNDP(VALSYM(LCNRL:LCNRR),DVAL,IERCNV)
              IF (IERCNV.EQ.0)
     &          THEN
                  RELSYM(ISYM) = DVAL
                ELSE
                  IER = -1
              ENDIF
C
            ELSE
              IER = -2
          ENDIF
          RETURN
        ENDIF
      ENDDO
C
      RETURN
      END
