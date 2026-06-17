

      SUBROUTINE EVSYMC(SYMBOL,VALSYM,NSYM,SYMC,CHRSYM,FOUND,IER)
C
C ======================================================================
C     EVSYMC===>evsymc   J.T. BETTS
C ======================================================================
C
C
C     PURPOSE:  Find SYMBOL in the list SYMC and load the value in
C               VALSYM string into CHRSYM.  WARNING: The length
C               of the strings SYMBOL and SYMC should be the same.
C
C     INPUT:
C
C        SYMBOL    CHARACTER STRING CONTAINING SYMBOL.
C        VALSYM    CHARACTER STRING CONTAINING VALUE.
C        SYMC      LIST OF VALID NAMES.
C
C     OUTPUT:
C
C        CHRSYM    ARRAY OF VALUES TO BE UPDATED.
C        FOUND     LOGICAL SPECIFYING WHETHER SYMBOL WAS FOUND IN SYMI.
C        IER       ERROR RETURN FLAG.
C
C     ------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER(LEN=*)  SYMBOL,VALSYM
      CHARACTER(LEN=6)  SYMC(NSYM)
      CHARACTER(LEN=1) TYPE
      CHARACTER(LEN=*)  CHRSYM(NSYM)
      LOGICAL  FOUND
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
        IF(SYMBOL.EQ.SYMC(ISYM)) THEN
          FOUND = .TRUE.
C
C             NOW FIND THE CHARACTER VALUE IN VALSYM.
C             AFTER THE CALL THE VALUE IS FOUND
C             IN VALSYM(LCNRL:LCNRR), AND ITS TYPE IS 
C             GIVEN BY THE VARIABLE TYPE.
C
          CALL HHFCNR(VALSYM,' ',LCNRL,LCNRR,TYPE)
C
C             CONVERT THE CHARACTER VARIABLE TO ITS NUMERIC
C             REPRESENTATION AS INDICATED BY THE TYPE FLAG.
C
          IF (TYPE.EQ.'N')
     &      THEN
C
C         COMPRESS THE BLANKS OUT OF THE STRING TO IDENTIFY THE
C         SYMBOL.
C
              LNVAL = LEN(VALSYM)
              CALL HHSDEL(VALSYM,' ',' ',LNVAL,ND,IERCNV)
C
              CHRSYM(ISYM) = VALSYM(1:LNVAL-ND)
            ELSE
              IER = -2
          ENDIF
          RETURN
        ENDIF
      ENDDO
C
      RETURN
      END
