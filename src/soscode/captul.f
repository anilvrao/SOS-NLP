
      SUBROUTINE CAPTUL(STRING,LNSTRN,IER)
C
C ======================================================================
C     CAPTUL===>captul   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  CONVERT THE INPUT CHARACTER STRING
C                   OF LENGTH LNSTRN TO CAPITAL LETTERS
C         Ref: FORTAN 90 Programming, 
C         T.M.R. Ellis, Ivor R. Philips, Thomas M. Lahey, p. 150
C
C         INPUT:
C
C            STRING CHARACTER VARIABLE OF LENGTH LNSTRN
C
C         OUTPUT:
C
C            STRING CAPITALIZED VARIABLE STRING
C            IER    ERROR FLAG
C                   NONZERO FOR REPLACEMENT ERROR
C
      CHARACTER(LEN=*) STRING
C     --local constant
      INTEGER, PARAMETER :: UP2LOW = IACHAR("a")-IACHAR("A")
C
      IER = 0
C
      DO KK = 1,LNSTRN
C
C         Check if argument is lower case alphabetic, upper case 
C         alphabetic, or non-alphabetic.  Do NOT change the
C         case of the following statement
C
        IF("a"<=STRING(KK:KK) .AND. STRING(KK:KK)<="z") THEN
C         lower case ---> convert to upper case
          STRING(KK:KK) = ACHAR(IACHAR(STRING(KK:KK))-UP2LOW)
        ENDIF
C
      ENDDO
C
      RETURN
      END
