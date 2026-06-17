      INTEGER FUNCTION JHSN74(IDUM)
C
C ======================================================================
C     JHSN74===>innset   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  THIS ROUTINE IS USED IN LIEU OF BLOCK DATA TO 
C                   INITIALIZE A SINGLE INTEGER VARIABLE.  IT IS USED
C                   TO INITIATE THE DEFAULT PROCEDURES FOR THE SPRNLP
C                   COMMON 
C
C         ARGUMENTS:  NONE (DUMMY ONLY)
C
      SAVE INTJER
      DATA INTJER /-1/
      JHSN74 = INTJER
      INTJER = 0
C
      RETURN
      END
