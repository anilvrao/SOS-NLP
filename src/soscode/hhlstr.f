      SUBROUTINE HHLSTR ( A, B, MODE, IPOSA )
C
C----------------------------------------------------------------------
C
C ... PURPOSE     LOCATE THE CHARACTER STRING B IN THE CHARACTER
C                 STRING A.
C
C ... INPUT       A        THE CHARACTER VARIABLE TO BE SCANNED.
C
C                 B        CHARACTER VARIABLE CONTAINING THE STRING
C                          TO BE LOCATED.
C
C                 MODE     DIRECTION OF SEARCH MODE.
C                          MODE .GE. 0    SEARCH LEFT TO RIGHT IN A
C                          MODE .LT. 0    SEARCH RIGHT TO LEFT IN A
C
C ... OUTPUT      IPOSA    FIRST POSITION IN A WHERE THE STRING IN B
C                          BEGINS.  IPOSA = 0 IF B WAS NOT FOUND IN A
C                          OR IF THE LENGTH OF B WAS GREATER THAN THE
C                          LENGTH OF A.
C
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C
C ... GLOBAL VARIABLES
C
         CHARACTER(LEN=*)  A, B
C
         INTEGER           MODE, IPOSA
C
C----------------------------------------------------------------------
C
C ... LOCAL VARIABLES
C
         INTEGER           IA, IAEND, IBGN, IEND, INCRMT, LENA, LENB,
     A                     LENBM1
C
C----------------------------------------------------------------------
C
C
C ... INITIALIZE
C
         LENA = LEN ( A )
         LENB = LEN ( B )
         IF ( LENB .GT. LENA ) THEN
C
C----------------------------------------------------------------------
C
C ...... SPECIAL RETURN FOR CASE - LENB .GT. LENA.
C
            IPOSA = 0
            GO TO 900
         ENDIF
C
C----------------------------------------------------------------------
C
C ... USE FORTRAN-77 INTRINSIC IF FORWARD SEARCH
C
         IF ( MODE .GE. 0 ) THEN
C
            IPOSA = INDEX ( A, B )
            GO TO 900
         ENDIF
C
C----------------------------------------------------------------------
C
C ... SET THE DO LOOP PARAMETERS FOR THE BACKWARD SEARCH THROUGH A
C
         LENBM1 =  LENB - 1
         IBGN   =  LENA - LENBM1
         IEND   =  1
         INCRMT = -1
C
C ... START OF THE LOOP THROUGH A
C
         DO IA = IBGN, IEND, INCRMT
C
            IAEND = IA + LENBM1
            IF (  B .EQ. A(IA:IAEND)  )  GO TO 400
C
         ENDDO
C
C ...... STRING B NOT FOUND IN A
C
            IPOSA  = 0
            GO TO 900
C
C ...... STRING B FOUND IN A
C
  400     CONTINUE
            IPOSA  = IA
C
C----------------------------------------------------------------------
C
C ... END OF HHLSTR
C
  900    CONTINUE
      RETURN
      END
