
      SUBROUTINE HHISCN ( A, B, MODE, IPOSA, IPOSB )
C
C----------------------------------------------------------------------
C
C ... PURPOSE     SCAN THE CHARACTER VARIABLE A FOR THE FIRST OCCURENCE
C                 OF ANY CHARACTER IN THE CHARACTER VARIABLE B.
C
C ... INPUT       A        THE CHARACTER VARIABLE TO BE SCANNED.
C
C                 B        CHARACTER VARIABLE CONTAINING THE POOL OF
C                          CHARACTERS TO BE SCANNED FOR.
C
C                 MODE     DIRECTION OF SEARCH MODE.
C                          MODE .GE. 0    SEARCH LEFT TO RIGHT IN A
C                          MODE .LT. 0    SEARCH RIGHT TO LEFT IN A
C
C ... OUTPUT      IPOSA    FIRST POSITION IN A OF THE CHARACTER IN B
C                          THAT WAS FOUND IN A.  IPOSA = 0 IF NONE
C                          FOUND.
C
C                 IPOSB    POSITION OF THE CHARACTER IN B WHICH WAS
C                          FOUND FIRST IN A.  IPOSB = 0 IF NONE
C                          FOUND.
C
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C
C ... GLOBAL VARIABLES
C
         CHARACTER(LEN=*)  A, B
C
         INTEGER           MODE, IPOSA, IPOSB
C
C----------------------------------------------------------------------
C
C ... LOCAL VARIABLES
C
         INTEGER           IA, IBGN, IEND, INCRMT, JB, LENA, LENB
C
C----------------------------------------------------------------------
C
C ... INITIALIZE
C
         LENA     = LEN ( A )
         LENB     = LEN ( B )
C
C ... SET THE DO LOOP PARAMETERS FOR THE SEARCH THROUGH A
C
         IF ( MODE .GE. 0 ) THEN
            IBGN     = 1
            IEND     = LENA
            INCRMT   = 1
C
         ELSE
            IBGN     = LENA
            IEND     = 1
            INCRMT   = -1
         ENDIF
C
C ... START OF THE LOOP THROUGH A
C
            DO IA = IBGN, IEND, INCRMT
C
C ...... LOOP THROUGH B
C
               DO JB = 1, LENB
                  IF (  B(JB:JB)  .EQ.  A(IA:IA)  )  THEN
C
C ... CHARACTER IN B FOUND IN A
C
                    IPOSA  = IA
                    IPOSB  = JB
                    GO TO 900
                  ENDIF
               ENDDO
C
            ENDDO
C
C ... NO CHARACTER IN B FOUND IN A
C
            IPOSA  = 0
            IPOSB  = 0
C
C ... END OF HHISCN
C
  900    CONTINUE
         RETURN
      END
