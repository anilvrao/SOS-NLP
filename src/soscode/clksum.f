
      SUBROUTINE CLKSUM(NCLOCK)
C
C ======================================================================
C     CLKSUM===>clksum   J.T. BETTS
C ======================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C       PURPOSE:  ACCUMULATE THE SUM FOR CLOCK NUMBER NCLOCK
C                 THE CALL TO THIS ROUTINE MUST BE PRECEDED BY
C                 (1) A CALL TO CLKBEG, AND
C                 (2) A CALL TO CLKSET WITH NCLOCK > 0
C
C       INPUT:
C
C         NCLOCK  NUMBER OF THE CLOCK TO BE ACCUMULATED
C
C                 NOTE IT IS ASSUMED THAT
C                 0 < NCLOCK LEQ MAXCLK
C                 THIS CONDITION IS NOT CHECKED!
C
      PARAMETER (MAXCLK=40)
C
      DOUBLE PRECISION INTCLK,INTOLD
      COMMON /CLKCOM/ CLKOLD(MAXCLK),CLKVAL(MAXCLK),INTCLK(MAXCLK),
     $    INTOLD(MAXCLK),IERCLK(MAXCLK)
C
      IF(IERCLK(NCLOCK).EQ.+1) THEN
        IERCLK(NCLOCK) = -1
        CALL XDSLT1(TYME,WLBEGN)
        DELCLK = TYME - CLKOLD(NCLOCK)
        CLKVAL(NCLOCK) = CLKVAL(NCLOCK) + DELCLK
        INTCLK(NCLOCK) = INTCLK(NCLOCK) + 1
      ELSE
        CLKVAL(NCLOCK) = -NCLOCK
        INTCLK(NCLOCK) = -NCLOCK
      ENDIF
C
      RETURN
      END
