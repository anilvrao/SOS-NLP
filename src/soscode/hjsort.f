      SUBROUTINE HJSORT ( A, N, IER )
C
C----------------------------------------------------------------------
C
C  PURPOSE  HJSORT SORTS AN INTEGER ARRAY OF N ELEMENTS INTO THE
C           STANDARD COLLATING SEQUENCE.
C
C  METHOD   HJSORT USES A QUICKSORT ALGORITHM IN THE STYLE OF THE
C           CACM PAPER BY BOB SEDGEWICK, OCTOBER 1978
C
C  INPUT    N        NUMBER OF ELEMENTS TO BE SORTED.
C
C  INPUT/   A        ARRAY THAT IS TO BE SORTED.
C    OUTPUT
C
C  OUTPUT   IER      SUCCESS/ERROR FLAG.
C
C              IER =  0    NORMAL RETURN
C              IER = -1    N .LE. 0
C              IER = -5    INTERNAL FAILURE
C
C----------------------------------------------------------------------
C
C
C----------------------------------------------------------------------
C
C ... GLOBAL VARIABLES
C
      INTEGER        N, IER
C
      INTEGER        A(*)
C
C----------------------------------------------------------------------
C
C ... LOCAL VARIABLES
C
      INTEGER     I, IP1, J, JM1, LEFT, LLEN, RIGHT, RLEN,
     1            STACK(50), STKLEN, TINY, TOP
C
C
      INTEGER     ATEMP(2)
      CHARACTER(LEN=8)  NAME
C
      DATA        STKLEN / 50 /,  TINY / 9 /,
     A            NAME / 'HJSORT'/
C
C----------------------------------------------------------------------
C
C ... CHECK INPUT
C
      IER = 0
C
      IF ( N .LE. 0 ) THEN
C
         IER = -1
         CALL HHERR ( 1, NAME, IER, 0 )
         GO TO 9000
      ENDIF
C
C----------------------------------------------------------------------
C
C ... PROGRAM IS A DIRECT TRANSLATION INTO FORTRAN OF SEDGEWICK'S
C     PROGRAM 2, WHICH IS NON-RECURSIVE, IGNORES FILES OF LENGTH
C     LESS THAN 'TINY' DURING PARTITIONING, AND USES MEDIAN OF THREE
C     PARTITIONING.
C
      IF (N .EQ. 1)  GO TO 9000
C
      TOP = 1
      LEFT = 1
      RIGHT = N
      IF ( N .LE. TINY ) GO TO 2000
C
C     ===========================================================
C     QUICKSORT -- PARTITION THE FILE UNTIL NO SUBFILE REMAINS OF
C     LENGTH GREATER THAN 'TINY'
C     ===========================================================
C
C     ... WHILE NOT DONE DO ...
C
C
C         ... FIND MEDIAN OF LEFT, RIGHT AND MIDDLE ELEMENTS OF CURRENT
C             SUBFILE, WHICH IS  A(LEFT), ..., A(RIGHT)
C             (CORRECTION TO CACM ARTICLE INCORPORATED)
C
  300     CONTINUE
          I        = ( LEFT + RIGHT ) / 2
          ATEMP(1) = A (I)
          A (I)    = A (LEFT)
          A (LEFT) = ATEMP(1)
C
          IF (  A(LEFT+1) .GT. A(LEFT)  ) THEN
              ATEMP(1)       = A (LEFT+1)
              A (LEFT+1)     = A (LEFT)
              A (LEFT)       = ATEMP(1)
          ENDIF
C
          IF (  A(LEFT) .GT. A(RIGHT)  ) THEN
              ATEMP(1)      = A (LEFT)
              A (LEFT)      = A (RIGHT)
              A (RIGHT)     = ATEMP(1)
C
              IF (  A (LEFT+1) .GT. A (LEFT)  )  THEN
                  ATEMP(1)       = A (LEFT+1)
                  A (LEFT+1)     = A (LEFT)
                  A (LEFT)       = ATEMP(1)
              ENDIF
C
          ENDIF
          ATEMP(2)   = A (LEFT)
C
C         ... ATEMP(2) IS NOW THE MEDIAN VALUE OF THE THREE A-S.  NOW
C             MOVE FROM THE LEFT AND RIGHT ENDS SIMULTANEOUSLY,
C             EXCHANGING A-S UNTIL ALL A-S LESS THAN  ATEMP(2) ARE
C             PACKED TO THE LEFT, ALL A-S LARGER THAN  ATEMP(2) ARE
C             PACKED TO THE RIGHT.
C
          I = LEFT+1
          J = RIGHT
C
C         LOOP
C             REPEAT I = I+1 UNTIL A(I) >= ATEMP(2)  ;
C             REPEAT J = J-1 UNTIL A(J) <= ATEMP(2)  ;
C         EXIT IF J < I;
C             << EXCHANGE AS I AND J >>
C         END
C
  700     CONTINUE
  800         CONTINUE
              I  = I + 1
              IF ( A(I) .LT. ATEMP(2)   )  GO TO 800
C
  900         CONTINUE
              J = J - 1
              IF ( A(J) .GT. ATEMP(2)   )  GO TO 900
C
          IF (J .GE. I)  THEN
              ATEMP(1)  = A (I)
              A (I)     = A (J)
              A (J)     = ATEMP(1)
              GO TO 700
          ENDIF
C
          ATEMP(1)     = A (LEFT)
          A (LEFT)     = A (J)
          A (J)        = ATEMP(1)
C
C
C         ... WE HAVE NOW PARTITIONED THE FILE INTO TWO SUBFILES,
C             ONE IS (LEFT ... J-1)  AND THE OTHER IS (I...RIGHT).
C             PROCESS THE SMALLER NEXT.  STACK THE LARGER ONE.
C
          LLEN = J-LEFT
          RLEN = RIGHT - I + 1
          IF ( RLEN .GT. LLEN ) GO TO 1100
C
C         ... LEFT SUBFILE ( LEFT ... J-1 ) IS THE LARGE SUBFILE.
C             TEST IF IT IS SMALL ENOUGH NOT TO PROCESS FURTHER.
C
              IF ( LLEN .LE. TINY ) GO TO 1200
C
C             ... IF RIGHT SUBFILE IS SMALL THEN PROCESS THE LEFT.
C                 ELSE STACK THE LEFT AND PROCESS THE RIGHT.
C
                  IF ( RLEN .LE. TINY ) THEN
                       RIGHT = J - 1
                       GO TO 300
                  ENDIF
C
                       IF ( TOP .GE. STKLEN ) GO TO 8000
                           STACK ( TOP    ) = LEFT
                           STACK ( TOP+1  ) = J - 1
                           TOP              = TOP + 2
                           LEFT             = I
                           GO TO 300
C
C        ... RIGHT SUBFILE ( I ... RIGHT ) IS THE LARGE SUBFILE.
C            TEST IF IT IS SMALL ENOUGH NOT TO PROCESS FURTHER.
C
 1100        CONTINUE
             IF ( RLEN .LE. TINY ) GO TO 1200
C
C            ... IF LEFT SUBFILE IS SMALL THEN PROCESS THE RIGHT.
C                ELSE STACK THE RIGHT AND PROCESS THE LEFT.
C
                 IF ( LLEN .LE. TINY ) THEN
                     LEFT = I
                     GO TO 300
                 ENDIF
C
                     IF ( TOP .GE. STKLEN ) GO TO 8000
                         STACK ( TOP     ) = I
                         STACK ( TOP + 1 ) = RIGHT
                         TOP               = TOP + 2
                         RIGHT             = J - 1
                         GO TO 300
C
C        ... BOTH LEFT AND RIGHT SUBFILE ARE SMALL.  POP THE STACK
C            TO GET THE NEXT SUBFILE TO PROCESS.
C
 1200        CONTINUE
             IF ( TOP .NE. 1 ) THEN
                 TOP    = TOP - 2
                 LEFT   = STACK ( TOP )
                 RIGHT  = STACK ( TOP + 1 )
                 GO TO 300
             ENDIF
C
C     ------------------------------------------------------------
C     INSERTION SORT THE ENTIRE FILE, WHICH CONSISTS OF A LIST
C     OF 'TINY' SUBFILES, LOCALLY OUT OF ORDER, GLOBALLY IN ORDER.
C     ------------------------------------------------------------
C
C     ... INSERTION SORT ... FOR I := N-1 STEP -1 TO 1 DO ...
C
 2000 CONTINUE
      I   = N - 1
      IP1 = N
C
 2100     CONTINUE
          IF (  A (I) .GT. A (IP1)  )  THEN
C
C             ... OUT OF ORDER ... MOVE UP TO CORRECT PLACE
C
              ATEMP(1)  = A (I)
              J         = IP1
              JM1       = I
C
C             ... REPEAT ... UNTIL 'CORRECT PLACE FOR ATEMP(1)   FOUND'
C
 2200             CONTINUE
                  A (JM1)     = A (J)
                  JM1         = J
                  J           = J + 1
                  IF ( J .LE. N ) THEN
                    IF (  A (J) .LT. ATEMP(1)    )  GO TO 2200
                  ENDIF
C
              A (JM1)     = ATEMP(1)
          ENDIF
C
          IP1 = I
          I   = I - 1
          IF ( I .GT. 0 )  GO TO 2100
C
      GO TO 9000
C
C---------------------------------------------------------------------
C
 8000 CONTINUE
      IER = -5
      CALL HHERR ( 3, NAME, IER, 0 )
C
C----------------------------------------------------------------------
C
C ... END OF HJSORT
C
 9000    CONTINUE
      RETURN
      END
