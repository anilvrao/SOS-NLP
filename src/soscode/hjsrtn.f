      SUBROUTINE HJSRTN ( A, N, IOPT, ISTBLE, INDX, IER )
C
C----------------------------------------------------------------------
C
C  PURPOSE  HJSRTN SORTS A INTEGER ARRAY OF N ELEMENTS INTO THE
C           ASCENDING ORDER AND GENERATES THE PERMUTATION
C           ARRAY INDX.
C
C  METHOD   HJSRTN USES A QUICKSORT ALGORITHM IN THE STYLE OF THE
C           CACM PAPER BY BOB SEDGEWICK, OCTOBER 1978
C
C  INPUT    N        NUMBER OF ELEMENTS TO BE SORTED.
C
C           IOPT     SORTING OPTION.  IF IOPT = 0, THE ARRAY A IS
C                    RETURNED IN SORTED ORDER.  OTHERWISE, A IS
C                    RETURNED IN ORIGINAL ORDER.
C
C           ISTBLE   STABILIZATION OPTION.  IF ISTBLE = 0, THE
C                    SORT IS STABILIZED.  OTHERWISE, IT IS NOT.
C
C  INPUT/   A        ARRAY THAT IS TO BE SORTED.  IF IOPT = 0
C    OUTPUT          AND IER = 0 A IS RETURNED IN SORTED ORDER
C                    IF IER = -5, A IS PARTIALLY SORTED.  OTHERWISE
C                    A IS LEFT UNCHANGED.
C
C  OUTPUT   INDX      PERMUTATION ARRAY.
C
C           IER      SUCCESS/ERROR FLAG.
C
C              IER =  0    NORMAL RETURN
C              IER = -1    N .LE. 0
C              IER = -5    INTERNAL FAILURE
C
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C
C ... GLOBAL VARIABLES
C
      INTEGER        N, INDX(*), IER, IOPT, ISTBLE
C
      INTEGER        A(*)
C
C----------------------------------------------------------------------
C
C ... LOCAL VARIABLES
C
      INTEGER     I, IP1, J, JM1, K, L, LEFT, LLEN, RIGHT,
     A            RLEN, STACK(50), TOP,STKLEN,TINY
C
C
      CHARACTER(LEN=8)  NAME
      INTEGER     ATEMP(2)
C
      DATA        STKLEN, TINY / 50, 9 /,
     A            NAME / 'HJSRTN' /
C
C----------------------------------------------------------------------
C
C ... CHECK INPUT
C
      IER = 0
C
      IF ( N .LE. 0 ) IER = - 1
C
      IF ( IER .NE. 0 ) THEN
C
         INDX(1) = 0
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
C ... SET INDX TO THE INDENTITY PERMUTATION
C
      DO I = 1, N
         INDX(I) = I
      ENDDO
C
      IF (N .EQ. 1)  GO TO 9000
C
      TOP = 1
      LEFT = 1
      RIGHT = N
      IF ( N .LE. TINY) GO TO 2000
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
          I            = ( LEFT + RIGHT ) / 2
C
          ATEMP ( 1 )  = A(I)
          A(I)         = A ( LEFT )
          A ( LEFT )   = ATEMP ( 1 )
C
          K            = INDX (I)
          INDX (I)      = INDX (LEFT)
          INDX (LEFT)   = K
C
          IF ( A(LEFT+1) .GT. A(LEFT) ) THEN
C
              ATEMP ( 1 )  = A ( LEFT+1 )
              A ( LEFT+1 ) = A ( LEFT )
              A ( LEFT )   = ATEMP ( 1 )
C
              K            = INDX (LEFT+1)
              INDX (LEFT+1) = INDX (LEFT)
              INDX (LEFT)  = K
          ENDIF
C
          IF ( A(LEFT) .GT. A(RIGHT) ) THEN
C
              ATEMP ( 1 ) = A ( LEFT )
              A ( LEFT )  = A ( RIGHT )
              A ( RIGHT ) = ATEMP ( 1 )
C
              K           = INDX (LEFT)
              INDX (LEFT)  = INDX (RIGHT)
              INDX (RIGHT) = K
C
              IF ( A(LEFT+1) .GT. A(LEFT) ) THEN
C
                  ATEMP ( 1 )   = A ( LEFT+1 )
                  A ( LEFT+1 )  = A ( LEFT )
                  A ( LEFT )    = ATEMP ( 1 )
C
                  K            = INDX (LEFT+1)
                  INDX (LEFT+1) = INDX (LEFT)
                  INDX (LEFT)   = K
              ENDIF
          ENDIF
C
          ATEMP(2) = A (LEFT)
C
C       ... ATEMP(2) IS NOW THE MEDIAN VALUE OF THE THREE A-S.  NOW MOVE
C           FROM THE LEFT AND RIGHT ENDS SIMULTANEOUSLY, EXCHANGING
C           INDX UNTIL ALL A-S LESS THAN  ATEMP(2)  ARE SYMBOLLICALLY
C           PACKED TO THE LEFT, ALL A-S LARGER THAN  ATEMP(2)  ARE
C           PACKED TO THE RIGHT.
C
          I = LEFT+1
          J = RIGHT
C
C         LOOP
C             REPEAT I = I+1 UNTIL A(I) >= ATEMP;
C             REPEAT J = J-1 UNTIL A(J) <= ATEMP;
C         EXIT IF J < I;
C             << EXCHANGE AS I AND J >>
C         END
C
  700     CONTINUE
  800     CONTINUE
              I  = I + 1
              IF ( A(I) .LT. ATEMP(2) )  GO TO 800
C
  900     CONTINUE
              J  = J - 1
              IF ( A(J) .GT. ATEMP(2) )  GO TO 900
C
          IF (J .GE. I)  THEN
C
              ATEMP(1)  = A(I)
              A(I)      = A(J)
              A(J)      = ATEMP(1)
C
              K         = INDX (I)
              INDX (I)   = INDX (J)
              INDX (J)   = K
              GO TO 700
          ENDIF
C
          ATEMP(1)   = A (LEFT)
          A(LEFT)    = A (J)
          A(J)       = ATEMP(1)
C
          K          = INDX (LEFT)
          INDX (LEFT) = INDX (J)
          INDX (J)    = K
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
          IF (  A (I) .LE. A (IP1) )  GO TO 2400
C
C             ... OUT OF ORDER ... MOVE UP TO CORRECT PLACE
C
              ATEMP(1)  = A   (I)
              K         = INDX (I)
              J         = IP1
              JM1       = I
C
C             ... REPEAT ... UNTIL 'CORRECT PLACE FOR ATEMP FOUND'
C
 2200         CONTINUE
                  A   (JM1) = A   (J)
                  INDX (JM1) = INDX (J)
                  JM1       = J
                  J         = J + 1
                  IF ( J .GT. N )  GO TO 2300
                  IF (  A (J) .LT. ATEMP(1) )  GO TO 2200
C
 2300         CONTINUE
              A   (JM1) = ATEMP(1)
              INDX (JM1) = K
C
 2400     CONTINUE
          IP1 = I
          I   = I - 1
          IF ( I .GT. 0 )  GO TO 2100
C
C     -----------------------------------------------------------
C     FILE IS NOW SORTED AND INDX IS THE POINTER.
C     STABILIZE THE SORT IF REQUESTED.
C     -----------------------------------------------------------
C
      IF ( ISTBLE .NE. 0 ) GO TO 4000
C
C ... SEARCH FOR A SEQUENCE OF IDENTICAL ENTRIES
C
      I = 1
C
 3100 CONTINUE
      IP1 = I + 1
      IF ( A(I) .NE. A(IP1) ) THEN
         I = IP1
         IF ( I .GE. N ) GO TO 4000
         GO TO 3100
      ENDIF
C
C ... 2 ENTRIES ARE IDENTICAL.  FIND END OF SEQUENCE
C
      J = IP1 + 1
C
 3300 CONTINUE
      IF ( J .LE. N ) THEN
         IF ( A(I) .NE. A(J) ) GO TO 3400
            J = J + 1
            GO TO 3300
      ENDIF
C
C ... SORT INDX FOR IDENTICAL ENTRIES
C
 3400 CONTINUE
      L = J - I
      CALL HJSORT ( INDX(I), L, IER )
      IF ( IER .NE. 0 ) THEN
C
C ... UNEXPECTED ERROR RETURN FROM HJSORT
C
        IER = -6
        GO TO 8300
      ENDIF
      I = J
      IF ( I .GE. N ) GO TO 4000
      GO TO 3100
C
C--------------------------------------------------------------------
C
C     -----------------------------------------
C     IF REQUESTED, RETURN A IN ORIGINAL ORDER.
C     -----------------------------------------
C
 4000 CONTINUE
      IF ( IOPT .EQ. 0 ) GO TO 9000
C
      CALL HJPRMY ( A, N, INDX, IER )
      IF( IER .NE. 0 ) THEN
C
C ... UNEXPECTED ERROR RETURN FROM HJPRMX
C
        IER = -7
        GO TO 8300
      ENDIF
C
      GO TO 9000
C
C---------------------------------------------------------------------
C
C ... ERROR TRAPS
C
C ... INTERNAL STACK OVERFLOW
C
 8000 CONTINUE
      IER = -5
      INDX(1) = 0
      CALL HHERR ( 3, NAME, IER, 0 )
      GO TO 9000
 8300 CONTINUE
      INDX(1) = 0
      CALL HHERR ( 5, NAME, IER, 0 )
C
C----------------------------------------------------------------------
C
C ... END OF HJSRTN
C
 9000    CONTINUE
      RETURN
      END
