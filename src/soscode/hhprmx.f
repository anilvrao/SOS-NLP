      SUBROUTINE HHPRMX ( A, N, INDX, CHWORK, IER )
C
C----------------------------------------------------------------------
C
C ... PURPOSE  HHPRMX APPLIES THE PERMUTATION IN INDX TO THE
C              CHARACTER ARRAY A.  A(INDX(I)) IS MOVED TO A(I).
C
C ... INPUT    N        NUMBER OF RECORDS IN A TO BE INTERCHANGED.
C
C              INDX      PERMUTATION ARRAY OF LENGTH N.
C
C ... INPUT/   A        CHARACTER ARRAY TO BE REORDERED.
C      OUTPUT
C
C ... WORK     CHWORK   TEMPORARY CHARACTER VARIABLE.
C
C ... OUTPUT   IER      SUCCESS/ERROR FLAG.
C
C                       IER = 0     NORMAL RETURN
C                           = -1    N .LE. 0
C                           = -2    INDX(1) .LE. 0
C                           = -3    LEN(CHWORK) .LT. LEN(A(1))
C
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C
C ... GLOBAL VARIABLES
C
         CHARACTER(LEN=*)  A(*), CHWORK
C
         INTEGER           N, INDX(*), IER
C
C----------------------------------------------------------------------
C
C ... LOCAL VARIABLES
C
         CHARACTER(LEN=8)  NAME
C
         INTEGER           I, NEXT, NOW
C
         DATA              NAME / 'HHPRMX' /
C
C----------------------------------------------------------------------
C
C ... CHECK INPUT
C
         IER = 0
            IF ( LEN( CHWORK ) .LT. LEN( A(1) ) ) IER = -3
            IF ( INDX(1) .LE. 0 )                  IER = -2
            IF ( N .LE. 0 )                       IER = -1
            IF ( IER .NE. 0 ) THEN
C
C ...... INPUT ERROR DETECTED.
C
               CALL HHERR ( 1, NAME, IER, 0 )
               GO TO 900
            ENDIF
C
C ... APPLY THE PERMUTATION
C
         IF ( N .EQ. 1 ) GO TO 900
C
C ...... SEARCH FOR THE FIRST ENTRY NOT PERMUTED WHICH IS INDICATED BY
C        A NONNEGATIVE VALUE IN INDX.
C
            DO I = 1, N
               IF ( INDX(I) .LE. 0 ) CYCLE
C
C ...... INITIALIZE TO FOLLOW THE CURRENT CHAIN OF PERMUTATIONS
C
               NOW         =  I
               NEXT        =  INDX(NOW)
               INDX(NOW)    = -NEXT
               IF ( NOW .EQ. NEXT ) CYCLE
               CHWORK      =  A(NOW)
C
C ...... FOLLOW THE CHAIN - PERMUTE AS YOU GO UNTIL THE CHAIN ENDS
C
  200          CONTINUE
               IF ( INDX(NEXT) .GT. 0 ) THEN
C
                  A(NOW)      =  A(NEXT)
                  NOW         =  NEXT
                  NEXT        =  INDX(NOW)
                  INDX(NOW)    = -NEXT
                  GO TO 200
               ENDIF
C
C ...... END OF THE CHAIN
C
               A(NOW) = CHWORK
C
C ... END OF SEARCH LOOP FOR THE NEXT CHAIN
C
            ENDDO
C
C----------------------------------------------------------------------
C
C ... PERMUTATION NOW FINISHED.  RESTORE THE INDX ARRAY.
C
            DO I = 1, N
               INDX(I) = - INDX(I)
            ENDDO
C
C----------------------------------------------------------------------
C
C ... END OF HHPRMX
C
  900    CONTINUE
      RETURN
      END
