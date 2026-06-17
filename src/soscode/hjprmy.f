      SUBROUTINE HJPRMY ( A, N, INDX, IER )
C
C----------------------------------------------------------------------
C
C ... PURPOSE  HJPRMY APPLIES THE INVERSE OF THE PERMUTATION IN INDX
C              TO THE INTEGER ARRAY A.  A(I) IS MOVED TO A(INDX(I)).
C
C ... INPUT    N        NUMBER OF RECORDS IN A TO BE INTERCHANGED.
C
C              INDX      PERMUTATION ARRAY OF LENGTH N.
C
C ... INPUT/   A        INTEGER ARRAY TO BE REORDERED.
C      OUTPUT
C
C ... OUTPUT   IER      SUCCESS/ERROR FLAG.
C
C                       IER = 0     NORMAL RETURN
C                           = -1    N .LE. 0
C                           = -2    INDX(1) .LE. 0
C
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C
C ... GLOBAL VARIABLES
C
         INTEGER           A(*)
C
         INTEGER           N, INDX(*), IER
C
C----------------------------------------------------------------------
C
C ... LOCAL VARIABLES
C
         CHARACTER(LEN=8)  NAME
C
         INTEGER           ATEMP1, ATEMP2
C
         INTEGER           I, NEXT, NOW
C
         DATA              NAME / 'HJPRMY' /
C
C----------------------------------------------------------------------
C
C ... CHECK INPUT
C
         IER = 0
            IF ( INDX(1) .LE. 0 )                     IER = -2
            IF ( N .LE. 0 )                          IER = -1
            IF ( IER .NE. 0 ) THEN
C
C ...... INPUT ERROR DETECTED.
C
               CALL HHERR ( 1, NAME, IER, 0 )
               GO TO 900
            ENDIF
C
C ... APPLY THE INVERSE OF THE PERMUTATION
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
               NEXT        = INDX(I)
               ATEMP1      = A(I)
C
C ...... FOLLOW THE CHAIN - PERMUTE AS YOU GO UNTIL THE CHAIN ENDS
C
  200          CONTINUE
               IF ( INDX(NEXT) .GT. 0 ) THEN
C
                  ATEMP2      =  ATEMP1
                  ATEMP1      =  A(NEXT)
                  A(NEXT)     =  ATEMP2
C
                  NOW         =  NEXT
                  NEXT        =  INDX(NOW)
                  INDX(NOW)    = -NEXT
                  GO TO 200
               ENDIF
C
C ...... END OF THE CHAIN AND END OF SEARCH LOOP FOR THE NEXT CHAIN
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
C ... END OF HJPRMY
C
  900    CONTINUE
         RETURN
      END
