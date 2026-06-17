
      SUBROUTINE ROWCHK(IROW,JCOLST,IVSTAT,IWORK,NDIM,MCON,NZERO,
     $    IPC,IER)
C
C ======================================================================
C     ROWCHK===>rowchk   J.T. BETTS
C ======================================================================
C
C     ===================================================================
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C         PURPOSE:  CHECK FOR A ZERO ROW IN THE JACOBIAN
C
      DIMENSION IROW(NZERO),JCOLST(NDIM+1),IVSTAT(MCON+1+NDIM),
     $    IWORK(MCON)
C
      IER = 0
C
C         COUNT THE NUMBER OF NONZEROS IN ROW I AND STORE
C         THE RESULT IN IWORK(I)
C
      IWORK(1:MCON) = 0
C
      DO J = 1,NDIM
        IV = MCON + 1 + J
        IF(IVSTAT(IV).NE.3) THEN
          DO I = JCOLST(J),JCOLST(J+1)-1
            IR = IROW(I)
            IWORK(IR) = IWORK(IR) + 1
          enddo
        ENDIF
      enddo
C
C         CHECK THAT EACH ROW HAS AT LEAST ONE ELEMENT
C
      DO I = 1,MCON
        IF(IWORK(I).LE.0) THEN
          IER = -1016
          IF (IPC.GE.100)
     &      WRITE(*,*) 'ROW  =',I,'  NO. NONZEROS =',IWORK(I)
        ENDIF
      enddo
C     ----------
C     ... RETURN
C     ----------
C
      RETURN
      END
