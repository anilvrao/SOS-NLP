      SUBROUTINE DSCTR ( M, X, INDX, Y )
C
C--------------------------------------------------------------------
C
C     DSCTR SCATTERS THE ENTRIES OF A PACKED FORM OF A VECTOR STORED
C     IN X AND INDX AND STORES THEM INTO A SCATTERED FORM IN THE
C     VECTOR Y.
C
C
C--------------------------------------------------------------------
C
         INTEGER     M, INDX(*), I, J
C
         DOUBLE PRECISION        Y(*), X(*)
C
C--------------------------------------------------------------------
C
         DO I = 1, M
            J    = INDX(I)
            Y(J) = X(I)
         ENDDO
C
      RETURN
      END
