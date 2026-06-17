      SUBROUTINE DGTHR ( M, Y, X, INDX )
C
C--------------------------------------------------------------------
C
C     DGTHR GATHERS THE ENTRIES FROM A SPARSE VECTOR REPRESENTED BY
C     Y AND STORES THEM INTO A PACKED FORM IN THE VECTOR Y.  THE
C     INDICIES GATHERED ARE GIVEN BY THE ARRAY INDX.
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
            X(I) = Y(J)
         ENDDO
C
      RETURN
      END
