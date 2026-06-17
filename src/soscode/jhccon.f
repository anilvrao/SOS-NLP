
      INTEGER FUNCTION jhccon( i )
      
      IMPLICIT         none
      INTEGER          i, mcon, idum
      REAL             rdum
      DOUBLE PRECISION ddum
      
        select case (i)
          case (1)                   ! no longer used
            mcon = huge(1)
          case (2)                   ! REAL RECL factor
            inquire(IOLENGTH=mcon) rdum
          case (3)                   ! max record size
            mcon = huge(1)
          case (4)                   ! INTEGER RECL factor
            inquire(IOLENGTH=mcon) idum
          case (5)                   ! DOUBLE PRECISION RECL factor
            inquire(IOLENGTH=mcon) ddum
          case default               ! Out of range
            mcon = huge(1)
        end select
      
        jhccon = mcon
      end
