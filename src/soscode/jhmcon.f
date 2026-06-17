

      INTEGER FUNCTION jhmcon( i )
      
      IMPLICIT none
      INTEGER  i, mcon
      
        select case (i)
          case (1)                   ! Clobber constant
            mcon = huge(1)
          case (2)                   ! Symmetric range
            mcon = huge(1)
          case (3)                   ! Overflow threshold
            mcon = huge(1)
          case (4)                   ! Standard error
            mcon = 6
          case (5)                   ! Standard in
            mcon = 5
          case (6)                   ! Standard out
            mcon = 6
          case (7)                   ! Standard in length
            mcon = 80
          case (8)                   ! Standard out length
            mcon = 132
          case (9)                   ! Decimals in integer
            mcon = int( log10(real(huge(1))) )
          case (10)                  ! Decimals in real
            mcon = precision(1.0)+1
          case (11)                  ! Decimals in double precision
            mcon = precision(1.0d0)+1
          case (12)                  ! Page size in ints
            mcon = 1
          case (13)                  ! Characters in int
            mcon = nint( log(real(huge(1))) / log(256.0) )
            if ( (mcon /= 4) .and. (mcon /= 8) ) mcon = -99999
          case (14)                  ! Ints in real
            mcon = nint( real(digits(1.0)) / real(digits(1)) )
            if ( mcon == 3 ) then
               mcon = 4
            elseif ( (mcon < 1) .or. (mcon > 4) ) then
               mcon = -99999
            endif
          case (15)                  ! Ints in double
            mcon = nint( real(digits(1.0d0)) / real(digits(1)) )
            if ( mcon == 3 ) then
               mcon = 4
            elseif ( (mcon < 1) .or. (mcon > 4) ) then
               mcon = -99999
            endif
          case default               ! Out of range
            mcon = huge(1)
        end select
      
        jhmcon = mcon
      return
      end
