
      DOUBLE PRECISION FUNCTION hdmcon( i )
      
      IMPLICIT         none
      INTEGER          i, onei, bigi
      DOUBLE PRECISION mcon, small, big, temp, e, pi, eulers, one
      PARAMETER        (
     &   onei   = 1  ,
     &   one    = 1.0D0,
     &   e      = 2.71828182845904523536028747135266249775724709D+0,
     &   pi     = 3.14159265358979323846264338327950288419716940D+0,
     &   eulers = 5.77215664901532860606512090082402431042159336D-1 )
      
        select case (i)
          case (1)                        ! Clobber constant
            mcon = huge(one)
          case (2)                        ! Symmetric range
            small = tiny(one)
            big   = huge(one)
            if (small*big > one) then
              mcon = one/small
            else
              mcon = big
            end if
          case (3)                        ! Overflow threshold
            mcon = huge(one)
          case (4)                        ! Underflow threshold
            mcon = tiny(one)
          case (5)                        ! Relative spacing
            mcon = epsilon(one)
          case (6)                        ! Relative precision
            mcon = epsilon(one) 
          case (7)                        ! Radix
            mcon = radix(one)
          case (8)                        ! Mantissa length
            mcon = digits(one)
          case (9)                        ! Exponential overflow threshold
            mcon = aint(log(huge(one)))
          case (10)                       ! Exponential underflow threshold
            mcon = aint(log(tiny(one)))
          case (11)                       ! Max floatable integer
            bigi = huge(onei) 
            temp = radix(one)/epsilon(one) - one
            if (bigi < temp) then
              mcon = bigi
            else
              mcon = temp
            end if
          case (12)                       ! The constant pi
            mcon = pi
          case (13)                       ! The constant e
            mcon = e
          case (14)                       ! Euler's constant
            mcon = eulers
          case (15)                       ! radians in degree
            mcon = pi/180
          case (16)                       ! degrees in radian
            mcon = 180/pi
          case default                    ! Out of range
            mcon = huge(one)
        end select
      
        hdmcon = mcon
      return
      end
