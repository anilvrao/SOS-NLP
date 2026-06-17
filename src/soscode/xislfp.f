      subroutine   xislfp   ( xpboxs, psboxs, inpexp, inpsiz, inzsiz, 
     1                        izfail, rtpexp, rtpsiz, rtzsiz, rzfail,
     2                        output )

 
c
c     purpose -- print summary statistics on "panel pivot failures",
c                which occur when the LDL^T factorization is unable
c                to find even a single pivot from a particular panel
c                within a front
c
c     created -- 04-0ct-2001  jgl & dkw
c
c     variables  
c
c         inpexp -- count of the number of times a pivot failure
c                   would not have occurred, had the pivot tolerance
c                   be set to 10**(-e) with decimal exponent  e
c                   equal to or smaller than the corresponding index in 
c                   inpexp
c         inpsiz -- count of the number of times a pivot failure that
c                   could have been avoided for some nonzero pivot
c                   tolerance, for a panel whose size is K, where
c                   ceiling (log_2 (K)) is the corresponding index in
c                   inpsiz
c         inzsiz -- as inpsiz, but a count of cases in which the
c                   failure is due to an identically zero diagonal
c                   block
c         izfail -- count of the total number of times the pivot 
c                   failure is due to an identically zero diagonal 
c                   block
c         rtpexp -- count of the number of times a pivot failure
c                   would not have occurred, had the pivot tolerance
c                   be set to 10**(-e) with decimal exponent  e
c                   equal to or smaller than the corresponding index in 
c                   rtpexp
c         rtpsiz -- count of the number of times a pivot failure that
c                   could have been avoided for some nonzero pivot
c                   tolerance, for a panel whose size is K, where
c                   ceiling (log_2 (K)) is the corresponding index in
c                   rtpsiz
c         rtzsiz -- as rtpsiz, but a count of cases in which the
c                   failure is due to an identically zero diagonal
c                   block
c         rzfail -- count of the total number of times the pivot 
c                   failure is due to an identically zero diagonal 
c                   block

      integer           xpboxs, psboxs,  izfail, rzfail, 
     1                  inpexp (xpboxs), inpsiz (psboxs), 
     2                  inzsiz (psboxs), rtpexp (xpboxs), 
     3                  rtpsiz (psboxs), rtzsiz (psboxs),
     4                  output

      integer   i, infail, rtfail

c     ------------------------------------------------------------------
      
      infail = izfail
      rtfail = rzfail
      do i = 1, xpboxs
         infail = infail + inpexp (i)
         rtfail = rtfail + rtpexp (i)
      end do

      if  ( max (infail, rtfail) .gt. 0 )  then

         write ( output, * ) ' Factorization encountered pivot ', 
     1                  'failures for at least one panel'
         write ( output, * )

         if  ( infail .gt. 0 )  then

c           ... at least one interior panel failure

            write ( output, * ) ' number of interior panels for which'
            write ( output, * ) '                 no pivot would work:', 
     1                   izfail      
            write ( output, * ) 'a smaller pivot tolerance would work:', 
     2                  infail - izfail 
            write ( output, * )
     
            write ( output, * ) ' interior pivot failure panel size',
     1                           ' summary'
            write ( output, * ) '     panel size counts'
            write ( output, * ) '     size   no pivot    smaller pivot'
            write ( output, * ) '           acceptable     acceptable'
            do i = 1, psboxs
               if ( max (inpsiz (i), inzsiz (i)) .gt. 0 )  then
                  write ( output, '(i9, i10, i14)') i, inzsiz(i), 
     1                                         inpsiz(i)
               end if
            end do

            write ( output, * )

            if  ( infail - izfail .gt. 0 )  then
               write ( output, * ) ' acceptable pivot exponents'
               write ( output, '(16i5)' ) (-i, i = 1, xpboxs), 
     1                               (inpexp (i), i = 1, xpboxs)
               write ( output, * )
            end if

         end if

         if  ( rtfail .gt. 0 )  then

c           ... at least one root panel failure

            write ( output, * ) ' number of root panels for which:'
            write ( output, * )  '                no pivot would work:', 
     1                   rzfail      
            write ( output, * ) 'a smaller pivot tolerance would work:', 
     2                  rtfail - rzfail 
            write ( output, * )
     
            write ( output, * ) ' root pivot failure panel size ',
     1                           'summary'
            write ( output, * ) '     panel size counts'
            write ( output, * ) '     size   no pivot    smaller pivot'
            write ( output, * ) '           acceptable     acceptable'
            do i = 1, psboxs
               if ( max (rtpsiz (i), rtzsiz (i)) .gt. 0 )  then
                  write ( output, '(i9, i10, i14)') i, rtzsiz(i), 
     1                                         rtpsiz(i)
               end if
            end do

            write ( output, * )

            if  ( rtfail - rzfail .gt. 0 )  then
               write ( output, * ) ' acceptable pivot exponents'
               write ( output, '(16i5)' ) (-i, i = 1, xpboxs), 
     1                               (rtpexp (i), i = 1, xpboxs)
               write ( output, * )
            end if

         end if

      else

         write ( output, * ) 'factorization completed with no pivot',
     1                        ' failures'
         write ( output, * )

      end if

      end
