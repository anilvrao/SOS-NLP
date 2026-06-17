      subroutine   QPERPT   ( msglvl, iPrint, msg   , inform, prbtyp,
     1                        nZ    , nZr   , maxnZ , nerror, numinf,
     2                        obj   , minsum, errmax )
c     ==================================================================
c     ==================================================================
c     ====  QPERPT / 
c     ====  qperpt -- set error flag and print concluding messages  ====
c     ==================================================================
c     ==================================================================

c     ... parameters
      
      integer           msglvl, iPrint, inform, nZ    , nZr   , maxnZ ,
     1                  nerror, numinf

      logical           minsum
      
      double precision  obj   , errmax
      
      character(len=6)  msg

      character(len=2)  prbtyp
      
c     ... local variables

      logical           prnt
      
c     ===================================================================

c     ... last modification -- 25-march-1996

c     ===================================================================

      prnt = msglvl .gt. 0  .and.  iPrint .gt. 0
      
      if  (      msg .eq. 'feasbl' )  then
         
         inform = 0
         if  ( prnt )  write (iPrint, 2001)
         
      else
     1if  ( msg .eq. 'optiml' )  then
         
         inform = 0
         if  ( prnt )  write (iPrint, 2002) prbtyp
         
      else
     1if  ( msg .eq. 'deadpt'  .or.  msg .eq. 'weak  ' )  then
         
         inform = 1
         if  ( prnt )  then
            if  ( prbtyp .eq. 'QP' )  then
               write (iPrint, 2010)
               if  ( nZ .gt. nZr) write (iPrint, 2015) nZ-nZr
            else
               write (iPrint, 2011)
            end if
         end if
         
      else
     1if  ( msg .eq. 'unbndd' )  then

         inform = 2
         if  ( prnt )  write (iPrint, 2020) prbtyp

      else
     1if  ( msg .eq. 'infeas' )  then

         inform = 3
         if  ( prnt )  write (iPrint, 2030)

      else
     1if  ( msg .eq. 'rowerr' )  then

         inform = 3
         if  ( prnt )  write (iPrint, 2035)

      else
     1if  ( msg .eq. 'itnlim' )  then

         inform = 4
         if  ( prnt )  write (iPrint, 2040)

      else
     1if  ( msg .eq. 'Rz2big' )  then

         inform = 5
         if  ( prnt )  write (iPrint, 2050) maxnZ

      else
     1if  ( msg .eq. 'errors' )  then

         inform = 6
         if  ( prnt )  write (iPrint, 2060) nerror

      else

         inform = 7
         if  ( prnt ) write (iprint, 2070) 

      end if

      if  ( prnt )  then
         
         if  ( inform .lt.   5 )  then
            
            if      (numinf .eq. 0 )  then
               if  ( prbtyp .ne. 'FP') write (iPrint, 3000) prbtyp, obj

            else
     1      if  ( inform .eq. 3 )  then

               if  ( msg .eq. 'infeas' )  then
                  if  ( .not. minsum )  then
                     write (iPrint, 3010) obj
                  else
                     write (iPrint, 3011) obj
                  end if
               else if  ( msg .eq. 'rowerr' )  then
                  write (iPrint, 3015) errmax
               end if
            else
               write (iPrint, 3020) obj
            end if
            
         end if
         
      end if

      return
       

 2001 format(/ ' Exit QPOPT  - Feasible point found.     ')
      
 2002 format(/ ' Exit QPOPT  - Optimal ', a2, ' solution.')
      
 2010 format(/ ' Exit QPOPT  - Iterations terminated at a dead-point',
     1         ' (check the optimality conditions).     ')
      
 2011 format(/ ' Exit QPOPT  - Optimal solution is not unique.' )
      
 2015 format(  '            - Artificial constraints in working set = ',
     1         i4 )
      
 2020 format(/ ' Exit QPOPT  - ', a2, ' solution is unbounded.' )
      
 2030 format(/ ' Exit QPOPT  - No feasible point for the linear',
     1         ' constraints.')
      
 2035 format(/ ' Exit QPOPT  - Cannot satisfy the constraints to the',
     1         ' accuracy requested.')
      
 2040 format(/ ' Exit QPOPT  - Too many iterations.')
      
 2050 format(/ ' Exit QPOPT  - Reduced Hessian exceeds assigned',
     1         ' dimension.   maxnZ = ', i4
     2       / ' Problem abandoned.' )
      
 2060 format(/ ' Exit QPOPT  - ', i10, ' errors found in the input',
     1         ' parameters.  Problem abandoned.' )
      
 2070 format(/ ' Exit QPOPT  - Internal Failure Detected.  ',
     1         ' Problem abandoned.'  )
      
 3000 format(/ ' Final ', a2, ' objective value =', g16.7 )

 3010 format(/ ' Sum of infeasibilities =',         g16.7 )

 3011 format(/ ' Minimum sum of infeasibilities =', g16.7 )

 3015 format(/ ' Maximum row error =',              g16.7 )

 3020 format(/ ' Final sum of infeasibilities =',   g16.7 )

c     end of QPERPT (qperpt)
      end
