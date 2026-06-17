      subroutine   QPPRND   ( prbtyp, header, Rset,
     1                        msglvl, iter,
     2                        isdel, jdel, jadd,
     3                        n, nclin, nactiv,
     4                        nfree, nZ, nZr,
     5                        ldR, ldT, istate,
     6                        iPrint, iSumm , lines1, lines2,
     7                        alfa, condRz, condT,
     8                        Dzz, gZrnrm,
     9                        numinf, suminf, notOpt, objqp, trusml,
     A                        Ax, R, T, x, work )

c     ==================================================================
c     ==================================================================
c     ====  QPPRND /                                                ====
c     ====  qpprnd -- qp specialized output routine                 ====
c     ==================================================================
c     ==================================================================

      integer            msglvl, iter  , isdel , jdel  , jadd  ,
     1                   n     , nclin , nactiv, nfree , nZ    ,
     2                   nZr   , ldR   , ldT   , iPrint, iSumm ,
     3                   lines1, lines2, numinf, notOpt

      logical            header, Rset

      character(len=2)   prbtyp

      double precision   alfa  , condRz, condT , Dzz   , gZrnrm, suminf,
     1                   objqp , trusml

      integer            istate (*)

      double precision   Ax (*), R (ldR,*), T (ldT,*), x (n), work (n)

c     ==================================================================
c     QPPRND / qpprnd prints various levels of output for
c     qpcore (QPCORE).
c
c     derived from qpopt version 1.0
c     last modification -- 25-March-1996
c
c           Original version of qpprnd written by PEG, 31-October-1984.
c          This version of  qpprnd  dated  23-Dec-92.
c
c           msg        cumulative result
c           ---        -----------------
c
c       .le.  0        no output.
c
c       .eq.  1        nothing now (but full output later).
c
c       .eq.  5        one terse line of output.
c
c       .ge. 10        same as 5 (but full output later).
c
c       .ge. 20        constraint status,  x  and  Ax.
c
c       .ge. 30        diagonals of  T  and  R.
c
c     ==================================================================

      integer            Itn   , j     , k     , kadd  , kdel  , ndf

      logical            newSet, prtHdr

      double precision   obj

      character(len=2)   ladd, ldel

      character(len=9)   lcondR, lDzz

      character(len=15)  lmchar

      integer            mLine1, mLine2

      double precision   zero, one

      parameter        ( mLine1 = 40 ,
     1                   mLine2 = 5  ,
     2                   zero   = 0.0d0,
     3                   one    = 1.0d0 )

      character(len=2)   lstate(0:5)

      data               lstate(0), lstate(1), lstate(2)
     1                  /'  '     , 'L '     , 'U '     /
      data               lstate(3), lstate(4), lstate(5)
     1                  /'E '     , 'F '     , 'A '     /

c     ==================================================================

      if  ( msglvl .ge. 15  .and.  iPrint .gt. 0)  then
         write (iPrint, 1000) iter, prbtyp
      endif

      if  ( msglvl .ge. 5 )  then

c        ---------------------------------------------------------------
c        Some printing required.  Set up information for the terse line.
c        ---------------------------------------------------------------

         Itn = mod( iter, 1000  )
         ndf = mod( nZr , 10000 )

         if  ( jdel .ne. 0 )  then

            if  ( notOpt .gt. 0 )  then
               write ( lmchar, '( i5, 1p,e10.2 )' ) notOpt, trusml
            else
               write ( lmchar, '( 5x, 1p,e10.2 )' )         trusml
            end if

            if  ( jdel .gt. 0 )  then
               kdel = isdel
            else
     1      if  ( jdel .lt. 0 )  then
               jdel = nZ - nZr + 1
               kdel = 5
            end if

         else

            jdel = 0
            kdel = 0
            lmchar = '               '

         end if

         lDzz   = '         '
         lcondR = '         '
         if  ( Rset  .and.  nZr .gt. 0 )  then
            write( lcondR, '( 1p,e9.1 )' ) condRz
            if  ( Dzz .ne. one)  then
               write( lDzz, '( 1p,e9.1 )' ) Dzz
            endif
         end if

         if  ( jadd .gt. 0 )  then
            kadd = istate(jadd)
         else
            kadd = 0
         end if

         ldel   = lstate (kdel)
         ladd   = lstate (kadd)

         if  ( numinf .gt. 0 )  then
            obj = suminf
         else
            obj = objqp
         end if

c        -----------------------------------
c        If necessary, print a header. 
c        Print a single line of information.
c        -----------------------------------

         if  ( iPrint .gt. 0 )  then

c           ------------------------------
c           Terse line for the Print file.
c           ------------------------------

            newSet = lines1 .ge. mLine1
            prtHdr = msglvl .ge. 15  .or.  header
     1                               .or.  newSet

            if  ( prtHdr )  then
               if  ( prbtyp .eq. 'QP' )  then
                  write (iPrint, 1300)
               else
                  write (iPrint, 1200)
               end if
               lines1 = 0
            end if

            write (iPrint, 1700) Itn, jdel, ldel, jadd, ladd,
     1                          alfa, numinf, obj,
     2                          gZrnrm, ndf, nZ-nZr, 
     3                          n-nfree, nactiv, lmchar, condT,
     4                          lcondR, lDzz
            lines1 = lines1 + 1

         end if

         if  ( iSumm .gt. 0 )  then

c           --------------------------------
c           Terse line for the Summary file.
c           --------------------------------

            newSet = lines2 .ge. mLine2
            prtHdr =                      header 
     1                              .or.  newSet
            if  ( prtHdr )  then
               write (iSumm , 1100)
               lines2 = 0
            end if
            
            write (iSumm , 1700) Itn, jdel, ldel, jadd, ladd,
     1                          alfa, numinf, obj,
     2                          gZrnrm, ndf, nZ-nZr
            lines2 = lines2 + 1

         end if

         if  ( msglvl .ge. 20  .and.  iPrint .gt. 0 )  then

            write (iPrint, 2000) prbtyp
            write (iPrint, 2100) (x(j), istate(j),  j = 1, n)

            if  ( nclin .gt. 0)  then
               write (iPrint, 2200) (Ax(k), istate(n+k), k = 1, nclin)
            endif

            if  ( msglvl .ge. 30 )  then

c              ----------------------------------
c              Print the diagonals of  T  and  R.
c              ----------------------------------

               if  ( nactiv .gt. 0 )  then
                  call dcopy ( nactiv, T(1,nZ+1), ldT+1, work, 1 )
                  write (iPrint, 3000) prbtyp, (work(j), j=1,nactiv)
               end if

               if  ( Rset  .and.  nZr .gt. 0)  then
                  write (iPrint, 3100) prbtyp, (R(j,j) , j=1,nZr )
               endif
               
            end if

            write (iPrint, 5000)
            
         end if

      end if

      header = .false.
      jdel  = 0
      jadd  = 0
      alfa  = zero

      return

 1000 format (/// ' ', ' iteration', i5, ' :: ', a2 
     1          / ' =====================' )

 1100 format (// ' Itn Jdel  Jadd     Step Ninf  Sinf/Objective',
     1           ' Norm gZ   Zr  Art' )

 1200 format (// ' Itn Jdel  Jadd     Step Ninf  Sinf/Objective',
     1           ' Norm gZ   Zr  Art  Bnd  Lin NOpt    Min Lm  Cond T' )

 1300 format (// ' Itn Jdel  Jadd     Step Ninf  Sinf/Objective',
     1           ' Norm gZ   Zr  Art  Bnd  Lin NOpt    Min Lm  Cond T',
     2           '  Cond Rz     Rzz' )

 1700 format (    i4, i5, a1, i5, a1, 1p, e8.1, i5, e16.8, 
     1            e8.1, 2i5, 2i5, a15, e8.0, 2a9 )

 2000 format (/ ' Values and status of the ', a2, ' constraints'
     1        / ' ---------------------------------------' )

 2100 format (/ ' Variables...'                 /  (1x,5(1p,e15.6, i5)))

 2200 format (/ ' General linear constraints...'/  (1x,5(1p,e15.6, i5)))

 3000 format (/ ' Diagonals of ' , a2,' working set factor T'
     1        / (1p, 5e15.6))

 3100 format (/ ' Diagonals of ' , a2, ' triangle Rz        '
     1        / (1p, 5e15.6))

 5000 format (/// ' ---------------------------------------------------'
     1            ,'--------------------------------------------' )

c     end of QPPRND (qpprnd)
      
      end
