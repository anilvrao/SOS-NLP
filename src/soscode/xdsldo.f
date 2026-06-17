      subroutine xdsldo (pvttol, lfront, ncol  , pnlsiz, updsiz, matsiz, 
     1                   pospon, nodoff, zpcntl, npcntl, invp  , pvtblk, 
     2                   panell, panelu, swap  , temp  , fnzlf , npanel, 
     3                   xpanel, mpanel, fdiag , wafil5, watrn5, 
     4                   iopos1, wafil1, watrn1, wafil4, watrn4,
     5                   iopos2, wafil2, walen2, watrn2,
     6                   ocbufr, locbfr, fctops, slvops,
     7                   ppfmon, xpboxs, psboxs, inpexp, inpsiz, inzsiz,
     8                   izfail, rtpexp, rtpsiz, rtzsiz, rzfail, error )

c
c  purpose -- top level routine to factor the current front.
c             out-of-core unsymmetric version
c
c  created            -- 25-jul-98, rgg, derived from xdsleo and xdsldi
c  last modifications -- 31-jul-98, rgg, added wafil5 for out-of-core
c                                        processing of front
c                        01-sep-98, rgg, 32 bit integer mods
c                        14-apr-99, rgg, corrected passing of postponed 
c                                        columns between panels 
c                        30-jan-01, jgl  corrected error handling
c                                        for block elimination step
c
c  input variables --
c
c      pvttol -- pivoting tolerance
c      lfront -- size of the frontal matrix
c      ncol   -- number of columns corresponding to nodes in the
c                supernode
c      pnlsiz -- width of the panel used during the factorization
c      updsiz -- width of the panel used during the factorization
c                update
c      matsiz -- size of the lower tri. frontal matrix
c      pospon -- number of postponed columns in this front
c      nodoff -- number of nodes eliminated so far
c      zpcntl -- array for zero pivot controls
c      npcntl -- array for negative pivot controls
c
c  working storage --
c
c      panell -- temporary array of pnlsiz by lfront - lower tri
c      panelu -- temporary array of pnlsiz by lfront - upper tri
c      swap   -- temporary integer array of length lfront
c      temp   -- temporary array of size lfront by updsiz
c
c  input/output variable --
c
c      fnzlf  -- number of nonzeroes stored in lnz so far
c      npanel -- number of panels eliminated so far
c      xpanel -- pointer to first node of each panel
c      mpanel -- max. size of a panel
c      wafil5 -- i/o file for storing current front
c      watrn5 -- count of i/o transfer for i/o file wafil5
c      iopos1 -- current position on wafil1 for storing l on disk
c      wafil1 -- i/o file for storing l
c      watrn1 -- count of i/o transfer for i/o file wafil1
c      wafil4 -- i/o file for storing u
c      watrn4 -- count of i/o transfer for i/o file wafil4
c      iopos2 -- current position on wafil2 for storing update
c                matrices on disk
c      wafil2 -- i/o file for storing update matrices
c                on wafil2
c      walen2 -- length of i/o file wafil2
c      watrn2 -- count of i/o transfer for i/o file wafil2
c      ocbufr -- real buffer for i/o when storing l
c      locbfr -- length of real buffer
c      fctops -- factor operation count
c      slvops -- solve operation count
c
c  output variable --
c
c      ncol   -- the number of columns that were eliminated.  the
c                difference between the input value and the output
c                value of ncol is equal to the number of columns
c                postponed.
c      invp   -- inverse permutation array.
c      pvtblk -- integer vector indicates if diagonal element is a
c                1 x 1 or 2 x 2 block pivot.
c      fdiag  -- local section of diagonal of the factorization
c      error  -- error return,
c                if error =  0, success,
c                         = -1, exactly singular matrix (when pivoting)
c                         = -2, i/o error on wafil1 (factor entries)
c                         = -3, i/o error on wafil2 (stack entries)
c                         = -4, i/o error on wafil4 (upper entries)
c                         = -5, i/o error on wafil5 (single front)
c                         = -6, exact zero diagonal (when not pivoting)
c                         = -7, negative diagonal (when not pivoting
c                               and under user prohibition of negative
c                               diagonal)
c                         = -8, failure in panel pivot scheme (when
c                               pivoting)
c                        = -99, unknown error from lower level routine
c
c     ------------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer            lfront, ncol  , pnlsiz, matsiz, pospon, 
     1                   nodoff, npanel, wafil1, wafil4, wafil5,
     2                   iopos2, wafil2, locbfr, error , updsiz
 
      integer            invp  (*),      pvtblk (*),     xpanel(*),  
     1                   swap  (*),      iopos1(2)

      double precision   pvttol, fctops, slvops, watrn1, walen2, 
     1                   watrn2, watrn4, watrn5, fnzlf
 
      double precision   zpcntl(*),      npcntl(*),
     1                   panell(*),      panelu(*),
     1                   temp(lfront,*),
     1                   fdiag (*),      ocbufr(*)

c     ------------------------------------
c     ... global variables for diagnostics
c     ------------------------------------

      integer           xpboxs, psboxs,  izfail, rzfail,
     1                  inpexp (xpboxs), inpsiz (psboxs),
     2                  inzsiz (psboxs), rtpexp (xpboxs),
     3                  rtpsiz (psboxs), rtzsiz (psboxs)

c     -------------------
c     ... local variables
c     -------------------
 
      integer            cmajor, colbgn, colend, frtoff, 
     1                   iopsvu, j1    , j2    ,
     1                   ipanel, kcol  , loclfr, lstpnl, ncolf , 
     2                   pnlbgn, pnlcol, pnlrow, updoff, interr

      integer            iopos5, iops5f, iops5u, k     , l     ,
     1                   len   , mpanel, ncolpp, updlen, nswap ,
     2                   skip

      integer            iopsav(2), iopos(2), iops1w(2), idummy(1)

      logical            qzpvt1, qnpvt1, qroot, ppfmon

      double precision   ftemp
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      external            xdsle2, xdsld4, xdslea, xdsld5, xdsleb, 
     1                    xdsldh
 
c     ------------------------------------------------------------------

      colbgn = 1
      cmajor = 1
      matsiz = lfront * ( lfront + 1 ) / 2
      iopos5 = 1
      pnlbgn = npanel + 1

      iops1w(1) = iopos1(1)
      iops1w(2) = iopos1(2)
      iopsav(1) = iopos1(1)
      iopsav(2) = iopos1(2)

      qzpvt1 = .true.
      if ( zpcntl(2) .ne. 0. ) qzpvt1 = .false.
      qnpvt1 = .true.
      if ( npcntl(2) .ne. 0. ) qnpvt1 = .false.
      qroot = .false.
      if ( ncol .eq. lfront  ) qroot  = .true.
c.debug
c     write(6,'("in xdsldo - pnlsiz, ncol, pvttol = ", 
c    1            3i8,1pd15.5)') pnlsiz, ncol, pvttol
c.debug

c     -----------------------------------
c     ... start processing the next panel
c     -----------------------------------

  100 continue
      colend = min ( colbgn + pnlsiz - 1, ncol )

      pnlcol = colend - colbgn + 1
      pnlrow = lfront - colbgn + 1
      loclfr = lfront - colend

      npanel = npanel + 1
      xpanel(npanel) = colbgn + nodoff
c.debug
c     if ( colbgn .eq. 1 .and. colend .ne. ncol ) then
c         write(6,'("front split into more than 1 panel - npanel = ",
c    1                i8)') npanel
c     end if
c.debug

      kcol   = colbgn - 1
      frtoff = lfront * kcol - ( kcol*(kcol-1) ) / 2
      iops5f = iopos5 + frtoff

      kcol   = colend
      updoff = lfront * kcol - ( kcol*(kcol-1) ) / 2
      iops5u = iopos5 + updoff
c.debug
c     write(6,'(/,"npanel, colbgn, colend, pnlcol, pnlrow = ", 
c    1            5i8)')
c    1            npanel, colbgn, colend, pnlcol, pnlrow
c     write(6,'("frtoff, updoff, fnzlf, iops5f, iops5u   = ", 
c    1            5i8)')
c    1            frtoff, updoff, fnzlf, iops5f, iops5u
c.debug

c     --------------------------------------------------
c     ... extract panel from out-of-core data structure.
c         the columns of the front are stored as either
c         rows in the panel.
c     --------------------------------------------------

      call xdsleg ( cmajor, colbgn, colend, pnlrow, 
     1              iops5f, wafil5, watrn5, locbfr, ocbufr,
     2              pnlcol, pnlrow, panell, interr )
c.debug
c     write(6,'("after reading in the panell - interr = ",i8)') interr
c     write(6,'("ncol, watrn5 = ", i8,f12.0)') 
c    1            colend-colbgn+1, watrn5
c.debug

      if ( interr .ne. 0 ) then
         if  ( interr .eq. -1 )  then
            error = -5
         else
            error = -99
         end if
         return
      end if

c.debug
c     write(6,'("after reading in the panell ")')
c     if ( pnlrow .lt. 128 ) then
c     call xdslp5('panel after reading in', pnlcol*pnlrow, panell, 6 )
c     end if
c.debug

      call xdsleg ( cmajor, colbgn, colend, pnlrow, 
     1              iops5f+matsiz,  wafil5, watrn5, locbfr, ocbufr,
     2              pnlcol, pnlrow, panelu, interr )
c.debug
c     write(6,'("after reading in the panelu - interr = ",i8)') interr
c     write(6,'("watrn5 = ", f12.0)') watrn5
c.debug

      if ( interr .ne. 0 ) then
         if  ( interr .eq. -1 )  then
            error = -5
         else
            error = -99
         end if
         return
      end if

      j2 = 1
      do j1 = 1, pnlcol
          panelu(j2) = 0.0
          j2 = j2 + pnlrow + 1
      enddo

c.debug
c     write(6,'("after reading in the panelu ")')
c     if ( pnlrow .lt. 128 ) then
c     call xdslp5('panel after reading in', pnlcol*pnlrow, panelu, 6 )
c     end if
c.debug

c     ----------------------------------------------------
c     ... factor the columns in the panel
c         note:  on output, ncolf is the number of columns
c                actually factored
c     ----------------------------------------------------
c.debug
c     if ( qswtch ) then
c     write(6,'(/,"pnlcol, pnlrow         = ", 3i8)')
c    1              pnlcol, pnlrow 
c     write(6,'(  "pvttol                 = ", 1pd15.5)') pvttol
c     end if
c.debug

      if ( pvttol .eq. 0. ) then

          call xdsld4 ( zpcntl, npcntl, pnlcol, pnlrow, panell,
     1                  panelu, pvtblk(colbgn), 
     2                  fctops, slvops, interr )

          if ( interr .ne. 0 )  then
             if  ( interr .eq. -1 )  then
                error = -6
             else
     1       if  ( interr .eq. -2 )  then
                error = -7
             else
                error = -99
             end if
             return
          end if

          ncolf = pnlcol

      else

c.debug
c     if ( qswtch ) then
c     write(6,'("before xdsld5")')
c     end if
c.debug

          call xdsld5 ( pvttol, zpcntl, qroot,
     1                  pnlrow, pnlcol, pospon, colbgn-1,
     2                  invp(colbgn),   panell, panelu, pvtblk(colbgn), 
     3                  temp(1,1), temp(1,2), swap(colbgn),
     4                  ncolf , fctops, slvops,
     5                  ppfmon, xpboxs, psboxs, inpexp,
     6                  inpsiz, inzsiz, izfail, rtpexp,
     7                  rtpsiz, rtzsiz, rzfail, interr )


          if ( interr .ne. 0 )  then
             if  ( interr .eq. -1 )  then
                error = -1
             else
     1       if  ( interr .eq. -2 )  then
                error = -8
             else
                error = -99
             end if
             return
          end if

      end if

      if ( qzpvt1 .and. zpcntl(3) .ne. 0. ) then
          zpcntl(3) = zpcntl(3) + colbgn - 1
          qzpvt1    = .false.
      end if

      if ( qnpvt1 .and. npcntl(3) .ne. 0. ) then
          npcntl(3) = npcntl(3) + colbgn - 1
          qnpvt1    = .false.
      end if
c.debug
c     if ( qswtch ) then
c     write(6,'("after panel factor - interr = ", i8)') interr
c     end if
c.debug

c.debug
c     if ( qswtch ) then
c     write(6,'("after panel factor-ncolf, pnlcol, pnlrow = ", 3i8)') 
c    1                               ncolf, pnlcol, pnlrow
c     call xdslp5('panell after factor ', pnlcol*pnlrow, panell, 6 )
c     call xdslp5('panelu after factor ', pnlcol*pnlrow, panelu, 6 )
c     call xislp3('pvtblk after factor', ncolf, pvtblk(colbgn), 6 )
c     call xislp3('invp after factor'  , ncolf, invp  (colbgn), 6 )
c     call xislp3('swap after factor'  , ncolf, swap  (colbgn), 6 )
c     end if
c.debug

      xpanel(npanel+1) = xpanel(npanel) + ncolf

      if ( ncolf .eq. 0 ) go to 200 

      mpanel = max ( mpanel, pnlrow*ncolf - ( ncolf*(ncolf+1) ) / 2 )

c     ----------------------------------------------
c     ... update the rest of the front with factored 
c         columns from the panel
c     ----------------------------------------------

      if ( loclfr .eq. 0 ) go to 150

c     ------------------------------------------------------
c     ... carve up temp into an array of ncolf by updsiz and
c         an array of loclfr by updsiz
c     ------------------------------------------------------

      j2 = ( updsiz * ncolf + 1 ) / lfront
      j1 = ( updsiz * ncolf + 1 ) - j2 * lfront

      if ( j2 * lfront .eq. updsiz * ncolf + 1 ) then
         j1 = lfront
      else
         j2 = j2 + 1
      end if
c.debug
c     write(6,'("xdsldh - ncolf, loclfr, lfront, j1, j2 = ", 5i8)')
c    1                     ncolf, loclfr, lfront, j1, j2 
c     write(6,'("xdsldh - updsiz                        = ", 5i8)')
c    1                     updsiz                        
c     call xdslp5('panel before update ', pnlcol*pnlrow, panel, 6 )
c.debug

      call xdsldh ( ncolf , pnlcol, pnlrow, updsiz,
     1              panell, panelu, pvtblk(colbgn),
     1              loclfr, locbfr, ocbufr, matsiz, iops5u, wafil5,
     2              watrn5, temp, temp(j1,j2), fctops, interr )
 
      if ( interr .ne. 0 ) then
         if  ( ( interr .eq. -2 ) .or. ( interr .eq. -3 ) )  then
            error = -5
         else
            error = -99
         end if
         return
      end if

c.debug
c     write(6,'("after panel update")')
c     write(6,'("watrn5 = ", f12.0)') watrn5
c     call xdslp5('panel after update ', pnlcol*pnlrow, panel, 6 )
c.debug

c     -----------------------------------------------------------
c     ... put any postponed columns back into the front on wafil5
c     -----------------------------------------------------------

  150 continue
      if ( ncolf .gt. 0 .and. ncolf .ne. pnlcol ) then

          k      = colbgn - 1 + ncolf
          updoff = lfront * k - ( k*(k-1) ) / 2
          iops5u = iopos5 + updoff

c.debug
c     write(6,'("before stuffing postponed columns back into front")')
c     call xdslp5('postponed columns in panell', pnlcol*pnlrow, 
c    1            panell, 6 )
c     call xdslp5('postponed columns in panelu', pnlcol*pnlrow, 
c    1            panelu, 6 )
c.debug

          call xdslem ( cmajor, ncolf ,      1, pnlcol, pnlrow, 
     1                  pnlcol, pnlrow, panell, locbfr, ocbufr,
     2                  iops5u, wafil5, watrn5, interr )

          if ( interr .ne. 0 ) then
             if  ( interr .eq. -2 )  then
                error = -5
             else
                error = -99
             end if
             return
          end if

          call xdslem ( cmajor, ncolf ,      1, pnlcol, pnlrow, 
     1                  pnlcol, pnlrow, panelu, locbfr, ocbufr,
     2                  iops5u+matsiz , wafil5, watrn5, interr )

          if ( interr .ne. 0 ) then
             if  ( interr .eq. -2 )  then
                error = -5
             else
                error = -99
             end if
             return
          end if

      end if

c     -------------------------------------------------------
c     ... copy the panel into final storage for fdiag and lnz
c     -------------------------------------------------------

c.debug
c     write(6,'("before xdsle3")') 
c     write(6,'("cmajor, ncolf , pnlcol, pnlrow, fnzlf = ", 5i8)') 
c    1            cmajor, ncolf , pnlcol, pnlrow, fnzlf  
c     write(6,'("iops1w, wafil1, locbfr                = ", 3i8)') 
c    1            iops1w, wafil1, locbfr    
c.debug

      call xdsld8 ( ncolf , pnlcol, pnlrow, panell, panelu,
     1              fnzlf , fdiag(colbgn), matsiz,
     2              iops1w, wafil1, wafil4, ocbufr, locbfr, interr )

      if ( interr .ne. 0 ) then
         if  ( interr .eq. -2 )  then
            error = -2
         else if ( interr .eq. -4 )  then
            error = -4
         else
            error = -99
         end if
         return
      end if
  
c.debug
c     write(6,'("after xdsld8 - fnzlf, iops1w = ", 3i8)') 
c    1                           fnzlf, iops1w
c     call xdslp5('factor diagonal', ncolf, fdiag (colbgn), 6 )
c     call xislp3('pivot block    ', ncolf, pvtblk(colbgn), 6 )
c.debug

c     ----------------------------
c     ... loop back for next panel
c     ----------------------------

      if ( colend .ne. ncol ) then
          colbgn = colbgn + ncolf 
          go to 100
      end if

c     ------------------------------------------------------------------

c     ----------------------------------------------------------
c     ... factorization complete.  move update matrix to wafil2.
c         preserve column numbers of postponed columns at the 
c         start of the update matrix.  also move the upper tri.
c         of the factorization to its final position.
c     ----------------------------------------------------------

  200 continue
      iopos1(1) = iops1w(1)
      iopos1(2) = iops1w(2)

      ncolf  = colbgn + ncolf - 1
      ncolpp = ncol - ncolf
      ncol   = ncolf

      updoff = lfront * ncolf - ( ncolf*(ncolf-1) ) / 2
      iops5u = iopos5 + updoff
      iopsvu = iops5u + matsiz

c     ----------------------------------------------
c     ... put labels for postponed columns in ocbufr
c     ----------------------------------------------

      call icopy ( ncolpp, invp(ncolf+1), 1, ocbufr, 1 )

      k = ncolpp

c     ------------------------------------------------------
c     ... copy columns of lower tri. update matrix to ocbufr
c         and dump ocbufr
c     ------------------------------------------------------

      l      = lfront - ncolf
      updlen = l * ( l + 1 ) / 2

  210 continue
      len = min ( k + updlen, locbfr )
      if ( len .eq. 0 ) go to 240
c.debug
c     write(6,'("moving update matrix - k, len, iops5u, iopos2 = ",
c    1   4i8)') k, len, iops5u, iopos2
c.debug

      call xdslw1 ( wafil5, 2, idummy, idummy, ocbufr(k+1),
     1              iops5u, len - k, interr )

      if ( interr .ne. 0 ) then
         if  ( interr .eq. -1 )  then
            error = -5
         else
            error = -99
         end if
         return
      end if

      call xdslw2 ( wafil2, 2, idummy, idummy, ocbufr,
     1              iopos2, len, interr )

      if ( interr .ne. 0 ) then
         if  ( interr .eq. -1 )  then
            error = -3
         else
            error = -99
         end if
         return
      end if

      iops5u = iops5u + ( len - k ) 
      iopos2 = iopos2 + len
      updlen = updlen - ( len - k )
c.debug 
c     write(6,'("after moving lower update matrix - watrn5 = ", 
c    1           f12.0)')                            watrn5
c     call xdslp5('lower update matrix', len-k, ocbufr(k+1), 6 )
c.debug

      ftemp  = iopos2 - 1
      walen2 = max ( walen2, ftemp ) 
      watrn2 = watrn2 + len
      watrn5 = watrn5 + ( len - k )

      k      = 0

      go to 210

c     ------------------------------------------------------
c     ... copy columns of upper tri. update matrix to ocbufr
c         and dump ocbufr
c     ------------------------------------------------------

  240 continue
      iops5u = iopsvu 
      l      = lfront - ncolf
      updlen = l * ( l + 1 ) / 2

  250 continue
      len = min ( updlen, locbfr )
      if ( len .eq. 0 ) go to 260
c.debug
c     write(6,'("moving update matrix - len, iops5u, iopos2 = ",
c    1   4i8)') len, iops5u, iopos2
c.debug

      call xdslw1 ( wafil5, 2, idummy, idummy, ocbufr,
     1              iops5u, len, interr )

      if ( interr .ne. 0 ) then
         if  ( interr .eq. -1 )  then
            error = -5
         else
            error = -99
         end if
         return
      end if

      call xdslw2 ( wafil2, 2, idummy, idummy, ocbufr,
     1              iopos2, len, interr )

      if ( interr .ne. 0 ) then
         if  ( interr .eq. -1 )  then
            error = -3
         else
            error = -99
         end if
         return
      end if

      iops5u = iops5u + len
      iopos2 = iopos2 + len
      updlen = updlen - len 
c.debug 
c     write(6,'("after moving upper update matrix - watrn5 = ", 
c    1           f12.0)')                            watrn5
c     call xdslp5('upper update matrix', len, ocbufr, 6 )
c.debug

      ftemp  = iopos2 - 1
      walen2 = max ( walen2, ftemp ) 
      watrn2 = watrn2 + len
      watrn5 = watrn5 + len

      go to 250

c     --------------------------------------------------------
c     ... determine if any postprocessing of the factorization
c         of this front is required
c     --------------------------------------------------------

  260 continue
      if ( npanel - pnlbgn .ge. 1 .and. pvttol .gt. 0. ) then

          colbgn = ncol + 1
c.debug
c     write(6,'("pnlbgn, npanel                      = ", 5i8)')     
c    1            pnlbgn, npanel 
c     write(6,'("colbgn, ncolf, ncol                 = ", 5i8)')     
c    1            colbgn, ncolf, ncol                 
c.debug
          do ipanel = npanel, pnlbgn+1, -1
              colend = colbgn - 1
              colbgn = xpanel(ipanel) - nodoff
c.debug
c     write(6,'("ipanel, colbgn, colend              = ", 5i8)')     
c    1            ipanel, colbgn, colend              
c.debug
              do kcol = colbgn, colend
c.debug
c     write(6,'("kcol, swap(kcol)                    = ", 5i8)')     
c    1            kcol, swap(kcol)                    
c.debug
                  if ( swap(kcol) .ne. kcol ) then
                      lstpnl = ipanel
                      go to 300
                  end if

              enddo

          enddo

      end if

c     -------------------------------
c     ... no post processing required
c     -------------------------------

      return

c     ------------------------------------------------
c     ... post processing of the factorization of this
c         front is required
c     ------------------------------------------------

  300 continue
      iopos(1) = iopsav(1)
      iopos(2) = iopsav(2)
c.debug
c     write(6,'("postprocessing of factor - lstpnl   = ", 5i8)')     
c    1            lstpnl                              
c.debug

      nswap = colend

      do 310 ipanel = pnlbgn, lstpnl - 1

          colbgn = xpanel(ipanel)   - nodoff
          colend = xpanel(ipanel+1) - nodoff - 1
          pnlcol = colend - colbgn + 1
          pnlrow = lfront - colend 

          skip   = pnlcol * ( pnlcol - 1 ) / 2
          call xdslw9 ( iopos, skip )
          len    = pnlcol * pnlrow
c.debug
c     write(6,'(///)')     
c     write(6,'("ipanel, colbgn, colend, nodoff      = ", 5i8)')     
c    1            ipanel, colbgn, colend, nodoff      
c     write(6,'("pnlcol, pnlrow, iopos , len         = ", 5i8)')     
c    1            pnlcol, pnlrow, iopos , len         
c     write(6,'("nswap                               = ", 5i8)')     
c    1            nswap                               
c.debug

c.debug
c     write(6,'("reading in lnz - iopos, len = ",3i8)') 
c    1                             iopos, len
c.debug

          call xdslw5 ( wafil1, panell, iopos , len, interr )

          if ( interr .ne. 0 ) then
             if  ( interr .eq. -1 )  then
                error = -2
             else
                error = -99
             end if
             return
          end if

          call xdslel ( colend, pnlcol, pnlrow, nswap-colend,
     1                  swap(colend+1), panell )

          call xdslw6 ( wafil1, panell, iopos , len, interr )

          if ( interr .ne. 0 ) then
             if  ( interr .eq. -1 )  then
                error = -2
             else
                error = -99
             end if
             return
          end if

c         ----------------------------------
c         ... now do panel of upper triangle
c         ----------------------------------

          call xdslw5 ( wafil4, panelu, iopos , len, interr )

          if ( interr .ne. 0 ) then
             if  ( interr .eq. -1 )  then
                error = -4
             else
                error = -99
             end if
             return
          end if

          call xdslel ( colend, pnlcol, pnlrow, nswap-colend,
     1                  swap(colend+1), panelu )

          call xdslw6 ( wafil4, panelu, iopos , len, interr )

          if ( interr .ne. 0 ) then
             if  ( interr .eq. -1 )  then
                error = -4
             else
                error = -99
             end if
             return
          end if

          call xdslw9 ( iopos, len )
     
  310 continue
      
c     ------------------------------------------------------------------
 
      return
      end 
