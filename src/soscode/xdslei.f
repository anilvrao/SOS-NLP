      subroutine xdslei (pvttol, cmajor, lfront, ncol  , pnlsiz, updsiz,
     1                   pospon, nodoff, zpcntl, npcntl, invp  , front ,
     2                   pvtblk, panel , swap  , temp  , fnzlf , npanel,
     3                   xpanel, mpanel, qincor, fdiag , lnzloc, iopos1,
     4                   wafil1, ocbufr, lbuffr, inrtia, fctops, slvops,
     5                   ppfmon, xpboxs, psboxs, inpexp, inpsiz, inzsiz,
     6                   izfail, rtpexp, rtpsiz, rtzsiz, rzfail, llnzl ,
     7                   error )

 
c
c  purpose -- top level routine to factor the current front.
c
c  created            -- 07-feb-97, rgg
c  last modifications -- 28-oct-97, rgg, added zero and negative pivot
c                                        controls
c                         9-mar-98, rgg, added updsiz
c                         1-sep-98, rgg, 32 bit integer mods
c                        30-jan-01, jgl  corrected error handling
c                                        for block elimination step
c
c  input variables --
c
c      pvttol -- pivoting tolerance
c      cmajor -- column major flag
c                .eq. 0 row    major form for panel is used
c                .ne. 0 column major form for panel is used
c      lfront -- size of the frontal matrix
c      ncol   -- number of columns corresponding to nodes in the
c                supernode
c      pnlsiz -- width of the panel used during the factorization
c      updsiz -- width of the panel used during the factorization
c                update
c      pospon -- number of postponed columns in this front
c      nodoff -- number of nodes eliminated so far
c      zpcntl -- array for zero pivot control
c      npcntl -- array for negative pivot control
c      front  -- the symmetric frontal matrix in packed storage
c      npanel -- number of panels eliminated so far
c      qincor -- logical flag on whether l is to be stored in core 
c                or on disk
c
c  working storage --
c
c      panel  -- temporary array of pnlsiz by lfront
c      swap   -- temporary integer array of length lfront
c      temp   -- temporary array of size lfront by updsiz
c
c  input/output variable --
c
c      fnzlf  -- number of nonzeroes stored in lnz so far
c      npanel -- number of panels
c      xpanel -- pointer to first node of each panel
c      mpanel -- max. size of a panel
c      iopos1 -- current position on wafil1 for storing l on disk
c      wafil1 -- i/o file for storing l
c      ocbufr -- real buffer for i/o when storing l
c      lbuffr -- length of real buffer
c      inrtia -- matrix inertia
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
c      lnzloc -- local section of in-core storage for the factorization
c      error  -- error return,
c                if error = 0, success,
c                        = -1 exactly singular matrix (when pivoting)
c                        = -2 i/o error on file wafil1 (factor entries)
c                        = -3 failure in panel pivot scheme (when 
c                             pivoting)
c                        = -4 exact zero diagonal (when not pivoting)
c                        = -5 negative diagonal (when not pivoting
c                             and under user prohibition of negative
c                             diagonal)
c                        = -6 not enough memory in lnzloc
c                        = -99 unknown error from lower level routine
c
c     ------------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer            cmajor, lfront, ncol  , pnlsiz, pospon, 
     1                   nodoff, npanel, mpanel, llnzl ,
     2                   wafil1, lbuffr, error , updsiz
 
      integer            invp  (*),      pvtblk (*),     xpanel(*),  
     1                   swap  (*),      inrtia (3),     iopos1(2)

      logical            qincor, ppfmon
 
      double precision   pvttol, fctops, slvops, fnzlf
 
      double precision   zpcntl(*),      npcntl(*),
     1                   front (*),      panel (*),      
     1                   temp(lfront,*)      ,
     1                   fdiag (*),      lnzloc(llnzl),  ocbufr(*)

c.debug
c     logical qswtch
c     common /rggdbg/ qswtch
c.debug

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
 
      integer            colbgn, colend, frtbeg, interr,  j1    , 
     1                   j2    , kcol  , kloc  , loclfr, lstpnl, 
     2                   ncolf , pnlbgn, pnlcol, pnlrow, updbeg, 
     3                   len   , ipanel, nswap , skip  , need

      integer            iopsav(2), iopos(2)

      logical            qzpvt1, qnpvt1, qroot

      double precision   fnewnz
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      external            xdsle2, xdsle4, xdslea, xdsle5, xdsleb, 
     1                    xdsle9
 
c     ------------------------------------------------------------------

      colbgn = 1
      kloc   = 1
      pnlbgn = npanel + 1

      iopsav(1) = iopos1(1)
      iopsav(2) = iopos1(2)

      qzpvt1 = .true.
      if ( zpcntl(2) .ne. 0. ) qzpvt1 = .false.
      qnpvt1 = .true.
      if ( npcntl(2) .ne. 0. ) qnpvt1 = .false.
      qroot  = .false.
      if ( ncol .eq. lfront )  qroot  = .true.
c.debug
c     if ( qswtch ) then
c     write(6,'("in xdslei - cmajor, pnlsiz, ncol, pvttol = ", 
c    1            3i8,1pd15.5)') cmajor, pnlsiz, ncol, pvttol
c     write(6,'("in xdslei - wafil1                       = ", 
c    1            3i8,1pd15.5)') wafil1                      
c     end if
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
c         write(6,'("colbgn, colend, ncol = ", 3i8)')
c    1                colbgn, colend, ncol 
c     end if
c.debug

      kcol   = colbgn - 1
      frtbeg = lfront * kcol - ( kcol*(kcol-1) ) / 2 + 1

      kcol   = colend
      updbeg = lfront * kcol - ( kcol*(kcol-1) ) / 2 + 1
c.debug
c     if ( qswtch ) then
c     write(6,'(/,"npanel, colbgn, colend, pnlcol, pnlrow = ", 
c    1            5i8)')
c    1            npanel, colbgn, colend, pnlcol, pnlrow
c     write(6,'("frtbeg, updbeg, fnzlf                  = ", 
c    1            5i8)')
c    1            frtbeg, updbeg, fnzlf
c     call xdslp5('full front', lfront*(lfront+1)/2, front, 6 )
c     end if
c.debug

c     ----------------------------------------------
c     ... extract panel from in-core data structure.
c         the columns of the front are stored as 
c         either rows or columns in the panel.
c     ----------------------------------------------

      call xdsle2 ( cmajor, colbgn, colend, pnlrow, front(frtbeg),
     1              pnlcol, pnlrow, panel )
c.debug
c     if ( qswtch ) then
c     write(6,'("after panel extract")')
c     call xdslp5('panel after extract', pnlcol*pnlrow, panel, 6 )
c     end if
c.debug

c     ----------------------------------------------------
c     ... factor the columns in the panel
c         note:  on output, ncolf is the number of columns
c                actually factored
c     ----------------------------------------------------
c.debug
c     if ( qswtch ) then
c     write(6,'(/,"cmajor, pnlcol, pnlrow = ", 3i8)')
c    1              cmajor, pnlcol, pnlrow 
c     write(6,'(  "pvttol                 = ", 1pd15.5)') pvttol
c     end if
c.debug

      if ( pvttol .eq. 0. ) then

         if ( cmajor .eq. 0 ) then

            call xdslea ( zpcntl, npcntl, pnlcol, pnlrow, panel ,
     1                    pvtblk(colbgn), inrtia, fctops, slvops, 
     2                    interr )

         else

           call xdsle4 ( zpcntl, npcntl, pnlcol, pnlrow, panel ,
     1                   pvtblk(colbgn), inrtia, fctops, slvops, 
     2                   interr )
 
         end if

         if ( interr .ne. 0 )  then
            if       ( interr .eq. -1 )  then
               error = -4
            else if  ( interr .eq. -2 )  then
               error = -5
            else
               error = -99
            end if
            return
         end if

         ncolf = pnlcol

      else
         
         if ( cmajor .eq. 0 ) then
c.debug
c     if ( qswtch ) then
c     write(6,'("before xdsleb")')
c     end if
c.debug

              call xdsleb ( pvttol, zpcntl, qroot,
     1                      pnlrow, pnlcol, pospon, colbgn-1,
     1                      invp(colbgn), panel , pvtblk(colbgn), 
     2                      temp(1,1), temp(1,2), swap(colbgn),
     3                      ncolf , inrtia, fctops, slvops, interr)

          else
c.debug
c     if ( qswtch ) then
c     write(6,'("before xdsle5")')
c     end if
c.debug

              call xdsle5 ( pvttol, zpcntl, qroot,
     1                      pnlrow, pnlcol, pospon, colbgn-1,
     2                      invp(colbgn),   panel , pvtblk(colbgn), 
     3                      temp(1,1), temp(1,2), swap(colbgn),
     4                      ncolf , inrtia, fctops, slvops,
     5                      ppfmon, xpboxs, psboxs, inpexp,
     6                      inpsiz, inzsiz, izfail, rtpexp,
     7                      rtpsiz, rtzsiz, rzfail, interr )

          end if

c.debug
c$$$          if ( ncolf .eq. 0 )  then
c$$$             if  (qroot)  then
c$$$                write (*,*) 'root front panel pivot failure -- xdslei'
c$$$             else
c$$$                write (*,*) 'interior front panel pivot failure', 
c$$$     1                      ' -- xdslei'
c$$$             end if
c$$$             write (*,*) ' front begins at column:', nodoff
c$$$             write (*,*) '  panel range:', colbgn, colend
c$$$             if  ( colend .lt. ncol ) then
c$$$                write (*,*) '  front column dim.: ', ncol
c$$$             end if
c$$$          end if
c.debug

          if ( interr .ne. 0 )  then
             if       ( interr .eq. -1 )  then
                error = -1
             else if  ( interr .eq. -2 )  then
                error = -3
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
c     write(6,'("after panel factor - ncolf = ", i8)') ncolf
c     call xdslp5('panel after factor ', pnlcol*pnlrow, panel, 6 )
c     call xislp3('pvtblk after factor', pnlcol, pvtblk(colbgn), 6 )
c     call xislp3('invp after factor'  , pnlcol, invp  (colbgn), 6 )
c     call xislp3('swap after factor'  , pnlcol, swap  (colbgn), 6 )
c     end if
c.debug

      xpanel(npanel+1) = xpanel(npanel) + ncolf

      if ( ncolf .eq. 0 )  then
         go to 200 
      end if

      mpanel = max ( mpanel, pnlrow*ncolf - ( ncolf*(ncolf+1) ) / 2 )

c     ----------------------------------------------
c     ... update the rest of the front with factored 
c         columns from the panel
c     ----------------------------------------------
c.debug
c     if ( qswtch ) then
c     call xdslp5('front before update', loclfr*(loclfr+1)/2, 
c    1            front(updbeg), 6 )
c     end if
c.debug

      if ( loclfr .eq. 0 ) go to 150

      if ( cmajor .eq. 0 ) then

          call xdslef ( ncolf, pnlcol, pnlrow, panel , 
     1                  pvtblk(colbgn), loclfr, front(updbeg), 
     2                  temp, fctops  )

      else

c         ------------------------------------------------------
c         ... carve up temp into an array of ncolf by updsiz and
c             an array of loclfr by updsiz
c         -------------------------------------------------

          j2 = ( updsiz * ncolf + 1 ) / lfront
          j1 = ( updsiz * ncolf + 1 ) - j2 * lfront

          if ( j2 * lfront .eq. updsiz * ncolf + 1 ) then
             j1 = lfront
          else
             j2 = j2 + 1
          end if
c.debug
c     if ( qswtch ) then
c     write(6,'("xdsle9 - ncolf, loclfr, lfront, j1, j2 = ", 5i8)')
c    1                     ncolf, loclfr, lfront, j1, j2 
c     end if
c.debug

          call xdsle9 ( ncolf, pnlcol, pnlrow, updsiz, panel ,
     1                  pvtblk(colbgn), loclfr, front(updbeg), 
     2                  temp, temp(j1,j2), fctops  )

 
      end if

c.debug
c     if ( qswtch ) then
c     write(6,'("after panel update")')
c     call xdslp5('panel after update ', pnlcol*pnlrow, panel, 6 )
c     call xdslp5('front after update', loclfr*(loclfr+1)/2, 
c    1            front(updbeg), 6 )
c     end if
c.debug

c     -------------------------------------------------
c     ... put any postponed columns back into the front
c     -------------------------------------------------

  150 continue
      if ( ncolf .gt. 0 .and. ncolf .ne. pnlcol ) then

c.debug
c     write(6,'("before stuffing postponed columns back into front")')
c     call xdslp5('postponed columns in panel', pnlcol*pnlrow, 
c    1            panel, 6 )
c.debug

          call xdslek ( cmajor, ncolf, colbgn, colend, pnlrow, 
     1                  front(frtbeg), pnlcol, pnlrow, panel )

c.debug
c     call xdslp5('full front after stuffing', lfront*(lfront+1)/2, 
c    1            front, 6 )
c.debug


      end if

c     -------------------------------------------------------
c     ... copy the panel into final storage for fdiag and lnz
c     -------------------------------------------------------

      fnewnz = 0.
c.debug
c     if ( qswtch ) then
c     write(6,'("before xdsle3 - nzlf, iopos1 = ", 2i8)') 
c    1                            nzlf, iopos1
c     end if
c.debug

      if (qincor) then
        need = (ncolf*(ncolf-1))/2 + (pnlrow - ncolf)*ncolf
        if (kloc + need .gt. llnzl) then
          error = -6
          return
        endif
      endif
      call xdsle3 ( cmajor, ncolf , pnlcol, pnlrow, panel , 
     1              fnewnz, qincor, fdiag(colbgn),  lnzloc(kloc),
     2              iopos1, wafil1, ocbufr, lbuffr, interr  )

      if ( interr .ne. 0 ) then
         if  ( interr .eq. -2 )  then
            error = -2
         else
            error = -99
         end if
         return
      end if
  
c.debug
c     if ( qswtch ) then
c     write(6,'("after xdsle3 - kloc, fnewnz = ", 2i8)') 
c    1                           kloc, fnewnz
c     call xdslp5('factor diagonal', ncolf, fdiag(colbgn), 6 )
c     if ( qincor ) then
c         k = fnewnz
c         call xdslp5('lnzloc', k, lnzloc(kloc), 6 )
c     end if
c     end if
c.debug

      kloc  = kloc  + fnewnz
      fnzlf = fnzlf + fnewnz

c     ----------------------------
c     ... loop back for next panel
c     ----------------------------

      if ( colend .ne. ncol ) then
          colbgn = colbgn + ncolf 
          go to 100
      end if
      
c     ------------------------------------------------------------------

c     --------------------------
c     ... factorization complete
c     --------------------------

  200 continue
      ncol = colbgn + ncolf - 1

c     --------------------------------------------------------
c     ... determine if any postprocessing of the factorization
c         of this front is required
c     --------------------------------------------------------
c.debug
c     write(6,'("pnlbgn, npanel                      = ", 5i8)')     
c    1            pnlbgn, npanel 
c.debug

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
      kloc     = 1
      iopos(1) = iopsav(1)
      iopos(2) = iopsav(2)
c.debug
c     write(6,'("pnlbgn, lstpnl                      = ", 5i8)')     
c    1            pnlbgn, lstpnl                      
c.debug

      nswap = colend

      do 310 ipanel = pnlbgn, lstpnl - 1

          colbgn = xpanel(ipanel)   - nodoff
          colend = xpanel(ipanel+1) - nodoff - 1
          pnlcol = colend - colbgn + 1
          pnlrow = lfront - colend 
          skip   = pnlcol * ( pnlcol - 1 ) / 2
c.debug
c     write(6,'(///)')     
c     write(6,'("ipanel, colbgn, colend, nodoff      = ", 5i8)')     
c    1            ipanel, colbgn, colend, nodoff      
c     write(6,'("pnlcol, pnlrow, kloc  , iopsav      = ", 5i8)')     
c    1            pnlcol, pnlrow, kloc  , iopsav      
c     write(6,'("nswap                               = ", 5i8)')     
c    1            nswap                               
c     write(6,'("qincor                              = ", 5l8)')     
c    1            qincor                              
c.debug
          if ( qincor ) then

c         -----------------------------------------
c         ... factorization is being held in memory
c         -----------------------------------------

              kloc = kloc + skip

c.debug
c             call xdslp5('pre-lnzloc', pnlrow*pnlcol, lnzloc(kloc), 6 )
c.debug

              call xdslel ( colend, pnlcol, pnlrow, nswap-colend,
     1                      swap(colend+1), lnzloc(kloc) )

c.debug
c             call xdslp5('post-lnzloc',pnlrow*pnlcol, lnzloc(kloc), 6 )
c.debug

              kloc = kloc + pnlrow * pnlcol

c         ---------------------------------------
c         ... factorization is being held on disk  
c         ---------------------------------------

          else                         

              call xdslw9 ( iopos, skip )

              len = pnlcol * pnlrow
c.debug
c     write(6,'("reading in lnz")')
c.debug

              call xdslw5 ( wafil1, panel, iopos, len, interr )

              if ( interr .ne. 0 ) then
                 if  ( interr .eq. -1 )  then
                    error = -2
                 else
                    error = -99
                 end if
                 return
              end if

c.debug
c             call xdslp5('pre-lnzloc', pnlrow*pnlcol, panel, 6 )
c.debug
     
              call xdslel ( colend, pnlcol, pnlrow, nswap-colend,
     1                      swap(colend+1), panel )

c.debug
c             call xdslp5('post-lnzloc', pnlrow*pnlcol, panel, 6 )
c.debug

c.debug
c             write(6,'("writing out lnz")')
c.debug

              call xdslw6 ( wafil1, panel, iopos, len, interr )

              if ( interr .ne. 0 ) then
                 if  ( interr .eq. -1 )  then
                    error = -2
                 else
                    error = -99
                 end if
                 return
              end if

              call xdslw9 ( iopos, len )

          end if
     
  310 continue
      
c     ------------------------------------------------------------------
 
      return
      end 
