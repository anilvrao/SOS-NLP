      subroutine xdsldi (pvttol, lfront, ncol  , pnlsiz, updsiz,
     1                   pospon, nodoff, zpcntl, npcntl, invp  , 
     2                   frontl, frontu, pvtblk, panell, panelu,
     3                   swap  , temp  , fnzlf , npanel, xpanel, 
     4                   mpanel, qincor, fdiag , lnzloc, iopos1, 
     5                   wafil1, wafil4, ocbufr, lbuffr, fctops, slvops,
     6                   ppfmon, xpboxs, psboxs, inpexp, inpsiz, inzsiz,
     7                   izfail, rtpexp, rtpsiz, rtzsiz, rzfail, error )

c
c  purpose -- top level routine to factor the current front.
c             unsymmetric version, in memory
c
c  created            -- 21-may-98, rgg, derived from xdslei
c  last modifications -- 31-jul-98, rgg, added wafil4 to hold upper
c                                        tri. factor
c                        01-sep-98, rgg, 32 bit integer mods
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
c      pospon -- number of postponed columns in this front
c      nodoff -- number of nodes eliminated so far
c      zpcntl -- array for zero pivot control
c      npcntl -- array for negative pivot control
c      frontl -- the symmetric frontal matrix in packed storage
c                lower tri
c      frontu -- the symmetric frontal matrix in packed storage
c                upper tri
c      npanel -- number of panels eliminated so far
c      qincor -- logical flag on whether l is to be stored in core 
c                or on disk
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
c      npanel -- number of panels
c      xpanel -- pointer to first node of each panel
c      mpanel -- max. size of a panel
c      iopos1 -- current position on wafil1 for storing l on disk
c      wafil1 -- i/o file for storing l
c      wafil1 -- i/o file for storing u
c      ocbufr -- real buffer for i/o when storing l
c      lbuffr -- length of real buffer
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
c                        = -3 i/o error on file wafil4 (upper entries)
c                        = -4 exact zero diagonal (when not pivoting)
c                        = -5 negative diagonal (when not pivoting
c                             and under user prohibition of negative
c                             diagonal)
c                        = -6 failure in panel pivot scheme (when
c                             pivoting)
c                        = -99 unknown error from lower level routine
c
c     ------------------------------------------------------------------
 
c     --------------------
c     ... global variables
c     --------------------
 
      integer            lfront, ncol  , pnlsiz, pospon, 
     1                   nodoff, npanel, mpanel,
     2                   wafil1, wafil4, lbuffr, error , updsiz
 
      integer            invp  (*),      pvtblk (*),     xpanel(*),  
     1                   swap  (*),      iopos1(2)

      logical            qincor, ppfmon
 
      double precision   pvttol, fctops, slvops, fnzlf
 
      double precision   zpcntl(*),      npcntl(*),
     1                   frontl(*),      frontu(*),      
     1                   panell(*),      panelu(*),      
     1                   temp(lfront,*)           ,
     1                   fdiag (*),      lnzloc(*),      ocbufr(*)

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
 
      integer            colbgn, colend, frtbeg, j1    , j2    ,
     1                   kcol  , kloc  , loclfr, lstpnl, ncolf , 
     2                   pnlbgn, pnlcol, pnlrow, updbeg, ipanel, 
     3                   len   , nswap , cmajor, matsiz, skip  ,
     4                   interr

      integer            iopos(2), iopsav(2)

      logical            qzpvt1, qnpvt1, qroot
 
      double precision   fnewnz
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      external            xdsle2, xdsld4, xdsld5, xdsld9, xdsld3
 
c     ------------------------------------------------------------------

      cmajor = 1
      colbgn = 1
      kloc   = 1
      matsiz = lfront * ( lfront + 1 ) / 2

      pnlbgn = npanel + 1

      iopsav(1) = iopos1(1)
      iopsav(2) = iopos1(2)

      qzpvt1 = .true.
      if ( zpcntl(2) .ne. 0. ) qzpvt1 = .false.
      qnpvt1 = .true.
      if ( npcntl(2) .ne. 0. ) qnpvt1 = .false.
      qroot = .false.
      if ( ncol .eq. lfront  ) qroot  = .true.
c.debug
c     write(6,'("in xdsldi - pnlsiz, ncol, pvttol = ", 
c    1            2i8,1pd15.5)') pnlsiz, ncol, pvttol
c     write(6,'("in xdsldi - wafil1, wafil4               = ", 
c    1            3i8        )') wafil1, wafil4              
c     write(6,'("in xdsldi - qzpvt1, qnpvt1, qroot        = ", 
c    1            3l8        )') qzpvt1, qnpvt1, qroot       
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
c     call xdslp5('lower front', lfront*(lfront+1)/2, frontl, 6 )
c     call xdslp5('upper front', lfront*(lfront+1)/2, frontu, 6 )
c     end if
c.debug

c     ----------------------------------------------
c     ... extract panel from in-core data structure.
c         the columns of the front are stored as 
c         columns in the panel.
c     ----------------------------------------------

      call xdsle2 ( cmajor, colbgn, colend, pnlrow, frontl(frtbeg),
     1              pnlcol, pnlrow, panell )

      call xdsle2 ( cmajor, colbgn, colend, pnlrow, frontu(frtbeg),
     1              pnlcol, pnlrow, panelu )

      j2 = 1
      do j1 = 1, pnlcol
          panelu(j2) = 0.0
          j2 = j2 + pnlrow + 1
      enddo

c.debug
c     if ( ncol .eq. lfront ) then
c     write(6,'("after panel extract")')
c     call xdslp5('panell after extract', pnlcol*pnlrow, panell, 6 )
c     call xdslp5('panelu after extract', pnlcol*pnlrow, panelu, 6 )
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
     1                 panelu, pvtblk(colbgn), 
     2                 fctops, slvops, interr )

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
             if       ( interr .eq. -1 )  then
                error = -1
             else if  ( interr .eq. -2 )  then
                error = -6
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
c     if ( ncol .eq. lfront ) then
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
c.debug
c     if ( qswtch ) then
c     call xdslp5('front before update', loclfr*(loclfr+1)/2, 
c    1            frontl(updbeg), 6 )
c     call xdslp5('front before update', loclfr*(loclfr+1)/2, 
c    1            frontu(updbeg), 6 )
c     end if
c.debug

      if ( loclfr .eq. 0 ) go to 150

c     ------------------------------------------------------
c     ... carve up temp into an array of ncolf by updsiz and
c         an array of loclfr by updsiz
c     -------------------------------------------------

      j2 = ( updsiz * ncolf + 1 ) / lfront
      j1 = ( updsiz * ncolf + 1 ) - j2 * lfront

      if ( j2 * lfront .eq. updsiz * ncolf + 1 ) then
         j1 = lfront
      else
         j2 = j2 + 1
      end if
c.debug
c     if ( qswtch ) then
c     write(6,'("xdsld9 - ncolf, loclfr, lfront, j1, j2 = ", 5i8)')
c    1                     ncolf, loclfr, lfront, j1, j2 
c     end if
c.debug

      call xdsld9 ( ncolf , pnlcol, pnlrow, updsiz, panell,
     1              panelu, pvtblk(colbgn), loclfr, 
     2              frontl(updbeg),         frontu(updbeg),
     2              temp, temp(j1,j2), fctops  )

c.debug
c     if ( ncol .eq. lfront ) then
c     write(6,'("after panel update")')
c     call xdslp5('panell after update ', pnlcol*pnlrow, panell, 6 )
c     call xdslp5('panelu after update ', pnlcol*pnlrow, panelu, 6 )
c     call xdslp5('lower front after update', loclfr*(loclfr+1)/2, 
c    1            frontl(updbeg), 6 )
c     call xdslp5('upper front after update', loclfr*(loclfr+1)/2, 
c    1            frontu(updbeg), 6 )
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
c    1            panell, 6 )
c     call xdslp5('postponed columns in panel', pnlcol*pnlrow, 
c    1            panelu, 6 )
c     call xdslp5('full front before stuffing', lfront*(lfront+1)/2, 
c    1            frontl, 6 )
c     call xdslp5('full front before stuffing', lfront*(lfront+1)/2, 
c    1            frontu, 6 )
c.debug

      call xdslek ( cmajor, ncolf, colbgn, colend, pnlrow, 
     1              frontl(frtbeg), pnlcol, pnlrow, panell )

      call xdslek ( cmajor, ncolf, colbgn, colend, pnlrow, 
     1              frontu(frtbeg), pnlcol, pnlrow, panelu)

c.debug
c     call xdslp5('full front after stuffing', lfront*(lfront+1)/2, 
c    1            frontl, 6 )
c     call xdslp5('full front after stuffing', lfront*(lfront+1)/2, 
c    1            frontu, 6 )
c.debug


      end if

c     ------------------------------------------------------------
c     ... copy the panel into final storage for fdiag, lnz and unz
c     ------------------------------------------------------------

      fnewnz = 0.
c.debug
c     if ( qswtch ) then
c     write(6,'("before xdsld3 - fnzlf, iopos1, kloc = ", 3i8)') 
c    1                            fnzlf, iopos1, kloc
c     write(6,'("                qincor              = ", 3l8)') 
c    1                            qincor            
c     end if
c.debug

      call xdsld3 ( ncolf , pnlcol, pnlrow, panell, panelu,
     1              fnewnz, qincor, fdiag(colbgn),  lnzloc(kloc),
     2              iopos1, wafil1, wafil4, ocbufr, lbuffr, interr )

      if ( interr .ne. 0 ) then
         if  ( interr .eq. -2 )  then
            error = -2
         else
     1   if  ( interr .eq. -3 )  then
            error = -3
         else
            error = -99
         end if
         return
      end if
  
c.debug
c     if ( qswtch ) then
c     write(6,'("after xdsld3 - kloc, fnzlf = ", 2i8)') kloc, fnzlf
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
c.debug
c     write(6,'("kcol, swap(kcol)                    = ", 5i8)')     
c    1            kcol, swap(kcol)                    
c.debug
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
c     write(6,'("post processing of the factors      = ")')     
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
          matsiz = pnlrow * pnlcol + skip
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
c             call xdslp5('pre-lnzloc', pnlrow*pnlcol, 
c    1                    lnzloc(kloc), 6 )
c             call xdslp5('pre-lnzloc', pnlrow*pnlcol, 
c    1                    lnzloc(kloc+matsiz), 6 )
c.debug

              call xdslel ( colend, pnlcol, pnlrow, nswap-colend,
     1                      swap(colend+1), lnzloc(kloc) )

              call xdslel ( colend, pnlcol, pnlrow, nswap-colend,
     1                      swap(colend+1), lnzloc(kloc+matsiz) )

c.debug
c             call xdslp5('post-lnzloc',pnlrow*pnlcol, 
c   1                     lnzloc(kloc), 6 )
c             call xdslp5('post-lnzloc',pnlrow*pnlcol, 
c   1                     lnzloc(kloc+matsiz), 6 )
c.debug

              kloc = kloc + matsiz + pnlrow*pnlcol

c         ---------------------------------------
c         ... factorization is being held on disk  
c         ---------------------------------------

          else                         

              call xdslw9 ( iopos, skip )

              len = pnlcol * pnlrow
c.debug
c     write(6,'("reading in lnz")')
c.debug

              call xdslw5 ( wafil1, panell, iopos, len, interr )

              if ( interr .ne. 0 ) then
                 if  ( interr .eq. -1 )  then
                    error = -2
                 else
                    error = -99
                 end if
                 return
              end if

c.debug
c             call xdslp5('pre-lnzloc', pnlrow*pnlcol, panell, 6 )
c.debug
     
              call xdslel ( colend, pnlcol, pnlrow, nswap-colend,
     1                      swap(colend+1), panell )

c.debug
c             call xdslp5('post-lnzloc', pnlrow*pnlcol, panell, 6 )
c.debug

c.debug
c             write(6,'("writing out lnz")')
c.debug

              call xdslw6 ( wafil1, panell, iopos, len, interr )

              if ( interr .ne. 0 ) then
                 if  ( interr .eq. -1 )  then
                    error = -2
                 else
                    error = -99
                 end if
                 return
              end if

c             -------------------------
c             ... now do upper triangle
c             -------------------------

c.debug
c     write(6,'("reading in unz")')
c.debug

              call xdslw5 ( wafil4, panelu, iopos, len, interr )

              if ( interr .ne. 0 ) then
                 if  ( interr .eq. -1 )  then
                    error = -3
                 else
                    error = -99
                 end if
                 return
              end if

c.debug
c             call xdslp5('pre-lnzloc', pnlrow*pnlcol, panelu, 6 )
c.debug
     
              call xdslel ( colend, pnlcol, pnlrow, nswap-colend,
     1                      swap(colend+1), panelu )

c.debug
c             call xdslp5('post-lnzloc', pnlrow*pnlcol, panelu, 6 )
c.debug

c.debug
c             write(6,'("writing out unz")')
c.debug

              call xdslw6 ( wafil4, panelu, iopos, len, interr )

              if ( interr .ne. 0 ) then
                 if  ( interr .eq. -1 )  then
                    error = -3
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
