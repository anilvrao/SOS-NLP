      subroutine xdsls4 ( unsym , wafil1, wafil4, neqns , nrhs  , ndimb,
     1                    rhs   , nsuper, xsup  , xpanel, xlindx,
     2                    lindxg, diag,   pvtblk, slvbsz, temp1 , 
     3                    lwork , work  , ierr    )
 
c
c     created            -- 22-jun-87, rgg
c     last modification  -- 24-feb-89, rgg
c                           09-nov-95, rgg, added temp1 and temp2
c                                           to reflect mods to xdsls5
c                                           and xdsls6
c                           27-feb-97, rgg, converted to panels
c                                           switched to dgemv and dgemm
c                           10-mar-97, rgg, reconstructed i/o for
c                                           out-of-core solves
c                           27-mar-98, rgg, converted to processing 
c                                           slvbsz rhs at one time
c                           30-jul-98, rgg, added wafil4 for unsym.
c                                           problems
c                           01-sep-98, rgg, 32 bit integer mods
c
c     purpose
c     -------
c
c         xdsls4 performs the numeric solution of ax = b after
c         the matrix a has been factored with out-of-core multi-
c         frontal factorization.
c
c     method
c     ------
c
c         the matrix a has been decomposed into l * d * l' and
c         written to unit wafil1.  xdsls4 reads in the factorization
c         at the panel level and performs the standard forward
c         elimination, diagonal scaling and back substitution.
c         for unsym. problems, the upper tri. of the factorization
c         is on wafil4.
c
c     input variables
c     ---------------
c
c         unsym   l   unsymemtric flag
c         wafil1  i   word addressable logical unit containing the
c                     lower tri. of the factorization of a.
c         wafil4  i   word addressable logical unit containing the
c                     upper tri. of the factorization of a (unsym).
c         neqns   i   number of equations.
c         nrhs    i   number of right hand sides
c         ndimb   i   leading dimension of b.
c         nsuper  i   number of super nodes.
c         xsup    i   array giving the starting panel number
c                     for each super node.
c         xpanel  i   pointers into the panel partition
c         xlindx  i   array giving the starting pointer for the row
c                     indices in each column of l.
c         lindxg  i   array containing the global row indices for
c                     each column of l.
c         diag    d   array containing the diagonal of factorization
c         pvtblk  i   pivot block information
c         slvbsz  i   number of rhs vectors to process at one time
c         lwork   i   length of work
c
c     input/output variables
c     ----------------------
c
c         rhs     d   two dimensional array containing the right
c                     hand sides to be solved.
c
c     working storage
c     ---------------
c
c         temp1   d   work vector of length slvbsz*mxlfrt
c         work    d   temporary work storage of length lwork
c
c     output variables
c     ----------------
c
c         ierr    i   error flag.
c                     =  0    successful completion
c                     = -1    i/o error encountered on wafil1
c                     = -2    i/o error encountered on wafil4
c
c     subprograms used
c     ----------------
c
c         xdsls5  perform the forward elimination step for a given
c                 super node.
c         xdsls6  perform the back substitution step for a given
c                 super node.
c
c     ================================================================
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      logical             unsym
 
      integer             wafil1, wafil4, neqns , nrhs  , ndimb ,
     1                    nsuper, slvbsz, lwork , ierr
 
      integer             xsup(*),        xpanel(*),      xlindx(*),
     1                    lindxg(*),      pvtblk(*)
 
      double precision    rhs(ndimb,*),   diag(*),        work(*),
     1                    temp1(*)
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------
 
      integer             indbgn, istart, isuper,
     1                    l,      ldtemp, nodext, 
     2                    ipanel, ipbgn , ipend , n1    , n2

      integer             jpanel, jpnext, jsuper,
     1                    len, l1, l2, l3, lstpnl, npanel, pnlptr

      integer             ipos(2)
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      external            xdslw5,  xdsls6, xdsls5
 
c     ================================================================
 
c     -----------------------------------------------------
c     ... perform the forward elimination step and diagonal
c         scaling for each supernode
c     -----------------------------------------------------
 
      ierr   = 0
      lstpnl = 0
      npanel = xsup(nsuper+1) - 1

      ipos(1) = 1
      ipos(2) = 1
c.debug
c     write(6,'("entry to xdsls4 - nsuper, npanel, lwork = ",
c    1    4i8)')                    nsuper, npanel, lwork
c.debug
 
      do 100 isuper = 1, nsuper

          ipbgn  = xsup(isuper)
          ipend  = xsup(isuper+1) - 1
          n2     = xpanel(ipend+1) - xpanel(ipbgn)
          
          indbgn = xlindx ( isuper )
          nodext = xlindx ( isuper+1 ) - indbgn

          do ipanel = ipbgn, ipend

              if ( ipanel .gt. lstpnl ) then

c                 ---------------------------------------------
c                 ... decide how many panels to read in
c                     note:  in unsymmetric case the panels for
c                            l and u are intertwined.  they can 
c                            only be read one at a time.
c                 ---------------------------------------------

                  len    = 0
                  jsuper = isuper
                  jpnext = ipend + 1
                  l3     = nodext

                  do jpanel = ipanel, npanel
     
                      if ( jpanel .eq. jpnext ) then
                          jsuper = jsuper + 1
                          jpnext = xsup(jsuper+1)
                          l3     = xlindx(jsuper+1) - xlindx(jsuper)
                      end if
c.debug
c     write(6,'("in do 10 - jpanel, jsuper, jpnext, l3     = ", 4i8)')
c    1                       jpanel, jsuper, jpnext, l3     
c.debug

                      l1 = xpanel(jpanel+1) - xpanel(jpanel  )
                      l2 = xpanel(jpnext  ) - xpanel(jpanel+1)
                      l  = l1*(l1+l2+l3) - (l1*(l1+1))/2
c.debug
c     write(6,'("in do 10 - l1    , l2    , l     , len    = ", 4i8)')
c    1                       l1    , l2    , l     , len    
c     write(6,'("in do 10 - lwork                          = ", 4i8)')
c    1                       lwork                          
c.debug

                      if ( l + len .gt. lwork ) exit

                      len    = len + l
                      lstpnl = jpanel

                      if ( unsym ) exit

                  enddo

c                 ----------------------------------
c                 ... read in the next set of panels
c                 ----------------------------------

c.debug
c     write(6,'("after 20 - lstpnl, len   , ipos           = ", 4i8)')
c    1                       lstpnl, len   , ipos           
c.debug

                  call xdslw5 ( wafil1, work, ipos, len, ierr )
                  if ( ierr .ne. 0 ) go to 8000
 
                  call xdslw9 ( ipos, len )

                  pnlptr = 1

              end if

c             ------------------------------------------------
c             ... perform triangluar solve with diagonal block
c                 of super node, update rhs with just computed
c                 values, and scale with diagonal.
c             ------------------------------------------------
 
              istart = xpanel ( ipanel )
              n1     = xpanel ( ipanel+1 ) - istart

              n2     = n2 - n1
c.debug
c     write(6,'("before xdsls5 - n1, n2, istart, nodext    = ", 4i8)')
c    1                            n1, n2, istart, nodext
c     write(6,'("before xdsls5 - pnlptr                    = ", 4i8)')
c    1                            pnlptr                 
c.debug

              ldtemp = max ( 1, n2 + nodext )
 
              call xdsls5 ( istart-1, n1, n2, nodext, ndimb, nrhs,
     1                      pvtblk, diag, work(pnlptr), rhs,
     2                      lindxg(indbgn), ldtemp, slvbsz, temp1 )
c.debug
c     write(6,'("in xdsls4 after xdsls5 - ipanel = ",i8)') ipanel
c.debug

              l      = n1*(n1+n2+nodext) - (n1*(n1+1))/2
              pnlptr = pnlptr + l
 
c             --------------------------------------------------
c             ... end of forward elimination loop for this panel
c             --------------------------------------------------
 
          enddo
 
c         ------------------------------------------------------
c         ... end of forward elimination loop for this supernode
c         ------------------------------------------------------
 
  100 continue
 
c     ================================================================
 
c     ---------------------------------------------------------
c     ... perform the back substitution step for each supernode
c     ---------------------------------------------------------

      lstpnl = npanel + 1
 
      do 300 isuper = nsuper, 1, -1
 
          ipbgn  = xsup(isuper)
          ipend  = xsup(isuper+1) - 1
          n2     = 0
          
          indbgn = xlindx ( isuper )
          nodext = xlindx ( isuper+1 ) - indbgn

          do ipanel = ipend, ipbgn, -1

              if ( ipanel .lt. lstpnl ) then

c                 ---------------------------------------------
c                 ... decide how many panels to read in
c                     note:  in unsymmetric case the panels for
c                            l and u are intertwined.  they can 
c                            only be read one at a time.
c                 ---------------------------------------------

                  len    = 0
                  jsuper = isuper
                  jpnext = ipbgn - 1
                  l3     = nodext

                  do jpanel = ipanel, 1, -1
     
                      if ( jpanel .eq. jpnext ) then
                          jsuper = jsuper - 1
                          jpnext = xsup(jsuper) - 1
                          l3     = xlindx(jsuper+1) - xlindx(jsuper)
                      end if
c.debug
c     write(6,'("in do210 - jpanel, jsuper, jpnext, l3     = ", 4i8)')
c    1                       jpanel, jsuper, jpnext, l3     
c.debug

                      l1 = xpanel(jpanel+1)       - xpanel(jpanel  )
                      l2 = xpanel(xsup(jsuper+1)) - xpanel(jpanel+1)
                      l  = l1*(l1+l2+l3) - (l1*(l1+1))/2
c.debug
c     write(6,'("in do210 - l1    , l2    , l     , len    = ", 4i8)')
c    1                       l1    , l2    , l     , len    
c     write(6,'("in do210 - lwork                          = ", 4i8)')
c    1                       lwork                          
c.debug

                      if ( l + len .gt. lwork ) exit

                      len    = len + l
                      lstpnl = jpanel

                  enddo

c                 ----------------------------------
c                 ... read in the next set of panels
c                 ----------------------------------

                  call xdslw9 ( ipos, -len )
c.debug
c     write(6,'("after220 - lstpnl, len   , ipos           = ", 4i8)')
c    1                       lstpnl, len   , ipos           
c.debug

                  if ( unsym ) then 
                      call xdslw5 ( wafil4, work, ipos, len, ierr )
                      if ( ierr .ne. 0 ) go to 8100
                  else
                      call xdslw5 ( wafil1, work, ipos, len, ierr )
                      if ( ierr .ne. 0 ) go to 8000
                  end if

                  pnlptr = len + 1

              end if
 
c             ------------------------------------------------------
c             ... update rhs vectors with previously computed values
c                 and perform back substitution with diagonal block.
c             ------------------------------------------------------
 
              istart = xpanel ( ipanel )
              n1     = xpanel ( ipanel+1 ) - istart

              l      = n1 * ( n1 + n2 + nodext ) 
     1               - ( n1 * ( n1 + 1 ) ) / 2
              pnlptr = pnlptr - l
c.debug
c     write(6,'("before xdsls6 - n1, n2, istart, nodext    = ", 4i8)')
c    1                            n1, n2, istart, nodext
c     write(6,'("before xdsls6 - pnlptr                    = ", 4i8)')
c    1                            pnlptr                 
c     write(6,'("pre-s6 - isuper, ipanel, n1, n2, nodext, pnlptr = ", 
c    1          6i8)')     isuper, ipanel, n1, n2, nodext, pnlptr 
c     call xdslp5('start of i/o buffer', 5, work, 6 )
c.debug

              ldtemp = max ( 1, n2 + nodext )

              call xdsls6 ( istart-1, n1, n2, nodext, ndimb, nrhs,
     1                      pvtblk, work(pnlptr), rhs, lindxg(indbgn),
     2                      ldtemp, slvbsz, temp1 )
c.debug
c     write(6,'("in xdsls4 after xdsls6 - ipanel = ",i8)') ipanel
c.debug
 
c             ------------------------------------------------
c             ... end of back substitution loop for this panel
c             ------------------------------------------------

              n2 = n2 + n1
 
      enddo
 
c         ----------------------------------------------------
c         ... end of back substitution loop for this supernode
c         ----------------------------------------------------
 
  300 continue
 
      go to 9000
 
c     ================================================================
 
c     -----------------------------
c     ... i/o error trap for wafil1
c     -----------------------------
 
 8000 continue
      ierr = -1
      go to 9000
 
c     -----------------------------
c     ... i/o error trap for wafil4
c     -----------------------------
 
 8100 continue
      ierr = -2
      go to 9000
 
c     ================================================================
 
c     ------------------------
c     ... end of module xdsls4
c     ------------------------
 
 9000 continue
      return
      end
