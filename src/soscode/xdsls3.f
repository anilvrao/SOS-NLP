      subroutine xdsls3 (unsym , neqns , nsuper, nrhs  , ldrhs ,
     1                   xsup  , xpanel, xlindx, lindxg, pvtblk, diag  ,
     2                   lnz   , rhs   , slvbsz, temp1 )

c
c  purpose -- to solve a supernodal general sparse factorization
c             for multiple right hand sides
c
c  created            -- 07-jul-87, cca
c  last modifications -- 07-jul-87, cca
c                        09-nov-95, rgg, added temp1 and temp2
c                                        to reflect mods to xdsls5
c                                        and xdsls6
c                        27-feb-97, rgg, converted to panels
c                                        switched to dgemv and dgemm
c                        27-mar-98, rgg, converted to processing 
c                                        slvbsz rhs at one time
c                        01-sep-98, rgg, 32 bit integer mods
c
c  input variables --
c
c      unsym  -- unsymmetric flag
c      neqns  -- number of equations
c      nsuper -- number of supernodes
c      nrhs   -- number of right hand sides
c      ldrhs  -- leading dimension of the right hand side array
c      xsup   -- pointers into the supernode partition
c      xpanel -- pointers into the panel partition
c      xlindx -- pointers into the supernode indices
c      lindxg -- global supernodal indices
c      pvtblk -- size of the diagonal blocks in the factorization
c      diag   -- diagonal entries
c      lnz    -- lower triangular entries
c      rhs    -- right hand sides
c      slvbsz -- number of rhs vectors to process at one time
c      temp1  -- work vector of length slvbsz*mxlfrt 
c
c  output variables --
c
c      rhs    -- solutions
c
c  subprograms called --
c
c      none
c
c     ------------------------------------------------------------------
 
c     -----------------------
c     ...  passed parameters
c     -----------------------
 
      logical            unsym
 
      integer            neqns , nsuper, nrhs  , ldrhs ,
     1                   slvbsz

      integer            xsup(*),   xpanel(*), xlindx(*), 
     1                   lindxg(*), pvtblk (*)
 
      double precision   diag(*),   lnz(*),    rhs(ldrhs,*), 
     1                   temp1(*)
 
c     ----------------------
c     ...  local parameters
c     ----------------------
 
      integer            ierr,   indbgn, istart, isuper, kmtx,
     1                   l,      ldtemp, nodext, 
     2                   ipanel, ipbgn , ipend , n1    , n2
 
c     ------------------------------------------------------------------
 
c     -----------------------------------------------------
c     ... perform the forward elimination step and diagonal
c         scaling for each supernode
c     -----------------------------------------------------
 
      kmtx = 1
      ierr = 0
 
      do 100 isuper = 1, nsuper

          ipbgn  = xsup(isuper)
          ipend  = xsup(isuper+1) - 1
          n2     = xpanel(ipend+1) - xpanel(ipbgn)
          
          indbgn = xlindx ( isuper )
          nodext = xlindx ( isuper+1 ) - indbgn
c.debug
c     write(6,'("ipbgn, ipend, n2, indbgn, nodext = ", 5i8)')
c    1            ipbgn, ipend, n2, indbgn, nodext
c.debug

          do ipanel = ipbgn, ipend
 
              istart = xpanel ( ipanel )
              n1     = xpanel ( ipanel+1 ) - istart

              n2     = n2 - n1
c.debug
c     write(6,'("ipanel, istart, n1, n2           = ", 5i8)')
c    1            ipanel, istart, n1, n2          
c     call xdslp5 ( 'rhs before xdsls5', neqns, rhs, 6 )
c.debug
 
c             ------------------------------------------------
c             ... perform triangluar solve with diagonal block
c                 of super node, update rhs with just computed
c                 values, and scale with diagonal.
c             ------------------------------------------------

              ldtemp = max ( 1, n2 + nodext )

              call xdsls5 ( istart-1, n1, n2, nodext, ldrhs, nrhs,
     1                      pvtblk, diag, lnz(kmtx), rhs,
     2                      lindxg(indbgn), ldtemp, slvbsz, temp1 )
c.debug
c     write(6,'("in xdsls3 after xdsls5 - ipanel = ",i8)') ipanel
c     call xdslp5 ( 'rhs after  xdsls5', neqns, rhs, 6 )
c.debug

c             -------------------------------------------------
c             ... adjust pointer to the start of the next front
c             -------------------------------------------------
 
              l    = n1 * ( n1 + n2 + nodext ) 
     1             - ( n1 * ( n1 + 1 ) ) / 2
              kmtx = kmtx + l
    
              if ( unsym ) kmtx = kmtx + l
c.debug
c     write(6,'("l     , kmtx                     = ", 5i8)')
c    1            l     , kmtx                    
c.debug
 
c             -----------------------------------------------------
c             ... end of forward elimination loop for this panel
c             -----------------------------------------------------
 
          enddo
 
c         ------------------------------------------------------
c         ... end of forward elimination loop for this supernode
c         ------------------------------------------------------
 
  100 continue
 
c     ================================================================
 
c     ---------------------------------------------------------
c     ... perform the back substitution step for each supernode
c     ---------------------------------------------------------
 
      do 300 isuper = nsuper, 1, -1

          ipbgn  = xsup(isuper)
          ipend  = xsup(isuper+1) - 1
          n2     = 0
          
          indbgn = xlindx ( isuper )
          nodext = xlindx ( isuper+1 ) - indbgn
c.debug
c     write(6,'("ipbgn, ipend, n2, indbgn, nodext = ", 5i8)')
c    1            ipbgn, ipend, n2, indbgn, nodext
c.debug

          do ipanel = ipend, ipbgn, -1
 
              istart = xpanel ( ipanel )
              n1     = xpanel ( ipanel+1 ) - istart
c.debug
c     write(6,'("ipanel, istart, n1, n2           = ", 5i8)')
c    1            ipanel, istart, n1, n2          
c.debug
 
c             ----------------------------------------------------
c             ... adjust pointer to the start of the current front
c             ----------------------------------------------------
 
              l      = n1 * ( n1 + n2 + nodext ) 
     1               - ( n1 * ( n1 + 1 ) ) / 2
              kmtx   = kmtx - l
c.debug
c     write(6,'("l     , kmtx                     = ", 5i8)')
c    1            l     , kmtx                    
c.debug
 
c             ------------------------------------------------------
c             ... update rhs vectors with previously computed values
c                 and perform back substitution with diagonal block.
c             ------------------------------------------------------

              ldtemp = max ( 1, n2 + nodext )
 
              call xdsls6 ( istart-1, n1, n2, nodext, ldrhs, nrhs,
     1                      pvtblk, lnz(kmtx), rhs, lindxg(indbgn),
     2                      ldtemp, slvbsz, temp1 )
c.debug
c     write(6,'("in xdsls3 after xdsls6 - ipanel = ",i8)') ipanel
c     call xdslp5 ( 'rhs after  xdsls6', neqns, rhs, 6 )
c.debug

              if ( unsym ) kmtx = kmtx - l
    
              n2 = n2 + n1
 
c             ------------------------------------------------
c             ... end of back substitution loop for this panel
c             ------------------------------------------------
 
          enddo
 
c         ----------------------------------------------------
c         ... end of back substitution loop for this supernode
c         ----------------------------------------------------
 
  300 continue
 
c     ------------------------------------------------------------------
 
      return
      end
