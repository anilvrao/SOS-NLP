      subroutine xdslc2( unsym , neqns , nsuper, xsup  ,
     1                   xpanel, xlindx, lindxg, pvtblk, diag  , 
     2                   lnz   , tnorma, z,      temp1,  cndnum   )
 
c
c  purpose -- to compute the condition number of the factorization
c             performed by xdslf2.
c
c  created            -- 07-jul-89, cca
c  last modifications -- 31-aug-98, rgg, removed use of nzl 
c
c  input variables --
c
c      unsym  -- unsymmetric flag
c      neqns  -- number of equations
c      nsuper -- number of supernodes
c      xsup   -- pointers into the supernode partition
c      xpanel -- pointer into panel partition
c      xlindx -- pointers into the supernode indices
c      lindxg -- global supernodal indices
c      pvtblk -- size of the diagonal blocks in the factorization
c      diag   -- diagonal entries
c      lnz    -- lower triangular entries
c      tnorma -- 1 norm of orginial matrix
c
c  working storage
c
c      z      -- work array of length neqns which holds the
c                approximate null vector.
c      temp1  -- work array of length mxlfrt for use by 
c                solution subroutines
c
c  output variables --
c
c      cndnum -- estimate of the 1-norm condition number
c
c  subprograms called --
c
c      xdsls3
c
c     ------------------------------------------------------------------
 
c     -----------------------
c     ...  passed parameters
c     -----------------------
 
      logical            unsym
 
      integer            neqns , nsuper
 
      integer            xsup(*), xpanel(*), xlindx(*), 
     1                   lindxg(*), pvtblk (*)
 
      double precision   tnorma, cndnum
 
      double precision   diag(*), lnz(*), z(*), temp1(*)
 
c     ----------------------
c     ...  local parameters
c     ----------------------
 
      integer            fstind, fstnod, i     , isuper, j     ,
     1                   jj    , k1    , k2    , 
     2                   kk1   , kk2   , l     , lstind, 
     3                   lstnod, node   

      integer            indbgn, ipanel, ipbgn , ipend , istart,
     1                   ldtemp, n1    , n2    , nodext
 
      double precision   e     , offdia, temp  ,
     1                   x1    , x2    , zasum
 
c     ---------------------
c     ...  subprograms used
c     ---------------------
 
      external           xdsls3
 
c     ------------------------------------------------------------------
 
c     ---------------
c     initialization
c     ---------------
 
      k1     = 1
      e      = 1.0
 
      z(1:neqns) = 0.0d0
 
c     --------------------------------------------------------------
c     forward solve -- pick rhs to cause large growth in solution of
c     l*d*z = b
c     --------------------------------------------------------------
 
      do 80 isuper = 1,nsuper
c.cebug
c     write(6,'("xdslc2 - 80 - isuper = ", i8)') isuper
c.debug

          ipbgn  = xsup(isuper)
          ipend  = xsup(isuper+1) - 1
          n2     = xpanel(ipend+1) - xpanel(ipbgn)

          fstind = xlindx(isuper)
          lstind = xlindx(isuper+1) - 1
          nodext = lstind - fstind + 1
c.debug
c     write(6,'("ipbgn, ipend, n2, fstind, lstind, nodext = ", 6i8)')
c    1            ipbgn, ipend, n2, fstind, lstind, nodext 
c.debug

          do 79 ipanel = ipbgn, ipend

              fstnod = xpanel ( ipanel )
              lstnod = xpanel ( ipanel+1 ) - 1
 
              n1     = lstnod - fstnod + 1
              k2     = k1 + n1 * ( n1 - 1 ) / 2

              n2     = n2 - n1
c.debug
c     write(6,'("ipanel, fstnod, lstnod, n1, k2, n2 = ", 6i8)')
c    1            ipanel, fstnod, lstnod, n1, k2, n2
c.debug
 
              node   = fstnod
 
   20         continue
              if ( node .le. lstnod ) then
 
                if ( pvtblk(node) .eq. 1 ) then
 
                   if ( z(node) .ne. 0. ) e = sign ( e, -z(node) )
                   x1      = e - z(node)
                   z(node) = x1 / diag(node)
                   kk1     = k1
                   kk2     = k2
c.debug
c     write(6,'("node, z(node) = ", i8,1pd15.5)') node, z(node)
c.debug
 
cdir$ ivdep
                   do j = node+1, lstnod
                       z(j) = z(j) + lnz(kk1) * x1
                       kk1  = kk1 + 1
                   enddo
c.debug
c     call xdslp5 ( 'z - after 30', neqns, z, 6 )
c.debug
 
cdir$ ivdep
                   do j = lstnod+1, lstnod+n2
                       z(j) = z(j) + lnz(kk2) * x1
                       kk2  = kk2 + 1
                   enddo
c.debug
c     call xdslp5 ( 'z - after 35', neqns, z, 6 )
c     call xislp3 ( 'lindxg', lstind-fstind+1, lindxg(fstind), 6 )
c.debug
 
cdir$ ivdep
                   do j = fstind, lstind
                       jj    = lindxg(j)
c.debug
c     write(6,'("j, jj, kk2 = ", 3i8)') j, jj, kk2
c.debug
                       z(jj) = z(jj) + lnz(kk2) * x1
                       kk2   = kk2 + 1
                   enddo
c.debug
c     call xdslp5 ( 'z - after 40', neqns, z, 6 )
c.debug
 
                   node = node + 1
                   k1   = kk1
                   k2   = kk2
                   go to 20
 
                else
 
                   if ( z(node) .ne. 0. ) e = sign ( e, -z(node) )
                   x1 = e - z(node)
                   if ( z(node+1) .ne. 0. ) e = sign ( e, -z(node+1) )
                   x2 = e - z(node+1)
 
                   offdia    = lnz(k1)
                   temp      = diag(node) * diag(node+1) - offdia ** 2
                   z(node  ) = (  diag(node+1)*x1 - offdia    *x2 )
     1                       / temp
                   z(node+1) = ( -offdia      *x1 + diag(node)*x2 )
     1                       / temp
 
                   kk1 = k1 + 1
                   kk2 = k2
 
                   l   = lstnod - node - 1
 
cdir$ ivdep
                   do j = node+2, lstnod
                       z(j) = z(j) + lnz(kk1)*x1 + lnz(kk1+l)*x2
                       kk1  = kk1 + 1
                   enddo
 
                   kk1 = kk1 + l
 
                   l = n2 + nodext
 
cdir$ ivdep
                   do j = lstnod+1, lstnod+n2
                       z(j) = z(j) + lnz(kk2)*x1 + lnz(kk2+l)*x2
                       kk2  = kk2 + 1
                   enddo
 
cdir$ ivdep
                   do j = fstind, lstind
                       jj    = lindxg(j)
                       z(jj) = z(jj) + lnz(kk2)*x1 + lnz(kk2+l)*x2
                       kk2   = kk2 + 1
                   enddo
 
                   kk2 = kk2 + l
 
                   node = node + 2
                   k1   = kk1
                   k2   = kk2
                   go to 20
 
                end if
 
              end if
 
              if ( unsym ) then
                  k1 = k2 + n1*(n1-1)/2 + n1 * ( n2 + nodext )
              else
                  k1 = k2
              end if
 
   79     continue
c.debug
c     write(6,'("after 79 loop")')
c.debug
 
   80 continue
c.debug
c     write(6,'("after 80 loop")')
c.debug
 
c     -------------------------
c     scale z to have unit norm
c     -------------------------
 
      zasum = 0.d0
      do i=1,neqns
        zasum = zasum + abs(z(i))
      enddo
      temp = 1.0 / zasum
c.debug
c     write(6,'("after 80 - temp = ", 1pd15.5)') temp
c     call xdslp5 ( 'z - after 80', neqns, z, 6 )
c.debug
      z(1:neqns) = temp*z(1:neqns)
 
c     --------------
c     backward solve
c     --------------
 
      do 140 isuper = nsuper, 1, -1
 
          ipbgn  = xsup(isuper)
          ipend  = xsup(isuper+1) - 1
          n2     = 0

          indbgn = xlindx ( isuper )
          nodext = xlindx ( isuper+1 ) - indbgn
c.debug
c     write(6,'("ipbgn, ipend, n2, indbgn, nodext = ", 6i8)')
c    1            ipbgn, ipend, n2, indbgn, nodext 
c     call xislp3 ( 'indices', nodext, lindxg(indbgn), 6 )
c.debug

          do ipanel = ipend, ipbgn, -1

              istart = xpanel ( ipanel )
              n1     = xpanel ( ipanel+1 ) - istart
c.debug
c     write(6,'("ipanel, istart, n1 = ", 6i8)')
c    1            ipanel, istart, n1
c.debug

c             ----------------------------------------------------
c             ... adjust pointer to the start of the current front
c             ----------------------------------------------------

              l  = n1 * ( n1 + n2 + nodext ) - ( n1 * ( n1 + 1 ) ) / 2
              k1 = k1 - l
c.debug
c     write(6,'("k1, l = ", 2i8)') k1, l
c     call xdslp5 ( 'subset of lnz', l, lnz(k1), 6 )
c.debug

c             ------------------------------------------------------
c             ... update rhs vectors with previously computed values
c                 and perform back substitution with diagonal block.
c             ------------------------------------------------------

              ldtemp = max ( 1, n2 + nodext )

              call xdsls6 ( istart-1, n1, n2, nodext, neqns, 1,
     1                      pvtblk, lnz(k1), z, lindxg(indbgn),
     2                      ldtemp, 1, temp1 )
c.debug
c     write(6,'("in xdslc2 after xdsls6 - ipanel = ",i8)') ipanel
c     call xdslp5 ( 'z - after xdsls6', neqns, z, 6 )
c.debug
      
              if ( unsym ) k1 = k1 - l
    
              n2 = n2 + n1
 
c             ------------------------------------------------
c             ... end of back substitution loop for this panel
c             ------------------------------------------------
 
          enddo
c.debug
c     write(6,'("in new xdslc2 after 139 - isuper = ",i8)') isuper
c     call xdslp5 ( 'z - after 139', neqns, z, 6 )
c.debug
 
c         ----------------------------------------------------
c         ... end of back substitution loop for this supernode
c         ----------------------------------------------------
 
  140 continue
 
c     -------------------------
c     scale z to have unit norm
c     -------------------------
 
      zasum = 0.d0
      do i=1,neqns
        zasum = zasum + abs(z(i))
      enddo
      temp = 1.0 / zasum
c.debug
c     write(6,'("after 140 - temp = ", 1pd15.5)') temp
c     call xdslp5 ( 'z - after 140', neqns, z, 6 )
c.debug
      z(1:neqns) = temp*z(1:neqns)
 
c     -------------------------
c     ... compute full solution
c     -------------------------
 
c.debug
c     write(6,'("in xdslc2 before xdsls3")')
c.debug
      call xdsls3( unsym,  neqns , nsuper, 1     , neqns ,
     1             xsup ,  xpanel, xlindx, lindxg, pvtblk, diag  ,
     2             lnz  ,  z    ,  1     , temp1   )
c.debug
c     write(6,'("in xdslc2 after  xdsls3")')
c.debug
 
c     -----------------------------------------------
c     ... compute norm of z and then condition number
c     -----------------------------------------------
 
      zasum = 0.d0
      do i=1,neqns
        zasum = zasum + abs(z(i))
      enddo
      temp = zasum
c.debug
c     write(6,'("after xdsls3 - temp = ", 1pd15.5)') temp
c     call xdslp5 ( 'z - after xdsls3', neqns, z, 6 )
c.debug
 
      cndnum = max ( tnorma * temp , 1.0d0 )
c     ------------------------------------------------------------------
c
      return
      end
