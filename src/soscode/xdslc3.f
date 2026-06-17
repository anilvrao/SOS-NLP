      subroutine xdslc3( unsym , neqns , nsuper, xsup  ,
     1                   xpanel, xlindx, lindxg, pvtblk, diag  , 
     2                   wafil1, wafil4, lmatrx, matrix, tnorma, 
     3                   z     , temp1 , cndnum, error   )
 
c
c  purpose -- to compute the condition number of the factorization
c             performed by xdslf2.
c
c  created            -- 07-jul-89, cca
c  last modifications -- 15-nov-95, rgg, adjusted for mods to solve
c                                        subroutines
c                        30-jul-98, rgg, added wafil4 for upper tri. 
c                                        factor
c                        31-aug-98, rgg, removed use of nzl and
c                                        modified i/o
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
c      wafil1 -- word addressable i/o file holding lower triangular
c                entries
c      wafil4 -- word addressable i/o file holding upper triangular
c                entries
c      lmatrx -- length of matrix array
c      tnorma -- 1 norm of orginial matrix
c
c  working storage
c
c      matrix -- work array to hold current front
c      z      -- work array of length neqns which holds the
c                approximate null vector.
c      temp1  -- work array of length mxlfrt for use by solve
c                subroutines
c
c  output variables --
c
c      cndnum -- estimate of the 1-norm condition number
c      error  -- error return
c                = 0 normal
c                =-1 i/o error on wafil1
c                =-2 i/o error on wafil4
c
c  subprograms called --
c
c      xdsls4, xdsls6
c
c     ------------------------------------------------------------------
 
c     -----------------------
c     ...  passed parameters
c     -----------------------
 
      logical            unsym
 
      integer            neqns , nsuper, wafil1, wafil4, 
     1                   lmatrx, error
 
      integer            xsup(*), xpanel(*), xlindx(*), 
     1                   lindxg(*), pvtblk (*)
 
      double precision   tnorma, cndnum
 
      double precision   diag(*), matrix(*), z(*), temp1(*)
 
c     ----------------------
c     ...  local parameters
c     ----------------------
 
      integer            fstind, fstnod, i     , isuper, j     ,
     1                   jj    , k1    , k2    , kk1   ,
     2                   kk2   , l     , lstind, lstnod, 
     3                   node  , length, istart, indbgn,
     4                   nodext

      integer            ipanel, ipbgn , ipend , 
     1                   ldtemp, n1    , n2    

      integer            iopos(2)
 
      double precision   e     , offdia, temp  , x1    , x2  , zasum
 
c     ---------------------
c     ...  subprograms used
c     ---------------------
 
      external           xdsls4, xdsls6
 
c     ------------------------------------------------------------------
 
c     ---------------
c     initialization
c     ---------------
 
      error    = 0
      iopos(1) = 1
      iopos(2) = 1
 
      e      = 1.0
 
      z(1:neqns) = 0.0d0
 
c     --------------------------------------------------------------
c     forward solve -- pick rhs to cause large growth in solution of
c     l*d*z = b
c     --------------------------------------------------------------
 
      do 80 isuper = 1,nsuper
c.debug
c     write(6,'("xdslc3 - 80 - isuper = ", i8)') isuper
c.debug

          ipbgn  = xsup(isuper)
          ipend  = xsup(isuper+1) - 1
          n2     = xpanel(ipend+1) - xpanel(ipbgn)

          fstind = xlindx(isuper)
          lstind = xlindx(isuper+1) - 1
          nodext = lstind - fstind + 1

c.debug
c     write(6,'("ipbgn, ipend, n2, fstind, lstind, nodext, k1 = ", 
c    1             6i8)')
c    1            ipbgn, ipend, n2, fstind, lstind, nodext, k1
c.debug

          do 79 ipanel = ipbgn, ipend

              fstnod = xpanel ( ipanel )
              lstnod = xpanel ( ipanel+1 ) - 1
 
              n1     = lstnod - fstnod + 1

              n2     = n2 - n1
c.debug
c     write(6,'("ipanel, fstnod, lstnod, n1, k2, n2 = ", 6i8)')
c    1            ipanel, fstnod, lstnod, n1, k2, n2
c.debug
 
              k1     = 1
              k2     = k1 + n1 * ( n1 - 1 ) / 2
              length = k2 + ( n2 + lstind - fstind + 1 ) * n1 - 1
 
              call xdslw5 ( wafil1, matrix, iopos, length, error )
              if ( error .ne. 0 ) go to 8000
c.debug
c     write(6,'("iopos, length = ", 3i8)') iopos, length
c     call xdslp5 ( 'subset of matrix', length, matrix, 6 )
c.debug
 
              call xdslw9 ( iopos, length )
 
              node   = fstnod
 
   20         continue
              if ( node .le. lstnod ) then
c.debug
c     write(6,'("in 20 loop - node = ", i8)') node
c.debug
 
                if ( pvtblk(node) .eq. 1 ) then
 
                   if ( z(node) .ne. 0. ) e = sign ( e, -z(node) )
                   x1      = e - z(node)
                   z(node) = x1 / diag(node)
                   kk1     = k1
                   kk2     = k2
 
cdir$ ivdep
                   do j = node+1, lstnod
                       z(j) = z(j) + matrix(kk1) * x1
                       kk1  = kk1 + 1
                   enddo
c.debug
c     call xdslp5 ( 'z - after 30', neqns, z, 6 )
c     write(6,'("after 30 - z(411) = ", 1pd15.5)') z(411)
c.debug


cdir$ ivdep                 
                   do j = lstnod+1, lstnod+n2
                       z(j) = z(j) + matrix(kk2) * x1
                       kk2  = kk2 + 1
                   enddo
c.debug
c     call xdslp5 ( 'z - after 35', neqns, z, 6 )
c     call xislp3 ( 'lindxg', lstind-fstind+1, lindxg(fstind), 6 )
c     write(6,'("after 35 - z(411) = ", 1pd15.5)') z(411)
c.debug
 
cdir$ ivdep
                   do j = fstind, lstind
                       jj    = lindxg(j)
c.debug
c     write(6,'("j, jj, kk2 = ", 3i8)') j, jj, kk2
c.debug
                       z(jj) = z(jj) + matrix(kk2) * x1
                       kk2   = kk2 + 1
                   enddo
c.debug
c     call xdslp5 ( 'z - after 40', neqns, z, 6 )
c     write(6,'("after 40 - z(411) = ", 1pd15.5)') z(411)
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
 
                   offdia    = matrix(k1)
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
                       z(j) = z(j) + matrix(kk1)*x1 + matrix(kk1+l)*x2
                       kk1  = kk1 + 1
                   enddo
 
                   kk1 = kk1 + l
 
                   l = n2 + nodext

cdir$ ivdep                        
                   do j = lstnod+1, lstnod+n2
                       z(j) = z(j) + matrix(kk2)*x1 + matrix(kk2+l)*x2
                       kk2  = kk2 + 1
                   enddo

cdir$ ivdep
                   do j = fstind, lstind
                       jj    = lindxg(j)
                       z(jj) = z(jj) + matrix(kk2)*x1 + matrix(kk2+l)*x2
                       kk2   = kk2 + 1
                   enddo
 
                   kk2 = kk2 + l
 
                   node = node + 2
                   k1   = kk1
                   k2   = kk2
                   go to 20
 
                end if
 
              end if
 
   79     continue   
c.debug
c     zasum = 0.d0
c     do i=1,neqns
c       zasum = zasum + abs(z(i))
c     enddo
c     temp = zasum
c     write(6,'("in 80 loop - isuper, temp = ", i8, 1pd15.5)') 
c    1                         isuper, temp
c.debug

   80 continue
 
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
 
c     ---------------------------------------------------------
c     ... perform the back substitution step for each supernode
c     ---------------------------------------------------------
 
      do 140 isuper = nsuper, 1, -1
c.debug
c     write(6,'("xdslc3 - 140 - isuper = ", i8)') isuper
c.debug

          ipbgn  = xsup(isuper)
          ipend  = xsup(isuper+1) - 1
          n2     = 0
         
          indbgn = xlindx ( isuper )
          nodext = xlindx ( isuper+1 ) - indbgn
c.debug
c     write(6,'("ipbgn, ipend, n2, fstind, lstind, nodext = ", 6i8)')
c    1            ipbgn, ipend, n2, fstind, lstind, nodext
c.debug

          do ipanel = ipend, ipbgn, -1
 
              istart = xpanel ( ipanel )
              n1     = xpanel ( ipanel+1 ) - istart
c.debug
c     write(6,'("ipanel, istart, n1, n2 = ", 6i8)')
c    1            ipanel, istart, n1, n2
c.debug
 
c             ----------------------------------------
c             ... read in columns of the factorization
c             ----------------------------------------
 
              l      = n1 * ( n1 + n2 + nodext ) 
     1               - ( n1 * ( n1 + 1 ) ) / 2
              call xdslw9 ( iopos, -l )
 
              if ( unsym ) then 
                  call xdslw5 ( wafil4, matrix, iopos, l, error )
                  if ( error .ne. 0 ) go to 8100
              else
                  call xdslw5 ( wafil1, matrix, iopos, l, error )
                  if ( error .ne. 0 ) go to 8000
              end if
c.debug
c     write(6,'("k1, l = ", 3i8)') iopos, l
c     call xdslp5 ( 'subset of matrix', l, matrix, 6 )
c.debug
     
c             ------------------------------------------------------
c             ... update rhs vectors with previously computed values
c                 and perform back substitution with diagonal block.
c             ------------------------------------------------------
 
              ldtemp = max ( 1, n2 + nodext )
 
              call xdsls6 ( istart-1, n1, n2, nodext, neqns, 1,
     1                      pvtblk, matrix, z, lindxg(indbgn),
     2                      ldtemp, 1, temp1 )
c.debug
c     call xdslp5 ( 'z - after xdsls6', neqns, z, 6 )
c     write(6,'("in xdslc3 after xdsls6 - ipanel = ",i8)') ipanel
c.debug

              n2 = n2 + n1

c             ------------------------------------------------
c             ... end of back substitution loop for this panel
c             ------------------------------------------------

           enddo
c.debug
c     call xdslp5 ( 'z - after 139 - xdsls6', neqns, z, 6 )
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
 
      call xdsls4( unsym , wafil1, wafil4, neqns , 1, neqns, z, 
     1             nsuper, xsup  , xpanel, xlindx, lindxg, diag,
     2             pvtblk, 1     , temp1 , lmatrx, matrix, error )
c.debug
c     write(6,'("in xdslc3 after xdsls4")') 
c.debug
 
      if ( error .ne. 0 ) go to 8000
 
c     -----------------------------------------------
c     ... compute norm of z and then condition number
c     -----------------------------------------------
 
      zasum = 0.d0
      do i=1,neqns
        zasum = zasum + abs(z(i))
      enddo
      temp = zasum
c.debug
c     write(6,'("at end of xdslc3 - temp = ", 1pd15.5)') temp
c     call xdslp5 ( 'z - after 80', neqns, z, 6 )
c.debug
 
      cndnum = max ( tnorma * temp , 1.0d0 )
 
      return
 
c     ------------------------------------------------------------------
 
c     ----------------------------
c     ... i/o error trap on wafil1
c     ----------------------------
 
 8000 continue
      error = -1
      return
 
c     ----------------------------
c     ... i/o error trap on wafil4
c     ----------------------------
 
 8100 continue
      error = -2
      return
      end
