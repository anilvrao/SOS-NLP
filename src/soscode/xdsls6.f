      subroutine xdsls6 ( nodoff, n1   ,  n2    , nodext, ndimb,
     1                    nrhs,   pvtblk, matrix, rhs,    lindx,
     2                    ldtemp, slvbsz, temp1  )
c
c     created            -- 22-jun-87, rgg
c     last modified      -- 24-feb-89, rgg
c                           12-sep-89, rgg, to allow for pivoting
c                           09-nov-95, rgg, added temp1 and temp2
c                                           to do gather/scatter for
c                                           rhs, unrolled loops
c                           27-feb-97, rgg, converted to panels
c                                           switched to dgemv and dgemm
c                           27-mar-98, rgg, converted to processing 
c                                           slvbsz rhs at one time
c
c     purpose
c     -------
c
c         xdsls6 performs the back substitution step for a
c         given panel in a super node.
c
c     input variables
c     ---------------
c
c         nodoff  i   number of nodes already eliminated.
c         n1      i   number of nodes in the triangle of this panel
c         n2      i   number of nodes in the rectangle of this panel
c                     but still in the interior of the front
c         nodext  i   number of nodes in the exterior of the front
c         ndimb   i   leading dimension of the array b.
c         nrhs    i   number of right hand sides
c         pvtblk  i   pivoting information array
c         matrix  d   array holding the columns of the factorization
c                     associated with the super node
c         lindx   i   array holding the row indices for the entries
c                     in the column outside the diagonal block.
c         ldtemp  i   leading dimension of array temp1
c         slvbsz  i   number of columns in array temp1
c
c     working storage
c     ---------------
c
c         temp1   d   work vector of size ( mxlfrt, slvbsz )
c
c
c     input/output variables
c     ----------------------
c
c         rhs     d   array holding the right-hand/solution vectors
c
c     ================================================================
 
c     ------------------------
c     ... variable declaration
c     ------------------------
 
      integer             nodoff, n1,     n2,     nodext, ndimb,  
     1                    nrhs  , ldtemp, slvbsz
 
      integer             pvtblk(*),      lindx(*)
 
      double precision    matrix(*),      rhs(ndimb, *),
     1                    temp1(ldtemp,*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             fstnod, i,      irhs,   lstnod, node,   
     1                    jjrhs,  jrhs,
     2                    k1,     k2,     k3,     k4,     kcol,
     3                    k1sav,  k2sav,  k3sav,  krhs,   n3
 
      double precision    mone,   one,
     1                    s11,    s12,    s13,    s14,   
     1                    s21,    s22,    s23,    s24,
     2                    t1i,    t2i,
     3                    mk1,    mk2,    mk3,    mk4

      integer             mtxind, i1,     i2

      mtxind(i1,i2) = 1 + (i2-1)*i1 - (i2-1)*i2/2
 
c     ================================================================
 
c     ------------------------------------------------
c     ... update the nodes for this supernode with the
c         nodes previously computed.
c     ------------------------------------------------
 
      fstnod = nodoff + 1
      lstnod = nodoff + n1

      mone   = -1.0
      one    =  1.0

      n3     = n2 + nodext
      k1     = 1 + n1 * ( n1 - 1 ) / 2 
 
      if ( n3 .gt. 0 ) then
 
      do 130 irhs = 1, nrhs, slvbsz
 
          krhs = min ( slvbsz, nrhs - irhs + 1 )
c.debug
c     write(6,'("in xdsls6 - fstnod, lstnod, n1, n2 = ", 4i8)')
c    1                        fstnod, lstnod, n1, n2 
c     write(6,'("in xdsls6 - nodext, n3, irhs, krhs = ", 4i8)')
c    1                        nodext, n3, irhs, krhs
c.debug
c.debug
c     write(6,'("in xdsls5 - irhs, nrhs, slvbsz, krhs = ", 4i8)')
c    1                        irhs, nrhs, slvbsz, krhs 
c.debug

c         ----------------------------------------------
c         ... gather up info in rhs into temp1 and temp2
c         ----------------------------------------------

          do jjrhs = 1, krhs

              jrhs = irhs + jjrhs - 1

              temp1(1:n2,jjrhs) = rhs(lstnod+1:lstnod+n2,jrhs)

              call dgthr ( nodext, rhs(1,jrhs), 
     1                     temp1(n2+1,jjrhs), lindx )

          enddo

c         -------------------------------------------------
c         ... perform back elimination step
c             if krhs is small, use multiple calls to dgemv.
c             if krhs is more moderate, use one call to 
c             dgemm.  the breakpoint of 4 was chosen from
c             performance testing on a SGI Origin 2000 in
c             april 1998.
c         -------------------------------------------------

          if ( krhs .le. 4 ) then

              do jjrhs = 1, krhs

                  jrhs = irhs + jjrhs - 1 
 
                  call dgemv ( 't', n3, n1, mone, matrix(k1), n3,
     1                         temp1(1,jjrhs), 1,
     2                         one, rhs(fstnod,jrhs), 1 )

              enddo

          else
 
              call dgemm ( 't', 'n', n1, krhs, n3, mone, matrix(k1), n3,
     1                     temp1, ldtemp, one, rhs(fstnod,irhs), ndimb )

          end if
 
  130 continue
 
      endif

c-------------------------------------------------------------------
 
c     ------------------------------------------------------------
c     ... perform back substitution with the transpose of the unit
c         lower triangle of the frontal matrix
c     ------------------------------------------------------------
 
      do 240 irhs = 1, nrhs, 2
 
          krhs = min ( 2, nrhs - irhs + 1 )
c.debug
c     write(6,'("in 240  - irhs, krhs, lstnod     = ", 4i8)')
c    1                      irhs, krhs, lstnod     
c.debug

          if ( krhs .eq. 1 ) then 

c             -----------------------------
c             ... process one rhs at a time
c             -----------------------------
 
              node = lstnod + 1

c             -----------------------------------------------------
c             ... loop over the nodes by blocks of 3 or 4 until the 
c                 nodes are exhausted
c             -----------------------------------------------------

  210         continue
              if ( node .eq. fstnod ) go to 219

                  kcol = min ( 3, node - fstnod )

                  if ( pvtblk(node-kcol) .lt. 0 ) kcol = kcol + 1

                  node = node - kcol
c.debug
c     if ( nrhs .gt. 1 ) then
c     write(6,'("in xdsls6 - irhs, krhs, node, kcol = ", 4i8)')
c    1                        irhs, krhs, node, kcol 
c     end if
c.debug

                  if ( kcol .eq. 1 ) then

c                     -------------------------------------------
c                     ... handle a single column.  this must be a
c                         1x1 pivot.
c                     -------------------------------------------

                      s11 = 0.
                      k1  = mtxind(n1,node-nodoff  )

c.debug
c     write(6,'("before 211 loop - node, lstnod = ", 2i8)') 
c    1                              node, lstnod
c.debug
cdir$ ivdep
                      do i = node+1, lstnod
                          t1i = rhs(i,irhs  )
                          mk1 = matrix(k1)
                          s11 = s11 + t1i*mk1
                          k1  = k1 + 1
                      enddo

                      rhs(node,irhs) = rhs(node,irhs) - s11
c.debug
c     if ( nrhs .gt. 1 )then
c     write(6,'("node, irhs, rhs(node,irhs) = ", 2i8,1pd15.5)')
c    1            node, irhs, rhs(node,irhs)
c     end if
c.debug

                  else if ( kcol .eq. 2 ) then

c                     --------------------------------------------
c                     ... handle two columns.  this must be either
c                         2 1x1 pivots or 1 2x2 pivot.
c                     --------------------------------------------

                      s11   = 0.
                      s12   = 0.
                      k1sav = mtxind(n1,node-nodoff  ) 
                      k1    = k1sav + 1
                      k2    = mtxind(n1,node-nodoff+1)

c.debug
c     write(6,'("before 212 loop - node, lstnod = ", 2i8)') 
c    1                              node, lstnod
c.debug
cdir$ ivdep
                      do i = node+2, lstnod
                          t1i = rhs(i,irhs  )
                          mk1 = matrix(k1)
                          mk2 = matrix(k2)
                          s11 = s11 + t1i*mk1
                          s12 = s12 + t1i*mk2
                          k1            = k1 + 1
                          k2            = k2 + 1
                      enddo

                      rhs(node+1,irhs) = rhs(node+1,irhs) - s12

                      if ( pvtblk(node) .eq. 1 ) then
                          s11 = s11 + rhs(node+1,irhs)*matrix(k1sav)
                      end if

                      rhs(node,irhs) = rhs(node,irhs) - s11
c.debug 
c     if ( nrhs .gt. 1 ) then
c     write(6,'("node, irhs, rhs(node,irhs) = ", 2i8,1pd15.5)')
c    1            node, irhs, rhs(node,irhs)
c     write(6,'("            rhs(node+1,irhs) = ", 16x,1pd15.5)')
c    1                        rhs(node+1,irhs)
c     end if
c.debug

                  else if ( kcol .eq. 3 ) then

c                     --------------------------------------------
c                     ... handle three columns.  this must be either
c                         3 1x1 pivots or 1 1x1 and 1 2x2 pivot.
c                     --------------------------------------------

                      s11   = 0.
                      s12   = 0.
                      s13   = 0.
                      k1sav = mtxind(n1,node-nodoff  ) 
                      k2sav = mtxind(n1,node-nodoff+1) 
                      k1    = k1sav + 2
                      k2    = k2sav + 1
                      k3    = mtxind(n1,node-nodoff+2) 
c.debug
c     if ( nrhs .gt. 1 ) then
c     write(6,'("n1, node, nodoff = ", 3i8)') n1, node, nodoff
c     write(6,'("k1sav, k2sav, k3 = ", 3i8)') k1sav, k2sav, k3
c     write(6,'("k1, k2, k3   = ", 3i8)') k1, k2, k3
c     write(6,'("node, lstnod = ", 3i8)') node, lstnod
c     end if
c.debug

c.debug
c     write(6,'("before 213 loop - node, lstnod = ", 2i8)') 
c    1                              node, lstnod
c.debug
cdir$ ivdep
                      do i = node+3, lstnod
                          t1i = rhs(i,irhs  )
                          mk1 = matrix(k1)
                          mk2 = matrix(k2)
                          mk3 = matrix(k3)
                          s11 = s11 + t1i*mk1
                          s12 = s12 + t1i*mk2
                          s13 = s13 + t1i*mk3
                          k1  = k1 + 1
                          k2  = k2 + 1
                          k3  = k3 + 1
                      enddo
c.debug
c     if ( nrhs .gt. 1 ) then
c     write(6,'("after 213 s11, s12, s13 = ", 1p3d15.5)') 
c    1                      s11, s12, s13
c     end if
c.debug

                      rhs(node+2,irhs) = rhs(node+2,irhs) - s13

                      if ( pvtblk(node+1) .ne. 2 ) then
c.debug
c     if ( nrhs .gt. 1 ) then
c     write(6,'("node, irhs, k2sav = ", 3i8)')
c    1            node, irhs, k2sav
c     write(6,'("rhs(node+2,irhs), matrix(k2sav), s12 = ", 1p3d15.5)')
c    1            rhs(node+2,irhs), matrix(k2sav), s12
c     end if
c.debug
                          s12 = s12 + rhs(node+2,irhs)*matrix(k2sav)
                      end if

                      rhs(node+1,irhs) = rhs(node+1,irhs) - s12

                      if ( pvtblk(node) .eq. 1 ) then
                          s11 = s11 + rhs(node+1,irhs)*matrix(k1sav)
     1                              + rhs(node+2,irhs)*matrix(k1sav+1)
                      else
                          s11 = s11 + rhs(node+2,irhs)*matrix(k1sav+1)
                      end if

                      rhs(node,irhs) = rhs(node,irhs) - s11
c.debug
c     if ( nrhs .gt. 1 ) then
c     write(6,'("node, irhs, rhs(node,irhs) = ", 2i8,1pd15.5)')
c    1            node, irhs, rhs(node,irhs)
c     write(6,'("            rhs(node+1,irhs) = ", 16x,1pd15.5)')
c    1                        rhs(node+1,irhs)
c     write(6,'("            rhs(node+2,irhs) = ", 16x,1pd15.5)')
c    1                        rhs(node+2,irhs)
c     end if
c.debug

                  else 

c                     -----------------------------------------------
c                     ... handle four columns.  the first two columns
c                         are a 2x2 pivot.  the last two columns are
c                         either 2 1x2 pivots or 1 2x2 pivot.
c                     -----------------------------------------------

                      s11   = 0.
                      s12   = 0.
                      s13   = 0.
                      s14   = 0.
                      k1sav = mtxind(n1,node-nodoff  ) 
                      k2sav = mtxind(n1,node-nodoff+1) 
                      k3sav = mtxind(n1,node-nodoff+2) 
                      k1    = k1sav + 3
                      k2    = k2sav + 2
                      k3    = k3sav + 1
                      k4    = mtxind(n1,node-nodoff+3) 
c.debug
c     write(6,'("before 214 loop - node, lstnod = ", 2i8)') 
c    1                              node, lstnod
c.debug

cdir$ ivdep
                      do i = node+4, lstnod
                          t1i = rhs(i,irhs  )
                          mk1 = matrix(k1)
                          mk2 = matrix(k2)
                          mk3 = matrix(k3)
                          mk4 = matrix(k4)
                          s11 = s11 + t1i*mk1
                          s12 = s12 + t1i*mk2
                          s13 = s13 + t1i*mk3
                          s14 = s14 + t1i*mk4
                          k1  = k1 + 1
                          k2  = k2 + 1
                          k3  = k3 + 1
                          k4  = k4 + 1
                      enddo

                      rhs(node+3,irhs) = rhs(node+3,irhs) - s14

                      if ( pvtblk(node+2) .ne. 2 ) then
                          s13 = s13 + rhs(node+3,irhs)*matrix(k3sav)
                      end if

                      rhs(node+2,irhs) = rhs(node+2,irhs) - s13

                      s11 = s11 + rhs(node+2,irhs)*matrix(k1sav+1)
     1                          + rhs(node+3,irhs)*matrix(k1sav+2)
                      s12 = s12 + rhs(node+2,irhs)*matrix(k2sav  )
     1                          + rhs(node+3,irhs)*matrix(k2sav+1)

                      rhs(node+1,irhs) = rhs(node+1,irhs) - s12

                      if ( pvtblk(node) .eq. 1 ) then
                          s11 = s11 + rhs(node+1,irhs)*matrix(k1sav)
                      end if

                      rhs(node,irhs) = rhs(node,irhs) - s11
c.debug
c     if ( nrhs .gt. 1 ) then
c     write(6,'("node, irhs, rhs(node,irhs) = ", 2i8,1pd15.5)')
c    1            node, irhs, rhs(node,irhs)
c     write(6,'("            rhs(node+1,irhs) = ", 16x,1pd15.5)')
c    1                        rhs(node+1,irhs)
c     write(6,'("            rhs(node+2,irhs) = ", 16x,1pd15.5)')
c    1                        rhs(node+2,irhs)
c     write(6,'("            rhs(node+3,irhs) = ", 16x,1pd15.5)')
c    1                        rhs(node+3,irhs)
c     end if
c.debug

                  end if

                  go to 210
 
  219         continue

          else 

c             -----------------------------
c             ... process two rhs at a time
c             -----------------------------
 
              node = lstnod + 1

c             -----------------------------------------------------
c             ... loop over the nodes by blocks of 3 or 4 until the 
c                 nodes are exhausted
c             -----------------------------------------------------

  230         continue
              if ( node .eq. fstnod ) go to 239

                  kcol = min ( 3, node - fstnod )

                  if ( pvtblk(node-kcol) .lt. 0 ) kcol = kcol + 1

                  node = node - kcol

                  if ( kcol .eq. 1 ) then

c                     -------------------------------------------
c                     ... handle a single column.  this must be a
c                         1x1 pivot.
c                     -------------------------------------------

                      s11 = 0.
                      s21 = 0.
                      k1  = mtxind(n1,node-nodoff  )

cdir$ ivdep
                      do i = node+1, lstnod
                          t1i = rhs(i,irhs  )
                          t2i = rhs(i,irhs+1)
                          mk1 = matrix(k1)
                          s11 = s11 + t1i*mk1
                          s21 = s21 + t2i*mk1
                          k1  = k1 + 1
                      enddo

                      rhs(node,irhs  ) = rhs(node,irhs  ) - s11
                      rhs(node,irhs+1) = rhs(node,irhs+1) - s21

                  else if ( kcol .eq. 2 ) then

c                     --------------------------------------------
c                     ... handle two columns.  this must be either
c                         2 1x1 pivots or 1 2x2 pivot.
c                     --------------------------------------------

                      s11   = 0.
                      s12   = 0.
                      s21   = 0.
                      s22   = 0.
                      k1sav = mtxind(n1,node-nodoff  ) 
                      k1    = k1sav + 1
                      k2    = mtxind(n1,node-nodoff+1)

cdir$ ivdep
                      do i = node+2, lstnod
                          t1i = rhs(i,irhs  )
                          t2i = rhs(i,irhs+1)
                          mk1 = matrix(k1)
                          mk2 = matrix(k2)
                          s11 = s11 + t1i*mk1
                          s12 = s12 + t1i*mk2
                          s21 = s21 + t2i*mk1
                          s22 = s22 + t2i*mk2
                          k1            = k1 + 1
                          k2            = k2 + 1
                      enddo

                      rhs(node+1,irhs  ) = rhs(node+1,irhs  ) - s12
                      rhs(node+1,irhs+1) = rhs(node+1,irhs+1) - s22

                      if ( pvtblk(node) .eq. 1 ) then
                          s11 = s11 + rhs(node+1,irhs  )*matrix(k1sav)
                          s21 = s21 + rhs(node+1,irhs+1)*matrix(k1sav)
                      end if

                      rhs(node,irhs  ) = rhs(node,irhs  ) - s11
                      rhs(node,irhs+1) = rhs(node,irhs+1) - s21

                  else if ( kcol .eq. 3 ) then

c                     --------------------------------------------
c                     ... handle three columns.  this must be either
c                         3 1x1 pivots or 1 1x1 and 1 2x2 pivot.
c                     --------------------------------------------

                      s11   = 0.
                      s12   = 0.
                      s13   = 0.
                      s21   = 0.
                      s22   = 0.
                      s23   = 0.
                      k1sav = mtxind(n1,node-nodoff  ) 
                      k2sav = mtxind(n1,node-nodoff+1) 
                      k1    = k1sav + 2
                      k2    = k2sav + 1
                      k3    = mtxind(n1,node-nodoff+2) 

cdir$ ivdep
                      do i = node+3, lstnod
                          t1i = rhs(i,irhs  )
                          t2i = rhs(i,irhs+1)
                          mk1 = matrix(k1)
                          mk2 = matrix(k2)
                          mk3 = matrix(k3)
                          s11 = s11 + t1i*mk1
                          s12 = s12 + t1i*mk2
                          s13 = s13 + t1i*mk3
                          s21 = s21 + t2i*mk1
                          s22 = s22 + t2i*mk2
                          s23 = s23 + t2i*mk3
                          k1            = k1 + 1
                          k2            = k2 + 1
                          k3            = k3 + 1
                      enddo

                      rhs(node+2,irhs  ) = rhs(node+2,irhs  ) - s13
                      rhs(node+2,irhs+1) = rhs(node+2,irhs+1) - s23

                      if ( pvtblk(node+1) .ne. 2 ) then
                          s12 = s12 + rhs(node+2,irhs  )*matrix(k2sav)
                          s22 = s22 + rhs(node+2,irhs+1)*matrix(k2sav)
                      end if

                      rhs(node+1,irhs  ) = rhs(node+1,irhs  ) - s12
                      rhs(node+1,irhs+1) = rhs(node+1,irhs+1) - s22

                      if ( pvtblk(node) .eq. 1 ) then
                          s11 = s11 + rhs(node+1,irhs  )*matrix(k1sav)
     1                              + rhs(node+2,irhs  )*matrix(k1sav+1)
                          s21 = s21 + rhs(node+1,irhs+1)*matrix(k1sav)
     1                              + rhs(node+2,irhs+1)*matrix(k1sav+1)
                      else
                          s11 = s11 + rhs(node+2,irhs  )*matrix(k1sav+1)
                          s21 = s21 + rhs(node+2,irhs+1)*matrix(k1sav+1)
                      end if

                      rhs(node,irhs  ) = rhs(node,irhs  ) - s11
                      rhs(node,irhs+1) = rhs(node,irhs+1) - s21

                  else 

c                     -----------------------------------------------
c                     ... handle four columns.  the first two columns
c                         are a 2x2 pivot.  the last two columns are
c                         either 2 1x2 pivots or 1 2x2 pivot.
c                     -----------------------------------------------

                      s11   = 0.
                      s12   = 0.
                      s13   = 0.
                      s14   = 0.
                      s21   = 0.
                      s22   = 0.
                      s23   = 0.
                      s24   = 0.
                      k1sav = mtxind(n1,node-nodoff  ) 
                      k2sav = mtxind(n1,node-nodoff+1) 
                      k3sav = mtxind(n1,node-nodoff+2) 
                      k1    = k1sav + 3
                      k2    = k2sav + 2
                      k3    = k3sav + 1
                      k4    = mtxind(n1,node-nodoff+3) 

cdir$ ivdep
                      do i = node+4, lstnod
                          t1i = rhs(i,irhs  )
                          t2i = rhs(i,irhs+1)
                          mk1 = matrix(k1)
                          mk2 = matrix(k2)
                          mk3 = matrix(k3)
                          mk4 = matrix(k4)
                          s11 = s11 + t1i*mk1
                          s12 = s12 + t1i*mk2
                          s13 = s13 + t1i*mk3
                          s14 = s14 + t1i*mk4
                          s21 = s21 + t2i*mk1
                          s22 = s22 + t2i*mk2
                          s23 = s23 + t2i*mk3
                          s24 = s24 + t2i*mk4
                          k1            = k1 + 1
                          k2            = k2 + 1
                          k3            = k3 + 1
                          k4            = k4 + 1
                      enddo

                      rhs(node+3,irhs  ) = rhs(node+3,irhs  ) - s14
                      rhs(node+3,irhs+1) = rhs(node+3,irhs+1) - s24

                      if ( pvtblk(node+2) .ne. 2 ) then
                          s13 = s13 + rhs(node+3,irhs  )*matrix(k3sav)
                          s23 = s23 + rhs(node+3,irhs+1)*matrix(k3sav)
                      end if

                      rhs(node+2,irhs  ) = rhs(node+2,irhs  ) - s13
                      rhs(node+2,irhs+1) = rhs(node+2,irhs+1) - s23

                      s11 = s11 + rhs(node+2,irhs  )*matrix(k1sav+1)
     1                          + rhs(node+3,irhs  )*matrix(k1sav+2)
                      s12 = s12 + rhs(node+2,irhs  )*matrix(k2sav  )
     1                          + rhs(node+3,irhs  )*matrix(k2sav+1)
                      s21 = s21 + rhs(node+2,irhs+1)*matrix(k1sav+1)
     1                          + rhs(node+3,irhs+1)*matrix(k1sav+2)
                      s22 = s22 + rhs(node+2,irhs+1)*matrix(k2sav  )
     1                          + rhs(node+3,irhs+1)*matrix(k2sav+1)

                      rhs(node+1,irhs  ) = rhs(node+1,irhs  ) - s12
                      rhs(node+1,irhs+1) = rhs(node+1,irhs+1) - s22

                      if ( pvtblk(node) .eq. 1 ) then
                          s11 = s11 + rhs(node+1,irhs  )*matrix(k1sav)
                          s21 = s21 + rhs(node+1,irhs+1)*matrix(k1sav)
                      end if

                      rhs(node,irhs  ) = rhs(node,irhs  ) - s11
                      rhs(node,irhs+1) = rhs(node,irhs+1) - s21

                  end if

                  go to 230
 
  239         continue

          end if
 
  240 continue
 
c     ================================================================
 
c     ------------------------
c     ... end of module xdsls6
c     ------------------------
 
c.debug
c     write(6,'("at end of xdsls6 ")')
c.debug
      return
      end
