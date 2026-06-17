      subroutine xdsls5 ( nodoff, n1,     n2,     nodext, ndimb,  
     1                    nrhs,   pvtblk, diag,   matrix,
     2                    rhs,    lindx , ldtemp, slvbsz, temp1  )
 
c
c     created              -- 22-jun-87, rgg
c     last modification    -- 24-feb-89, rgg
c                             12-sep-89, rgg, to allow for pivoting
c                             09-nov-95, rgg, added temp1 and temp2
c                                             to do gather/scatter for
c                                             rhs, unrolled loops
c                             27-feb-97, rgg, converted to panels
c                                             switched to dgemv and dgemm
c                             27-mar-98, rgg, converted to processing 
c                                             slvbsz rhs at one time
c                             10-jun-01, dkw, changed temp1(ldtemp,2)
c                                             to temp1(ldtemp,*)
c
c     purpose
c     -------
c
c         xdsls5 performs the forward elimination step for a
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
c         pvtblk  i   pivot block information
c         diag    d   diagonal array
c         matrix  d   array holding the columns of the factorization
c                     associated with the super node
c         lindx   i   array holding the row indices for the entries
c                     in the column outside the diagonal block.
c         ldtemp  i   leading dimension of temp1
c         slvbsz  i   number of columns in temp1
c
c     working storage
c     ---------------

c         temp1   d   work vector of size ( nodext, slvbsz )
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
 
      double precision    diag(*),  matrix(*), rhs(ndimb, *),
     1                    temp1(ldtemp,*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      integer             fstnod, i,      irhs,   kcol,   jjrhs,
     1                    jrhs  , k1,     k2,     k3,     k4,
     1                    krhs  , lstnod, node,   n3

      logical             q2by2
 
      double precision    det,    diag1,  diag2,  mone,  
     1                    offdia, one,    temp,   
     2                    s11,    s12,    s13,    s14,
     3                    s21,    s22,    s23,    s24

c     ================================================================
 
c     ------------------------------------------------------------
c     ... perform forward elimination with the unit lower triangle
c         of the frontal matrix
c     ------------------------------------------------------------
 
      fstnod = nodoff + 1
      lstnod = nodoff + n1
      q2by2  = .false.

      one    =  1.0
      mone   = -1.0
 
      do 40 irhs = 1, nrhs, 2
 
          krhs = min ( 2, nrhs - irhs + 1 )
c.debug
c     write(6,'("irhs, krhs, fstnod, lstnod = ", 4i8)')
c    1            irhs, krhs, fstnod, lstnod 
c.debug

          if ( krhs .eq. 1 ) then 

c             -----------------------------
c             ... process one rhs at a time
c             -----------------------------
 
              k1   = 1
              node = fstnod

c             -----------------------------------------------------
c             ... loop over the nodes by blocks of 3 or 4 until the 
c                 nodes are exhausted
c             -----------------------------------------------------

   10         continue
              if ( node .gt. lstnod - 1 ) go to 19

                  kcol = min ( 3, lstnod - node )
                  if ( pvtblk(node+kcol-1) .eq. 2 ) kcol = kcol + 1
c.debug
c     write(6,'("after 10 - k1, node, kcol = ", 3i8)')
c    1                       k1, node, kcol 
c.debug

                  if ( kcol .eq. 1 ) then

c                     -------------------------------------------
c                     ... handle a single column.  this must be a
c                         1x1 pivot.
c                     -------------------------------------------

                      s11 = -rhs(node  ,irhs  )
c.debug
c     write(6,'("before 11 - s11, node, lstnod = ", 1pd15.5, 2i8)')
c    1                        s11, node, lstnod 
c.debug

cdir$ ivdep
                      do i = node+1, lstnod
                          rhs(i,irhs  ) = rhs(i,irhs  ) + s11*matrix(k1)
                          k1            = k1 + 1
                      enddo

                  else if ( kcol .eq. 2 ) then

c                     --------------------------------------------
c                     ... handle two columns.  this must be either
c                         2 1x1 pivots or 1 2x2 pivot.
c                     --------------------------------------------

                      s11 = -rhs(node  ,irhs  )

                      k2  = k1 + ( lstnod - node )

                      if ( pvtblk(node) .eq. 1 ) then
                          rhs(node+1,irhs  ) = rhs(node+1,irhs  ) 
     1                                       + s11*matrix(k1)                  
                      else
                          q2by2 = .true.
                      end if

                      s12 = -rhs(node+1,irhs  )

                      k1  = k1 + 1

cdir$ ivdep
                      do i = node+2, lstnod
                          rhs(i,irhs  ) = rhs(i,irhs  ) + s11*matrix(k1)
     1                                                  + s12*matrix(k2)
                          k1            = k1 + 1
                          k2            = k2 + 1
                      enddo

                      k1 = k2 

                  else if ( kcol .eq. 3 ) then

c                     ------------------------------------------------
c                     ... handle three columns.  this must be either
c                         3 1x1 pivots or 1 1x1 pivot and 1 2x2 pivot.
c                     ------------------------------------------------

                      s11 = -rhs(node  ,irhs  )

                      k2  = k1 + ( lstnod - node )
                      k3  = k2 + ( lstnod - node - 1 )
    
                      if ( pvtblk(node) .eq. 1 ) then

                          rhs(node+1,irhs  ) = rhs(node+1,irhs  ) 
     1                                       + s11*matrix(k1)                  
                          rhs(node+2,irhs  ) = rhs(node+2,irhs  ) 
     1                                       + s11*matrix(k1+1)                  

                          s12 = -rhs(node+1,irhs  )

                          if ( pvtblk(node+1) .eq. 1 ) then

                              rhs(node+2,irhs  ) = rhs(node+2,irhs  ) 
     1                                           + s12*matrix(k2)

                          else

                              q2by2 = .true.

                          end if

                      else

                          s12   = -rhs(node+1,irhs  )
                          q2by2 = .true.

                          rhs(node+2,irhs  ) = rhs(node+2,irhs  ) 
     1                                       + s11*matrix(k1+1)                  
     1                                       + s12*matrix(k2)                  

                      end if

                      s13 = -rhs(node+2,irhs  )

                      k1  = k1 + 2
                      k2  = k2 + 1

cdir$ ivdep
                      do i = node+3, lstnod
                          rhs(i,irhs  ) = rhs(i,irhs  ) + s11*matrix(k1)
     1                                                  + s12*matrix(k2)
     2                                                  + s13*matrix(k3)
                          k1            = k1 + 1
                          k2            = k2 + 1
                          k3            = k3 + 1
                      enddo

                      k1 = k3 

                  else 

c                     --------------------------------------------------
c                     ... handle four columns.  the last pair of column
c                         are from a 2x2 pivot.  the first pair may be 2
c                         1x1 pivots or 1 2x2 pivot.
c                     --------------------------------------------------

                      s11   = -rhs(node  ,irhs  )
                      q2by2 = .true.
    
                      k2  = k1 + ( lstnod - node )
                      k3  = k2 + ( lstnod - node - 1 )
                      k4  = k3 + ( lstnod - node - 2 )

                      if ( pvtblk(node) .eq. 1 ) then

                          rhs(node+1,irhs  ) = rhs(node+1,irhs  ) 
     1                                       + s11*matrix(k1)                  
    
                      end if

                      s12 = -rhs(node+1,irhs  )

                      rhs(node+2,irhs  ) = rhs(node+2,irhs  ) 
     1                                   + s11*matrix(k1+1)                  
     1                                   + s12*matrix(k2  )                  
                      rhs(node+3,irhs  ) = rhs(node+3,irhs  ) 
     1                                   + s11*matrix(k1+2)                  
     1                                   + s12*matrix(k2+1)                  

                      s13 = -rhs(node+2,irhs  )
                      s14 = -rhs(node+3,irhs  )

                      k1  = k1 + 3
                      k2  = k2 + 2
                      k3  = k3 + 1

cdir$ ivdep
                      do i = node+4, lstnod
                          rhs(i,irhs  ) = rhs(i,irhs  ) + s11*matrix(k1)
     1                                                  + s12*matrix(k2)
     2                                                  + s13*matrix(k3)
     3                                                  + s14*matrix(k4)
                          k1            = k1 + 1
                          k2            = k2 + 1
                          k3            = k3 + 1
                          k4            = k4 + 1
                      enddo

                      k1 = k4

                  end if

                  node = node + kcol
                  go to 10
 
   19         continue

          else

c             -----------------------------
c             ... process two rhs at a time
c             -----------------------------
 
              k1   = 1
              node = fstnod

c             -----------------------------------------------------
c             ... loop over the nodes by blocks of 3 or 4 until the 
c                 nodes are exhausted
c             -----------------------------------------------------

   30         continue
              if ( node .gt. lstnod - 1 ) go to 39

                  kcol = min ( 3, lstnod - node )
                  if ( pvtblk(node+kcol-1) .eq. 2 ) kcol = kcol + 1

                  if ( kcol .eq. 1 ) then

c                     -------------------------------------------
c                     ... handle a single column.  this must be a
c                         1x1 pivot.
c                     -------------------------------------------

                      s11 = -rhs(node  ,irhs  )
                      s21 = -rhs(node  ,irhs+1)

cdir$ ivdep
                      do i = node+1, lstnod
                          rhs(i,irhs  ) = rhs(i,irhs  ) + s11*matrix(k1)
                          rhs(i,irhs+1) = rhs(i,irhs+1) + s21*matrix(k1)
                          k1            = k1 + 1
                      enddo

                  else if ( kcol .eq. 2 ) then

c                     --------------------------------------------
c                     ... handle two columns.  this must be either
c                         2 1x1 pivots or 1 2x2 pivot.
c                     --------------------------------------------

                      s11 = -rhs(node  ,irhs  )
                      s21 = -rhs(node  ,irhs+1)

                      k2  = k1 + ( lstnod - node )

                      if ( pvtblk(node) .eq. 1 ) then
                          rhs(node+1,irhs  ) = rhs(node+1,irhs  ) 
     1                                       + s11*matrix(k1)                  
                          rhs(node+1,irhs+1) = rhs(node+1,irhs+1) 
     1                                       + s21*matrix(k1)                  
                      else
                          q2by2 = .true.
                      end if

                      s12 = -rhs(node+1,irhs  )
                      s22 = -rhs(node+1,irhs+1)

                      k1  = k1 + 1

cdir$ ivdep
                      do i = node+2, lstnod
                          rhs(i,irhs  ) = rhs(i,irhs  ) + s11*matrix(k1)
     1                                                  + s12*matrix(k2)
                          rhs(i,irhs+1) = rhs(i,irhs+1) + s21*matrix(k1)
     1                                                  + s22*matrix(k2)
                          k1            = k1 + 1
                          k2            = k2 + 1
                      enddo

                      k1 = k2 

                  else if ( kcol .eq. 3 ) then

c                     ------------------------------------------------
c                     ... handle three columns.  this must be either
c                         3 1x1 pivots or 1 1x1 pivot and 1 2x2 pivot.
c                     ------------------------------------------------

                      s11 = -rhs(node  ,irhs  )
                      s21 = -rhs(node  ,irhs+1)

                      k2  = k1 + ( lstnod - node )
                      k3  = k2 + ( lstnod - node - 1 )
    
                      if ( pvtblk(node) .eq. 1 ) then

                          rhs(node+1,irhs  ) = rhs(node+1,irhs  ) 
     1                                       + s11*matrix(k1)                  
                          rhs(node+2,irhs  ) = rhs(node+2,irhs  ) 
     1                                       + s11*matrix(k1+1)                  

                          rhs(node+1,irhs+1) = rhs(node+1,irhs+1) 
     1                                       + s21*matrix(k1)                  
                          rhs(node+2,irhs+1) = rhs(node+2,irhs+1) 
     1                                       + s21*matrix(k1+1)                  

                          s12 = -rhs(node+1,irhs  )
                          s22 = -rhs(node+1,irhs+1)

                          if ( pvtblk(node+1) .eq. 1 ) then

                              rhs(node+2,irhs  ) = rhs(node+2,irhs  ) 
     1                                           + s12*matrix(k2)
                              rhs(node+2,irhs+1) = rhs(node+2,irhs+1) 
     1                                           + s22*matrix(k2)

                          else

                              q2by2 = .true.

                          end if

                      else

                          s12   = -rhs(node+1,irhs  )
                          s22   = -rhs(node+1,irhs+1)
                          q2by2 = .true.

                          rhs(node+2,irhs  ) = rhs(node+2,irhs  ) 
     1                                       + s11*matrix(k1+1)                  
     1                                       + s12*matrix(k2)                  
                          rhs(node+2,irhs+1) = rhs(node+2,irhs+1) 
     1                                       + s21*matrix(k1+1)                  
     1                                       + s22*matrix(k2)                  

                      end if

                      s13 = -rhs(node+2,irhs  )
                      s23 = -rhs(node+2,irhs+1)

                      k1  = k1 + 2
                      k2  = k2 + 1

cdir$ ivdep
                      do i = node+3, lstnod
                          rhs(i,irhs  ) = rhs(i,irhs  ) + s11*matrix(k1)
     1                                                  + s12*matrix(k2)
     2                                                  + s13*matrix(k3)
                          rhs(i,irhs+1) = rhs(i,irhs+1) + s21*matrix(k1)
     1                                                  + s22*matrix(k2)
     2                                                  + s23*matrix(k3)
                          k1            = k1 + 1
                          k2            = k2 + 1
                          k3            = k3 + 1
                      enddo

                      k1 = k3 

                  else 

c                     --------------------------------------------------
c                     ... handle four columns.  the last pair of column
c                         are from a 2x2 pivot.  the first pair may be 2
c                         1x1 pivots or 1 2x2 pivot.
c                     --------------------------------------------------

                      s11   = -rhs(node  ,irhs  )
                      s21   = -rhs(node  ,irhs+1)
                      q2by2 = .true.
    
                      k2  = k1 + ( lstnod - node )
                      k3  = k2 + ( lstnod - node - 1 )
                      k4  = k3 + ( lstnod - node - 2 )

                      if ( pvtblk(node) .eq. 1 ) then

                          rhs(node+1,irhs  ) = rhs(node+1,irhs  ) 
     1                                       + s11*matrix(k1)                  
                          rhs(node+1,irhs+1) = rhs(node+1,irhs+1) 
     1                                       + s21*matrix(k1)                  
    
                      end if

                      s12 = -rhs(node+1,irhs  )
                      s22 = -rhs(node+1,irhs+1)

                      rhs(node+2,irhs  ) = rhs(node+2,irhs  ) 
     1                                   + s11*matrix(k1+1)                  
     1                                   + s12*matrix(k2  )                  
                      rhs(node+3,irhs  ) = rhs(node+3,irhs  ) 
     1                                   + s11*matrix(k1+2)                  
     1                                   + s12*matrix(k2+1)                  

                      rhs(node+2,irhs+1) = rhs(node+2,irhs+1) 
     1                                   + s21*matrix(k1+1)                  
     1                                   + s22*matrix(k2  )                  
                      rhs(node+3,irhs+1) = rhs(node+3,irhs+1) 
     1                                   + s21*matrix(k1+2)                  
     1                                   + s22*matrix(k2+1)                  

                      s13 = -rhs(node+2,irhs  )
                      s14 = -rhs(node+3,irhs  )
                      s23 = -rhs(node+2,irhs+1)
                      s24 = -rhs(node+3,irhs+1)

                      k1  = k1 + 3
                      k2  = k2 + 2
                      k3  = k3 + 1

cdir$ ivdep
                      do i = node+4, lstnod
                          rhs(i,irhs  ) = rhs(i,irhs  ) + s11*matrix(k1)
     1                                                  + s12*matrix(k2)
     2                                                  + s13*matrix(k3)
     3                                                  + s14*matrix(k4)
                          rhs(i,irhs+1) = rhs(i,irhs+1) + s21*matrix(k1)
     1                                                  + s22*matrix(k2)
     2                                                  + s23*matrix(k3)
     3                                                  + s24*matrix(k4)
                          k1            = k1 + 1
                          k2            = k2 + 1
                          k3            = k3 + 1
                          k4            = k4 + 1
                      enddo

                      k1 = k4

                  end if

                  node = node + kcol
                  go to 30
 
   39         continue

          end if
c.debug
c     write(6,'("at end of 40 loop")')
c.debug
 
   40 continue
 
c     ----------------------------------------------------
c     ... update the remaining nodes in the frontal matrix
c     ----------------------------------------------------
 
      if ( n2 + nodext .le. 0 ) go to 190

      n3 = n2 + nodext
      k1 = 1 + n1 * ( n1 - 1 ) / 2 
 
      do 130 irhs = 1, nrhs, slvbsz
 
          krhs = min ( slvbsz, nrhs - irhs + 1 )
c.debug
c     write(6,'("in xdsls5 - irhs, nrhs, slvbsz, krhs = ", 4i8)')
c    1                        irhs, nrhs, slvbsz, krhs 
c.debug

c         ------------------------------
c         ... gather rhs info into temp1 
c         ------------------------------

          do jjrhs = 1, krhs

              jrhs = irhs + jjrhs - 1 

              temp1(1:n2,jjrhs) = rhs(lstnod+1:lstnod+n2,jrhs)

              call dgthr ( nodext, rhs(1,jrhs  ), 
     1                             temp1(n2+1,jjrhs), lindx )

          enddo

c         -------------------------------------------------
c         ... perform update of info in temp1 and temp2
c             if krhs is small, use multiple calls to dgemv.
c             if krhs is more moderate, use one call to 
c             dgemm.  the breakpoint of 4 was chosen from
c             performance testing on a SGI Origin 2000 in
c             april 1998.
c         -------------------------------------------------

          if ( krhs .le. 4 ) then

              do jjrhs = 1, krhs

                  jrhs = irhs + jjrhs - 1 
 
                  call dgemv ( 'n', n3, n1, mone, matrix(k1), n3,
     1                         rhs(fstnod,jrhs), 1,
     2                         one, temp1(1,jjrhs), 1 )

              enddo

          else
 
              call dgemm ( 'n', 'n', n3, krhs, n1, mone, matrix(k1), n3,
     1                     rhs(fstnod,irhs), ndimb, one, temp1, ldtemp )

          end if

c         --------------------------------------------------
c         ... distribute info in temp1 and temp2 back in rhs
c         --------------------------------------------------

          do jjrhs = 1, krhs

              jrhs = irhs + jjrhs - 1 

              rhs(lstnod+1:lstnod+n2,jrhs) = temp1(1:n2,jjrhs)

              call dsctr ( nodext, temp1(n2+1,jjrhs), lindx, 
     1                         rhs(1,jrhs  ) )

          enddo
 
  130 continue

c     ================================================================
 
c     -------------------------------
c     ... scale by the diagonal block
c     -------------------------------
 
  190 continue
      node = fstnod

      if ( q2by2 ) then 
 
  200     continue
          if (node .le. lstnod ) then
 
              if ( pvtblk (node) .ne. 2 ) then
 
                 do irhs = 1, nrhs
                     rhs(node,irhs) = rhs(node,irhs) / diag(node)
                 enddo
 
                 node = node + 1
 
              else
 
                 i      = abs (pvtblk (node+1))
                 offdia = matrix(i)
 
                 diag1  = diag (node)
                 diag2  = diag (node+1)
                 det    = diag1 * diag2 - offdia ** 2
                 temp   = diag1
                 diag1  = diag2 / det
                 diag2  = temp  / det
                 offdia = -offdia / det
 
                 do irhs = 1, nrhs
                     temp   = rhs (node, irhs)
                     rhs (node, irhs)   = diag1  * temp
     1                                  + offdia * rhs (node+1, irhs)
                     rhs (node+1, irhs) = offdia * temp
     1                                  + diag2  * rhs (node+1, irhs)
                 enddo
 
                 node = node + 2
     
              endif
     
              go to 200

          end if

      else

c         -----------------------
c         ... strictly 1x1 blocks
c         -----------------------

          do irhs = 1, nrhs
              do node = fstnod, lstnod
                  rhs(node,irhs) = rhs(node,irhs) / diag(node)
              enddo
          enddo
 
      endif
c.debug
c     write(6,'("after diagonal scaling")')
c.debug
 
c     ================================================================
 
c     ------------------------
c     ... end of module xdsls5
c     ------------------------
 
      return
      end
