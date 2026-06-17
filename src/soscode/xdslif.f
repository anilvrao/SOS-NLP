      subroutine xdslif ( work,   lwork,  needs,  error )
 
c
c     purpose
c     -------
c
c     xdslif is the top level driver for finalization of structural
c     input of the matrix.
c
c     created         27-jan-89   -- rgg --
c     modified        08-feb-91   -- rgg -- mods to allow no i/o to
c                                           sqfile
c                     08-feb-91   -- rgg -- added use of xdslil
c                     04-dec-96   -- rgg -- complete rewrite to allow
c                                           out-of-core assembly
c                     09-jun-00   -- rgg -- corrected an error in use
c                                           of xdslil
c
c     input arguments
c     ---------------
c
c     lwork       i   length of work array.
c
c     input/output arguments
c     ----------------------
c
c     work        d   work array.  on input it holds the linked list
c                     representation of the assembled matrix structure
c                     built by xisli3 and xisli4.  on output it holds
c                     the full adjacency structure.
c
c     output arguments
c     ----------------
c
c     needs       i   amount of workspace required for the next stage.
c
c     error       i   error flag
c                     =    0  normal return
c                     = -100  incorrect processing path
c                     = -101  lwork not large enough.
c                     = -110  i/o error on wafil1 and wafil2
c
c     ------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      integer             lwork,  error,  needs
 
      double precision    work(*)
 
c     -------------------------------------
c     ... include global.CMNication area
c     -------------------------------------
 
      include '../commons/bcsext4.CMN'                                           
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------
 
      integer             acolls, adjncy, anrecd, anzero, arowls, 
     1                    cmpmap, colptr, inadj , ianrec,
     2                    incmp , incmp1, ineqns, ineqn1, inuse , 
     3                    inzcmp, itemp1, itemp2, itemp3, imxrec,
     4                    itemp4, itemp5, itemp6, itemp7, maxrec,
     5                    lstcol, maxnz , msglvl, mxused, nadj  ,
     6                    ncomp , neqns , rowlst, nzcomp, nzla  ,
     7                    output, sort  , sqfile, srtlst, stage ,
     8                    wafil1, wafil2, wkreqd, xadj  , xrowls

      integer             collst, colls2, colls3, colstr, ianzer,
     1                    ilen  , ilset , iwaln1, iwaln2, iwatr1,
     2                    iwatr2, k     , lset  , maxrow,
     3                    rowls2, rowls3, temp
 
      double precision    t1,     t2,     t3,
     1                    w1,     w2,     w3
 
c     --------------------
c     ... subprograms used
c     --------------------
 
      integer             xdslil, xdslni
 
      external            xislmv, xisliv, xisliw, xisliz, 
     1                    xdslil, xdslni, xdslt2, xislp1, xislp2
 
c---------------------------------------------------------------------
 
      error  = 0
 
      stage  = work ( qstage )
      msglvl = work ( qmsglv )
      output = work ( qoutpu )
      neqns  = work ( qneqns )

      ineqns = xdslni ( neqns )
      ineqn1 = xdslni ( neqns + 1 )

      maxrec = min ( neqns, mxrecd )
      imxrec = xdslni(maxrec)

      t1 = work ( qinptm )
      w1 = work ( qinpwl )
      call xdslt2 ( t1, w1, t2, w2 )
c.debug
c     write(6,'("at start of xdslif - accum. time = ", f15.6)') t2
c.debug
 
      if ( stage .eq. 1. ) then
 
c         ---------------------------
c         ... process diagonal matrix
c         ---------------------------
 
          ncomp  = neqns
          nzcomp = 0
          anzero = 0
          nzla   = 0

          incmp  = xdslni ( ncomp )
          incmp1 = xdslni ( ncomp + 1 )
          inzcmp = xdslni ( nzcomp )
 
          nadj    = 2 * nzcomp
          inadj   = xdslni ( nadj   )

          xrowls = lncomm + 1
          rowlst = xrowls + incmp1
          cmpmap = rowlst + ineqns
          xadj   = cmpmap + ineqns
          adjncy = xadj   + ineqn1

          mxused = adjncy + inadj  - 1

          call iramp ( neqns+1, work(xrowls), 1 )
          call iramp ( neqns  , work(rowlst), 1 )
          call iramp ( neqns  , work(cmpmap), 1 )
          call izero ( neqns+1, work(xadj  ), 1 )
 
      else if ( stage .eq. 2 ) then
 
c         ------------------------------
c         ... process nondiagonal matrix
c         ------------------------------

          maxnz  = work(qnzab )
          anzero = work(qnzla )
          anrecd = work(qtemp1)
          lstcol = work(qtemp2)
          sort   = work(qtemp3)
          
          arowls = work(qarowi)
          acolls = work(qacols) 
          srtlst = work(qtemp4)
 
          wafil1   = work(qwafl1)

c         ----------------------------------
c         ... allocate temporary work arrays 
c         ----------------------------------

          itemp1 = lncomm + 1
          itemp2 = itemp1 + imxrec
          itemp3 = itemp2 + imxrec

c.debug
          if ( itemp3 + ineqn1 .ne. acolls ) then
              write(6,'("storage oops in xdslif")')
              write(6,'("ineqn1, acolls = ", 3i8)')
     1                    ineqn1, acolls 
              stop
          end if
c.debug

          ilen = xdslil ( lwork - lncomm )

c         ---------------------------------------
c         ... sort and compress rowlst and collst
c         ---------------------------------------

c.debug
c     write(6,'(/ "in xdslif before xisliz")')
c     write(6,'("anzero, maxnz  = ", 2i8)') anzero, maxnz
c     write(6,'("anrecd, wafil1 = ", 2i8)') anrecd, wafil1
c     write(6,'("arowls, acolls = ", 2i8)') arowls, acolls
c
c     call xislp3 ( 'row    list', anzero, work(arowls), output )
c     call xislp3 ( 'column list', anzero, work(acolls), output )
c.debug
c.debug
c     t1 = work ( qinptm )
c     w1 = work ( qinpwl )
c     call xdslt2 ( t1, w1, t3, w3 )
c     write(6,'(/15x,"prior to xisliz        - accum. time = ", 
c    1          f15.6)') t3
c.debug

          call xisliz ( .true., neqns, anzero, maxnz, anrecd, 
     1                  maxrec, lstcol, sort, wafil1, ilen, 
     2                  work(arowls), work(acolls), work(srtlst),
     3                  work(itemp1), work(itemp2), work(itemp3), 
     4                  error )

          if ( error .ne. 0 ) go to 8900

c.debug
c     write(6,'(15x,"in xdslif after xisliz")')
c     write(6,'(15x,"anzero, maxnz , anrecd = ", 3i8)') 
c    1                anzero, maxnz , anrecd
c
c     call xislp3 ( 'row    list', anzero, work(arowls), output )
c     call xislp3 ( 'column list', anzero, work(acolls), output )
c     if ( anrecd .gt. 0 ) then
c     call xislp3 ( 'record pos.', anrecd, work(itemp1), output )
c     call xislp3 ( 'record len.', anrecd, work(itemp2), output )
c     end if
c.debug
c.debug
c     t1 = work ( qinptm )
c     w1 = work ( qinpwl )
c     call xdslt2 ( t1, w1, t3, w3 )
c     write(6,'(15x,"after    xisliz        - accum. time = ", 
c    1          f15.6)') t3
c.debug

c
c     ===============================================================
c

c         -------------------------------------------------
c         ... remainder of processing depends on whether
c             the rowlst/collst structure is in memory or 
c             not
c             in either case the result is a compressed
c             representation of the strict lower triangular
c             matrix representation
c         -------------------------------------------------

          if ( anrecd .eq. 0 ) then

c             ----------------------------------------
c             ... rowlst and collst are held in memory.
c                 rearrange for in-memory processing
c             ----------------------------------------

              ianzer = xdslni ( anzero )
              temp   = arowls
              arowls = lwork + 1 - ianzer

              call xislmv ( anzero, work, xdslil(temp  -1)+1, 
     1                                    xdslil(arowls-1)+1 )

              temp   = acolls
              acolls = arowls - ianzer

              call xislmv ( anzero, work, xdslil(temp  -1)+1, 
     1                                    xdslil(acolls-1)+1 )

              itemp2 = itemp1 + ineqn1
              itemp3 = itemp2 + ineqn1
              itemp4 = itemp3 + ineqn1
              itemp5 = itemp4 + ineqn1
              itemp6 = itemp5 + ineqn1
              itemp7 = itemp6 + ineqn1
              colptr = itemp7 + ineqn1

c.debug
              if ( colptr + ineqn1 .gt. acolls ) then
                  write(6,'("storage oops 2 in xdslif")')
                  write(6,'("colptr, ineqn1, acolls = ", 3i8)')
     1                        colptr, ineqn1, acolls 
                  stop
              end if
c.debug

c             ---------------------------------------------------
c             ... generate two random permutations to be used for
c                 computing 2 of the 3 checksums.  permutations
c                 will be stored in itemp1 and itemp2.
c             ---------------------------------------------------

              call xisliw ( neqns, work(itemp1), 
     1                      work(itemp2), work(itemp3) )

c             -----------------------------------------------------
c             ... find the graph compression and build a compressed
c                 column major representation of the strict lower
c                 triangle of the original matrix.
c             -----------------------------------------------------

              nzla   = anzero 
              mxused = lwork - 2*maxnz + 2*anzero

c.debug
c     write(6,'("in xdslif before in-core xisliv")')
c     write(6,'("anzero, maxnz = ", 2i8)') anzero, maxnz
c.debug

              call xisliv ( neqns , anzero, maxnz,
     1                      work(colptr), work(arowls), work(acolls),
     2                      work(itemp1), work(itemp2), work(itemp3),
     3                      work(itemp4), work(itemp5), work(itemp6),
     4                      work(itemp7), ncomp , nzcomp )
 
              if ( msglvl .ge. 1 ) then 
                  write ( output, 81000 ) neqns, anzero, ncomp, nzcomp
              end if
c.debug
c     t1 = work ( qinptm )
c     w1 = work ( qinpwl )
c     call xdslt2 ( t1, w1, t3, w3 )
c     write(6,'(15x,"after    xisliv        - accum. time = ", 
c    1         f15.6)') t3
c.debug

c             --------------------------------------
c             ... compressed xadj   is in itemp1 and
c                 compressed adjncy is in acolls
c             --------------------------------------

          else

c             --------------------------------------------
c             ... rowlst and collst are held out-of-memory
c                 note that the last call to xisliz flushed
c                 all of rowlst and collst to wafil1.
c                 work(itemp1) holds record position info
c                 work(itemp2) holds record length   info
c             --------------------------------------------

              mxused = lwork

c             ---------------------------------------------------
c             ... perform a sort merge of rowlst/collst which
c                 is held on file wafil1.  The sorted result will
c                 reside on the unit number stored in wafil2.
c                 the unit numbers stored in wafil1 and wafil2 
c                 will flip-flop.
c             ---------------------------------------------------

              wafil2 = work(qwafl2)
c.debug
c     write(6,'("before open in wafil2 = ", i8)') wafil2
c.debug
             
              call xislw4 ( wafil2, error )
              if ( error .ne. 0 ) go to 8900  

              ianrec = xdslni ( anrecd )

              k      = itemp2
              itemp2 = itemp1 + ianrec

              call xislmv ( anrecd, work, xdslil(k     -1)+1, 
     1                                    xdslil(itemp2-1)+1)

              itemp3 = itemp2 + ianrec
              itemp4 = itemp3 + ianrec
              itemp5 = itemp4 + ianrec

              lset  = xdslil ( ( lwork - itemp5 ) / 8 )
c.debug
c             lset  = min ( lset, 45 )
c.debug
              ilset = xdslni ( lset )

              rowlst = itemp5
              collst = rowlst + ilset
              rowls2 = collst + ilset
              colls2 = rowls2 + ilset
              rowls3 = colls2 + ilset
              colls3 = rowls3 + 2*ilset
c.debug
c     write(6,'(15x,"in xdslif before out-of-core xislis")')
c     write(6,'(15x,"anrecd, lset = ", 2i8)') anrecd, lset
c.debug

              call xislis ( anrecd, sort,
     1                      wafil1, work(itemp1), work(itemp2),
     2                      iwaln1, iwatr1, 
     3                      wafil2, work(itemp3), work(itemp4),
     4                      iwaln2, iwatr2, lset,
     5                      work(rowlst), work(collst), 
     6                      work(rowls2), work(colls2), 
     7                      work(rowls3), work(colls3), 
     8                      anzero, error )

              if ( error .ne. 0 ) go to 8900

              work(qwaln1) = xdslni ( iwaln1 )
              work(qwatr1) = xdslni ( iwatr1 )
              work(qwaln2) = xdslni ( iwaln2 )
              work(qwatr2) = xdslni ( iwatr2 )

              nzla   = anzero 
c.debug
c     t1 = work ( qinptm )
c     w1 = work ( qinpwl )
c     call xdslt2 ( t1, w1, t3, w3 )
c     write(6,'(15x,"after    xislis        - accum. time = ", 
c    1         f15.6)') t3
c.debug

c             ------------------------------------------------
c             ... repartition work space for graph compression
c                 phase (out-of-core processing)
c             ------------------------------------------------

              itemp2 = itemp1 + ineqn1
              itemp3 = itemp2 + ineqn1
              itemp4 = itemp3 + ineqn1
              itemp5 = itemp4 + ineqn1
              itemp6 = itemp5 + ineqn1
              itemp7 = itemp6 + ineqn1
              colptr = itemp7 + ineqn1

c             ---------------------------------------------------
c             ... generate two random permutations to be used for
c                 computing 2 of the 3 checksums.  permutations
c                 will be stored in itemp1 and itemp2.
c             ---------------------------------------------------

              call xisliw ( neqns, work(itemp1), 
     1                      work(itemp2), work(itemp3) )
c.debug
c     t1 = work ( qinptm )
c     w1 = work ( qinpwl )
c     call xdslt2 ( t1, w1, t3, w3 )
c     write(6,'(15x,"after    xisliw        - accum. time = ", 
c    1         f15.6)') t3
c.debug

c             -----------------------------------------------------
c             ... find the graph compression and build a compressed
c                 column major representation of the strict lower
c                 triangle of the original matrix.
c             -----------------------------------------------------

              acolls = colptr + ineqn1
              k      = lwork - acolls + 1
              maxnz  = xdslil ( k / 2 )

              if ( maxnz .lt. neqns ) then
                  wkreqd = 10 * ineqn1 + lncomm
                  go to 8100
              end if

              arowls = acolls + xdslni ( maxnz )
 
              call xisliu ( neqns , anzero, maxnz, wafil2, output,
     1                      work(colptr), work(arowls), work(acolls),
     1                      work(acolls),
     2                      work(itemp1), work(itemp2), work(itemp3),
     3                      work(itemp4), work(itemp5), work(itemp6),
     4                      work(itemp7), 
     5                      ncomp , nzcomp, error )
c.debug
c     t1 = work ( qinptm )
c     w1 = work ( qinpwl )
c     call xdslt2 ( t1, w1, t3, w3 )
c     write(6,'(15x,"after    xisliu        - accum. time = ", 
c    1         f15.6)') t3
c.debug

              if ( error .eq. -2 ) go to 8900
 
              if ( msglvl .ge. 1 ) then 
                  write ( output, 81000 ) neqns, anzero, ncomp, nzcomp
              end if

              if ( error .eq. -1 ) then
                  wkreqd = lncomm + 8*ineqn1 + xdslni(2*nzcomp)
                  go to 8100
              end if

c             -------------------
c             ... close i/o files
c             -------------------
c.debug
c     write(6,'("before closing files")')
c.debug

              call xislw3 ( wafil1, error )
              if ( error .ne. 0 ) go to 8900  

              call xislw3 ( wafil2, error )
              if ( error .ne. 0 ) go to 8900  
c.debug
c     write(6,'("after  closing files")')
c.debug

          end if

c
c     ===============================================================
c

c         ---------------------------------------------------------
c         ... the compressed representation of the matrix structure
c             of the lower triangle has now been built in memory
c
c             itemp2 has xrowls - pointer to rowlst, a list of column
c                                 numbers for each compressed node
c                                 ncomp+1 long
c             itemp3 has rowlst - list of column numbers for each
c                                 compress node, neqns long
c             itemp4 has cmpmap - a map from row numbers to
c                                 compressed nodes, neqns long
c             itemp5 has xadj   - column start pointer for each 
c                                 compressed node, neqns+1 long
c             collst has ladjnc - the compressed row indices of the
c                                 lower triangle nzcomp long
c         ---------------------------------------------------------
c.debug
c     write(6,'("before storage constants computation")')
c.debug

          incmp  = xdslni ( ncomp )
          incmp1 = xdslni ( ncomp + 1 )
          inzcmp = xdslni ( nzcomp )
 
          nadj   = 2*nzcomp
          inadj  = xdslni ( nadj   )

          xrowls = lncomm + 1
          rowlst = xrowls + incmp1
          cmpmap = rowlst + ineqns
          colstr = cmpmap + ineqns
          adjncy = colstr + ineqn1
          xadj   = adjncy + inadj
          temp   = xadj   + ineqn1
c.debug
c     write(6,'("ncomp , neqns , nzcomp, nadj   = ", 4i8)')
c    1            ncomp , neqns , nzcomp, nadj   
c     write(6,'("itemp2, itemp3, itemp4, itemp5 = ", 4i8)')
c    1            itemp2, itemp3, itemp4, itemp5 
c     write(6,'("acolls                         = ", 4i8)')
c    1            acolls                         
c     write(6,'("xrowls, rowlst, cmpmap, colstr = ", 4i8)')
c    1            xrowls, rowlst, cmpmap, colstr 
c     write(6,'("adjncy, xadj  , temp           = ", 4i8)')
c    1            adjncy, xadj  , temp           
c.debug

          wkreqd = temp + ineqns - 1

          if ( wkreqd .gt. lwork ) go to 8100

          mxused = max ( mxused, wkreqd )
          inuse  = wkreqd
c.debug
c         call xislp3 ( 'before data move - itemp2',
c    1                  ncomp+1, work(itemp2), 6 )
c         call xislp3 ( 'before data move - itemp3',
c    1                  neqns  , work(itemp3), 6 )
c         call xislp3 ( 'before data move - itemp4',
c    1                  neqns  , work(itemp4), 6 )
c         call xislp3 ( 'before data move - itemp5',
c    1                  neqns+1, work(itemp5), 6 )
c         call xislp3 ( 'before data move - acolls',
c    1                  nzcomp , work(acolls), 6 )
c.debug

          call xislmv ( ncomp+1, work, xdslil(itemp2-1)+1, 
     1                                 xdslil(xrowls-1)+1)
          call xislmv ( neqns  , work, xdslil(itemp3-1)+1, 
     1                                 xdslil(rowlst-1)+1)
          call xislmv ( neqns  , work, xdslil(itemp4-1)+1, 
     1                                 xdslil(cmpmap-1)+1)
          call xislmv ( neqns+1, work, xdslil(itemp5-1)+1, 
     1                                 xdslil(colstr-1)+1)
          call xislmv ( nzcomp , work, xdslil(acolls-1)+1, 
     1                                 xdslil(adjncy-1)+1)

c.debug
c         call xislp3 ( 'after  data move - xrowls',
c    1                  ncomp+1, work(xrowls), 6 )
c         call xislp3 ( 'after  data move - rowlst',
c    1                  neqns  , work(rowlst), 6 )
c         call xislp3 ( 'after  data move - cmpmap',
c    1                  neqns  , work(cmpmap), 6 )
c         call xislp3 ( 'after  data move - colstr',
c    1                  neqns+1, work(colstr), 6 )
c         call xislp3 ( 'after  data move - adjncy',
c    1                  nzcomp , work(adjncy), 6 )
c.debug
c         -----------------------------------------------
c         ... expand lower adjacency (current contents of 
c             (colstr,adjncy)) to full adjacency
c         -----------------------------------------------
c.debug
c     t1 = work ( qinptm )
c     w1 = work ( qinpwl )
c     call xdslt2 ( t1, w1, t3, w3 )
c     write(6,'(15x,"before   xisli7        - accum. time = ", 
c    1         f15.6)') t3
c.debug

          call xisli7 ( neqns, neqns+1, work(colstr),
     1                  work(adjncy), work(xadj), work(temp), 
     2                  maxrow )
c.debug
c     write(6,'("after xisli7")')
c     call xislp3 ( 'compressed adjacency pointer array',
c    1              neqns+1, work(xadj), 6 )
c     call xislp3 ( 'compressed adjacency array',
c    1              nzcomp, work(adjncy), 6 )
c.debug
c.debug
c     t1 = work ( qinptm )
c     w1 = work ( qinpwl )
c     call xdslt2 ( t1, w1, t3, w3 )
c     write(6,'(15x,"after    xisli7        - accum. time = ", 
c    1         f15.6)') t3
c.debug

c         -----------------------------------------------
c         ... move contents of work(xadj) to work(colstr)
c         -----------------------------------------------

          call xislmv ( neqns+1, work, xdslil(xadj  -1)+1, 
     1                                 xdslil(colstr-1)+1)

          inuse  = xadj - 1
          xadj   = colstr

c.debug
c         write(6,'("after xisli7 and data move")')
c         call xislp3 ( 'compressed adjacency pointer array',
c    1                  neqns+1, work(xadj), 6 )
c         call xislp3 ( 'compressed adjacency array',
c    1                  nzcomp, work(adjncy), 6 )
c.debug
 
      else
 
          go to 8000
 
      end if
 
c     ---------------------
c     ... finish processing
c     ---------------------
 
      t1 = work ( qinptm )
      w1 = work ( qinpwl )
 
      call xdslt2 ( t1, w1, t3, w3 )
 
      if ( msglvl .eq. 3 ) then
          call xislp2 ( 'compressed adjacency pointer array',
     1                  neqns+1, work(xadj), output )
      end if
 
      if ( msglvl .ge. 4 ) then
          call xislp1 ( 'compressed adj. structure of org. matrix',
     1                 neqns, work(xadj), work(adjncy), output )
      end if
 
c     ----------------------------------
c     ... print optional timing messages
c     ----------------------------------
 
      if ( msglvl .ge. 1 ) then 
          write ( output, 80000 ) t2, w2, t3, w3
      end if
 
c     ------------------------------------------------------------
c     ... compute workspace required for next stage
c
c         here we use the knowledge that the number of fundamental
c         supernodes (nsup1 in xdslor) must be less than the 
c         number of compressed nodes, ncomp.  So we use ncomp
c         as an estimate for nsup1.  
c
c         If nsup1 .le. ncomp .le. 3 * neqns / 8 then temporary
c         storage for ordering is dominated by the itemp1
c         term.
c     ------------------------------------------------------------

      itemp1 =      9 * neqns +  2 * ncomp + 1
      itemp2 =      7 * neqns +  8 * ncomp + 2
      itemp3 =      6 * neqns + 10 * ncomp + 2

      temp = max ( itemp1, itemp2, itemp3 )

c.debug
c     write(6,'("storage comp. for temp. storage in ordering")')
c     write(6,'("old, itemp1, itemp2, itemp3, temp = ", / 5i12)')
c    1   16*neqns, itemp1, itemp2, itemp3, temp
c.debug

      inuse = adjncy + xdslni ( 2 * nzcomp ) - 1

      needs = lncomm 
     1      + xdslni ( 2*nzcomp + ncomp + 3*neqns + temp + 40 )
c.debug
c     write(6,'("in xdslif - needs = ", i8)') needs
c.debug
 
      sqfile = work ( qsqfl1 )
      if ( sqfile .lt. 0 ) needs = needs 
     1                           + xdslni ( 2*nzcomp + neqns + 3 )
c.debug
c     write(6,'("in xdslif - adjusting for sqfil1 - needs = ", i8)')
c    1                                               needs
c     write(6,'("neqns, nzcomp, sqfile = ", 3i8)')
c    1            neqns, nzcomp, sqfile 
c.debug
 
c     ---------------------------------------------
c     ... store information into.CMNication area
c     ---------------------------------------------

c.debug
c         call xislp3 ( 'compressed adjacency pointer array-2',
c    1                  neqns+1, work(xadj), 6 )
c         call xislp3 ( 'compressed adjacency array-2',
c    1                  2*nzcomp, work(adjncy), 6 )
c.debug
 
      work ( qstage ) = 10
      work ( qinuse ) = inuse 
      work ( qneeds ) = needs
      work ( qmxuse ) = mxused
      work ( qmsglv ) = msglvl
      work ( qoutpu ) = output
      work ( qtemp1 ) = 0.
 
      work ( qneqns ) = neqns
      work ( qnzla  ) = nzla
      work ( qnzlb )  = min ( work ( qnzlb ), work ( qnzla ) ) 
      work ( qncomp ) = ncomp
      work ( qnzcmp ) = nzcomp
 
      work ( qxadj  ) = xadj
      work ( qadjnc ) = adjncy

c.debug
c     call xislp3 ( 'compressed adjacency pointer array-end',
c    1              neqns+1, work(xadj), 6 )
c     call xislp3 ( 'compressed adjacency array-end',
c    1              2*nzcomp, work(adjncy), 6 )
c.debug

      work ( qxrwls ) = xrowls
      work ( qrwlst ) = rowlst
      work ( qcmpmp ) = cmpmap
 
      work ( qinptm ) = t3
      work ( qinpwl ) = w3
 
      go to 9000
 
c---------------------------------------------------------------------
 
c     --------------
c     ... error trap
c     --------------
 
c     ---------------------------------
c     ... incorrect processing sequence
c     ---------------------------------
 
 8000 continue
      error = -100
      call hherr ( 3, 'xdslif', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88000 ) error, stage
      go to 8999
 
c     ---------------------------------------
c     ... insufficient storage for this stage
c     ---------------------------------------
 
 8100 continue
      error = -101
      call hherr ( 2, 'xdslif', error, wkreqd )
      if ( msglvl .gt. 0 ) write ( output, 88100 ) error, wkreqd, lwork
      needs = wkreqd
      go to 8999

c     ----------------------------
c     ... i/o error on file wafil1
c     ----------------------------
 
 8900 continue
      error = -110
      call hherr ( 3, 'xdslif', error, 0 )
      if ( msglvl .gt. 0 ) write ( output, 88900 ) error, wafil1, wafil2
      go to 8999

c     --------------------------------------
c     ... set stage to prevent further calls
c     --------------------------------------
 
 8999 continue
      work(qstage) = -1
 
c---------------------------------------------------------------------
 
c     ------------------------
c     ... end of module xdslif
c     ------------------------
 
 9000 continue
      return
 
c---------------------------------------------------------------------
 
c     -----------
c     ... formats
c     -----------
 
80000 format ( /5x, 'cpu  time for structural input          = ', f15.6
     1         /5x, 'wall time for structural input          = ', f15.6 
     2         /5x, 'cpu  time for struct. inp. processing   = ', f15.6 
     3         /5x, 'wall time for struct. inp. processing   = ', f15.6)
 
81000 format ( /5x, 'number of rows in matrix                = ', i15
     1         /5x, 'number of nonzeroes in full matrix      = ', i15
     2         /5x, 'number of compressed nodes in matrix    = ', i15
     3         /5x, 'number of compressed nonzeroes          = ', i15 )
 
88000 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslif executed in an'
     2         /5x, 'incorrect sequence.  current stage = ', i10,
     3          5x, 'should be either 1 or 2.')
 
88100 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslif requires ', i15
     2         /5x, 'words of workspace and has only ', i15,
     3              ' available.' )

88900 format ( /5x, '*** fatal error no. ', i5, ' *** subroutine ',
     1              'xdslif encountered i/o error'
     2         /5x, 'on i/o file no. ', i15,
     3          5x, 'or on i/o file no. ', i15 )
 
c---------------------------------------------------------------------
 
      end
