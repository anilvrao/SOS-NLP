      subroutine xislp1 ( title, length, pointr, ilist, msgunt )
 
c
c  purpose -- to print out a list in easily readable 80 column format
c
c  created            -- 15-jun-87, cca
c  last modifications -- 11-jul-89, rgg
c
c  input variables --
c
c      title  -- title of vector being listed
c      length -- number of partitions in the list
c      pointr -- pointer into the list
c      ilist  -- integer vector containing list items
c      msgunt -- message unit number
c
c  output variables --
c
c      none
c
c-----------------------------------------------------------------------
 
c     --------------
c     ... parameters
c     --------------
 
      character(len=*)    title
 
      integer             msgunt, length, pointr(*), ilist(*)
 
c     -------------------
c     ... local variables
c     -------------------
 
      character(len=75)   line
 
      integer             i     , j     , j1    , j2    , j3    ,
     1                    l     , n
 
c-----------------------------------------------------------------------
 
c     ---------------
c     ... write title
c     ---------------
 
      l = min ( len (title), 75 )
 
      do i = 1, l
          line(i:i) = '-'
      enddo
 
      do i = l+1, 75
          line(i:i) = ' '
      enddo
 
      write ( msgunt, 2000 ) title (1:l), line (1:l)
 
c     ------------------------
c     ... write out the vector
c     ------------------------
 
      j1 = 1
 
      do i = 1,length
 
          n  = pointr(i+1) - pointr(i)
          j3 = j1 + n - 1
          j2 = min ( j3, j1+9 )
 
          write ( msgunt, 2100 ) i, ( ilist(j), j = j1, j2 )
 
          if ( j3 .gt. j2 )
     1        write ( msgunt, 2200 ) ( ilist(j), j = j2+1, j3 )
 
          j1 = j1 + n
 
      enddo
 
      return
 
c-----------------------------------------------------------------------
 
c     -----------
c     ... formats
c     -----------
 
 2000 format ( /5x, a  /5x, a  / )
 
 2100 format ( 5x, i7, ' : ', 10i7 )
 
 2200 format ( 15x, 10i7 )
 
c-----------------------------------------------------------------------
 
      end
