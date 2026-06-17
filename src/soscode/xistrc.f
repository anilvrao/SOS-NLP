      logical function xistrc ( strnga, strngb )
 
c
c     purpose
c     -------
c
c     xistrc does a case insenstive comparison of two strings to       
c     determine if they are the same modulo case.
c
c     created         05-may-00   -- rgg --
c     last modified   
c
c     input arguments
c     ---------------
c
c     strnga      c   character string a
c     strngb      c   character string b
c
c     output arguments
c     ----------------
c
c     xistrc      l   .true. is string a is equal to string b 
c                     modulo case.  .false. otherwise.
c
c---------------------------------------------------------------------
 
c     -------------------------------
c     ... global variable declaration
c     -------------------------------
 
      character(len=*)    strnga, strngb  
 
c     ------------------------------
c     ... local variable declaration
c     ------------------------------
 
      integer             i     , lena  , lenb

c     ----------------------
c     ... external functions
c     ----------------------

      logical             lsame

      external            lsame
 
c---------------------------------------------------------------------

      xistrc = .false.
     
      lena = len ( strnga )
      lenb = len ( strngb )

      if ( lena .lt. lenb ) return

      do i = 1, lenb
          if ( .not. lsame ( strnga(i:i), strngb(i:i) ) ) return
      enddo
 
c---------------------------------------------------------------------
 
c     -----------------
c     ... strings match
c     -----------------

      xistrc = .true.
 
      return
      end
