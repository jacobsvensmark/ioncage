MODULE io_general_module

 IMPLICIT NONE

 CONTAINS

 SUBROUTINE file_lines(filename, NL)
 !// Given a filename and an integer NL, this program
 !// counts the number of lines the file and stores in NL. // JSV
  IMPLICIT NONE
  character(len=300) :: filename
  integer            :: NL, nlines, lunit
  nlines = 0
  lunit = 22
  OPEN (lunit, file = filename)
     DO
        READ (lunit,*, END=10)
        nlines = nlines + 1
     END DO
  10 CLOSE (lunit)
  NL = nlines
 END SUBROUTINE file_lines

 INTEGER FUNCTION findnl(s)
 !// Locates position of newline in string
 IMPLICIT NONE
     character(len=*) :: s
     integer :: i
     findnl = len(s)+1
     do i = 1, len(s)
       if (s(i:i) .eq. achar(10)) then
         findnl = i
         return
       end if
     end do
 END FUNCTION

END MODULE io_general_module
