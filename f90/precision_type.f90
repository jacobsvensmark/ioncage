MODULE precision_type
   IMPLICIT NONE
!  integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: dp = SELECTED_REAL_KIND (15,300)
END MODULE precision_type
