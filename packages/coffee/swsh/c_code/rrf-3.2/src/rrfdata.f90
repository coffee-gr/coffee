MODULE rrfdata

IMPLICIT NONE

PUBLIC

CHARACTER(LEN=6) :: version="3.2.02"
CHARACTER(LEN=18) :: date="18 June 2009"

!  Machine-dependent values
!  8-byte integers
INTEGER, PARAMETER :: long=SELECTED_INT_KIND(16)
!  If 8-byte integers are not available, use
!  INTEGER, PARAMETER :: long=SELECTED_INT_KIND(8)
INTEGER(KIND=long), PARAMETER :: max_long_int=huge(1_long)

#if defined(NAGF95) || defined(G95) || defined(NOR16)
!  16-byte reals not available
INTEGER, PARAMETER :: qp=SELECTED_REAL_KIND(12)
REAL(KIND=qp), PARAMETER :: max_real=huge(1.0_qp)
#else
!  16-byte reals
INTEGER, PARAMETER :: qp=SELECTED_REAL_KIND(24)
REAL(KIND=qp), PARAMETER :: max_real=huge(1.0_qp)
#endif

!  Characters
CHARACTER(LEN=1), PARAMETER ::                                         &
    arrow="^", space=" ", dot=".", plus="+", minus="-", uci="I",       &
    lci="i", slash="/", star="*", shriek="!", bra="(", ket=")"

CHARACTER(LEN=1), PARAMETER :: digit(0:10)=                            &
    (/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '0'/)

!  LIN and LOUT are the units used for READ and WRITE respectively.
INTEGER :: lin=5, lout=6

LOGICAL :: debug=.false.

END MODULE rrfdata
