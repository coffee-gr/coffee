MODULE factorials

USE rrfdata
USE rrflist
USE arithmetic
IMPLICIT NONE

PRIVATE

PUBLIC setf, multf, clearf, maxfct

INTEGER, PARAMETER :: MAXFCT=1000
TYPE(rrf), PUBLIC, SAVE :: factorial(0:MAXFCT)

CONTAINS

SUBROUTINE setf(n,k)
!  Check that factorial(N) exists, and calculate it if it does not,
!  where N=abs(I).  If argument K is present, set RRF K to
!  factorial(N).

INTEGER, INTENT(IN) :: n
TYPE(node), POINTER, OPTIONAL :: k
INTEGER :: m

if (n .eq. 0) then
  if (present(k)) call unit(k)
  return
endif

if (n .gt. maxfct) then
  call report('Reference to untabulated factorial')
  if (present(k)) call clear(k)
  return
endif

if (.not.associated(factorial(n)%rrf)) then
!  Factorial(N) is not tabulated; calculate it.
  m=n-1
  do while (.not. associated(factorial(m)%rrf))
    m=m-1
  end do
  do while (m .lt. n)
    m=m+1
    call fromi(int(m,long),factorial(m)%rrf)
    call mult(factorial(m)%rrf,factorial(m-1)%rrf,1)
  end do
!       write (lout,'(a,i3,a)') 'Factorial', N, ' tabulated'
endif
if (present(k)) call copy(factorial(n)%rrf,k)

END SUBROUTINE setf

!-----------------------------------------------------------------------

SUBROUTINE multf (k, n, mexp)
!  Multiply K by (factorial N)**MX

TYPE(node), POINTER :: k
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(IN), OPTIONAL ::  mexp

if (n .le. 1) return
if (n .gt. maxfct) then
  call report('Reference to untabulated factorial')
  call clear(k)
  return
endif
!  Check that the factorial is tabulated
call setf(n)
call mult (k, factorial(n)%rrf, mexp)

END SUBROUTINE multf

!-----------------------------------------------------------------------

SUBROUTINE clearf(n1,n2)
!  Clear the factorials from n1! to n2! inclusive

INTEGER, INTENT(IN) :: n1,n2
INTEGER :: m,m1,m2

m1=max(n1,2)
m2=min(n2,maxfct)

do m=m1,m2
  call clear( factorial(m)%rrf)
end do

END SUBROUTINE clearf

END MODULE factorials
