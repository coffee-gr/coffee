SUBROUTINE start(quiet,primes)
USE rrfdata
USE rrflist
USE factorials

IMPLICIT NONE
INTEGER, INTENT(IN), OPTIONAL :: primes
LOGICAL, INTENT(IN), OPTIONAL :: quiet
INTEGER(KIND=long) :: m, i, j, ok
LOGICAL :: verbose=.true.

if (present(quiet)) then
  verbose=.not. quiet
else
  verbose=.true.
end if
if (present(primes)) then 
  nprime=primes
else
  nprime=10000
endif
!  Hard errors
hard_errors=.true.

print "(20x,a/27x,a/20x,2a,3x,a)", "Root rational fraction program",   &
    "by Anthony Stone",                                                &
    "Version ", trim(version), trim(date)
if (verbose) then
  print "(/a)",    "Constructing table of primes ..."
end if

!  Set up table of first NPRIME primes
allocate(prime(nprime),stat=ok)
if (ok>0) then
  print "(a,i0,a)", "Can't allocate ", nprime, " primes"
end if
prime(1)=2
m=1
j=1
do
  m=m+2
  do i=1,j
    if (mod(m,prime(i)) .eq. 0) exit    ! not prime
    if (m .lt. prime(i)*prime(i)) then  ! prime found
      j=j+1
      prime(j)=m
      ! if (j .le. 50) print "(i5,i10)", j, prime(j)
      exit
    endif
  end do
  if (j .eq. nprime) exit
end do

!  Initialise factorial table (initially contains just 0! and 1!).
call unit(factorial(0)%rrf)
call unit(factorial(1)%rrf)

if (verbose) then
  print "(a,i0)",                                                      &
      "Prime table includes primes up to ", prime(NPRIME)
end if

END SUBROUTINE start
