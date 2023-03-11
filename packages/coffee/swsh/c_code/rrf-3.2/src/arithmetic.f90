MODULE arithmetic

USE rrflist
USE rrfdata

IMPLICIT NONE

PRIVATE

PUBLIC mult, power, add, accum, subtr, root, fromi, fromr, multi

CONTAINS

!-----------------------------------------------------------------------

SUBROUTINE mult(k1,k2,mexp)
!  Multiply k1 by k2**mexp. mexp is optional, default value 1.

TYPE(node), POINTER :: k1, k2
INTEGER, INTENT(IN), OPTIONAL :: mexp

TYPE(node), POINTER :: kq1
INTEGER :: mx
CHARACTER(LEN=80) :: string

if (present(mexp)) then
  mx=mexp
else
  mx=1
end if

!  Special cases
if (iszero(k1) .or. mx .eq. 0) return
if (associated(k1,k2)) then
  string="Invalid arguments to MULT"
  call report(string)
endif
if (iszero(k2)) then
  if (mx < 0) then
    call report('Division by zero')
    call clear(k1)
  else
    call zero(k1)
  endif
  return
endif

!  General case
k2=>k2%next
kq1=>k1
k1=>k1%next
do
  if (k1%prime .eq. k2%prime) then
!  Term present in both K1 and K2: form product term
    k1%x=k1%x+mx*k2%x
    if (k1%prime .lt. 0) exit
    if (k1%x .eq. 0) then
!  Product term has zero exponent -- remove it and return link to free list
      kq1%next=>k1%next
      k1%next=>free
      free=>k1
      k1=>kq1%next
    else
      kq1=>k1
      k1=>k1%next
    endif
    k2=>k2%next
  else if (k2%prime .gt. 0 .and.                                      &
          (k1%prime .lt. 0 .or. k1%prime .gt. k2%prime)) then
!  Term in K2 not in K1: insert it
    call link(kq1)
    kq1%prime=k2%prime
    kq1%x=mx*k2%x
    k2=>k2%next
    cycle
  else
!  Term in K1 not in K2: advance to next term
    kq1=>k1
    k1=>k1%next
  endif
end do

END SUBROUTINE mult

!-----------------------------------------------------------------------

SUBROUTINE power(k, m)
!  Replace K by K**M

TYPE(node), POINTER :: k
INTEGER, INTENT(IN) :: m

if (m .eq. 0) then
  call unit(k)
else if (iszero(k)) then
  if (m < 0) call report('Division by zero')
else
  do
    k=>k%next
    k%x=m*k%x
    if (k%prime < 0) exit
  end do
endif

END SUBROUTINE power

!-----------------------------------------------------------------------

SUBROUTINE add(k1,k2)
!  Add K2 to K1

TYPE(node), POINTER :: k1, k2
INTEGER(KIND=long) :: a1

a1=1
call accum(k1,a1, k2)
if (a1 == 1) return
call report('Integer overflow')
call clear(k1)

END SUBROUTINE add

!-----------------------------------------------------------------------

SUBROUTINE accum(k1,a1, k2)
!  Add K2 to K1*A1.
!  Normally A1 will contain 1 on entry and exit.  If however it proves
!  impossible to decompose the sum into primes small enough to occupy
!  one integer location, the oversize factor is returned in A1.  Thus
!  the sum is K1*A1.  Subsequent additions can start from this point;
!  at the end of a sequence of additions the result can usually be
!  factorised satisfactorily.

TYPE(node), POINTER :: k1, k2
INTEGER(KIND=long), INTENT(INOUT) :: a1

LOGICAL :: next1, next2
TYPE(node), POINTER :: k=>null(), k10=>null(), k20=>null()
INTEGER :: ld
INTEGER(KIND=long) :: a2, p

!  Special cases
if (iszero(k2)) then
!  K2=0:  force factorisation of A1 if possible
  if (a1 == 1) return

  call fromi(a1,k)
  call mult(k1,k,1)
  call clear(k)
  return
else if (iszero(k1)) then
  call copy (k2, k1)
  a1=1
  return
else if (associated(k1,k2)) then
  call report('Invalid arguments to ACCUM')
end if

call unit(k)
a2=1
next1=.true.
next2=.true.
k10=>k1
k20=>k2

!  The next section constructs the highest common factor of K1 and
!  K2 and stores it as K.  The multipliers are A1 and A2, so that
!  K1=A1*K and K2=A2*K.

do
  if (next1) k1=>k1%next
  if (next2) k2=>k2%next
  if (k1%prime == k2%prime) then
    !  Entry in both K1 and K2
    k%prime=k1%prime
    k%x=min0(k1%x,k2%x)
    ld=k2%x-k1%x
    next1=.true.
    next2=.true.
  else if ((0 < k2%prime .and. k2%prime < k1%prime) .or. k1%prime .lt. 0) then
    !  Entry in K2 only
    k%prime=k2%prime
    k%x=min0(0,k2%x)
    ld=k2%x
    next2=.true.
    next1=.false.
  else if ((0 < k1%prime .and. k1%prime < k2%prime) .or. k2%prime .lt. 0) then
    !  Entry in K1 only
    k%prime=k1%prime
    k%x=min0(k1%x,0)
    ld=-k1%x
    next2=.false.
    next1=.true.
  end if
  ! print "(3i3,2f8.0)", k%prime, k%x, ld, a1, a2

  if (ld .ne. 0) then
    if (mod(ld,2) .ne. 0) then
      call report('Incompatible summands')
      call list(k10)
      call list(k20)
      stop
    endif
    ld=ld/2
    p=k%prime
    if (ld .gt. 0) then
      if (p<0 .or. a2 < max_long_int/p**ld) then
        a2=a2*p**ld
      else
        print "(a,i0,a,i0,a,i0)", "a2 = ", a2, "  p = ", p, "  ld = ", ld
        call report("Overflow in ADD")
      end if
      else
        if (p<0 .or. a1 < max_long_int/p**(-ld)) then
        a1=a1*p**(-ld)
      else
        print "(a,i0,a,i0,a,i0)", "a1 = ", a1, "  p = ", p, "  ld = ", ld
        call report("Overflow in ADD")
      end if
    end if
  end if
  if (k%prime < 0) then
    exit
  else
    !  Get new link for next entry -- but if current link has zero
    !  exponent it is re-used
    if (k%x .ne. 0) call link(k)
  end if

!  End of decomposition into h.c.f. and multipliers
end do

! print "(a,f10.0,a,f10.0)", "accum: a1 = ", a1, ", a2 = ", a2
call clear(k1)
if (real(abs(a1+a2),qp) > real(max_long_int,qp) ) then
  call report('overflow in add')
  call clear(k)
  call clear(k1)
  k2=>k20
endif
a1=a1+a2
call fromi(a1,k1)
! if (a1 > 1)                                                       &
!     print '(a,i0)', 'Intermediate value not factorized: ', A1
a1=1
call mult(k1,k,1)
call clear(k)

END SUBROUTINE accum

!-----------------------------------------------------------------------

SUBROUTINE subtr(k1,k2)
!  Subtract K2 from K1

TYPE(node), POINTER :: k1, k2

call chsign(k2)
call add(k1,k2)
call chsign(k2)

END SUBROUTINE subtr

!-----------------------------------------------------------------------

SUBROUTINE root(k)

TYPE(node), POINTER :: k

LOGICAL :: ok

if (iszero(k)) return

k%x=modulo(k%x,4)
ok=.true.
do
  k=>k%next
  if (mod(k%x,2) .eq. 0) then
    k%x=k%x/2
  else
    ok=.false.
  endif
  if (k%prime .lt. 0) exit
end do

if (.not. ok) then
  call report('Square root argument invalid')
  ! call clear(k)
end if

END SUBROUTINE root

!-----------------------------------------------------------------------

SUBROUTINE fromi(n,k)
!  Express integer N as PP K

INTEGER(KIND=long), INTENT(IN) :: n
TYPE(node), POINTER :: k

INTEGER(KIND=long) :: m
INTEGER :: p
LOGICAL :: new

call zero(k)
if (n .eq. 0) return

if (debug) print "(a,i0)", "Factorising ", n
m=n
call unit(k)
if (n .lt. 0) then
  m=-m
  k%x=2
endif
if (m .eq. 1) return

if (debug) print "(2i5,3x,i0)", k%prime, k%x, m
do p=1,nprime
  new=.true.
  do while (modulo(m,prime(p)) .eq. 0)
    if (new) then
      call link(k)
      k%x=0
      k%prime=prime(p)
      new=.false.
    end if
    k%x=k%x+2
    m=m/prime(p)
  end do
  if (debug .and. .not. new .and. k%x>0) print "(2i5,3x,i0)", k%prime, k%x, m
  if (m .eq. 1) exit
  if (m < prime(p)*prime(p)) then
    !  Remaining factor must be prime
    call link(k)
    k%x=2
    k%prime=m
    if (debug) print "(2i5,3x,i0)", k%prime, k%x, m
    m=1
    exit
  end if
end do
if (m>1) then
  !  Treat remaining factor as prime
  call link(k)
  k%x=2
  k%prime=m
  if (debug) print "(2i5,3x,i0)", k%prime, k%x, m
end if

k=>k%next
if (debug) call list(k)

END SUBROUTINE fromi

!-----------------------------------------------------------------------

SUBROUTINE fromr(a, k)
!  Express A (assumed integer value) as RRF K.
!  If the factor remaining after all tabulated prime factors have been
!  extracted is too large to express as an integer value, it is
!  returned in A.  Otherwise A is set to 1.  A>0 on exit in any case,
!  though it need not be positive on entry.

REAL(KIND=qp), INTENT(INOUT) :: a
TYPE(node), POINTER :: k

REAL(KIND=qp) :: pr, b
LOGICAL :: new
INTEGER :: p
INTEGER, PARAMETER ::  dp=kind(1d0)

call zero(k)
b=a
a=1d0
if (abs(b) .lt. 0.5d0) return

call unit(k)
if (b .lt. 0d0) then
  b=-b
  k%x=2
endif
if (b .gt. max_real) then
  call report('Integer value too large to store as real_qp')
  call clear(k)
  return
endif
if (b .lt. 1.5d0) return

b=floor(b,long)
do p=1,nprime
pr=prime(p)
new=.true.
do while (mod(b,pr) .eq. 0d0)
  if (new) then
    call link(k)
    k%x=0
    k%prime=prime(p)
    new=.false.
  endif
  k%x=k%x+2
  b=floor(b/pr)
end do
if (b .lt. 1.5d0) exit
if (b .lt. pr*pr) then
  !  Remaining factor assumed prime
  if (b .gt. real(max_long_int,qp)) then
    a=b
  else
    call link(k)
    k%x=2
    k%prime=b
  end if
  exit
end if
end do
k=>k%next

END SUBROUTINE fromr

!-----------------------------------------------------------------------

SUBROUTINE multi(k,i,m)
!  Multiply K by I**M, I, M integers

TYPE(node), POINTER :: k
INTEGER(KIND=long), INTENT(IN) :: i
INTEGER, INTENT(IN) :: m
TYPE(node), POINTER :: kz=>null()

call fromi(i,kz)
call mult(k,kz,m)
call clear(kz)

END SUBROUTINE multi

END MODULE arithmetic
