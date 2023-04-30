MODULE conversions

USE rrflist
USE rrfdata
USE arithmetic
IMPLICIT NONE

INTERFACE number
  MODULE PROCEDURE number4, number_long
END INTERFACE

PRIVATE

PUBLIC fromch, number, display, to4i, imult, toreal, char4i, putch,    &
    tochar, char

CONTAINS

!-----------------------------------------------------------------------

!  Conversion from character representations

SUBROUTINE fromch(m,l, np, k)

!  Read an RRF from a character buffer M, starting at position L+1.
!  The number must be in the form

!  sign  x1 x2 ...  xNP (.) n(!)(^exp) (.) n(!)(^exp) ...

!  where sign is +, (+)i, -, or -i, exp is (+)n(/2) or -n(/2),
!  and n is an unsigned integer.
!  In all the above expressions, parentheses denote optional items and
!  do not appear in the actual input. Spaces are required as shown, but
!  two or more spaces may be replaced by a single space.
!  x1 ... xNP are the exponents of the first NP primes. (NP may be zero)
!  Subsequent integer and factorial terms may occur in any order and may
!  be interspersed with each other.
!  The ^exp may be omitted if exp is unity; the ^ may be omitted if
!  the next character is + or -.

USE factorials

CHARACTER(LEN=*), INTENT(IN) :: m
INTEGER, INTENT(INOUT) :: l
INTEGER, INTENT(IN) :: np
TYPE(node), POINTER :: k

INTEGER(KIND=long) :: nn
INTEGER :: n, x, p
LOGICAL :: sign
CHARACTER :: c
TYPE(node), POINTER :: kz=>null()

call unit(k)
sign=.false.

do
  if (m(l:) == "") then
    if (sign) then
      return
    else
      call report('Empty power-of-primes expression')
      call clear(k)
      return
    endif
  end if
  l=l+1

  c=m(l:l)
  ! print "(3a)", '"', c, '"'
  select case(c)
  case(space)
!  Space
    if (sign) exit
  case(plus,minus)
    if (sign) then
      call report('Syntax error in input')
      call clear(k)
      return
    endif
    if (c == minus) then
      k%x=k%x+2
    endif
    sign=.true.
  case(uci,lci)
!  I or i
    k%x=k%x+1
    sign=.true.
!   exit
  case(digit(0))
!  0
    call clear(k)
    return
!  case(digit(1))
!  1
!    if (m(l+1:l+1) .eq. space .or. m(l+1:l+1) .eq. dot) then
!      l=l+1
!    endif
!    exit
  case default
    call report('Sign factor missing in power-of-primes expression')
    exit
  end select
end do

!  Exponents of first NP primes
do p=1,np
  x=2
  call number(m,l, nn,x)
  if (x .eq. 0) exit
  if (n .eq. 0) cycle
  call link(k)
  k%prime=prime(p)
  k%x=(nn+nn)/x
end do
k=>k%next

!  Explicit n^exponent or n!^exponent
! l=l-1
do
  if (m(l:) == "") exit
  do while (m(l:l) .eq. space .or. m(l:l) .eq. dot)
    l=l+1
  end do
  x=1
  call number (m,l, nn,x)
  if (x .eq. 0) return
  if (m(l:l) .eq. shriek) then
    n=nn
    call setf(n,kz)
    l=l+1
  else
    call fromi(nn,kz)
  endif
  select case(m(l:l))
  case(arrow,plus,minus)
    if (m(l:l) == arrow) l=l+1
    x=2
    call number(m,l, nn,x)
    if (x .eq. 0) call report('Syntax error in power-of-primes expression')
    if (x .eq. 2) call root(kz)
    n=nn
    call mult(k,kz,n)
  case(space)
    call mult(k,kz,1)
  case default
    call report('Syntax error in power-of-primes expression')
  end select
end do

call clear(kz)

END SUBROUTINE fromch

!-----------------------------------------------------------------------

SUBROUTINE number_long(m,l, n,x)
!  If X=1 on entry, read an integer (-)nnn... terminated by any
!  non-digit;  if X has any other value on entry, read an integer or
!  half-integer (-)nnn...(/2).
!  In either case, the number is read from the buffer M, starting at
!  position L and terminated by space, non-numerical character or end of
!  buffer. L is updated to point at the next character after the number.
!  If no digit is found before a terminator, or if the digit is missing
!  after /, X is set to zero; otherwise X is set to m if /m is
!  present and expected and to 1 otherwise.

CHARACTER(LEN=*), INTENT(IN) :: m
INTEGER(KIND=long), INTENT(OUT) :: n
INTEGER, INTENT(INOUT) :: l, x

INTEGER :: i, sign
LOGICAL :: found, whole

whole=x .eq. 1
n=0
x=1
sign=1
found=.false.


if (debug) print "(a,a,i3/a,i2)", trim(m), "  length ", len(m), "Number: l=", l
l=l-1
chars: do
  l=l+1
  if (debug) print "(i2,1x,a,1x,i0)", l, m(l:l), n
  ! if (l .gt. len(trim(m)) .and. found) exit

  do i=0,9
    if (m(l:l) .eq. digit(i)) then
      !  Digit found
      if (x .eq. 0) then
        x=i
        l=l+1
        exit
      else if (n .gt. max_long_int/10) then
        call report('Integer overflow')
        x=0
        exit
      else
        found=.true.
        n=n*10+i
        cycle chars
      endif
    end if
  end do

  !  Non-digit
  select case (m(l:l))
  case(space)
    if (x .eq. 0) call report('Denominator missing')
    if (found) exit
  case(plus)
    if (found) exit
    sign=1
  case(minus)
    if (found) exit
    sign=-1
  case(slash)
    if (.not. found) then
      call report('Numerator missing')
      x=0
      exit
    else if (whole)then
      exit
    else
      x=0
      if (l .ge. len(trim(m))) then
        call report('Denominator missing')
        exit
      end if
      cycle
    end if
  case default
    if (found) then
      exit
    else
      call report('Non-digit while attempting to read number')
      x=0
      exit
    endif
  end select
end do chars

n=n*sign
if (debug) print "(a,i0)", "n = ", n

END SUBROUTINE number_long

!-----------------------------------------------------------------------

SUBROUTINE number4(m,l, n,x)

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: m
INTEGER, INTENT(OUT) :: n
INTEGER, INTENT(INOUT) :: l, x
INTEGER(KIND=long) :: nn

call number_long(m,l, nn,x)
n=int(nn)

END SUBROUTINE number4

!-----------------------------------------------------------------------

!  Conversion to other representations

SUBROUTINE display(k,mode)

!  Display: mode = 1: power-of prime
!                  2: real
!                  3: (i1/i2)*sqrt(i3/i4)

TYPE(node), POINTER :: k
INTEGER, INTENT(IN) :: mode

CHARACTER :: c
CHARACTER(LEN=1024) :: iout
REAL(KIND=qp) :: a
INTEGER :: j, m

select case(mode)
case(1)
  m=1
  call tochar(k, iout, m, 0)
  if (.not. error) write (lout,"(a)") iout(1:m)
  
case(2)
  call toreal(k,a,j)
  if (nonzero(k) .and. a==0.0_qp) then
    write(lout,"(8x,a)") "Real value too big"
  else
    if (j .eq. 1) then
      c = lci
    else
      c = space
    endif
    write(lout,"(g30.16,1x,a1)") a, c
  end if
  
case(3)
  m=1
  call char4i(k, iout, m,132)
  if (.not. error) write (lout,"(a)") iout(1:m)
end select

error=.false.

END SUBROUTINE display

!-----------------------------------------------------------------------

SUBROUTINE to4i(k0,i1,i2,i3,i4)
!  Express PP K0 as (I1/I2)sqrt(I3/I4).  I2, I4 >0; I1, I3 may be
!  negative.
!  If integer overflow occurs, I2 is set to zero and I1 and I3 contain
!  rubbish.

TYPE(node), POINTER :: k0
INTEGER(KIND=long), INTENT(OUT) :: i1,i2,i3,i4

TYPE(node), POINTER :: k
INTEGER :: n
LOGICAL :: ovflow=.false.

k=>k0
i1=1
i2=1
i3=1
i4=1
if (iszero(k)) then
  i1=0
  return
end if

!  Sign
k%x=mod(k%x,4)
if (k%x .lt. 0) k%x=k%x+4
if (mod(k%x,2) .ne. 0) i3=-1
if (k%x .ge. 2) i1=-1

ovflow=.false.
do
  k=>k%next
  if (k%prime .lt. 0) exit
  n=k%x
  if (n .gt. 0) then
    if (mod(n,2) .ne. 0) then
      call imult(i3,k%prime,1)
      if (i3 .eq. 0) ovflow=.true.
      n=n-1
    endif
    n=n/2
    call imult (i1,k%prime,n)
    if (i1 .eq. 0) ovflow=.true.
  else
    n=-n
    if (mod(n,2) .ne. 0) then
      call imult(i4,k%prime,1)
      if (i4 .eq. 0) ovflow=.true.
      n=n-1
    endif
    n=n/2
    call imult(i2,k%prime,n)
    if (i2 .eq. 0) ovflow=.true.
  endif
end do

if (ovflow) then
  i2=0
  call report('Integer overflow')
end if

END SUBROUTINE to4i

!-----------------------------------------------------------------------

SUBROUTINE imult (i, n, l)
!  Multiply I by N**L, checking for overflow.  L >= 0.
!  If overflow occurs, I is set to zero.

INTEGER(KIND=long), INTENT(INOUT) :: i
INTEGER(KIND=long), INTENT(IN) :: n
INTEGER, INTENT(IN) :: l
INTEGER(KIND=long) :: j, m

if (l .eq. 0) return
m=max_long_int/n
do j=1,l
  if (abs(i) .gt. m) then
    i=0
    exit
  end if
  i=i*n
end do

END SUBROUTINE imult

!-----------------------------------------------------------------------

SUBROUTINE toreal (k,a,i)
!  Express K as A*i**I, where I=0 or 1 and A is real
!  If the real is going to be too big, A is set to zero and an error
!  report is generated.

TYPE(node), POINTER :: k
INTEGER, INTENT(OUT) :: i
REAL(KIND=qp), INTENT(OUT) :: a

REAL(KIND=qp) :: an

i=0
if (iszero(k)) then
  a=0d0
  return
endif

a=1d0

do
  k=>k%next
  if (k%prime .gt. 0) then
    an=k%prime
    if (k%x > 0 .and. max_real/(an**k%x) < a) then
      call report("Real overflow")
      a=0
      return
    else
      a=a*an**k%x
    end if
  else
    a=sqrt(a)
    k%x=mod(k%x,4)
    if (k%x .lt. 0) k%x=k%x+4
    if (mod(k%x,2) .ne. 0) i=1
    if (k%x .ge. 2) a=-a
    exit
  end if
end do

END SUBROUTINE toreal

!-----------------------------------------------------------------------

SUBROUTINE char4i(k, m, m1,m2)
!  Convert K to 4-integer form (I1/I2)*sqrt(I3/I4) as characters
!  in the string M, starting at M1.
!  The pointer M1 is incremented to point at the next available
!  position in the string.

TYPE(node), POINTER :: k
INTEGER, INTENT(IN) ::  m2
CHARACTER(LEN=*), INTENT(INOUT) :: m
INTEGER, INTENT(INOUT) :: m1

INTEGER(KIND=long) :: i1,i2,i3,i4, modi1
LOGICAL :: brackt

call to4i (k, i1,i2,i3,i4)
!  Overflow?
if (error) return

!  Clear buffer
m(m1:m2) = ""

!  Check for zero
if (i1 .eq. 0 .or. i3 .eq. 0) then
  call putch(m, m1,m2, "0")
  return
endif

!  Attach sign
if (i1 .lt. 0) call putch(m, m1,m2, minus)
modi1=abs(i1)
!  Check for anything before square root
if (modi1 .ne. 1 .or. i2 .ne. 1) then
  !  Brackets only if I2, and I3 or I4, not equal to 1
  brackt=(i2 .ne. 1 .and. (i3 .ne. 1 .or. i4 .ne. 1))
  if (brackt) call putch(m, m1,m2, bra)
  call char(modi1,1, m, m1,m2, -1)
  !  Put in slash if I2 not equal to 1
  if(i2.ne.1) then
    call putch(m, m1,m2, slash)
    call char(i2,1, m,m1,m2, -1)
    if (brackt) call putch(m, m1,m2, ket)
  endif
  if (i3 .eq. 1 .and. i4 .eq. 1) return
  call putch(m,m1,m2, star)
endif
if (i3 .ne. 1 .or. i4 .ne. 1) then
  m(m1:m1+4)="sqrt("
  m1=m1+5
  call char(i3,1, m,m1,m2, -1)
  if (i4 .ne. 1) then
    call putch(m, m1,m2, slash)
    call char(i4,1, m, m1,m2, -1)
  endif
  call putch(m,m1,m2, ")")
  return
endif
call putch(m, m1,m2, digit(1))

END SUBROUTINE char4i

!-----------------------------------------------------------------------

SUBROUTINE putch(m, m1,m2, ch)

CHARACTER(LEN=*), INTENT(INOUT) :: m
CHARACTER(LEN=1), INTENT(IN) :: ch
INTEGER, INTENT(INOUT) :: m1
INTEGER, INTENT(IN) :: m2

if (m1 > m2) call report('Buffer overflow')
m(m1:m1)=ch
m1=m1+1

END SUBROUTINE putch

!-----------------------------------------------------------------------

SUBROUTINE tochar(k, m, m1, mp)
!  Converts K to character form in array M(MAX), starting at position
!  M1.  M1 is updated to point at the next available position
!  in the buffer.  The exponents only are output for the first
!  MP primes; subsequent entries are output as  . prime^exponent.

TYPE(node), POINTER :: k
CHARACTER(LEN=*), INTENT(INOUT) :: m
INTEGER, INTENT(INOUT) :: m1
INTEGER, INTENT(IN) :: mp

INTEGER :: i, p, ml

ml=len(m)
if (ml .lt. mp*3+2) call report('Buffer overflow')
do i=m1,ml
  m(i:i)=space
end do
if (iszero(k)) then
  m(m1+1:m1+1)=digit(0)
  m1=m1+2
  return
endif

!  Sign
k%x=mod(k%x,4)
if (k%x .lt. 0) k%x=k%x+4
m(m1:m1)=plus
if (k%x .ge. 2) m(m1:m1)=minus
i=mod(k%x,2)
if (i .eq. 1) m(m1+1:m1+1)=lci
if (associated(k%next,k)) then  ! Sign factor only
  if (i .eq. 0) m(m1+1:m1+1)=digit(1)
  m1=m1+3
  return
end if
k=>k%next
m1=m1+3

!  Powers of first MP primes
p=1
do while (p .le. mp)
  if (prime(p) .eq. k%prime) then
    call char(int(k%x,long), 2, m, m1, m1+4, 1)
    k=>k%next
    if (k%prime .lt. 0) return
  else
    call char(0_long, 2, m, m1, m1+5, 1)
  endif
  p=p+1
end do

do
  if (ml .lt. m1+5) call report('Buffer overflow')

  call char(k%prime,1,m,m1,ml,-1)
  if (k%x .ne. 2) then
    m(m1:m1)=arrow
    if (ml .lt. m1+2) call report('Buffer overflow')
    m1=m1+1
    call char(int(k%x,long), 2, m, m1, ml, -1)
  end if
  k=>k%next
  if (k%prime .lt. 0) exit
  m(m1+1:m1+1)=dot
  m1=m1+3
end do

END SUBROUTINE tochar

!-----------------------------------------------------------------------

SUBROUTINE char(n, x, m, m1, m2, jr)
!  Converts N/X (N, X integers) into characters in string M, starting at
!  M1 and finishing at M2, left or right justified according
!  as JR <0 or >0.  M1 is reset on exit to the next available
!  element of M.  M is assumed to be cleared to spaces on entry.
!  The value of X must be in the range 1 to 9 (it will normally be
!  1 or 2)

INTEGER(KIND=long), INTENT(IN) :: n
INTEGER, INTENT(IN) :: x, m2, jr
INTEGER, INTENT(INOUT) :: m1
CHARACTER(LEN=*), INTENT(INOUT) :: m

INTEGER :: i, j, ns, l, md
INTEGER(KIND=long) :: nn

if (m1 .gt. len(m)) call report('Buffer overflow')
ns=m1
nn=n
l=m2+1
md=0
if (n .lt. 0) then
  nn=-n
  m(m1:m1)=minus
  ns=m1+1
end if

if (x .eq. 2) then
  if (mod(nn,2_long) .eq. 0) then
    nn=nn/2_long
    md=2
  else
    m(m2:m2)=digit(x)
    m(m2-1:m2-1)=slash
  end if
  l=l-2
end if

do
  j=mod(nn, 10_long)
  l=l-1
  if (l .lt. ns) then
    do i=m1,m2
      m(i:i)=star
    end do
    m1=m2+1
    return
  end if
  m(l:l)=digit(j)
  nn=nn/10_long
  if (nn .eq. 0) exit
end do

if (l .eq. ns) then
  m1=m2+1
else if (jr .gt. 0) then
  !  Right justify (move sign)
  m(l-1:l-1)=m(m1:m1)
  m(m1:m1)=space
  m1=m2+1
else
  !  Left justify
  do i=l,m2
    m(ns:ns)=m(i:i)
    ns=ns+1
  end do
  l=max0(ns,l)
  do i=l,m2
    m(i:i)=space
  end do
  m1=ns-md
end if

END SUBROUTINE char

END MODULE conversions
