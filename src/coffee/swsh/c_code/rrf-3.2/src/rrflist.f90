MODULE rrflist

USE rrfdata, ONLY : long
IMPLICIT NONE

PRIVATE

INTEGER, PUBLIC :: nprime
INTEGER(KIND=long), ALLOCATABLE, PUBLIC :: prime(:)

LOGICAL, PUBLIC :: hard_errors=.true., error=.false.

TYPE node
  INTEGER :: x
  INTEGER(KIND=long) :: prime
  TYPE(node), POINTER :: next=>null()
END TYPE node
TYPE rrf
  TYPE(node), POINTER :: rrf=>null()
END TYPE rrf

!  free points to a linked list of free pointers, or is null
TYPE(node), POINTER, PUBLIC, SAVE :: free=>null()

PUBLIC node, rrf, get, link, clear, copy, rename, unit, zero,          &
    iszero, nonzero, equal, exch, chsign, conjugate, list, report, stri

CONTAINS

SUBROUTINE get(k)
!  Get a single link and point K at it

TYPE(node), POINTER :: k
INTEGER :: ok

if (associated(free)) then
  k=>free
  free=>free%next
else
  allocate(k,stat=ok)
  if (ok>0) print "(a)", "Couldn't allocate node"
end if

END SUBROUTINE get

!-----------------------------------------------------------------------

SUBROUTINE link(k)
!  Insert a new link after K and point K at it

TYPE(node), POINTER :: k
TYPE(node), POINTER :: kz=>null()

call get(kz)
kz%next=>k%next
k%next=>kz
k=>kz
kz=>null()

END SUBROUTINE link

!-----------------------------------------------------------------------

SUBROUTINE clear(k)
!  Free chain K and nullify K

TYPE(node), POINTER :: k
TYPE(node), POINTER :: i=>null()

if (associated(k)) then
  if (associated(free)) then
    i=>free
    free=>k%next
    k%next=>i
    i=>null()
  else
    free=>k%next
    k%next=>null()
  end if
  k=>null()
endif

END SUBROUTINE clear

!-----------------------------------------------------------------------

SUBROUTINE copy(k,kz)
!  Copy chain K to KZ

TYPE(node), POINTER :: k, kz

if (.not. associated(k)) call report("Source for COPY is undefined")
if (associated(k,kz)) call report("Source and target of COPY are identical")

call clear(kz)
if (iszero(k)) then
  call zero(kz)
else
  call unit(kz)
  kz%x=k%x
  do
    k=>k%next
    if (k%prime .lt. 0) exit
    call link(kz)
    kz%prime=k%prime
    kz%x=k%x
  end do
  !  Step to -1 entry of kz
  kz=>kz%next 
end if

END SUBROUTINE copy

!-----------------------------------------------------------------------

SUBROUTINE rename(k1,k2)

TYPE(node), POINTER :: k1, k2

call clear(k2)
k2=>k1
k1=>null()

END SUBROUTINE rename

!-----------------------------------------------------------------------

SUBROUTINE unit(k)
!  Set K to unity

TYPE(node), POINTER :: k

call clear(k)
call get(k)
k%next=>k
k%prime=-1
k%x=0

END SUBROUTINE unit

!-----------------------------------------------------------------------

SUBROUTINE zero(k)
!  Set K to zero

TYPE(node), POINTER :: k

call clear(k)
call get(k)
k%next=>k
k%prime=0
k%x=0

END SUBROUTINE zero

!-----------------------------------------------------------------------

LOGICAL FUNCTION iszero(k)

TYPE(node), POINTER :: k

if (.not. associated(k)) then
  call report("Undefined rrf")
else
  iszero=(k%prime==0)
end if

END FUNCTION iszero

!-----------------------------------------------------------------------

LOGICAL FUNCTION nonzero(k)

TYPE(node), POINTER :: k

if (.not. associated(k)) then
  call report("Undefined rrf")
else
  nonzero=(k%prime.ne.0)
end if

END FUNCTION nonzero

!-----------------------------------------------------------------------

LOGICAL FUNCTION equal(k1,k2)

TYPE(node), POINTER :: k1, k2
TYPE(node), POINTER :: kz1, kz2

if (k1%prime==0 .and. k2%prime==0) then
  equal=.true.
else if (k1%prime==0 .or. k2%prime==0) then
  !  One is null, one not
  equal=.false.
else
  kz1=>k1
  kz2=>k2
  !  Standardise signs
  k1%x=modulo(k1%x,4)
  k2%x=modulo(k2%x,4)
  do
    kz1=>kz1%next
    kz2=>kz2%next
    if (kz1%prime .ne. kz2%prime .or. kz1%x .ne. kz2%x) then
      equal=.false.
      exit
    else if (kz1%prime .le. 0) then
      equal=.true.
      exit
    end if
  end do
end if

END FUNCTION equal

!-----------------------------------------------------------------------

SUBROUTINE exch(k1,k2)

TYPE(node), POINTER :: k1, k2
TYPE(node), POINTER :: i=>null()
i=>k1
k1=>k2
k2=>i
i=>null()
END SUBROUTINE exch

!-----------------------------------------------------------------------

SUBROUTINE chsign(k)

TYPE(node), POINTER :: k

if (k%prime .eq. 0) return
k%x=modulo(k%x+2,4)

END SUBROUTINE chsign

!-----------------------------------------------------------------------

SUBROUTINE conjugate(k)

TYPE(node), POINTER :: k

if (k%prime .eq. 0) return
k%x=modulo(4-k%x,4)

END SUBROUTINE conjugate

!-----------------------------------------------------------------------

SUBROUTINE list(k,string)
!  Print chain K

USE rrfdata

TYPE(node), POINTER :: k
CHARACTER(LEN=*), OPTIONAL :: string

TYPE(node), POINTER :: p

if (present(string)) then
  print "(a)", string
end if
if (.not. associated(k)) then
  print "(a)", '  undefined'
! else if (k%prime==0) then
!   print "(a)", '  zero'
else
  p=>k
  do
    print "(i10,i5)", p%prime, p%x
    p=>p%next
    if (associated(p)) then
      if (associated(p,k)) then
        !  Circular list -- back to start
        print "(1x)"
        exit
      endif
    else
      print "(a)", "         => null"
      exit
    end if
  end do
endif

END SUBROUTINE list

!-----------------------------------------------------------------------

SUBROUTINE report(m)

!  Print message m; then stop if hard_errors or error is true.
!  Otherwise set error true.

USE rrfdata

CHARACTER(LEN=*), INTENT(IN) :: m

write (lout,"(/a)") m

if (hard_errors .or. error) then
  write (lout,'(/1x,a)')                                                  &
      '**********       Program terminated       **********'
  stop
else
  error=.true.
endif

END SUBROUTINE report

!-----------------------------------------------------------------------

FUNCTION stri(i)
!  Return a left-justified character string representation of the
!  integer I

INTEGER, INTENT(IN) :: i
CHARACTER(LEN=20) :: stri, buffer

write (unit=buffer,fmt="(i12)") i
stri=adjustl(buffer)

END FUNCTION stri

END MODULE rrflist

