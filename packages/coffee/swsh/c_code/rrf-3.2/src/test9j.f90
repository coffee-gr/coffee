PROGRAM test9j

USE rrfdata
USE rrflist
USE factorials
USE arithmetic
USE conversions, ONLY : number, char4i
USE wigner
USE input
IMPLICIT NONE

INTERFACE
  SUBROUTINE start(quiet)
  LOGICAL, INTENT(IN), OPTIONAL :: quiet
  END SUBROUTINE start
END INTERFACE


INTEGER :: a, b, c, d, e, f, g, h, i, m, n, p, x, minj=0, maxj=12, rule
TYPE(node), POINTER :: k=>null(), k1=>null(), k2=>null(), sum=>null(), &
    z=>null(), one=>null()
LOGICAL :: eof, verbose=.false., printall=.false., allok, single=.false.
CHARACTER(LEN=16) :: key

call start(.true.)
call unit(one)

data: do
  call read_line(eof)
  if (eof) exit
  call readu(key)
  select case(key)
    case("VERBOSE")
      verbose=.true.
    case("QUIET")
      verbose=.false.
    case ("PRINT")
      call readu(key)
      select case(key)
      case("ALL")
        printall=.true.
      case("FAILED")
        printall=.false.
      end select
    case("MIN")
      call readi(minj)
      minj=2*minj
    case("MAX")
      call readi(maxj)
      maxj=2*maxj


    case("SUM-RULE","SUMRULE","SUM", "TEST")
      rule=1
      do while (item<nitems)
        call readu(key)
        select case(key)
        case("RULE")
          call readi(rule)
        case("SINGLE")
          single=.true.
        case default
          call reread(-1)
          call readi(rule)
        end select
      end do
      select case(rule)

      case(1)
        !  Sum rule 1 (orthogonality)
        x=2
        allok=.true.
        print "(/a)", "Sum rule 1"
        if (printall) print "(/a)", " a    b    d    e    g    h    i    m    n"
        if (single) then
          do
            call read_line(eof)
            if (eof) exit data
            if (nitems==0) cycle
            call readargs(a,b,d,e,g,h,i,m,n)
            call check(a,b, d,e, g,h, i, m,n, .true.)
          end do
        else
          do a=minj,maxj
            do b=minj,a
              do d=minj,a
                do e=minj,a
                  do g=max(minj,abs(a-d)),min(a+d,maxj),2
                    do h=max(minj,abs(b-e)),min(b+e,maxj),2
                      do i=max(minj,abs(g-h)),min(g+h,maxj),2
                        do m=max(minj,abs(a-d)),min(a+d,maxj),2
                          do n=max(minj,abs(b-e)),min(b+e,maxj),2
                            call check(a,b, d,e, g,h, i, m,n, printall)
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end if
        call clearall

    end select
  end select
end do data
if (allok) print "(a)", "Test completed -- no errors"

CONTAINS

SUBROUTINE check(a,b, d,e, g,h, i, m,n, print)

INTEGER, INTENT(IN) :: a,b,d,e,g,h,i,m,n
LOGICAL, INTENT(IN) :: print

INTEGER :: c, f
LOGICAL ok
CHARACTER(LEN=45) :: args
CHARACTER(LEN=8) :: verdict
CHARACTER(LEN=40) :: string

call zero(sum)
do c=abs(a-b),a+b,2
  do f=abs(d-e),d+e,2
    call ninej(a, b, c, d, e, f, g, h, i, x, k1)
    call ninej(a, b, c, d, e, f, m, n, i, x, k2)
    call fromi(int((c+1)*(f+1)*(g+1)*(h+1),long),k)
    call mult(k,k1)
    call mult(k,k2)
    call add(sum,k)
  end do
end do
if (g==m .and. h==n) then
  ok=equal(sum,one)
else
  ok=iszero(sum)
end if
if (ok) then
  verdict="  O.K."
  ok=.true.
else
  verdict=" Wrong!"
  if (allok .and. .not. printall)            &
      print "(/a)", " a    b    d    e    g    h    i    m    n"
  allok=.false.
end if
if (print .or. .not. ok) then
  p=1
  args=""
  call put(a, args, p)
  call put(b, args, p)
  call put(d, args, p)
  call put(e, args, p)
  call put(g, args, p)
  call put(h, args, p)
  call put(i, args, p)
  call put(m, args, p)
  call put(n, args, p)
  p=1
  string=""
  call char4i(sum, string, p, 40)
  print "(4a)", args, verdict, "  Sum = ", trim(string)
end if

END SUBROUTINE check

SUBROUTINE readargs(a,b,d,e,g,h,i,m,n)

INTEGER, INTENT(OUT) :: a,b,d,e,g,h,i,m,n

call readarg(a)
call readarg(b)
call readarg(d)
call readarg(e)
call readarg(g)
call readarg(h)
call readarg(i)
call readarg(m)
call readarg(n)

END SUBROUTINE readargs

SUBROUTINE readarg(a)

INTEGER, INTENT(OUT) :: a
INTEGER :: p, x
CHARACTER(LEN=10) :: buff

call reada(buff)
p=1
x=2
call number(buff,p,a,x)
if (x==1) a=2*a

END SUBROUTINE readarg

SUBROUTINE put(j, s,p)

CHARACTER(LEN=*) :: s
INTEGER, INTENT(IN) :: j
INTEGER, INTENT(INOUT) :: p

if (mod(j,2) == 0) then
  write (unit=s(p:p+1),fmt="(i2)") j/2
else
  write (unit=s(p:p+3),fmt="(i2,a2)") j,"/2"
endif
p=p+5

END SUBROUTINE put

SUBROUTINE clearall

call clear(k)
call clear(k1)
call clear(k2)
call clear(sum)
call clear(z)
call clear(one)

END SUBROUTINE clearall

END PROGRAM test9j
