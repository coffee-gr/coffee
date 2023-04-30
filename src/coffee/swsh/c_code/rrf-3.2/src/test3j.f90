PROGRAM test3j

!  Test the RRF program, and specifically the 3j symbol routine, using the
!  sum rules in Brink & Satchler, Appendix II.

USE rrfdata
USE rrflist
USE factorials
USE arithmetic
USE conversions
USE wigner
USE input
IMPLICIT NONE

INTERFACE
  SUBROUTINE start(quiet)
  LOGICAL, INTENT(IN), OPTIONAL :: quiet
  END SUBROUTINE start
END INTERFACE


INTEGER :: a, b, c, d, e, f, p, x, minj=0, maxj=12, test
INTEGER :: ma, mb, mc, md, me, mf
TYPE(node), POINTER :: k=>null(), k1=>null(), k2=>null(), k3=>null(),  &
    k4=>null(), sum=>null(), z=>null(), one=>null()
LOGICAL :: eof, verbose=.false., ok, allok, printall=.false.
CHARACTER(LEN=16) :: key
CHARACTER(LEN=30) :: args
CHARACTER(LEN=8) :: verdict
CHARACTER(LEN=40) :: string

call start(.true.)

do
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
      test=1
      do while (item<nitems)
        call readu(key)
        select case(key)
        case("RULE")
          call readi(test)
        case default
          call reread(-1)
          call readi(test)
        end select
      end do
      select case(test)

      case(1)
        !  Sum rule 1
        call unit(one)
        x=2
        print "(/a)", "Sum rule 1"
        if (printall) print "(/a)", " a    b    c    d    e    f"
        allok=.true.
        do a=minj,maxj
          do b=a,maxj
            do e=abs(a-b),min(a+b,maxj),2
              do d=minj,maxj
                do c=abs(d-e),min(d+e,maxj),2
                  do f=abs(a-c),min(a+c,maxj),2
                    call zero(sum)
                    mc=c
                    do ma=-a,a,2
                      do mb=-b,b,2
                        me=ma+mb
                        if (me<-e .or. me>e) cycle
                        call threej(a,b,e, ma,mb,-me,x, k1)
                        md=-mc-me
                        if (md<-d .or. md>d) cycle
                        call threej(d,c,e, md,mc,me, x, k2)
                        mf=mb+md
                        if (mf<-f .or. mf>f) cycle
                        call threej(b,d,f, mb,md,-mf, x, k3)
                        call threej(c,a,f, mc,ma,mf, x, k4)
                        call fromi(int(c+1,long),k)
                        call mult(k,k1)
                        call mult(k,k2)
                        call mult(k,k3)
                        call mult(k,k4)
                        if (modulo((a+b+c+d+f-e-ma-md)/2,2) .ne. 0) call chsign(k)
                        call add(sum,k)
                      end do
                    end do
                    call sixj(a, b, e, d, c, f, x, z)
                    p=1
                    args=""
                    call put(a, args, p)
                    call put(b, args, p)
                    call put(e, args, p)
                    call put(c, args, p)
                    call put(d, args, p)
                    call put(f, args, p)
                    if (equal(sum,z)) then
                      verdict="  O.K."
                      ok=.true.
                    else
                      verdict=" Wrong!"
                      ok=.false.
                      if (allok .and. .not. printall) print "(/a)", " a    b    c    d    e    f"
                      allok=.false.
                    end if
                    if (printall .or. .not. ok) then
                      p=1
                      string=""
                      call char4i(sum, string, p, 40)
                      print "(4a)", args, verdict, "  Sum = ", trim(string)
                    end if
                  end do
                end do
              end do
            end do
          end do
        end do
        if (allok) print "(a)", "Test completed -- no errors"

    end select
  end select
end do

CONTAINS

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

END PROGRAM test3j
