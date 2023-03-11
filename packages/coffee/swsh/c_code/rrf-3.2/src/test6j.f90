PROGRAM test6j

!  Test the RRF program, and specifically the 6j symbol routine, using the
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


INTEGER :: a, b, c, d, e, f, g, p, x, minj=0, maxj=12, test
TYPE(node), POINTER :: k=>null(), k1=>null(), k2=>null(), sum=>null(), &
    z=>null(), one=>null()
LOGICAL :: eof, verbose=.false., printall, allok, ok
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
    case("LIST")

      !  Tabulate
      !  a - f are double the actual j value
      do a=minj,maxj
        do b=a,maxj
          do d=0,maxj
            do e=abs(a-b),min(a+b,maxj),2
              do c=abs(d-e),min(d+e,maxj),2
                do f=abs(a-c),min(a+c,maxj),2
                  x=2
                  call sixj(a, b, c, d, e, f, x, k)
                  if ( nonzero(k) ) then
                    p=1
                    args=""
                    call put(a, args,p)
                    call put(b, args,p)
                    call put(c, args,p)
                    call put(d, args,p)
                    call put(e, args,p)
                    call put(f, args,p)
                    p=1
                    string=""
                    call char4i(k, string, p, 40)
                    print "(a,2x,a)", args, trim(string)
                  end if
                end do
              end do
            end do
          end do
        end do
      end do

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
        if (printall) print "(/a)", " a    b    f"
        allok=.true.
        do a=minj,maxj
          do b=a,maxj
            d=a
            e=b
            do f=abs(a-b),min(a+b,maxj),2
              call zero(sum)
              do c=abs(a-b),abs(a+b),2
                call sixj(a, b, c, d, e, f, x, k)
                if (verbose) then
                  p=1
                  args=""
                  call put(a, args,p)
                  call put(b, args,p)
                  call put(c, args,p)
                  call put(d, args,p)
                  call put(e, args,p)
                  call put(f, args,p)
                  p=1
                  string=""
                  call char4i(k, string, p, 40)
                  print "(a,2x,a)", args, trim(string)
                end if
                call fromi(int(c+1,long),z)
                call mult(k,z)
                call add(sum,k)
              end do
              if (mod(a+b,2) .ne. 0) call chsign(sum)
              if (equal(sum,one)) then
                verdict="  O.K."
                ok=.true.
              else
                verdict=" Wrong!"
                ok=.false.
                if (allok .and. .not. printall) print "(/a)", " a    b    f"
                allok=.false.
              end if
              if (printall .or. .not. ok) then
                p=1
                args=""
                call put(a, args, p)
                call put(b, args, p)
                call put(f, args, p)
                p=1
                string=""
                call char4i(sum, string, p, 40)
                print "(4a)", args(1:15), verdict, "  Sum = ", trim(string)
              end if
            end do
          end do
        end do
        if (allok) print "(a)", "Test completed -- no errors"

      case(2)
        !  Sum rule 2
        call clear(sum)
        call clear(z)
        x=2
        print "(/a)", "Sum rule 2"
        if (printall) print "(/a)", " a    b    f"
        do a=minj,maxj
          do b=a,maxj
            d=b
            e=a
            do f=abs(a-b),min(a+b,maxj),2
              call zero(sum)
              do c=abs(a-b),abs(a+b),2
                call sixj(a, b, c, d, e, f, x, k)
                if (verbose) then
                  p=1
                  args=""
                  call put(a, args,p)
                  call put(b, args,p)
                  call put(c, args,p)
                  call put(d, args,p)
                  call put(e, args,p)
                  call put(f, args,p)
                  p=1
                  string=""
                  call char4i(k, string, p, 40)
                  print "(a,2x,a)", args, trim(string)
                end if
                call fromi(int(c+1,long),z)
                call mult(k,z)
                if (mod((a+b+c)/2,2) .ne. 0) call chsign(k)
                call add(sum,k)
              end do
              if (f==0) then
                call fromi(int((a+1)*(b+1),long),z)
                call root(z)
              else
                call zero(z)
              end if
              if (equal(sum,z)) then
                verdict="  O.K."
                ok=.true.
              else
                verdict=" Wrong!"
                ok=.false.
                if (allok .and. .not. printall) print "(/a)", " a    b    f"
                allok=.false.
              end if
              if (printall .or. .not. ok) then
                args=""
                p=1
                call put(a, args,p)
                call put(b, args,p)
                call put(f, args,p)
                p=1
                string=""
                call char4i(sum, string, p, 40)
                print "(4a)", args(1:15), verdict, "Sum = ", trim(string)
              end if
            end do
          end do
        end do
        if (allok) print "(a)", "Test completed -- no errors"

      case(3)
        print "(/a)", "Sum rule 3"
        if (printall) print "(/a)", " a    b    d    e    f    g"
        call clear(sum)
        call clear(z)
        call unit(one)
        x=2
        allok=.true.
        !  a - g are double the actual j value
        do a=minj,maxj
          do b=a,maxj
            do d=0,maxj
              do e=0,maxj
                if (mod(d+e,2) .ne. mod(a+b,2)) cycle
                do f=max(abs(a-e),abs(b-d)),min(a+e,b+d,maxj),2
                  do g=f,min(a+e,b+d,maxj),2
                    call zero(sum)
                    do c=a+b,abs(a-b),-2
                      x=2
                      call sixj(a, b, c, d, e, f, x, k1)
                      call sixj(a, b, c, d, e, g, x, k2)
                      call fromi(int((c+1)*(f+1),long),k)
                      call mult(k,k1)
                      call mult(k,k2)
                      call add(sum,k)
                    end do
                    if (f .ne. g .and. iszero(sum)) then
                      verdict="  O.K."
                    else if (f==g) then
                      if (equal(sum,one)) verdict="  O.K."
                      ok=.true.
                    else
                      verdict=" Wrong!"
                      ok=.false.
                      if (allok .and. .not. printall)                  &
                          print "(/a)", " a    b    d    e    f    g"
                      allok=.false.
                    end if
                    if (printall .or. .not. ok) then
                      args=""
                      p=1
                      call put(a, args,p)
                      call put(b, args,p)
                      call put(d, args,p)
                      call put(e, args,p)
                      call put(f, args,p)
                      call put(g, args,p)
                      p=1
                      string=""
                      call char4i(sum, string, p, 40)
                      print "(4a)", args, verdict, "Sum = ", trim(string)
                    end if
                  end do
                end do
              end do
            end do
          end do
        end do
        call clearall
        if (allok) print "(a)", "Test completed -- no errors"

      case(4)
        print "(/a)", "Sum rule 4"
        if (printall) print "(/a)", " a    b    d    e    f    g"
        call clear(sum)
        call clear(z)
        call unit(one)
        x=2
        allok=.true.
        !  a - g are double the actual j  value
        do a=minj,maxj
          do d=a,maxj
            do b=0,maxj
              do e=0,maxj
                if (mod(b+e,2) .ne. mod(a+d,2)) cycle
                do c=max(abs(a-b),abs(d-e)),min(a+b,d+e,maxj),2
                  do f=max(abs(a-b),abs(d-e)),min(a+b,d+e,maxj),2
                    call zero(sum)
                    do g=abs(a-d),a+d,2
                      x=2
                      call sixj(a, d, g, e, b, c, x, k1)
                      call sixj(a, d, g, b, e, f, x, k2)
                      call fromi(int(g+1,long),k)
                      call mult(k,k1)
                      call mult(k,k2)
                      if (mod((c+f+g)/2,2) .ne. 0) call chsign(k)
                      call add(sum,k)
                    end do
                    call sixj(a, b, c, d, e, f, x, z)
                    if (equal(sum,z)) then
                      verdict="  O.K."
                      ok=.true.
                    else
                      verdict=" Wrong!"
                      ok=.false.
                      if (allok .and. .not. printall)                  &
                          print "(/a)", " a    b    d    e    f    g"
                      allok=.false.
                    end if
                    if (printall .or. .not. ok) then
                      args=""
                      p=1
                      call put(a, args,p)
                      call put(b, args,p)
                      call put(c, args,p)
                      call put(d, args,p)
                      call put(e, args,p)
                      call put(f, args,p)
                      p=1
                      string=""
                      call char4i(sum, string, p, 40)
                      print "(4a)", args, verdict, "Sum = ", trim(string)
                    end if
                  end do
                end do
              end do
            end do
          end do
        end do
        call clearall
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

END PROGRAM test6j
