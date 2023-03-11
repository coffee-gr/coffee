PROGRAM rrfcalc

!  This program carries out rrf calculations in the manner of a
!  Hewlett-Packard calculator, using Reverse Polish Notation for
!  the commands, which are read from the data stream a line at a
!  time.

USE rrfdata
USE rrflist
USE factorials
USE arithmetic
USE conversions
USE wigner
IMPLICIT NONE

INTERFACE
  SUBROUTINE start(quiet,primes)
  LOGICAL, INTENT(IN), OPTIONAL :: quiet
  INTEGER, INTENT(IN), OPTIONAL :: primes
  END SUBROUTINE start
END INTERFACE

INTEGER, PARAMETER :: TOP=25
LOGICAL :: verify, echo
CHARACTER(LEN=1) :: ibuf(80), c
CHARACTER(LEN=6) :: word
CHARACTER(LEN=20) :: prompt=':'
CHARACTER(LEN=80) :: buffer

INTEGER :: k(9),j(9), npr=1, vmode=3, ok
INTEGER :: nw=80
INTEGER :: i, l, m, n, x
INTEGER(KIND=long) :: nn
EQUIVALENCE (buffer, ibuf(1))

TYPE(node), POINTER :: ka=>null(), kb=>null()
!  Stack and memory
TYPE(rrf) :: stack(TOP), s(TOP)
!  Stack pointer
INTEGER :: ipt

call start(quiet=.true.)
!  Soft errors
hard_errors=.false.
!  debug=.true.

ipt=0
vmode=3
!  Clear stack and memory
do i=1,TOP
  s(i)%rrf=>null()
  stack(i)%rrf=>null()
end do
echo=.false.
verify=.true.


lines: do
!  Prompt for input if a prompt string is set
if (npr .gt. 0) write (lout,'(a)', advance="no") trim(prompt)//" "
!  Read new input line
  read (lin,"(a)", iostat=ok) buffer
  if ( ok<0 ) exit lines  ! End-of-file

  error=.false.
!  Copy input to output if echo switch is on
  if (echo) write (lout,"(/a)") trim(buffer)

  l=1
  commands: do

    if (buffer(l:) == "") exit
    !  Check for 1-character commands and space
    c=buffer(l:l)
    if (debug) print "(2a)", "Char: ", c
    select case(c)
    case(space)
      l=l+1
      cycle
    case(bra)
      !  Ignore rest of line
      exit
    case(plus,minus,star,slash)
      l=l+1
      if (ipt .le. 1) then
!  Too few items in stack
        call abort('Too few items in stack')
        exit
      else
        select case(c)
        case(plus,minus)
!  Add or subtract
          call rename(stack(ipt)%rrf,ka)
          ipt=ipt-1
          if (c .eq. minus) call chsign(ka)
          call copy(stack(ipt)%rrf,kb)
          call add(kb,ka)
          call clear(ka)
          if (error) then
            call abort("Failure in addition")
            exit commands
          endif
          call rename(kb,stack(ipt)%rrf)
          verify=.true.
        case(star,slash)
          if (c .eq. star) then
            !  Multiply
            i = 1
          else
            !  Divide
            i = -1
            if (iszero(stack(ipt)%rrf)) then
              call abort('Attempt to divide by zero')
              exit
            end if
          endif
          call mult(stack(ipt-1)%rrf,stack(ipt)%rrf,i)
          call clear(stack(ipt)%rrf)
          ipt = ipt-1
          verify=.true.
        end select
        cycle
      endif
      cycle
    end select

    !  No recognised single character. Look for command word terminated
    !  by space.
    m=index(buffer(l:)," ")
    if (m == 0) then
      word=buffer(l:)
      l=len(buffer)
    else
      word=buffer(l:l+m-2)
      l=l+m-1
    endif
    if (debug) print "(2a)", "Word: ", word
    call upcase(word)
    select case(word)
    case("")
      exit
    case("POP")
      !  Pop off stack
      if (ipt .gt. 0) then
        call clear(stack(ipt)%rrf)
        ipt = ipt-1
      end if
      verify=.true.
    case("SQRT")  
      !  Take square root
      call copy(stack(ipt)%rrf,kb)
      call root(kb)
      if (error) then
        call abort("Current input line abandoned")
        call clear(kb)
      else
        call rename(kb,stack(ipt)%rrf)
      end if
      verify=.true.
    case("POWER","PWR")
      !  Raise to power
      x=2
      call number(buffer,l, nn,x)
      i=nn
      if (x.le.0 .or. x .gt. 2) then
        call abort("Syntax error")
        exit
      end if
      if (iszero(stack(ipt)%rrf).and.i.lt.0) then
        call abort("Division by zero")
        exit
      end if
      call copy(stack(ipt)%rrf,kb)
      call power(kb,i)
      if (x.eq.2) call root(kb)
      call rename(kb,stack(ipt)%rrf)
      verify=.true.
    case("CHS")
      !  Change sign
      call chsign(stack(ipt)%rrf)
      verify=.true.
    case("PUSH","DUP")
      !  Push onto stack (duplicate)
      if (ipt .gt. 0) call copy(stack(ipt)%rrf,stack(ipt+1)%rrf)
      ipt = ipt+1
      verify=.true.
    case("SWAP","SWOP","EXCH")
      !  Swop
      ka=>stack(ipt)%rrf
      stack(ipt)%rrf=>stack(ipt-1)%rrf
      stack(ipt-1)%rrf=>ka
      verify=.true.
    case("CLR","CLEAR")
      !  Clear stack
      do i=1,TOP
        call clear(stack(i)%rrf)
      end do
      ipt = 0
      verify=.true.
    case default
      !  Read an integer from the buffer
      l=l-m+1
      x=1
      call number(buffer,l, nn,x)
      if (x.le.0) then
        call abort("Malformed number or unrecognised command")
        exit
      end if
      if (ipt==TOP) then
        call abort("Stack full")
      else
        ipt=ipt+1
        if (buffer(l:l) == shriek) then
          l=l+1
          call setf(int(nn),stack(ipt)%rrf)
        else
          call fromi(nn,stack(ipt)%rrf)
        end if
      end if
      verify=.true.
    case("PP")
!  Read an rrf number in power-of-prime form
      call fromch(buffer, l, 0, stack(ipt+1)%rrf)
      verify=verify .or. .not. error
      if (error) then
        call abort("Malformed power-of-primes")
        exit
      end if
      ipt = ipt+1
      verify=.true.
    case("VP")
      !  Verify top of stack
      call show(ipt,1)
      verify=.false.
    case("VF")
      !  Verify top of stack
      call show(ipt,2)
      verify=.false.
    case("VI")
      !  Verify top of stack
      call show(ipt,3)
      verify=.false.
    case("0J","3J","6J","9J","CG")
      !  Wigner 3n-j symbols and Clebsch-Gordan coefficient
      select case(word)
      case("0J")
        n=3
      case("CG","3J","6J")
        n=6
      case("9J")
        n=9
      end select
      !  Read arguments
      x = 1
      do i=1,n
        k(i)=2
        call number(buffer,l, j(i),k(i))
        if (k(i).le.0) then
          call abort("Malformed number")
          exit commands
        endif
        if (k(i).eq.2) x = 2
      end do
      if (x .eq. 2) then
        do i=1,n
          if (k(i) .eq. 1) j(i) = 2*j(i)
        end do
      end if
      select case(word)
      case("0J")
        !  Wigner 3j symbol, all m values zero
        if (x .ne. 1) then
          call report('Half-integer argument in 0j')
        else
          call three0(j(1),j(2),j(3), stack(ipt+1)%rrf)
        endif
      case("CG")
        !  Clebsch-Gordan coefficient
          call CG(j(1),j(2),j(3),j(4),j(5),j(6),x,stack(ipt+1)%rrf)
      case("3J")
        !  Wigner 3j symbol
          call threej(j(1),j(2),j(3),j(4),j(5),j(6),x,stack(ipt+1)%rrf)
      case("6J")
        !  Wigner 6j symbol
        call sixj(j(1),j(2),j(3),j(4),j(5),j(6),x,stack(ipt+1)%rrf)
      case("9J")
        !  Wigner 9j symbol
        call ninej(j(1),j(2),j(3),j(4),j(5),j(6),j(7),j(8),j(9),       &
            x,stack(ipt+1)%rrf)
      end select
      if (error) then
        call abort("Error in n-j symbol")
        exit
      end if
      ipt = ipt+1
      verify=.true.
    case("VMP")
!  Change verification mode to primes
      vmode = 1
      verify=.true.
    case("VMF")
!  Change verification mode to floating-point
      vmode = 2
      verify=.true.
    case("VMI")
!  Change verification mode to integer
      vmode = 3
      verify=.true.
    case("LIST")
!  List items in stack
      do m=ipt,1,-1
        call show(m,vmode)
      end do
      verify=.false.
      ! call list(-90)
    case("V")
!  Display top-of-stack
      call show(ipt,vmode)
      verify=.false.
    case("STO","STORE")
!  Store
      x=1
      call number(buffer,l, i,x)
      if (x .ne. 1 .or. i .lt. 1 .or. i .gt. 10) then
        call abort("Invalid memory number")
        exit
      else
        call copy(stack(ipt)%rrf,s(i)%rrf)
      endif
    case("RCL","RECALL")
!  Recall
      x=1
      call number(buffer,l, i,x)
      if (x .ne. 1 .or. i .lt. 1 .or. i .gt. 10) then
        call abort("Invalid memory number")
        exit
      else
        ipt=ipt+1
        call copy(s(i)%rrf,stack(ipt)%rrf)
      endif
      verify=.true.
    case("Q","QUIT","STOP")
      STOP
    case("ECHO+","ECHO")
      echo=.true.
    case("ECHO-","QUIET")
      echo=.false.
    case("PROMPT")
!  Set prompt string.  Delimited by a pair of apostrophes, or of
!  any other character which does not appear in the string.  Can
!  be null.  Truncated if longer than 20 characters.
      npr=0
      do while (buffer(l:l) .eq. space)
        l=l+1
      end do
      if (l .lt. nw) then
        c=buffer(l:1)
        l=l+1
        do while (buffer(l:l) .ne. c .and. l .le. nw)
          npr=npr+1
          if (npr .le. 20) prompt(npr:npr)=buffer(l:1)
          l=l+1
        end do
        l=l+1
      endif
    case("CLF")
!  Clear factorial table
      call clearf(2,100)
    case("HELP","?")
!  Display command list
      print "(a)",                                                      &
          'integer      push value onto stack',                         &
          '-            subtract top of stack from next',               &
          '+            add top two items in stack',                    &
          '*            multiply top two items in stack',               &
          '/            divide top of stack into next',                 &
          'CHS          change sign',                                   &
          'SQRT         square root',                                   &
          'POWER,PWR n  raise to power',                                &
          'POP          pop stack (remove top item)',                   &
          'PUSH         push stack (copy top item)',                    &
          'SWOP,SWAP    swop top two items in stack',                   &
          'PP           read rrf in power-of-prime form'
      print "(a)",                                                      &
          '3J,6J,9J     Wigner 3n-j symbol',                            &
          '0J j1 j2 j3 = 3J j1 j2 j3 0 0 0',                            &
          'CG j1 j2 j3 m1 m2 m3',                                       &
          '             Clebsch-Gordan coefficient <j1j2m1m2|j3m3>',    &
          'STO n        store in memory n (1 <= n <= 10)',              &
          'RCL n        recall from memory n',                          &
          'CLEAR,CLR    clear stack',                                   &
          'CLF          clear factorial table',                         &
          'Q,QUIT       stop',                                          &
          '( or !       ignore rest of line',                           &
          'V            verify top of stack',                           &
          'VP           verify in power-of-prime form',                 &
          'VF           verify as floating-point',                      &
          'VI           verify as (I1/I2)*sqrt(I3/I4)'
      print "(2a)",                                                     &
          'VMF          switch default verification mode',              &
                                             ' to floating-point',      &
          'VMI             "      "          "       "  ',              &
                                             ' to (I1/I2)*sqrt(I3/I4)', &
          'VMP             "      "          "       "  ',              &
                                             ' to power-of-prime'
      print "(a)",                                                      &
          'LIST         list stack',                                    &
          'PROMPT "ss"  set prompt string to ss',                       &
          'ECHO         reflect command line to output',                &
          'QUIET        do not reflect commands (default)',             &
          'HELP,?       type list of commands'
    end select
  end do commands
!  Before reading a new line of input, the top-of-stack is verified if
!  there has been any change to the stack or the verify mode since the
!  last verification.
  if (verify) then
    call show(ipt,vmode)
  endif
  verify=.false.
end do lines

stop

CONTAINS

SUBROUTINE abort(string)
!  Abandon current line
CHARACTER(LEN=*) :: string
INTEGER :: j
error=.false.
write (lout,"(a)") string
write (lout,'(132a1)') ibuf
write (lout,'(132a1)') (space, j=3,l), star
END SUBROUTINE abort

!-----------------------------------------------------------------------

SUBROUTINE show(n,mode)

INTEGER :: n, mode

if (n .ge. 1) then
  call display(stack(n)%rrf,mode)
else
  write (lout,"(a)") "Stack empty"
endif

END SUBROUTINE show

!-----------------------------------------------------------------------

SUBROUTINE upcase(word)
IMPLICIT NONE
CHARACTER(LEN=*) :: word

CHARACTER(LEN=26) :: uc="ABCDEFGHIJKLMNOPQRSTUVWXYZ",            &
    lc="abcdefghijklmnopqrstuvwxyz"
INTEGER :: i, k

do i=1,len(word)
  k=index(lc,word(i:i))
  if (k .ne. 0) word(i:i)=uc(k:k)
end do

END SUBROUTINE upcase

END PROGRAM rrfcalc
