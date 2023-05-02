MODULE wigner
USE rrflist
USE rrfdata
USE factorials
USE arithmetic
#if defined(TEST9J) || defined(TEST6J) || defined(TEST3J)
USE conversions
#endif
IMPLICIT NONE

PRIVATE
PUBLIC threej, cg, three0, sixj, ninej

CONTAINS

!-----------------------------------------------------------------------

SUBROUTINE threej (j1,j2,j3,m1,m2,m3, x, k)

!  Wigner 3-j symbol:  k = (j1/x  j2/x  j3/x)
!                          (m1/x  m2/x  m3/x)

INTEGER, INTENT(IN) :: j1,j2,j3,m1,m2,m3, x
TYPE(node), POINTER :: k

call clear(k)
if (x .eq. 1) then
  call regge3(j1+m1, j2+m2, j3+m3,                                &
      j1-m1, j2-m2, j3-m3, k)
else
  if (mod(j1+j2+j3,x) .ne. 0                                      &
      .or. mod(j1+m1,x) .ne. 0                                    &
      .or. mod(j2+m2,x) .ne. 0                                    &
      .or. mod(j3+m3,x) .ne. 0) then
!  Sum of j's is not an integer, or a j value and its m value are not
!  both integer or both half-odd.
    call report('Invalid arguments for 3j symbol')
    return
  endif
  call regge3((j1+m1)/x, (j2+m2)/x, (j3+m3)/x,                    &
      (j1-m1)/x, (j2-m2)/x, (j3-m3)/x, k)
end if

END SUBROUTINE threej

!-----------------------------------------------------------------------

SUBROUTINE CG(j1,j2,j3, m1,m2,m3, x, k)

!  Clebsch-Gordan coefficient: k = <j1 j2 m1 m2 | j3 m3>

INTEGER, INTENT(IN) :: j1,j2,j3,m1,m2,m3, x
TYPE(node), POINTER :: k
TYPE(node), POINTER :: kz=>null()

call clear(k)
call threej(j1,j2,j3,m1,m2,-m3, x, k)
if (x .eq. 1) then
  call fromi(int(2*j3+1,long),kz)
  call root(kz)
#ifdef TESTCG
write (unit=lout,fmt="(a,t20)",advance="no") "sqrt(2j_3+1) = "
call display(kz, 3)
#endif
  call mult(k, kz, 1)
  if (modulo(j1-j2+m3,2) == 1) call chsign(k)
else
  call fromi(int(j3+1,long),kz)
  call root(kz)
  call mult(k, kz, 1)
  if (modulo((j1-j2+m3)/2,2) == 1) call chsign(k)
endif
call clear(kz)

END SUBROUTINE CG

!-----------------------------------------------------------------------

SUBROUTINE three0(j1,j2,j3, k)

!  Wigner 3-j symbol:  K = (J1  J2  J3)
!                          ( 0   0   0)

INTEGER, INTENT(IN) :: j1,j2,j3
TYPE(node), POINTER :: k
! LOGICAL :: odd
INTEGER :: n, ng, l1,l2,l3

#ifdef TEST3J
print "(/a,3(i0,1x),a)", "Evaluating 3j symbol ", j1, j2, j3, "  0 0 0
#endif
call zero(k)
n=j1+j2+j3
if (mod(n,2) .ne. 0) return
ng=n/2
l1=-j1+j2+j3
l2=j1-j2+j3
l3=j1+j2-j3
if (l1 .lt. 0 .or. l2 .lt. 0 .or. l3 .lt. 0) return
if (n .ge. maxfct) then
  call report('Reference to untabulated factorial')
  return
endif
call setf(l1,k)
call multf(k,l2,1)
call multf(k,l3,1)
call multf(k,n+1,-1)
call root(k)
call multf(k,ng,1)
call multf(k,ng-j1,-1)
call multf(k,ng-j2,-1)
call multf(k,ng-j3,-1)
if (mod(ng,2) .ne. 0) call chsign(k)
#ifdef TEST
write (unit=lout,fmt="(a,t20)",advance="no") "3j = "
call display(k, 3)
#endif

END SUBROUTINE three0

!-----------------------------------------------------------------------

SUBROUTINE regge3 (jp1,jp2,jp3,jm1,jm2,jm3, k)

!  Wigner 3-j symbol        k = (j1 j2 j3)
!                               (m1 m2 m3)

!  and JP1=J1+M1, JM1=J1-M1, etc., so that all arguments are
!  integers

INTEGER, INTENT(IN) :: jp1,jp2,jp3,jm1,jm2,jm3
TYPE(node), POINTER :: k
TYPE(node), POINTER :: knu=>null(), kz=>null()
INTEGER :: k1,k2,k3, l, l1,l2,l3, nn, nu, numin, numax

call zero(k)
l=jp1+jp2+jp3
!  Result is zero if M1+M2+M3~=0
if (l .ne. jm1+jm2+jm3) return
l1=jp1-jm2
l2=jp2-jm3
l3=jp3-jm1
!  K1=-J1+J2+J3, etc.
k1=l-jp1-jm1
k2=l-jp2-jm2
k3=l-jp3-jm3
numax=min(k3,jm1,jp2)
numin=max(0,l2,-l3)
!  If NUMIN>NUMAX the triangle conditions are not satisfied, or
!  an M value is out of range
if (numin .gt. numax) return
if (l .ge. maxfct) then
  call report('Reference to untabulated factorial')
  return
end if
nn=(-1)**(numin+l1)
do nu=numin,numax
  call unit(knu)
  call multf(knu,nu,-1)
  call multf(knu,k3-nu,-1)
  call multf(knu,jm1-nu,-1)
  call multf(knu,jp2-nu,-1)
  call multf(knu,nu-l2,-1)
  call multf(knu,nu+l3,-1)
  if (nn .lt. 0) call chsign(knu)
  call add(k,knu)
  if (error) then !  Failure in addition
    call clear(knu)
    call clear(kz)
    return
  endif
  nn=-nn
end do

call unit(kz)
call multf(kz,k1, 1)
call multf(kz,k2, 1)
call multf(kz,k3, 1)
call multf(kz,l+1,-1)
call multf(kz,jm1,1)
call multf(kz,jp1,1)
call multf(kz,jp2,1)
call multf(kz,jm2,1)
call multf(kz,jp3,1)
call multf(kz,jm3,1)
call root(kz)
call mult(k,kz,1)
call clear(knu)
call clear(kz)

END SUBROUTINE regge3

!-----------------------------------------------------------------------

SUBROUTINE delta(a,b,c, x, k)

!  Delta coefficient sqrt[(a+b-c)!(a+c-b)!(b+c-a)!/(a+b+c+1)!]

INTEGER, INTENT(IN) :: a, b, c, x
TYPE(node), POINTER :: k

call clear(k)
call setf((a+b-c)/x, k)
call multf(k, (a+c-b)/x, 1)
call multf(k, (b+c-a)/x, 1)
call multf(k, (a+b+c)/x+1, -1)
call root(k)

END SUBROUTINE delta

!-----------------------------------------------------------------------

SUBROUTINE sixj(j1, j2, j3, l1, l2, l3, x, k)

!  Wigner 6-j symbol:  k = <j1/x  j2/x  j3/x>
!                          <l1/x  l2/x  l3/x>

!                      where x = 1 or 2

INTEGER, INTENT(IN) :: j1, j2, j3, l1, l2, l3, x
TYPE(node), POINTER :: k

TYPE(node), POINTER :: kz=>null(), knu=>null(), kd=>null()
INTEGER :: n, ns, nu, numin, numax

#ifdef TEST6J
print "(/a,2(3(i0,1x),1x),a,i1,a)",                                   &
    "Evaluating 6j symbol ", j1, j2, j3, l1, l2, l3, " (X=", x, ")"
#endif
call zero(k)
n=j1+j2+j3
if (x .eq. 2) then
  if (odd(n) .or. odd(j1+l2+l3) .or. odd(l1+j2+l3) .or. odd(l1+l2+j3)) return
  n=n/2
endif
numax=min((j1+j2+l1+l2)/x,(j2+j3+l2+l3)/x,(j3+j1+l3+l1)/x)
numin=max(n,(j1+l2+l3)/x,(l1+j2+l3)/x,(l1+l2+j3)/x)
#ifdef TEST6J
  print "(a,i2,a,i2)", "numin = ", numin, "   numax = ", numax
#endif
!  If NUMIN>NUMAX a triangle condition is not satisfied
if (numin .gt. numax) then
#ifdef TEST6J
  print "(a)", "6j is zero"
#endif
  return
endif
if (numax .ge. maxfct) then
  call report('Reference to untabulated factorial')
  return
endif

ns=1
if (odd(numin)) ns=-1
do nu=numin,numax
  !       write (lout,'(a,i3,a)') 'Factorial', nu+1, ' needed'
  call setf(nu+1,knu)
  call multf(knu,(j1+j2+l1+l2)/x-nu,-1)
  call multf(knu,(j2+j3+l2+l3)/x-nu,-1)
  call multf(knu,(j3+j1+l3+l1)/x-nu,-1)
  call multf(knu,nu-n,-1)
  call multf(knu,nu-(j1+l2+l3)/x,-1)
  call multf(knu,nu-(l1+j2+l3)/x,-1)
  call multf(knu,nu-(l1+l2+j3)/x,-1)
  if (ns .lt. 0) call chsign(knu)
#ifdef TEST6J
  write (unit=lout,fmt="(a,i0,a,t20)",advance="no") "nu = ", nu, "  term = "
  call display(knu, 1)
#endif
  call add(k,knu)
#ifdef TEST6J
  write (unit=lout,fmt="(a,t20)",advance="no") "Subtotal = "
  call display(k, 1)
#endif
  if (error) then   !  Failure in addition
    call clear(knu)
    return
  end if
  ns=-ns
end do
#ifdef TEST6J
write (unit=lout,fmt="(a,t20)",advance="no") "Sum ="
call display(k, 1)
#endif

call delta(j1,j2,j3, x, kz)
call delta(j1,l2,l3, x, kd)
call mult(kz,kd,1)
call delta(l1,j2,l3, x, kd)
call mult(kz,kd,1)
call delta(l1,l2,j3, x, kd)
call mult(kz,kd,1)
call mult(k,kz,1)
#ifdef TEST6J
write (unit=lout,fmt="(a,t20)",advance="no") "Coefficient"
call display(kz, 1)
write (unit=lout,fmt="(a,t20)",advance="no") "6J ="
call display(k, 1)
#endif

call clear(knu)
call clear(kd)
call clear(kz)

CONTAINS

LOGICAL FUNCTION odd(m)

INTEGER, INTENT(IN) :: m

odd=mod(m,2) .ne. 0

END FUNCTION odd

END SUBROUTINE sixj

!-----------------------------------------------------------------------

SUBROUTINE ninej (a,b,c, d,e,f, g,h,i, x, k)

!  Wigner 9-j symbol:  k = (a/x  b/x  c/x)
!                          (d/x  e/x  f/x)
!                          (g/x  h/x  i/x)

!                      where x = 1 or 2

INTEGER, INTENT(IN)  :: a,b,c,d,e,f,g,h,i, x
TYPE(node), POINTER :: k
TYPE(node), POINTER :: k1=>null(), k2=>null(), k3=>null()
INTEGER :: kk, kmax, kmin
INTEGER(KIND=long) :: z

kmax=min(a+i, d+h, b+f)
kmin=max(abs(a-i), abs(d-h), abs(b-f))
call zero(k)
!  If kmin>kmax a triangle condition is not satisfied
if (kmin .gt. kmax) return

if (max(min(a+i+h+g, a+kmax+h+g, i+kmax+d+g),                          &
    min(b+f+d+h, b+kmax+d+e, f+kmax+h+e),                              &
    min(a+i+f+b, a+kmax+f+c, i+kmax+b+c))                              &
    .ge. x*maxfct) then
  call report('Reference to untabulated factorial')
  return
endif

#ifdef TEST9J
  print "(/a,3(3(i0,1x),1x),a,i1,a)",                                   &
      "Evaluating 9j symbol ", a, b, c, d, e, f, g, h, i, " (X=", x, ")"
#endif
z=kmin*2/x+1
do kk=kmin,kmax,x
  call sixj(a,i,kk, h,d,g, x, k1)
#ifdef TEST9J
    write (unit=6,fmt="(a,6i2,a,t20)",advance="no") "6j: ", a,i,kk, h,d,g, " = "
    call display(k1,3)
#endif
  call sixj(b,f,kk, d,h,e, x, k2)
#ifdef TEST9J
    write (unit=6,fmt="(a,6i2,a,t20)",advance="no") "6j: ", b,f,kk, d,h,e, " = "
    call display(k2,3)
#endif
  call sixj(a,i,kk, f,b,c, x, k3)
#ifdef TEST9J
    write (unit=6,fmt="(a,6i2,a,t20)",advance="no") "6j: ", a,i,kk, f,b,c, " = "
    call display(k3,3)
#endif
  call mult(k1,k2,1)
  call mult(k1,k3,1)
  call multi(k1,z,1)
#ifdef TEST9J
  write (unit=6,fmt="(a,t20,i0/a,t35)",advance="no") "z = ", z, "Product = "
  call display(k1,3)
#endif
  call add(k,k1)
  z=z+2
end do
#ifdef TEST9J
    write (unit=6,fmt="(a,t35)",advance="no") "Sum = "
    call display(k,3)
#endif
if (mod(z,2_long) .eq. 0) call chsign(k)
#ifdef TEST9J
    write (unit=6,fmt="(a,t35)",advance="no") "9j  = "
    call display(k,3)
#endif

call clear(k1)
call clear(k2)
call clear(k3)

END SUBROUTINE ninej

END MODULE wigner
