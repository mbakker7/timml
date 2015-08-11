! FORTRAN implementation of higher order Bessel line-elements
! Accurate for orders 0 - 7. Give numerical problems at edge for orders higher than 7
! This is not checked; any order can be specified and may result in bogus values
! Details of formulation in maqelements.tex
! (c) Mark Bakker, 2003

subroutine potbeslsho(x,y,x1in,y1in,x2in,y2in,Naquifers,lambdain,order,rv)

    ! Input:
    !   x,y: Point where potential is computed
    !   x1in,y1in: Begin point of line-sink
    !   x2in,y2in: End point of line-sink
    !   Naquifers: Number of aquifers
    !   lambdain(Naquifers): Array with zero in first spot and lambda's in remaining spots
    !   order: Order of the line-sink
    !   rv(Naquifers): Array to store return value
    ! Output:
    !   rv(Naquifers): Value of potentials with Laplace value in first spot
    !                  and mod.Helmholtz potentials in remaining spots

    implicit none
    ! Input / Output
    real*8, intent(in) :: x,y,x1in,y1in,x2in,y2in
    integer, intent(in) :: Naquifers, order
    real*8, intent(in) :: lambdain(*)
    real*8,intent(inout) :: rv(*)

    ! Locals
    integer, parameter :: Nterms = 8

    real*8, parameter, dimension(0:8) :: rRange = (/0.d0,1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,7.d0,8.d0/)
    real*8, parameter :: rZERO=0.d0, rONE=1.d0, rTWO=2.d0, rFOUR=4.d0
    real*8, parameter :: pi=3.1415926535897931d0, tiny = 1.d-8

    integer :: power, n, m, i
    real*8 :: Rconv, Lin, pot, biglab
    complex*16 :: zin, z1in, z2in, z, zplus1, zmin1
    complex*16 :: pcor, comega

    real*8, dimension(0:Nterms) :: ac,bc
    real*8, dimension(Naquifers) :: lambda
    ! integer, parameter :: Nmax = 50
    ! real*8, dimension(0:Nmax) :: rRange  ! Maximum order set to 50 (like you ever want that!)
    real*8, dimension(0:Nterms,0:Nterms) :: rbinom

    ! Coefficients of Table 1
    ac(0) = -0.500004211065677e0;   bc(0) = 0.115956920789028e0
    ac(1) = -0.124989431448700e0;   bc(1) = 0.278919134951974e0
    ac(2) = -0.781685695640975e-2;  bc(2) = 0.252752008621167e-1
    ac(3) = -0.216324415010288e-3;  bc(3) = 0.841879407543506e-3
    ac(4) = -0.344525393452639e-5;  bc(4) = 0.152425102734818e-4
    ac(5) = -0.315133836828774e-7;  bc(5) = 0.148292488399579e-6
    ac(6) = -0.296636186265427e-9;  bc(6) = 0.157622547107156e-8
    ac(7) = -0.313689942474032e-12; bc(7) = 0.117975437124933e-11
    ac(8) = -0.112031912249579e-13; bc(8) = 0.655345107753534e-13

    ! Radius of convergence
    Rconv = 7.0

    ! First lambda is expected to be zero, and is not used
    if ( Naquifers > 1 ) then
        lambda(1:Naquifers-1) = lambdain(2:Naquifers)
    end if

    ! Precompute binomials
    do n=0,Nterms
        do m=0,Nterms
            rbinom(n,m) = product(rRange(m+1:n)) / product(rRange(1:n-m))
        end do
    end do

    zin = cmplx(x,y,8); z1in = cmplx(x1in,y1in,8); z2in = cmplx(x2in,y2in,8)
    Lin = abs(z2in-z1in)
    z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in)
    zplus1 = z + rONE; zmin1 = z - rONE

    ! If at cornerpoint, move slightly
    if  (abs(zplus1) < tiny*rTWO/Lin) zplus1 = zplus1 + tiny
    if  (abs(zmin1) < tiny*rTWO/Lin) zmin1 = zmin1 + tiny

    ! Laplace linesink

    power = order + 1 ;
    pcor = cmplx(0.d0,0.d0) ;
    do n=1,(power+1)/2
        pcor = pcor + z**(power-2*n+1) / (2*n-1) ;
    end do
    pcor = rTWO * pcor ;

    comega = z**power * log( (zmin1) / (zplus1) ) + pcor - log(zmin1) + (-rONE)**power*log(zplus1) ;
    comega = -comega*Lin/(rFOUR*pi*power)

    rv(1) = real(comega)

    ! N-1 leakage factors
    do i=1,Naquifers-1
        pot = 0.0

        ! Check whether entire linesink is outside radius of convergence
        ! Outside if |z-zc|>L/2+7lab, and thus |Z|>1+7lab*2/L, or |zeta|>1/biglab+7 (zeta is called z here)
        biglab = rTWO * lambda(i) / Lin
        z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in) / biglab

        if ( abs(z) < (Rconv + rONE/biglab) ) then
            call IntegralF(zin,z1in,z2in,Lin,lambda(i),order,Nterms,ac,bc,Rconv,rbinom,pot)
            rv(i+1) = -Lin/rTWO * pot
        else
            rv(i+1) = rZERO
        end if

    end do

end subroutine potbeslsho

subroutine disbeslsho(x,y,x1in,y1in,x2in,y2in,Naquifers,lambdain,order,rvx,rvy)

    ! Input:
    !   x,y: Point where discharge is computed
    !   x1in,y1in: Begin point of line-sink
    !   x2in,y2in: End point of line-sink
    !   Naquifers: Number of aquifers
    !   lambdain(Naquifers): Array with zero in first spot and lambda's in remaining spots
    !   order: Order of the line-sink
    !   rvx(Naquifers),rvy(Naquifers): Arrays to store return values
    ! Output:
    !   rvx(Naquifers),rvy(Naquifers): Values of Qx and Qy with Laplace value in first spot
    !                   and mod.Helmholtz potentials in remaining spots

    implicit none
    ! Input / Output
    real*8, intent(in) :: x,y,x1in,y1in,x2in,y2in
    integer, intent(in) :: Naquifers, order
    real*8, intent(in) :: lambdain(*)
    real*8, intent(inout) :: rvx(*), rvy(*)

    ! Locals
    integer, parameter :: Nterms = 8

    real*8, parameter, dimension(0:8) :: rRange = (/0.d0,1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,7.d0,8.d0/)
    real*8, parameter :: rZERO=0.d0, rONE=1.d0, rTWO=2.d0, rFOUR=4.d0
    real*8, parameter :: pi=3.1415926535897931d0, tiny = 1.d-8
    complex*16, parameter :: cZERO=(0.d0,0.d0)

    integer :: n, m, i
    real*8 :: Rconv, Lin, biglab
    complex*16 :: zin, z1in, z2in, z, zplus1, zmin1
    complex*16 :: pcor, wdis, cdum

    real*8, dimension(0:Nterms) :: ac,bc
    real*8, dimension(Naquifers) :: lambda
    real*8, dimension(0:Nterms,0:Nterms) :: rbinom

    ! Coefficients of Table 1
    ac(0) = -0.500004211065677e0;   bc(0) = 0.115956920789028e0
    ac(1) = -0.124989431448700e0;   bc(1) = 0.278919134951974e0
    ac(2) = -0.781685695640975e-2;  bc(2) = 0.252752008621167e-1
    ac(3) = -0.216324415010288e-3;  bc(3) = 0.841879407543506e-3
    ac(4) = -0.344525393452639e-5;  bc(4) = 0.152425102734818e-4
    ac(5) = -0.315133836828774e-7;  bc(5) = 0.148292488399579e-6
    ac(6) = -0.296636186265427e-9;  bc(6) = 0.157622547107156e-8
    ac(7) = -0.313689942474032e-12; bc(7) = 0.117975437124933e-11
    ac(8) = -0.112031912249579e-13; bc(8) = 0.655345107753534e-13

    ! Radius of convergence
    Rconv = 7.0

    ! First lambda is expected to be zero, and is not used
    if ( Naquifers > 1 ) then
        lambda(1:Naquifers-1) = lambdain(2:Naquifers)
    end if

    ! Precompute binomials
    do n=0,Nterms
        do m=0,Nterms
            rbinom(n,m) = product(rRange(m+1:n)) / product(rRange(1:n-m))
        end do
    end do

    zin = cmplx(x,y,8); z1in = cmplx(x1in,y1in,8); z2in = cmplx(x2in,y2in,8)
    Lin = abs(z2in-z1in)
    z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in)
    zplus1 = z + rONE; zmin1 = z - rONE

    ! If at cornerpoint, move slightly
    if  (abs(zplus1) < tiny*rTWO/Lin) zplus1 = zplus1 + tiny
    if  (abs(zmin1) < tiny*rTWO/Lin) zmin1 = zmin1 + tiny

    ! Laplace linesink

    pcor = cZERO ;
    do n=1,(order+1)/2
        pcor = pcor + float(order-2*n+2) * z**(order+1-2*n) / float(2*n-1) ;
    end do
    pcor = rTWO * pcor ;

    cdum = 1.d0 / (order+1)  ! Without this intermediate statement it didn't seem to work
    wdis = float(order+1) * z**order * log( (zmin1) / (zplus1) ) + pcor
    
    wdis = wdis + (z**(order+1) - rONE) / zmin1 - (z**(order+1) - (-rONE)**(order+1)) / zplus1

    wdis = wdis*Lin/rTWO/(z2in-z1in)/pi * cdum;

    rvx(1) = real(wdis);
    rvy(1) = -aimag(wdis)

    ! N-1 leakage factors
    do i=1,Naquifers-1
        wdis = cZERO

        ! Check whether entire linesink is outside radius of convergence
        ! Outside if |z-zc|>L/2+7lab, and thus |Z|>1+7lab*2/L, or |zeta|>1/biglab+7 (zeta is called z here)
        biglab = rTWO * lambda(i) / Lin
        z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in) / biglab

        if ( abs(z) < (Rconv + rONE/biglab) ) then
            call IntegralG(zin,z1in,z2in,Lin,lambda(i),order,Nterms,ac,bc,Rconv,rbinom,wdis)
            wdis = rTWO * Lin / (z2in-z1in) / biglab * wdis

            rvx(i+1) = real(wdis)
            rvy(i+1) = -aimag(wdis)
        else
            rvx(i+1) = rZERO; rvy(i+1) = rZERO
        end if

    end do

end subroutine disbeslsho

subroutine potbesonlylsho(x,y,x1in,y1in,x2in,y2in,Naquifers,lambda,order,rv)

    ! Input:
    !   x,y: Point where potential is computed
    !   x1in,y1in: Begin point of line-sink
    !   x2in,y2in: End point of line-sink
    !   Naquifers: Number of aquifers
    !   lambda(Naquifers-1): Array lambda's (requires confined system)
    !   order: Order of the line-sink
    !   rv(Naquifers-1): Array to store return value
    ! Output:
    !   rv(Naquifers-1): Value of mod.Helmholtz potentials

    implicit none

    ! Input / Output
    real*8, intent(in) :: x,y,x1in,y1in,x2in,y2in
    integer, intent(in) :: Naquifers, order
    real*8, intent(in) :: lambda(*)
    real*8, intent(inout) :: rv(*)

    ! Locals
    integer, parameter :: Nterms = 8

    real*8, parameter, dimension(0:8) :: rRange = (/0.d0,1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,7.d0,8.d0/)
    real*8, parameter :: rZERO=0.d0, rONE=1.d0, rTWO=2.d0, rFOUR=4.d0
    real*8, parameter :: pi=3.1415926535897931d0, tiny = 1.d-8

    integer :: power, n, m, i
    real*8 :: Rconv, Lin, pot, biglab
    complex*16 :: zin, z1in, z2in, z, zplus1, zmin1
    complex*16 :: pcor, comega

    real*8, dimension(0:Nterms) :: ac,bc
    real*8, dimension(0:Nterms,0:Nterms) :: rbinom

    ! Coefficients of Table 1
    ac(0) = -0.500004211065677e0;   bc(0) = 0.115956920789028e0
    ac(1) = -0.124989431448700e0;   bc(1) = 0.278919134951974e0
    ac(2) = -0.781685695640975e-2;  bc(2) = 0.252752008621167e-1
    ac(3) = -0.216324415010288e-3;  bc(3) = 0.841879407543506e-3
    ac(4) = -0.344525393452639e-5;  bc(4) = 0.152425102734818e-4
    ac(5) = -0.315133836828774e-7;  bc(5) = 0.148292488399579e-6
    ac(6) = -0.296636186265427e-9;  bc(6) = 0.157622547107156e-8
    ac(7) = -0.313689942474032e-12; bc(7) = 0.117975437124933e-11
    ac(8) = -0.112031912249579e-13; bc(8) = 0.655345107753534e-13

    ! Radius of convergence
    Rconv = 7.0

    ! Precompute binomials
    do n=0,Nterms
        do m=0,Nterms
            rbinom(n,m) = product(rRange(m+1:n)) / product(rRange(1:n-m))
        end do
    end do

    zin = cmplx(x,y,8); z1in = cmplx(x1in,y1in,8); z2in = cmplx(x2in,y2in,8)
    Lin = abs(z2in-z1in)
    z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in)
    zplus1 = z + rONE; zmin1 = z - rONE

    ! If at cornerpoint, move slightly
    if  (abs(zplus1) < tiny*rTWO/Lin) zplus1 = zplus1 + tiny
    if  (abs(zmin1) < tiny*rTWO/Lin) zmin1 = zmin1 + tiny

    do i=1,Naquifers-1
        pot = 0.0

        ! Check whether entire linesink is outside radius of convergence
        ! Outside if |z-zc|>L/2+7lab, and thus |Z|>1+7lab*2/L, or |zeta|>1/biglab+7 (zeta is called z here)
        biglab = rTWO * lambda(i) / Lin
        z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in) / biglab

        if ( abs(z) < (Rconv + rONE/biglab) ) then
            call IntegralF(zin,z1in,z2in,Lin,lambda(i),order,Nterms,ac,bc,Rconv,rbinom,pot)
            rv(i) = -Lin/rTWO * pot
        else
            rv(i) = rZERO
        end if

    end do

end subroutine potbesonlylsho

subroutine disbesonlylsho(x,y,x1in,y1in,x2in,y2in,Naquifers,lambda,order,rvx,rvy)

    ! Input:
    !   x,y: Point where discharge is computed
    !   x1in,y1in: Begin point of line-sink
    !   x2in,y2in: End point of line-sink
    !   Naquifers: Number of aquifers
    !   lambda(Naquifers-1): Array with lambda's (requires confined system)
    !   order: Order of the line-sink
    !   rvx(Naquifers-1),rvy(Naquifers-1): Arrays to store return values
    ! Output:
    !   rvx(Naquifers-1),rvy(Naquifers-1): Values of Qx and Qy with mod.Helmholtz discharges

    implicit none

    ! Input / Output
    real*8, intent(in) :: x,y,x1in,y1in,x2in,y2in
    integer, intent(in) :: Naquifers, order
    real*8, intent(in) :: lambda(*)
    real*8,intent(inout) :: rvx(*), rvy(*)

    ! Locals
    integer, parameter :: Nterms = 8

    real*8, parameter, dimension(0:8) :: rRange = (/0.d0,1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,7.d0,8.d0/)
    real*8, parameter :: rZERO=0.d0, rONE=1.d0, rTWO=2.d0, rFOUR=4.d0
    real*8, parameter :: pi=3.1415926535897931d0, tiny = 1.d-8
    complex*16, parameter :: cZERO=(0.d0,0.d0)

    integer :: n, m, i
    real*8 :: Rconv, Lin, biglab
    complex*16 :: zin, z1in, z2in, z, zplus1, zmin1
    complex*16 :: pcor, wdis, cdum

    real*8, dimension(0:Nterms) :: ac,bc
    real*8, dimension(0:Nterms,0:Nterms) :: rbinom

    ! Coefficients of Table 1
    ac(0) = -0.500004211065677e0;   bc(0) = 0.115956920789028e0
    ac(1) = -0.124989431448700e0;   bc(1) = 0.278919134951974e0
    ac(2) = -0.781685695640975e-2;  bc(2) = 0.252752008621167e-1
    ac(3) = -0.216324415010288e-3;  bc(3) = 0.841879407543506e-3
    ac(4) = -0.344525393452639e-5;  bc(4) = 0.152425102734818e-4
    ac(5) = -0.315133836828774e-7;  bc(5) = 0.148292488399579e-6
    ac(6) = -0.296636186265427e-9;  bc(6) = 0.157622547107156e-8
    ac(7) = -0.313689942474032e-12; bc(7) = 0.117975437124933e-11
    ac(8) = -0.112031912249579e-13; bc(8) = 0.655345107753534e-13

    ! Radius of convergence
    Rconv = 7.0

    ! Precompute binomials
    do n=0,Nterms
        do m=0,Nterms
            rbinom(n,m) = product(rRange(m+1:n)) / product(rRange(1:n-m))
        end do
    end do

    zin = cmplx(x,y,8); z1in = cmplx(x1in,y1in,8); z2in = cmplx(x2in,y2in,8)
    Lin = abs(z2in-z1in)
    z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in)
    zplus1 = z + rONE; zmin1 = z - rONE

    ! If at cornerpoint, move slightly
    if  (abs(zplus1) < tiny*rTWO/Lin) zplus1 = zplus1 + tiny
    if  (abs(zmin1) < tiny*rTWO/Lin) zmin1 = zmin1 + tiny

    do i=1,Naquifers-1
        wdis = cZERO

        ! Check whether entire linesink is outside radius of convergence
        ! Outside if |z-zc|>L/2+7lab, and thus |Z|>1+7lab*2/L, or |zeta|>1/biglab+7 (zeta is called z here)
        biglab = rTWO * lambda(i) / Lin
        z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in) / biglab

        if ( abs(z) < (Rconv + rONE/biglab) ) then
            call IntegralG(zin,z1in,z2in,Lin,lambda(i),order,Nterms,ac,bc,Rconv,rbinom,wdis)
            wdis = rTWO * Lin / (z2in-z1in) / biglab * wdis

            rvx(i) = real(wdis)
            rvy(i) = -aimag(wdis)
        else
            rvx(i) = rZERO; rvy(i) = rZERO
        end if

    end do

end subroutine disbesonlylsho

subroutine potlaplsho(x,y,x1in,y1in,x2in,y2in,order,rv)

    ! Input:
    !   x,y: Point where potential is computed
    !   x1in,y1in: Begin point of line-sink
    !   x2in,y2in: End point of line-sink
    !   order: Order of the line-sink
    !   rv(order+1): Array to store return value
    ! Output:
    !   rv(order+1): Return values

    implicit none

    ! Input / Output
    real*8, intent(in) :: x,y,x1in,y1in,x2in,y2in
    integer, intent(in) :: order
    real*8,intent(inout) :: rv(*)

    ! Locals

    real*8, parameter :: rONE=1.d0, rTWO=2.d0, rFOUR=4.d0
    real*8, parameter :: pi=3.1415926535897931d0, tiny = 1.d-8
    complex*16, parameter :: cZERO=(0.d0,0.d0), cTWO = (2.d0,0.d0)

    integer :: n, m, i
    real*8 :: Lin
    complex*16 :: zin, z1in, z2in, z, zplus1, zmin1, zpower
    complex*16 :: clog1, clog2, clog3
    complex*16, dimension(order+1) :: comega
    complex*16, dimension(order+2) :: qm, qmnew

    zin = cmplx(x,y,8); z1in = cmplx(x1in,y1in,8); z2in = cmplx(x2in,y2in,8)
    Lin = abs(z2in-z1in)
    z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in)
    zplus1 = z + rONE; zmin1 = z - rONE

    ! If at cornerpoint, move slightly
    if  (abs(zplus1) < tiny*rTWO/Lin) zplus1 = zplus1 + tiny
    if  (abs(zmin1) < tiny*rTWO/Lin) zmin1 = zmin1 + tiny

    ! Laplace linesink

    !qm(1) = cZERO
    !do m=1,order+1
    !   qm(m+1) = cZERO
    !   do n=1,(m+1)/2
    !       qm(m+1) = qm(m+1) + z**(m-rTWO*float(n)+rONE) / (rTWO*float(n) - rONE)
    !   end do
    !end do
    !qm = rTWO * qm
    !write(*,*) 'qm ',qm
    
    qmnew(1) = cZERO
    qmnew(2) = cTWO
    do m=3,order+1,2
        qmnew(m+1) = qmnew(m-1) * z * z + cTWO / float(m)
    end do
    do m=2,order+1,2
        qmnew(m+1) = qmnew(m) * z 
    end do
    !write(*,*) 'qmnew ',qmnew

    clog1 = log( (zmin1) / (zplus1) )
    clog2 = log(zmin1)
    clog3 = log(zplus1)

    zpower = rONE
    do i = 1, order+1
        zpower = zpower * z
        comega(i) = -Lin/(rFOUR*pi*i) * ( zpower * clog1 + qmnew(i+1) - clog2 + (-rONE)**i * clog3 )
        rv(i) = real(comega(i))
    end do
   
end subroutine potlaplsho

subroutine omegalaplsho(x,y,x1in,y1in,x2in,y2in,order,rv,rvpsi)

    ! Input:
    !   x,y: Point where potential is computed
    !   x1in,y1in: Begin point of line-sink
    !   x2in,y2in: End point of line-sink
    !   order: Order of the line-sink
    !   rv(order+1): Array to store return value
    ! Output:
    !   rv(order+1): Return values

    implicit none

    ! Input / Output
    real*8, intent(in) :: x,y,x1in,y1in,x2in,y2in
    integer, intent(in) :: order
    real*8,intent(inout) :: rv(*), rvpsi(*)

    ! Locals

    real*8, parameter :: rONE=1.d0, rTWO=2.d0, rFOUR=4.d0
    real*8, parameter :: pi=3.1415926535897931d0, tiny = 1.d-8
    complex*16, parameter :: cZERO=(0.d0,0.d0), cTWO = (2.d0,0.d0)

    integer :: n, m, i
    real*8 :: Lin
    complex*16 :: zin, z1in, z2in, z, zplus1, zmin1, zpower
    complex*16 :: clog1, clog2, clog3
    complex*16, dimension(order+1) :: comega
    complex*16, dimension(order+2) :: qm, qmnew

    zin = cmplx(x,y,8); z1in = cmplx(x1in,y1in,8); z2in = cmplx(x2in,y2in,8)
    Lin = abs(z2in-z1in)
    z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in)
    zplus1 = z + rONE; zmin1 = z - rONE

    ! If at cornerpoint, move slightly
    if  (abs(zplus1) < tiny*rTWO/Lin) zplus1 = zplus1 + tiny
    if  (abs(zmin1) < tiny*rTWO/Lin) zmin1 = zmin1 + tiny

    ! Laplace linesink

    !qm(1) = cZERO
    !do m=1,order+1
    !   qm(m+1) = cZERO
    !   do n=1,(m+1)/2
    !       qm(m+1) = qm(m+1) + z**(m-rTWO*float(n)+rONE) / (rTWO*float(n) - rONE)
    !   end do
    !end do
    !qm = rTWO * qm
    !write(*,*) 'qm ',qm
    
    qmnew(1) = cZERO
    qmnew(2) = cTWO
    do m=3,order+1,2
        qmnew(m+1) = qmnew(m-1) * z * z + cTWO / float(m)
    end do
    do m=2,order+1,2
        qmnew(m+1) = qmnew(m) * z 
    end do
    !write(*,*) 'qmnew ',qmnew

    clog1 = log( (zmin1) / (zplus1) )
    clog2 = log(zmin1)
    clog3 = log(zplus1)

    zpower = rONE
    do i = 1, order+1
        zpower = zpower * z
        comega(i) = -Lin/(rFOUR*pi*i) * ( zpower * clog1 + qmnew(i+1) - clog2 + (-rONE)**i * clog3 )
        rv(i) = real(comega(i))
        rvpsi(i) = aimag(comega(i))
    end do
   
end subroutine omegalaplsho

subroutine dislaplsho(x,y,x1in,y1in,x2in,y2in,order,rvx,rvy)

    ! Input:
    !   x,y: Point where discharge is computed
    !   x1in,y1in: Begin point of line-sink
    !   x2in,y2in: End point of line-sink
    !   order: Order of the line-sink
    !   rvx(order+1),rvy(order+1): Arrays to store return values
    ! Output:
    !   rvx(order+1),rvy(order+1): Return values

    implicit none

    ! Input / Output
    real*8, intent(in) :: x,y,x1in,y1in,x2in,y2in
    integer, intent(in) :: order
    real*8,intent(inout) :: rvx(*), rvy(*)

    ! Locals
    real*8, parameter :: rZERO=0.d0, rONE=1.d0, rTWO=2.d0, rFOUR=4.d0
    real*8, parameter :: pi=3.1415926535897931d0, tiny = 1.d-8
    complex*16, parameter :: cZERO=(0.d0,0.d0)

    integer :: n, m, i
    real*8 :: Lin
    complex*16 :: zin, z1in, z2in, z, zplus1, zmin1
    complex*16 :: logterm, zterm1, zterm2, zterm3
    complex*16 :: qm, cdum
    complex*16, dimension(order+1) :: wdis

    zin = cmplx(x,y,8); z1in = cmplx(x1in,y1in,8); z2in = cmplx(x2in,y2in,8)
    Lin = abs(z2in-z1in)
    z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in)
    zplus1 = z + rONE; zmin1 = z - rONE

    ! If at cornerpoint, move slightly
    if  (abs(zplus1) < tiny*rTWO/Lin) zplus1 = zplus1 + tiny
    if  (abs(zmin1) < tiny*rTWO/Lin) zmin1 = zmin1 + tiny

    ! Laplace linesink

    ! Laplace line-doublet

    zterm1 = rONE / zmin1
    zterm2 = rONE / zplus1
    zterm3 = zterm1 - zterm2
    logterm = log(zmin1/zplus1)

    do i=1,order+1
        wdis(i) = float(i) * z**(i-1) * logterm
        wdis(i) = wdis(i) + z**i * zterm3 - zterm1 + (-rONE)**i * zterm2
        qm = cZERO
        do n=1,(i+1)/2
            qm = qm + float(i-2*n+1) * z**(i-2*n) / float(2*n-1)
        end do
        wdis(i) = wdis(i) + rTWO * qm
        cdum = 1.d0 / i  ! Without this intermediate statement it didn't seem to work previously, so I leave it in
        wdis(i) = wdis(i) * Lin / rTWO / (z2in-z1in) / pi * cdum;
        rvx(i) = real(wdis(i));
        rvy(i) = -aimag(wdis(i))
    end do

end subroutine dislaplsho

subroutine potbesldho(x,y,x1in,y1in,x2in,y2in,Naquifers,lambdain,order,rv)

    ! Input:
    !   x,y: Point where potential is computed
    !   x1in,y1in: Begin point of line-doublet
    !   x2in,y2in: End point of line-doublet
    !   Naquifers: Number of aquifers
    !   lambdain(Naquifers): Array with zero in first spot and lambda's in remaining spots
    !   order: Order of the line-doublet
    !   rv(Naquifers): Array to store return value
    ! Output:
    !   rv(Naquifers): Value of potentials with Laplace value in first spot
    !                  and mod.Helmholtz potentials in remaining spots

    implicit none

    ! Input / Output
    real*8, intent(in) :: x,y,x1in,y1in,x2in,y2in
    integer, intent(in) :: Naquifers, order
    real*8, intent(in) :: lambdain(*)
    real*8,intent(inout) :: rv(*)

    ! Locals
    integer, parameter :: Nterms = 8

    real*8, parameter, dimension(0:8) :: rRange = (/0.d0,1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,7.d0,8.d0/)
    real*8, parameter :: rZERO=0.d0, rONE=1.d0, rTWO=2.d0, rFOUR=4.d0
    real*8, parameter :: pi=3.1415926535897931d0, tiny = 1.d-8
    complex*16, parameter :: cZERO=(0.d0,0.d0), ci = (0.d0,1.d0)

    integer :: power, n, m, i, m1, m2, NLS
    real*8 :: Rconv, Lin, pot, biglab, del0, ra
    complex*16 :: zin, z1in, z2in, z, zplus1, zmin1
    complex*16 :: pcor, comega, qm, z1, z2

    real*8, dimension(0:Nterms) :: ac,bc
    real*8, dimension(Naquifers) :: lambda
    real*8, dimension(0:Nterms,0:Nterms) :: rbinom

    ! Coefficients of K1
    ac(0) = 0.250000197208863e0;   bc(0) = -0.307966963840932e0  
    ac(1) = 0.312495047266388e-1;  bc(1) = -0.853676915840295e-1 
    ac(2) = 0.130228768540005e-2;  bc(2) = -0.464343185899275e-2 
    ac(3) = 0.270943632576982e-4;  bc(3) = -0.112338260301993e-3 
    ac(4) = 0.341642640876988e-6;  bc(4) = -0.157491472277549e-5 
    ac(5) = 0.271286480571077e-8;  bc(5) = -0.133414321160211e-7 
    ac(6) = 0.197096143802491e-10; bc(6) = -0.106342159633141e-9 
    ac(7) = 0.329351103353054e-13; bc(7) = -0.159531523882074e-12
    ac(8) = 0.573031034976631e-15; bc(8) = -0.340195779923156e-14

    ! Radius of convergence
    Rconv = 7.0

    ! First lambda is expected to be zero, and is not used
    if ( Naquifers > 1 ) then
        lambda(1:Naquifers-1) = lambdain(2:Naquifers)
    end if

    ! Precompute binomials
    do n=0,Nterms
        do m=0,Nterms
            rbinom(n,m) = product(rRange(m+1:n)) / product(rRange(1:n-m))
        end do
    end do

    zin = cmplx(x,y,8); z1in = cmplx(x1in,y1in,8); z2in = cmplx(x2in,y2in,8)
    Lin = abs(z2in-z1in)
    z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in)
    zplus1 = z + rONE; zmin1 = z - rONE

    ! If at cornerpoint, move slightly
    if  (abs(zplus1) < tiny*rTWO/Lin) zplus1 = zplus1 + tiny
    if  (abs(zmin1) < tiny*rTWO/Lin) zmin1 = zmin1 + tiny

    ! Laplace line-doublet
    comega = z**order * log(zmin1/zplus1)
    qm = cZERO
    do n=1,(order+1)/2
        qm = qm + z**(order-rTWO*float(n)+rONE) / (rTWO*float(n) - rONE)
    end do
    comega = rONE/(rTWO*pi*ci) * ( comega + rTWO * qm )

    rv(1) = real(comega)


    ! N-1 leakage factors
    do i=1,Naquifers-1
        pot = rZERO

        ! Check whether entire linedoublet is outside radius of convergence
        ! Outside if |z-zc|>L/2+7lab, and thus |Z|>1+7lab*2/L, or |zeta|>1/biglab+7 (zeta is called z here)
        biglab = rTWO * lambda(i) / Lin
        z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in) / biglab

        if ( abs(z) < (Rconv + rONE/biglab) ) then
            call findm1m2(zin,z1in,z2in,Lin,lambda(i),Rconv,m1,m2,NLS)
            comega = cZERO
            if (m1 > 0) then   ! Otherwise outside radius of convergence
                z1 = z1in + float(m1-1)/float(NLS) * (z2in-z1in)
                z2 = z1in + float(m2)/float(NLS)   * (z2in-z1in)
                del0 = float(1-m1-m2+NLS)/float(1-m1+m2)
                ra = float(NLS) / float(1+m2-m1)
                call IntegralLapLineDipole(zin,z1,z2,del0,ra,order,rbinom,Nterms,comega)
            end if
            call IntegralF(zin,z1in,z2in,Lin,lambda(i),order,Nterms,ac,bc,Rconv,rbinom,pot)
            rv(i+1) = real(comega/ci) + aimag(z) / biglab * pot  ! Note that z is really zeta in analysis
        else
            rv(i+1) = rZERO
        end if

    end do

end subroutine potbesldho

subroutine disbesldho(x,y,x1in,y1in,x2in,y2in,Naquifers,lambdain,order,rvx,rvy)

    ! Input:
    !   x,y: Point where discharge is computed
    !   x1in,y1in: Begin point of line-doublet
    !   x2in,y2in: End point of line-doublet
    !   Naquifers: Number of aquifers
    !   lambdain(Naquifers): Array with zero in first spot and lambda's in remaining spots
    !   order: Order of the line-doublet
    !   rvx(Naquifers),rvy(Naquifers): Arrays to store return values
    ! Output:
    !   rvx(Naquifers),rvy(Naquifers): Values of Qx and Qy with Laplace value in first spot
    !                   and mod.Helmholtz potentials in remaining spots

    implicit none

    ! Input / Output
    real*8, intent(in) :: x,y,x1in,y1in,x2in,y2in
    integer, intent(in) :: Naquifers, order
    real*8, intent(in) :: lambdain(*)
    real*8,intent(inout) :: rvx(*), rvy(*)

    ! Locals
    integer, parameter :: Nterms = 8

    real*8, parameter, dimension(0:8) :: rRange = (/0.d0,1.d0,2.d0,3.d0,4.d0,5.d0,6.d0,7.d0,8.d0/)
    real*8, parameter :: rZERO=0.d0, rONE=1.d0, rTWO=2.d0, rFOUR=4.d0
    real*8, parameter :: pi=3.1415926535897931d0, tiny = 1.d-8
    complex*16, parameter :: cZERO=(0.d0,0.d0), ci = (0.d0,1.d0)

    integer :: n, m, i, m1, m2, NLS
    real*8 :: Rconv, Lin, biglab, del0, ra, pot
    complex*16 :: zin, z1in, z2in, z, z1, z2, zplus1, zmin1
    complex*16 :: qm, wdis, wdis1, wdis2, wdis3

    real*8, dimension(0:Nterms) :: ac,bc
    real*8, dimension(Naquifers) :: lambda
    real*8, dimension(0:Nterms,0:Nterms) :: rbinom

    ! Coefficients of K1
    ac(0) = 0.250000197208863e0;   bc(0) = -0.307966963840932e0  
    ac(1) = 0.312495047266388e-1;  bc(1) = -0.853676915840295e-1 
    ac(2) = 0.130228768540005e-2;  bc(2) = -0.464343185899275e-2 
    ac(3) = 0.270943632576982e-4;  bc(3) = -0.112338260301993e-3 
    ac(4) = 0.341642640876988e-6;  bc(4) = -0.157491472277549e-5 
    ac(5) = 0.271286480571077e-8;  bc(5) = -0.133414321160211e-7 
    ac(6) = 0.197096143802491e-10; bc(6) = -0.106342159633141e-9 
    ac(7) = 0.329351103353054e-13; bc(7) = -0.159531523882074e-12
    ac(8) = 0.573031034976631e-15; bc(8) = -0.340195779923156e-14

    ! Radius of convergence
    Rconv = 7.0

    ! First lambda is expected to be zero, and is not used
    if ( Naquifers > 1 ) then
        lambda(1:Naquifers-1) = lambdain(2:Naquifers)
    end if

    ! Precompute binomials
    do n=0,Nterms
        do m=0,Nterms
            rbinom(n,m) = product(rRange(m+1:n)) / product(rRange(1:n-m))
        end do
    end do

    zin = cmplx(x,y,8); z1in = cmplx(x1in,y1in,8); z2in = cmplx(x2in,y2in,8)
    Lin = abs(z2in-z1in)
    z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in)
    zplus1 = z + rONE; zmin1 = z - rONE

    ! If at cornerpoint, move slightly
    if  (abs(zplus1) < tiny*rTWO/Lin) zplus1 = zplus1 + tiny
    if  (abs(zmin1) < tiny*rTWO/Lin) zmin1 = zmin1 + tiny

    ! Laplace line-doublet

    if (order == 0) then
        wdis = - ( rONE / zmin1 - rONE / zplus1 ) / (pi * ci * (z2in - z1in) )
    else
        wdis = float(order) * z**(order-1) * log(zmin1/zplus1)
        wdis = wdis + z**order * ( rONE / zmin1 - rONE / zplus1 )
        qm = cZERO
        if (order > 1) then ! To avoid a possible problem of 0 * 0^(-1)
            do n=1,order/2
                qm = qm + float(order-2*n+1) * z**(order-2*n) / float(2*n-1)
            end do
        end if
        wdis = - ( wdis + rTWO * qm ) / (pi*ci*(z2in-z1in)) ;
    end if

    rvx(1) = real(wdis); rvy(1) = -aimag(wdis)

    ! N-1 leakage factors
    do i=1,Naquifers-1
        wdis = cZERO

        ! Check whether entire line-doublet is outside radius of convergence
        ! Outside if |z-zc|>L/2+7lab, and thus |Z|>1+7lab*2/L, or |zeta|>1/biglab+7 (zeta is called z here)
        biglab = rTWO * lambda(i) / Lin
        z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in) / biglab

        if ( abs(z) < (Rconv + rONE/biglab) ) then
            call findm1m2(zin,z1in,z2in,Lin,lambda(i),Rconv,m1,m2,NLS)
            wdis1 = cZERO
            if (m1 > 0) then
                z1 = z1in + float(m1-1)/float(NLS) * (z2in-z1in)
                z2 = z1in + float(m2)/float(NLS)   * (z2in-z1in)
                del0 = float(1-m1-m2+NLS)/float(1-m1+m2)
                ra = float(NLS) / float(1+m2-m1)
                call IntegralLapLineDipoleDis(zin,z1,z2,del0,ra,order,rbinom,Nterms,wdis1)
                wdis1 = -rTWO  * wdis1 / (ci*(z2-z1))
            end if
            call IntegralF(zin,z1in,z2in,Lin,lambda(i),order,Nterms,ac,bc,Rconv,rbinom,pot)
            call IntegralG(zin,z1in,z2in,Lin,lambda(i),order,Nterms,ac,bc,Rconv,rbinom,wdis2)
            wdis3 =  pot / (rTWO * ci) + wdis2 * aimag(z)
            
            wdis = wdis1 - rFOUR * wdis3 / ( biglab**2 * (z2in-z1in) )
            
        end if 

        rvx(i+1) = real(wdis); rvy(i+1) = -aimag(wdis)

    end do

end subroutine disbesldho

subroutine potlapldho(x,y,x1in,y1in,x2in,y2in,order,rv)

    ! Input:
    !   x,y: Point where potential is computed
    !   x1in,y1in: Begin point of line-doublet
    !   x2in,y2in: End point of line-doublet
    !   order: Order of the line-doublet
    !   rv(order+1): Array to store return value
    ! Output:
    !   rv(order+1): Return values

    implicit none

    ! Input / Output
    real*8, intent(in) :: x,y,x1in,y1in,x2in,y2in
    integer, intent(in) :: order
    real*8,intent(inout) :: rv(*)

    ! Locals
    real*8, parameter :: rZERO=0.d0, rONE=1.d0, rTWO=2.d0, rFOUR=4.d0
    real*8, parameter :: pi=3.1415926535897931d0, tiny = 1.d-8
    integer, parameter :: Nfarfield = 8 ! means minimally 8/2 = 4 non-zero terms
    complex*16, parameter :: cZERO=(0.d0,0.d0), cTWO = (2.d0,0.d0), ci = (0.d0,1.d0)

    integer :: power, n, m, i
    real*8 :: Lin, pot, biglab, del0, ra, Rfarfield
    complex*16 :: zin, z1in, z2in, z, zplus1, zmin1
    complex*16 :: z1, z2, oneoverZsq, oneoverZpower, Zpower
    complex*16, dimension(order+1) :: comega, qm, qmnew
    complex*16, dimension(0:order+Nfarfield) :: series

    Rfarfield = 5.0d0 
    zin = cmplx(x,y,8); z1in = cmplx(x1in,y1in,8); z2in = cmplx(x2in,y2in,8)
    z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in)

    if ( abs(z) < Rfarfield ) then

        zplus1 = z + rONE; zmin1 = z - rONE

        ! If at cornerpoint, move slightly
        Lin = abs(z2in-z1in)
        if  (abs(zplus1) < tiny*rTWO/Lin) zplus1 = zplus1 + tiny
        if  (abs(zmin1) < tiny*rTWO/Lin) zmin1 = zmin1 + tiny

        comega(1) = log(zmin1/zplus1)
        do n=1,order
            comega(n+1) = z * comega(n)
        end do
        
        qmnew(1) = cZERO
        if (order > 0) qmnew(2) = cTWO
        do m=3,order,2
            qmnew(m+1) = qmnew(m-1) * z * z + cTWO / float(m)
        end do
        do m=2,order,2
            qmnew(m+1) = qmnew(m) * z 
        end do

        comega = rONE/(rTWO*pi*ci) * ( comega + qmnew )

    else

        series = cZERO
        oneoverZsq = 1.0d0/ ( z * z )
        oneoverZpower = 1.0 / z
        do i = 1, order + Nfarfield, 2
            series(i) = oneoverZpower / float(i)
            oneoverZpower = oneoverZpower * oneoverZsq
        end do
        Zpower = 1.0
        do i = 0, order
            comega(i+1) = sum( series(i:order+Nfarfield) ) * Zpower
            Zpower = Zpower * z
        end do
        comega = -rONE / ( pi * ci ) * comega

    end if

    do i = 1, order+1
        rv(i) = real( comega(i) )
    end do

end subroutine potlapldho

subroutine dislapldho(x,y,x1in,y1in,x2in,y2in,order,rvx,rvy)

    ! Input:
    !   x,y: Point where discharge is computed
    !   x1in,y1in: Begin point of line-doublet
    !   x2in,y2in: End point of line-doublet
    !   order: Order of the line-doublet
    !   rvx(order+1),rvy(order+1): Arrays to store return values
    ! Output:
    !   rvx(order+1),rvy(order+1): Return values

    implicit none

    ! Input / Output
    real*8, intent(in) :: x,y,x1in,y1in,x2in,y2in
    integer, intent(in) :: order
    real*8,intent(inout) :: rvx(*), rvy(*)

    real*8, parameter :: rZERO=0.d0, rONE=1.d0, rTWO=2.d0, rFOUR=4.d0
    real*8, parameter :: pi=3.1415926535897931d0, tiny = 1.d-8
    integer, parameter :: Nfarfield = 8 ! means minimally 8/2 = 4 non-zero terms
    complex*16, parameter :: cZERO=(0.d0,0.d0), ci = (0.d0,1.d0)

    integer :: n, m, i
    real*8 :: Lin, Rfarfield
    complex*16 :: zin, z1in, z2in, z, z1, z2, zplus1, zmin1, logterm, zterm, qm
    complex*16 :: oneoverZsq, oneoverZpower, Zpower
    complex*16, dimension(order+1) :: wdis
    complex*16, dimension(0:order+Nfarfield) :: series, nseries

    Rfarfield = 5.0d0 
    zin = cmplx(x,y,8); z1in = cmplx(x1in,y1in,8); z2in = cmplx(x2in,y2in,8)
    z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in)

    if ( abs(z) < Rfarfield ) then

        zplus1 = z + rONE; zmin1 = z - rONE

        ! If at cornerpoint, move slightly
        Lin = abs(z2in-z1in)
        if  (abs(zplus1) < tiny*rTWO/Lin) zplus1 = zplus1 + tiny
        if  (abs(zmin1) < tiny*rTWO/Lin) zmin1 = zmin1 + tiny

        zterm = rONE / zmin1 - rONE / zplus1 
        wdis(1) = - zterm / (pi * ci * (z2in - z1in) )
        if (order > 0) then
            logterm = log(zmin1/zplus1)
        end if
        do i=1,order
            wdis(i+1) = float(i) * z**(i-1) * logterm
            wdis(i+1) = wdis(i+1) + z**i * zterm
            qm = cZERO
            if (i > 1) then ! To avoid a possible problem of 0 * 0^(-1)
                do n=1,i/2
                    qm = qm + float(i-2*n+1) * z**(i-2*n) / float(2*n-1)
                end do
            end if
            wdis(i+1) = - ( wdis(i+1) + rTWO * qm ) / (pi*ci*(z2in-z1in)) ;
        end do

    else

        series = cZERO
        nseries = cZERO
        oneoverZsq = 1.0d0/ ( z * z )
        oneoverZpower = 1.0 / z
        do i = 1, order + Nfarfield, 2
            series(i) = oneoverZpower / float(i)
            nseries(i) = series(i) * float(i)
            oneoverZpower = oneoverZpower * oneoverZsq
        end do
        Zpower = 1.0 / z
        do i = 0, order
            wdis(i+1) = sum( nseries(i:order+Nfarfield) - float(i)*series(i:order+Nfarfield) ) * Zpower
            Zpower = Zpower * z
        end do
        wdis = -rONE / ( pi * ci ) * wdis * rTWO / (z2in - z1in)

    endif

    do i = 1, order+1
        rvx(i) = real( wdis(i) )
        rvy(i) = -aimag( wdis(i) )
    end do

end subroutine dislapldho

subroutine IntegralF(zin,z1in,z2in,Lin,lambda,order,Nterms,ac,bc,Rconv,rbinom,pot)

    implicit none
    ! Input
    integer :: order, Nterms
    real*8 :: Lin, lambda, Rconv
    complex*16 :: zin, z1in, z2in
    real*8, dimension(0:Nterms) :: ac,bc
    real*8, dimension(0:Nterms,0:Nterms) :: rbinom

    ! Out
    real*8 :: pot

    ! Local
    integer :: NLS,m1,m2,j, m, n
    real*8 :: L, biglab, del1, del2
    complex*16 :: z, zbar, z1, z2, cInt, cd1minz, cd2minz, cln1, cln2
    complex*16, dimension(0:Nterms) :: czmzbarp
    complex*16, dimension(0:Nterms,0:Nterms) :: cgamma
    complex*16, dimension(0:2*Nterms) :: calphat, cbetat
    complex*16, dimension(0:order+1) :: cc
    complex*16, dimension(0:2*Nterms+order) :: calpha, cbeta
    real*8, parameter :: rZERO=0.d0, rONE=1.d0, rTWO=2.d0, rFOUR=4.d0
    real*8, parameter :: pi=3.1415926535897931d0, tiny = 1.d-8

    call findm1m2(zin,z1in,z2in,Lin,lambda,Rconv,m1,m2,NLS)
    if ( m1 == 0 ) then
        pot = rZERO
        return
    end if

    ! Compute zeta (called z here). This is the regular value of the entire element
    L = abs(z2in-z1in)
    biglab = rTWO * lambda / L
    z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in) / biglab; zbar=conjg(z);

    ! Coefficients gamma(n,m), Eq. 21
    ! Store coefficents in matrix.
    do n = 0, Nterms
        czmzbarp(n) = (z-zbar)**n
    end do
    do n = 0, Nterms
        do m = 0, n
            cgamma(n,m) = rbinom(n,m) * czmzbarp(n-m)
        end do
    end do

    ! Eq. 23; These coefficients should be modified for a higher order linesink
    do n = 0, 2*Nterms
        calphat(n) = cmplx(0.d0,0.d0); cbetat(n)=cmplx(0.d0,0.d0)
        do m = max(0,n-Nterms), n/2
            calphat(n) = calphat(n) + ac(n-m) * cgamma(n-m,m)
            cbetat(n) = cbetat(n) + bc(n-m) * cgamma(n-m,m)
        end do
    end do

    ! Compute coefficients of delta^p
    do m = 0,order
        cc(m) = rbinom(order,m) * z**(order-m) * biglab**order
    end do
    if ( order > 0 ) then
        do n=0,2*Nterms+order
            calpha(n) = cmplx(0.d0,0.d0); cbeta(n) = cmplx(0.d0,0.d0)
            do m = max(0,n-2*Nterms), min(n,order)
                calpha(n) = calpha(n) + cc(m) * calphat(n-m)
                cbeta(n) = cbeta(n) + cc(m) * cbetat(n-m)
            end do
        end do
    else
        calpha = calphat
        cbeta = cbetat
    end if

    ! Evaluation of integral, Eq. 25
    cInt = cmplx(0.d0,0.d0)
    del1 = -rONE + rTWO * (float(m1)-rONE) / float(NLS)
    del2 = -rONE + rTWO * float(m2) / float(NLS)
    cd1minz = del1 / biglab - z;  cd2minz = del2 / biglab - z
    if ( abs(cd1minz) < tiny/lambda) cd1minz = cd1minz + tiny
    if ( abs(cd2minz) < tiny/lambda) cd2minz = cd2minz + tiny
    cln1 = log(cd1minz);  cln2 = log(cd2minz)
    do n=0,2*Nterms + order
        cInt = cInt + ( rTWO * calpha(n) * cln2 - rTWO * calpha(n)/(n+1) + cbeta(n) ) * (cd2minz)**(n+1) / float(n+1)
        cInt = cInt - ( rTWO * calpha(n) * cln1 - rTWO * calpha(n)/(n+1) + cbeta(n) ) * (cd1minz)**(n+1) / float(n+1)
    end do
    pot = real(cInt) * biglab / (rTWO*pi)

end subroutine IntegralF

subroutine IntegralG(zin,z1in,z2in,Lin,lambda,order,Nterms,ac,bc,Rconv,rbinom,wdis)

    implicit none
    ! Input
    integer :: order, Nterms
    real*8 :: Lin, lambda, Rconv
    complex*16 :: zin, z1in, z2in
    real*8, dimension(0:Nterms) :: ac,bc
    real*8, dimension(0:Nterms,0:Nterms) :: rbinom

    ! Out
    complex*16 :: wdis

    ! Local
    integer :: NLS,m1,m2,j, m, n
    real*8 :: L, biglab, del1, del2, del0, ra, biglabin
    complex*16 :: z, zbar, z1, z2, cInt, cd1minz, cd2minz, cln1, cln2
    complex*16 :: g1, g2, g3, comega
    complex*16, dimension(0:Nterms) :: czmzbarp
    complex*16, dimension(0:Nterms,0:Nterms) :: cgamma
    complex*16, dimension(0:2*(Nterms-1)) :: cahat, cbhat
    complex*16, dimension(0:2*Nterms) :: calphat, cbetat
    complex*16, dimension(0:order+1) :: cc
    complex*16, dimension(0:2*Nterms+order) :: calpha, cbeta
    real*8, parameter :: rZERO=0.d0, rONE=1.d0, rTWO=2.d0, rFOUR=4.d0
    real*8, parameter :: pi=3.1415926535897931d0, tiny = 1.d-8
    complex*16, parameter :: cZERO=(0.d0,0.d0)

    biglabin = rTWO * lambda / Lin

    call findm1m2(zin,z1in,z2in,Lin,lambda,Rconv,m1,m2,NLS)
    if ( m1 == 0 ) then
        wdis = cZERO
        return
    end if
!    write(*,*)'m1,m2 ',m1,m2

    ! Compute zeta (called z here). This is the regular value of the entire element
    L = abs(z2in-z1in)
    biglab = rTWO * lambda / L
    z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in) / biglab; zbar=conjg(z);

    ! Coefficients gamma(n,m), Eq. 21
    ! Store coefficents in matrix.
    do n = 0, Nterms
        czmzbarp(n) = (z-zbar)**n
    end do
    do n = 0, Nterms
    do m = 0, n
        cgamma(n,m) = rbinom(n,m) * czmzbarp(n-m)
    end do
    end do

    ! Integral g1
    ! Implemented with different z, rather than Delta1 and Delta2

    z1 = z1in + float(m1-1)/float(NLS) * (z2in-z1in)
    z2 = z1in + float(m2)/float(NLS)   * (z2in-z1in)
    del0 = float(1-m1-m2+NLS)/float(1-m1+m2)
    ra = float(NLS) / float(1+m2-m1)
    comega = cZERO
    call IntegralLapLineDipole(zin,z1,z2,del0,ra,order,rbinom,Nterms,comega)

    g1 = -ac(0) * biglabin * comega

    ! Integral g2

    ! Compute hat coefficients
    do n = 0, Nterms-1
        cahat(n) = float(n+1) * ac(n+1)
        cbhat(n) = ac(n+1) + float(n+1) * bc(n+1)
    end do

    ! Eq. 23
    do n = 0, 2*Nterms-1
        calphat(n) = cZERO; cbetat(n) = cZERO   
        do m = max(0,n-Nterms+1),(n+1)/2
            calphat(n) = calphat(n) + cahat(n-m) * cgamma(n-m+1,m)
            cbetat(n) = cbetat(n) + cbhat(n-m) * cgamma(n-m+1,m)
        end do
    end do

    ! Compute coefficients of delta^p
    do m = 0,order
        cc(m) = rbinom(order,m) * z**(order-m) * biglab**order
    end do
    if ( order > 0 ) then
        do n = 0, 2*Nterms - 1 + order
            calpha(n) = cZERO; cbeta(n) = cZERO
            do m = max(0,n-2*Nterms+1), min(n,order)
                calpha(n) = calpha(n) + cc(m) * calphat(n-m)
                cbeta(n) = cbeta(n) + cc(m) * cbetat(n-m)
            end do
        end do
    else
        calpha = calphat
        cbeta = cbetat
    end if

    ! Computation of integral
    g2 = cZERO
    del1 = -rONE + rTWO * (float(m1)-rONE) / float(NLS)
    del2 = -rONE + rTWO * float(m2) / float(NLS)
    cd1minz = del1 / biglab - z;  cd2minz = del2 / biglab - z
    if ( abs(cd1minz) < tiny/lambda) cd1minz = cd1minz + tiny
    if ( abs(cd2minz) < tiny/lambda) cd2minz = cd2minz + tiny
    cln1 = log(cd1minz);  cln2 = log(cd2minz)
    do n = 0, 2*Nterms - 1 + order
        g2 = g2 - ( calpha(n) * cln2 - calpha(n)/(n+1) + cbeta(n) ) * (cd2minz)**(n+1) / float(n+1)
        g2 = g2 + ( calpha(n) * cln1 - calpha(n)/(n+1) + cbeta(n) ) * (cd1minz)**(n+1) / float(n+1)
    end do
    g2 = biglabin * g2 / (2*pi)

    ! Integral g3

    ! Eq. 23
    calphat(0) = cZERO
    do n = 1, 2*Nterms-1  ! Loop start at 1, because of bug in Digital Fortran
        calphat(n) = cZERO
        do m = max(0,n-Nterms), (n-1)/2
            calphat(n) = calphat(n) + cahat(n-m-1) * cgamma(n-m-1,m) * (-rONE)**(n-1-2*m)
        end do
    end do

    ! Compute coefficients of delta^p 
    do m = 0,order
        cc(m) = rbinom(order,m) * zbar**(order-m) * biglab**order
    end do
    if ( order > 0 ) then
        do n = 0, 2*Nterms-1+order
            calpha(n) = cZERO
            do m = max(0,n-2*Nterms+1), min(n,order)
                calpha(n) = calpha(n) + cc(m) * calphat(n-m)
            end do
        end do
    else
        calpha = calphat
    end if

    ! Computation of integral
    g3 = cZERO
    ! cd1minz = del1 / biglab - zbar;  cd2minz = del2 / biglab - zbar
    ! if ( abs(cd1minz) < tiny) cd1minz = cd1minz + tiny
    ! if ( abs(cd2minz) < tiny) cd2minz = cd2minz + tiny
    ! By definition log is conjugate of previous log; this avoids problems with signs along the line (and saves logs).
    cd1minz = conjg(cd1minz); cd2minz = conjg(cd2minz)
    cln1 = conjg(cln1);  cln2 = conjg(cln2)
    do n = 0, 2*Nterms - 1 + order
        g3 = g3 - ( calpha(n) * cln2 - calpha(n)/(n+1) ) * (cd2minz)**(n+1) / float(n+1)
        g3 = g3 + ( calpha(n) * cln1 - calpha(n)/(n+1) ) * (cd1minz)**(n+1) / float(n+1)
    end do
    g3 = biglabin * g3 / (rTWO*pi)

    wdis = g1 + g2 + g3

end subroutine IntegralG

subroutine findm1m2(zin,z1in,z2in,Lin,lambda,Rconv,m1,m2,NLS)

    implicit none
    ! Input
    integer :: m1, m2
    real*8 :: Lin, lambda, Rconv
    complex*16 :: zin, z1in, z2in

    ! Local
    integer :: NLS, j
    real*8 :: biglab, L
    complex*16 :: z, z1, z2
    real*8, parameter :: rTWO=2.d0
    
    ! Break integral up in sections of max one lambda
    ! and find first (m1) and last (m2) section within radius of convergence
    NLS = ceiling(Lin/lambda)
    m1 = 0; m2 = 0
    do j = 1, NLS
        z1 = z1in + float(j-1)/NLS * (z2in-z1in)
        z2 = z1 + (z2in-z1in)/NLS
        L = abs(z2-z1)
        biglab = rTWO * lambda / L
        z = ( rTWO*zin - (z1 + z2) ) / (z2 - z1) / biglab
        if ( m1 == 0 ) then
            if ( abs(z) < Rconv ) m1 = j
        else
            if ( abs(z) > Rconv ) then
                m2 = j-1
                exit
            end if
        end if
    end do
    if ( m2 == 0 ) m2 = NLS

end subroutine findm1m2

subroutine IntegralLapLineDipole(zin,z1,z2,del0,ra,order,rbinom,Nterms,comega)

    implicit none
    ! Input
    integer :: order, Nterms
    real*8 :: del0,ra
    real*8, dimension(0:Nterms,0:Nterms) :: rbinom
    complex*16 :: zin, z1in, z2in, comega

    ! Local
    integer :: m, n
    real*8 :: L
    complex*16 :: z, z1, z2, zmin1, zplus1, zterm, qm, qmtot
    complex*16, dimension(0:order+1) :: cg
    real*8, parameter :: rZERO=0.d0, rONE=1.d0, rTWO=2.d0
    real*8, parameter :: pi=3.1415926535897931d0, tiny = 1.d-8
    complex*16, parameter :: cZERO=(0.d0,0.d0)

    z = ( rTWO*zin - (z1 + z2) ) / (z2 - z1)
    zplus1 = z + rONE; zmin1 = z - rONE
    L = abs(z2-z1)

    ! If at cornerpoint, move slightly
    if  (abs(zplus1) < tiny*rTWO/L) zplus1 = zplus1 + tiny
    if  (abs(zmin1) < tiny*rTWO/L) zmin1 = zmin1 + tiny

    ! We don't always have to do this, so maybe put in condition?
    ! Determine coefficients of powers of Delta

    do m = 0,order
        cg(m) = rbinom(order,m) * (-del0 )**(order-m) /ra**order
    end do

    zterm = cZERO
    do n = 0,order
        zterm = zterm + cg(n) * z**n
    enddo
    
    qmtot = cZERO
    do m = 1,order
        qm = cZERO
        do n=1,(m+1)/2
            qm = qm + z**(m-2*n+1) / float(2*n - 1)
        end do
        qmtot = qmtot + rTWO * cg(m) * qm
    end do

    comega = ( zterm * log(zmin1/zplus1) + qmtot ) / (rTWO * pi)
    
end subroutine IntegralLapLineDipole

subroutine IntegralLapLineDipoleDis(zin,z1,z2,del0,ra,order,rbinom,Nterms,wdis)

    implicit none
    ! Input
    integer :: order, Nterms
    real*8 :: del0,ra
    real*8, dimension(0:Nterms,0:Nterms) :: rbinom
    complex*16 :: zin, z1in, z2in, wdis

    ! Local
    integer :: m, n
    real*8 :: L
    complex*16 :: z, z1, z2, zmin1, zplus1, zterm1, zterm2, qm, qmtot
    complex*16, dimension(0:order+1) :: cg
    real*8, parameter :: rZERO=0.d0, rONE=1.d0, rTWO=2.d0
    real*8, parameter :: pi=3.1415926535897931d0, tiny = 1.d-8
    complex*16, parameter :: cZERO=(0.d0,0.d0)

    z = ( rTWO*zin - (z1 + z2) ) / (z2 - z1)
    zplus1 = z + rONE; zmin1 = z - rONE
    L = abs(z2-z1)

    ! If at cornerpoint, move slightly
    if  (abs(zplus1) < tiny*rTWO/L) zplus1 = zplus1 + tiny
    if  (abs(zmin1) < tiny*rTWO/L) zmin1 = zmin1 + tiny

    ! Determine coefficients of powers of Delta for [ (Delta-Delta_0)/a ] ^ p 
    do m = 0,order
        cg(m) = rbinom(order,m) * (-del0 )**(order-m) /ra**order
    end do

    zterm1 = cZERO; zterm2 = cZERO
    do n = 1,order
        zterm1 = zterm1 + cg(n) * float(n) * z**(n-1)
    enddo
    do n = 0,order
        zterm2 = zterm2 + cg(n) * z**n
    end do
    
    qmtot = cZERO
    do m = 2,order
        qm = cZERO
        do n=1,m/2
            qm = qm + float(m-2*n+1) * z**(m-2*n) / float(2*n - 1)
        end do
        qmtot = qmtot + rTWO * cg(m) * qm
    end do

    wdis = ( zterm1 * log(zmin1/zplus1) + zterm2 * (rONE/zmin1 - rONE/zplus1) + qmtot ) / (rTWO * pi)
    
end subroutine IntegralLapLineDipoleDis

!!! Program used for testing
Program BesselLS
  integer :: Naquifers,order
  real*8 :: x,y,x1in,y1in,x2in,y2in,lambdain(5),xright,ytop,del,disxnum,disynum
  real*8 :: rvpot(5),rvdisx(5),rvdisy(5),rvpotright(5),rvpottop(5),disx(5),disy(5)
  real*8 :: rvbes(4),lambda(4),rvxbes(4),rvybes(4)
!!  write(*,*) 'Naquifers,del '
!!  read (*,*) Naquifers,del
!  Naquifers = 1
!  del = 0.0001d0
  write(*,*) 'hello'
  write(*,*) 'x,y,x1,y1,x2,y2,order '
  read (*,*) x,y,x1in,y1in,x2in,y2in,order
!  write(*,*) 'x,y,x1,y1,x2,y2,order '
!  read (*,*) x,y,x1in,y1in,x2in,y2in,order
!  lambdain(1)=0.0
!  xright = x + del
!  ytop = y + del
  Naquifers = 5
  lambdain(1) = 0.0
  lambdain(2) = 1.0
  lambdain(3) = 2.0
  lambdain(4) = 5.0
  lambdain(5) = 10.0
!  rvpot = lambdain * lambdain
!  write(*,*) 'rvpot ',rvpot
!  call dislapldho(x,y,x1in,y1in,x2in,y2in,order,rvdisx,rvdisy)
!  write(*,*) 'disx ',rvdisx
!  write(*,*) 'disy ',rvdisy
  lambda(1:4) = lambdain(2:5)
!  call potbesonlylsho(x,y,x1in,y1in,x2in,y2in,Naquifers,lambda,order,rvbes)
!  write(*,*) 'new ',rvbes
!  write(*,*) 'qx ',rvdisx
!  write(*,*) 'qy ',rvdisy
!  do i = 0,order
!      call disbeslsho(x,y,x1in,y1in,x2in,y2in,Naquifers,lambdain,i,rvdisx,rvdisy)
!      write(*,*)'i,qx,qy ',i,rvdisx,rvdisy
!      call disbesonlylsho(x,y,x1in,y1in,x2in,y2in,Naquifers,lambdain,i,rvxbes,rvybes)
!      write(*,*)'inew    ',i,rvxbes,rvybes
!  end do
!    write(*,*)' pot ',rvpot(1),rvpot(2)
!  Naquifers = 1
!  lambdain(1)=0.0
!  call disbeslsho(x,y,x1in,y1in,x2in,y2in,Naquifers,lambdain,order,rvdisx,rvdisy)
!  write(*,*) 'rvdisx,rvdisy   ',rvdisx(1),rvdisy(1)
!  call dislaplsho(x,y,x1in,y1in,x2in,y2in,order,disx,disy)
!  write(*,*) 'disx ',disx
!  write(*,*) 'disy ',disy
!  call potbeslsho(xright,y,x1in,y1in,x2in,y2in,Naquifers,lambdain,order,rvpotright)
!  call potbeslsho(x,ytop,x1in,y1in,x2in,y2in,Naquifers,lambdain,order,rvpottop)
!  call disbeslsho(x,y,x1in,y1in,x2in,y2in,Naquifers,lambdain,order,rvdisx,rvdisy)
!  disx = (rvpot - rvpotright) / del
!  disy = (rvpot - rvpottop) / del
!  write(*,*) 'qx    ',rvdisx(1),rvdisx(2)
!  write(*,*) 'qxnum ',disx(1),disx(2)
!  write(*,*) 'qy    ',rvdisy(1),rvdisy(2)
!  write(*,*) 'qynum ',disy(1),disy(2)
!  write(*,*) 'pot ',rvpot(1),rvpot(2)
end
