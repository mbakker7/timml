module besselaesnew

    integer :: Nterms
    real(kind=8), dimension(0:8) :: rRange
    real(kind=8) :: rZERO, rONE, rTWO, rFOUR, pi, tiny
    real(kind=8), dimension(0:8,0:8) :: rbinom
    real(kind=8), dimension(0:8) :: ac, bc, ac1, bc1
    complex(kind=8) :: cZERO, ci

contains

    subroutine initialize
    
        implicit none
        integer :: n, m
    
        Nterms = 8
        
        do n = 0, 8
            rRange(n) = real(n)
        end do
        
        rZERO=0.d0
        cZERO=(0.d0,0.d0)
        ci = (0.d0, 1.d0)
        rONE=1.d0
        rTWO=2.d0
        rFOUR=4.d0
        pi=3.1415926535897931d0
        tiny = 1.d-8
        
        do n=0,Nterms
            do m=0,Nterms
                rbinom(n,m) = product(rRange(m+1:n)) / product(rRange(1:n-m))
            end do
        end do
    
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
        
        ! Coefficients of K1
        ac1(0) = 0.250000197208863e0;   bc1(0) = -0.307966963840932e0  
        ac1(1) = 0.312495047266388e-1;  bc1(1) = -0.853676915840295e-1 
        ac1(2) = 0.130228768540005e-2;  bc1(2) = -0.464343185899275e-2 
        ac1(3) = 0.270943632576982e-4;  bc1(3) = -0.112338260301993e-3 
        ac1(4) = 0.341642640876988e-6;  bc1(4) = -0.157491472277549e-5 
        ac1(5) = 0.271286480571077e-8;  bc1(5) = -0.133414321160211e-7 
        ac1(6) = 0.197096143802491e-10; bc1(6) = -0.106342159633141e-9 
        ac1(7) = 0.329351103353054e-13; bc1(7) = -0.159531523882074e-12
        ac1(8) = 0.573031034976631e-15; bc1(8) = -0.340195779923156e-14
        
    end subroutine initialize
    
    function potbeslsho(x, y, z1, z2, lambda, order, ilap, naq) result (rv)
    
        ! Input:
        !   x,y: Point where potential is computed
        !   z1: Complex begin point of line-sink
        !   z2: Complex end point of line-sink
        !   lambda(naq): lambda's (zero for first lambda if Laplace)
        !   order: Order of the line-sink
        !   ilap: equals 1 when first value is Laplace line-sink and first lambda equals zero
        !   naq: Number of aquifers
        !   rv(naq): Array to store return value (must be pre-allocated)
        ! Output:
        !   rv(naq): Potentials. Fist spot is Laplace value if ilap=1
    
        implicit none
        ! Input / Output
        real(kind=8), intent(in) :: x, y
        complex(kind=8), intent(in) :: z1, z2
        integer, intent(in) :: naq, order, ilap
        real(kind=8), dimension(naq), intent(in) :: lambda
        real(kind=8), dimension(naq) :: rv
    
        ! Locals
        integer, parameter :: Nterms = 8
        integer :: power, n, i, lstype
        real(kind=8) :: Rconv, Lin, pot, biglab
        complex(kind=8) :: zin, z1in, z2in, z, zplus1, zmin1
        complex(kind=8) :: pcor, comega
    
        ! lstype = 1 means line-sink
        lstype = 1
    
        ! Radius of convergence
        Rconv = 7.0
        
        !if (ilap==1) then
        !    istart = 1
        !else
        !    istart = 0
        !end if
    
        zin = cmplx(x,y,kind=8); z1in = z1; z2in = z2
        Lin = abs(z2in-z1in)
        z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in)
        zplus1 = z + rONE; zmin1 = z - rONE
    
        ! If at cornerpoint, move slightly
        if  (abs(zplus1) < tiny*rTWO/Lin) zplus1 = zplus1 + tiny
        if  (abs(zmin1) < tiny*rTWO/Lin) zmin1 = zmin1 + tiny
    
        ! Laplace linesink

        if (ilap == 1) then
            power = order + 1 ;
            pcor = cmplx(0.d0,0.d0,kind=8) ;
            do n=1,(power+1)/2
                pcor = pcor + z**(power-2*n+1) / (2*n-1) ;
            end do
            pcor = rTWO * pcor ;
        
            comega = z**power * log( (zmin1) / (zplus1) ) + pcor - log(zmin1) + (-rONE)**power*log(zplus1) ;
            comega = -comega*Lin/(rFOUR*pi*power)
        
            rv(1) = real(comega)
        end if
    
        ! N-1 leakage factors
        do i=ilap+1, naq
    
            ! Check whether entire linesink is outside radius of convergence
            ! Outside if |z-zc|>L/2+7lab, and thus |Z|>1+7lab*2/L, or |zeta|>1/biglab+7 (zeta is called z here)
            biglab = rTWO * lambda(i) / Lin
            z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in) / biglab
    
            if ( abs(z) < (Rconv + rONE/biglab) ) then
                pot = IntegralF(zin,z1in,z2in,Lin,lambda(i),order,Rconv,lstype)
                rv(i) = -Lin/rTWO * pot
            else
                rv(i) = rZERO
            end if
        end do
    
    end function potbeslsho
    
    function potbeslsv(x, y, z1, z2, lab, order, ilap, naq) result(pot)
        implicit none
        real(kind=8), intent(in) :: x, y
        complex(kind=8), intent(in) :: z1, z2
        integer, intent(in) :: naq
        real(kind=8), dimension(naq), intent(in) :: lab
        integer, intent(in) :: order, ilap
        real(kind=8), dimension(order+1, naq) :: pot
        ! locals
        integer :: n
        ! Check if endpoints need to be adjusted using the largest lambda (the first one)
        do n = 0, order
            !print *, n*naq + 1, (n+1)*naq
            pot(n+1, 1:naq) = potbeslsho(x,y,z1,z2,lab,n,ilap,naq)
        end do
    end function potbeslsv
    
    function disbeslsho(x,y,z1,z2,lambda,order,ilap,naq) result(rv)
    
        ! Input:
        !   x,y: Point where discharge is computed
        !   z1: Complex begin point of line-sink
        !   z2: Complex end point of line-sink
        !   lambdain(Naquifers): Array with zero in first spot and lambda's in remaining spots
        !   order: Order of the line-sink
        !   Naquifers: Number of aquifers
        !   rvx(Naquifers),rvy(Naquifers): Arrays to store return values
        ! Output:
        !   rvx(Naquifers),rvy(Naquifers): Values of Qx and Qy with Laplace value in first spot
        !                   and mod.Helmholtz potentials in remaining spots
    
        implicit none
        ! Input / Output
        real(kind=8), intent(in) :: x, y
        complex(kind=8), intent(in) :: z1, z2
        integer, intent(in) :: naq, order, ilap
        real(kind=8), dimension(naq), intent(in) :: lambda
        real(kind=8), dimension(2,naq) :: rv        
    
        ! Locals
        integer :: n, i, lstype
        real(kind=8) :: Rconv, Lin, biglab
        complex(kind=8) :: zin, z1in, z2in, z, zplus1, zmin1
        complex(kind=8) :: pcor, wdis, cdum
    
        ! Radius of convergence
        Rconv = 7.0
    
        ! lstype = 1 means line-sink
        lstype = 1
    
        zin = cmplx(x,y,8); z1in = z1; z2in = z2
        Lin = abs(z2in-z1in)
        z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in)
        zplus1 = z + rONE; zmin1 = z - rONE
    
        ! If at cornerpoint, move slightly
        if  (abs(zplus1) < tiny*rTWO/Lin) zplus1 = zplus1 + tiny
        if  (abs(zmin1) < tiny*rTWO/Lin) zmin1 = zmin1 + tiny
    
        ! Laplace linesink
        if (ilap == 1) then
            pcor = cZERO ;
            do n=1,(order+1)/2
                pcor = pcor + float(order-2*n+2) * z**(order+1-2*n) / float(2*n-1) ;
            end do
            pcor = rTWO * pcor ;
        
            cdum = 1.d0 / (order+1)  ! Without this intermediate statement it didn't seem to work
            wdis = float(order+1) * z**order * log( (zmin1) / (zplus1) ) + pcor
            
            wdis = wdis + (z**(order+1) - rONE) / zmin1 - (z**(order+1) - (-rONE)**(order+1)) / zplus1
        
            wdis = wdis*Lin/rTWO/(z2in-z1in)/pi * cdum;
        
            rv(1,1) = real(wdis)
            rv(2,1) = -aimag(wdis)
        end if
    
        do i=ilap+1, naq
            wdis = cZERO
    
            ! Check whether entire linesink is outside radius of convergence
            ! Outside if |z-zc|>L/2+7lab, and thus |Z|>1+7lab*2/L, or |zeta|>1/biglab+7 (zeta is called z here)
            biglab = rTWO * lambda(i) / Lin
            z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in) / biglab
    
            if ( abs(z) < (Rconv + rONE/biglab) ) then
                wdis = IntegralG(zin,z1in,z2in,Lin,lambda(i),order,Rconv,lstype)
                wdis = rTWO * Lin / (z2in-z1in) / biglab * wdis
    
                rv(1,i) = real(wdis)
                rv(2,i) = -aimag(wdis)
            else
                rv(1,i) = rZERO; rv(2,i) = rZERO
            end if
    
        end do
    
    end function disbeslsho
    
    function disbeslsv(x, y, z1, z2, lab, order, ilap, naq) result(qxqy)
        implicit none
        real(kind=8), intent(in) :: x, y
        complex(kind=8), intent(in) :: z1, z2
        integer, intent(in) :: naq
        real(kind=8), dimension(naq), intent(in) :: lab
        integer, intent(in) :: order, ilap
        real(kind=8), dimension(2*(order+1), naq) :: qxqy
        ! locals
        integer :: n
        real(kind=8), dimension(2,naq) :: rv
        ! Check if endpoints need to be adjusted using the largest lambda (the first one)
        do n = 0, order
            rv = disbeslsho(x,y,z1,z2,lab,n,ilap,naq)
            qxqy(n+1, 1:naq) = rv(1,1:naq)
            qxqy(n+1+order+1, 1:naq) = rv(2,1:naq)
        end do
    end function disbeslsv
    
    function potbesldho(x, y, z1, z2, lambda, order, ilap, naq) result (rv)
    
        ! Input:
        !   x,y: Point where potential is computed
        !   z1: Complex begin point of line-doublet
        !   z2: Complex end point of line-doublet
        !   lambda(naq): lambda's (zero for first lambda if Laplace)
        !   order: Order of the line-doublet
        !   ilap: equals 1 when first value is Laplace line-doublet and first lambda equals zero
        !   naq: Number of aquifers
        !   rv(naq): Array to store return value (must be pre-allocated)
        ! Output:
        !   rv(naq): Potentials. Fist spot is Laplace value if ilap=1
    
        implicit none
    
        ! Input / Output
        real*8, intent(in) :: x, y
        complex(kind=8), intent(in) :: z1, z2
        integer, intent(in) :: naq, order, ilap
        real(kind=8), dimension(naq), intent(in) :: lambda
        real(kind=8), dimension(naq) :: rv
    
        ! Locals
        integer, parameter :: Nterms = 8
        integer :: n, i, lstype
        integer :: NLS, m1, m2
        real(kind=8) :: Rconv, Lin, pot, biglab, del0, ra
        complex(kind=8) :: zin, z1in, z2in, z, zplus1, zmin1, z1new, z2new
        complex(kind=8) :: comega, qm
    
        ! Radius of convergence
        Rconv = 7.0
        
        ! lstype=2 means line-doublet
        lstype = 2
    
        zin = cmplx(x, y, 8); z1in = z1; z2in = z2
        Lin = abs(z2in - z1in)
        z = (rTWO * zin - (z1in + z2in) ) / (z2in - z1in)
        zplus1 = z + rONE; zmin1 = z - rONE
    
        ! If at cornerpoint, move slightly
        if  (abs(zplus1) < tiny * rTWO / Lin) zplus1 = zplus1 + tiny
        if  (abs(zmin1) < tiny * rTWO / Lin) zmin1 = zmin1 + tiny
    
        ! Laplace line-doublet
        if (ilap == 1) then
            comega = z**order * log(zmin1 / zplus1)
            qm = cZERO
            do n=1,(order+1)/2
                qm = qm + z**(order-rTWO*float(n)+rONE) / (rTWO*float(n) - rONE)
            end do
            comega = rONE/(rTWO*pi*ci) * ( comega + rTWO * qm )
            rv(1) = real(comega)
        end if
            
        ! N-1 leakage factors
        do i=ilap+1, naq
            pot = rZERO
    
            ! Check whether entire linedoublet is outside radius of convergence
            ! Outside if |z-zc|>L/2+7lab, and thus |Z|>1+7lab*2/L, or |zeta|>1/biglab+7 (zeta is called z here)
            biglab = rTWO * lambda(i) / Lin
            z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in) / biglab
    
            if ( abs(z) < (Rconv + rONE/biglab) ) then
                call findm1m2(zin,z1in,z2in,Lin,lambda(i),Rconv,m1,m2,NLS)
                comega = cZERO
                if (m1 > 0) then   ! Otherwise outside radius of convergence
                    z1new = z1in + float(m1-1)/float(NLS) * (z2in-z1in)
                    z2new = z1in + float(m2)/float(NLS)   * (z2in-z1in)
                    del0 = float(1-m1-m2+NLS)/float(1-m1+m2)
                    ra = float(NLS) / float(1+m2-m1)
                    comega = IntegralLapLineDipole(zin,z1new,z2new,del0,ra,order)
                end if
                pot = IntegralF(zin,z1in,z2in,Lin,lambda(i),order,Rconv,lstype)
                rv(i) = real(comega/ci) + aimag(z) / biglab * pot  ! Note that z is really zeta in analysis
            else
                rv(i) = rZERO
            end if
    
        end do
        
    end function potbesldho
    
    function potbesldv(x, y, z1, z2, lab, order, ilap, naq) result(pot)
        implicit none
        real(kind=8), intent(in) :: x, y
        complex(kind=8), intent(in) :: z1, z2
        integer, intent(in) :: naq
        real(kind=8), dimension(naq), intent(in) :: lab
        integer, intent(in) :: order, ilap
        real(kind=8), dimension(order+1, naq) :: pot
        ! locals
        integer :: n
        ! Check if endpoints need to be adjusted using the largest lambda (the first one)
        do n = 0, order
            !print *, n*naq + 1, (n+1)*naq
            pot(n+1, 1:naq) = potbesldho(x,y,z1,z2,lab,n,ilap,naq)
        end do
    end function potbesldv
    
    function disbesldho(x,y,z1,z2,lambda,order,ilap,naq) result(rv)
    
        ! Input:
        !   x,y: Point where discharge is computed
        !   z1: Complex begin point of line-sink
        !   z2: Complex end point of line-sink
        !   lambdain(Naquifers): Array with zero in first spot and lambda's in remaining spots
        !   order: Order of the line-sink
        !   naq: Number of aquifers
        ! Output:
        !   rv(2, Naquifers),rvy(Naquifers): Values of Qx and Qy with Laplace value in first spot
        !                   and mod.Helmholtz potentials in remaining spots
    
        implicit none    
        ! Input / Output
        real(kind=8), intent(in) :: x, y
        complex(kind=8), intent(in) :: z1, z2
        integer, intent(in) :: naq, order, ilap
        real(kind=8), dimension(naq), intent(in) :: lambda
        real(kind=8), dimension(2,naq) :: rv
        
        ! Locals
        integer :: n, i, m1, m2, NLS, lstype
        real(kind=8) :: Rconv, Lin, biglab, del0, ra, pot
        complex(kind=8) :: zin, z1in, z2in, z, zplus1, zmin1, z1new, z2new
        complex(kind=8) :: qm, wdis, wdis1, wdis2, wdis3
    
        ! Locals
    
        !integer :: n, m, i, m1, m2, NLS
        !real*8 :: Rconv, Lin, biglab, del0, ra, pot
        !complex*16 :: zin, z1in, z2in, z, z1, z2, zplus1, zmin1
        !complex*16 :: qm, wdis, wdis1, wdis2, wdis3
        
        ! Radius of convergence
        Rconv = 7.0
        
        ! lstype=2 means line-doublet
        lstype = 2
    
        zin = cmplx(x,y,8); z1in = z1; z2in = z2
        Lin = abs(z2in-z1in)
        z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in)
        zplus1 = z + rONE; zmin1 = z - rONE
    
        ! If at cornerpoint, move slightly
        if  (abs(zplus1) < tiny*rTWO/Lin) zplus1 = zplus1 + tiny
        if  (abs(zmin1) < tiny*rTWO/Lin) zmin1 = zmin1 + tiny
    
        ! Laplace line-doublet
        if (ilap == 1) then
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
            rv(1,1) = real(wdis)
            rv(2,1) = -aimag(wdis)
        end if
    
        ! N-1 or N leakage factors
        do i = ilap + 1, naq
            wdis = cZERO
    
            ! Check whether entire line-doublet is outside radius of convergence
            ! Outside if |z-zc|>L/2+7lab, and thus |Z|>1+7lab*2/L, or |zeta|>1/biglab+7 (zeta is called z here)
            biglab = rTWO * lambda(i) / Lin
            z = ( rTWO*zin - (z1in + z2in) ) / (z2in - z1in) / biglab
    
            if ( abs(z) < (Rconv + rONE/biglab) ) then
                call findm1m2(zin,z1in,z2in,Lin,lambda(i),Rconv,m1,m2,NLS)
                wdis1 = cZERO
                if (m1 > 0) then
                    z1new = z1in + float(m1-1)/float(NLS) * (z2in-z1in)
                    z2new = z1in + float(m2)/float(NLS)   * (z2in-z1in)
                    del0 = float(1-m1-m2+NLS)/float(1-m1+m2)
                    ra = float(NLS) / float(1+m2-m1)
                    wdis1 = IntegralLapLineDipoleDis(zin,z1new,z2new,del0,ra,order)
                    wdis1 = -rTWO  * wdis1 / (ci*(z2new-z1new))
                end if
                pot = IntegralF(zin,z1in,z2in,Lin,lambda(i),order,Rconv,lstype)
                wdis2 = IntegralG(zin,z1in,z2in,Lin,lambda(i),order,Rconv,lstype)
                wdis3 =  pot / (rTWO * ci) + wdis2 * aimag(z)
                
                wdis = wdis1 - rFOUR * wdis3 / ( biglab**2 * (z2in-z1in) )
                
            end if 
    
            rv(1,i) = real(wdis)
            rv(2,i) = -aimag(wdis)
    
        end do
    
    end function disbesldho
    
    function disbesldv(x, y, z1, z2, lab, order, ilap, naq) result(qxqy)
        implicit none
        real(kind=8), intent(in) :: x, y
        complex(kind=8), intent(in) :: z1, z2
        integer, intent(in) :: naq
        real(kind=8), dimension(naq), intent(in) :: lab
        integer, intent(in) :: order, ilap
        real(kind=8), dimension(2*(order+1), naq) :: qxqy
        ! locals
        integer :: n
        real(kind=8), dimension(2,naq) :: rv
        ! Check if endpoints need to be adjusted using the largest lambda (the first one)
        do n = 0, order
            rv = disbesldho(x,y,z1,z2,lab,n,ilap,naq)
            qxqy(n+1, 1:naq) = rv(1,1:naq)
            qxqy(n+1+order+1, 1:naq) = rv(2,1:naq)
        end do
    end function disbesldv
    
    function IntegralF(zin,z1in,z2in,Lin,lambda,order,Rconv,lstype) result (pot)
    
        implicit none
        ! Input
        integer :: order, lstype
        real(kind=8) :: Lin, lambda, Rconv
        complex(kind=8) :: zin, z1in, z2in
    
        ! Out
        real(kind=8) :: pot
    
        ! Local
        integer :: NLS,m1,m2, m, n
        real(kind=8) :: L, biglab, del1, del2
        complex(kind=8) :: z, zbar, cInt, cd1minz, cd2minz, cln1, cln2
        complex(kind=8), dimension(0:Nterms) :: czmzbarp
        complex(kind=8), dimension(0:Nterms,0:Nterms) :: cgamma
        complex(kind=8), dimension(0:2*Nterms) :: calphat, cbetat
        complex(kind=8), dimension(0:order+1) :: cc
        complex(kind=8), dimension(0:2*Nterms+order) :: calpha, cbeta
    
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
            calphat(n) = cmplx(0.d0,0.d0,kind=8); cbetat(n)=cmplx(0.d0,0.d0,kind=8)
            do m = max(0,n-Nterms), n/2
                if (lstype == 1) then
                    calphat(n) = calphat(n) + ac(n-m) * cgamma(n-m,m)
                    cbetat(n) = cbetat(n) + bc(n-m) * cgamma(n-m,m)
                else
                    calphat(n) = calphat(n) + ac1(n-m) * cgamma(n-m,m)
                    cbetat(n) = cbetat(n) + bc1(n-m) * cgamma(n-m,m)
                end if
            end do
        end do
    
        ! Compute coefficients of delta^p
        do m = 0,order
            cc(m) = rbinom(order,m) * z**(order-m) * biglab**order
        end do
        if ( order > 0 ) then
            do n=0,2*Nterms+order
                calpha(n) = cmplx(0.d0,0.d0,kind=8); cbeta(n) = cmplx(0.d0,0.d0,kind=8)
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
        cInt = cmplx(0.d0,0.d0,kind=8)
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
    
    end function IntegralF
    
    function IntegralG(zin,z1in,z2in,Lin,lambda,order,Rconv,lstype) result(wdis)
    
        implicit none
        ! Input
        integer :: order, lstype
        real(kind=8) :: Lin, lambda, Rconv
        complex(kind=8) :: zin, z1in, z2in

        ! Out
        complex(kind=8) :: wdis
    
        ! Local
        integer :: NLS,m1,m2, m, n
        real(kind=8) :: L, biglab, del1, del2, del0, ra, biglabin
        complex(kind=8) :: z, zbar, z1, z2, cd1minz, cd2minz, cln1, cln2
        complex(kind=8) :: g1, g2, g3, comega
        complex(kind=8), dimension(0:Nterms) :: czmzbarp
        complex(kind=8), dimension(0:Nterms,0:Nterms) :: cgamma
        complex(kind=8), dimension(0:2*(Nterms-1)) :: cahat, cbhat
        complex(kind=8), dimension(0:2*Nterms) :: calphat, cbetat
        complex(kind=8), dimension(0:order+1) :: cc
        complex(kind=8), dimension(0:2*Nterms+order) :: calpha, cbeta
    
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
        !comega = cZERO
        comega = IntegralLapLineDipole(zin,z1,z2,del0,ra,order)
        
        if (lstype == 1) then
            g1 = -ac(0) * biglabin * comega
            ! Integral g2
            ! Compute hat coefficients
            do n = 0, Nterms-1
                cahat(n) = float(n+1) * ac(n+1)
                cbhat(n) = ac(n+1) + float(n+1) * bc(n+1)
            end do
        else
            g1 = -ac1(0) * biglabin * comega
            ! Integral g2
            ! Compute hat coefficients
            do n = 0, Nterms-1
                cahat(n) = float(n+1) * ac1(n+1)
                cbhat(n) = ac1(n+1) + float(n+1) * bc1(n+1)
            end do
        end if
    
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
    
    end function IntegralG
    
    function IntegralLapLineDipole(zin,z1,z2,del0,ra,order) result(comega)
    
        implicit none
        ! Input
        integer :: order
        real(kind=8) :: del0,ra
        complex(kind=8) :: zin, comega
    
        ! Local
        integer :: m, n
        real(kind=8) :: L
        complex(kind=8) :: z, z1, z2, zmin1, zplus1, zterm, qm, qmtot
        complex(kind=8), dimension(0:order+1) :: cg
    
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
        
    end function IntegralLapLineDipole
    
    function IntegralLapLineDipoleDis(zin,z1,z2,del0,ra,order) result(wdis)
    
        implicit none
        ! Input
        integer :: order
        real(kind=8) :: del0,ra
        complex(kind=8) :: zin, wdis
    
        ! Local
        integer :: m, n
        real(kind=8) :: L
        complex(kind=8) :: z, z1, z2, zmin1, zplus1, zterm1, zterm2, qm, qmtot
        complex(kind=8), dimension(0:order+1) :: cg
    
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
        
    end function IntegralLapLineDipoleDis
    
    subroutine findm1m2(zin,z1in,z2in,Lin,lambda,Rconv,m1,m2,NLS)
    
        implicit none
        ! Input
        integer :: m1, m2
        real(kind=8) :: Lin, lambda, Rconv
        complex(kind=8) :: zin, z1in, z2in
    
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
        
end module besselaesnew


!program besseltest
!
!    use besselaesnew
!    real(kind=8) :: x, y
!    complex(kind=8) :: z1, z2
!    real(kind=8), dimension(3) :: lab, pot
!    real(kind=8), dimension(12) :: potv
!    integer :: naq, ilap, order
!    
!    call initialize()
!    
!    x = 2.d0
!    y = 1.d0
!    z1 = dcmplx(-2.d0, -1.d0)
!    z2 = dcmplx(1.d0, 0.d0)
!    naq = 3
!    lab(1) = 10.d0
!    lab(2) = 1.d0
!    lab(3) = 5.d0
!    order = 3
!    ilap = 0
!    
!    do n = 0, order
!        pot = potbeslsho(x,y,z1,z2,lab,n,ilap,naq)
!        print *, 'pot: ', pot
!    end do
!    
!    potv = potbeslsv(x, y, z1, z2, lab, order, ilap, naq)
!    print *, 'potv: ', potv
!
!end


program besseltest

    use besselaesnew
    real(kind=8) :: x, y
    complex(kind=8) :: z1, z2
    real(kind=8), dimension(3) :: lab
    real(kind=8), dimension(3) :: pot
    real(kind=8), dimension(2,3) :: potv
    real(kind=8), dimension(2,3) :: qxqy
    real(kind=8), dimension(4,3) :: qxqyv
    integer :: naq, ilap, order
    
    call initialize()
    
    x = 2.d0
    y = 1.d0
    z1 = dcmplx(-3.d0, -1.d0)
    z2 = dcmplx(2.d0, 2.d0)
    naq = 3
    lab(1) = 0.d0
    lab(2) = 2.d0
    lab(3) = 11.d0
    order = 1
    ilap = 1
    
    do n = 0, order
        pot = potbesldho(x,y,z1,z2,lab,n,ilap,naq)
        print *, 'n, pot: ', n, pot
    end do
    !
    potv = potbesldv(x, y, z1, z2, lab, order, ilap, naq)
    do n = 0, order
        print *, 'n: ', potv(n+1,1:naq)
    end do
    
    qxqy = disbesldho(x,y,z1,z2,lab,0,ilap,naq)
    print *,qxqy(1,1:naq)
    print *,qxqy(2,1:naq)
    qxqy = disbesldho(x,y,z1,z2,lab,1,ilap,naq)
    print *,qxqy(1,1:naq)
    print *,qxqy(2,1:naq)
    qxqyv = disbesldv(x, y, z1, z2, lab, 1, ilap, naq)
    print *,qxqyv(1,1:naq)
    print *,qxqyv(1+order+1,1:naq)
    print *,qxqyv(2,1:naq)
    print *,qxqyv(2+order+1,1:naq)

end
