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
        integer :: power, n, i
        !integer :: istart
        real(kind=8) :: Rconv, Lin, pot, biglab
        complex(kind=8) :: zin, z1in, z2in, z, zplus1, zmin1
        complex(kind=8) :: pcor, comega
    
        ! Radius of convergence
        Rconv = 7.0
    
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
                pot = IntegralF(zin,z1in,z2in,Lin,lambda(i),order,Rconv)
                rv(i) = -Lin/rTWO * pot
            else
                rv(i) = rZERO
            end if
        end do
    
    end function potbeslsho