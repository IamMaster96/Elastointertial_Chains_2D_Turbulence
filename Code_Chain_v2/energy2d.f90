!!      energy2d.f90 
!! Calculates the spectra for a given velocity field, this is the
!! serial version. This should not change as we move from one 
!! serial machine to another
		subroutine energy2d
	use omp_lib
        use mod_2dflu
	implicit none
	double precision :: energy_norm,rk2,rk,energy,dissipation, &
                      one_by_nsqr, volume_element,enstrophy,entdiss
	integer :: i1,i2,i3,k1,k2,k3,ksqr,mshl,ireal,iimag 
!! --------------------------------------------------------------
	one_by_nsqr = 1.0d0/dfloat(n1*n2)
	energy_norm = one_by_nsqr*one_by_nsqr
!$omp parallel do  &
!$omp default(shared) &
!$omp private(i1,i2,ireal,iimag,k1,k2,ksqr,rk2,rk,mshl,energy,enstrophy,entdiss)
	do i2 = 1,n2
		k2 = (i2-1) - n2*(i2/(n1hf+1))	
		do i1 =1,n1hf
			k1 = i1 -1
			ksqr = k1*k1 + k2*k2  
			rk2 = factor*factor*dfloat(ksqr)
			rk = dsqrt(dble(ksqr))
			mshl = nint(rk)
			if(ksqr.gt.kasqr)then
			E_Omega(mshl,1) = 0.0d0
			E_Omega(mshl,2) = 0.0d0
			E_Omega(mshl,3) = 0.0d0
			else
			ireal = 2*i1 -1   !!Actually the read value is stored second by fftw
			iimag = 2*i1
		energy = energy_norm*(ukx(ireal,i2)**2 +  &
			ukx(iimag,i2)**2 + uky(ireal,i2)**2 + uky(iimag,i2)**2 )
		enstrophy = energy_norm*(omega(ireal,i2)**2 + & 
			omega(iimag,i2)**2) 
			entdiss = rk2*enstrophy
			!$omp critical
			if((k1.eq.0).or.(k1.eq.n1h)) then
			E_Omega(mshl,1) = E_Omega(mshl,1) + energy
			E_Omega(mshl,2) = E_Omega(mshl,2) + enstrophy 
			E_Omega(mshl,3) = E_Omega(mshl,3) + entdiss
			else
			E_Omega(mshl,1) = E_Omega(mshl,1) + 2.0d0*energy
			E_Omega(mshl,2) = E_Omega(mshl,2) + 2.0d0*enstrophy
			E_Omega(mshl,3) = E_Omega(mshl,3) + 2.0d0*entdiss
			endif
			!$omp end critical
			endif	
!! ----------------- 
		enddo
	enddo
!$omp end parallel do
!! --------------------------------------------------------
			end subroutine energy2d
