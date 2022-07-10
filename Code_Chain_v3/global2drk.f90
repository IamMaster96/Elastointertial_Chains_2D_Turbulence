!        global2d.f90
!! This subroutine evaluates the global arrays used in the serial 
!! spectral code. This doesnot use any system call thus should
!! remain unchanged from one serial machine to another. 
!! changed for slaved adam-bashforth scheme. 
!! -----------------------------------------------------------------
			subroutine global2d
	use mod_2dflu
	implicit none
	integer :: j,l,k1,k2,ksqr,mshl,imode
	double precision :: rk,rk2,nu
!! defination of other variables 
!! --the density of states calculations ----
!! --the formula below does the folding in from i to j 
!! where i varies from 1 to n1 and j should be from -n1/2 + 1 
!! to + n1/2 : 
!!          j = (i -1) - n1*(i/(n1hf+1))
!! this has been checked in the program modulo.f in home/res/fortran
!! in the theory machines.  
!! -----------------------------------------
	factor = 2.0d0*pi/length
!! ---
	do j=1,n2
		k2 = (j-1) - n2*(j/(n1hf+1))	
		do l=1,n1hf
			k1=l-1
!! -------------
			ksqr = k1*k1 + k2*k2 
			rk2 = factor*factor*dfloat(ksqr)
			rk = dsqrt(dfloat(ksqr))
			mshl = nint(rk)
			nu = vis + vis2*rk2*rk2*rk2
			if(mshl > kdrag)then
			time_increment(ksqr)= dexp(-delta*(nu*rk2+mu2)/2.0d0)
			else
			time_increment(ksqr)= dexp(-delta*(nu*rk2+mu)/2.0d0)
			endif
!			write(*,*) ksqr,mshl,time_increment(ksqr)
!! -------------
			if((k1.eq.0).or.(k1.eq.n1h)) then
			den_state(mshl)=den_state(mshl) + 1
			else
			den_state(mshl)=den_state(mshl) + 2
			endif
		enddo
	enddo
! --------------------------------------------
			end subroutine global2d
