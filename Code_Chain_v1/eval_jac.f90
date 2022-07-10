
!!         eval_jac.f90
!! This subroutine evaluates the Jacobean in the omega-psi
!! formulation of the 2-d navier-stokes equation. 
!! This is a serial code.This IS the part that contains
!! call to fourier transform subroutines. This SHOULD be changed
!! from one serial machine to another.
!! -----------------------------------------------------------------
subroutine eval_jac(flag) 
  use mod_2dflu
  use mod_part_interp
  use omp_lib   
!  use omp_lib
  implicit none
  double precision :: re_omega,im_omega,jplus,jminus,enx,eny,enkx,enky,ekx,eky
  double precision :: sumUx,sumUy, sumNoisex,sumNoisey, sumR
  integer ::i1,i2,i3,k1,k2,k3,ireal,iimag,ksqr
  integer:: id,tid,flag,ilink,ibead
  double precision,dimension(:,:),allocatable :: f
  double precision,allocatable,dimension(:) :: ran, fspring

!! --------------------------------------------------------------
  !! We recognise that the two derivative of psi appearing in the
  !! Jacobean are simply the velocities. Hence we first evaluate
  !! the velocities
 allocate(f(2,Nparticle))
 allocate(ran(2),fspring(Nbead-1)) 
 fspring = 0.0d0;sumR = 0.0d0

call getpsivel
  !! ------------*************************************---------------
  !! ------------using FFTW ----------
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,omega,0)
  !
  !$omp parallel workshare
  uktmp(1:n1+2,1:n2) = ukx(1:n1+2,1:n2)
  !$omp end parallel workshare
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,uktmp,0)
  !$omp parallel workshare
  ukx(1:n1,1:n2) = uktmp(1:n1,1:n2)
  !$omp end parallel workshare
  !
  !$omp parallel workshare
  uktmp(1:n1+2,1:n2) = uky(1:n1+2,1:n2)
  !$omp end parallel workshare
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,uktmp,0)
  !$omp parallel workshare
  uky(1:n1,1:n2) = uktmp(1:n1,1:n2)
  !$omp end parallel workshare
  !!
  !! -------normalise etc ------------------------------
  !	write(*,*) ' normalising ...'

  do i1=n1+1,n1+2
     ukx(i1,:) = 0.0d0
     uky(i1,:) = 0.0d0
     omega(i1,:) = 0.0d0
  enddo

  !$omp parallel do default (shared) private(i1,i2)
  do i2 = 1,n2
     do i1 =1,n1 
  ukx(i1,i2) = ukx(i1,i2)*scale
  uky(i1,i2) = uky(i1,i2)*scale
  omega(i1,i2) = omega(i1,i2)*scale
     enddo
  enddo
  !$omp end parallel do

  !!
  !! --Pad the velocities

  call apply_pbc

  if (particl) then
  !!Calculate bead locations
	call beadpositions
  !!Get fluid velocity at current bead locations
        call linear_interp

  	if(flag==1) then
	  do i1=1,Nparticle

	     do ibead = 1, Nbead
	     !!Gaussian noise from uniform random numbers
     		call random_number(ran)
     		noisex(i1,ibead)=dsqrt(-2.0d0*dlog(ran(1)))*dcos(2.0d0*pi*ran(2))
     		noisey(i1,ibead)=dsqrt(-2.0d0*dlog(ran(1)))*dsin(2.0d0*pi*ran(2))
	     enddo

	 enddo
	else
	endif

  do istokes = 1,nstokes
!$omp parallel do default (none) &
!$omp private(i1,ilink,sumR,fspring) &
!$omp shared(istokes,flag,noisex,noisey,Xc,Yc,Xc0,Yc0,Rx,Ry,Rx0,Ry0,vxc,vyc,vRx,vRy) &
!$omp shared(uxb,uyb,tau,Rmag,Rmax,Req,delta,length,Nparticle,Nbead,Lchain)
  do i1=1,Nparticle
 
  if (flag==1)then
  !! case of flag = 2: initial step with delta/2 advance

   !! Storing previous CM locations and separation vectors
     Xc0(i1,istokes) = Xc(i1,istokes) 
     Yc0(i1,istokes) = Yc(i1,istokes)

   do ilink = 1, Nbead-1
     Rx0(i1,ilink,istokes) = Rx(i1,ilink,istokes)
     Ry0(i1,ilink,istokes) = Ry(i1,ilink,istokes)
   enddo

  !! Center of mass evolution
     vxc(i1,istokes) = sum(uxb(i1,:,istokes))/Nbead
     vyc(i1,istokes) = sum(uyb(i1,:,istokes))/Nbead

     Xc(i1,istokes) = Xc0(i1,istokes) + (delta/2.0d0)*vxc(i1,istokes) + &
 			dsqrt((delta/2.0d0)*(Req*Req)/(2.0d0*tau(istokes)))*sum(noisex(i1,:))/Nbead
     Yc(i1,istokes) = Yc0(i1,istokes) + (delta/2.0d0)*vyc(i1,istokes) + &
			dsqrt((delta/2.0d0)*(Req*Req)/(2.0d0*tau(istokes)))*sum(noisey(i1,:))/Nbead
 !    Xp(:,istokes) = Xp(:,istokes) - int(Xp(:,istokes)/length)*length + 0.5d0*(1.0d0-Xp(:,istokes)/abs(Xp(:,istokes)))*length 
     Xc(i1,istokes) = modulo(Xc(i1,istokes),length)
     Yc(i1,istokes) = modulo(Yc(i1,istokes),length)

 !! Separation vectors (links) evolution

    !! Cal the spring forces
    do ilink = 1, Nbead-1
       	fspring(ilink) = 1.0d0/(1.0d0-(Rmag(i1,ilink,istokes)/Rmax)**2.0d0)
    enddo

    !! Cal the total deterministic force
    do ilink = 2, Nbead-2
     vRx(i1,ilink,istokes) = (uxb(i1,ilink+1,istokes) - uxb(i1,ilink,istokes)) &
		           - 1.0d0/(4.0d0*tau(istokes))*(2.0d0*fspring(ilink)*Rx(i1,ilink,istokes) &
				-fspring(ilink+1)*Rx(i1,ilink+1,istokes)-fspring(ilink-1)*Rx(i1,ilink-1,istokes))
     vRy(i1,ilink,istokes) = (uyb(i1,ilink+1,istokes) - uyb(i1,ilink,istokes)) &
		           - 1.0d0/(4.0d0*tau(istokes))*(2.0d0*fspring(ilink)*Ry(i1,ilink,istokes) &
				-fspring(ilink+1)*Ry(i1,ilink+1,istokes)-fspring(ilink-1)*Ry(i1,ilink-1,istokes))
    enddo
    ilink = 1
     vRx(i1,ilink,istokes) = (uxb(i1,ilink+1,istokes) - uxb(i1,ilink,istokes)) &
		           - 1.0d0/(4.0d0*tau(istokes))*(2.0d0*fspring(ilink)*Rx(i1,ilink,istokes) &
				-fspring(ilink+1)*Rx(i1,ilink+1,istokes))
     vRy(i1,ilink,istokes) = (uyb(i1,ilink+1,istokes) - uyb(i1,ilink,istokes)) &
		           - 1.0d0/(4.0d0*tau(istokes))*(2.0d0*fspring(ilink)*Ry(i1,ilink,istokes) &
				-fspring(ilink+1)*Ry(i1,ilink+1,istokes))

    ilink = Nbead-1
     vRx(i1,ilink,istokes) = (uxb(i1,ilink+1,istokes) - uxb(i1,ilink,istokes)) &
		           - 1.0d0/(4.0d0*tau(istokes))*(2.0d0*fspring(ilink)*Rx(i1,ilink,istokes) &
				-fspring(ilink-1)*Rx(i1,ilink-1,istokes))
     vRy(i1,ilink,istokes) = (uyb(i1,ilink+1,istokes) - uyb(i1,ilink,istokes)) &
		           - 1.0d0/(4.0d0*tau(istokes))*(2.0d0*fspring(ilink)*Ry(i1,ilink,istokes) &
				-fspring(ilink-1)*Ry(i1,ilink-1,istokes))
 
   !! Propagating the separation vectors

   do ilink = 1,Nbead-1 	
     Rx(i1,ilink,istokes) = Rx0(i1,ilink,istokes) + (delta/2.0d0)*vRx(i1,ilink,istokes) + &
			dsqrt((delta/2.0d0)*(Req*Req)/(2.0d0*tau(istokes)))*(noisex(i1,ilink+1) - noisex(i1,ilink))
     Ry(i1,ilink,istokes) = Ry0(i1,ilink,istokes) + (delta/2.0d0)*vRy(i1,ilink,istokes) + &
			dsqrt((delta/2.0d0)*(Req*Req)/(2.0d0*tau(istokes)))*(noisey(i1,ilink+1) - noisey(i1,ilink))
     Rmag(i1,ilink,istokes) = dsqrt(Rx(i1,ilink,istokes)**2.0d0 + Ry(i1,ilink,istokes)**2.0d0)

     if(Rmag(i1,ilink,istokes) .ge. (1.0d0-dsqrt(delta/tau(istokes)))*Rmax) then
	Rx(i1,ilink,istokes)=Rx0(i1,ilink,istokes); Ry(i1,ilink,istokes)=Ry0(i1,ilink,istokes);
        Rmag(i1,ilink,istokes) = dsqrt(Rx(i1,ilink,istokes)**2.0d0 + Ry(i1,ilink,istokes)**2.0d0)
     else
     endif

   enddo

  else
  !! case of flag = 2: second final step with delta advance

  !! Center of mass evolution
     vxc(i1,istokes) = sum(uxb(i1,:,istokes))/Nbead
     vyc(i1,istokes) = sum(uyb(i1,:,istokes))/Nbead

     Xc(i1,istokes) = Xc0(i1,istokes) + (delta)*vxc(i1,istokes) + &
 			dsqrt((delta)*(Req*Req)/(2.0d0*tau(istokes)))*sum(noisex(i1,:))/Nbead
     Yc(i1,istokes) = Yc0(i1,istokes) + (delta)*vyc(i1,istokes) + &
			dsqrt((delta)*(Req*Req)/(2.0d0*tau(istokes)))*sum(noisey(i1,:))/Nbead
 !    Xp(:,istokes) = Xp(:,istokes) - int(Xp(:,istokes)/length)*length + 0.5d0*(1.0d0-Xp(:,istokes)/abs(Xp(:,istokes)))*length 
     Xc(i1,istokes) = modulo(Xc(i1,istokes),length)
     Yc(i1,istokes) = modulo(Yc(i1,istokes),length)

 !! Separation vectors (links) evolution

    !! Cal the spring forces
    do ilink = 1, Nbead-1
       	fspring(ilink) = 1.0d0/(1.0d0-(Rmag(i1,ilink,istokes)/Rmax)**2.0d0)
    enddo

    !! Cal the total deterministic force
    do ilink = 2, Nbead-2
     vRx(i1,ilink,istokes) = (uxb(i1,ilink+1,istokes) - uxb(i1,ilink,istokes)) &
		           - 1.0d0/(4.0d0*tau(istokes))*(2.0d0*fspring(ilink)*Rx(i1,ilink,istokes) &
				-fspring(ilink+1)*Rx(i1,ilink+1,istokes)-fspring(ilink-1)*Rx(i1,ilink-1,istokes))
     vRy(i1,ilink,istokes) = (uyb(i1,ilink+1,istokes) - uyb(i1,ilink,istokes)) &
		           - 1.0d0/(4.0d0*tau(istokes))*(2.0d0*fspring(ilink)*Ry(i1,ilink,istokes) &
				-fspring(ilink+1)*Ry(i1,ilink+1,istokes)-fspring(ilink-1)*Ry(i1,ilink-1,istokes))
    enddo
    ilink = 1
     vRx(i1,ilink,istokes) = (uxb(i1,ilink+1,istokes) - uxb(i1,ilink,istokes)) &
		           - 1.0d0/(4.0d0*tau(istokes))*(2.0d0*fspring(ilink)*Rx(i1,ilink,istokes) &
				-fspring(ilink+1)*Rx(i1,ilink+1,istokes))
     vRy(i1,ilink,istokes) = (uyb(i1,ilink+1,istokes) - uyb(i1,ilink,istokes)) &
		           - 1.0d0/(4.0d0*tau(istokes))*(2.0d0*fspring(ilink)*Ry(i1,ilink,istokes) &
				-fspring(ilink+1)*Ry(i1,ilink+1,istokes))

    ilink = Nbead-1
     vRx(i1,ilink,istokes) = (uxb(i1,ilink+1,istokes) - uxb(i1,ilink,istokes)) &
		           - 1.0d0/(4.0d0*tau(istokes))*(2.0d0*fspring(ilink)*Rx(i1,ilink,istokes) &
				-fspring(ilink-1)*Rx(i1,ilink-1,istokes))
     vRy(i1,ilink,istokes) = (uyb(i1,ilink+1,istokes) - uyb(i1,ilink,istokes)) &
		           - 1.0d0/(4.0d0*tau(istokes))*(2.0d0*fspring(ilink)*Ry(i1,ilink,istokes) &
				-fspring(ilink-1)*Ry(i1,ilink-1,istokes))
 
   !! Propagating the separation vectors
   sumR = 0.0d0
   do ilink = 1,Nbead-1 	
     Rx(i1,ilink,istokes) = Rx0(i1,ilink,istokes) + (delta)*vRx(i1,ilink,istokes) + &
			dsqrt((delta)*(Req*Req)/(2.0d0*tau(istokes)))*(noisex(i1,ilink+1) - noisex(i1,ilink))
     Ry(i1,ilink,istokes) = Ry0(i1,ilink,istokes) + (delta)*vRy(i1,ilink,istokes) + &
			dsqrt((delta)*(Req*Req)/(2.0d0*tau(istokes)))*(noisey(i1,ilink+1) - noisey(i1,ilink))
     Rmag(i1,ilink,istokes) = dsqrt(Rx(i1,ilink,istokes)**2.0d0 + Ry(i1,ilink,istokes)**2.0d0)

     if(Rmag(i1,ilink,istokes) .ge. (1.0d0-dsqrt(delta/tau(istokes)))*Rmax) then
	Rx(i1,ilink,istokes)=Rx0(i1,ilink,istokes); Ry(i1,ilink,istokes)=Ry0(i1,ilink,istokes);
        Rmag(i1,ilink,istokes) = dsqrt(Rx(i1,ilink,istokes)**2.0d0 + Ry(i1,ilink,istokes)**2.0d0)
     else
     endif
     sumR = sumR + Rmag(i1,ilink,istokes)
   enddo
   Lchain(i1,istokes) = sumR  
  
  endif

 enddo
!$omp end parallel do

enddo

 if (mod(counter,navgpart)==0 .and. (flag .ne. 1)) then
        call beadpositions
        call eval_lam_interp
  endif 

else
endif

  !$omp parallel do default (shared) private(i1,i2)
  do i2=1,n2
     do i1=1,n1
        ukx(i1,i2) = ukx(i1,i2)*omega(i1,i2)
        uky(i1,i2) =  uky(i1,i2)*omega(i1,i2)
     enddo
  enddo
  !$omp end parallel do


  !! ------------using FFTW ----------
  !!       Forward Transform
  !$omp parallel workshare
  uktmp = 0.0d0
  !$omp end parallel workshare
  !$omp parallel workshare
  uktmp(1:n1,1:n2) = ukx(1:n1,1:n2)
  !$omp end parallel workshare
  call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,uktmp,0)
  !$omp parallel workshare
  ukx(1:n1+2,1:n2) = uktmp(1:n1+2,1:n2)
  !$omp end parallel workshare
  !%
  !$omp parallel workshare
  uktmp(1:n1,1:n2) = uky(1:n1,1:n2)
  !$omp end parallel workshare
  call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,uktmp,0)
  !$omp parallel workshare
  uky(1:n1+2,1:n2) = uktmp(1:n1+2,1:n2)
  !$omp end parallel workshare
  !! -----------------------------------------------

 deallocate(f)
 deallocate(ran,fspring) 
end subroutine eval_jac
