program flu2d 
  use omp_lib
  use mod_2dflu
  use mod_part_interp
  implicit none
  !! other defination of variables 
  integer ::fno, i1,i2,i3,ispectra,iouter,isnap,nouter,nalias,icall,ilink,ibead
  double precision :: En,Ent,Entdiss,epsilon,sumEkbyki,compare,dummy,wtime
  integer :: k1,k2,ksqr,ireal,iimag,mshli,vort_count
  character*80 :: fname1,fname2
!  character*500 :: cnpar,formp,formvort
!  character*1000 :: fnum

  call system('mkdir vorticity spectra')

  !! %---read flu.in ----
  open(unit=1,file='para_2d.in',status='old')
  read(1,*)
  read(1,*) nn,vis,vis2,delta,maxiter,nrun, nrunpart, pmax,navg
  read(1,*)
  read(1,*) kini,kf, famp, mu, kdrag, mu2, nthreads  
  read(1,*)
  read(1,*) particl, Nparticle, nstokes, Nbead, Req_chain, Rmax_chain, Rini_chain, navgpart
  close(1)
  !! %---------------------
  !! Obtain individual link parameters from chain inputs 
  Req=(Req_chain)/dble(Nbead-1)
  Rmax=(Rmax_chain)/dble(Nbead-1)
  Rini=Rini_chain/dble(Nbead-1)
  !! %---------------------
  if(particl) then
  	do i1 = 1,nstokes
   	 write(fname1,'(g8.0)')  i1
	 call system('mkdir tau'//trim(adjustl(fname1)))
	enddo
  else
  endif

!! % --- Set the number of threads
  call omp_set_num_threads(nthreads)
  open ( unit = 111, file = 'stat.out',status='unknown')
  write ( 111, '(a,i8)' ) &
    '  The number of processors available = ', omp_get_num_procs ( )
  write ( 111, '(a,i8)' ) &
    '  The number of threads available    = ', omp_get_max_threads ( )
  close (111)
  wtime = omp_get_wtime ( )
  !! %---------------------
  !! \section{Initialisation}  
  n1 = nn !the whole boxsize
  n2 = nn 
  n1h = n1/2
  n1hf = n1/2 + 1
  ksqr_max = 2*(nn/2)*(nn/2) 
  nshell = int(1.414*nn/2.0d0) + 1 	
  nalias = int(nn/3) 
  kasqr = 2*nalias*nalias
  	write(*,*) kasqr,nalias,ksqr_max
  !! nshell = \sqrt{ 2 nn}, the diagonal of the box in fourier space. 
  !	write(*,*) nshell
  !! The length of the box is fixed to be twice \pi.  
  pi = 4.0d0*datan(1.0d0)
  length = 2.0d0*pi
  factor = 2.0d0*pi/length  
  dx = length/dble(n1);
  dy = dx
  !! allocate and calculate the global arrays --------
  allocate(den_state(0:nshell))
  allocate(time_increment(0:ksqr_max),alpha(0:ksqr_max))
  allocate(diss_rate(maxiter),tot_energy(maxiter),energy_time(maxiter))
  allocate(energy_dissrate(maxiter),enstrphy_dissrate(maxiter))
  den_state = 0
  time_increment = 0.0d0
  call global2d
  !! %---and also energy and dissipation arrays ------
  allocate(E_Omega(0:nshell,3))
  E_Omega = 0.0d0
  !! %-stream-function and vorticity array are allocated  ----
  allocate(psi(n1+2,n2),omega(n1+2,n2),fomega(n1+2,n2))
  psi = 0.0d0
  omega = 0.0d0
  fomega = 0.0d0
  !! %-allocate Jacobean etc 
  !	allocate(jac(n1+2,n2))
  allocate(jac_old(n1+2,n2))
  !	jac = 0.0d0
  jac_old = 0.0d0
  !! -- allocate velocity in fourier space 
  !! We pad the velocity arrays so that in real space the velocity 
  !! can be easily made periodic
  allocate(ukx(-1:n1+2,-1:n2+2),uky(-1:n1+2,-1:n2+2)) 
  allocate(uktmp(n1+2,n2)) 
  !! --- Allocate arrays for particle tracking
  allocate(Xc(Nparticle,nstokes),Yc(Nparticle,nstokes),Xc0(Nparticle,nstokes),Yc0(Nparticle,nstokes))
  allocate(Xb(Nparticle,Nbead,nstokes),Yb(Nparticle,Nbead,nstokes))
  allocate(Rx(Nparticle,Nbead-1,nstokes),Ry(Nparticle,Nbead-1,nstokes),Rmag(Nparticle,Nbead-1,nstokes))
  allocate(Rx0(Nparticle,Nbead-1,nstokes),Ry0(Nparticle,Nbead-1,nstokes))
  allocate(vRx(Nparticle,Nbead-1,nstokes),vRy(Nparticle,Nbead-1,nstokes))
  allocate(uxb(Nparticle,Nbead,nstokes),uyb(Nparticle,Nbead,nstokes))
  allocate(vxc(Nparticle,nstokes),vyc(Nparticle,nstokes))
  allocate(lamc(Nparticle,nstokes),lamb(Nparticle,Nbead,nstokes))
  allocate(Lchain(Nparticle,nstokes))
!  allocate(dxux(Nparticle,nstokes),dxuy(Nparticle,nstokes),dyux(Nparticle,nstokes))
!  allocate(dxux1(Nparticle,nstokes),dxuy1(Nparticle,nstokes),dyux1(Nparticle,nstokes))
!  allocate(dxux2(Nparticle,nstokes),dxuy2(Nparticle,nstokes),dyux2(Nparticle,nstokes))
!  allocate(GXp(Nparticle,nstokes),GYp(Nparticle,nstokes))
!  allocate(GXp0(Nparticle,nstokes),GYp0(Nparticle,nstokes))
!  allocate(Gvxint0(Nparticle,nstokes),Gvyint0(Nparticle,nstokes))
!  allocate(Guxint(Nparticle,nstokes),Guyint(Nparticle,nstokes))
!  allocate(Gvxint(Nparticle,nstokes),Gvyint(Nparticle,nstokes))
!  allocate(Glamint(Nparticle,nstokes))
!  allocate(Gdxux(Nparticle,nstokes),Gdxuy(Nparticle,nstokes),Gdyux(Nparticle,nstokes))
  allocate(tau_chain(nstokes),tau(nstokes))

!$omp parallel workshare
  uxb = 0.0d0; uyb = 0.0d0;
  vxc = 0.0d0; vyc = 0.0d0; lamc = 0.0d0; lamb = 0.0d0; vRx = 0.0d0; vRy = 0.0d0;
  Xc = 0.0d0;Yc = 0.0d0; Xb = 0.0d0;Yb = 0.0d0; Rx = 0.0d0;Ry = 0.0d0;
  Xc0 = 0.0d0;Yc0 = 0.0d0;Rx0 = 0.0d0;Ry0 = 0.0d0;Lchain=0.0d0
!  dxux = 0.0d0;dyux = 0.0d0;dxuy = 0.0d0;
!  Guxint = 0.0d0; Guyint = 0.0d0
!  Gvxint = 0.0d0; Gvyint = 0.0d0
!$omp end parallel workshare
!!
  write(cnpar,'(g8.0)') Nparticle
  formp = '('//trim(adjustl(cnpar))//'es16.3E2)'
  write(cnpar,'(g8.0)') n1 
  formvort = '('//trim(adjustl(cnpar))//'es16.3E2)'
  !! -------------------------------------------------
  !! %          peculiar for FFTW               *
  !! % We set up the plan for FFTW here. At present we are not using 
  !! WISDOM*
  dim(1) = n1 
  dim(2) = n2  
  scale = 1.0d0/dfloat(n1*n2)
  !! %-------create plan for forward transform -------------------------
  call rfftwnd_f77_create_plan(pfor,2,dim,FFTW_REAL_TO_COMPLEX, & 
       FFTW_MEASURE + FFTW_IN_PLACE + FFTW_THREADSAFE)
  !! %------create plan for inverse transform -------------------------
  call rfftwnd_f77_create_plan(pinv,2,dim,FFTW_COMPLEX_TO_REAL, & 
       FFTW_MEASURE + FFTW_IN_PLACE + FFTW_THREADSAFE)
  !! %-------------------******************------------------------------
  !!
  call gen_force
  !! % ---------------------------
  !! %------------------------------------------
  !! we set an initial configuration.  
  if(nrun.eq.1) then
     call iniconf2d
     !  ------------- Initialize the particle positions
  
     if (particl)  call init_particles
     ! ---------------
     !! obtain stream function and velocity from vorticity
     call getpsivel
     call energy2d
     energy_time(1) = sum(E_Omega(:,1))	
     !	write(*,*) energy_time(1)
     energy_time = 0.0d0
     open(unit=110,file='initial_spectrum.out',status='unknown')
     do ispectra = 1,nshell
        write(110,*) dfloat(ispectra),E_Omega(ispectra,1)
     enddo
     close(110)
     call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,omega,0)
     do i1=n1+1,n1+2
        omega(i1,:) = 0.0d0
     enddo
     omega = omega*scale
     call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,omega,0)
     call getpsivel
     E_Omega = 0.0d0
     call energy2d
     open(unit=110,file='check_initial_spectrum.out',status='unknown')
     do ispectra = 1,nshell
        write(110,*) dfloat(ispectra),E_Omega(ispectra,1)
     enddo
     close(110)

     energy_time(1) = sum(E_Omega(:,1))
     !        write(*,*) energy_time(1)
     !! Initial spectra
  else
     !! Otherwise (i.e. if this is not the first run ) we read the
     !! configuration of vorticity
     open(unit=11,file='vortex.in',form='unformatted',status='old')
!!     read(11) ((omega(i1,i2),i1=1,n1+2),i2=1,n2)
     read(11) ((omega(i1,i2),i1=1,n1),i2=1,n2)
     close(11)
     do i1=n1+1,n1+2
        omega(i1,:) = 0.0d0
     enddo
     call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,omega,0)
     call getpsivel
     call energy2d
     energy_time(1) = sum(E_Omega(:,1))
     !        write(*,*) energy_time(1)
     energy_time = 0.0d0
     !! Initial spectra
     open(unit=110,file='initial_spectra.out',status='unknown')
     do ispectra = 1,nshell
        write(110,*) dfloat(ispectra),E_Omega(ispectra,1)
     enddo
     close(110)
     !! Initial particle position
     if (particl)  call init_particles
  endif


  !! \section{Time Marching}
  open(unit=15,file='energies.out',status='unknown')
counter = 1
  do icall = 1,maxiter/navg
     do iinner = 1,navg
        call rnkt
        call getpsivel
        E_Omega = 0.0d0
        call energy2d
        write(15,*) counter*delta,0.5d0*sum(E_Omega(:,1)),-vis*sum(E_Omega(:,2))
 !       write(*,*) counter*delta,0.5d0*sum(E_Omega(:,1)),-vis*sum(E_Omega(:,2))
    if (mod(counter,navgpart)==0) then
    if (particl) then
     do istokes = 1,nstokes
     fno = 1000*istokes + 30
     write(fno,formp) (Xc(i1,istokes),i1=1,Nparticle)
     fno = fno + 1
     write(fno,formp) (Yc(i1,istokes),i1=1,Nparticle)
     fno = fno + 1
     write(fno,formp) (lamc(i1,istokes),i1=1,Nparticle)
     fno = fno + 1
     write(fno,formp) (Lchain(i1,istokes),i1=1,Nparticle)

	do ilink = 1, Nbead-1
    	 fno = fno + 1
   	 write(fno,formp) (Rx(i1,ilink,istokes),i1=1,Nparticle)
    	 fno = fno + 1
    	 write(fno,formp) (Ry(i1,ilink,istokes),i1=1,Nparticle)
	enddo
	do ibead = 1, Nbead
    	 fno = fno + 1
    	 write(fno,formp) (lamb(i1,ibead,istokes),i1=1,Nparticle)
	enddo

     enddo
     else
     endif
     else
     endif
     counter = counter + 1
     enddo
     call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,omega,0)
     do i1=n1+1,n1+2
        omega(i1,:) = 0.0d0
     enddo
     omega = omega*scale
     write(fnum,'(g8.0)') icall 
   do i1=n1+1,n1+2
        omega(i1,:) = 0.0d0
    enddo
   open(unit=224,file='vorticity/vortex'//trim(adjustl(fnum))//'.out',status='unknown')
   do i1 = 1,n1
   write(224,formvort)(omega(i1,i2),i2=1,n2)
   enddo
    close(224) 
     call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,omega,0)
    open(unit=225,file='spectra/spectra'//trim(adjustl(fnum))//'.out',status='unknown')
    do ispectra= 1,nshell
      write(225,'(I9,3D16.6)') ispectra,E_Omega(ispectra,1),E_Omega(ispectra,2),E_Omega(ispectra,3)
    enddo
    close(225)
  enddo
  !!ENd of time marching

  call getpsivel
  E_Omega = 0.0d0
  call energy2d
  energy_time(1) = sum(E_Omega(:,1))
  write(*,*) 'Now here',energy_time(1)

  energy_time = 0.0d0
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,omega,0)
  do i1=n1+1,n1+2
     omega(i1,:) = 0.0d0
  enddo
  omega = omega*scale
  call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,omega,0)
  call getpsivel
  E_Omega = 0.0d0
  call energy2d
  energy_time(1) = sum(E_Omega(:,1))
  write(*,*) 'And here',energy_time(1)

  !! Save vorticity file
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,omega,0)
  omega = omega*scale
  do i1=n1+1,n1+2
     omega(i1,:) = 0.0d0
  enddo
  open(unit=110,file='vortex.out',form='unformatted',status='unknown')
  write(110) ((omega(i1,i2),i1=1,n1),i2=1,n2)
  close(110)
  !! Save the streamfunction
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,psi,0)
  omega = omega*scale
  do i1=n1+1,n1+2
     omega(i1,:) = 0.0d0
  enddo
  open(unit=110,file='stream.out',form='unformatted',status='unknown')
  write(110) ((psi(i1,i2),i1=1,n1),i2=1,n2)
  close(110)

!! Saving state of particles
  if(particl) then
    write(cnpar,'(g8.0)') Nparticle
    formp = '('//trim(adjustl(cnpar))//'es16.3E2)'

   do istokes = 1,nstokes
    write(fname1,'(g8.0)')  istokes
    fno = 1000*istokes + 70

    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/partstate.out',status='unknown')
    	write(fno,formp) (Xc(i1,istokes),i1=1,Nparticle)
	write(fno,formp) (Yc(i1,istokes),i1=1,Nparticle)

	do ilink = 1 ,Nbead-1
    		write(fno,formp) (Rx(i1,ilink,istokes),i1=1,Nparticle)
		write(fno,formp) (Ry(i1,ilink,istokes),i1=1,Nparticle)
	enddo

    close(fno)
   enddo
 else
 endif

  call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,omega,0)
  call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,psi,0)
  call getpsivel
  E_Omega = 0.0d0
  call energy2d
  energy_time(1) = sum(E_Omega(:,1))
  !        write(*,*) 'First',energy_time(1)
  !! Finally we deallocate the arrays allocated,  
  deallocate(den_state)
  deallocate(time_increment)
  deallocate(psi,omega,ukx,uky)
  deallocate(E_Omega)
  !! Destroy the plans created (for FFTW only ) 
  call rfftwnd_f77_destroy_plan(pfor)
  call rfftwnd_f77_destroy_plan(pinv)
  write(*,*) 'End'

  wtime = omp_get_wtime ( ) - wtime
  open ( unit = 111, file = 'stat.out',Access = 'append',Status='old')
  write ( 111, '(a,g14.6)' ) '  Wall clock time = ', wtime
  close (111)
  !! %-------------*********************-------------------------------
  !! and end our program 			
end program  flu2d
