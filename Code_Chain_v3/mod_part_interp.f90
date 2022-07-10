!! Module to be used if particle present
module mod_part_interp
  implicit none
  save
  !
  real*8,allocatable,dimension(:,:)::Xc,Yc,Xc0,Yc0,Lchain,vxc0,vyc0
  real*8,allocatable,dimension(:,:,:)::Xb,Yb
  real*8,allocatable,dimension(:,:,:)::Rx,Ry,Rx0,Ry0,Rmag
  real*8,allocatable,dimension(:,:)::vxc,vyc
  real*8,allocatable,dimension(:,:,:)::vRx,vRy,vRx0,vRy0
  real*8,allocatable,dimension(:,:)::lamc
  real*8,allocatable,dimension(:,:,:)::lamb
  real*8,allocatable,dimension(:,:,:)::uxb,uyb,vxb,vyb
! real*8,allocatable,dimension(:,:)::dxux,dxuy,dyux,dxux1,dxuy1,dyux1,dxux2,dxuy2,dyux2
!  real*8,allocatable,dimension(:,:)::GXp,GYp,GXp0,GYp0
!  real*8,allocatable,dimension(:,:)::Gvxint,Gvyint,Gvxint0,Gvyint0
!  real*8,allocatable,dimension(:,:)::Guxint,Guyint,Glamint
!  real*8,allocatable,dimension(:,:)::Gdxux,Gdxuy,Gdyux
  real*8,allocatable,dimension(:)::tau_chain,tau
        double precision :: taup
  double precision,allocatable,dimension(:,:) :: noisex,noisey
  character*500 :: cnpar,formp,formvort
  character*1000 :: fnum
  integer::Nparticle
  integer::nstokes,istokes
  logical:: particl
  !
contains
  !! ---------------SUBROUINES ---------------
  !%--
  subroutine init_particles
    use mod_2dflu
    implicit none
    integer::i1,istokes,fno,ilink, ibead
    real*8,allocatable,dimension(:) :: ran
    double precision :: sumR
    character*80 :: fname1,fname2
    allocate(ran(Nbead-1))
    allocate(noisex(Nparticle,Nbead),noisey(Nparticle,Nbead)) 
    ran=0.0d0


    open(unit = 12,file='stokes.in',status='old')
    do i1 = 1,nstokes
    	read(12,*) tau_chain(i1)
    	tau(i1)=6.0d0*tau_chain(i1)/dble(Nbead*(Nbead+1)) 
    enddo
    close(12)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!(mohit_code_editing)
    open(unit = 12,file='taup.in',status='old')
    
    	read(12,*) taup
    
    close(12)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!(mohit_code_editing)
    call random_seed()

    !! Initialize particle position
  if(nrunpart.eq.1) then

    do istokes = 1,nstokes
        call random_number(Xc(:,istokes)); call random_number(Yc(:,istokes));
      call random_number(vxc(:,istokes)); call random_number(vyc(:,istokes));


     	do i1=1,Nparticle

    		Xc(i1,istokes) = modulo(Xc(i1,istokes)*length,length); 
		Yc(i1,istokes) = modulo(Yc(i1,istokes)*length,length);
		sumR=0.0d0

	    do ilink = 1 ,Nbead-1
        	call random_number(ran)
    		Rx(i1,ilink,istokes) = Rini*(ran(ilink))/dsqrt((ran(ilink))**2.0d0 &
					+(1.0d0-ran(ilink))**2.0d0)
    		Ry(i1,ilink,istokes) = Rini*(1.0d0-ran(ilink))/dsqrt((ran(ilink))**2.0d0 &
					+(1.0d0-ran(ilink))**2.0d0)
    		Rmag(i1,ilink,istokes) = dsqrt(Rx(i1,ilink,istokes)**2.0d0 &
						+ Ry(i1,ilink,istokes)**2.0d0)
		sumR = sumR + Rmag(i1,ilink,istokes)
	   enddo

	  Lchain(i1,istokes) = sumR

    	enddo

    enddo
  else

    write(cnpar,'(g8.0)') Nparticle
    formp = '('//trim(adjustl(cnpar))//'es16.3E2)'
    do istokes = 1,nstokes
    	write(fname1,'(g8.0)')  istokes
    	fno = 1000*istokes + 80
	open(unit=fno,file='tau'//trim(adjustl(fname1))//'/partstate.in',status='unknown')
    	read(fno,formp) (Xc(i1,istokes),i1=1,Nparticle)
	read(fno,formp) (Yc(i1,istokes),i1=1,Nparticle)
        read(fno,formp) (vxc(i1,istokes),i1=1,Nparticle)
        read(fno,formp) (vyc(i1,istokes),i1=1,Nparticle)
	do ilink = 1 ,Nbead-1
    		read(fno,formp) (Rx(i1,ilink,istokes),i1=1,Nparticle)
		read(fno,formp) (Ry(i1,ilink,istokes),i1=1,Nparticle)
	enddo
	close(fno)
     	

        do i1=1,Nparticle

	        sumR=0.0d0
		do ilink = 1 ,Nbead-1
    		    Rmag(i1,ilink,istokes) = & 
			dsqrt(Rx(i1,ilink,istokes)**2.0d0 + Ry(i1,ilink,istokes)**2.0d0)
		    sumR = sumR + Rmag(i1,ilink,istokes)
    		enddo
		Lchain(i1,istokes) = sumR

	enddo

    enddo

  endif

     if(particl) then
     do istokes = 1,nstokes
    write(fname1,'(g8.0)')  istokes
    fno = 1000*istokes + 30
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/xctrack.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/yctrack.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/lamtrack.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/Lchaintrack.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/vxctrack.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/vyctrack.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/uxint.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/uyint.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dxux.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dxuy.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dyux.out',status='unknown')
	do ilink = 1, Nbead-1
	 write(fname2,'(g8.0)')  ilink
    	 fno = fno + 1
   	 open(unit=fno,file='tau'//trim(adjustl(fname1))//'/Rx' &
			//trim(adjustl(fname2))//'track.out',status='unknown')
    	 fno = fno + 1
    	 open(unit=fno,file='tau'//trim(adjustl(fname1))//'/Ry' &
			//trim(adjustl(fname2))//'track.out',status='unknown')
	enddo
	do ibead = 1, Nbead
	 write(fname2,'(g8.0)')  ibead
    	 fno = fno + 1
    	 open(unit=fno,file='tau'//trim(adjustl(fname1))//'/lamb' &
			//trim(adjustl(fname2))//'track.out',status='unknown')
	enddo
    enddo
   else
   endif


!    Xp = pi+0.10456432211924;Yp = pi+0.861819520612071;
  end subroutine init_particles
  !%--

 !!------Position of beads from CM and separation vectors------------
 subroutine beadpositions
    use mod_2dflu
    use omp_lib

    integer :: i1,ilink,ibead,k1,j
    double precision :: sumRx,sumRy

do istokes = 1,nstokes
!$omp parallel do default (none) &
!$omp private(i1,ilink,j,k,sumRx,sumRy) &
!$omp shared(istokes,Nbead,Nparticle,Xb,Yb,Xc,Yc,Rx,Ry,length)
    do i1 = 1,Nparticle

	sumRx=0.0d0; sumRy = 0.0d0;
!	do j = 2, Nbead
!	  do k1 = 1,j-1
!	    sumRx=sumRx+Rx(i1,k1,istokes)
!	    sumRy=sumRy+Ry(i1,k1,istokes)	
!	  enddo
!	enddo

      	Xb(i1,1,istokes)=modulo(Xc(i1,istokes),length)
      	Yb(i1,1,istokes)=modulo(Yc(i1,istokes),length)

	do ilink = 1, Nbead-1
	   Xb(i1,ilink+1,istokes)=modulo(Xb(i1,ilink,istokes)+Rx(i1,ilink,istokes),length)
	   Yb(i1,ilink+1,istokes)=modulo(Yb(i1,ilink,istokes)+Ry(i1,ilink,istokes),length)	   
	enddo

    enddo
!$omp end parallel do
 enddo

 end subroutine beadpositions
 
!subroutine beadvelocities
 !   use mod_2dflu
 !   use omp_lib
!
!    integer :: i1,ilink,ibead,k1,j
!    double precision :: sumVx,sumVy

!do istokes = 1,nstokes
!!$omp parallel do default (none) &
!!$omp private(i1,ilink,j,k,sumVx,sumVy) &
!!$omp shared(istokes,Nbead,Nparticle,vxb,vyb,vxc,vyc,vRx,vRy,length)
 !   do i1 = 1,Nparticle
!
!	sumVx=0.0d0; sumVy = 0.0d0;
!	do j = 2, Nbead
!	  do k1 = 1,j-1
!	    sumVx=sumVx+vRx(i1,k1,istokes)
!	    sumVy=sumVy+vRy(i1,k1,istokes)	
!	  enddo
!	enddo
!
  !    	vxb(i1,1,istokes)=vxc(i1,istokes)
 !     	vyb(i1,1,istokes)=vyc(i1,istokes)
!
!	do ilink = 1, Nbead-1
!	   vxb(i1,ilink+1,istokes)=vxb(i1,ilink,istokes)+vRx(i1,ilink,istokes)
!	   vyb(i1,ilink+1,istokes)=vyb(i1,ilink,istokes)+vRy(i1,ilink,istokes)
!	enddo
!
  !  enddo
!!$omp end parallel do
 !enddo

! end subroutine beadvelocities
!!-----------------------------------------------------------------
 !!------Length of polymer chain: could calculate during evolution of chains------------
! subroutine total_chainlength
!    use mod_2dflu
!    use omp_lib

!    integer::ip,ilink
!    double precision :: sumR

!do istokes = 1,nstokes
!!$omp parallel do default (none) &
!!$omp private(ip,ilink,sumR) &
!!$omp shared(istokes,Nbead,Nparticle,Rmag)
!    do ip = 1,Nparticle

!	sumR=0.0d0
!	do ilink = 1, Nbead-1
!	   sumR = sumR + Rmag(ip,ilink,istokes) 	   
!	enddo
!	
!	Lchain(ip,istokes)=sumR
!    enddo
!!$omp end parallel do
! enddo

! end subroutine total_chainlength
 !!-----------------------------------------------------------------
  !! --- Linear Interpolation of velocities --
  subroutine linear_interp
    use mod_2dflu
    use omp_lib
    implicit none
    integer::i2,i1,ibead
    integer::intx,inty,intx_nxt,inty_nxt
!    integer::Gintx,Ginty
    real*8::x,y,xi,yj,xip1,yjp1
!    real*8::Gx,Gy,Gxi,Gyj,Gxip1,Gyjp1

do istokes = 1,nstokes
!$omp parallel do default (none) &
!$omp private(i1,ibead,x,y,intx,inty,intx_nxt,inty_nxt,xi,yj,xip1,yjp1) &
!$omp shared(istokes,Nbead,n1,Xb,Yb,Rx,Ry,dx,dy,ukx,uky) &
!$omp shared(uxb,uyb,length,Nparticle)
    do i1 = 1,Nparticle

      do ibead = 1, Nbead

        x = modulo(Xb(i1,ibead,istokes),length)
        y = modulo(Yb(i1,ibead,istokes),length)
	intx=floor(x/dx)+1
	intx_nxt=(1-intx/n1)*intx+1
	inty=floor(y/dy)+1
	inty_nxt=(1-inty/n1)*inty+1
!  print*, istokes,i1,intx,intx_nxt,inty,inty_nxt
!  print*, istokes,i1,Xp(i1,istokes), Yp(i1,istokes), Rx(i1,istokes), Ry(i1,istokes)
        xi = dfloat(intx-1)*dx; yj = dfloat(inty-1)*dy;
        xip1 = xi + dx; yjp1 = yj + dy;
        uxb(i1,ibead,istokes) = (xip1 - x)*(yjp1 - y)*ukx(intx,inty) + &
	    (xip1 - x)*(y - yj)*ukx(intx,inty_nxt) + &
	    (x - xi)*(yjp1 - y)*ukx(intx_nxt,inty) + &
	    (x - xi)*(y - yj)*ukx(intx_nxt,inty_nxt)
        uyb(i1,ibead,istokes) = (xip1 - x)*(yjp1 - y)*uky(intx,inty) + &
	    (xip1 - x)*(y - yj)*uky(intx,inty_nxt) + &
	    (x - xi)*(yjp1 - y)*uky(intx_nxt,inty) + &
	    (x - xi)*(y - yj)*uky(intx_nxt,inty_nxt)
        uxb(i1,ibead,istokes) = uxb(i1,ibead,istokes)/(dx*dy)
        uyb(i1,ibead,istokes) = uyb(i1,ibead,istokes)/(dx*dy)

      enddo

    enddo
!$omp end parallel do
 enddo
  end subroutine linear_interp

  !! ---------------------------------------------------
  !!
  !! Method-II (see notebook for details) to calculate 
  !! $\Lambda$ along inertial particle trajectories.
  !!
  !! ---------------------------------------------------
  subroutine eval_lam_interp
    use mod_2dflu
    implicit none
    integer::i2,i,j,ip1,jp1,ip2,jp2,im1,jm1,i1,intx,inty,ibead
    real*8::x,y,xi,yj,xip1,yjp1,xip2,yjp2,xim1,yjm1
    real*8::df00,df01,df10,df11,dxux,dxuy,dyux
!    integer::Gip1,Gjp1,Gintx,Ginty,Gi,Gjy,Gj
!    real*8::Gx,Gy,Gxi,Gyj,Gxip1,Gyjp1
!    real*8::Gdf00,Gdf01,Gdf10,Gdf11


do istokes = 1,nstokes
!$omp parallel do default (none) &
!$omp private(i1,x,y,i,j,ip1,jp1,im1,jm1,ip2,jp2,xi,yj,xip1,yjp1,xip2,yjp2,xim1,yjm1) &
!$omp private(df00,df01,df10,df11,dxux,dxuy,dyux,ibead) &
!$omp shared(istokes,Nbead,n1,Xc,Yc,Xb,Yb,dx,dy,ukx,uky) &
!$omp shared(lamc,lamb,length,Nparticle)
    do i1 = 1,Nparticle

!! Center of mass of Dumbbell
       x = modulo(Xc(i1,istokes),length); y = modulo(Yc(i1,istokes),length); 
       i=floor(x/dx)+1; j=floor(y/dy)+1;
       xi = dble(i-1)*dx; yj = dble(j-1)*dy;
       xip1 = xi + dx; yjp1 = yj + dy; 
       ip1 = floor(modulo(xip1,length)/dx)+1; jp1 = floor(modulo(yjp1,length)/dy)+1;
       xip2 = xi + 2.0d0*dx; yjp2 = yj + 2.0d0*dy; 
       ip2 = floor(modulo(xip2,length)/dx)+1; jp2 = floor(modulo(yjp2,length)/dy)+1;
       xim1 = xi - dx; yjm1 = yj - dy; 
       im1 = floor(modulo(xim1,length)/dx)+1; jm1 = floor(modulo(yjm1,length)/dy)+1;

       df00 = (ukx(ip1,j) - ukx(im1,j))/(2.0d0*dx)
       df01 = (ukx(ip1,jp1) - ukx(im1,jp1))/(2.0d0*dx)
       df10 = (ukx(ip2,j) - ukx(i,j))/(2.0d0*dx)
       df11 = (ukx(ip2, jp1) - ukx(i, jp1))/(2.0d0*dx)
       dxux = (xip1 - x)*(yjp1 - y)*df00 + &
            (xip1 - x)*(y - yj)*df01 + &
            (x - xi)*(yjp1 - y)*df10 + &
            (x - xi)*(y - yj)*df11 

       df00 = (uky(ip1,j) - uky(im1,j))/(2.0d0*dx)
       df01 = (uky(ip1,jp1) - uky(im1,jp1))/(2.0d0*dx)
       df10 = (uky(ip2,j) - uky(i,j))/(2.0d0*dx)
       df11 = (uky(ip2, jp1) - uky(i, jp1))/(2.0d0*dx)
       dxuy = (xip1 - x)*(yjp1 - y)*df00 + &
            (xip1 - x)*(y - yj)*df01 + &
            (x - xi)*(yjp1 - y)*df10 + &
            (x - xi)*(y - yj)*df11

       df00 = (ukx(i,jp1) - ukx(i,jm1))/(2.0d0*dx)
       df01 = (ukx(i,jp2) - ukx(i,j))/(2.0d0*dx)
       df10 = (ukx(ip1,jp1) - ukx(ip1,jm1))/(2.0d0*dx)
       df11 = (ukx(ip1, jp2) - ukx(ip1, j))/(2.0d0*dx)
       dyux = (xip1 - x)*(yjp1 - y)*df00 + &
            (xip1 - x)*(y - yj)*df01 + &
            (x - xi)*(yjp1 - y)*df10 + &
            (x - xi)*(y - yj)*df11

    ! Determinant of strain tensor
    dxux = dxux/(dx*dy); 
    dyux = dyux/(dx*dy);
    dxuy = dxuy/(dx*dy);
    lamc(i1,istokes) = -(dxux**2.0d0) - dxuy*dyux; 
    !! Use incompressibility, dyuy = -dxux

!! Beads
    do ibead =1, Nbead
       x = modulo(Xb(i1,ibead,istokes),length); 
       y = modulo(Yb(i1,ibead,istokes),length); 
       i=floor(x/dx)+1; j=floor(y/dy)+1;
       xi = dble(i-1)*dx; yj = dble(j-1)*dy;
       xip1 = xi + dx; yjp1 = yj + dy; 
       ip1 = floor(modulo(xip1,length)/dx)+1; jp1 = floor(modulo(yjp1,length)/dy)+1;
       xip2 = xi + 2.0d0*dx; yjp2 = yj + 2.0d0*dy; 
       ip2 = floor(modulo(xip2,length)/dx)+1; jp2 = floor(modulo(yjp2,length)/dy)+1;
       xim1 = xi - dx; yjm1 = yj - dy; 
       im1 = floor(modulo(xim1,length)/dx)+1; jm1 = floor(modulo(yjm1,length)/dy)+1;

       df00 = (ukx(ip1,j) - ukx(im1,j))/(2.0d0*dx)
       df01 = (ukx(ip1,jp1) - ukx(im1,jp1))/(2.0d0*dx)
       df10 = (ukx(ip2,j) - ukx(i,j))/(2.0d0*dx)
       df11 = (ukx(ip2, jp1) - ukx(i, jp1))/(2.0d0*dx)
       dxux = (xip1 - x)*(yjp1 - y)*df00 + &
            (xip1 - x)*(y - yj)*df01 + &
            (x - xi)*(yjp1 - y)*df10 + &
            (x - xi)*(y - yj)*df11 

       df00 = (uky(ip1,j) - uky(im1,j))/(2.0d0*dx)
       df01 = (uky(ip1,jp1) - uky(im1,jp1))/(2.0d0*dx)
       df10 = (uky(ip2,j) - uky(i,j))/(2.0d0*dx)
       df11 = (uky(ip2, jp1) - uky(i, jp1))/(2.0d0*dx)
       dxuy = (xip1 - x)*(yjp1 - y)*df00 + &
            (xip1 - x)*(y - yj)*df01 + &
            (x - xi)*(yjp1 - y)*df10 + &
            (x - xi)*(y - yj)*df11

       df00 = (ukx(i,jp1) - ukx(i,jm1))/(2.0d0*dx)
       df01 = (ukx(i,jp2) - ukx(i,j))/(2.0d0*dx)
       df10 = (ukx(ip1,jp1) - ukx(ip1,jm1))/(2.0d0*dx)
       df11 = (ukx(ip1, jp2) - ukx(ip1, j))/(2.0d0*dx)
       dyux = (xip1 - x)*(yjp1 - y)*df00 + &
            (xip1 - x)*(y - yj)*df01 + &
            (x - xi)*(yjp1 - y)*df10 + &
            (x - xi)*(y - yj)*df11

    ! Determinant of strain tensor
    dxux = dxux/(dx*dy); 
    dyux = dyux/(dx*dy);
    dxuy = dxuy/(dx*dy);
    lamb(i1,ibead,istokes) = -(dxux**2.0d0) - dxuy*dyux; 
    !! Use incompressibility, dyuy = -dxux
   enddo

    enddo
!$omp end parallel do
 enddo

  end subroutine eval_lam_interp

end module mod_part_interp
