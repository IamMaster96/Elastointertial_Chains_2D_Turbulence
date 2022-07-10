!! set up the initial configuration for the serial spectral code
!! obeying the K41 spectra. And also divergenceless. 
!! ---------------------------------------------------------------
			subroutine iniconf2d
	use mod_2dflu
	implicit none
	double precision,dimension(3) :: ran
	double precision :: rk2,rk,ek1,vk,p11,p12,p31,p22,p33,p23,rk2inv
	double precision :: iniamp,omegak,n_sqr,rkini,mabs,prim,pert,x,y
	integer :: i1,i2,i3,k1,k2,k3,ksqr,mshl,ik,iseed,ireal,iimag,m1,m2
	integer,dimension(3) :: tarray
	integer,dimension(2) :: tofday
	integer :: ierr
!! -----------------------------------------------------------------

  double precision :: omega_amp=652d0,komega=4.0d0
  double precision :: center1,center2,x1,x2,y1,y2,rad;
!! -----------------------------------------------------------------
!! vorticity has a flat spectra with random phase
        n_sqr = dfloat(n1*n2)
!!        call gettimeofday(tofday,ierr)
!!        iseed = tofday(1)
!! ---------
	iseed = 5
	omega = 0.d0
        do i2 = 1,n2
                k2 = (i2-1) - n2*(i2/(n1hf+1))
                do i1 =1,n1hf
                        k1 = i1 -1
                        ksqr = k1*k1 + k2*k2
                        rk2 = factor*factor*dfloat(ksqr)
                        rk = dsqrt(rk2)
                        mshl = nint(rk)
                        rkini = factor*dfloat(kini)
                        iniamp = (rk2*rk2)*exp(-(rk2*rk2))/den_state(mshl)
                        omegak = dsqrt(iniamp)
                        omegak = omegak*n_sqr
!! -------------
                        call getrandom(1,ran(1),ran(2),ran(3),i1+i2-1,iseed)
                        do ik = 1,3
                        ran(ik)= -pi + 2.0d0*pi*ran(ik)
                        enddo
                        ireal=2*i1-1
                        iimag=2*i1
!! dealias
                        if(ksqr > kasqr)then
                        omega(ireal,i2) = 0.0d0
                        omega(iimag,i2) = 0.0d0
                        else
                        omega(ireal,i2) = omegak*dcos(ran(1))
                        omega(iimag,i2) = omegak*dsin(ran(1))
                        endif

!! ----------------
                enddo
        enddo
!! ------------------------------------------------------------


		end subroutine iniconf2d
