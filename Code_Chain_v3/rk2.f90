!!       rnkt2.f90
!! This is the code for serial runge-kutta iterations for 
!! the spectral code. This does not make any explicit call to
!! the fourier transform libraries, thus may go unchanged from
!! one serial machine to another. 
!! --------------------------------------------------------------
subroutine rnkt
  use mod_2dflu
  implicit none
  double precision,dimension(3) :: ran
  double precision :: rk2,rk2inv,delta1,tr,tc,temph,omre,omim
  double precision :: encheck,kx,ky,n_sqr
  integer :: i1,i2,k1,k2,ksqr,ireal,iimag,j,ik
  integer,dimension(3) :: tarray
  integer,dimension(2) :: tofday
  integer :: ierr
  !! ----------------
  factor = 2.0d0*pi/length 
  !! ---------------------------------------------------------
  !! first  store omega in temporary  array.  
  jac_old = omega	
  !! -evaluate the Jacobean ----------------------------
  call eval_jac(1) 
  !! -----------------------------------------------------------
  n_sqr = dfloat(n1*n2)
  encheck = 0.0d0
!$omp parallel do &
!$omp default(shared) &
!$omp private(i1,i2,ireal,iimag,k1,k2,ksqr,rk2,kx,ky,temph,tr,tc,omre,omim)
  do i2 = 1,n2
     k2 = (i2-1) - n2*(i2/(n1hf+1))	
     do i1 =1,n1hf
        k1 = i1 -1
        ireal = 2*i1-1
        iimag=2*i1
        ksqr = k1*k1 + k2*k2 
        if((ksqr.gt.kasqr).or.(ksqr.eq.0))then
           omega(ireal,i2) = 0.0d0
           omega(iimag,i2) = 0.0d0 
        else
           rk2 = factor*factor*dfloat(ksqr)
           kx = factor*dfloat(k1)
           ky = factor*dfloat(k2)
           temph = time_increment(ksqr)
           !! -------------			
           tr = (ky*uky(iimag,i2)+kx*ukx(iimag,i2))  
           tc = -(ky*uky(ireal,i2)+ kx*ukx(ireal,i2))
           omre = jac_old(ireal,i2)
           omim = jac_old(iimag,i2)
           !! ---time stepping  -----------------------------------
           !! ------------------------------------------------------------
           omega(ireal,i2)=temph*(omre+(delta/2.0d0)*(tr + fomega(ireal,i2))) 
           omega(iimag,i2)=temph*(omim+(delta/2.0d0)*(tc + fomega(iimag,i2))) 
           !! ------------------------
        endif
     enddo
  enddo
!$omp end parallel do
  !! -the second step ----------------------------------
  !! -------evaluating the non-linear part ----------------
  call eval_jac(2)
  !! --the second time-stepping -------------------------------------
  encheck = 0.0d0
!$omp parallel do &
!$omp default(shared) &
!$omp private(i1,i2,ireal,iimag,k1,k2,ksqr,rk2,kx,ky,temph,tr,tc,omre,omim)
  do i2 = 1,n2
     k2 = (i2-1) - n2*(i2/(n1hf+1))	
     do i1 =1,n1hf
        k1 = i1 -1
        ksqr = k1*k1 + k2*k2 
        rk2 = factor*factor*dfloat(ksqr)
        ireal = 2*i1-1
        iimag=2*i1
        if((ksqr.gt.kasqr).or.(ksqr.eq.0))then
           omega(ireal,i2) = 0.0d0
           omega(iimag,i2) = 0.0d0 
        else
           rk2 = factor*factor*dfloat(ksqr)
           kx = factor*dfloat(k1)
           ky = factor*dfloat(k2)
           temph = time_increment(ksqr)
           tr = (ky*uky(iimag,i2)+kx*ukx(iimag,i2))  
           tc = -(ky*uky(ireal,i2)+ kx*ukx(ireal,i2))
           omre = jac_old(ireal,i2)
           omim = jac_old(iimag,i2)
           !! ---time stepping  -----------------------------------
           omega(ireal,i2) =  temph*(temph*omre+delta*(tr + fomega(ireal,i2))) 
           omega(iimag,i2) =  temph*(temph*omim+delta*(tc + fomega(iimag,i2)))
           !! ---------------------------------------------------
        endif
     enddo
  enddo
  !$omp end parallel do

end subroutine rnkt
