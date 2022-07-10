subroutine getpsivel
  use omp_lib
  use mod_2dflu
  implicit none
  integer :: i1,i2,k1,k2,ireal,iimag,ksqr
  double precision :: rk2,kx,ky
  !! get psi and velocity from the omega
  factor = 2.0d0*pi/length
!$omp parallel do &
!$omp default(shared) &
!$omp private(i1,i2,ireal,iimag,k1,k2,ksqr,rk2,kx,ky)
  do i2 = 1,n2
     k2 = (i2-1) - n2*(i2/(n1hf+1))	
     do i1 =1,n1hf
        k1 = i1 -1
        ksqr = k1*k1 + k2*k2 
        rk2 = factor*factor*dfloat(ksqr)
        ireal = 2*i1 -1
        iimag = 2*i1
        kx = factor*dfloat(k1)
        ky = factor*dfloat(k2)
!! Dealiasing and setting the base value of psi to zero
        if((ksqr.gt.kasqr).or.(ksqr.eq.0))then
           omega(ireal,i2) = 0.0d0
           omega(iimag,i2) = 0.0d0
           psi(ireal,i2) = 0.0d0
           psi(iimag,i2) = 0.0d0
           ukx(ireal,i2) = 0.0d0
           ukx(iimag,i2) = 0.0d0
           uky(ireal,i2) = 0.0d0
           uky(iimag,i2) = 0.0d0
        else
           psi(ireal,i2) = -omega(ireal,i2)/rk2
           psi(iimag,i2) = -omega(iimag,i2)/rk2
           ukx(ireal,i2) = ky*psi(iimag,i2)
           ukx(iimag,i2) = -ky*psi(ireal,i2)
           uky(ireal,i2) = -kx*psi(iimag,i2)
           uky(iimag,i2) = kx*psi(ireal,i2)
        endif
     enddo
  enddo
!$omp end parallel do
  !! --------------------- 
end subroutine getpsivel
