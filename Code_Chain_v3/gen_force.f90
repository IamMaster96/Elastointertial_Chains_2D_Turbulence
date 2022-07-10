!!
subroutine gen_force
  use omp_lib
  use mod_2dflu

  implicit none

  integer:: i1,i2
  real*8::x,y

  ! Generate the forcing field in real space

  !$omp parallel default(shared) private(i1,i2,x,y)

  !$omp workshare
  fomega = 0.0d0 ! Initialize the force field
  !$omp end workshare
  ! Checking how many threads are running
  !$omp single
  open ( unit = 111, file = 'stat.out',Access = 'append',Status='old')
  write ( 111, '(a,i8)' ) &
    '  The number of running threads  = ', omp_get_num_threads ( )
  close (111)
  !$omp end single
  
  !$omp do
  do i2 = 1,n2
     y = i2*dy
     do i1 = 1,n1
        x = i1*dx
	fomega(i1,i2) = -famp*kf*dcos(kf*x)
     enddo
  enddo
  !$omp end do

  !$omp end parallel
  !
  ! Fourier tranform the force field
  !
  call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,fomega,0)
end subroutine gen_force
