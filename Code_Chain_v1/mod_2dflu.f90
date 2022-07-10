!!             mod_2dflu.f90
!! This is the module containing the defination and commons for the 
!! serial code. As all the serial ffts requal arrays of same type 
!! we expect this should not be changed as we go from one machine  
!! to another.
module mod_2dflu
  implicit none
  save
  !! ----velocites  -------------
  double precision,allocatable,dimension(:,:) :: psi,omega, & 
       jac_old,ukx,uky, &
       fomega,uktmp 
  double precision,allocatable,dimension(:,:) :: omega_tm
  !! ---parameters of simulations -----------------
  double precision,allocatable,dimension(:,:) :: E_Omega
  integer :: pmax
  !! -----common parameters ----------------------------
  integer :: thousand,nn,maxiter,nrun,nrunpart,n1,n2,n1h,n1hf, &
       nshell,navg,ksqr_max,mm,nthreads,iinner,counter,navgpart
  integer :: kini,kdrag,kasqr,kf,no_of_fmode,t_n,ip
  double precision :: vis,vis2,delta,pi,length,factor,famp,mu,mu2,finp
  double precision :: fixed_energy1,fixed_energy2,fixed_energy3
  real*8::dx,dy
  real*8::Req_chain,Rmax_chain, Rini_chain,Req,Rmax, Rini
  integer :: Nbead
  !! -----global arrays --------------
  double precision,allocatable,dimension(:) :: time_increment,alpha, &
	diss_rate,tot_energy,energy_time, &
	enstrphy_dissrate,energy_dissrate
  integer,allocatable,dimension(:) :: den_state
  integer,allocatable,dimension(:,:) :: fmode 
  !! ----------------------------------------------------------
  !! ---parameters of simulations -----------------
  double precision :: edr,tlr_micro_scale,dissipation_scale,vkrms,Rlambda, &
       tau_eddy,integral_scale,Rbox,Rintscale
  !! -----------For FFTW --------------------------------------
  integer*8 :: pfor,pinv
  integer,dimension(2) :: dim
  double precision :: scale
  include "fftw_f77.i"
  !! ---------------------------------------------------------
!! Dario's code
      INTEGER Nt,Nchi,Nstat,k
      REAL*8 D,mu_st,sum_st,random,var1,var2
      REAL*8 grad(2,2)
      PARAMETER(Nchi=100,Nstat=200) !Nchi number of bins, Nstat number of realizations
      REAL*8 theta(Nstat),chi(Nstat)
      REAL*8 pdf(Nchi)
      EXTERNAL integration,statistics

end module mod_2dflu
