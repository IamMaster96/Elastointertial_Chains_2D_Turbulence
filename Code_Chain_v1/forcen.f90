subroutine forcen(f)
use mod_2dflu
use mod_part_interp
implicit none
!double precision,dimension(:,:),allocatable:: Xp,Yp
double precision,dimension(1:2,1:10000):: f
integer :: nbox,o,p,z1,z2,z3,z4,z5,z6
double precision,dimension(:,:,:),allocatable :: z
double precision,dimension(:,:),allocatable :: q
double precision :: wi,wj,ri,rj,dx1,dy1,nbin,rf,meand,fnp
double precision,dimension(1:2) :: ftemp
nbox=100
allocate(z(0:Nparticle,1:nbox,1:nbox))
!allocate(Xp(Nparticle))
!allocate(f(1:2,1:Nparticle))
!allocate(Yp(Nparticle))
allocate(q(1:5,1:Nparticle))
z=0
f=0
meand=sqrt(tau(istokes)*9.0d0*vis/(2.0d0*1000d0))
nbin=dfloat(nbox)/length
!assign bin number to n particles
do o=1,Nparticle
wi=Xp(o,istokes)
wj=Yp(o,istokes)
z1=abs(int(wi*nbin))+1
z2=abs(int(wj*nbin))+1
q(1,o)=o
q(2,o)=Xp(o,istokes)
q(3,o)=Yp(o,istokes)
q(4,o)=z1
q(5,o)=z2
z(0,z1,z2)=z(0,z1,z2)+1
z(z(0,z1,z2),z1,z2)=o
end do
!force on ith particle
do o=1,Nparticle
wi=q(2,o)
wj=q(3,o)
do z3=-1,1
do z4=-1,1
z5=q(4,o)+z3
z6=q(5,o)+z4
if(z5.gt.nbox)then
z5=z5-nbox
else if(z5.lt.1)then
z5=nbox
else 
z5=z5
end if
if(z6.gt.nbox)then
z6=z6-nbox
else if(z6.lt.1)then
z6=nbox
else 
z6=z6
end if

do p=1,z(0,z5,z6)
ri=q(2,z(p,z5,z6))
rj=q(3,z(p,z5,z6))
dx1=wi-ri
dy1=wj-rj
if(abs(dx1).gt.(length*0.5d0))then
dx1=dx1-length
end if
if(abs(dy1).gt.(length*0.5d0))then
dy1=dy1-length
end if
if(dx1==0.and.dy1==0)then
f(1,o)=f(1,o)+0
f(2,o)=f(2,o)+0
end if
if((dx1.ne.0).or.(dy1.ne.0))then
rf=sqrt(dx1**2+dy1**2)
  if(rf.lt.meand)then
  fnp=(1.0d0-rf/meand)/meand
  ftemp(1)=fnp*(dx1/rf)
  ftemp(2)=fnp*(dy1/rf)
  else
  ftemp(1)=0
  ftemp(2)=0
  end if
f(1,o)=f(1,o)+ftemp(1)
f(2,o)=f(2,o)+ftemp(2)
else
f(1,o)=f(1,o)+0
f(2,o)=f(2,o)+0
end if
end do

end do
end do

end do
return
end subroutine forcen
