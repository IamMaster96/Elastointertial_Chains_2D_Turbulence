      subroutine trumbbell
	use mod_2dflu
	use mod_part_interp
      IMPLICIT NONE
	double precision :: dt
	integer :: j

! The parameters with an exact equivalent in the theory are D and tau. Here mu_st has a different
! meaning as compared to the draft of the paper. The stiffness parameter is A=1/(D*tau).
	dt = delta
	do ip = 1,Nparticle
         grad(1,1)=dxux(ip)   !the velocity gradient coming from NS should go here
         grad(1,2)=dxuy(ip)
         grad(2,1)=dyux(ip)
         grad(2,2)=-grad(1,1)

!        realizations
         
         DO j=1,Nstat

            var1=theta(j)
            var2=chi(j)
            
            CALL integration(var1,var2,D,mu_st,grad,dt)
            
            theta(j)=var1
            chi(j)=var2
         
            CALL statistics(chi(j),Nchi,pdf)

        enddo                  !end realizations
	enddo

         IF(iinner.EQ.navg)THEN !write pdf

            sum_st=0.d0
            DO k=1,Nchi
               sum_st=sum_st+pdf(k)*2.D0*pi/dble(Nchi)
            ENDDO
      
            OPEN(10,file='pdf-chi.dat')
            DO k=1,Nchi
               WRITE(10,*)(dble(k)-0.5d0)*2.D0*pi/Nchi,pdf(k)/sum_st
            ENDDO
            CLOSE(10)
         ENDIF

      end subroutine trumbbell



      SUBROUTINE integration(theta,chi,D,mu_st,grad,dt)

      IMPLICIT NONE

      REAL*8 pi,theta,chi,D,mu_st,grad(2,2),dt,twopi,x
      PARAMETER(pi=acos(-1.d0))
      REAL*8 diff(2,2),drift(2),w1,w2,coeffw
      REAL*8 noise,elastic,flow,coeff1,coeff2,numerator,denominator
      REAL*8 coeff11,coeff12,coeff21,coeff22

      w1=0.d0
      w2=0.d0
      twopi=2.d0*pi

      x=(chi/twopi-INT(chi/twopi))*twopi
      IF(x.LT.0.D0)x=twopi+x

      noise=-6.d0*sin(chi)/((2.d0-cos(chi))*(2.d0+cos(chi))**2.d0)

      elastic=-3.d0*(pi-x)/(2.d0+cos(chi))

      coeff11=2.d0*cos(chi+theta)**2.d0*sin(chi)+cos(theta)*(4.d0*&
          sin(theta)-cos(chi)*sin(chi+theta))
      coeff12=4.d0*sin(theta)**2.d0-(2.d0+cos(chi))*sin(theta)*&
          sin(chi+theta)+2.d0*cos(chi)*sin(chi+theta)**2.d0
      coeff21=-4.d0*cos(theta)**2.d0+(2.d0+cos(chi))*cos(theta)*&
          cos(chi+theta)-2.d0*cos(chi)*cos(chi+theta)**2.d0
      coeff22=2.d0*cos(theta)*(-2.d0*sin(theta)+sin(chi+theta))+&
          cos(chi)*(cos(chi+theta)*sin(theta)-sin(2.d0*(chi+theta)))

      flow=(coeff11*grad(1,1)+coeff12*grad(1,2)+coeff21*grad(2,1)&
          +coeff22*grad(2,2))/(-4.d0+cos(chi)**2.d0)

      drift(1)=D*noise+mu_st*elastic+flow

      noise=12.d0*sin(chi)/((2.d0-cos(chi))*(2.d0+cos(chi))**2.d0)

      elastic=6.d0*(pi-x)/(2.d0+cos(chi))
      
      flow=sin(chi)/(2.d0+cos(chi))*(-2.d0*cos(2.d0*theta+chi)*&
		(grad(1,1) -grad(2,2))-2.d0*sin(2.d0*theta+chi)*&
          (grad(1,2)+grad(2,1)))

      drift(2)=D*noise+mu_st*elastic+flow

      diff(1,1)=24.d0*D/(7.d0-cos(2.d0*chi))
      diff(1,2)=-6.d0*D/(2.d0+cos(chi))
      diff(2,1)=diff(1,2)
      diff(2,2)=12.d0*D/(2.d0+cos(chi))
  
      CALL RANDOM_NUMBER(w1)
      CALL RANDOM_NUMBER(w2)

      coeffw=sqrt(12.d0*dt)

      theta=theta+drift(1)*dt+sqrt(diff(1,1))*coeffw*(w1-0.5d0)

      chi=chi+drift(2)*dt+diff(1,2)/sqrt(diff(1,1))*coeffw*(w1-0.5d0) &
       +sqrt(3.d0*D)*coeffw*(w2-0.5d0)

      RETURN
      END
      
!*********************************************************************
      
      SUBROUTINE statistics(chi,Nchi,pdf)

      IMPLICIT NONE

      INTEGER j,Nchi
      REAL*8 chi,pdf(Nchi),pi,delta,x,twopi
      PARAMETER(pi=acos(-1.d0))

      twopi=2.d0*pi

      delta=twopi/Nchi

      x=(chi/twopi-INT(chi/twopi))*twopi
      IF(x.LT.0.D0)x=twopi+x
      j=ABS(x/delta)  
      IF((x.LT.0).OR.(x.GE.twopi))THEN
         PRINT*, 'NO',j,chi,x
	write(*,*) 1000
         STOP
      ENDIF
      pdf(j+1)=pdf(j+1)+0.001d0

      RETURN
      END
      
      
