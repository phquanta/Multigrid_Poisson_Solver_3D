MODULE LevelType
TYPE :: LevelStruct
 INTEGER:: ndim
 INTEGER:: downsweep
 INTEGER:: upsweep
 INTEGER:: PlaneMaxDim
 INTEGER:: LineMaxDim
 INTEGER:: BuffElements
 REAL*8 :: dVol
 REAL*8,pointer::x(:)
 REAL*8,pointer::y(:)
 REAL*8,pointer::z(:)
 REAL*8,pointer::rhs(:,:,:)
 REAL*8,pointer::phi(:,:,:)

 REAL*8,pointer::eps(:,:,:,:)
 REAL*8,pointer::err(:,:,:)


 REAL*8,pointer::PlaneBuffRecv(:,:)
 REAL*8,pointer::LineBuffRecV(:,:,:)
 REAL*8,pointer::PointBuffRecV(:,:,:)

 INTEGER, pointer:: nElements(:)
 INTEGER, pointer:: nb(:)
 INTEGER, pointer:: X1(:)
 INTEGER, pointer:: X2(:)
 INTEGER, pointer:: Y1(:)
 INTEGER, pointer:: Y2(:)
 INTEGER, pointer:: Z1(:)		
 INTEGER, pointer:: Z2(:)


 INTEGER, pointer:: RecV(:,:)
 INTEGER, pointer:: SenD(:,:)
 INTEGER, pointer:: PlanesSend(:,:)

 INTEGER, pointer:: X1Plane(:,:)
 INTEGER, pointer:: X2Plane(:,:)
 INTEGER, pointer:: Y1Plane(:,:)
 INTEGER, pointer:: Y2Plane(:,:)
 INTEGER, pointer:: Z1Plane(:,:)		
 INTEGER, pointer:: Z2Plane(:,:)



 INTEGER, pointer:: RecV1(:,:,:)
 INTEGER, pointer:: SenD1(:,:,:)
 INTEGER, pointer:: LinesSend(:,:,:)


 INTEGER, pointer:: RecV2(:,:,:)
 INTEGER, pointer:: SenD2(:,:,:)
 INTEGER, pointer:: PointsSend(:,:,:)

 INTEGER, pointer:: X1Line(:,:,:)
 INTEGER, pointer:: X2Line(:,:,:)
 INTEGER, pointer:: Y1Line(:,:,:)
 INTEGER, pointer:: Y2Line(:,:,:)
 INTEGER, pointer:: Z1Line(:,:,:)		
 INTEGER, pointer:: Z2Line(:,:,:)

 INTEGER, pointer:: XPoint(:,:,:)
 INTEGER, pointer:: YPoint(:,:,:)
 INTEGER, pointer:: ZPoint(:,:,:)		

END TYPE LevelStruct
END Module LevelType

PROGRAM POISSON_MGRID
Use LevelType
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER MYID,TOTPS,IERR,stat(MPI_STATUS_SIZE),NCPUS

!########## GAUSS UNITS IMPLIED ############
INTEGER I,J,K,L,M,P,NDIM,ITER,MAXITER,MAXDIM,Ip,Im,Jp,Jm,Kp,Km,NDIM2H,NDIM4H,NLEVELS,Jnow
INTEGER :: ndim1,i1,j1,k1,PlaneMaxDim,NBb,II,JJ,KK,MM
REAL*8 A2B,lamCoarse


REAL*8 Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,dhX,dhY,dhZ,summ_prev,summE,summEprev,summE1
REAL*8 coeff,d_map,const,R(6),summ,Xoffset,Yoffset,Zoffset,Rpol,summRho,vol,dvol
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: PHI,RHO
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: Rk,Rkt,Rk1,Rkt1,Rkp1,Pk,PkT,Rtemp,Rk2H,PHI2H,RHO2H,Rk4H,PHI4H,Rtemp2H

REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: EPS
REAL*8 summ_even,summ_odd,omega,coeff1,coeff2,delta,Pi,summHartree,summHartree1,Xtemp,Ytemp,lamN,Momega
REAL*8 eps1,eps2,eps3,eps4,eps5,eps6,Ztemp,temp1,temp2,temp3,temp,alphaK,betaK,dVol2H,summRhoTemp
REAL*8 lamN1,lamN2,summRes(1000)
REAL t_start(2),t_stop(2),result
REAL, DIMENSION(:), ALLOCATABLE :: X,Y,Z,Xmid,Ymid,Zmid,X2H,Y2H,Z2H
PARAMETER(A2B=1./0.529177249)
PARAMETER(NDIM=120)
PARAMETER(NDIM2H=NDIM/2)
PARAMETER(NDIM4H=NDIM2H/2)
PARAMETER(Xoffset=4.)
PARAMETER(Yoffset=12.)
PARAMETER(Zoffset=5.)
PARAMETER(Rpol=3.)
PARAMETER(MAXITER=22242)



TYPE (LevelStruct), DIMENSION(:), ALLOCATABLE  ::Level




 
open (UNIT=85, FILE='eps', STATUS='unknown')



open (UNIT=86, FILE='contour', STATUS='unknown')
open (UNIT=84, FILE='potential', STATUS='unknown')
open (UNIT=83, FILE='field_magniture', STATUS='unknown')
open (UNIT=101, FILE='field', STATUS='unknown')

open (UNIT=103, FILE='ratioHelm_par.dat', STATUS='unknown')
open (UNIT=104, FILE='residiual.dat', STATUS='unknown')


Pi=4.*Atan(1.)


 call mpi_init( IERR)
 call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)
 call MPI_COMM_SIZE(MPI_COMM_WORLD,NCPUS,IERR)


ALLOCATE(EPS(6,NDIM,NDIM,NDIM))
ALLOCATE(RHO(NDIM,NDIM,NDIM))
ALLOCATE(PHI(NDIM,NDIM,NDIM))



ALLOCATE(X(NDIM))
ALLOCATE(Y(NDIM))
ALLOCATE(Z(NDIM))





lamN1=0.99
lamN2=0.99
lamCoarse=0.92

!lamN1=0.993
!lamN2=0.993
!lamCoarse=0.92


  call GetNLevels(NDIM,nlevels)
!  nlevels=2
!   nlevels=3
  ALLOCATE(Level(nlevels))

  ndim1=ndim
  DO i=nlevels,1,-1
  write(*,*) "NLEVELS=",nlevels,ndim1
   Level(i)%ndim=ndim1
   ALLOCATE(level(i)%x(ndim1))
   ALLOCATE(level(i)%y(ndim1))
   ALLOCATE(level(i)%z(ndim1))
   ALLOCATE(level(i)%rhs(ndim1,ndim1,ndim1))
   ALLOCATE(level(i)%err(ndim1,ndim1,ndim1))
   ALLOCATE(level(i)%phi(ndim1,ndim1,ndim1))
   ALLOCATE(level(i)%eps(6,ndim1,ndim1,ndim1))
   ALLOCATE(level(i)%X1(NCPUS))
   ALLOCATE(level(i)%X2(NCPUS))
   ALLOCATE(level(i)%Y1(NCPUS))
   ALLOCATE(level(i)%Y2(NCPUS))
   ALLOCATE(level(i)%Z1(NCPUS))
   ALLOCATE(level(i)%Z2(NCPUS))
   ALLOCATE(level(i)%Nb(NCPUS))
   ALLOCATE(level(i)%Nelements(NCPUS))
   ALLOCATE(level(i)%RecV(NCPUS,6))
   ALLOCATE(level(i)%SenD(NCPUS,6))
   ALLOCATE(level(i)%PlanesSend(NCPUS,6))
   ALLOCATE(level(i)%X1Plane(NCPUS,6))
   ALLOCATE(level(i)%X2Plane(NCPUS,6))
   ALLOCATE(level(i)%Y1Plane(NCPUS,6))
   ALLOCATE(level(i)%Y2Plane(NCPUS,6))
   ALLOCATE(level(i)%Z1Plane(NCPUS,6))
   ALLOCATE(level(i)%Z2Plane(NCPUS,6))
   ALLOCATE(level(i)%RecV1(NCPUS,6,4))
   ALLOCATE(level(i)%SenD1(NCPUS,6,4))
   ALLOCATE(level(i)%LinesSend(NCPUS,6,4))
   ALLOCATE(level(i)%PointsSend(NCPUS,6,4))
   ALLOCATE(level(i)%X1Line(4,NCPUS,6))
   ALLOCATE(level(i)%X2Line(4,NCPUS,6))
   ALLOCATE(level(i)%Y1Line(4,NCPUS,6))
   ALLOCATE(level(i)%Y2Line(4,NCPUS,6))
   ALLOCATE(level(i)%Z1Line(4,NCPUS,6))
   ALLOCATE(level(i)%Z2Line(4,NCPUS,6))
   ALLOCATE(level(i)%RecV2(NCPUS,6,4))
   ALLOCATE(level(i)%SenD2(NCPUS,6,4))
   ALLOCATE(level(i)%PointsSend(NCPUS,6,4))
   ALLOCATE(level(i)%XPoint(4,NCPUS,6))
   ALLOCATE(level(i)%YPoint(4,NCPUS,6))
   ALLOCATE(level(i)%ZPoint(4,NCPUS,6))
   level(i)%phi=0.
   level(i)%err=0.
   level(i)%rhs=0.
   level(i)%downsweep=3
   level(i)%upsweep=4
   ndim1=ndim1/2
  ENDDO
 

  call INITIALIZE(MYID,NCPUS,LEVEL(nlevels),LEVEL(nlevels),nlevels,nlevels)
  DO i=nlevels-1,1,-1
   call INITIALIZE(MYID,NCPUS,LEVEL(i),LEVEL(i+1),i,nlevels)
  ENDDO



  DO i=nlevels,1,-1
   call GetNelements(NCPUS,LEVEL(i))
   ALLOCATE(level(i)%PlaneBuffRecV(level(i)%PlaneMaxDim+7,level(i)%NB(MYID+1)))
   ALLOCATE(level(i)%LineBuffRecV(level(i)%LineMaxDim+7,level(i)%NB(MYID+1),4))
   ALLOCATE(level(i)%PointBuffRecV(5,level(i)%NB(MYID+1),4))
   level(i)%BUFFElements=MAXVAL(level(i)%NElements)
   IF(MYID.eq.0) write(*,*) "PlaneMaxDim,levelN ", Level(i)%PlaneMaxDim,level(i)%BUFFELEMENTS,i
  ENDDO

  call ReadData(vol,LEVEL(nlevels))

  DO i=nlevels-1,1,-1
    call propagateXYZ(LEVEL(I+1),LEVEL(I))	
  ENDDO


  DO i=nlevels,1,-1
    call fillEps(level(i),Xoffset,Yoffset,Zoffset)
  ENDDO



summ=0.
summ_prev=0.

level(nlevels)%rhs(:,:,:)=level(nlevels)%rhs(:,:,:)*level(nlevels)%eps(2,:,:,:)
level(nlevels)%phi=0.



!DO I=1,NDIM
! DO J=1,NDIM
!  DO K=1,NDIM
!	IF(MYID.eq.0) write(84,*)  level(nlevels)%rhs(i,j,k)
!ENDDO
!ENDDO
!ENDDO

123 continue	
write(*,*) "**********************************"
  call etime (t_start,result)
Jnow=0
!     call  MPI_BARRIER(MPI_COMM_WORLD,IERR)
DO J=1,90
  call  MPI_BARRIER(MPI_COMM_WORLD,IERR)
  call GaussSeidelSmoothSOR(LEVEL(nlevels),level(nlevels)%downsweep,Vol,1,1,lamN1,NCPUS,MYID,MPI_COMM_WORLD)
  call computeDefectpar(MYID,NCPUS,LeveL(nlevels),MPI_COMM_WORLD,summE,1)
  call MPI_ALLREDUCE(summE,summE1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
  summRes(J)=summE1
  IF(DSQRT(summRes(J)).le.1.d-5)  goto 122
 IF(MYID.eq.0) write(*,*) DSQRT(summE1),J

 DO I=Nlevels-1,1,-1
     level(i)%phi=0.
    if(i.eq.1) then
     call restrictDefectpar(MYID,Level(i+1),Level(i),MPI_COMM_WORLD)
     call GaussSeidelSmoothSOR(LEVEL(i),22,Vol,0,1,lamCoarse,NCPUS,MYID,MPI_COMM_WORLD)
    else
     call restrictDefectpar(MYID,Level(i+1),Level(i),MPI_COMM_WORLD)
     call GaussSeidelSmoothSOR(LEVEL(i),LEVEL(i)%downsweep,Vol,1,1,lamN1,NCPUS,MYID,MPI_COMM_WORLD)
     call computeDefectpar(MYID,NCPUS,LeveL(i),MPI_COMM_WORLD,summE,0)
    endif
  ENDDO


  DO I=1,Nlevels-1
	call prolongateDefectpar(MYID,Level(i+1),Level(i),MPI_COMM_WORLD)
        call GaussSeidelSmoothSOR(LEVEL(i+1),LEVEL(i+1)%upsweep,Vol,0,1,lamN2,NCPUS,MYID,MPI_COMM_WORLD)
  ENDDO
ENDDO



122 continue
 call etime (t_stop,result )
 IF(MYID.eq.0) write(*,*) "Elapsed time", t_stop(1)-t_start(1),t_stop(2)-t_start(2)
 dhX=ABS(level(nlevels)%x(2)-level(nlevels)%x(1))
 summ=0.
 DO K=LEVEL(nlevels)%Z1(MYID+1),LEVEL(nlevels)%Z2(MYID+1)
  DO J=LEVEL(nlevels)%Y1(MYID+1),LEVEL(nlevels)%Y2(MYID+1)
   DO I=LEVEL(nlevels)%X1(MYID+1),LEVEL(nlevels)%X2(MYID+1)
      summ=summ+level(nlevels)%phi(i,j,k)*level(nlevels)%eps(2,i,j,k)*level(nlevels)%rhs(i,j,k)/(4.*PI)*dhX*0.5
   ENDDO
  ENDDO
 ENDDO
 call MPI_ALLREDUCE(summ,summE,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 IF(MYID.eq.0) write(*,*) "Energy=",summE
 IF(MYID.eq.0) write(*,*) "I'm here"




DEALLOCATE(EPS)
DEALLOCATE(PHI)
DEALLOCATE(X)
DEALLOCATE(Y)
DEALLOCATE(Z)
DEALLOCATE(RHO)

  DO i=nlevels,1,-1


   DEALLOCATE(level(i)%x)
   DEALLOCATE(level(i)%y)
   DEALLOCATE(level(i)%z)
   DEALLOCATE(level(i)%rhs)
   DEALLOCATE(level(i)%err)
   DEALLOCATE(level(i)%phi)
   DEALLOCATE(level(i)%eps)


   DEALLOCATE(level(i)%X1)
   DEALLOCATE(level(i)%X2)
   DEALLOCATE(level(i)%Y1)
   DEALLOCATE(level(i)%Y2)
   DEALLOCATE(level(i)%Z1)
   DEALLOCATE(level(i)%Z2)

   DEALLOCATE(level(i)%Nb)
   DEALLOCATE(level(i)%Nelements)

   DEALLOCATE(level(i)%RecV)
   DEALLOCATE(level(i)%SenD)
   DEALLOCATE(level(i)%PlanesSend)

   DEALLOCATE(level(i)%X1Plane)
   DEALLOCATE(level(i)%X2Plane)
   DEALLOCATE(level(i)%Y1Plane)
   DEALLOCATE(level(i)%Y2Plane)
   DEALLOCATE(level(i)%Z1Plane)
   DEALLOCATE(level(i)%Z2Plane)


   DEALLOCATE(level(i)%RecV1)
   DEALLOCATE(level(i)%SenD1)
   DEALLOCATE(level(i)%LinesSend)

   DEALLOCATE(level(i)%X1Line)
   DEALLOCATE(level(i)%X2Line)
   DEALLOCATE(level(i)%Y1Line)
   DEALLOCATE(level(i)%Y2Line)
   DEALLOCATE(level(i)%Z1Line)
   DEALLOCATE(level(i)%Z2Line)


   DEALLOCATE(level(i)%RecV2)
   DEALLOCATE(level(i)%SenD2)
   DEALLOCATE(level(i)%PointsSend)

   DEALLOCATE(level(i)%XPoint)
   DEALLOCATE(level(i)%YPoint)
   DEALLOCATE(level(i)%ZPoint)




   DEALLOCATE(level(i)%PlaneBuffRecV)
   DEALLOCATE(level(i)%LineBuffRecV)
   DEALLOCATE(level(i)%PointBuffRecV)

  ENDDO
   DEALLOCATE(level)
END PROGRAM POISSON_MGRID






SUBROUTINE computeDefectpar(MYID,NCPUS,LeveL,COMM,summE,flag)
Use LevelType
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER :: MYID,IERR,stat(MPI_STATUS_SIZE),COMM,JJ,II,KK,Df1,Df2,PNT,CNT,Df3
INTEGER :: NCPUS,PlaneMaxDim,NBI,FLAG1,NN,FLAG2
REAL*8  :: phiIm,phiIp,phiJm,phiJp,phiKm,phiKp,summE
INTEGER :: CpuN
INTEGER :: NDIM,ITER,MAXITER,I,J,K,IP,IM,JP,JM,KP,KM,P,NDIM2,CNT2,CNT1,flag,M,O,CNTTT

REAL*8  :: SUMMHARTREE,EPS1,D_MAP,const,diag,offdiag,tmp
TYPE (LevelStruct)::Level



 CpuN=MYID+1

NDIM=LEVEL%NDIM


! call UpdateBoundary(MYID,Level,COMM,0)
summE=0.
IF(FLAG.eq.1) THEN
DO K=LEVEL%Z1(CpuN),LEVEL%Z2(CpuN)
  DO J=LEVEL%Y1(CpuN),LEVEL%Y2(CpuN)
	DO I=LEVEL%X1(CpuN),LEVEL%X2(CpuN)
	 ip=I+1
	 im=I-1
         jp=J+1
	 jm=J-1
	 kp=K+1
	 km=K-1
	 IF(I.EQ.NDIM) ip=1
	 IF(J.EQ.NDIM) jp=1
	 IF(K.EQ.NDIM) kp=1
	 IF(I.EQ.1) im=NDIM
	 IF(J.EQ.1) jm=NDIM
	 IF(K.EQ.1) km=NDIM	

         phiIm=LEVEL%phi(im,j,k)
	 phiIp=LEVEL%phi(ip,j,k)
	 phiJm=LEVEL%phi(i,jm,k)		
	 phiJp=LEVEL%phi(i,jp,k)		
	 phiKm=LEVEL%phi(i,j,km)
	 phiKp=LEVEL%phi(i,j,kp)
	 eps1=LEVEL%eps(5,i,j,k)
         const=(6.+eps1)
	 offdiag=-(phiIm+phiIp+phiJm+phiJp+phiKm+phiKp)
	 diag=const*LEVEL%phi(i,j,k)
	 tmp=LEVEL%rhs(I,J,K)-(diag+offdiag)	
	 summE=summE+tmp*tmp
	 LEVEL%err(I,J,K)=tmp
    ENDDO
   ENDDO
 ENDDO

ELSE
DO K=LEVEL%Z1(CpuN),LEVEL%Z2(CpuN)
  DO J=LEVEL%Y1(CpuN),LEVEL%Y2(CpuN)
	DO I=LEVEL%X1(CpuN),LEVEL%X2(CpuN)
	  ip=I+1
	  im=I-1
          jp=J+1
	  jm=J-1
	  kp=K+1
	  km=K-1
	  IF(I.EQ.NDIM) ip=1
	  IF(J.EQ.NDIM) jp=1
	  IF(K.EQ.NDIM) kp=1
	  IF(I.EQ.1) im=NDIM
	  IF(J.EQ.1) jm=NDIM
	  IF(K.EQ.1) km=NDIM	
          phiIm=LEVEL%phi(im,j,k)
	  phiIp=LEVEL%phi(ip,j,k)
	  phiJm=LEVEL%phi(i,jm,k)		
	  phiJp=LEVEL%phi(i,jp,k)		
	  phiKm=LEVEL%phi(i,j,km)
	  phiKp=LEVEL%phi(i,j,kp)
	  eps1=LEVEL%eps(5,i,j,k)
          const=(6.+eps1)
	  offdiag=-(phiIm+phiIp+phiJm+phiJp+phiKm+phiKp)
	  diag=const*LEVEL%phi(i,j,k)
  	  LEVEL%err(I,J,K)=LEVEL%rhs(I,J,K)-(diag+offdiag)	
    ENDDO
   ENDDO
 ENDDO
ENDIF

END SUBROUTINE computeDefectpar



SUBROUTINE restrictDefectpar(MYID,Level,Level1,COMM)
Use LevelType
IMPLICIT NONE
INCLUDE 'mpif.h'
REAL*8 summE,summE1
INTEGER NDIM,NDIM2H,I,J,K,Ih,Jh,Kh,Im,Jm,Km,Ip,Jp,Kp,MYID,CpuN,COMM
REAL*8 temp,tmp1
TYPE (LevelStruct)::Level,Level1

NDIM=LEVEL%NDIM
NDIM2H=LEVEL1%NDIM
 CpuN=MYID+1

 call UPDATE(MYID,LEVEL,COMM,1)
 call UpdateBoundary(MYID,Level,COMM,1)

summE=0.
temp=1./64.

 DO K=LEVEL1%Z1(CpuN),LEVEL1%Z2(CpuN)
   DO J=LEVEL1%Y1(CpuN),LEVEL1%Y2(CpuN)
	DO I=LEVEL1%X1(CpuN),LEVEL1%X2(CpuN)
	Ih=2*I-1
	Jh=2*J-1
	Kh=2*K-1
    	Ip=MOD(Ih,NDIM)+1
	Im=Ih-1
        Jp=MOD(Jh,NDIM)+1
        Jm=Jh-1
        Kp=MOD(Kh,NDIM)+1
        Km=Kh-1
        IF(Ih.EQ.1) Im=NDIM
        IF(Jh.EQ.1) Jm=NDIM
        IF(Kh.EQ.1) Km=NDIM


tmp1=(4.*LEVEL%ERR(Ih,Jh,Km)+2.*LEVEL%ERR(Im,Jh,Km)+2.*LEVEL%ERR(Ip,Jh,Km)+2.*LEVEL%ERR(Ih,Jm,Km)+2.*LEVEL%ERR(Ih,Jp,Km)+LEVEL%ERR(Im,Jp,Km)&
&+LEVEL%ERR(Im,Jm,Km)+LEVEL%ERR(Ip,Jm,Km)+LEVEL%ERR(Ip,Jp,Km))
tmp1=(8.*LEVEL%ERR(Ih,Jh,Kh)+4.*LEVEL%ERR(Im,Jh,Kh)+4.*LEVEL%ERR(Ip,Jh,Kh)+4.*LEVEL%ERR(Ih,Jm,Kh)+4.*LEVEL%ERR(Ih,Jp,Kh)+2.*LEVEL%ERR(Im,Jp,Km)&
& +2.*LEVEL%ERR(Im,Jm,Km)+2.*LEVEL%ERR(Ip,Jm,Km)+2.*LEVEL%ERR(Ip,Jp,Km))+tmp1

tmp1=(4.*LEVEL%ERR(Ih,Jh,Kp)+2.*LEVEL%ERR(Im,Jh,Kp)+2.*LEVEL%ERR(Ip,Jh,Kp)+2.*LEVEL%ERR(Ih,Jm,Kp)+2.*LEVEL%ERR(Ih,Jp,Kp)+LEVEL%ERR(Im,Jp,Kp)&
&+LEVEL%ERR(Im,Jm,Kp)+LEVEL%ERR(Ip,Jm,Kp)+LEVEL%ERR(Ip,Jp,Kp))+tmp1
LEVEL1%RHS(I,J,K)=tmp1*temp
!LEVEL1%RHS(I,J,K)=1.
     ENDDO
   ENDDO
ENDDO


END SUBROUTINE restrictDefectpar





SUBROUTINE prolongateDefectpar(MYID,Level,Level1,COMM)
Use LevelType

IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER NDIM,NDIM2H,I1,J1,K1,Ih,Jh,Kh,Im,Jm,Km,Ip,Jp,Kp,I,J,K,MYID,CpuN,COMM

TYPE (LevelStruct)::Level,Level1
REAL*8 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,PhiIp,PhiJp,PhiKp,PhiIpJp,PhiIpKp,PhiJpKp
REAL*8 Half,Fourth,Eighth
 CpuN=MYID+1

NDIM2H=Level1%NDIM
NDIM=Level%NDIM

Half=1./2.
Fourth=1./4.
Eighth=1./8.
 DO K=LEVEL1%Z1(CpuN),LEVEL1%Z2(CpuN)
   DO J=LEVEL1%Y1(CpuN),LEVEL1%Y2(CpuN)
	DO I=LEVEL1%X1(CpuN),LEVEL1%X2(CpuN)

	Ih=2*I-1
	Jh=2*J-1
	Kh=2*K-1
	
    	Ip=MOD(I,NDIM2H)+1

        Jp=MOD(J,NDIM2H)+1

        Kp=MOD(K,NDIM2H)+1

	PhiIp=LEVEL1%phi(Ip,J,K)
	PhiJp=LEVEL1%phi(I,Jp,K)
	PhiKp=LEVEL1%phi(I,J,Kp)
	PhiJpKp=LEVEL1%phi(I,Jp,Kp)
	PhiIpKp=LEVEL1%phi(Ip,J,Kp)
	PhiIpJp=LEVEL1%phi(Ip,Jp,K)
	tmp1=LEVEL1%phi(I,J,K)
	tmp2=Half*(tmp1+PhiIp)
	tmp3=Half*(tmp1+PhiJp)
	tmp4=Half*(tmp1+PhiKp)

	tmp5=Fourth*(tmp1+PhiKp+PhiJp+PhiJpKp)
	tmp6=Fourth*(tmp1+PhiIp+PhiJp+PhiIpJp)
	tmp7=Fourth*(tmp1+PhiIp+PhiKp+PhiIpKp)

	tmp8=Eighth*(tmp1+PhiIp+PhiJp+PhiKp+PhiIpJp+PhiIpKp+PhiJpKp+LEVEL1%phi(Ip,Jp,Kp))


	LEVEL%phi(Ih,Jh,Kh)=LEVEL%phi(Ih,Jh,Kh)+tmp1
	LEVEL%phi(Ih+1,Jh,Kh)=LEVEL%phi(Ih+1,Jh,Kh)+tmp2
	LEVEL%phi(Ih,Jh+1,Kh)=LEVEL%phi(Ih,Jh+1,Kh)+tmp3
	LEVEL%phi(Ih,Jh,Kh+1)=LEVEL%phi(Ih,Jh,Kh+1)+tmp4

	LEVEL%phi(Ih,Jh+1,Kh+1)=LEVEL%phi(Ih,Jh+1,Kh+1)+tmp5
	LEVEL%phi(Ih+1,Jh+1,Kh)=LEVEL%phi(Ih+1,Jh+1,Kh)+tmp6
	LEVEL%phi(Ih+1,Jh,Kh+1)=LEVEL%phi(Ih+1,Jh,Kh+1)+tmp7

	LEVEL%phi(Ih+1,Jh+1,Kh+1)=LEVEL%phi(Ih+1,Jh+1,Kh+1)+tmp8

     ENDDO
   ENDDO
ENDDO


END SUBROUTINE prolongateDefectpar








SUBROUTINE GaussSeidelSmoothSOR(LEVEL,MAXITER,Vol,flag,flag1,lamN,NCPUS,MYID,COMM)
Use LevelType
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER MYID,IERR,stat(MPI_STATUS_SIZE),COMM,JJ,II,KK,Df1,Df2,PNT,CNT,Df3
INTEGER NCPUS,PlaneMaxDim,NBI,FLAG1,NN,FLAG2
REAL*8  :: phiIm,phiIp,phiJm,phiJp,phiKm,phiKp,summE
INTEGER:: XrecV(7,6),CpuN,Nelements(NCPUS),firstTimeflagX,firstTimeflagY,firstTimeflagZ,X1s,X2s,Y1s,Y2s,Z1s,Z2s
INTEGER NDIM,ITER,MAXITER,I,J,K,IP,IM,JP,JM,KP,KM,P,NDIM2,CNT2,CNT1,flag,M,O,CNTTT
REAL*8 SUMMHARTREE,EPS1,EPS2,EPS3,EPS4,EPS5,EPS6,D_MAP,summ,summ_prev,error
REAL*8 Vol,coeff2,omega,coeff1,summ1,IND(7),lamN,Momega,summHRecv,summRecV,summE1
INTEGER,SAVE:: FirstTime=1
LOGICAL FlagRB
TYPE (LevelStruct)::Level

NDIM=LEVEL%NDIM
omega=2./(1.+SQRT(1.-lamN))
Momega=1.-omega
 CpuN=MYID+1
summRecV=0.
summHRecV=0.

summE=0.
 DO K=LEVEL%Z1(CpuN),LEVEL%Z2(CpuN)
   DO J=LEVEL%Y1(CpuN),LEVEL%Y2(CpuN)
	DO I=LEVEL%X1(CpuN),LEVEL%X2(CpuN)
	summE=summE+abs(LEVEL%rhs(i,j,k))
 ENDDO
ENDDO
ENDDO

  call MPI_ALLREDUCE(summE,summE1,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM,IERR)
IF(MYID.eq.0) write(*,*) "SMOTRI=ABS(PHI)=",summE1

DO ITER=1,MAXITER
! call  MPI_BARRIER(COMM,IERR)
 call UPDATE(MYID,LEVEL,COMM,0)
 summHartree=0.
 summ_prev=summRecV
 summ=0.
 FlagRB=.TRUE.
 call ITERATE(LEVEL,FlagRB,lamN,MYID,COMM,0)
 FlagRB=.FALSE.
 call  MPI_BARRIER(COMM,IERR)
 call UPDATE(MYID,LEVEL,COMM,0)
 call ITERATE(LEVEL,FlagRB,lamN,MYID,COMM,0)
  DO K=LEVEL%Z1(CpuN),LEVEL%Z2(CpuN)
    DO J=LEVEL%Y1(CpuN),LEVEL%Y2(CpuN)
      DO I=LEVEL%X1(CpuN),LEVEL%X2(CpuN)
	coeff2=LEVEL%phi(i,j,k)
 	summ=abs(LEVEL%phi(i,j,k))+summ
	summHartree=summHartree+coeff2*LEVEL%eps(2,i,j,k)
	ENDDO
    ENDDO
ENDDO

 summHartree=summHartree*LEVEL%dVol
 summHartree=summHartree/Vol

 call MPI_ALLREDUCE(summHartree,summHRecv,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM,IERR)
 call MPI_ALLREDUCE(summ,summRecv,1,MPI_DOUBLE_PRECISION,MPI_SUM,COMM,IERR)


IF(MYID.eq.0) write(*,*) "summHreCv=",summHRecV,LEVEL%dVol,Vol

IF(FLAG1.eq.1) THEN
DO K=LEVEL%Z1(CpuN),LEVEL%Z2(CpuN)
  DO J=LEVEL%Y1(CpuN),LEVEL%Y2(CpuN)
	DO I=LEVEL%X1(CpuN),LEVEL%X2(CpuN)
	LEVEL%PHI(I,J,K)=LEVEL%PHI(I,J,K)*LEVEL%eps(2,I,J,K)	
	LEVEL%PHI(I,J,K)=LEVEL%PHI(I,J,K)-summHRecV
	LEVEL%PHI(I,J,K)=LEVEL%PHI(I,J,K)/LEVEL%eps(2,I,J,K)	
 	ENDDO
  ENDDO
ENDDO

ENDIF


IF(MYID.eq.0) then
	IF(flag.eq.1) WRITE(*,*) "Error",ABS(summRecV-summ_prev)/summRecV*100,ITER
ENDIF

ENDDO


 call UPDATE(MYID,LEVEL,COMM,0)
 call UpdateBoundary(MYID,Level,COMM,0)


END SUBROUTINE GaussSeidelSmoothSOR






SUBROUTINE GetNLevels(NDIM,nlevels)
IMPLICIT NONE
INTEGER I,J,K,NDIM,nlevels,NDIM1
nlevels=1


NDIM1=NDIM
 DO WHILE(1.GT.0)
	IF(MOD(NDIM1,2).EQ.0.AND.NDIM1/2.GE.10) THEN
!	IF(MOD(NDIM,2).EQ.0.AND.NDIM.GE.8) THEN
!	IF(MOD(NDIM,2).EQ.0) THEN
	 nlevels=nlevels+1
	 NDIM1=NDIM1/2
	 ELSE
	 goto 122
	ENDIF

 ENDDO
122 continue
END SUBROUTINE GetNLevels


SUBROUTINE ReadData(vol,LEVEL)
USE levelType
IMPLICIT NONE
INTEGER I,J,K
REAL*8 summRho
REAL*8 Xtemp,Ytemp,Ztemp,A2B,dhX,dhY,dhZ,Xmin,Ymin,Zmin,Xmax,Ymax,Zmax,dVol,Pi,vol
TYPE (LevelStruct)::Level



open (UNIT=82, FILE='sim_rho_r.dat', STATUS='OLD')
!open (UNIT=102, FILE='hartree.dat', STATUS='OLD')



write(*,*) LEVEL%ndim
!pause
A2B=1./0.529177249
Pi=4.*Atan(1.d0)
summRho=0.
DO K=1,LEVEL%NDIM
  DO J=1,LEVEL%NDIM
    DO I=1,LEVEL%NDIM
	READ(82,*) LEVEL%X(I),LEVEL%Y(J),LEVEL%Z(K),LEVEL%RHS(i,j,k)
	summRho=summRho+LEVEL%RHS(i,j,k)
!	READ(102,*) Xtemp,Ytemp,Ztemp,vHartree(I,J,K)
    ENDDO
  ENDDO
ENDDO




LEVEL%X=LEVEL%X*A2B
LEVEL%Y=LEVEL%Y*A2B
LEVEL%Z=LEVEL%Z*A2B
write(*,*) LEVEL%X(1),LEVEL%X(LEVEL%NDIM-1)
write(*,*) LEVEL%X(1),LEVEL%X(LEVEL%NDIM-1)

dhX=ABS(LEVEL%X(2)-LEVEL%X(1))
dhY=ABS(LEVEL%Y(2)-LEVEL%Y(1))
dhZ=ABS(LEVEL%Z(2)-LEVEL%Z(1))


Xmin=MINVAL(LEVEL%X)
Ymin=MINVAL(LEVEL%Y)
Zmin=MINVAL(LEVEL%Z)
Xmax=MAXVAL(LEVEL%X)
Ymax=MAXVAL(LEVEL%Y)
Zmax=MAXVAL(LEVEL%Z)


LEVEL%dvol=dhX*dhY*dhZ
vol=(Xmax-Xmin+dhX)*(Ymax-Ymin+dhY)*(Zmax-Zmin+dhZ)



write(*,*) "vol,dvol  =",vol,LEVEL%dVol,dhX
write(*,*) "SumRho  before   =",summRho*dVol

write(*,*) "X(0), X(NDIM)   =",LEVEL%X(1)/A2B, LEVEL%X(LEVEL%NDIM)/A2B

LEVEL%RHS=LEVEL%RHS-summRho/Vol*LEVEL%dVol


LEVEL%RHS=LEVEL%RHS*4.*Pi*dhX*dhX

summRho=0.



DO K=1,LEVEL%NDIM
  DO J=1,LEVEL%NDIM
    DO I=1,LEVEL%NDIM
	summRho=summRho+LEVEL%RHS(i,j,k)*LEVEL%dVol
    ENDDO
  ENDDO
ENDDO


write(*,*) "SumRho    after =",summRho*LEVEL%dVol
write(*,*) Xmin,Ymin,Zmin
write(*,*) Xmax,Ymax,Zmax

 close (82)
 close (102)



END SUBROUTINE ReadData




!SUBROUTINE    propagateXYZ(Xb,Yb,Zb,X,Y,Z,NDIMb,NDIM)	
SUBROUTINE    propagateXYZ(LEVEL1,LEVEL)	
USE levelType
IMPLICIT NONE
INTEGER I,Ism
TYPE (LevelStruct)::LEVEL,LEVEL1


DO I=1,LEVEL%NDIM
	LEVEL%X(I)=LEVEL1%X(2*I-1)
	LEVEL%Y(I)=LEVEL1%Y(2*I-1)
	LEVEL%Z(I)=LEVEL1%Z(2*I-1)
ENDDO


END SUBROUTINE propagateXYZ


SUBROUTINE    fillEps(LEVEL,Xoffset,Yoffset,Zoffset)	
Use LevelType
IMPLICIT NONE
INTEGER I,J,K,M,NDIM,Ip,Im,Jp,Jm,Kp,Km,P
REAL*8 Xoffset,Yoffset,Zoffset
REAL*8 R(6),dhX,dhY,dhZ,b,v,a1,a2,a3,a,d12,d32,temp1,temp2,temp,F
TYPE (LevelStruct)::LEVEL

NDIM=LEVEL%NDIM
!dhX=ABS(LEVEL%X(2)-LEVEL%X(1))
!dhY=ABS(LEVEL%Y(2)-LEVEL%Y(1))
!dhZ=ABS(LEVEL%Z(2)-LEVEL%Z(1))
b=12.
v=0.8
a=1.
DO K=1,NDIM
  DO J=1,NDIM
	DO I=1,NDIM
	R(1)=SQRT((LEVEL%X(I)-Xoffset)**2.+(LEVEL%Y(J)-Yoffset)**2.+(LEVEL%Z(K)-Zoffset)**2.)
	ip=MOD(I,NDIM)+1
	im=I-1
        jp=MOD(J,NDIM)+1
	jm=J-1
	kp=MOD(K,NDIM)+1
	km=K-1
        IF(I.EQ.1) im=NDIM
	IF(J.EQ.1) jm=NDIM
	IF(K.EQ.1) km=NDIM

!	LEVEL%EPS(1,I,J,K)=a+b*DExp(-v*(R(1))**2.)
	LEVEL%EPS(1,I,J,K)=1.
	LEVEL%EPS(2,I,J,K)=LEVEL%EPS(1,I,J,K)**(-0.5D0)
	LEVEL%EPS(3,I,J,K)=LEVEL%EPS(1,I,J,K)**(0.5D0)
	ENDDO
  ENDDO
ENDDO


DO K=1,NDIM
  DO J=1,NDIM
	DO I=1,NDIM
	R(1)=SQRT((LEVEL%X(I)-Xoffset)**2.+(LEVEL%Y(J)-Yoffset)**2.+(LEVEL%Z(K)-Zoffset)**2.)
	ip=MOD(I,NDIM)+1
	im=I-1
        jp=MOD(J,NDIM)+1
	jm=J-1
	kp=MOD(K,NDIM)+1
	km=K-1
        IF(I.EQ.1) im=NDIM
	IF(J.EQ.1) jm=NDIM
	IF(K.EQ.1) km=NDIM
        LEVEL%EPS(4,I,J,K)=LEVEL%EPS(3,Ip,J,K)+LEVEL%EPS(3,Im,J,K)+LEVEL%EPS(3,I,Jp,K)+LEVEL%EPS(3,I,Jm,K)+LEVEL%EPS(3,I,J,Kp)+&
        & LEVEL%EPS(3,I,J,Km)-6.*LEVEL%EPS(3,I,J,K)
	LEVEL%EPS(5,I,J,K)=LEVEL%EPS(4,I,J,K)*LEVEL%EPS(2,I,J,K)
	ENDDO
  ENDDO
ENDDO

END SUBROUTINE    fillEps






SUBROUTINE INITIALIZE(MYID,NCPUS,LEVEL,Level1,levelN,nlevels)
USE LevelType
IMPLICIT NONE
INTEGER NDIM,NDIM3,NSUB,NCPUS,MYID,NDIV,CNT,I,J,K,CNT1(3),MULT,X_int,Y_int,Z_int,Nx,Ny,Nz,NSUB1,P,INCX,INCY,INCZ,MODX,MODY,MODZ
INTEGER levelN,nlevels,KK
INTEGER, DIMENSION(:), ALLOCATABLE ::NDIM_MAX,NDIM_MAX_LINE
TYPE (LevelStruct)::Level,Level1


NDIM=LEVEL%NDIM


NDIM3=NDIM*NDIM*NDIM
NSUB=NDIM3/NCPUS

NDIV=INT(NCPUS**(1./3.))+1

!NDIV=NCPUS-1

NSUB1=NCPUS
 CNT=0
121 continue
 DO I=NDIV,1,-1
	IF(MOD(NSUB1,I).EQ.0) THEN
	NSUB1=NSUB1/I
	CNT=CNT+1
	CNT1(CNT)=I
	goto 122
	ENDIF
 ENDDO

 122 continue
 IF(CNT.eq.2) THEN 
	 goto 123 
     ELSE
         goto 121	
     ENDIF

 123 continue
   CNT=CNT+1
  CNT1(CNT)=NSUB1
 MULT=1
! IF(MYID.eq.0) WRITE(*,*) "######################################"



 call SSORT(CNT1,3)



 IF(MYID.eq.0) THEN
  WRITE(*,*) "######################################"
  WRITE(*,*) CNT1(1),CNT1(2),CNT1(3),NDIM
 ENDIF

 DO I=1,CNT
  MULT=MULT*CNT1(I)
 ENDDO


 CNT=0
DO K=1,CNT1(3)
 DO J=1,CNT1(2)
  DO I=1,CNT1(1)	
   CNT=CNT+1 
   LEVEL%X1(CNT)=NDIM/CNT1(1)*(I-1)+1
   LEVEL%X2(CNT)=NDIM/CNT1(1)*I
   LEVEL%Y1(CNT)=NDIM/CNT1(2)*(J-1)+1
   LEVEL%Y2(CNT)=NDIM/CNT1(2)*J
   LEVEL%Z1(CNT)=NDIM/CNT1(3)*(K-1)+1
   LEVEL%Z2(CNT)=NDIM/CNT1(3)*K

   IF(I.eq.CNT1(1)) THEN
    LEVEL%X2(CNT)=NDIM
   ENDIF

   IF(J.eq.CNT1(2)) THEN
    LEVEL%Y2(CNT)=NDIM
   ENDIF

   IF(K.eq.CNT1(3)) THEN
    LEVEL%Z2(CNT)=NDIM
   ENDIF

 IF(levelN.le.nlevels-1) THEN
  LEVEL%X1(CNT)=(LEVEL1%X1(CNT)+1)/2
  LEVEL%Y1(CNT)=(LEVEL1%Y1(CNT)+1)/2
  LEVEL%Z1(CNT)=(LEVEL1%Z1(CNT)+1)/2
  LEVEL%X2(CNT)=(LEVEL1%X2(CNT)+1)/2
  LEVEL%Y2(CNT)=(LEVEL1%Y2(CNT)+1)/2
  LEVEL%Z2(CNT)=(LEVEL1%Z2(CNT)+1)/2

   DO KK=1,CNT-1
	 IF(CNT.gt.1.and.(LEVEL%X1(CNT).eq.LEVEL%X2(KK))) THEN
!	   LEVEL%X1(CNT)=LEVEL%X1(CNT)+1
	 ENDIF

	 IF(CNT.gt.1.and.(LEVEL%Y1(CNT).eq.LEVEL%Y2(KK))) THEN
!	   LEVEL%Y1(CNT)= LEVEL%Y1(CNT)+1
	 ENDIF
	 IF(CNT.gt.1.and.(LEVEL%Z1(CNT).eq.LEVEL%Z2(KK))) THEN
!	   LEVEL%Z1(CNT)= LEVEL%Z1(CNT)+1
	 ENDIF
   ENDDO

   IF(I.eq.CNT1(1)) THEN
    LEVEL%X2(CNT)=NDIM
   ENDIF

   IF(J.eq.CNT1(2)) THEN
    LEVEL%Y2(CNT)=NDIM
   ENDIF

   IF(K.eq.CNT1(3)) THEN
    LEVEL%Z2(CNT)=NDIM
   ENDIF

 ENDIF


 IF(MYID.EQ.0) THEN
 WRITE(*,'(A1,I3,A1,I3,A1,A5,A1,I3,A1,I3,A1,A5,A1,I3,A1,I3,A1,I3)') "(",LEVEL%X1(CNT),",",LEVEL%X2(CNT),")","     ","(",LEVEL%Y1(CNT),",",LEVEL%Y2(CNT)&
,")","     ","(",LEVEL%Z1(CNT),",",LEVEL%Z2(CNT),")",CNT
ENDIF


  ENDDO
 ENDDO
ENDDO

LEVEL%Nb=0
DO I=1,NCPUS
  X_int=LEVEL%X1(I)-1
   IF(LEVEL%X1(I).EQ.1) X_int=NDIM
   DO J=1,NCPUS
    IF((X_int.GE.LEVEL%X1(J)).AND.(X_int.LE.LEVEL%X2(J)).AND.(LEVEL%Y1(I).EQ.LEVEL%Y1(J)).AND.(LEVEL%Y2(I).EQ.LEVEL%Y2(J))&
    &.AND.(LEVEL%Z1(I).EQ.LEVEL%Z1(J)).AND.(LEVEL%Z2(I).EQ.LEVEL%Z2(J)).AND.(I.NE.J)) THEN 
    LEVEL%NB(I)=LEVEL%NB(I)+1
    LEVEL%RecV(I,LEVEL%NB(I))=J
    LEVEL%X1Plane(I,LEVEL%NB(I))=X_int
    LEVEL%X2Plane(I,LEVEL%NB(I))=X_int
    LEVEL%Y1Plane(I,LEVEL%NB(I))=LEVEL%Y1(I)
    LEVEL%Y2Plane(I,LEVEL%NB(I))=LEVEL%Y2(I)
    LEVEL%Z1Plane(I,LEVEL%NB(I))=LEVEL%Z1(I)
    LEVEL%Z2Plane(I,LEVEL%NB(I))=LEVEL%Z2(I)
    ENDIF
   ENDDO	
   X_int=LEVEL%X2(I)+1
   IF(LEVEL%X2(I).EQ.NDIM) X_int=1
   DO J=1,NCPUS
    IF((X_int.GE.LEVEL%X1(J)).AND.(X_int.LE.LEVEL%X2(J)).AND.(LEVEL%Y1(I).EQ.LEVEL%Y1(J)).AND.(LEVEL%Y2(I).EQ.LEVEL%Y2(J))&
    &.AND.(LEVEL%Z1(I).EQ.LEVEL%Z1(J)).AND.(LEVEL%Z2(I).EQ.LEVEL%Z2(J)).AND.(I.NE.J)) THEN 
    LEVEL%NB(I)=LEVEL%NB(I)+1
    LEVEL%RecV(I,LEVEL%NB(I))=J
    LEVEL%X1Plane(I,LEVEL%NB(I))=X_int
    LEVEL%X2Plane(I,LEVEL%NB(I))=X_int
    LEVEL%Y1Plane(I,LEVEL%NB(I))=LEVEL%Y1(I)
    LEVEL%Y2Plane(I,LEVEL%NB(I))=LEVEL%Y2(I)
    LEVEL%Z1Plane(I,LEVEL%NB(I))=LEVEL%Z1(I)
    LEVEL%Z2Plane(I,LEVEL%NB(I))=LEVEL%Z2(I)

    ENDIF
   ENDDO	

  Y_int=LEVEL%Y1(I)-1
   IF(LEVEL%Y1(I).EQ.1) Y_int=NDIM
   DO J=1,NCPUS
    IF((Y_int.GE.LEVEL%Y1(J)).AND.(Y_int.LE.LEVEL%Y2(J)).AND.(LEVEL%X1(I).EQ.LEVEL%X1(J)).AND.(LEVEL%X2(I).EQ.LEVEL%X2(J))&
   & .AND.(LEVEL%Z1(I).EQ.LEVEL%Z1(J)).AND.(LEVEL%Z2(I).EQ.LEVEL%Z2(J)).AND.(I.NE.J)) THEN 
    LEVEL%NB(I)=LEVEL%NB(I)+1
    LEVEL%RecV(I,LEVEL%NB(I))=J
    LEVEL%X1Plane(I,LEVEL%NB(I))=LEVEL%X1(I)
    LEVEL%X2Plane(I,LEVEL%NB(I))=LEVEL%X2(I)
    LEVEL%Y1Plane(I,LEVEL%NB(I))=Y_int
    LEVEL%Y2Plane(I,LEVEL%NB(I))=Y_int
    LEVEL%Z1Plane(I,LEVEL%NB(I))=LEVEL%Z1(I)
    LEVEL%Z2Plane(I,LEVEL%NB(I))=LEVEL%Z2(I)
    ENDIF
   ENDDO	
   Y_int=LEVEL%Y2(I)+1
   IF(LEVEL%Y2(I).EQ.NDIM) Y_int=1
   DO J=1,NCPUS
    IF((Y_int.GE.LEVEL%Y1(J)).AND.(Y_int.LE.LEVEL%Y2(J)).AND.(LEVEL%X1(I).EQ.LEVEL%X1(J)).AND.(LEVEL%X2(I).EQ.LEVEL%X2(J)).AND.(LEVEL%Z1(I).EQ.LEVEL%Z1(J))&
   &.AND.(LEVEL%Z2(I).EQ.LEVEL%Z2(J)).AND.(I.NE.J)) THEN 
    LEVEL%NB(I)=LEVEL%NB(I)+1
    LEVEL%RecV(I,LEVEL%NB(I))=J
    LEVEL%X1Plane(I,LEVEL%NB(I))=LEVEL%X1(I)
    LEVEL%X2Plane(I,LEVEL%NB(I))=LEVEL%X2(I)
    LEVEL%Y1Plane(I,LEVEL%NB(I))=Y_int
    LEVEL%Y2Plane(I,LEVEL%NB(I))=Y_int
    LEVEL%Z1Plane(I,LEVEL%NB(I))=LEVEL%Z1(I)
    LEVEL%Z2Plane(I,LEVEL%NB(I))=LEVEL%Z2(I)
    ENDIF
   ENDDO	


  Z_int=LEVEL%Z1(I)-1
   IF(LEVEL%Z1(I).EQ.1) Z_int=NDIM
   DO J=1,NCPUS
    IF((Z_int.GE.LEVEL%Z1(J)).AND.(Z_int.LE.LEVEL%Z2(J)).AND.(LEVEL%Y1(I).EQ.LEVEL%Y1(J)).AND.(LEVEL%Y2(I).EQ.LEVEL%Y2(J)).AND.(LEVEL%X1(I).EQ.LEVEL%X1(J))&
    &.AND.(LEVEL%X2(I).EQ.LEVEL%X2(J)).AND.(I.NE.J)) THEN 
    LEVEL%NB(I)=LEVEL%NB(I)+1
    LEVEL%RecV(I,LEVEL%NB(I))=J
    LEVEL%X1Plane(I,LEVEL%NB(I))=LEVEL%X1(I)
    LEVEL%X2Plane(I,LEVEL%NB(I))=LEVEL%X2(I)
    LEVEL%Y1Plane(I,LEVEL%NB(I))=LEVEL%Y1(I)
    LEVEL%Y2Plane(I,LEVEL%NB(I))=LEVEL%Y2(I)
    LEVEL%Z1Plane(I,LEVEL%NB(I))=Z_int
    LEVEL%Z2Plane(I,LEVEL%NB(I))=Z_int
    ENDIF
   ENDDO	

   Z_int=LEVEL%Z2(I)+1
   IF(LEVEL%Z2(I).EQ.NDIM) Z_int=1
   DO J=1,NCPUS
    IF((Z_int.GE.LEVEL%Z1(J)).AND.(Z_int.LE.LEVEL%Z2(J)).AND.(LEVEL%Y1(I).EQ.LEVEL%Y1(J)).AND.(LEVEL%Y2(I).EQ.LEVEL%Y2(J)).AND.(LEVEL%X1(I).EQ.LEVEL%X1(J))&
    &.AND.(LEVEL%X2(I).EQ.LEVEL%X2(J)).AND.(I.NE.J)) THEN 
    LEVEL%NB(I)=LEVEL%NB(I)+1
    LEVEL%RecV(I,LEVEL%NB(I))=J
    LEVEL%X1Plane(I,LEVEL%NB(I))=LEVEL%X1(I)
    LEVEL%X2Plane(I,LEVEL%NB(I))=LEVEL%X2(I)
    LEVEL%Y1Plane(I,LEVEL%NB(I))=LEVEL%Y1(I)
    LEVEL%Y2Plane(I,LEVEL%NB(I))=LEVEL%Y2(I)
    LEVEL%Z1Plane(I,LEVEL%NB(I))=Z_int
    LEVEL%Z2Plane(I,LEVEL%NB(I))=Z_int
    ENDIF
   ENDDO	

ENDDO



 DO I=1,NCPUS
  DO K=1,NCPUS
   DO J=1,LEVEL%NB(K)
    IF(LEVEL%RecV(K,J).eq.I) THEN 
	 LEVEL%SenD(I,J)=K
 	 LEVEL%PlanesSend(I,J)=J
    ENDIF
   ENDDO
  ENDDO
 ENDDO



DO I=1,NCPUS
   DO J=1,LEVEL%NB(I)
    IF((J.eq.1).or.(J.eq.2)) then
	DO P=1,4
	 LEVEL%X1Line(P,I,J)=LEVEL%X1Plane(I,J)
 	 LEVEL%X2Line(P,I,J)=LEVEL%X2Plane(I,J)
	IF(P.eq.1) then
	 LEVEL%Y1Line(P,I,J)=LEVEL%Y1Plane(I,J)-1
	 LEVEL%Y2Line(P,I,J)=LEVEL%Y1Plane(I,J)-1
	 LEVEL%Z1Line(P,I,J)=LEVEL%Z1Plane(I,J)
	 LEVEL%Z2Line(P,I,J)=LEVEL%Z2Plane(I,J)
	   IF(LEVEL%Y1Plane(I,J).eq.1) then 
		 LEVEL%Y1Line(P,I,J)=NDIM
		 LEVEL%Y2Line(P,I,J)=NDIM
	   endif
 	endif

	IF(P.eq.2) then
	 LEVEL%Y1Line(P,I,J)=LEVEL%Y2Plane(I,J)+1
	 LEVEL%Y2Line(P,I,J)=LEVEL%Y2Plane(I,J)+1
	 LEVEL%Z1Line(P,I,J)=LEVEL%Z1Plane(I,J)
	 LEVEL%Z2Line(P,I,J)=LEVEL%Z2Plane(I,J)
	   IF(LEVEL%Y2Plane(I,J).eq.NDIM) then 
		 LEVEL%Y1Line(P,I,J)=1
		 LEVEL%Y2Line(P,I,J)=1
	   endif
 	endif
   
	IF(P.eq.3) then
	 LEVEL%Y1Line(P,I,J)=LEVEL%Y1Plane(I,J)
	 LEVEL%Y2Line(P,I,J)=LEVEL%Y2Plane(I,J)
	 LEVEL%Z1Line(P,I,J)=LEVEL%Z1Plane(I,J)-1
	 LEVEL%Z2Line(P,I,J)=LEVEL%Z1Plane(I,J)-1
	   IF(LEVEL%Z1Plane(I,J).eq.1) then 
		 LEVEL%Z1Line(P,I,J)=NDIM
		 LEVEL%Z2Line(P,I,J)=NDIM
	   endif
 	endif
   
	IF(P.eq.4) then
	 LEVEL%Y1Line(P,I,J)=LEVEL%Y1Plane(I,J)
	 LEVEL%Y2Line(P,I,J)=LEVEL%Y2Plane(I,J)
	 LEVEL%Z1Line(P,I,J)=LEVEL%Z2Plane(I,J)+1
	 LEVEL%Z2Line(P,I,J)=LEVEL%Z2Plane(I,J)+1
	   IF(LEVEL%Z2Plane(I,J).eq.NDIM) then 
		 LEVEL%Z1Line(P,I,J)=1
		 LEVEL%Z2Line(P,I,J)=1
	   endif
 	endif

       ENDDO	
      endif	




    IF((J.eq.3).or.J.eq.4) then
	DO P=1,4
	 LEVEL%Y1Line(P,I,J)=LEVEL%Y1Plane(I,J)
 	 LEVEL%Y2Line(P,I,J)=LEVEL%Y2Plane(I,J)
	IF(P.eq.1) then
	 LEVEL%X1Line(P,I,J)=LEVEL%X1Plane(I,J)-1
	 LEVEL%X2Line(P,I,J)=LEVEL%X1Plane(I,J)-1
	 LEVEL%Z1Line(P,I,J)=LEVEL%Z1Plane(I,J)
	 LEVEL%Z2Line(P,I,J)=LEVEL%Z2Plane(I,J)
	   IF(LEVEL%X1Plane(I,J).eq.1) then 
		 LEVEL%X1Line(P,I,J)=NDIM
		 LEVEL%X2Line(P,I,J)=NDIM
	   endif
 	endif

	IF(P.eq.2) then
	 LEVEL%X1Line(P,I,J)=LEVEL%X2Plane(I,J)+1
	 LEVEL%X2Line(P,I,J)=LEVEL%X2Plane(I,J)+1
	 LEVEL%Z1Line(P,I,J)=LEVEL%Z1Plane(I,J)
	 LEVEL%Z2Line(P,I,J)=LEVEL%Z2Plane(I,J)
	   IF(LEVEL%X2Plane(I,J).eq.NDIM) then 
		 LEVEL%X1Line(P,I,J)=1
		 LEVEL%X2Line(P,I,J)=1
	   endif
 	endif
   
	IF(P.eq.3) then
	 LEVEL%X1Line(P,I,J)=LEVEL%X1Plane(I,J)
	 LEVEL%X2Line(P,I,J)=LEVEL%X2Plane(I,J)
	 LEVEL%Z1Line(P,I,J)=LEVEL%Z1Plane(I,J)-1
	 LEVEL%Z2Line(P,I,J)=LEVEL%Z1Plane(I,J)-1
	   IF(LEVEL%Z1Plane(I,J).eq.1) then 
		 LEVEL%Z1Line(P,I,J)=NDIM
		 LEVEL%Z2Line(P,I,J)=NDIM
	   endif
 	endif
   
	IF(P.eq.4) then
	 LEVEL%X1Line(P,I,J)=LEVEL%X1Plane(I,J)
	 LEVEL%X2Line(P,I,J)=LEVEL%X2Plane(I,J)
	 LEVEL%Z1Line(P,I,J)=LEVEL%Z2Plane(I,J)+1
	 LEVEL%Z2Line(P,I,J)=LEVEL%Z2Plane(I,J)+1
	   IF(LEVEL%Z2Plane(I,J).eq.NDIM) then 
		 LEVEL%Z1Line(P,I,J)=1
		 LEVEL%Z2Line(P,I,J)=1
	   endif
 	endif

       ENDDO	
    endif	


    IF((J.eq.5).or.J.eq.6) then

	DO P=1,4
	 LEVEL%Z1Line(P,I,J)=LEVEL%Z1Plane(I,J)
 	 LEVEL%Z2Line(P,I,J)=LEVEL%Z2Plane(I,J)
	IF(P.eq.1) then
	 LEVEL%X1Line(P,I,J)=LEVEL%X1Plane(I,J)-1
	 LEVEL%X2Line(P,I,J)=LEVEL%X1Plane(I,J)-1
	 LEVEL%Y1Line(P,I,J)=LEVEL%Y1Plane(I,J)
	 LEVEL%Y2Line(P,I,J)=LEVEL%Y2Plane(I,J)
	   IF(LEVEL%X1Plane(I,J).eq.1) then 
		 LEVEL%X1Line(P,I,J)=NDIM
		 LEVEL%X2Line(P,I,J)=NDIM
	   endif
 	endif

	IF(P.eq.2) then
	 LEVEL%X1Line(P,I,J)=LEVEL%X2Plane(I,J)+1
	 LEVEL%X2Line(P,I,J)=LEVEL%X2Plane(I,J)+1
	 LEVEL%Y1Line(P,I,J)=LEVEL%Y1Plane(I,J)
	 LEVEL%Y2Line(P,I,J)=LEVEL%Y2Plane(I,J)
	   IF(LEVEL%X2Plane(I,J).eq.NDIM) then 
		 LEVEL%X1Line(P,I,J)=1
		 LEVEL%X2Line(P,I,J)=1
	   endif
 	endif
   
	IF(P.eq.3) then
	 LEVEL%X1Line(P,I,J)=LEVEL%X1Plane(I,J)
	 LEVEL%X2Line(P,I,J)=LEVEL%X2Plane(I,J)
	 LEVEL%Y1Line(P,I,J)=LEVEL%Y1Plane(I,J)-1
	 LEVEL%Y2Line(P,I,J)=LEVEL%Y1Plane(I,J)-1
	   IF(LEVEL%Y1Plane(I,J).eq.1) then 
		 LEVEL%Y1Line(P,I,J)=NDIM
		 LEVEL%Y2Line(P,I,J)=NDIM
	   endif
 	endif
   
	IF(P.eq.4) then
	 LEVEL%X1Line(P,I,J)=LEVEL%X1Plane(I,J)
	 LEVEL%X2Line(P,I,J)=LEVEL%X2Plane(I,J)
	 LEVEL%Y1Line(P,I,J)=LEVEL%Y2Plane(I,J)+1
	 LEVEL%Y2Line(P,I,J)=LEVEL%Y2Plane(I,J)+1
	   IF(LEVEL%Y2Plane(I,J).eq.NDIM) then 
		 LEVEL%Y1Line(P,I,J)=1
		 LEVEL%Y2Line(P,I,J)=1
	   endif
 	endif

       ENDDO	
    endif	
! IF(MYID.eq.0.and.NDIM.eq.120.and.J.eq.5) write(*,*) "Z1%%,Z2= ",LEVEL%Z1Line(3,2,5),LEVEL%Z2Line(3,2,5),P,I,J


   ENDDO	
ENDDO


DO I=1,NCPUS
   DO J=1,LEVEL%NB(I)
     IF((J.eq.1).or.(J.eq.2)) then
       DO P=1,4	
	 LEVEL%XPoint(P,I,J)=LEVEL%X1Line(P,I,J)
	IF(P.eq.1) then
	 LEVEL%YPoint(P,I,J)=LEVEL%Y1(I)-1
	 LEVEL%ZPoint(P,I,J)=LEVEL%Z1(I)-1
	   IF(LEVEL%Y1(I).eq.1) then 
		LEVEL%YPoint(P,I,J)=NDIM
	   endif
	   IF(LEVEL%Z1(I).eq.1) then 
		LEVEL%ZPoint(P,I,J)=NDIM
	   endif
        ENDIF

	IF(P.eq.2) then
	 LEVEL%YPoint(P,I,J)=LEVEL%Y2(I)+1
	 LEVEL%ZPoint(P,I,J)=LEVEL%Z1(I)-1
	   IF(LEVEL%Y2(I).eq.NDIM) then 
		LEVEL%YPoint(P,I,J)=1
	   endif
	   IF(LEVEL%Z1(I).eq.1) then 
		LEVEL%ZPoint(P,I,J)=NDIM
	   endif
        ENDIF


	IF(P.eq.3) then
	 LEVEL%YPoint(P,I,J)=LEVEL%Y1(I)-1
	 LEVEL%ZPoint(P,I,J)=LEVEL%Z2(I)+1
	   IF(LEVEL%Y1(I).eq.1) then 
		LEVEL%YPoint(P,I,J)=NDIM
	   endif
	   IF(LEVEL%Z2(I).eq.NDIM) then 
		LEVEL%ZPoint(P,I,J)=1
	   endif
        ENDIF

	IF(P.eq.4) then
	 LEVEL%YPoint(P,I,J)=LEVEL%Y2(I)+1
	 LEVEL%ZPoint(P,I,J)=LEVEL%Z2(I)+1
	   IF(LEVEL%Y2(I).eq.NDIM) then 
		LEVEL%YPoint(P,I,J)=1
	   endif
	   IF(LEVEL%Z2(I).eq.NDIM) then 
		LEVEL%ZPoint(P,I,J)=1
	   endif
        ENDIF
     ENDDO
    ENDIF


     IF((J.eq.3).or.(J.eq.4)) then
       DO P=1,4	
	 LEVEL%YPoint(P,I,J)=LEVEL%Y1Line(P,I,J)
	IF(P.eq.1) then
	 LEVEL%XPoint(P,I,J)=LEVEL%X1(I)-1
	 LEVEL%ZPoint(P,I,J)=LEVEL%Z1(I)-1
	   IF(LEVEL%X1(I).eq.1) then 
		LEVEL%XPoint(P,I,J)=NDIM
	   endif
	   IF(LEVEL%Z1(I).eq.1) then 
		LEVEL%ZPoint(P,I,J)=NDIM
	   endif
        ENDIF

	IF(P.eq.2) then
	 LEVEL%XPoint(P,I,J)=LEVEL%X2(I)+1
	 LEVEL%ZPoint(P,I,J)=LEVEL%Z1(I)-1
	   IF(LEVEL%X2(I).eq.NDIM) then 
		LEVEL%XPoint(P,I,J)=1
	   endif
	   IF(LEVEL%Z1(I).eq.1) then 
		LEVEL%ZPoint(P,I,J)=NDIM
	   endif
        ENDIF


	IF(P.eq.3) then
	 LEVEL%XPoint(P,I,J)=LEVEL%X1(I)-1
	 LEVEL%ZPoint(P,I,J)=LEVEL%Z2(I)+1
	   IF(LEVEL%X1(I).eq.1) then 
		LEVEL%XPoint(P,I,J)=NDIM
	   endif
	   IF(LEVEL%Z2(I).eq.NDIM) then 
		LEVEL%ZPoint(P,I,J)=1
	   endif
        ENDIF

	IF(P.eq.4) then
	 LEVEL%XPoint(P,I,J)=LEVEL%X2(I)+1
	 LEVEL%ZPoint(P,I,J)=LEVEL%Z2(I)+1
	   IF(LEVEL%X2(I).eq.NDIM) then 
		LEVEL%XPoint(P,I,J)=1
	   endif
	   IF(LEVEL%Z2(I).eq.NDIM) then 
		LEVEL%ZPoint(P,I,J)=1
	   endif
        ENDIF
     ENDDO
    ENDIF


     IF((J.eq.5).or.(J.eq.6)) then
       DO P=1,4	
	 LEVEL%ZPoint(P,I,J)=LEVEL%Z1Line(P,I,J)
	IF(P.eq.1) then
	 LEVEL%XPoint(P,I,J)=LEVEL%X1(I)-1
	 LEVEL%YPoint(P,I,J)=LEVEL%Y1(I)-1
	   IF(LEVEL%Y1(I).eq.1) then 
		LEVEL%YPoint(P,I,J)=NDIM
	   endif
	   IF(LEVEL%X1(I).eq.1) then 
		LEVEL%XPoint(P,I,J)=NDIM
	   endif
        ENDIF

	IF(P.eq.2) then
	 LEVEL%XPoint(P,I,J)=LEVEL%X1(I)-1
	 LEVEL%YPoint(P,I,J)=LEVEL%Y2(I)+1

	   IF(LEVEL%Y2(I).eq.NDIM) then 
		LEVEL%YPoint(P,I,J)=1
	   endif
	   IF(LEVEL%X1(I).eq.1) then 
		LEVEL%XPoint(P,I,J)=NDIM
	   endif
        ENDIF


	IF(P.eq.3) then
	 LEVEL%XPoint(P,I,J)=LEVEL%X2(I)+1
	 LEVEL%YPoint(P,I,J)=LEVEL%Y1(I)-1
	   IF(LEVEL%Y1(I).eq.1) then 
		LEVEL%YPoint(P,I,J)=NDIM
	   endif
	   IF(LEVEL%X2(I).eq.NDIM) then 
		LEVEL%XPoint(P,I,J)=1
	   endif
        ENDIF

	IF(P.eq.4) then
	 LEVEL%XPoint(P,I,J)=LEVEL%X2(I)+1
	 LEVEL%YPoint(P,I,J)=LEVEL%Y2(I)+1
	   IF(LEVEL%Y2(I).eq.NDIM) then 
		LEVEL%YPoint(P,I,J)=1
	   endif
	   IF(LEVEL%X2(I).eq.NDIM) then 
		LEVEL%XPoint(P,I,J)=1
	   endif
        ENDIF
     ENDDO
    ENDIF

 ENDDO
ENDDO



!pause
IF(MYID.eq.0) THEN
DO I=1,NCPUS
write(*,*) "## LINES  LINES  LINES  LINES  LINES   LINES ##"

 DO J=1,LEVEL%NB(I)
 DO P=1,4
IF(MYID.eq.0) then
WRITE(*,'(A1,I3,A1,I3,A1,A5,A1,I3,A1,I3,A1,A5,A1,I3,A1,I3,A1,I3)') "(",LEVEL%X1Line(P,I,J),",",LEVEL%X2Line(P,I,J),")","     ","(",&
     & LEVEL%Y1Line(P,I,J),",",LEVEL%Y2Line(P,I,J),")","      ","(",LEVEL%Z1Line(P,I,J),",",LEVEL%Z2Line(P,I,J),")",I
ENDIF

 ENDDO
WRITE(*,*) "************"
 ENDDO
ENDDO

ENDIF


DO I=1,NCPUS
 DO K=1,NCPUS
  DO J=1,LEVEL%NB(I)
   DO P=1,4
     IF((LEVEL%X1Line(P,I,J).GE.LEVEL%X1(K)).AND.(LEVEL%X2Line(P,I,J).LE.LEVEL%X2(K)).AND.(LEVEL%Y1Line(P,I,J).GE.LEVEL%Y1(K))&
      &.AND.(LEVEL%Y2Line(P,I,J).LE.LEVEL%Y2(K)).AND.(LEVEL%Z1Line(P,I,J).GE.LEVEL%Z1(K)).AND.(LEVEL%Z2Line(P,I,J).LE.LEVEL%Z2(K)).AND.(I.NE.K)) then
	LEVEL%RecV1(I,J,P)=K
    endif
   ENDDO
  ENDDO
 ENDDO
ENDDO





DO I=1,NCPUS
 DO K=1,NCPUS
  DO J=1,LEVEL%NB(I)
   DO P=1,4
     IF((LEVEL%XPoint(P,I,J).GE.LEVEL%X1(K)).AND.(LEVEL%XPoint(P,I,J).LE.LEVEL%X2(K)).AND.(LEVEL%YPoint(P,I,J).GE.LEVEL%Y1(K))&
      &.AND.(LEVEL%YPoint(P,I,J).LE.LEVEL%Y2(K)).AND.(LEVEL%ZPoint(P,I,J).GE.LEVEL%Z1(K)).AND.(LEVEL%ZPoint(P,I,J).LE.LEVEL%Z2(K)).AND.(I.NE.K)) then
	LEVEL%RecV2(I,J,P)=K
    endif
   ENDDO
  ENDDO
 ENDDO
ENDDO







LEVEL%SenD1=0
DO I=1,NCPUS
 DO K=1,NCPUS
  DO J=1,LEVEL%NB(I)
   DO P=1,4   	
    IF(LEVEL%RecV1(K,J,P).eq.I) THEN 
	LEVEL%SenD1(I,J,P)=K
	LEVEL%LinesSend(I,J,P)=P
   ENDIF
   ENDDO
  ENDDO
 ENDDO
ENDDO



LEVEL%SenD2=0
DO I=1,NCPUS
 DO K=1,NCPUS
  DO J=1,LEVEL%NB(I)
   DO P=1,4   	
    IF(LEVEL%RecV2(K,J,P).eq.I) THEN 
	LEVEL%SenD2(I,J,P)=K
	LEVEL%PointsSend(I,J,P)=P
   ENDIF
   ENDDO
  ENDDO
 ENDDO
ENDDO




 Allocate(NDIM_MAX(MAXVAL(LEVEL%NB)*NCPUS))
 NDIM_MAX=0
 CNT=0
 DO I=1,NCPUS
  DO J=1,LEVEL%NB(I)
   Nx=INT(ABS(LEVEL%X1Plane(I,J)-LEVEL%X2Plane(I,J)))+1
   Ny=INT(ABS(LEVEL%Y1Plane(I,J)-LEVEL%Y2Plane(I,J)))+1
   Nz=INT(ABS(LEVEL%Z1Plane(I,J)-LEVEL%Z2Plane(I,J)))+1
   CNT=CNT+1
   NDIM_MAX(CNT)=Nx*Ny*Nz
  ENDDO
 ENDDO

 LEVEL%PlaneMaxDim=MAXVAL(NDIM_MAX)
 LEVEL%RECV=LEVEL%RECV-1
 LEVEL%SEND=LEVEL%SEND-1

 Allocate(NDIM_MAX_LINE(MAXVAL(LEVEL%NB)*4*NCPUS))
 NDIM_MAX_LINE=0
 CNT=0
 DO I=1,NCPUS
  DO J=1,LEVEL%NB(I)
    DO P=1,4   	
    Nx=INT(ABS(LEVEL%X1Line(P,I,J)-LEVEL%X2Line(P,I,J)))+1
    Ny=INT(ABS(LEVEL%Y1Line(P,I,J)-LEVEL%Y2Line(P,I,J)))+1
    Nz=INT(ABS(LEVEL%Z1Line(P,I,J)-LEVEL%Z2Line(P,I,J)))+1
    CNT=CNT+1
    NDIM_MAX_LINE(CNT)=Nx*Ny*Nz
   ENDDO
  ENDDO
 ENDDO


 LEVEL%RECV1=LEVEL%RECV1-1
 LEVEL%SEND1=LEVEL%SEND1-1
 LEVEL%LineMaxDim=MAXVAL(NDIM_MAX_LINE)


 LEVEL%RECV2=LEVEL%RECV2-1
 LEVEL%SEND2=LEVEL%SEND2-1

IF(MYID.eq.0) write(*,*) "PLANE_LINE_DIM=",LEVEL%LineMaxDIM

 DEAllocate(NDIM_MAX)
 DEAllocate(NDIM_MAX_LINE)
END SUBROUTINE INITIALIZE

SUBROUTINE GetNelements(NCPUS,LEVEL)
USE LevelType
INTEGER I,J,K,Nx,Ny,Nz
TYPE (LevelStruct)::Level

 DO I=1,NCPUS
   Nx=INT(ABS(LEVEL%X2(I)-LEVEL%X1(I)))+1
   Ny=INT(ABS(LEVEL%Y2(I)-LEVEL%Y1(I)))+1
   Nz=INT(ABS(LEVEL%Z2(I)-LEVEL%Z1(I)))+1
   LEVEL%NElements(I)=Nx*Ny*Nz
 ENDDO



END SUBROUTINE GetNelements



SUBROUTINE UPDATE(MYID,LEVEL,COMM,FLAG1)
USE LevelType
IMPLICIT NONE
INCLUDE 'mpif.h'
TYPE (LevelStruct)::Level
INTEGER :: MYID,TOTPS,IERR,stat(MPI_STATUS_SIZE),COMM,JJ,II,KK,P
INTEGER :: NDIM,I,J,K,NCPUS,CNT,FLAG,FLAG1
REAL*8::PlaneBuffSend(LEVEL%PlaneMaxDim+7,LEVEL%NB(MYID+1))
REAL*8::PlaneBuffSendLines(LEVEL%LineMaxDim+7,LEVEL%NB(MYID+1),4)
REAL*8::PlaneBuffSendPoints(5,LEVEL%NB(MYID+1),4)

 I=MYID+1
 PlaneBuffSend=0.
 LEVEL%PlaneBuffRecv=0.
 DO J=1,LEVEL%NB(I)
    PlaneBuffSend(1,J)=DBLE(LEVEL%PlanesSend(I,J))
    PlaneBuffSend(2,J)=DBLE(LEVEL%X1Plane(LEVEL%SenD(I,J)+1,LEVEL%PlanesSend(I,J)))
    PlaneBuffSend(3,J)=DBLE(LEVEL%X2Plane(LEVEL%SenD(I,J)+1,LEVEL%PlanesSend(I,J)))
    PlaneBuffSend(4,J)=DBLE(LEVEL%Y1Plane(LEVEL%SenD(I,J)+1,LEVEL%PlanesSend(I,J)))
    PlaneBuffSend(5,J)=DBLE(LEVEL%Y2Plane(LEVEL%SenD(I,J)+1,LEVEL%PlanesSend(I,J)))
    PlaneBuffSend(6,J)=DBLE(LEVEL%Z1Plane(LEVEL%SenD(I,J)+1,LEVEL%PlanesSend(I,J)))
    PlaneBuffSend(7,J)=DBLE(LEVEL%Z2Plane(LEVEL%SenD(I,J)+1,LEVEL%PlanesSend(I,J)))
!	DO P1

 ENDDO



 DO J=1,LEVEL%NB(I)
  CNT=7
  DO KK=INT(PlaneBuffSend(6,J)),INT(PlaneBuffSend(7,J))
    DO JJ=INT(PlaneBuffSend(4,J)),INT(PlaneBuffSend(5,J))
     DO II=INT(PlaneBuffSend(2,J)),INT(PlaneBuffSend(3,J))
	CNT=CNT+1 	
   IF(FLAG1.eq.0) PlaneBuffSend(CNT,J)=LEVEL%PHI(II,JJ,KK)
   IF(FLAG1.eq.1) PlaneBuffSend(CNT,J)=LEVEL%ERR(II,JJ,KK)
    ENDDO
   ENDDO
  ENDDO
 ENDDO



 DO J=1,LEVEL%NB(I)
   call  MPI_SendRecv(PlaneBuffSend(:,J),LEVEL%PlaneMaxDim+7,MPI_DOUBLE_PRECISION,LEVEL%SEND(I,J),MPI_ANY_TAG,LEVEL%PlaneBuffRecv(:,J),&
   & LEVEL%PlaneMaxDim+7,MPI_DOUBLE_PRECISION,LEVEL%RECV(I,J),MPI_ANY_TAG,COMM,stat,ierr)
 ENDDO



 PlaneBuffSendLines=0.
 LEVEL%LineBuffRecV=0.
  DO P=1,4	
   DO J=1,LEVEL%NB(I)
    PlaneBuffSendLines(1,J,P)=DBLE(LEVEL%LinesSend(I,J,P))
    PlaneBuffSendLines(2,J,P)=DBLE(LEVEL%X1Line(P,LEVEL%SenD1(I,J,P)+1,LEVEL%PlanesSend(I,J)))
    PlaneBuffSendLines(3,J,P)=DBLE(LEVEL%X2Line(P,LEVEL%SenD1(I,J,P)+1,LEVEL%PlanesSend(I,J)))
    PlaneBuffSendLines(4,J,P)=DBLE(LEVEL%Y1Line(P,LEVEL%SenD1(I,J,P)+1,LEVEL%PlanesSend(I,J)))
    PlaneBuffSendLines(5,J,P)=DBLE(LEVEL%Y2Line(P,LEVEL%SenD1(I,J,P)+1,LEVEL%PlanesSend(I,J)))
    PlaneBuffSendLines(6,J,P)=DBLE(LEVEL%Z1Line(P,LEVEL%SenD1(I,J,P)+1,LEVEL%PlanesSend(I,J)))
    PlaneBuffSendLines(7,J,P)=DBLE(LEVEL%Z2Line(P,LEVEL%SenD1(I,J,P)+1,LEVEL%PlanesSend(I,J)))
  ENDDO
 ENDDO



DO P=1,4	
 DO J=1,LEVEL%NB(I)
  CNT=7
   DO KK=INT(PlaneBuffSendLines(6,J,P)),INT(PlaneBuffSendLines(7,J,P))
    DO JJ=INT(PlaneBuffSendLines(4,J,P)),INT(PlaneBuffSendLines(5,J,P))
     DO II=INT(PlaneBuffSendLines(2,J,P)),INT(PlaneBuffSendLines(3,J,P))
	CNT=CNT+1 	
	   IF(FLAG1.eq.0) PlaneBuffSendLines(CNT,J,P)=LEVEL%PHI(II,JJ,KK)
	   IF(FLAG1.eq.1) PlaneBuffSendLines(CNT,J,P)=LEVEL%ERR(II,JJ,KK)
    ENDDO
   ENDDO
  ENDDO
  ENDDO
 ENDDO

DO P=1,4
 DO J=1,LEVEL%NB(I)
   call  MPI_SendRecv(PlaneBuffSendLines(:,J,P),LEVEL%LineMaxDim+7,MPI_DOUBLE_PRECISION,LEVEL%SEND1(I,J,P),MPI_ANY_TAG,LEVEL%LineBuffRecV(:,J,P),&
   & LEVEL%LineMaxDim+7,MPI_DOUBLE_PRECISION,LEVEL%RECV1(I,J,P),MPI_ANY_TAG,COMM,stat,ierr)
 ENDDO
ENDDO




 PlaneBuffSendPoints=0.
 LEVEL%PointBuffRecV=0.
  DO P=1,4	
   DO J=1,LEVEL%NB(I)
    PlaneBuffSendPoints(1,J,P)=DBLE(LEVEL%PointsSend(I,J,P))
    PlaneBuffSendPoints(2,J,P)=DBLE(LEVEL%XPoint(P,LEVEL%SenD2(I,J,P)+1,LEVEL%PlanesSend(I,J)))
    PlaneBuffSendPoints(3,J,P)=DBLE(LEVEL%YPoint(P,LEVEL%SenD2(I,J,P)+1,LEVEL%PlanesSend(I,J)))
    PlaneBuffSendPoints(4,J,P)=DBLE(LEVEL%ZPoint(P,LEVEL%SenD2(I,J,P)+1,LEVEL%PlanesSend(I,J)))
  ENDDO
 ENDDO



DO P=1,4	
 DO J=1,LEVEL%NB(I)
   KK=INT(PlaneBuffSendPoints(4,J,P))
   JJ=INT(PlaneBuffSendPoints(3,J,P))
   II=INT(PlaneBuffSendPoints(2,J,P))
   IF(FLAG1.eq.0)   PlaneBuffSendPoints(5,J,P)=LEVEL%PHI(II,JJ,KK)
   IF(FLAG1.eq.1)   PlaneBuffSendPoints(5,J,P)=LEVEL%ERR(II,JJ,KK)
 
  ENDDO
ENDDO

DO P=1,4
 DO J=1,LEVEL%NB(I)
   call  MPI_SendRecv(PlaneBuffSendPoints(:,J,P),5,MPI_DOUBLE_PRECISION,LEVEL%SEND2(I,J,P),MPI_ANY_TAG,LEVEL%PointBuffRecV(:,J,P),&
   & 5,MPI_DOUBLE_PRECISION,LEVEL%RECV2(I,J,P),MPI_ANY_TAG,COMM,stat,ierr)
 ENDDO
ENDDO


END SUBROUTINE UPDATE


SUBROUTINE SSORT (X, N)
      IMPLICIT NONE
      INTEGER N
      INTEGER X(N)
      INTEGER TEMP
      INTEGER I, J

	DO I=1,N-1
            DO J=I,N-1
              IF(X(I).LT.X(J+1)) THEN
                TEMP=X(I)
                X(I)=X(J+1)
                X(J+1)=TEMP
	      ENDIF
	    ENDDO
      ENDDO
END SUBROUTINE SSORT


SUBROUTINE ITERATE(LEVEL,flag,lamN,MYID,COMM,FLAG1)
USE LevelType
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER MYID,IERR,stat(MPI_STATUS_SIZE),COMM,JJ,II,KK,Df1,Df2,PNT,CNT,Df3,PlaneMaxDim
INTEGER NBI,FLAG1,NN
REAL*8  :: phiIm,phiIp,phiJm,phiJp,phiKm,phiKp
INTEGER:: CpuN
INTEGER NDIM,ITER,MAXITER,I,J,K,IP,IM,JP,JM,KP,KM,P,NDIM2,CNT2,CNT1,M,O,CNTTT
REAL*8 SUMMHARTREE,EPS1,EPS2,EPS3,EPS4,EPS5,EPS6,D_MAP,s8umm,summ_prev,error
REAL*8 coeff2,omega,coeff1,summ1,IND(7),lamN,Momega
LOGICAL TEST,FLAG
TYPE (LevelStruct)::Level

omega=2./(1.+SQRT(1.-lamN))
Momega=1.-omega



 NDIM=LEVEL%NDIM
 call UpdateBoundary(MYID,Level,COMM,FLAG1)
 CpuN=MYID+1
DO K=LEVEL%Z1(CpuN),LEVEL%Z2(CpuN)
  DO J=LEVEL%Y1(CpuN),LEVEL%Y2(CpuN)
	DO I=LEVEL%X1(CpuN),LEVEL%X2(CpuN)
	TEST=MOD((I+J+K),2).EQ.0
	IF(TEST.EQ.FLAG) THEN
	  ip=I+1
	  im=I-1
          jp=J+1
	  jm=J-1
	  kp=K+1
	  km=K-1
         

	  IF(I.EQ.NDIM) ip=1
	  IF(J.EQ.NDIM) jp=1
	  IF(K.EQ.NDIM) kp=1
	  IF(I.EQ.1) im=NDIM
	  IF(J.EQ.1) jm=NDIM
	  IF(K.EQ.1) km=NDIM	

          phiIm=LEVEL%phi(im,j,k)
	  phiIp=LEVEL%phi(ip,j,k)
	  phiJm=LEVEL%phi(i,jm,k)		
	  phiJp=LEVEL%phi(i,jp,k)		
	  phiKm=LEVEL%phi(i,j,km)
	  phiKp=LEVEL%phi(i,j,kp)
   	  eps1=LEVEL%eps(5,i,j,k)

	  d_map=1./(6.+eps1)
	  coeff1=LEVEL%phi(i,j,k)
	  coeff2=(phiIm+phiIp+phiJm+phiJp+phiKm+phiKp)*d_map
	  coeff2=coeff2*omega+Momega*coeff1+omega*LEVEL%rhs(i,j,k)*d_map
	  LEVEL%phi(i,j,k)=coeff2

	ENDIF
    ENDDO
   ENDDO
 ENDDO

END SUBROUTINE ITERATE


SUBROUTINE UpdateBoundary(MYID,LeveL,COMM,FLAG1)
USE LevelType
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER MYID,IERR,stat(MPI_STATUS_SIZE),COMM,JJ,II,KK,Df1,Df2,PNT,CNT,Df3
INTEGER NBI,FLAG1,NN
INTEGER::X1,X2,Y1,Y2,Z1,Z2,L
REAL*8  :: phiIm,phiIp,phiJm,phiJp,phiKm,phiKp
INTEGER:: CpuN
INTEGER NDIM,ITER,MAXITER,I,J,K,IP,IM,JP,JM,KP,KM,P,NDIM2,CNT2,CNT1,M,O,CNTTT
REAL*8 SUMMHARTREE,EPS1,EPS2,EPS3,EPS4,EPS5,EPS6,D_MAP,summ,summ_prev,error
REAL*8 coeff2,omega,coeff1,summ1,IND(7),lamN,Momega
TYPE (LevelStruct)::Level

 NDIM=LEVEL%NDIM
 CpuN=MYID+1

DO P=1,LEVEL%NB(CpuN)
  DO K=INT(LEVEL%PlaneBuffRecV(6,P)),INT(LEVEL%PlaneBuffRecV(7,P))
    DO J=INT(LEVEL%PlaneBuffRecV(4,P)),INT(LEVEL%PlaneBuffRecV(5,P))
     DO I=INT(LEVEL%PlaneBuffRecV(2,P)),INT(LEVEL%PlaneBuffRecV(3,P))

	  IF(INT(LEVEL%PlaneBuffRecv(1,P)).eq.1) then	
			Df1=ABS(INT(LEVEL%PlaneBuffRecV(4,P))-INT(LEVEL%PlaneBuffRecV(5,P)))+1
			PNT=(J-INT(LEVEL%PlaneBuffRecV(4,P)))+Df1*(K-INT(LEVEL%PlaneBuffRecV(6,P)))+8
	 	        phiIm=LEVEL%PlaneBuffRecV(PNT,1)
		        IF(FLAG1.eq.0)   LEVEL%phi(i,j,k)=phiIm
		        IF(FLAG1.eq.1)   LEVEL%err(i,j,k)=phiIm
	  ENDIF	
	
	  IF(INT(LEVEL%PlaneBuffRecv(1,P)).eq.2) then	
			Df1=ABS(INT(LEVEL%PlaneBuffRecV(4,P))-INT(LEVEL%PlaneBuffRecV(5,P)))+1
			PNT=(J-INT(LEVEL%PlaneBuffRecV(4,P)))+Df1*(K-INT(LEVEL%PlaneBuffRecV(6,P)))+8
		        phiIp=LEVEL%PlaneBuffRecV(PNT,2)
			IF(FLAG1.eq.0)  LEVEL%phi(i,j,k)=phiIp
			IF(FLAG1.eq.1)  LEVEL%err(i,j,k)=phiIp
	  ENDIF	

	  IF(INT(LEVEL%PlaneBuffRecv(1,P)).eq.3) then	
			Df3=ABS(INT(LEVEL%PlaneBuffRecV(2,P))-INT(LEVEL%PlaneBuffRecV(3,P)))+1
			PNT=(I-INT(LEVEL%PlaneBuffRecV(2,P)))+Df3*(K-INT(LEVEL%PlaneBuffRecV(6,P)))+8
		        phiJm=LEVEL%PlaneBuffRecV(PNT,3)
		        IF(FLAG1.eq.0) LEVEL%phi(i,j,k)=phiJm
			IF(FLAG1.eq.1) LEVEL%err(i,j,k)=phiJm
	  ENDIF	

	  IF(INT(LEVEL%PlaneBuffRecv(1,P)).eq.4) then	
			Df3=ABS(INT(LEVEL%PlaneBuffRecV(2,P))-INT(LEVEL%PlaneBuffRecV(3,P)))+1
			PNT=(I-INT(LEVEL%PlaneBuffRecV(2,P)))+Df3*(K-INT(LEVEL%PlaneBuffRecV(6,P)))+8
			phiJp=LEVEL%PlaneBuffRecV(PNT,4)
		        IF(FLAG1.eq.0) LEVEL%phi(i,j,k)=phiJp
			IF(FLAG1.eq.1) LEVEL%err(i,j,k)=phiJp
	  ENDIF	

	  IF(INT(LEVEL%PlaneBuffRecv(1,P)).eq.5) then	
			Df3=ABS(INT(LEVEL%PlaneBuffRecV(2,P))-INT(LEVEL%PlaneBuffRecV(3,P)))+1
			PNT=(I-INT(LEVEL%PlaneBuffRecV(2,P)))+Df3*(J-INT(LEVEL%PlaneBuffRecV(4,P)))+8
			phiKm=LEVEL%PlaneBuffRecV(PNT,5)
			IF(FLAG1.eq.0)  LEVEL%phi(i,j,k)=phiKm
			IF(FLAG1.eq.1)  LEVEL%err(i,j,k)=phiKm
	  ENDIF	

	  IF(INT(LEVEL%PlaneBuffRecv(1,P)).eq.6) then	
			Df3=ABS(INT(LEVEL%PlaneBuffRecV(2,P))-INT(LEVEL%PlaneBuffRecV(3,P)))+1
			PNT=(I-INT(LEVEL%PlaneBuffRecV(2,P)))+Df3*(J-INT(LEVEL%PlaneBuffRecV(4,P)))+8
		        phiKp=LEVEL%PlaneBuffRecV(PNT,6)
			IF(FLAG1.eq.0) LEVEL%phi(i,j,k)=phiKp
			IF(FLAG1.eq.1) LEVEL%err(i,j,k)=phiKp
	  ENDIF	
    ENDDO
   ENDDO
  ENDDO
 ENDDO

 DO J=1,LEVEL%NB(CpuN)
  DO P=1,4
   CNT=7	
   KK=INT(LEVEL%PointBuffRecV(4,J,P))
   JJ=INT(LEVEL%PointBuffRecV(3,J,P))
   II=INT(LEVEL%PointBuffRecV(2,J,P))
     IF(FLAG1.eq.0)   LEVEL%phi(II,JJ,KK)=LEVEL%PointBuffRecV(5,J,P)
     IF(FLAG1.eq.1)   LEVEL%err(II,JJ,KK)=LEVEL%PointBuffRecV(5,J,P)
   DO KK=INT(LEVEL%LineBuffRecV(6,J,P)),INT(LEVEL%LineBuffRecV(7,J,P))
     DO JJ=INT(LEVEL%LineBuffRecV(4,J,P)),INT(LEVEL%LineBuffRecV(5,J,P))
      DO II=INT(LEVEL%LineBuffRecV(2,J,P)),INT(LEVEL%LineBuffRecV(3,J,P))
	CNT=CNT+1
	IF(FLAG1.eq.0)	LEVEL%phi(II,JJ,KK)=LEVEL%LineBuffRecV(CNT,J,P)
	IF(FLAG1.eq.1)	LEVEL%err(II,JJ,KK)=LEVEL%LineBuffRecV(CNT,J,P)	
   ENDDO
   ENDDO
  ENDDO
 ENDDO
ENDDO

END SUBROUTINE UpdateBoundary
