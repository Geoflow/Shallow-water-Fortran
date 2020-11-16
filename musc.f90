Module MUSC


use numerics
use flow
use rusanov
use omp_lib


contains

subroutine MS_STATE(U,wg,wd)

 real(rp) , dimension(:,:), intent(in) :: U
 real(rp) , dimension(Ns,3), intent(out) :: wg, wd
 real(rp) , dimension(Ns,3) :: sigma

 integer :: i
 
call minmod_sys(u,sigma)
 
!$OMP PARALLEL   

  DO i=2,Ns-1
    
    wg(i,:)=u(i,:)-0.5_rp*dx_sys*sigma(i,:)*(u(i,:)-u(i-1,:))
    wd(i,:)=u(i,:)+0.5_rp*dx_sys*sigma(i,:)*(u(i+1,:)-u(i,:))
 
 END DO
!$OMP END PARALLEL   
  
     wg(1,:)=u(1,:)-0.5_rp*dx_sys*sigma(Ns-1,:)*(u(2,:)-u(1,:))
     wd(1,:)=u(1,:)-0.5_rp*dx_sys*sigma(1,:)*(u(2,:)-u(1,:))  

     wg(Ns,:)=u(Ns,:)-0.5_rp*dx_sys*sigma(Ns,:)
     wd(Ns,:)=u(Ns,:)+0.5_rp*dx_sys*sigma(Ns,:)
    
end subroutine MS_STATE




Subroutine MUSCL_Solve() ! Solveur MUSCL


!-----------------------------VARIABLES --------------------------------
    integer :: i,step, cpt
    real(rp) , dimension(:), allocatable::  tmp
    real(rp) , dimension(:,:), allocatable:: Fluxg,Fluxd, wl, wr

    real(rp) ::  t
    character(len=80) :: FileName

    
!-------------------------ALLOCATION------------------------------------
    
    allocate(Fluxg(Ns,3), Fluxd(Ns,3), wl(Ns,3), wr(Ns,3))
    allocate(tmp(Ns))

!--------------------------INITIALISATION-------------------------------
    
    
    t=zero
    cpt=0
   wl=U_prevs
   wr=U_prevs
    

!-------------------------BOUCLE EN TEMPS-------------------------------
step=1
DO WHILE (t<Tm .AND. step < 10000)

!----------------VECTEUR DE FLUX----------------------------------------
    
    !Call HLL_flow(wr(2,:),wl(1,:),Fluxg(1,:))
   ! Call HLL_flow(wr(1,:),wl(1,:),Fluxd(1,:))
    
 !$OMP PARALLEL   
    DO i=2,Ns-1
        
        Call HLL_flow(wr(i-1,:), wl(i,:),Fluxg(i,:))
        Call HLL_flow(wr(i,:),wl(i+1,:),Fluxd(i,:))
        tmp(i) =eigen_sys(U_prevs(i,:), U_prevs(i+1,:) )   
                 
    END DO 
 !$OMP END PARALLEL
     
   !call HLL_flow(U_prevs(Ns-1,:),U_prevs(Ns,:),Fluxg(Ns,:))
   !call HLL_flow(U_prevs(Ns,:),U_prevs(Ns,:),Fluxd(Ns,:))
   
!----------------------------BORD---------------------------------------
   
   Call HLL_flow(U_prevs(Ns,:),U_prevs(1,:),Fluxg(1,:))
   Call HLL_flow(U_prevs(1,:),U_prevs(2,:),Fluxd(1,:))
    tmp(1) =eigen_sys(U_prevs(Ns,:), U_prevs(1,:) )
    
     
   call HLL_flow(U_prevs(Ns-1,:),U_prevs(Ns,:),Fluxg(Ns,:))
   tmp(Ns) =eigen_sys(U_prevs(Ns-1,:), U_prevs(Ns,:))
   
!---------------Dt EN FONCTION DE LA CFL-------------------------------
   
    dt_sys=0.25_rp*dx_sys/maxval(tmp)
  
!----------------------ITERATION EN ESPACE ------------------------------

!$OMP PARALLEL
    DO i=2,Ns-1
       
        U_nexts(i,:)=U_prevs(i,:)-dt_sys/dx_sys*(Fluxd(i,:)-Fluxg(i,:))
       
    END DO 
!$OMP END PARALLEL   

    U_nexts(1,:)=U_prevs(1,:)-dt_sys/dx_sys*(Fluxd(1,:)-Fluxg(1,:))
    U_nexts(Ns,:)=U_prevs(Ns,:)-dt_sys/dx_sys*(Fluxg(Ns,:)-Fluxd(Ns-1,:))

    U_prevs=U_nexts
    CALL MS_STATE(U_prevs,wl,wr)

    t=t+dt_sys

    step=step+1
    
    !if (mod(step,1000)==0)then
    ! write(FileName,'(A,I3.3,A)') 'resultat/hauteur_ms',cpt,'.txt'
    ! CALL save_file(X_sys,U_nexts(:,1),FileName)
    ! cpt=cpt+1
    !end if

END DO

!-----------------------Sauvegarde-------------------------------------

     !write(FileName,'(A,I4.4,A)') 'resultat/hauteur_hll',cpt,'.txt'
     !CALL save_file(X_sys,U_nexts(:,1),FileName)
CALL save_file(X_sys,U_nexts(:,1),'hauteur_ms.txt')
CALL save_file(X_sys,U_nexts(:,2),'vitesse_ms.txt')
CALL save_file(X_sys,U_nexts(:,3),'temperature_ms.txt')


end Subroutine MUSCL_Solve




!-----------------------LIMITEUR DE PENTE------------------------------


subroutine MINMOD_SYS(u,minmod)
  real(rp) ,dimension(:,:),intent(in):: u
  real(rp) ,dimension(Ns,3), intent(out) ::  minmod
  integer :: i
  

IF( (u(1,1)-u(Ns,1))/dx_sys>zero .AND. (u(2,1)-u(1,1))/dx_sys>zero)THEN
 
		minmod(1,1)=min((u(1,1)-u(Ns,1))/dx_sys, (u(2,1)-u(1,1))/dx_sys)
		
		ELSE IF(  (u(1,1)-u(Ns,1))/dx_sys<zero .AND. (u(2,1)-u(1,1))/dx_sys<zero )THEN
		
		minmod(1,1)=max((u(1,1)-u(Ns,1))/dx_sys, (u(2,1)-u(1,1))/dx_sys)
		
		ELSE
		minmod(1,1)=zero
		END IF      
 !------------------------------------------------------        
        IF( (u(1,2)-u(Ns,2))/dx_sys>zero .AND. (u(2,2)-u(1,2))/dx_sys>zero)THEN
        
		minmod(1,2)=min((u(1,2)-u(Ns,2))/dx_sys, (u(2,2)-u(1,2))/dx_sys)
		
		ELSE IF(  (u(1,2)-u(Ns,2))/dx_sys<zero .AND. (u(2,2)-u(1,2))/dx_sys<zero )THEN
		
		minmod(1,2)=max( (u(1,2)-u(Ns,2))/dx_sys, (u(2,2)-u(1,2))/dx_sys)
		
		ELSE
		minmod(1,2)=zero
		END IF
!--------------------------------------------------------		
		IF( (u(1,3)-u(Ns,3))/dx_sys>zero .AND. (u(2,3)-u(1,3))/dx_sys>zero)THEN
		
		minmod(1,3)=min((u(1,3)-u(Ns,3))/dx_sys, (u(2,3)-u(1,3))/dx_sys)
		
		ELSE IF(  (u(1,3)-u(Ns,3))/dx_sys<zero .AND. (u(2,3)-u(1,3))/dx_sys<zero )THEN
		
		minmod(1,3)=max((u(1,3)-u(Ns,3))/dx_sys, (u(2,3)-u(1,3))/dx_sys)
		
		ELSE
		minmod(1,3)=zero
		END IF
 !$OMP PARALLEL   

  DO i=2,Ns-1 
  
		IF( (u(i,1)-u(i-1,1))/dx_sys>zero .AND. (u(i+1,1)-u(i,1))/dx_sys>zero)THEN
		minmod(i,1)=min((u(i,1)-u(i-1,1))/dx_sys, (u(i+1,1)-u(i,1))/dx_sys)
		ELSE IF(  (u(i,1)-u(i-1,1))/dx_sys<zero .AND. (u(i+1,1)-u(i,1))/dx_sys<zero )THEN
		
		minmod(i,1)=max((u(i,1)-u(i-1,1))/dx_sys, (u(i+1,1)-u(i,1))/dx_sys)
		ELSE
		minmod(i,1)=zero
		END IF
         
         
        IF( (u(i,2)-u(i-1,2))/dx_sys>zero .AND. (u(i+1,2)-u(i,2))/dx_sys>zero)THEN
		minmod(i,2)=min((u(i,2)-u(i-1,2))/dx_sys, (u(i+1,2)-u(i,2))/dx_sys)
		ELSE IF(  (u(i,2)-u(i-1,2))/dx_sys<zero .AND. (u(i+1,2)-u(i,2))/dx_sys<zero )THEN
		
		minmod(i,2)=max( (u(i,2)-u(i-1,2))/dx_sys, (u(i+1,2)-u(i,2))/dx_sys)
		ELSE
		minmod(i,2)=zero
		END IF
		
		IF( (u(i,3)-u(i-1,3))/dx_sys>zero .AND. (u(i+1,3)-u(i,3))/dx_sys>zero)THEN
		minmod(i,3)=min((u(i,3)-u(i-1,3))/dx_sys, (u(i+1,3)-u(i,3))/dx_sys)
		ELSE IF(  (u(i,3)-u(i-1,3))/dx_sys<zero .AND. (u(i+1,3)-u(i,3))/dx_sys<zero )THEN
		
		minmod(i,3)=max((u(i,3)-u(i-1,3))/dx_sys, (u(i+1,3)-u(i,3))/dx_sys)
		ELSE
		minmod(i,3)=zero
		END IF
         
         
         
         
         
         
  END DO
 !$OMP END PARALLEL   
 
 IF( (u(Ns,1)-u(Ns-1,1))/dx_sys>zero .AND. (u(1,1)-u(Ns,1))/dx_sys>zero)THEN
 
		minmod(Ns,1)=min((u(Ns,1)-u(Ns-1,1))/dx_sys, (u(1,1)-u(Ns,1))/dx_sys)
		
		ELSE IF(  (u(Ns,1)-u(Ns-1,1))/dx_sys<zero .AND. (u(1,1)-u(Ns,1))/dx_sys<zero )THEN
		
		minmod(Ns,1)=max((u(Ns,1)-u(Ns-1,1))/dx_sys, (u(1,1)-u(Ns,1))/dx_sys)
		ELSE
		minmod(Ns,1)=zero
		END IF      
 !------------------------------------------------------        
        IF( (u(Ns,2)-u(Ns-1,2))/dx_sys>zero .AND. (u(1,2)-u(Ns,2))/dx_sys>zero)THEN
		minmod(Ns,2)=min((u(Ns,2)-u(Ns-1,2))/dx_sys, (u(1,2)-u(Ns,2))/dx_sys)
		
		ELSE IF(  (u(Ns,2)-u(Ns-1,2))/dx_sys<zero .AND. (u(1,2)-u(Ns,2))/dx_sys<zero )THEN
		
		minmod(Ns,2)=max( (u(Ns,2)-u(Ns-1,2))/dx_sys, (u(1,2)-u(Ns,2))/dx_sys)
		ELSE
		minmod(Ns,2)=zero
		END IF
!--------------------------------------------------------		
		IF( (u(Ns,3)-u(Ns-1,3))/dx_sys>zero .AND. (u(1,3)-u(Ns,3))/dx_sys>zero)THEN
		minmod(Ns,3)=min((u(Ns,3)-u(Ns-1,3))/dx_sys, (u(1,3)-u(Ns,3))/dx_sys)
		ELSE IF(  (u(Ns,3)-u(Ns-1,3))/dx_sys<zero .AND. (u(1,3)-u(Ns,3))/dx_sys<zero )THEN
		
		minmod(Ns,3)=max((u(Ns,3)-u(Ns-1,3))/dx_sys, (u(1,3)-u(Ns,3))/dx_sys)
		ELSE
		minmod(Ns,3)=zero
		END IF
return
end subroutine MINMOD_SYS


end module musc
