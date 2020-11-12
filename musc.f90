Module MUSC


use numerics
use flow
use rusanov

contains

subroutine MS_flow(U,wg,wd)

 real(rp) , dimension(:,:), intent(in) :: U
 real(rp) , dimension(Ns,3), intent(out) :: wg, wd
 real(rp) , dimension(Ns,3) :: sigma

 integer :: i
 
call minmod_sys(u,sigma)
 
  DO i=1,Ns-1
    
    wg(i,:)=u(i,:)+0.5_rp*dx_sys*sigma(i,:)
    wd(i,:)=u(i+1,:)-0.5_rp*dx_sys*sigma(i+1,:)
 
 END DO
      
    wg(Ns,:)=u(Ns,:)-0.5_rp*dx_sys*sigma(Ns,:)
    wd(Ns,:)=u(Ns,:)
end subroutine MS_flow




Subroutine MUSCL_Solve() ! Solveur MUSCL


!-----------------------------VARIABLES ----------------------------
    integer :: i,step, cpt
    real(rp) , dimension(:), allocatable::  tmp
    real(rp) , dimension(:,:), allocatable:: Fluxg,Fluxd, wl, wr

    real(rp) ::  t
    character(len=80) :: FileName
!----------------INITIALISATION-----------------------------
    
    
    t=zero
    cpt=0

    
    
!-----------------ALLOCATION--------------------------------------------
    
    allocate(Fluxg(Ns,3), Fluxd(Ns,3), wl(Ns,3), wr(Ns,3))
    allocate(tmp(Ns))


!----------------BOUCLE EN TEMPS----------------------------------------
step=1
DO WHILE (t<Tm .AND. step < 10000)

!----------------VECTEUR DE FLUX----------------------------------------
    
    CALL MS_flow(U_prevs,wl,wr)
    Call HLL_flow(wr(2,:),wl(1,:),Fluxg(1,:))
    Call HLL_flow(wr(1,:),wl(1,:),Fluxd(1,:))
    
    tmp(1) =eigen_sys(U_prevs(1,:), U_prevs(2,:) )    
    
    DO i=2,Ns-1
        
        Call HLL_flow(wr(i,:), wl(i-1,:),Fluxg(i,:))
        Call HLL_flow(wr(i-1,:),wl(i,:),Fluxd(i,:))
        tmp(i) =eigen_sys(U_prevs(i,:), U_prevs(i+1,:) )   
                 
    END DO 
 
     
   call HLL_flow(U_prevs(Ns-1,:),U_prevs(Ns,:),Fluxg(Ns,:))
   call HLL_flow(U_prevs(Ns,:),U_prevs(Ns,:),Fluxd(Ns,:))
   
   tmp(Ns) =eigen_sys(U_prevs(Ns-1,:), U_prevs(Ns,:))
   
!---------------Dt EN FONCTION DE LA CFL-------------------------------
   
    dt_sys=0.5_rp*dx_sys/maxval(tmp)
  
!----------------------Iteration en temps ------------------------------


    DO i=2,Ns-1
       
        U_nexts(i,:)=U_prevs(i,:)-dt_sys/dx_sys*(Fluxd(i+1,:)-Fluxd(i,:))
       
    END DO 
    
    U_nexts(1,:)=U_nexts(2,:)
    U_nexts(Ns,:)=zero

    U_prevs=U_nexts

    t=t+dt_sys

    step=step+1
    
    if (mod(step,10)==0)then
     write(FileName,'(A,I3.3,A)') 'resultat/hauteur_ms',cpt,'.txt'
     CALL save_file(X_sys,U_nexts(:,1),FileName)
     cpt=cpt+1
    end if

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
  

  minmod(1,1)=max(zero,min(un,(u(1,1)/dx_sys), (u(2,1)-u(1,1))/dx_sys))
  minmod(1,2)=max(zero,min(un,(u(1,2)/dx_sys), (u(2,2)-u(1,2))/dx_sys))
  minmod(1,3)=max(zero,min(un,(u(1,3)/dx_sys), (u(2,3)-u(1,3))/dx_sys))

  DO i=2,Ns-1 
         minmod(i,3)=max(zero,min(un,(u(i,1)-u(i-1,1))/dx_sys), (u(i+1,1)-u(i,1))/dx_sys)
         minmod(i,3)=max(zero,min(un,(u(i,2)-u(i-1,2))/dx_sys), (u(i+1,2)-u(i,2))/dx_sys)
         minmod(i,3)=max(zero,min(un,(u(i,3)-u(i-1,3))/dx_sys), (u(i+1,3)-u(i,3))/dx_sys)
  END DO
  
  minmod(Ns,1)=max(zero,min(un,(u(Ns,1)-u(Ns-1,1))/dx_sys), u(Ns,1)/dx_sys)
  minmod(Ns,2)=max(zero,min(un,(u(Ns,2)-u(Ns-1,2))/dx_sys), u(Ns,2)/dx_sys)
  minmod(Ns,3)=max(zero,min(un,(u(Ns,3)-u(Ns-1,3))/dx_sys), u(Ns,3)/dx_sys)

return
end subroutine MINMOD_SYS


end module musc
