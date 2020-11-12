module rusanov

use numerics
use flow
use omp_lib

!------------------------------Systeme----------------------------------
     
      REAL(rp), dimension(:), allocatable ::  X_sys
      REAL(rp), dimension(:,:), allocatable :: u_prevs, u_nexts
      REAL(rp) :: a,b, dx_sys, tm, dt_sys
      INTEGER ::  ns
     
!-----------------------------------------------------------------------



contains 

!------------------------Initialisation---------------------------------

subroutine init_syst(space_number,final_time)

    real(rp) , intent(in) :: final_time
    integer, intent(in) :: space_number
    
      Ns=space_number
       Tm=final_time
!--------------------------Allocation----------------------------------    
    
        allocate(u_prevs(Ns,3), u_nexts(Ns,3), X_sys(Ns))

!----------------------------Constantes--------------------------------


        a=-1._rp
        b=1._rp
        dx_sys=(b-a)/Ns
       


!-----------------Initialisation donnée intiale-------------------------

    Do i=1,Ns

        X_sys(i)=a+i*dx_sys  ! espace

        u_prevs(i,1)=init_dam(x_sys(i))! hauteur
        u_prevs(i,2)=0._rp!un/(un+abs(x_sys(i)))! u vitesse
        u_prevs(i,3)=5._rp ! temperature

    END DO 


end subroutine init_syst


function init_dam(xx) !donnee initiale
  real(rp) , intent(in):: xx
  real(rp) :: init_dam

if ( xx < 0._rp ) then
    
    init_dam=un
    
else
     init_dam=zero
     
end if  

return
end function init_dam

!---------------------------Valeurs propres-----------------------------

function eigen_min(Ug,Ud)


real(rp) , dimension(3), intent(in) :: Ug, Ud
real(rp) ::lam_g , lam_d, eigen_min



lam_g=min( Ug(1)-sqrt(grv*Ug(3)*Ug(1)),Ug(2), Ug(2)+sqrt(grv*Ug(3)*Ug(1)))
lam_d=min( Ud(1)-sqrt(grv*Ud(3)*Ud(1)),Ud(2), Ud(2)+sqrt(grv*Ud(3)*Ud(1)))

eigen_min=min(lam_g , lam_d)

return
end function eigen_min

function eigen_max(Ug,Ud)


real(rp) , dimension(3), intent(in) :: Ug, Ud
real(rp) ::lam_g , lam_d, eigen_max




lam_g=max( Ug(1)-sqrt(grv*Ug(3)*Ug(1)),Ug(2), Ug(2)+sqrt(grv*Ug(3)*Ug(1)))
lam_d=max( Ud(1)-sqrt(grv*Ud(3)*Ud(1)),Ud(2), Ud(2)+sqrt(grv*Ud(3)*Ud(1)))

eigen_max=max(lam_g , lam_d)

return
end function eigen_max


function eigen_sys(Ug,Ud)


real(rp) , dimension(3), intent(in) :: Ug, Ud
real(rp) ::lam_g , lam_d, eigen_sys




lam_g=max( abs(Ug(1)-sqrt(grv*Ug(3)*Ug(1))),abs(Ug(2)), abs(Ug(2)+sqrt(grv*Ug(3)*Ug(1))))
lam_d=max( abs(Ud(1)-sqrt(grv*Ud(3)*Ud(1))),abs(Ud(2)), abs(Ud(2)+sqrt(grv*Ud(3)*Ud(1))))

eigen_sys=max(lam_g , lam_d)

return
end function eigen_sys

!-----------------------------Flux--------------------------------------

subroutine F_sys(X,Y)

real(rp) , dimension(3), intent(in) :: X
real(rp) , dimension(3), intent(out) :: Y




Y(1)=X(2)*X(1)
Y(2)=X(2)**2+0.5*g*X(2)*X(3)*X(1)**2
Y(3)=X(2)*X(3)*X(1)




end subroutine F_sys


subroutine rusanov_flow(Ug,Ud,F)


real(rp) , dimension(3), intent(in) :: Ug, Ud
real(rp) , dimension(3), intent(out) :: F
real(rp) , dimension(3) :: Fg , Fd
real(rp) :: cc


cc=eigen_sys(Ug,Ud)

call F_sys(Ug,Fg)
call F_sys(Ud,Fd)





F(1)= 0.5_rp*(Fg(1)+Fd(1))-cc*( Ud(1)-Ug(1))*0.5_rp
F(2)= 0.5_rp*(Fg(2)+Fd(2)) -cc*( Ud(2)-Ug(2))*0.5_rp
F(3)= 0.5_rp*(Fg(3)+Fd(3))-cc*(Ud(3)-Ug(3))*0.5_rp




end subroutine

!----------------------------------------------------------------------



subroutine HLL_flow(Ug,Ud,F)


real(rp) , dimension(3), intent(in) :: Ug, Ud
real(rp) , dimension(3), intent(out) :: F
real(rp) , dimension(3) :: Fg , Fd
real(rp) :: lam_min,lam_max


!--------------------min/max des valeurs propres------------------------
lam_min=eigen_min(Ug,Ud)
lam_max=eigen_max(Ug,Ud)
!--------------------Calcul de F---------------------------------------

call F_sys(Ug,Fg)
call F_sys(Ud,Fd)

!----------------Flux en fonction Val. propores-------------------------

    if (lam_min < zero .and. zero <lam_max ) then

    
            F=(lam_max*Fg-lam_min*Fd +(Ud-Ug)*(lam_min*lam_max))/(lam_max-lam_min)
    
       
    else if (lam_max < zero) then
            
            F=Fd
            
    else if ( zero  < lam_min) then
            
            F=Fg
       
    end if 
    
     


end subroutine HLL_flow



!--------------------------Solveurs-------------------------------------

Subroutine RUS_Solve() ! Solveur basé sur Rusanov


!-----------------------------VARIABLES --------------------------------
    integer :: i, cpt, step
    real(rp) , dimension(:), allocatable::  tmp
    real(rp) , dimension(:,:), allocatable:: Fluxx
    real(rp) ::    t
    character(len=80) :: FileName

  
!-------------------------INITIALISATION--------------------------------
    
    
    t=zero
!----------------------------ALLOCATION---------------------------------
    allocate(Fluxx(Ns,3), tmp(Ns))


!------------------------BOUCLE EN TEMPS--------------------------------
cpt=0
step=1
DO WHILE (t<Tm )

!--------------------VECTEUR DE FLUX-----------------------------------
    
     
   
   
!$OMP PARALLEL
    
    DO i=2,Ns-1
        
        Call rusanov_flow(U_prevs(i,:),U_prevs(i+1,:),Fluxx(i,:))
        tmp(i) =eigen_sys(U_prevs(i,:), U_prevs(i+1,:) )
            
    END DO 
!$OMP END PARALLEL
!----------------------------BORD---------------------------------------
   
    Call rusanov_flow(U_prevs(Ns,:),U_prevs(1,:),Fluxx(1,:))
    tmp(1) =eigen_sys(U_prevs(Ns,:), U_prevs(1,:) )
     

   
   
   
   
   
   
   
   
!-------------------Dt EN FONCTION DE LA CFL----------------------------
   
    dt_sys=0.5_rp*dx_sys/maxval(tmp)
  
!----------------------ITERATION EN ESPACE------------------------------

!$OMP PARALLEL
    DO i=2,Ns-1
       
        U_nexts(i,:)=U_prevs(i,:)-dt_sys/dx_sys*(Fluxx(i,:)-Fluxx(i-1,:))
       
    END DO 
!$OMP END PARALLEL  
   U_nexts(1,:)=U_prevs(1,:)-dt_sys/dx_sys*(Fluxx(2,:)-Fluxx(1,:))
   U_nexts(Ns,:)=U_prevs(Ns,:)-dt_sys/dx_sys*(Fluxx(1,:)-Fluxx(Ns-1,:))

    U_prevs=U_nexts
    
!----------------------------------------------------------------------
    t=t+dt_sys
    step=step+1
    
    if (mod(step,1000)==0)then
     write(FileName,'(A,I4.4,A)') 'resultat/hauteur_rus',cpt,'.txt'
     CALL save_file(X_sys,U_nexts(:,1),FileName)
     cpt=cpt+1
  
   end if


END DO

!------------------------ Sauvegarde------------------------------------

CALL save_file(X_sys,U_nexts(:,1),'hauteur_rus.txt')
CALL save_file(X_sys,U_nexts(:,2),'vitesse_rus.txt')
CALL save_file(X_sys,U_nexts(:,3),'temperature_rus.txt')


end Subroutine RUS_Solve



Subroutine HLL_Solve() ! Solveur basé sur HLL


!-----------------------------VARIABLES ----------------------------
    integer :: i,step, cpt
    real(rp) , dimension(:), allocatable::  tmp
    real(rp) , dimension(:,:), allocatable:: Fluxx
    real(rp) ::  t
    character(len=80) :: FileName
!----------------INITIALISATION-----------------------------
    
    
    t=zero
    cpt=0

    
    
!-----------------ALLOCATION--------------------------------------------
    
    allocate(Fluxx(Ns,3),tmp(Ns))


!----------------BOUCLE EN TEMPS----------------------------------------
step=1

 write(FileName,'(A,I4.4,A)') 'resultat/hauteur_hll',cpt,'.txt'
 CALL save_file(X_sys,U_prevs(:,1),FileName)

cpt=1

DO WHILE (t<Tm )

!----------------VECTEUR DE FLUX----------------------------------------
    
   
!$OMP PARALLEL 
    DO i=2,Ns-1
        
        Call HLL_flow(U_prevs(i-1,:),U_prevs(i,:),Fluxx(i-1,:))
        tmp(i) =eigen_sys(U_prevs(i-1,:), U_prevs(i+1,:) )
            
    END DO 
!$OMP END PARALLEL
!----------------------------BORD---------------------------------------
   
   Call HLL_flow(U_prevs(Ns,:),U_prevs(1,:),Fluxx(1,:))
    tmp(1) =eigen_sys(U_prevs(Ns,:), U_prevs(1,:) )
    
     
   call HLL_flow(U_prevs(Ns-1,:),U_prevs(Ns,:),Fluxx(Ns-1,:))
   tmp(Ns) =eigen_sys(U_prevs(Ns-1,:), U_prevs(Ns,:))
   
!---------------Dt EN FONCTION DE LA CFL-------------------------------
   
    dt_sys=0.25_rp*dx_sys/maxval(tmp)
  
!----------------------ITERATION EN ESPACE ------------------------------

    
!$OMP PARALLEL
    DO i=2,Ns-1
       
           U_nexts(i,:)=U_prevs(i,:)-dt_sys/dx_sys*(Fluxx(i,:)-Fluxx(i-1,:))
       
    END DO 
!$OMP END PARALLEL
    
   U_nexts(1,:)=U_prevs(1,:)-dt_sys/dx_sys*(Fluxx(2,:)-Fluxx(1,:))
  U_nexts(Ns,:)=U_prevs(Ns,:)-dt_sys/dx_sys*(Fluxx(1,:)-Fluxx(Ns-1,:))

    

    U_prevs=U_nexts

    t=t+dt_sys

    step=step+1
    
   ! if (mod(step,1000)==0)then
    ! write(FileName,'(A,I4.4,A)') 'resultat/hauteur_hll',cpt,'.txt'
  !   CALL save_file(X_sys,U_nexts(:,1),FileName)
   !  cpt=cpt+1
  
  ! end if

END DO

!-----------------------Sauvegarde-------------------------------------

     
     !CALL save_file(X_sys,U_nexts(:,1),FileName)
CALL save_file(X_sys,U_nexts(:,1),'hauteur_hll.txt')
!CALL save_file(X_sys,U_nexts(:,2),'vitesse_hll.txt')
!CALL save_file(X_sys,U_nexts(:,3),'temperature_hll.txt')


end Subroutine HLL_Solve






end module rusanov
