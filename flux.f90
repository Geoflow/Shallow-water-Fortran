module flow

use numerics
implicit none




contains


function finit(x)
  real(rp) , intent(in):: x
  real(rp) :: finit

if ( x < 0._rp ) then
    
    finit=un
    
else
     finit=-un
     
end if  


return
end function finit

function finit_bis(x)
  real(rp) , intent(in):: x
  real(rp) :: finit_bis

if ( x <=0._rp ) then

    finit_bis=0.5_rp
    
   else
    
    finit_bis=zero
     
end if  


return
end function finit_bis



subroutine finit_save()
!---------------INITIALISATION DONNEE INITIALE-------------------------
   integer :: i
   
    
    do i=1, N
    
       u_prev(i)=finit(x(i))  
       
        
    end do
    
!-------------------------STOCKAGE U0----------------------------------

   call save_file(x,u_prev,'u_init.txt' ) 
    
    


end subroutine finit_save



function f(u)
    
    real(rp) , intent(in):: u
    real(rp) :: f
   
   
    f=u**3 - u

  return 
   
end function f


function fp(x)

real(rp) , intent(in):: x
real(rp) :: fp
   
   
   fp=3._rp*x**2 - un

return
 
end function fp



function minmod(ug,uc,ud,dx)
  real(rp) , intent(in):: ug,uc,ud, dx
  real(rp) :: init, sr, sl, minmod
   
  sr=(uc-ud)/dx
  sl=(ug-uc)/dx
  
  if( sl>zero  .and. sr >zero)then 
  
      minmod=min(sl,sr)
  else if( sl<zero  .and. sr <zero) then
   
       minmod=max(sl,sr)
  else
  
   minmod=zero
  end if
   
return
end function minmod


subroutine LF(ul,ur,dx,dt,g) ! Flux LAX FRIEDRIECH



real(rp), intent(in) :: ul, ur, dx , dt
real(rp), intent(out) :: g


g=0.5*(f(ul)+f(ur))-0.5*(dx/dt)*(ur-ul)  ! Lax-Friedriech flow

return

end subroutine LF


subroutine Godu(ul,ur,i, g) !! Flux Godunov

real(rp), intent(in) :: ul, ur
integer , intent(in) :: i
real(rp), intent(out) :: g
real(rp) :: sigma, u_star



!-------------------------f convexe-------------------------------------

    
if(ul== ur)then

g=f(ul)

return
end if




if(ul > zero .and. ur > zero ) then 

    
    if(ul < ur)then
    
    

        if(i*dx < t*fp(ul)) then

            g=f(ul)
            
         else if(  t*fp(ul)< i*dx  .and. i*dx< t*fp(ur) ) then
         
            g=sqrt(((i*dx/t)+un)/3._rp)
              
        
        else if (  t*fp(ur) < i*dx )then
            
            g=f(ur)
        
        end if
        
    else if (ur < ul)then 
    
        sigma=ur**2+ur*ul+ul**2 -un
    
        if(i*dx< t*sigma)then
        
            g=f(ul)
            
        else if(t*sigma < i*dx)then
        
            g=f(ur)
        
        
        end if  
        
    end if
   
!-----------------------------f concave---------------------------------   
    
else if (ul < zero .and. ur < zero ) then 


    if(ul >  ur)then
       
        
    
    
        if(i*dx< t*fp(ul))then !detente
        
            g=f(ul)
            
        else if( t*fp(ul)< i*dx .and. i*dx < t*fp(ur) )then
        
            g=f(sqrt(((i*dx/t)+un)/3._rp))
        else
            
            g=f(ur)
        end if  
     
     
       
        
    else if (ul < ur)then  ! choc 
    
        sigma=ur**2+ur*ul+ul**2 -un
    
        if(i*dx< t*sigma)then

            g=f(ul)
        
        else if (t*sigma < i*dx) then
        
            g=f(ur)
        
        end if
        
    end if
  
  
  
    
!----------------------f ni concave ni convexe--------------------------    

else if (ul < zero .and. ur > zero ) then 

           
        
            u_star=-0.5_rp*ur 
            sigma=(f(u_star)-f(ul))/(u_star-ul)
       
            
            if( i*dx <  fp(ul)*t )then
            
                g=f(ul)
                
            else if( fp(ul)*t <  i*dx  .and.  i*dx <  t*sigma  )then
            
                g=f(u_star)
        
            else if (  sigma*t <  i*dx .and.  i*dx < fp(ur)*t ) then 
            
                g=f(sqrt(((i*dx/t)+un)/3._rp))
            
           
            else if  (  fp(ur)*t  <  i*dx) then
           
                g=f(ur)
            
            end if
        
        
       
            
            
            
        
else if(  zero  < ul   .and.  ur <  zero ) then
        
   
       
            
            u_star=-0.5_rp*ur
            
            sigma=(f(u_star)-f(ul))/(u_star-ul)
       
            if(i*dx< t*f(ul) )then 
        
                g=f(ul)
             else if(t*f(ul)< i*dx .and. i*dx < t*sigma)   then
             
				g=f(u_star)
                
            
            else if( t*sigma< i*dx  .and. i*dx<t*fp(ur) )then
        
                g=f(sqrt(((i*dx/t)+un)/3._rp))
            
            else if( t*fp(ur)< i*dx  )then
        
                g=f(ur)
                
            end if
            
 


end if




end subroutine Godu



SUBROUTINE   LF_solve()
    
    INTEGER :: i
    
    
    t=zero
    
    
     DO i=1,N
          
         u_prev(i)=finit(g+i*dx)    
    END DO
    
    
    
       
DO WHILE (t<Tmax)

   
!---------------Dt EN FONCTION DE LA CFL-----------------------------
    
    DO i=1,N
         tmp(i)=fp(u_prev(i))       
    END DO
    
    dt=0.5_rp*dx/abs(maxval(tmp))
    
!--------------------VECTEUR DE FLUX------------------------------------

    flux(1)=f(u_prev(1)) ! BORD GAUCHE 
    
    DO i=2, N
       
        CALL LF(u_prev(i-1),u_prev(i),dx,dt,flux(i)) ! FLUX NUMERIQUE 
   
    END DO
    
    !flux(N)=f(u_prev(N)) ! BORD DROIT
    
!---------------------MISE A JOUR---------------------------------------
    DO i=1, N-1
    
       u_next(i)=u_prev(i)-(dt/dx)*(flux(i+1)-flux(i))
       
         
    END DO
    
    
    u_next(N)=u_prev(N)-(dt/dx)*(flux(N)-flux(N-1))
    
    u_prev=u_next
    t=t+dt
       
END DO

!-------------------------Sauvegarde ----------------------------------   
    
    
    CALL save_file(x,u_next, 'u_approx_LF.txt') 
    
    
END SUBROUTINE   LF_solve



SUBROUTINE   Godu_solve()
    
    INTEGER :: k
    
    
    t=zero
    
     DO k=1,N
         u_prev(k)=0.5_rp*finit(g+k*dx)   
         
    END DO
    
    
DO WHILE (t<Tmax)
     
      !write(*,*) t
    
    
!---------------Dt EN FONCTION DE LA CFL-----------------------------
    
    DO k=1,N
          
         tmp(k)=abs(fp(u_prev(k)) )      
    END DO
       
         
    dt=0.5_rp*dx/maxval(tmp)
        
!--------------------VECTEUR DE FLUX------------------------------------

    
    
    ! BORD GAUCHE 
    flux(0)=f(u_prev(1))
    CALL Godu(u_prev(1),u_prev(2),k,flux(1))
    
    DO k=2, N-1
       
        CALL Godu(u_prev(k-1),u_prev(k),k,flux(k)) ! FLUX NUMERIQUE 
        
    END DO
    
    CALL Godu(u_prev(N-1),u_prev(N),N,flux(N)) ! BORD DROIT
    
!---------------------MISE A JOUR---------------------------------------
    
    !write(*,*) maxval(flux)
    u_next(1)=u_prev(1)-(dt/dx)*(flux(1)-flux(0))
    
    DO k=2, N
    
       u_next(k)=u_prev(k)-(dt/dx)*(flux(k)-flux(k-1))
       
         
    END DO
    
    
    u_prev=u_next
    t=t+dt
      
END DO

!-------------------------Sauvegarde ----------------------------------   
    
    
    CALL save_file(x,u_next, 'u_approx_Godunov.txt') 
    
    
END SUBROUTINE   Godu_solve


subroutine Muscl()  

   
    real(rp), dimension(:), allocatable  :: u, uplus, umoin, tmp
    real(rp) ::  sigma, gg , dd 
    integer :: i
    
    
    
    
    allocate(uplus(n))
    allocate(umoin(n))
    allocate(u(n))
    allocate(tmp(n))
    DO i=1,N
          
         u(i)=finit(g+i*dx)    
    END DO
    
  Do While(t<Tmax) 
  
  
  !---------------Dt EN FONCTION DE LA CFL-----------------------------
    
    DO i=1,N
          
         tmp(i)=fp(u(i))       
    END DO
    
    dt=0.5_rp*dx/abs(maxval(tmp)) 
!-----------------------------------------------------------------------    

   umoin(1)=u(1)-0.5_rp*dx*fp(u(1))
   uplus(1)=u(1)+0.5_rp*dx*fp(u(1))

    do i=2,n-1
        sigma=minmod(u(i-1),u(i),u(i+1),dx)
        umoin(i)=u(i)-0.5_rp*dx*sigma
        uplus(i)=u(i)+0.5_rp*dx*sigma
    end do  
    
    umoin(n)=u(n)-0.5_rp*dx*sigma
    uplus(n)=u(n)+0.5_rp*dx*sigma

!--------------BOUCLE EN ESPACE----------------------------------------
    
    
  
    do i=2,n-1

        CALL LF(uplus(i),umoin(i+1),dx,dt,gg)
        CALL LF(uplus(i-1),umoin(i),dx,dt,dd)
        
            u(i)=u(i)-(dt/dx)*(gg-dd)  
    end do 
    


   CALL LF(uplus(n),umoin(n),dx,dt,gg)
   CALL LF(uplus(n-1),umoin(n),dx,dt,dd)
   
   u(n)=u(n)-(dt/dx)*(gg-dd)
   
   
   t=t+dt
   
  END DO 
   
 
   Call save_file(x,u, 'Muscl.txt')
   
   
   
end subroutine



end module flow
