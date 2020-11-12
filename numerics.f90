
 MODULE numerics
  
!-------------------------Parametres------------------------------------
    
      INTEGER,  PARAMETER  ::  rp = 8
      REAL(rp),  PARAMETER  ::  un   = 1.0_rp
      REAL(rp),  PARAMETER  ::  zero = 0.0_rp
      REAL(rp),  PARAMETER  ::   grv=9.81_rp
      
!---------------------------Scalaire------------------------------------
      REAL(rp), dimension(:), allocatable :: x, u0,tmp, u_prev, u_next, flux, up, um
      REAL(rp), dimension(:,:), allocatable :: uex
      REAL(rp) ::  d, g, dx, dt, Tmax, t
      INTEGER :: n
     
!-----------------------------------------------------------------------      
contains



subroutine init_space()


integer :: i

!-----------------------------------------------------------------------

 
    n=510 ! nombre points en espace
    g=-1._rp
    d=1._rp
    Tmax=2.5_rp  !temps maximal
    
    dx=(d-g)/(n-1)  ! pas d'espace
    
    t=0._rp
    dt=0._rp !pas de temps






!---------------------------Allocation----------------------------------
    
    allocate(x(n))
    allocate(u_prev(n))
    allocate(u_next(n)) 
    allocate(up(n))
    allocate(um(n))
    allocate(tmp(n))
    allocate(flux(0:n))

!---------------------INITIALISATION ESPACE----------------------------
    
    do i=1, N
    
       x(i)=g+i*dx
      
     
    end do
    



end subroutine init_space



SUBROUTINE save_file(x,u, fichier)

 real(rp), intent(in), dimension(:) :: u, x
 character(len=*), intent(in) :: fichier
! stockage dans un fichier
    
    open(unit=13, file=fichier, status='REPLACE')
    
    
        
        DO i=1, SIZE(X)
                
                WRITE(13,*) X(i), u(i)
        
      
        END DO
    
    close(13)
    


END SUBROUTINE 

END MODULE numerics
