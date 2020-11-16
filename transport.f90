    PROGRAM TRANSPORT
    use numerics
    use rusanov
    use musc
    use omp_lib
    
    IMPLICIT NONE
    
   

    Call init_space()
    Call init_syst(100,.5_rp)
    

    
    CALL HLL_Solve()  ! OK
    CALL RUS_Solve()  ! FONCTIONNE BIEN 
    CALL MUSCL_Solve() 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    END PROGRAM TRANSPORT
    
    
