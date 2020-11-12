    PROGRAM TRANSPORT
    use numerics
    use flow
    use rusanov
    use musc
    use omp_lib
    
    IMPLICIT NONE
    
   

    Call init_space()
    Call init_syst(101,.5_rp)
    

    
    !CALL Muscl()
    !CALL LF_solve()
    !CALL Godu_Solve()
    !CALL HLL_Solve()  ! OK
    CALL RUS_Solve()  ! FONCTIONNE BIEN 
    
    !CALL MUSCL_Solve() 
    
    END PROGRAM TRANSPORT
    
    
