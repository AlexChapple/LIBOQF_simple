!
! Fortran version of Little Boxes Simple 
! Mostly for making the plots much quicker for my dissertation 
!
! THIS FILE IS MAKING THE  g2 PLOTS SPECIFICALLY

program main 

    IMPLICIT NONE 

    ! Declare general variables and parameters 
    REAL (KIND=8), PARAMETER :: dt = 0.05 
    INTEGER (KIND=8), PARAMETER :: segment_count = 200d0 
    INTEGER (KIND=8), PARAMETER :: iterations = 500000d0 
    REAL (KIND=8), PARAMETER :: time_interval = 10.0d0 
    REAL (KIND=8), PARAMETER :: Omega = 1.0d0 
    REAL (kind=8), PARAMETER :: pi = 3.14159265358979323846d0
    INTEGER (KIND=8) :: index, t
    REAL (KIND=8) :: total, rand_num

    REAL (KIND=8), DIMENSION(segment_count) :: rand_list, g2_list 

    COMPLEX (KIND=8), DIMENSION(2) :: coeffs
    COMPLEX (KIND=8) :: g_0_new, e_0_new, g_1_new, e_1_new, e_up, cur_e_up, g2

    ! Initialise 
    g_0_new = 0.0d0; g_1_new = 0.0d0; e_0_new = 0.0d0
    e_1_new = 0.0d0; e_up = 0.0d0; g2_list = 0.0d0
    rand_list = 0.0d0; g2 = 0.0d0; cur_e_up = 0.0d0 
    
    ! -------------------------------------------------
    !       Start running code 
    ! -------------------------------------------------

    do index = 1, iterations ! Iterates over time intervals 

        ! At the start of the time interval we force the system to go back 
        ! to its initial state 

        e_up = SQRT(modulo_func(e_0_new)**2 + & ! Keep track expectation from the last time
                    modulo_func(e_1_new)**2) 

        coeffs = 0.0d0; coeffs(1) = 1.0d0 
        g_0_new = 0.0d0; e_0_new = 0.0d0; g_1_new = 0.0d0; e_1_new = 0

        call RANDOM_NUMBER(rand_list)

        do t = 1, segment_count

            g_0_new = coeffs(1) - ((Omega/2)*coeffs(2)*dt)
            e_0_new = (coeffs(2)*(1-dt/2)) + (Omega/2)*coeffs(1)*dt
            g_1_new = CMPLX(0,-1) * SQRT(dt/2) * coeffs(2)
            e_1_new = CMPLX(0,-1) * SQRT(dt/2) * coeffs(2)

            total = SQRT(modulo_func(e_0_new)**2 + &
                    modulo_func(g_0_new)**2 + modulo_func(e_1_new)**2 &
                    + modulo_func(g_1_new)**2)

            g_0_new = g_0_new / total 
            e_0_new = e_0_new / total 
            g_1_new = g_1_new / total 
            e_1_new = e_1_new / total 

            rand_num = rand_list(t)

            if (rand_num <= 2 * modulo_func(g_1_new)**2) then 

                cur_e_up = SQRT(modulo_func(e_0_new)**2 + & 
                    modulo_func(e_1_new)**2)

                g2 = (e_up * cur_e_up) 
                g2_list(t) = g2_list(t) + modulo_func(g2)

                coeffs = 0.0d0 
                coeffs(1) = 1.0d0 

            else 

                total = SQRT(modulo_func(e_0_new)**2 + &
                    modulo_func(g_0_new)**2)    

                g_0_new = g_0_new / total 
                e_0_new = e_0_new / total 

                coeffs = 0.0d0 
                coeffs(1) = g_0_new
                coeffs(2) = e_0_new             


            end if 


        end do 

        IF (mod(index, 1000) == 0) THEN 
            PRINT *, index
        END IF 


    end do 

    g2_list = g2_list / (iterations / segment_count)

    open(1, file="results_g2/g2.txt", status="replace")

    do index = 1, segment_count
        write(1,*) g2_list(index)
    end do  

    close(1)


    contains 

    function modulo_func(z) result(c)

        ! Takes in a complex number z and returns its modulus

        implicit none

        ! Declare var types
        real (kind=8) :: a, b, c
        complex (kind=8), intent(in) :: z

        a = real(z)
        b = aimag(z)

        c = sqrt(a**2 + b**2)

    end function 








end program main 