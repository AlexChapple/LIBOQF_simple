!
! Fortran version of Little Boxes Simple 
! Mostly for making the plots much quicker for my dissertation 
!
! THIS FILE IS MAKING THE  g2 PLOTS SPECIFICALLY

program main 

    IMPLICIT NONE 

    ! Declare general variables and parameters 
    REAL (KIND=8), PARAMETER :: dt = 0.002 
    INTEGER (KIND=8), PARAMETER :: segment_count = 5000d0 
    INTEGER (KIND=8), PARAMETER :: time_list_size = 5000d0 
    INTEGER (KIND=8), PARAMETER :: iterations = 5000d0 
    REAL (KIND=8), PARAMETER :: time_interval = 10.0d0 
    REAL (KIND=8), PARAMETER :: Omega = 2.3d0 
    REAL (kind=8), PARAMETER :: pi = 3.14159265358979323846d0
    INTEGER (KIND=8) :: index, t
    REAL (KIND=8) :: total, total2, rand_num

    REAL (KIND=8), DIMENSION(time_list_size) :: rand_list
    REAL (KIND=8), DIMENSION(segment_count) :: g2_list 

    COMPLEX (KIND=8), DIMENSION(2) :: coeffs
    COMPLEX (KIND=8) :: g_0_new, e_0_new, g_1_new, e_1_new, e_up, cur_e, g2
    COMPLEX (KIND=8) :: cur_g, g_down

    INTEGER (KIND=8) :: emission_count 

    ! Initialise 
    g_0_new = 0.0d0; g_1_new = 0.0d0; e_0_new = 0.0d0
    e_1_new = 0.0d0; e_up = 0.0d0; g2_list = 0.0d0
    rand_list = 0.0d0; g2 = 0.0d0; cur_e = 0.0d0 
    cur_g = 0.0d0; g_down = 0.0d0; total2 = 0.0d0 

    emission_count = 0d0 
    
    ! -------------------------------------------------
    !       Start running code 
    ! -------------------------------------------------

    coeffs = 0.0d0; coeffs(1) = 1.0d0 

    do index = 1, iterations ! Iterates over time intervals 

        ! At the start of the time interval we force the system to go back 
        ! to its initial state 

        e_up = (modulo_func(coeffs(2))**2)

        coeffs = 0.0d0; coeffs(1) = 1.0d0 
        g_0_new = 0.0d0; e_0_new = 0.0d0; g_1_new = 0.0d0; e_1_new = 0.0d0

        call RANDOM_NUMBER(rand_list)

        do t = 1, time_list_size

            g_0_new = coeffs(1) - ((Omega/2)*coeffs(2)*dt)
            e_0_new = (coeffs(2)*(1-dt/2)) + (Omega/2)*coeffs(1)*dt
            g_1_new = CMPLX(0,-1) * SQRT(dt/2) * coeffs(2)
            e_1_new = CMPLX(0,-1) * SQRT(dt/2) * coeffs(2)

            total = SQRT(modulo_func(e_0_new)**2 + &
                    modulo_func(g_0_new)**2 + (modulo_func(e_1_new)**2) &
                    + (modulo_func(g_1_new)**2))

            g_0_new = g_0_new / total 
            e_0_new = e_0_new / total 
            g_1_new = g_1_new / total 
            e_1_new = e_1_new / total 

            rand_num = rand_list(t)

            if (rand_num <= 2 * modulo_func(g_1_new)**2) then 

                cur_e = (modulo_func(e_0_new)**2 + 2*(modulo_func(e_1_new)**2))
                ! cur_e = (modulo_func(e_0_new)**2)

                g2 = (e_up * cur_e)
                g2_list(t) = g2_list(t) + (g2)

                coeffs = 0.0d0 
                coeffs(1) = 1.0d0 

                emission_count = emission_count + 1 

            else 

                cur_e = (modulo_func(e_0_new)**2 + 2*(modulo_func(e_1_new)**2))
                ! cur_e = (modulo_func(e_0_new)**2)

                g2 = (e_up * cur_e)
                g2_list(t) = g2_list(t) + (g2)

                total = SQRT(modulo_func(e_0_new)**2 + &
                    modulo_func(g_0_new)**2)    

                g_0_new = g_0_new / total 
                e_0_new = e_0_new / total 

                coeffs = 0.0d0 
                coeffs(1) = g_0_new
                coeffs(2) = e_0_new             


            end if 


        end do 

        ! cur_e = SQRT(modulo_func(e_0_new)**2 + & 
        !             modulo_func(e_1_new)**2)
        ! g2_list(index) = g2_list(index) + (e_up * cur_e)

        IF (mod(index, 1000) == 0) THEN 
            PRINT *, index
        END IF 


    end do 

    print *, emission_count

    g2_list = g2_list / iterations

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