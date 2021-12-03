! 
! Revised version of the Little Boxes revised 
!

program main 

    implicit none 

    ! Declare variables 
    real (kind=8), parameter :: start_time = 0.0d0 
    complex (kind=8), dimension(6) :: coeffs 
    integer, parameter :: end_time = 60d0 
    integer (kind=8), parameter :: time_steps = 60000d0 
    real (kind=8), parameter :: dt = 0.01d0 
    real (kind=8), dimension(time_steps) :: time_list
    real (kind=8), dimension(time_steps) :: sigma_z_list, sigma_L_list, sigma_R_list 
    integer (kind=8), parameter :: num_of_simulations = 10000d0
    complex (kind=8), parameter :: Omega = 2.3d0 
    real (kind=8), parameter :: pi = 3.14159265358979323846d0 
    integer (kind = 8), parameter :: period = 5d0 
    real (kind=8) :: Delta_t , rand_num, prob, total 
    real (kind=8), dimension(time_steps) :: rand_list 

    integer (kind=8) :: sim, t, index, index2 
    integer (kind=8) :: beginning, ended_time, rate 

    complex (kind=8) :: g_0, e_0, g_L, e_L, g_R, e_R 

    ! Photon emission tracking 
    integer (kind=8) :: emission_index, total_emisison_counter 
    integer (kind=8), parameter :: tracking_bin_width = 1000
    real (kind=8), dimension(num_of_simulations, tracking_bin_width) :: emission_tracking_list



    Delta_t = period * dt
    
    sigma_z_list = 0.0d0; sigma_L_list = 0.0d0; sigma_R_list = 0.0d0 
    time_list = 0.0d0 
    rand_list = 0.0d0 
    coeffs = 0.0d0 
    
    ! Experimental 
    total_emisison_counter = 0d0 

    ! ---------------------------------------------------
    !           Start running code 
    ! ---------------------------------------------------

    ! Program execution time tracking 
    call system_clock(beginning, rate)

    call linspace(start=start_time, end=end_time, list=time_list)

    do sim = 1, num_of_simulations

        ! Initialise the random number list 
        call random_number(rand_list)

        coeffs = 0.0d0; coeffs(1) = 1.0d0 
        g_0 = 0.0d0; e_0 = 0.0d0; g_L = 0.0d0; e_L = 0.0d0; g_R = 0.0d0; e_R = 0.0d0 
        emission_index = 1d0 

        do t = 1, size(time_list)

            ! Update the coefficients one time step 
            ! g_0 = coeffs(1) + ((Omega/2) * coeffs(2) * dt) 

            ! e_0 = coeffs(2) + ((((-Omega/2) * coeffs(1)) &
            !         + (cmplx(0,-1) * sqrt(1/(2 * Delta_t)) * coeffs(3)) &
            !         + (cmplx(0,-1) * sqrt(1/(2 * Delta_t)) * coeffs(5))) * dt)

            ! g_L = coeffs(3) + (((-Omega/2)*coeffs(4) + cmplx(0,-1)*sqrt(1/(2 * Delta_t)) * coeffs(2))*dt)

            ! e_L = coeffs(4) + (((Omega/2)*coeffs(3)) * dt)

            ! g_R = coeffs(5) + (((-Omega/2)*coeffs(6) + cmplx(0,-1)*sqrt(1/(2 * Delta_t)) * coeffs(2))*dt)

            ! e_R = coeffs(6) + (((Omega/2)*coeffs(5)) * dt)

            g_0 = coeffs(1) - ((Omega/2)*coeffs(2)*dt)
            e_0 = (coeffs(2)*(1-dt/2)) + (Omega/2)*coeffs(1)*dt
            g_L = cmplx(0,-1) * sqrt(dt/2) * coeffs(2)
            e_R = cmplx(0,-1) * sqrt(dt/2) * coeffs(2)

            ! Normalise 
            total = sqrt(modulo_func(g_0)**2 + modulo_func(e_0)**2 &
                    + 2*modulo_func(g_L)**2 + 2*modulo_func(e_L)**2) 

            g_0 = g_0 / total; e_0 = e_0 / total 
            g_L = g_L / total; e_L = e_L / total 
            g_R = g_R / total; e_R = e_R / total

            if (mod(t,period) == 0) then ! Check the boxes now 

                rand_num = rand_list(t)

                prob = 2 * modulo_func(g_L)**2

                if (rand_num <= prob) then 

                    ! A photon has been emitted

                    coeffs = 0.0d0 
                    coeffs(1) = 1.0d0 

                    sigma_z_list(t) = sigma_z_list(t) + modulo_func(coeffs(2))**2 - modulo_func(coeffs(1))**2 
                    sigma_L_list(t) = sigma_L_list(t) - (conjg(coeffs(1)) * coeffs(2)) ! negative sign from the phase shift of e^{-i \pi} 
                    sigma_R_list(t) = sigma_R_list(t) + (conjg(coeffs(2)) * coeffs(1))


                    ! Store the emission data 
                    emission_tracking_list(sim, emission_index) = time_list(t)
                    emission_index = emission_index + 1
                    total_emisison_counter = total_emisison_counter + 1 

                else 

                    coeffs = 0.0d0 

                    total = sqrt(modulo_func(g_0)**2 + modulo_func(e_0)**2)

                    g_0 = g_0 / total 
                    e_0 = e_0 / total 

                    coeffs(1) = g_0; coeffs(2) = e_0 

                    sigma_z_list(t) = sigma_z_list(t) + modulo_func(coeffs(2))**2 - modulo_func(coeffs(1))**2 
                    sigma_L_list(t) = sigma_L_list(t) - (conjg(coeffs(1)) * coeffs(2))
                    sigma_R_list(t) = sigma_R_list(t) + (conjg(coeffs(2)) * coeffs(1))



                end if

            else 

                coeffs(1) = g_0; coeffs(2) = e_0
                coeffs(3) = 0.0d0; coeffs(4) = 0.0d0
                coeffs(5) = 0.0d0; coeffs(6) = 0.0d0 

                sigma_z_list(t) = sigma_z_list(t) + modulo_func(coeffs(2))**2 - modulo_func(coeffs(1))**2 
                sigma_L_list(t) = sigma_L_list(t) - (conjg(coeffs(1)) * coeffs(2))
                sigma_R_list(t) = sigma_R_list(t) + (conjg(coeffs(2)) * coeffs(1))
                

            end if 




        end do 

        if (mod(sim,1000) == 0) then 
            print *, sim 
        end if 


    end do 

    print *, total_emisison_counter

    sigma_z_list = sigma_z_list / num_of_simulations
    sigma_L_list = sigma_L_list / num_of_simulations
    sigma_R_list = sigma_R_list / num_of_simulations

    open(1, file="testing/sigma_z_test.txt", status="replace")
    open(2, file="testing/sigma_L_test.txt", status="replace")
    open(3, file="testing/sigma_R_test.txt", status="replace")
    open(10, file="testing/emission_tracking.txt", status="replace")

    do index = 1, size(time_list)

        write(1,*) time_list(index), sigma_z_list(index)
        write(2,*) time_list(index), sigma_L_list(index)
        write(3,*) time_list(index), sigma_R_list(index)

    end do 

    do index = 1, num_of_simulations 
        do index2 = 1, tracking_bin_width
            if (emission_tracking_list(index, index2) == 0d0) then 
                exit
            else 
                write(10,*) emission_tracking_list(index, index2)
            end if 
        end do 

        write(10,*) end_time + 50d0 

    end do 

    close(1); close(2); close(3); close(10)


    call system_clock(ended_time)

    print *, "All simulations completed. Execution time: ", real(ended_time - beginning) / real(rate), " seconds."


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

    subroutine linspace(start, end, list) 

        real (kind=8), intent(in) :: start
        integer, intent(in) :: end 
        real (kind=8), intent(out) :: list(:)
        real (kind=8) :: range
        integer :: n, i

        n = size(time_list)
        range = end - start
        
        do i = 1,n
            time_list(i) = start + (range * (i - 1) / (n - 1))
        end do

    end subroutine






end program main 