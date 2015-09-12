PROGRAM RM_effect_mc_V05

!   Purpose: To determine the best values for the spin-orbit alignment and rotational velocity of the host star from a 
!        given data set of the Radial velocities during the RM effect. This program accepts input data from IDL
!        (see input data), does all the theoretical model calculations, and outputs data back to IDL for plotting
!        (see output data). This code was written in FORTRAN to increase speed and computational efficiency. The best
!        fit Keplarian is determined through a chi squared analysis of the parameter grid space using specified step
!        sizes determined by the user. While this may not be the most efficient method for calculating vsini and
!        lambda from the RV data, it has the benefit of calculating the chi squared values over the whole
!        parameter space and avoiding local minimums.  The parameters which the model iterates over
!        (i.e. parameters the user seeks) includes the rotational velocity of the host star (vsini) and the spin-
!        orbit alignment (lambda). Note: This model has been updated to calculate the precise location of the 
!            planet in it's orbit at specific data points by solving Kepler's equation 
!            (E - e*sin(E) = 2*pi*((time - JD_time_peri)/period)) for the eccentric anomaly using a Bessel function.
!            Previous calculations for the planet's orbital location were done through a numeric method of calculating
!            the planet's position at small time steps and moving the planet a distance determined through the orbital
!            angular velocity at each step. This method isn't as accurate or as quick as using a more analytical
!            method with Bessel functions.


!   Input parameters: omega_arg_periastron, Ecc 
!            Mp, Orbital_period 
!            JD_time_peri_interval, JD_time_peri_begin, JD_time_peri_end, Ms, Number_fit, long_periapsis, 
!                     Number_iterations, datafilelength, True_anomaly_start, True_anomaly, Number_orbits, Data_array, input_RV_filename,
!                     directory_location, Planet_name, output_temp_data, Jupiter_Earth_units, Phase_shift,

!   Output parameters: Chi_squared_array, Reduced_Chi_Squared_array, Parameter_array, omega_arg_periastron_array, Ecc_array, 
!                      Mp_array, Rorb_array, JD_time_peri_array, Number_omega_arg_periastron, Number_Ecc, Number_Mp,
!                      Number_Rorb, Number_JD_time_peri
use ran_mod
implicit none

!   Constant and variable decleration

double precision, parameter :: Rss = 6.9634D8                !Radius of our Sun (in meters).
double precision, parameter :: Rj = 7.1492D7                 !Radius of Jupiter (in meters).
double precision, parameter :: RE = 6.3781D6                 !Radius of Earth (in meters if planet is given in Earth radii).
double precision, parameter :: AU = 1.4960D11                !One astronomical unit in meters.
double precision, parameter :: pi = 3.1415926535           !Constant value of pi.
double precision, parameter :: G = 6.67384D-11                !Gravitation constant.
double precision, parameter :: Mss = 1.9891D30               !Mass of the sun in kilograms.
double precision, parameter :: Me = 5.9722D24                !Mass of the Earth in kg.
double precision, parameter :: Mj = 1.8986D27                !Mass of Jupiter in kg.
double precision, parameter :: day_sec = 86400.0D0         !Length of a day in seconds.

double precision :: a, b, c         !variables used to read in data
double precision :: Albedo
double precision :: Alpha1
double precision :: another_min
double precision :: Aplan
double precision :: Area_pixel
double precision :: average_delta_chi
double precision :: Bessel_function_exit
double precision :: Bessel_value
double precision :: best_Ecc_fit
double precision :: best_Ecc_total
double precision :: best_Ecc_nettotal
double precision :: best_Inc_fit
double precision :: best_Inc_total
double precision :: best_Inc_nettotal
double precision :: best_JD_time_mid_transit_fit
double precision :: best_JD_time_mid_transit_total
double precision :: best_JD_time_mid_transit_nettotal
double precision :: best_Mp_fit
double precision :: best_Mp_total
double precision :: best_Mp_nettotal
double precision :: best_omega_arg_periastron_fit
double precision :: best_omega_arg_periastron_total
double precision :: best_omega_arg_periastron_nettotal
double precision :: best_orbital_period_fit
double precision :: best_orbital_period_total
double precision :: best_orbital_period_nettotal
double precision :: best_Rp_fit
double precision :: best_Rp_total
double precision :: best_Rp_nettotal
double precision :: best_RV_offset_datasets_total
double precision :: best_RV_offset_datasets_nettotal
double precision :: best_RV_offset_datasets_fit
double precision :: best_RV_zero_offset_fit
double precision :: best_RV_zero_offset_total
double precision :: best_RV_zero_offset_nettotal
double precision :: best_spin_orbit_angle
double precision :: best_spin_orbit_fit
double precision :: best_spin_orbit_fit_prior
double precision :: best_spin_orbit_total
double precision :: best_spin_orbit_nettotal
double precision :: best_vsini
double precision :: best_vsini_fit
double precision :: best_vsini_fit_prior
double precision :: best_vsini_total
double precision :: best_vsini_nettotal
double precision :: Beta1
double precision :: bininterval
double precision :: Center_of_planet_x
double precision :: Center_of_planet_y
double precision :: chi_squared_change
double precision :: chi_squared_change_fit
double precision :: chi_squared_change_mc
double precision :: chi_square_Ecc_minus_fit
double precision :: chi_square_Ecc_minus_total
double precision :: chi_square_Ecc_minus_nettotal
double precision :: chi_square_Ecc_plus_fit
double precision :: chi_square_Ecc_plus_total
double precision :: chi_square_Ecc_plus_nettotal
double precision :: chi_square_Inc_minus_fit
double precision :: chi_square_Inc_minus_total
double precision :: chi_square_Inc_minus_nettotal
double precision :: chi_square_Inc_plus_fit
double precision :: chi_square_Inc_plus_total
double precision :: chi_square_Inc_plus_nettotal
double precision :: chi_square_JD_time_mid_transit_minus_fit
double precision :: chi_square_JD_time_mid_transit_minus_total
double precision :: chi_square_JD_time_mid_transit_minus_nettotal
double precision :: chi_square_JD_time_mid_transit_plus_fit
double precision :: chi_square_JD_time_mid_transit_plus_total
double precision :: chi_square_JD_time_mid_transit_plus_nettotal
double precision :: chi_square_Mp_minus_fit
double precision :: chi_square_Mp_minus_total
double precision :: chi_square_Mp_minus_nettotal
double precision :: chi_square_Mp_plus_fit
double precision :: chi_square_Mp_plus_total
double precision :: chi_square_Mp_plus_nettotal
double precision :: chi_square_omega_arg_periastron_minus_fit
double precision :: chi_square_omega_arg_periastron_minus_total
double precision :: chi_square_omega_arg_periastron_minus_nettotal
double precision :: chi_square_omega_arg_periastron_plus_fit
double precision :: chi_square_omega_arg_periastron_plus_total
double precision :: chi_square_omega_arg_periastron_plus_nettotal
double precision :: chi_square_orbital_period_minus_fit
double precision :: chi_square_orbital_period_minus_total
double precision :: chi_square_orbital_period_minus_nettotal
double precision :: chi_square_orbital_period_plus_fit
double precision :: chi_square_orbital_period_plus_total
double precision :: chi_square_orbital_period_plus_nettotal
double precision :: chi_square_Rp_minus_fit
double precision :: chi_square_Rp_minus_total
double precision :: chi_square_Rp_minus_nettotal
double precision :: chi_square_Rp_plus_fit
double precision :: chi_square_Rp_plus_total
double precision :: chi_square_Rp_plus_nettotal
double precision :: chi_square_RV_offset_datasets_minus_fit
double precision :: chi_square_RV_offset_datasets_minus_total
double precision :: chi_square_RV_offset_datasets_minus_nettotal
double precision :: chi_square_RV_offset_datasets_plus_fit
double precision :: chi_square_RV_offset_datasets_plus_total
double precision :: chi_square_RV_offset_datasets_plus_nettotal
double precision :: chi_square_RV_zero_offset_minus_fit
double precision :: chi_square_RV_zero_offset_minus_total
double precision :: chi_square_RV_zero_offset_minus_nettotal
double precision :: chi_square_RV_zero_offset_plus_fit
double precision :: chi_square_RV_zero_offset_plus_total
double precision :: chi_square_RV_zero_offset_plus_nettotal
double precision :: chi_square_spin_orbit_minus
double precision :: chi_square_spin_orbit_minus_fit
double precision :: chi_square_spin_orbit_minus_total
double precision :: chi_square_spin_orbit_minus_mc
double precision :: chi_square_spin_orbit_minus_nettotal
double precision :: chi_square_spin_orbit_plus
double precision :: chi_square_spin_orbit_plus_fit
double precision :: chi_square_spin_orbit_plus_total
double precision :: chi_square_spin_orbit_plus_mc
double precision :: chi_square_spin_orbit_plus_nettotal
double precision :: chi_square_vsini_minus
double precision :: chi_square_vsini_minus_total
double precision :: chi_square_vsini_minus_mc
double precision :: chi_square_vsini_minus_nettotal
double precision :: chi_square_vsini_minus_fit
double precision :: chi_square_vsini_plus
double precision :: chi_square_vsini_plus_total
double precision :: chi_square_vsini_plus_mc
double precision :: chi_square_vsini_plus_nettotal
double precision :: chi_square_vsini_plus_fit
double precision :: chi_2
double precision :: chi2
double precision :: constant
double precision :: cpu_time_1
double precision :: cpu_time_2
double precision :: cpu_time_3
double precision :: cpu_time_4
double precision :: cpu_time_it1
double precision :: cpu_time_it2
double precision :: cpu_time_it3
double precision :: cpu_time_it4
double precision :: data_plot_model_interval
double precision :: delta_chi_square
double precision :: Dist2
double precision :: Distance_center
double precision :: Dist_center_pixel
double precision :: dist_cent1_int
double precision :: Dist_planet_pixel
double precision :: Ecc
double precision :: Ecc_prior
double precision :: Ecc_end
double precision :: Ecc_begin
double precision :: Ecc_interval
double precision :: Ecc_1sigerr
double precision :: ecc_anomaly
double precision :: ecc_anomaly_after
double precision :: ecc_anomaly_before
double precision :: ecc_anomaly_start
double precision :: ecc_anomaly_transit
double precision :: Ecc_minus_error_fit
double precision :: Ecc_minus_error_total = 0.0D0
double precision :: Ecc_minus_error_nettotal = 0.0D0
double precision :: Ecc_plus_error_fit
double precision :: Ecc_plus_error_total = 0.0D0
double precision :: Ecc_plus_error_nettotal = 0.0D0
double precision :: fdpo
double precision :: fdpp
double precision :: impact_prior
double precision :: Inc_prior
double precision :: Inc_end
double precision :: Inc_begin
double precision :: Inc_interval
double precision :: Inc_1sigerr
double precision :: Inc
double precision :: Inc_minus_error_fit
double precision :: Inc_minus_error_total = 0.0D0
double precision :: Inc_minus_error_nettotal = 0.0D0
double precision :: Inc_plus_error_fit
double precision :: Inc_plus_error_total = 0.0D0
double precision :: Inc_plus_error_nettotal = 0.0D0
double precision :: Io
double precision :: Io_Pixel
double precision :: Io_planet
double precision :: JD_time_peri
double precision :: JD_time_peri_model
double precision :: JD_time_mid_transit
double precision :: JD_time_mid_transit_1
double precision :: JD_time_mid_transit_prior
double precision :: JD_time_mid_transit_end
double precision :: JD_time_mid_transit_begin
double precision :: JD_time_mid_transit_interval
double precision :: JD_time_mid_transit_1sigerr
double precision :: JD_time_mid_transit_minus_error_total = 0.0D0
double precision :: JD_time_mid_transit_minus_error_nettotal = 0.0D0
double precision :: JD_time_mid_transit_minus_error_fit
double precision :: JD_time_mid_transit_plus_error_total = 0.0D0
double precision :: JD_time_mid_transit_plus_error_nettotal = 0.0D0
double precision :: JD_time_mid_transit_plus_error_fit
double precision :: L_total_star
double precision :: last_model_time
double precision :: Lblocked
double precision :: Lblocked2
double precision :: Ldpo
double precision :: Ldpp
double precision :: Length
double precision :: length_time_compare
double precision :: long_periapsis
double precision :: Mean_anomaly_transit
double precision :: min_chi_squared
double precision :: min_chi_squared_total
double precision :: min_chi_squared_nettotal
double precision :: min_chi_squared_fit
double precision :: min_chi_squared_fit_prior
double precision :: min_reduced_chi_squared
double precision :: min_reduced_chi_squared_total
double precision :: min_reduced_chi_squared_nettotal
double precision :: min_reduced_chi_squared_fit
double precision :: min_reduced_chi_squared_fit_prior
double precision :: model_data_difference
double precision :: model_data_difference1
double precision :: model_data_difference_total2
double precision :: Mp
double precision :: Mpp
double precision :: Mp_prior
double precision :: Mp_end
double precision :: Mp_begin
double precision :: Mp_interval
double precision :: Mp_1sigerr
double precision :: Mp_minus_error_fit
double precision :: Mp_minus_error_total = 0.0D0
double precision :: Mp_minus_error_nettotal = 0.0D0
double precision :: Mp_plus_error_fit
double precision :: Mp_plus_error_total = 0.0D0
double precision :: Mp_plus_error_nettotal = 0.0D0
double precision :: Ms
double precision :: Ms_solar
double precision :: Number_orbits
double precision :: net_chi_squared_change
double precision :: omega_arg_periastron
double precision :: omega_arg_periastron_prior
double precision :: omega_arg_periastron_end
double precision :: omega_arg_periastron_begin
double precision :: omega_arg_periastron_interval
double precision :: omega_arg_periastron_1sigerr
double precision :: omega_arg_periastron_minus_error_fit
double precision :: omega_arg_periastron_minus_error_total = 0.0D0
double precision :: omega_arg_periastron_minus_error_nettotal = 0.0D0
double precision :: omega_arg_periastron_plus_error_fit
double precision :: omega_arg_periastron_plus_error_total = 0.0D0
double precision :: omega_arg_periastron_plus_error_nettotal = 0.0D0
double precision :: Orbital_period
double precision :: Orbital_period_1
double precision :: Orbital_period_2
double precision :: Orbital_period_prior
double precision :: Orbital_period_end
double precision :: Orbital_period_begin
double precision :: Orbital_period_interval
double precision :: Orbital_period_1sigerr
double precision :: orbital_period_minus_error_total = 0.0D0
double precision :: orbital_period_minus_error_nettotal = 0.0D0
double precision :: orbital_period_minus_error_fit
double precision :: orbital_period_plus_error_total = 0.0D0
double precision :: orbital_period_plus_error_nettotal = 0.0D0
double precision :: orbital_period_plus_error_fit
double precision :: Phase_angle
double precision :: Phase_angle_observed
double precision :: Phase_angle_start
double precision :: Phase_orbit_n
double precision :: phase_time
double precision :: Pixel
double precision :: Planet_star_distance
double precision :: q_1
double precision :: q_2
double precision :: r_chi
double precision :: r_Chi_2
double precision :: r_transit
double precision :: Radius_planet_array
double precision :: radius_ratio_prior
double precision :: Rorb
double precision :: Rorb_prior
double precision :: Rorb_star
double precision :: Rorb_star_prior
double precision :: Rp
double precision :: Rpp
double precision :: Rp_prior
double precision :: Rp_end
double precision :: Rp_begin
double precision :: Rp_interval
double precision :: Rp_1sigerr
double precision :: Rp2
double precision :: Rp_minus_error_fit
double precision :: Rp_minus_error_total = 0.0D0
double precision :: Rp_minus_error_nettotal = 0.0D0
double precision :: Rp_plus_error_fit
double precision :: Rp_plus_error_total = 0.0D0
double precision :: Rp_plus_error_nettotal = 0.0D0
double precision :: Rs
double precision :: Rs_solar
double precision :: Rs2
double precision :: RV
double precision :: RVamplitude
double precision :: RV_offset_datasets
double precision :: RV_offset_datasets_begin
double precision :: RV_offset_datasets_end
double precision :: RV_offset_datasets_interval
double precision :: RV_offset_datasets_prior
double precision :: RV_offset_datasets_1sigerr
double precision :: RV_offset_datasets_minus_error_total = 0.0D0
double precision :: RV_offset_datasets_minus_error_nettotal = 0.0D0
double precision :: RV_offset_datasets_minus_error_fit
double precision :: RV_offset_datasets_plus_error_total = 0.0D0
double precision :: RV_offset_datasets_plus_error_nettotal = 0.0D0
double precision :: RV_offset_datasets_plus_error_fit
double precision :: RV_zero_offset
double precision :: RV_zero_offset_prior
double precision :: RV_zero_offset_end
double precision :: RV_zero_offset_begin
double precision :: RV_zero_offset_interval
double precision :: RV_zero_offset_1sigerr
double precision :: RV_zero_offset_minus_error_fit
double precision :: RV_zero_offset_minus_error_total = 0.0D0
double precision :: RV_zero_offset_minus_error_nettotal = 0.0D0
double precision :: RV_zero_offset_plus_error_fit
double precision :: RV_zero_offset_plus_error_total = 0.0D0
double precision :: RV_zero_offset_plus_error_nettotal = 0.0D0
double precision :: set_distance_center
double precision :: slope_min
double precision :: spin_orbit_delta_chi
double precision :: spin_orbit_error_mc
double precision :: spin_orbit_minus_error
double precision :: spin_orbit_minus_error_total = 0.0D0
double precision :: spin_orbit_minus_error_mc = 0.0D0
double precision :: spin_orbit_minus_error_mc1 = 0.0D0
double precision :: spin_orbit_minus_error_mc2 = 0.0D0
double precision :: spin_orbit_minus_error_nettotal = 0.0D0
double precision :: spin_orbit_minus_error_fit
double precision :: spin_orbit_plus_error
double precision :: spin_orbit_plus_error_total = 0.0D0
double precision :: spin_orbit_plus_error_mc = 0.0D0
double precision :: spin_orbit_plus_error_mc1 = 0.0D0
double precision :: spin_orbit_plus_error_mc2 = 0.0D0
double precision :: spin_orbit_plus_error_nettotal = 0.0D0
double precision :: spin_orbit_plus_error_fit
double precision :: stellar_rotation_angle_begin
double precision :: stellar_rotation_angle_end
double precision :: stellar_rotation_angle_interval
double precision :: stellar_rotation_angle_prior
double precision :: stellar_rotation_angle_1sigerr
double precision :: stellar_rotation_angle
double precision :: Sub_planet_velocity
double precision :: sum_ecc_anomaly
double precision :: Transit_length_prior
double precision :: Time
double precision :: Time_bef_aft
double precision :: Time_check
double precision :: Time_compare_vel
double precision :: Time_mid_occultation
double precision :: Time_mid_transit
double precision :: Time_occultation_end
double precision :: Time_occultation_start
double precision :: time_peri_passage
double precision :: Time_ref
double precision :: time_remaining = 0.0D0
double precision :: Time_start
double precision :: Transit_start_no_inc
double precision :: Time_transit
double precision :: Time_transit_end
double precision :: Time_transit_start
double precision :: Theta_inside
double precision :: Transit_end_no_inc
double precision :: True_anomaly
double precision :: True_anomaly_start
double precision :: True_anomaly_transit
double precision :: True_phase
double precision :: total_chi_change
double precision :: Total_L
double precision :: Total_RM
double precision :: u
double precision :: v_rm
double precision :: variance_spin_orbit
double precision :: variance_spin_orbit_chi_squared
double precision :: variance_vsini
double precision :: variance_vsini_chi_squared
double precision :: vmacro
double precision :: vmicro
double precision :: vsini
double precision :: vsini_begin
double precision :: vsini_delta_chi
double precision :: vsini_end
double precision :: vsini_interval
double precision :: vsini_prior
double precision :: vsini_error_mc
double precision :: vsini_1sigerr
double precision :: vsini_minus_error
double precision :: vsini_minus_error_total = 0.0D0
double precision :: vsini_minus_error_mc = 0.0D0
double precision :: vsini_minus_error_mc1 = 0.0D0
double precision :: vsini_minus_error_mc2 = 0.0D0
double precision :: vsini_minus_error_nettotal = 0.0D0
double precision :: vsini_minus_error_fit
double precision :: vsini_plus_error
double precision :: vsini_plus_error_total = 0.0D0
double precision :: vsini_plus_error_mc = 0.0D0
double precision :: vsini_plus_error_mc1 = 0.0D0
double precision :: vsini_plus_error_mc2 = 0.0D0
double precision :: vsini_plus_error_nettotal = 0.0D0
double precision :: vsini_plus_error_fit
double precision :: vturb
double precision :: X_int_1
double precision :: X_int_2
double precision :: X_pixel
double precision :: X_pixel2
double precision :: X_pixel_prime
double precision :: X_prime
double precision :: x_prime_distance
double precision :: X_start_position
double precision :: X_transit
double precision :: Xpos
double precision :: Xpos_in
double precision :: Xpos_inside
double precision :: XpXp
double precision :: Y_int_1
double precision :: Y_int_2
double precision :: Y_pixel
double precision :: Y_planet_transit
double precision :: Y_start_position
double precision :: Ypos
double precision :: Ypos2
double precision :: Ypos_inside
double precision :: Z_start_position
double precision :: Zpos

DOUBLE PRECISION :: BESSJ
integer(4) :: counter_error_loop = 0
integer(4) :: counter_chi_loop
integer(4) :: counter_chi2_loop
integer(4) :: flag = 0
integer(4) :: Iteration_number = 0
integer(4), dimension(1) :: loc_min_chi_squared_total, loc_min_chi_squared_fit, loc_min_reduced_chi_squared_total, loc_min_reduced_chi_squared_fit
integer(4), dimension(1) :: loc_vsini_plus_error_total, loc_vsini_minus_error_total, loc_spin_orbit_plus_error_total, loc_spin_orbit_minus_error_total
integer(4), dimension(1) :: loc_RV_offset_datasets_plus_error_total, loc_RV_offset_datasets_minus_error_total, loc_orbital_period_plus_error_total
integer(4), dimension(1) :: loc_orbital_period_minus_error_total, loc_JD_time_mid_transit_plus_error_total, loc_JD_time_mid_transit_minus_error_total
integer(4), dimension(1) :: loc_Mp_plus_error_total, loc_Mp_minus_error_total, loc_Rp_plus_error_total, loc_Rp_minus_error_total, loc_Ecc_plus_error_total
integer(4), dimension(1) :: loc_Ecc_minus_error_total, loc_Inc_plus_error_total, loc_Inc_minus_error_total, loc_omega_arg_periastron_plus_error_total
integer(4), dimension(1) :: loc_omega_arg_periastron_minus_error_total, loc_RV_zero_offset_plus_error_total, loc_RV_zero_offset_minus_error_total
integer(4), dimension(1) :: loc_min_chi_squared_fit_prior, loc_min_reduced_chi_squared_fit_prior

integer(4), dimension(1) :: loc_vsini_plus_error_mc, loc_vsini_minus_error_mc, loc_spin_orbit_plus_error_mc, loc_spin_orbit_minus_error_mc

integer(4), dimension(1) :: loc_vsini_plus_error_nettotal, loc_vsini_minus_error_nettotal, loc_spin_orbit_plus_error_nettotal, loc_spin_orbit_minus_error_nettotal
integer(4), dimension(1) :: loc_RV_offset_datasets_plus_error_nettotal, loc_RV_offset_datasets_minus_error_nettotal, loc_orbital_period_plus_error_nettotal
integer(4), dimension(1) :: loc_orbital_period_minus_error_nettotal, loc_JD_time_mid_transit_plus_error_nettotal, loc_JD_time_mid_transit_minus_error_nettotal
integer(4), dimension(1) :: loc_Mp_plus_error_nettotal, loc_Mp_minus_error_nettotal, loc_Rp_plus_error_nettotal, loc_Rp_minus_error_nettotal, loc_Ecc_plus_error_nettotal
integer(4), dimension(1) :: loc_Ecc_minus_error_nettotal, loc_Inc_plus_error_nettotal, loc_Inc_minus_error_nettotal, loc_omega_arg_periastron_plus_error_nettotal
integer(4), dimension(1) :: loc_omega_arg_periastron_minus_error_nettotal, loc_RV_zero_offset_plus_error_nettotal, loc_RV_zero_offset_minus_error_nettotal

integer(4), dimension(1) :: loc_vsini_plus_error_fit, loc_vsini_minus_error_fit, loc_spin_orbit_plus_error_fit, loc_spin_orbit_minus_error_fit
integer(4), dimension(1) :: loc_RV_offset_datasets_plus_error_fit, loc_RV_offset_datasets_minus_error_fit, loc_orbital_period_plus_error_fit
integer(4), dimension(1) :: loc_orbital_period_minus_error_fit, loc_JD_time_mid_transit_plus_error_fit, loc_JD_time_mid_transit_minus_error_fit
integer(4), dimension(1) :: loc_Mp_plus_error_fit, loc_Mp_minus_error_fit, loc_Rp_plus_error_fit, loc_Rp_minus_error_fit, loc_Ecc_plus_error_fit
integer(4), dimension(1) :: loc_Ecc_minus_error_fit, loc_Inc_plus_error_fit, loc_Inc_minus_error_fit, loc_omega_arg_periastron_plus_error_fit
integer(4), dimension(1) :: loc_omega_arg_periastron_minus_error_fit, loc_RV_zero_offset_plus_error_fit, loc_RV_zero_offset_minus_error_fit
integer(4), dimension(1) :: loc_min_chi_squared_nettotal, loc_min_reduced_chi_squared_nettotal
integer(4) :: cutoff_position
integer(4) :: num_compare = 0
integer(4) :: num_compare_out
integer(4) :: Number_iterations_left_total
integer(4) :: Number = 0              !Counter.
integer(4) :: Number_points = 1                            !The number of points that are used from the predicted light curve.
integer(4) :: Number_points1 = 1                           !The number of points that are used from the predicted light curve.
integer(4) :: Number_vsini = 1
integer(4) :: Number_vsini_mc
integer(4) :: Number_stellar_rotation_angle = 1
integer(4) :: Number_stellar_rotation_angle_mc
integer(4) :: Number_Orbital_period = 1
integer(4) :: Number_JD_time_mid_transit = 1
integer(4) :: Number_Mp = 1
integer(4) :: Number_Rp = 1
integer(4) :: Number_Ecc = 1
integer(4) :: Number_Inc = 1
integer(4) :: Number_omega_arg_periastron = 1
integer(4) :: Number_RV_zero_offset = 1
integer(4) :: Number_RV_offset_datasets = 1
integer(4) :: num_my_rv
integer(4) :: num_other_rv
integer(4) :: num_degrees_freedom_fit
integer(4) :: num_degrees_freedom_total
integer(4) :: net_num_degrees_freedom_total
integer(4) :: Number_iterations
integer(4) :: Number_iterations_fit
integer(4) :: Number_iterations_total
integer(4) :: Number_iterations_nettotal
integer(4) :: Number_fit
integer(4) :: fit_count = 0
integer(4) :: occultation_End_position
integer(4) :: occultation_mid_position
integer(4) :: Occultation_start_position
integer(4) :: order
integer(4) :: percent_count=1
integer :: percent_int, exitstat, cmdstat
REAL :: percent_processed
integer(4) :: Pixel_interval
integer(4) :: Position_of_min
integer(4) :: Start_position
integer(4) :: Time_loop
integer(4) :: transit_End_position
integer(4) :: transit_End_position_no_inc
integer(4) :: transit_mid_position
integer(4) :: transit_start_position
integer(4) :: Transit_start_no_inc_position
!Loop variables
integer(4) :: i, ii, l, ll, j, jj, s, p, M_pixels, e, f, k, kk, vsini_loop, Num_vsini, stellar_rotation_angle_loop
integer(4) :: Num_stellar_rotation_angle, stellar_rotation_angle_begin_loop, stellar_rotation_angle_end_loop
integer(4) :: stellar_rotation_angle_interval_loop, vsini_end_loop, vsini_begin_loop, vsini_interval_loop
integer(4) :: percent_processed_loop, datafilelength, string_bar_num = 1, datastep, total_interval, Num_RV_offset_datasets
integer(4) :: Num_Inc, Num_omega_arg_periastron, Num_Ecc, Num_Mp, Num_Rp, Num_Orbital_period, Num_JD_time_mid_transit
integer(4) :: Num_RV_zero_offset, Mp_loop, Rp_loop, Orbital_period_loop, Ecc_loop, Inc_loop, JD_time_mid_transit_loop
integer(4) :: omega_arg_periastron_loop, RV_zero_offset_loop, aa, bb, c_c, dd, ee, ff, gg, hh, fit_counter, Monte_carlo_loop
integer(4) :: RV_offset_datasets_loop, zz, count_chi_total, count_chi_fit, chi_counter, k1, mc_sample_size, count_chi_nettotal

character(LEN=100) :: cmd1
character(len=1), Allocatable :: cmdmsg(:)
character(LEN=1) :: add_RV_zero_offset
character(LEN=1) :: Inc_fit
character(LEN=1) :: omega_arg_periastron_fit
character(LEN=1) :: Ecc_fit
character(LEN=1) :: Mp_fit
character(LEN=1) :: Rp_fit
character(LEN=1) :: vsini_fit
character(LEN=1) :: spin_orbit_fit
character(LEN=1) :: stellar_rotation_angle_fit
character(LEN=1) :: Orbital_period_fit
character(LEN=1) :: JD_time_mid_transit_fit
character(LEN=1) :: RV_zero_offset_fit
character(LEN=1) :: RV_offset_datasets_fit
character(LEN=25) :: string_bar
character(LEN=1) :: Jupiter_Earth_units
character(LEN=1) :: Phase_shift
character(LEN=1) :: transit_flag
character(LEN=1) :: time_compare_flag = 'N'
character(LEN=1) :: stellar_rotation_angle_flag = 'N'
character(LEN=1) :: occultation_flag
character(LEN=1) :: other_RV_files
character(LEN=1) :: subtract_JD_mid_time
character(LEN=1) :: flag_exit = 'N'
character(LEN=1) :: fit_flag = 'N', fit_flag1 = 'N', fit_flag2 = 'N', fit_flag3 = 'N', fit_flag4 = 'N', fit_flag5 = 'N', fit_flag6 = 'N'
character(LEN=1) :: fit_flag7 = 'N', fit_flag8 = 'N', fit_flag9 = 'N', fit_flag10 = 'N', array_resorted = 'N', fit_flag0 = 'N'
character(LEN=60) :: input_data_filename = '/Volumes/Seagate Backup Plus Drive/Cyclops/RM/parameters.txt'
character(LEN=60) :: python_program_location = '/Users/z3345211/Dropbox/RV_RM_transit_program/'
character(LEN=150) :: input_RM_filename
character(LEN=150) :: input_my_data_filename
character(LEN=150) :: input_other_data_filename
character(LEN=20) :: Planet_name
character(LEN=150) :: output_temp_data
character(LEN=150) :: directory_location
character(LEN=150) :: mc_parameter_array_input
character(LEN=150) :: output_chi_squared_array_filename
character(LEN=150) :: output_chi_squared_array_total_filename
character(LEN=150) :: output_chi_squared_array_total_filename_mc
character(LEN=150) :: output_chi_squared_array_nettotal_filename
character(LEN=150) :: output_chi_squared_array_1sigcut_filename_mc
character(LEN=150) :: output_chi_squared_array_fit_filename
character(LEN=150) :: output_delta_chi_square_mc_filename
character(LEN=150) :: output_delta_chi_square_mc_filename2
character(LEN=150) :: output_data_fit_array_filename
character(LEN=150) :: output_Ecc_array_filename
character(LEN=150) :: output_Ecc_chi_diff_filename
character(LEN=150) :: output_Inc_array_filename
character(LEN=150) :: output_Inc_chi_diff_filename
character(LEN=150) :: output_JD_time_mid_transit_array_filename
character(LEN=150) :: output_JD_time_mid_transit_chi_diff_filename
character(LEN=150) :: output_Mp_array_filename
character(LEN=150) :: output_Mp_chi_diff_filename
character(LEN=150) :: output_omega_arg_periastron_array_filename
character(LEN=150) :: output_omega_arg_periastron_chi_diff_filename
character(LEN=150) :: output_orbital_period_array_filename
character(LEN=150) :: output_orbital_period_chi_diff_filename
character(LEN=150) :: output_Rp_array_filename
character(LEN=150) :: output_Rp_chi_diff_filename
character(LEN=150) :: output_RV_zero_offset_array_filename
character(LEN=150) :: output_RV_zero_offset_chi_diff_filename
character(LEN=150) :: output_RV_offset_datasets_array_filename
character(LEN=150) :: output_RV_offset_datasets_chi_diff_filename
character(LEN=150) :: output_reduced_chi_array_filename
character(LEN=150) :: output_reduced_chi_array_total_filename
character(LEN=150) :: output_reduced_chi_array_total_filename_mc
character(LEN=150) :: output_reduced_chi_array_nettotal_filename
character(LEN=150) :: output_reduced_chi_array_1sigcut_filename_mc
character(LEN=150) :: output_reduced_chi_array_fit_filename
character(LEN=150) :: output_parameter_array_filename
character(LEN=150) :: output_parameter_array_total_filename
character(LEN=150) :: output_parameter_array_total_mc_filename
character(LEN=150) :: output_parameter_array_nettotal_filename
character(LEN=150) :: output_parameter_array_fit_filename
character(LEN=150) :: output_parameter_idl_filename
character(LEN=150) :: output_parameter_fit_idl_filename
character(LEN=150) :: output_parameter_total_idl_filename
character(LEN=150) :: output_parameter_total_mc_idl_filename
character(LEN=150) :: output_parameter_cutoff_mc_idl_filename
character(LEN=150) :: output_parameter_nettotal_idl_filename
character(LEN=150) :: output_vsini_array_filename
character(LEN=150) :: output_vsini_array_mc_filename
character(LEN=150) :: output_1sigma_vsini_mc_array_filename
character(LEN=150) :: output_stellar_rotation_angle_array_filename
character(LEN=150) :: output_stellar_rotation_angle_array_mc_filename
character(LEN=150) :: output_1sigma_spin_orbit_mc_array_filename
character(LEN=150) :: output_number_parameters_filename
character(LEN=150) :: output_best_parameters_filename
character(LEN=150) :: output_best_parameters_total_filename
character(LEN=150) :: output_best_parameters_mc_filename
character(LEN=150) :: output_best_vsini_mc_array_filename
character(LEN=150) :: output_best_spin_orbit_mc_array_filename
character(LEN=150) :: output_best_parameters_nettotal_filename
character(LEN=150) :: output_best_parameters_fit_filename
character(LEN=150) :: output_best_parameters_netfit_filename
character(LEN=150) :: output_RV_theory_fit_array_filename
character(LEN=150) :: output_RV_theory_all_array_filename
character(LEN=150) :: output_residual_all_array_filename
character(LEN=150) :: output_residual_array_filename
character(LEN=150) :: my_data_file_name
character(LEN=150) :: output_fit_array_indeces_filename
character(LEN=150) :: output_net_degrees_freedom_filename
character(LEN=150) :: input_net_delta_chi_square_filename
character(LEN=150) :: python_command
character(LEN=150) :: output_mc_uncert_vsini_lambda_filename
character(LEN=150) :: output_mc_parm_uncert_vsini_lambda_filename
character(LEN=150) :: output_mc_cut_array_indeces_filename
character(LEN=150) :: my_data_plotting_symbols
character(LEN=150) :: other_data_plotting_symbols
character(LEN=1) :: linear_quadratic
character(LEN=1) :: RV_offset_datasets_free_parameter='Y'
character(LEN=1) :: orbital_period_free_parameter='Y'
character(LEN=1) :: JD_time_mid_transit_free_parameter='Y'
character(LEN=1) :: Mp_free_parameter='Y'
character(LEN=1) :: Rp_free_parameter='Y'
character(LEN=1) :: Ecc_free_parameter='Y'
character(LEN=1) :: Inc_free_parameter='Y'
character(LEN=1) :: omega_arg_periastron_free_parameter='Y'
character(LEN=1) :: RV_zero_offset_free_parameter='Y'

DOUBLE PRECISION, Allocatable :: Data(:,:), Data_my_rv_offset(:,:), Data_my_rv(:,:), Data_other_rv(:,:), Data_adjusted(:,:)
DOUBLE PRECISION, Allocatable :: RV_theory(:,:), time_data_array(:), sorted_data_array(:,:), RV_offset_data_array(:,:)
DOUBLE PRECISION, Allocatable :: RV_synthetic_data_array(:,:), min_chi_squared_mc(:), min_reduced_chi_squared_mc(:), mc_parameter_array(:,:)
DOUBLE PRECISION, Allocatable :: JD_time_mid_transit_array(:), Parameter_array_total(:,:), Parameter_array_fit(:,:)
DOUBLE PRECISION, Allocatable :: mc_best_JD_time_mid_transit_array(:,:), sorted_mc_best_JD_time_mid_transit_array(:,:), cut_mc_best_JD_time_mid_transit_array(:,:)
REAL, Allocatable :: Parameter_chi_fit_array(:,:)
REAL, Allocatable :: Chi_squared_array_fit(:), Chi_squared_array_total(:), Parameter_chi_total_array(:,:)
REAL, Allocatable :: Reduced_Chi_Squared_array_fit(:), Reduced_Chi_Squared_array_total(:), Reduced_Chi_Squared_array_nettotal(:)
REAL, Allocatable :: vsini_array(:), stellar_rotation_angle_array(:), Inc_array(:), Parameter_chi_nettotal_array(:,:)
REAL, Allocatable :: RV_offset_datasets_array(:), chi_square_temp_array(:)
REAL, Allocatable :: omega_arg_periastron_array(:), Ecc_array(:), Mp_array(:), Rp_array(:)
REAL, Allocatable :: Orbital_period_array(:), RV_zero_offset_array(:)
REAL, Allocatable :: mc_best_lambda_array(:,:), mc_best_vsini_array(:,:)
REAL, Allocatable :: mc_best_RV_offset_datasets_array(:,:), mc_best_orbital_period_array(:,:)
REAL, Allocatable :: mc_best_Mp_array(:,:), mc_best_Rp_array(:,:), mc_best_Ecc_array(:,:), mc_best_Inc_array(:,:)
REAL, Allocatable :: mc_best_omega_arg_periastron_array(:,:), mc_best_RV_zero_offset_array(:,:), sorted_min_chi_squared_mc(:)
REAL, Allocatable :: sorted_min_reduced_chi_squared_mc(:), sorted_mc_best_vsini_array(:,:), sorted_mc_best_lambda_array(:,:)
REAL, Allocatable :: sorted_mc_best_RV_offset_datasets_array(:,:), sorted_mc_best_orbital_period_array(:,:)
REAL, Allocatable :: sorted_mc_best_Mp_array(:,:), sorted_mc_best_Rp_array(:,:)
REAL, Allocatable :: sorted_mc_best_Ecc_array(:,:), sorted_mc_best_Inc_array(:,:), sorted_mc_best_omega_arg_periastron_array(:,:)
REAL, Allocatable :: sorted_mc_best_RV_zero_offset_array(:,:), cut_min_chi_squared_mc(:)
REAL, Allocatable :: cut_min_reduced_chi_squared_mc(:), cut_mc_best_vsini_array(:,:), cut_mc_best_lambda_array(:,:)
REAL, Allocatable :: cut_mc_best_RV_offset_datasets_array(:,:), cut_mc_best_orbital_period_array(:,:)
REAL, Allocatable :: cut_mc_best_Mp_array(:,:), cut_mc_best_Rp_array(:,:)
REAL, Allocatable :: cut_mc_best_Ecc_array(:,:), cut_mc_best_Inc_array(:,:)
REAL, Allocatable :: cut_mc_best_omega_arg_periastron_array(:,:), cut_mc_best_RV_zero_offset_array(:,:)

INTEGER(4), Allocatable :: fit_array_indeces(:), Number_array_total(:), chi_position_array(:)
INTEGER(4), DIMENSION(:), ALLOCATABLE :: index_array
CHARACTER(len=1), Allocatable :: data_points_compare(:), check_array_nettotal(:)

!Read in output data from IDL and assaign values to parameters and variables.
OPEN(unit=99, FILE=input_data_filename, status='old', action='read')
!Now read line by line of data file and assaign values to parameters and variable.
READ(99,49,end=50) input_RM_filename, input_my_data_filename, input_other_data_filename, num_my_rv, num_other_rv, &
directory_location, Planet_name, my_data_file_name, output_temp_data, other_RV_files, &
Jupiter_Earth_units, Phase_shift, subtract_JD_mid_time, add_RV_zero_offset, data_plot_model_interval, &
Phase_angle_start, Ms_solar, Inc_prior, Inc_end, Inc_begin, &
Inc_interval, Inc_1sigerr, Inc_fit, Number_fit, Number_orbits, &
Number_iterations_fit, Number_iterations_total, datafilelength, omega_arg_periastron_prior, omega_arg_periastron_end, &
omega_arg_periastron_begin, omega_arg_periastron_interval, omega_arg_periastron_1sigerr, omega_arg_periastron_fit, Ecc_prior, &
Ecc_end, Ecc_begin, Ecc_interval, Ecc_1sigerr, Ecc_fit, &
Mp_prior, Mp_end, Mp_begin, Mp_interval, Mp_1sigerr, &
Mp_fit, Rp_prior, Rp_end, Rp_begin, Rp_interval, &
Rp_1sigerr, Rp_fit, Rs_solar, Orbital_period_prior, Orbital_period_end, &
Orbital_period_begin, Orbital_period_interval, Orbital_period_1sigerr, Orbital_period_fit, JD_time_mid_transit_prior, &
JD_time_mid_transit_end, JD_time_mid_transit_begin, JD_time_mid_transit_interval, JD_time_mid_transit_1sigerr, JD_time_mid_transit_fit, &
vsini_prior, vsini_end, vsini_begin, vsini_interval, vsini_1sigerr, &
vsini_fit, Number_vsini_mc, stellar_rotation_angle_prior, stellar_rotation_angle_end, stellar_rotation_angle_begin, &
stellar_rotation_angle_interval, stellar_rotation_angle_1sigerr, stellar_rotation_angle_fit, Number_stellar_rotation_angle_mc, vmacro, &
vmicro, RV_zero_offset_prior, RV_zero_offset_end, RV_zero_offset_begin, RV_zero_offset_interval, &
RV_zero_offset_1sigerr, RV_zero_offset_fit, RV_offset_datasets_prior, RV_offset_datasets_end, RV_offset_datasets_begin, &
RV_offset_datasets_interval, RV_offset_datasets_1sigerr, RV_offset_datasets_fit, Time_bef_aft, Time_compare_vel, &
M_pixels, Albedo, linear_quadratic, u, q_1, &
q_2, Bessel_function_exit, chi_squared_change, chi_squared_change_fit, num_degrees_freedom_fit, &
num_degrees_freedom_total, mc_sample_size, my_data_plotting_symbols, other_data_plotting_symbols
49 FORMAT(A150, /, A150, /, A150, /, I20, /, I20, /, &
          A150, /, A150, /, A150, /, A150, /, A1, /, &
          A1, /, A1, /, A1, /, A1, /, F50.10, /, &
          F50.5, /, F50.5, /, F50.5, /, F50.5, /, F50.5, /, &
          F50.5, /, F50.5, /, A1, /, I20, /, F10.5, /, &
          I20, /, I20, /, I20, /, F50.5, /, F50.5, /, &
          F50.5, /, F50.5, /, F50.5, /, A1, /, F50.5, /, &
          F50.5, /, F50.5, /, F50.5, /, F50.5, /, A1, /, &
          F50.5, /, F50.5, /, F50.5, /, F50.5, /, F50.5, /, &
          A1, /, F50.5, /, F50.5, /, F50.5, /, F50.5, /, &
          F50.5, /, A1, /, F50.5, /, F50.10, /, F50.10, /, &
          F50.10, /, F50.10, /, F50.10, /, A1, /, F50.10, /, &
          F50.10, /, F50.10, /, F50.10, /, F50.10, /, A1, /, &
          F50.5, /, F50.5, /, F50.5, /, F50.5, /, F50.5, /, &
          A1, /, I20, /, F50.5, /, F50.5, /, F50.5, /, &
          F50.5, /, F50.5, /, A1, /, I20, /, F50.5, /, &
          F50.5, /, F50.5, /, F50.5, /, F50.5, /, F50.5, /, &
          F50.5, /, A1, /, F50.5, /, F50.5, /, F50.5, /, &
          F50.5, /, F50.5, /, A1, /, F50.10, /, F50.10, /, &
          I20, /, F50.5, /, A1, /, F50.5, /, F50.5, /, &
          F50.5, /, E50.10, /, F50.5, /, F50.5, /, I20, /, &
          I20, /, I20, /, A150, /, A150)
50 CLOSE(99)

input_RM_filename = ADJUSTL(input_RM_filename)
input_RM_filename = TRIM(input_RM_filename)
PRINT *, 'input_RM_filename = ', input_RM_filename

input_my_data_filename = ADJUSTL(input_my_data_filename)
input_my_data_filename = TRIM(input_my_data_filename)
PRINT *, 'input_my_data_filename = ', input_my_data_filename

input_other_data_filename = ADJUSTL(input_other_data_filename)
input_other_data_filename = TRIM(input_other_data_filename)
PRINT *, 'input_other_data_filename = ', input_other_data_filename

directory_location = ADJUSTL(directory_location)
directory_location = TRIM(directory_location)
PRINT *, 'directory_location = ', directory_location
Planet_name = ADJUSTL(Planet_name)
Planet_name = TRIM(Planet_name)
PRINT *, 'Planet_name = ', Planet_name
output_temp_data = ADJUSTL(output_temp_data)
output_temp_data = TRIM(output_temp_data)
PRINT *, 'output_temp_data = ', output_temp_data

mc_parameter_array_input = TRIM(TRIM(output_temp_data) // 'mc_parameter_array.txt')
print *, 'mc_parameter_array_input: ', mc_parameter_array_input

output_chi_squared_array_total_filename = TRIM(TRIM(output_temp_data) // 'chi_squared_array_total.txt')
output_chi_squared_array_total_filename_mc = TRIM(TRIM(output_temp_data) // 'chi_squared_array_total_mc.txt')
output_chi_squared_array_nettotal_filename = TRIM(TRIM(output_temp_data) // 'chi_squared_array_nettotal.txt')
output_chi_squared_array_1sigcut_filename_mc = TRIM(TRIM(output_temp_data) // 'chi_squared_array_1sigcut_mc.txt')
output_chi_squared_array_fit_filename = TRIM(TRIM(output_temp_data) // 'chi_squared_array_fit.txt')
output_delta_chi_square_mc_filename = TRIM(TRIM(output_temp_data) // 'delta_chi_square_mc.txt')
output_delta_chi_square_mc_filename2 = TRIM(TRIM(output_temp_data) // 'delta_chi_square_mc2.txt')
output_reduced_chi_array_total_filename = TRIM(TRIM(output_temp_data) // 'reduced_chi_squared_array_total.txt')
output_reduced_chi_array_total_filename_mc = TRIM(TRIM(output_temp_data) // 'reduced_chi_squared_array_total_mc.txt')
output_reduced_chi_array_nettotal_filename = TRIM(TRIM(output_temp_data) // 'reduced_chi_squared_array_nettotal.txt')
output_reduced_chi_array_1sigcut_filename_mc = TRIM(TRIM(output_temp_data) // 'reduced_chi_squared_array_1sigcut_mc.txt')
output_reduced_chi_array_fit_filename = TRIM(TRIM(output_temp_data) // 'reduced_chi_squared_array_fit.txt')
output_vsini_array_filename = TRIM(TRIM(output_temp_data) // 'vsini_array.txt')
output_vsini_array_mc_filename = TRIM(TRIM(output_temp_data) // 'vsini_array_mc.txt')
output_best_vsini_mc_array_filename = TRIM(TRIM(output_temp_data) // 'best_vsini_mc_array.txt')
output_1sigma_vsini_mc_array_filename = TRIM(TRIM(output_temp_data) // '1sigma_vsini_mc_array.txt')
output_best_spin_orbit_mc_array_filename = TRIM(TRIM(output_temp_data) // 'best_stellar_rotation_angle_mc_array.txt')
output_1sigma_spin_orbit_mc_array_filename = TRIM(TRIM(output_temp_data) // '1sigma_stellar_rotation_angle_mc_array.txt')
output_stellar_rotation_angle_array_filename = TRIM(TRIM(output_temp_data) // 'stellar_rotation_angle_array.txt')
output_stellar_rotation_angle_array_mc_filename = TRIM(TRIM(output_temp_data) // 'stellar_rotation_angle_array_mc.txt')
output_orbital_period_array_filename = TRIM(TRIM(output_temp_data) // 'orbital_period_array.txt')
output_JD_time_mid_transit_array_filename = TRIM(TRIM(output_temp_data) // 'JD_time_mid_transit_array.txt')
output_Mp_array_filename = TRIM(TRIM(output_temp_data) // 'Mp_array.txt')
output_Rp_array_filename = TRIM(TRIM(output_temp_data) // 'Rp_array.txt')
output_Ecc_array_filename = TRIM(TRIM(output_temp_data) // 'Ecc_array.txt')
output_Inc_array_filename = TRIM(TRIM(output_temp_data) // 'Inc_array.txt')
output_omega_arg_periastron_array_filename = TRIM(TRIM(output_temp_data) // 'omega_arg_periastron_array.txt')
output_RV_zero_offset_array_filename = TRIM(TRIM(output_temp_data) // 'RV_zero_offset_array.txt')
output_RV_offset_datasets_array_filename = TRIM(TRIM(output_temp_data) // 'RV_offset_datasets_array.txt')
output_fit_array_indeces_filename = TRIM(TRIM(output_temp_data) // 'fit_array_indeces.txt')
output_parameter_array_total_filename = TRIM(TRIM(output_temp_data) // 'parameter_total_array.txt')
output_parameter_array_total_mc_filename = TRIM(TRIM(output_temp_data) // 'parameter_total_mc_array.txt')
output_parameter_array_nettotal_filename = TRIM(TRIM(output_temp_data) // 'parameter_nettotal_array.txt')
output_parameter_array_fit_filename = TRIM(TRIM(output_temp_data) // 'parameter_fit_array.txt')
output_residual_array_filename = TRIM(TRIM(output_temp_data) // 'residual_array.txt')
output_data_fit_array_filename = TRIM(TRIM(output_temp_data) // 'data_fit_array.txt')
output_RV_theory_fit_array_filename = TRIM(TRIM(output_temp_data) // 'RV_Theory_fit_array.txt')
output_residual_all_array_filename = TRIM(TRIM(output_temp_data) // 'residual_all_array.txt')
output_RV_theory_all_array_filename = TRIM(TRIM(output_temp_data) // 'RV_Theory_all_array.txt')
output_parameter_fit_idl_filename = TRIM(TRIM(output_temp_data) // 'parameter_fit_idl.txt')
output_parameter_total_idl_filename = TRIM(TRIM(output_temp_data) // 'parameter_total_idl.txt')
output_parameter_total_mc_idl_filename = TRIM(TRIM(output_temp_data) // 'parameter_total_mc_idl.txt')
output_parameter_cutoff_mc_idl_filename = TRIM(TRIM(output_temp_data) // 'parameter_cutoff_mc_idl.txt')
output_parameter_nettotal_idl_filename = TRIM(TRIM(output_temp_data) // 'parameter_nettotal_idl.txt')
output_number_parameters_filename = TRIM(TRIM(output_temp_data) // 'number_parameters.txt')
output_best_parameters_filename = TRIM(TRIM(output_temp_data) // 'best_parameters.txt')
output_best_parameters_total_filename = TRIM(TRIM(output_temp_data) // 'best_parameters_total.txt')
output_best_parameters_mc_filename = TRIM(TRIM(output_temp_data) // 'best_parameters_mc.txt')
output_best_parameters_nettotal_filename = TRIM(TRIM(output_temp_data) // 'best_parameters_nettotal.txt')
output_best_parameters_fit_filename = TRIM(TRIM(output_temp_data) // 'best_parameters_fit.txt')
output_best_parameters_netfit_filename = TRIM(TRIM(output_temp_data) // 'best_parameters_fit_nettotal.txt')
output_RV_offset_datasets_chi_diff_filename = TRIM(TRIM(output_temp_data) // 'RV_offset_datasets_chi_diff_array.txt')
output_orbital_period_chi_diff_filename = TRIM(TRIM(output_temp_data) // 'orbital_period_chi_diff_array.txt')
output_JD_time_mid_transit_chi_diff_filename = TRIM(TRIM(output_temp_data) // 'JD_time_mid_transit_chi_diff_array.txt')
output_Mp_chi_diff_filename = TRIM(TRIM(output_temp_data) // 'Mp_chi_diff_array.txt')
output_Rp_chi_diff_filename = TRIM(TRIM(output_temp_data) // 'Rp_chi_diff_array.txt')
output_Ecc_chi_diff_filename = TRIM(TRIM(output_temp_data) // 'Ecc_chi_diff_array.txt')
output_Inc_chi_diff_filename = TRIM(TRIM(output_temp_data) // 'Inc_chi_diff_array.txt')
output_omega_arg_periastron_chi_diff_filename = TRIM(TRIM(output_temp_data) // 'omega_arg_periastron_chi_diff_array.txt')
output_RV_zero_offset_chi_diff_filename = TRIM(TRIM(output_temp_data) // 'RV_zero_offset_chi_diff_array.txt')
output_net_degrees_freedom_filename = TRIM(TRIM(output_temp_data) // 'net_degrees_freedom.txt')
input_net_delta_chi_square_filename = TRIM(TRIM(output_temp_data) // 'net_delta_chi_square.txt')
output_mc_uncert_vsini_lambda_filename = TRIM(TRIM(output_temp_data) // 'mc_uncertainty_vsini_lambda.txt')
output_mc_parm_uncert_vsini_lambda_filename = TRIM(TRIM(output_temp_data) // 'mc_parm_uncertainty_vsini_lambda.txt')
output_mc_cut_array_indeces_filename = TRIM(TRIM(output_temp_data) // 'mc_cut_array_indeces.txt')
python_command = TRIM('python ' // TRIM(python_program_location) // 'critical_chi_square_value.py')

Num_vsini = NINT(((vsini_end - vsini_begin)/vsini_interval) + 1)               !Number of interations for the vsini loop.
!Number of interations for the stellar_rotation_angle loop.
Num_stellar_rotation_angle = NINT(((stellar_rotation_angle_end - stellar_rotation_angle_begin) &
                                 /stellar_rotation_angle_interval) + 1)

!Increase mc_sample_size by 1. Last iteration of MC is with the psuedo free parameters fixed to their prior values.
mc_sample_size = mc_sample_size + 1
!Number_iterations_total refers to the total number of iterations using the grid search. We will
!redifine it here to mean the number of Monte Carlo iterations times the number of vsini and
!times the number of lambda.
Number_iterations_total = Number_vsini_mc*Number_stellar_rotation_angle_mc*mc_sample_size

!Number_iterations_fit refers to the total number of iterations from parameters that are being fitted using
!the grid search. We will redifine it here to mean the number of vsini times the number of lambda.
Number_iterations_fit = Number_vsini_mc*Number_stellar_rotation_angle_mc

!Allocate the monte carlo parameter array as well as the other parameter arrays.
Allocate(mc_parameter_array(mc_sample_size,9))
Allocate(Data_my_rv(num_my_rv,3))
Allocate(Data_other_rv(num_other_rv,3))
Allocate(min_chi_squared_mc(mc_sample_size))
Allocate(min_reduced_chi_squared_mc(mc_sample_size))
Allocate(mc_best_vsini_array(mc_sample_size,6))
Allocate(mc_best_lambda_array(mc_sample_size,6))
Allocate(mc_best_RV_offset_datasets_array(mc_sample_size,3))
Allocate(mc_best_orbital_period_array(mc_sample_size,3))
Allocate(mc_best_JD_time_mid_transit_array(mc_sample_size,3))
Allocate(mc_best_Mp_array(mc_sample_size,3))
Allocate(mc_best_Rp_array(mc_sample_size,3))
Allocate(mc_best_Ecc_array(mc_sample_size,3))
Allocate(mc_best_Inc_array(mc_sample_size,3))
Allocate(mc_best_omega_arg_periastron_array(mc_sample_size,3))
Allocate(mc_best_RV_zero_offset_array(mc_sample_size,3))

!Will determine the size of these other arrays soon. Placed them here for the time being.
Allocate(Number_array_total(Number_iterations_total))
!Allocate(check_array_nettotal(Number_iterations_total))

!The size of the following array are based on the mc_sample_size.
Allocate(Inc_array(mc_sample_size))
Allocate(omega_arg_periastron_array(mc_sample_size))
Allocate(Ecc_array(mc_sample_size))
Allocate(Rp_array(mc_sample_size))
Allocate(Mp_array(mc_sample_size))
Allocate(Orbital_period_array(mc_sample_size))
Allocate(JD_time_mid_transit_array(mc_sample_size))
Allocate(RV_zero_offset_array(mc_sample_size))
Allocate(RV_offset_datasets_array(mc_sample_size))
Allocate(Chi_squared_array_total(Number_iterations_total))
Allocate(Reduced_Chi_Squared_array_total(Number_iterations_total))
Allocate(Parameter_array_total(Number_iterations_total,11))
!Read in the monte carlo parameter set.
OPEN(unit=99, FILE=mc_parameter_array_input, status='old', action='read')
DO i=1, mc_sample_size - 1
   READ(99,48) mc_parameter_array(i,1), mc_parameter_array(i,2), mc_parameter_array(i,3), mc_parameter_array(i,4), mc_parameter_array(i,5), &
   mc_parameter_array(i,6), mc_parameter_array(i,7), mc_parameter_array(i,8), mc_parameter_array(i,9)
   48 FORMAT(F20.8, 2X, F20.8, 2X, F20.8, 2X, F20.8, 2X, F20.8, 2X, F20.8, 2X, F20.8, 2X, F20.8, 2X, F20.8)
   RV_offset_datasets_array(i) = mc_parameter_array(i,1)
   Orbital_period_array(i) = mc_parameter_array(i,2)
   JD_time_mid_transit_array(i) = mc_parameter_array(i,3)
   Mp_array(i) = mc_parameter_array(i,4)
   Rp_array(i) = mc_parameter_array(i,5)
   Ecc_array(i) = mc_parameter_array(i,6)
   Inc_array(i) = mc_parameter_array(i,7)
   omega_arg_periastron_array(i) = mc_parameter_array(i,8)
   RV_zero_offset_array(i) = mc_parameter_array(i,9)
END DO
CLOSE(99)

RV_offset_datasets_array(mc_sample_size) = RV_offset_datasets_prior
Orbital_period_array(mc_sample_size) = Orbital_period_prior
JD_time_mid_transit_array(mc_sample_size) = JD_time_mid_transit_prior
Mp_array(mc_sample_size) = Mp_prior
Rp_array(mc_sample_size) = Rp_prior
Ecc_array(mc_sample_size) = Ecc_prior
Inc_array(mc_sample_size) = Inc_prior
omega_arg_periastron_array(mc_sample_size) = omega_arg_periastron_prior
RV_zero_offset_array(mc_sample_size) = RV_zero_offset_prior
!We can go ahead and populate the grid of values for the parameters vsini and lambda.
!Using do concurrent loop to allow parallization.
!DO CONCURRENT (vsini_loop = 1:Num_vsini)
!   vsini_array(vsini_loop) = vsini_end - (vsini_interval * (vsini_loop - 1))
!END DO

!DO CONCURRENT (stellar_rotation_angle_loop = 1:Num_stellar_rotation_angle)
!   stellar_rotation_angle_array(stellar_rotation_angle_loop) = stellar_rotation_angle_end - (stellar_rotation_angle_interval * (stellar_rotation_angle_loop - 1))
!END DO

!Now let's read in the data.
!Read in output data from IDL and assaign values to the data array.
OPEN(unit=99, FILE=input_my_data_filename, status='old', action='read')
!Now read line by line of data file and dump values into the data array.

DO l = 1, num_my_rv
   READ (99, 51) a, b, c
   51 FORMAT(E50.10, 5X, F50.5, 5X, F50.5)
   Data_my_rv(l,1) = a
   Data_my_rv(l,2) = b
   Data_my_rv(l,3) = c
END DO 

CLOSE(99)

IF (other_RV_files == 'Y') THEN
   !Read in output data from IDL and assaign values to the data array.
   OPEN(unit=99, FILE=input_other_data_filename, status='old', action='read')
   !Now read line by line of data file and dump values into the data array.

   DO l = 1, num_other_rv
      READ (99, 52) a, b, c
      52 FORMAT(E50.10, 5X, F50.5, 5X, F50.5)
      Data_other_rv(l,1) = a
      Data_other_rv(l,2) = b
      Data_other_rv(l,3) = c
   END DO 

   CLOSE(99)
END IF

e = 1
k = 1
ii = 1

Num_Inc = mc_sample_size
Num_omega_arg_periastron = mc_sample_size
Num_Ecc = mc_sample_size
Num_Mp = mc_sample_size
Num_Rp = mc_sample_size
Num_Orbital_period = mc_sample_size
Num_JD_time_mid_transit = mc_sample_size
Num_RV_zero_offset = mc_sample_size
Num_RV_offset_datasets = mc_sample_size

Number_iterations_total = mc_sample_size*Number_vsini_mc*Number_stellar_rotation_angle_mc

Number_iterations_left_total = Number_iterations_total

PRINT *, "Number_vsini_mc: ", Number_vsini_mc
PRINT *, "Number_stellar_rotation_angle_mc ", Number_stellar_rotation_angle_mc
PRINT *, "Num_Inc: ", Num_Inc
PRINT *, "Num_omega_arg_periastron: ", Num_omega_arg_periastron
PRINT *, "Num_Ecc: ", Num_Ecc
PRINT *, "Num_Mp: ", Num_Mp
PRINT *, "Num_Rp: ", Num_Rp
PRINT *, "Num_Orbital_period: ", Num_Orbital_period
PRINT *, "Num_JD_time_mid_transit: ", Num_JD_time_mid_transit
PRINT *, "Num_RV_zero_offset: ", Num_RV_zero_offset
PRINT *, "Num_RV_offset_datasets: ", Num_RV_offset_datasets
PRINT *, "Number_iterations_total: ", Number_iterations_total

Rs = Rss * Rs_solar
Ms = Ms_solar * Mss

total_interval = datafilelength

Rs2 = Rs**2.0D0                              !Square the radius of star to speed up calculations.

IF (linear_quadratic == 'q') THEN
   Io = 6.0D0 / (pi*Rs2*(6.0D0 - (2.0D0*q_1) - q_2))     !The initial light intensity equation with limb darkening (quadratic)
ELSE
   Io = 1.0D0 / (pi*Rs2*(1.0D0-(u/3.0D0)))      !The initial light intensity equation with limb darkening (linear)
                                                !(normalize Io such that total star luminosity is 1). 
END IF

vturb = SQRT(vmicro**2.0D0 + vmacro**2.0D0)            !Velocity width of spectral line due to mechanisms other than rotation (i.e. micro and macro turbulence).

!Based on given priors, determine the approximate time in the simulation to start and finish comparing RV's with model.
IF (Jupiter_Earth_units == 'Y') THEN
   Mp = Mp_prior*Mj
   Rp = Rp_prior*Rj
END IF

IF (Jupiter_Earth_units == 'N') THEN
   Mp = Mp_prior*Me
   Rp = Rp_prior*RE
END IF
Rorb_prior = (((Orbital_period_prior*day_sec)**2.0d0*G*(Ms + Mp))/(4.0d0*pi**2.0d0))**(1.0d0/3.0d0)     !Semi-major axis.
Rorb_star_prior = (((Orbital_period_prior*day_sec)**2.0D0*G*((Mp**(3.0D0))/(Ms + Mp)**(2.0D0)))/(4.0d0*pi**2.0D0))**(1.0d0/3.0d0)     !Semi-major of stellar axis.
radius_ratio_prior = Rp/Rs
impact_prior = ((Rorb_prior*cos(Inc_prior*(pi/180.0D0)))/Rs)*((1.0D0 - Ecc_prior**2.0D0)/(1.0D0 + (Ecc_prior*sin(omega_arg_periastron_prior*(pi/180.0D0)))))
Transit_length_prior = ((Orbital_period_prior*day_sec)/pi) * asin((Rs/Rorb_prior)*(sqrt((1.0D0 + radius_ratio_prior)**2.0D0 - impact_prior**2.0D0)/sin(Inc_prior*(pi/180.0D0)))) * &
                       (sqrt(1.0D0 - Ecc_prior**2.0D0)/(1.0D0 + (Ecc_prior*sin(omega_arg_periastron_prior*(pi/180.0D0)))))
PRINT *, "Orbital_period_prior: ", Orbital_period_prior
PRINT *, "Orbital_period_prior * day_sec: ", Orbital_period_prior*day_sec
PRINT *, "Rorb_prior: ", Rorb_prior
PRINT *, "Rorb_star_prior: ", Rorb_star_prior
PRINT *, "Inc_prior: ", Inc_prior
PRINT *, "Rp_prior: ", Rp_prior
PRINT *, "Mp_prior: ", Mp_prior
PRINT *, "Ecc_prior: ", Ecc_prior
PRINT *, "omega_arg_periastron_prior: ", omega_arg_periastron_prior
PRINT *, "Rs: ", Rs
PRINT *, "Impact parameter prior: ", impact_prior
PRINT *, "Radius ratio prior: ", radius_ratio_prior
PRINT *, "Length of transit with given priors: ", Transit_length_prior
length_time_compare = (Transit_length_prior/2.0D0) + Time_compare_vel
PRINT *, "Length of time to compare: ", length_time_compare

call cpu_time(cpu_time_1)
call cpu_time(cpu_time_it1)
call cpu_time(cpu_time_it2)

!STOP

!************************************************************************************************************************************************
!Start the Monte Carlo simulation.
DO Monte_carlo_loop = ii, mc_sample_size
e = 1
k = 1
!reset the fit_count
fit_count = 0

RV_offset_datasets = RV_offset_datasets_array(Monte_carlo_loop)
Orbital_period_1 = Orbital_period_array(Monte_carlo_loop)
Orbital_period_2 = (Orbital_period_1 - Orbital_period_prior)*10000.0D0
Orbital_period = Orbital_period_1*day_sec               !Convert orbital period in days to seconds.
JD_time_mid_transit = JD_time_mid_transit_array(Monte_carlo_loop)
JD_time_mid_transit_1 = (JD_time_mid_transit - JD_time_mid_transit_prior)*100.0D0
Mpp = Mp_array(Monte_carlo_loop)
Rpp = Rp_array(Monte_carlo_loop)
Ecc = Ecc_array(Monte_carlo_loop)
Inc = Inc_array(Monte_carlo_loop)
omega_arg_periastron = omega_arg_periastron_array(Monte_carlo_loop)
RV_zero_offset = RV_zero_offset_array(Monte_carlo_loop)

!Now apply the RV offset to my data.
Allocate(Data_my_rv_offset(num_my_rv,3), data(datafilelength,3))

DO l = 1, num_my_rv
   Data_my_rv_offset(l,1) = Data_my_rv(l,1)
   Data_my_rv_offset(l,2) = Data_my_rv(l,2) + RV_offset_datasets
   Data_my_rv_offset(l,3) = Data_my_rv(l,3)
END DO

IF (other_RV_files == 'Y') THEN
   DO l = 1, datafilelength
      IF (l <= num_other_rv) THEN
         Data(l,1) = Data_other_rv(l,1)
         Data(l,2) = Data_other_rv(l,2)
         Data(l,3) = Data_other_rv(l,3)
      ELSE
         zz = l - num_other_rv
         Data(l,1) = Data_my_rv_offset(zz,1)
         Data(l,2) = Data_my_rv_offset(zz,2)
         Data(l,3) = Data_my_rv_offset(zz,3)
      END IF
   END DO
ELSE
   DO l = 1, datafilelength
      Data(l,1) = Data_my_rv_offset(l,1)
      Data(l,2) = Data_my_rv_offset(l,2)
      Data(l,3) = Data_my_rv_offset(l,3)
   END DO
END IF

Allocate(Data_adjusted(datafilelength,3), sorted_data_array(datafilelength,3))
!Adjust data time relative to mid transit time based on the new orbital period and JD time mid transit.
DO l = 1, datafilelength
   IF (((Data(l,1) - (JD_time_mid_transit - JD_time_mid_transit_prior))/Orbital_period_1) - (FLOOR((Data(l,1) - &
      (JD_time_mid_transit - JD_time_mid_transit_prior))/Orbital_period_1)) < 0.5) THEN

      phase_time = ((Data(l,1) - (JD_time_mid_transit - JD_time_mid_transit_prior))/Orbital_period_1) - &
                   (FLOOR((Data(l,1) - (JD_time_mid_transit - JD_time_mid_transit_prior))/Orbital_period_1)) + 1.0D0

   ELSE
      phase_time = ((Data(l,1) - (JD_time_mid_transit - JD_time_mid_transit_prior))/Orbital_period_1) - &
                   (FLOOR((Data(l,1) - (JD_time_mid_transit - JD_time_mid_transit_prior))/Orbital_period_1))
   END IF

   Data_adjusted(l,1) = (Orbital_period_1*phase_time) - Orbital_period_1
   Data_adjusted(l,2) = Data(l,2)
   Data_adjusted(l,3) = Data(l,3)
END DO

!Verify that data array is ascending order as a function of time
CALL Selection_sort(Data_adjusted(:,1),index_array)
  
DO l = 1, datafilelength
   sorted_data_array(l,1) = Data_adjusted(index_array(l), 1)
   sorted_data_array(l,2) = Data_adjusted(index_array(l), 2)
   sorted_data_array(l,3) = Data_adjusted(index_array(l), 3)
!   PRINT *, "Data(l,1): ", Data(l,1)
!   PRINT *, "Data_adjusted(l,1): ", Data_adjusted(l,1)
!   PRINT *, "sorted_data_array(l,1): ", sorted_data_array(l,1)
   IF (index_array(l) /= l) THEN
      !PRINT *, "index_array(l): ", index_array(l)
      !PRINT *, "l: ", l
      array_resorted = 'Y'
   ELSE
      array_resorted = 'N'
   END IF
END DO

IF (Jupiter_Earth_units == 'Y') THEN
   Mp = Mpp*Mj
END IF

IF (Jupiter_Earth_units == 'N') THEN
   Mp = Mpp*Me
END IF
Rorb = ((Orbital_period**2.0d0*G*(Ms + Mp))/(4.0d0*pi**2.0d0))**(1.0d0/3.0d0)     !Semi-major axis.
Rorb_star = ((Orbital_period**2.0D0*G*((Mp**(3.0D0))/(Ms + Mp)**(2.0D0)))/(4.0d0*pi**2.0D0))**(1.0d0/3.0d0)     !Semi-major of stellar axis.

IF (Jupiter_Earth_units == 'Y') THEN
   Rp = Rpp*Rj
END IF

IF (Jupiter_Earth_units == 'N') THEN
   Rp = Rpp*RE
END IF

Aplan = pi*Rp**2.0D0                         !Surface area of the planet.

Rp2 = Rp**2.0D0                              !Square the radius of planet to speed up calculations.

IF (Ecc == 0) THEN 
   !Maximum amplitude caused by the exoplanet in a circular orbit.
   !RVamplitude = (Mp/Ms)*sqrt((G*(Ms + Mp))/Rorb)*sin(Inc*(pi/180.0D0))
   RVamplitude = SQRT((G*(Mp**(3.0D0))*(sin(Inc*(pi/180.0D0)))**(3.0D0))/(Rorb_star*sin(Inc*(pi/180.0D0))*(Ms + Mp)**(2.0D0)))
ELSE
   !Maximum amplitude caused by the exoplanet in an eccentric orbit.
   !RVamplitude = (1.0D0/sqrt(1.0D0 - Ecc**2.0D0))*(Mp/Ms)*sqrt((G*(Ms + Mp))/Rorb)*sin(Inc*(pi/180.0D0))
   RVamplitude = SQRT((G*(Mp**(3.0D0))*(sin(Inc*(pi/180.0D0)))**(3.0D0))/((1.0D0 - Ecc**2.0D0)*Rorb_star*sin(Inc*(pi/180.0D0))*(Ms + Mp)**(2.0D0)))
END IF

True_anomaly_start = pi - (omega_arg_periastron*(pi/180.0D0))
!PRINT *, "True_anomaly_start", True_anomaly_start*(180.0D0/pi)

IF (True_anomaly_start >= pi) THEN
   ecc_anomaly_start = 2.0D0*atan(tan(True_anomaly_start/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc))) + 2.0D0*pi 
ELSE IF (True_anomaly_start <= -pi) THEN
   ecc_anomaly_start = 2.0D0*atan(tan(True_anomaly_start/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc))) - 2.0D0*pi
ELSE
   ecc_anomaly_start = 2.0D0*atan(tan(True_anomaly_start/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc)))
END IF

!Calculate the amount of time (in seconds) after the passage of periastron has occurred when the planet is at mid transit.
!This is determined from the transit mid time and the argument of periastron.
IF (omega_arg_periastron > 180.0D0) THEN
   Mean_anomaly_transit = (2.0D0*pi) - (omega_arg_periastron*(pi/180.0D0)) + pi
ELSE
   Mean_anomaly_transit = pi - (omega_arg_periastron*(pi/180.0D0))
END IF

JD_time_peri = (JD_time_mid_transit*day_sec) - ((Mean_anomaly_transit*Orbital_period)/(2.0D0*pi))

time_peri_passage = (JD_time_mid_transit*day_sec) - JD_time_peri

True_anomaly_transit = ((3.0D0*pi)/2.0D0) - (omega_arg_periastron*(pi/180.0D0))
!PRINT *, "True_anomaly_transit", True_anomaly_transit*(180.0D0/pi)

IF (True_anomaly_transit >= pi) THEN
   ecc_anomaly_transit = 2.0D0*atan(tan(True_anomaly_transit/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc))) + 2.0D0*pi 
ELSE IF (True_anomaly_transit <= -pi) THEN
   ecc_anomaly_transit = 2.0D0*atan(tan(True_anomaly_transit/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc))) - 2.0D0*pi
ELSE
   ecc_anomaly_transit = 2.0D0*atan(tan(True_anomaly_transit/2.0D0)*sqrt((1.0D0-Ecc)/(1.0D0+Ecc)))
END IF

IF (Ecc == 0) THEN
   Time_transit = ((ecc_anomaly_transit*Orbital_period)/(2.0D0*pi)) + time_peri_passage
ELSE
   Time_transit = (((ecc_anomaly_transit - (Ecc*sin(ecc_anomaly_transit)))*Orbital_period)/(2.0D0*pi)) + time_peri_passage
END IF

!Create a new time array for the data which has been created for the JD time of mid transit and ecc anomalies.
Allocate(time_data_array(datafilelength))
DO l = 1, datafilelength
   time_data_array(l) = (sorted_data_array(l,1)*day_sec) + Time_transit
!   PRINT *, "time_data_array(l): ", time_data_array(l)
END DO

!Now apply the change in the RV zero offset from the prior offset to the data.
Allocate(RV_offset_data_array(datafilelength,3))
DO l = 1, datafilelength
   RV_offset_data_array(l,1) = sorted_data_array(l,1)
   RV_offset_data_array(l,2) = sorted_data_array(l,2) + (RV_zero_offset - RV_zero_offset_prior)
   RV_offset_data_array(l,3) = sorted_data_array(l,3)
!   PRINT *, "sorted_data_array(l,2): ", sorted_data_array(l,2)
END DO

!Let's create a synthetic data set, using the actual data, that is normally distributed about the 1 sigma uncertainty of each data point.
!I have included the required subroutines and functions for doing this at the end of this program.
Allocate(RV_synthetic_data_array(datafilelength,3))
DO l = 1, datafilelength
   IF (Monte_carlo_loop == mc_sample_size) THEN
      RV_synthetic_data_array(l,1) = RV_offset_data_array(l,1)
      RV_synthetic_data_array(l,2) = RV_offset_data_array(l,2)
      RV_synthetic_data_array(l,3) = RV_offset_data_array(l,3)
   ELSE   
      RV_synthetic_data_array(l,1) = RV_offset_data_array(l,1)
      RV_synthetic_data_array(l,2) = normal(RV_offset_data_array(l,2), RV_offset_data_array(l,3))
      RV_synthetic_data_array(l,3) = RV_offset_data_array(l,3)
   END IF
!   PRINT *, '************************************'
!   PRINT *, "RV_offset_data_array(l,1): ", RV_offset_data_array(l,1)
!   PRINT *, "RV_offset_data_array(l,2): ", RV_offset_data_array(l,2)
!   PRINT *, "RV_offset_data_array(l,3): ", RV_offset_data_array(l,3)
!   PRINT *, ''
!   PRINT *, "RV_synthetic_data_array(l,1): ", RV_synthetic_data_array(l,1)
!   PRINT *, "RV_synthetic_data_array(l,2): ", RV_synthetic_data_array(l,2)
!   PRINT *, "RV_synthetic_data_array(l,3): ", RV_synthetic_data_array(l,3)
!   PRINT *, '************************************'
END DO

!STOP

!Allocate the parameter arrays for just iterating over vsini and lambda.
Allocate(Parameter_array_fit(Number_iterations_fit,11))
Allocate(Chi_squared_array_fit(Number_iterations_fit))
Allocate(Reduced_Chi_Squared_array_fit(Number_iterations_fit))
Allocate(fit_array_indeces(Number_iterations_fit))

!Allocate and then generate a random evenly distributed values for vsini and lambda for each MC iteration.
Allocate(vsini_array(Number_vsini_mc))
Allocate(stellar_rotation_angle_array(Number_stellar_rotation_angle_mc))

DO l = 1, Number_vsini_mc
   vsini_array(l) = spread_rand(vsini_begin, vsini_end)
!   PRINT *, 'grid vsini_array: ', vsini_begin + (l*vsini_interval)
!   PRINT *, 'random vsini_array: ', vsini_array(l)
!   PRINT *, '*************************************'
END DO

!STOP

DO l = 1, Number_stellar_rotation_angle_mc
   stellar_rotation_angle_array(l) = spread_rand(stellar_rotation_angle_begin, stellar_rotation_angle_end)
!   PRINT *, 'grid stellar_rotation_angle_array: ', stellar_rotation_angle_begin + (l*stellar_rotation_angle_interval)
!   PRINT *, 'random stellar_rotation_angle_array: ', stellar_rotation_angle_array(l)
!   PRINT *, '*************************************'
END DO

!STOP

!PRINT *, 'Made it to line 1271'

!MCMC routines














!****************************************************************************************************************************************************
!I tried to construct these nested loops such that the inner loops have the most iterations while the outer loops have less.
!Doing this should help to increase performance.
DO vsini_loop = e, Number_vsini_mc
   k = 1
   vsini = vsini_array(vsini_loop)
   
   Number_vsini = Number_vsini + 1     !Increase Number_vsini by 1 each time this do loop runs.
   
   !probably won't be using this anymore.
   IF ((vsini < vsini_prior + (vsini_interval/10.0D0)) .AND. (vsini > vsini_prior - (vsini_interval/10.0D0)) .OR. (vsini_fit == 'Y')) THEN
      fit_flag9 = 'Y'
   ELSE
      fit_flag9 = 'N'
   END IF

!**********************************************************************************************************************************************************
   !$OMP PARALLEL DO
   DO CONCURRENT (stellar_rotation_angle_loop = k:Number_stellar_rotation_angle_mc)
   !DO stellar_rotation_angle_loop = k, Num_stellar_rotation_angle
      stellar_rotation_angle = stellar_rotation_angle_array(stellar_rotation_angle_loop)

      !probably won't be using this anymore.
      IF ((stellar_rotation_angle < stellar_rotation_angle_prior + (stellar_rotation_angle_interval/10.0D0)) &
         .AND. (stellar_rotation_angle > stellar_rotation_angle_prior - (stellar_rotation_angle_interval/10.0D0)) &
         .OR. (stellar_rotation_angle_fit == 'Y')) THEN
         fit_flag10 = 'Y'
         !PRINT *, "fit_flag Mp: ", fit_flag4
      ELSE
         fit_flag10 = 'N'
      END IF

      !probably won't be using this anymore.
      IF ((fit_flag9 == 'Y') .AND. (fit_flag10 == 'Y')) THEN
         fit_flag = 'Y'
      ELSE
         fit_flag = 'N'
      END IF

      Parameter_array_total(Number + stellar_rotation_angle_loop,1) = RV_offset_datasets
      Parameter_array_total(Number + stellar_rotation_angle_loop,2) = (Orbital_period_1 - Orbital_period_prior)*10000.0D0
      Parameter_array_total(Number + stellar_rotation_angle_loop,3) = (JD_time_mid_transit - JD_time_mid_transit_prior)*100.0D0
      Parameter_array_total(Number + stellar_rotation_angle_loop,4) = Mpp
      Parameter_array_total(Number + stellar_rotation_angle_loop,5) = Rpp
      Parameter_array_total(Number + stellar_rotation_angle_loop,6) = Ecc
      Parameter_array_total(Number + stellar_rotation_angle_loop,7) = Inc
      Parameter_array_total(Number + stellar_rotation_angle_loop,8) = omega_arg_periastron
      Parameter_array_total(Number + stellar_rotation_angle_loop,9) = RV_zero_offset
      !Vsini values at each iteration are placed in the Parameter_array at (Number,1).
      Parameter_array_total(Number + stellar_rotation_angle_loop,10) = vsini      
      !stellar_rotation_angle values at each iteration are placed in the Parameter_array at (number,2).
      Parameter_array_total(Number + stellar_rotation_angle_loop,11) = stellar_rotation_angle

      !For the Monte Carlo simulation, Parameter_array_fit is just for the vsini and lambda iterations.
      Parameter_array_fit(fit_count + stellar_rotation_angle_loop,1) = RV_offset_datasets
      Parameter_array_fit(fit_count + stellar_rotation_angle_loop,2) = (Orbital_period_1 - Orbital_period_prior)*10000.0D0
      Parameter_array_fit(fit_count + stellar_rotation_angle_loop,3) = (JD_time_mid_transit - JD_time_mid_transit_prior)*100.0D0
      Parameter_array_fit(fit_count + stellar_rotation_angle_loop,4) = Mpp
      Parameter_array_fit(fit_count + stellar_rotation_angle_loop,5) = Rpp
      Parameter_array_fit(fit_count + stellar_rotation_angle_loop,6) = Ecc
      Parameter_array_fit(fit_count + stellar_rotation_angle_loop,7) = Inc
      Parameter_array_fit(fit_count + stellar_rotation_angle_loop,8) = omega_arg_periastron
      Parameter_array_fit(fit_count + stellar_rotation_angle_loop,9) = RV_zero_offset
      Parameter_array_fit(fit_count + stellar_rotation_angle_loop,10) = vsini  
      Parameter_array_fit(fit_count + stellar_rotation_angle_loop,11) = stellar_rotation_angle
!        Keep track of which iterations will be used to determine best fits (not all parameters are used for the fit
!        as some parameters are iterated to propagate the errors in the parameters which are fitted).
      fit_array_indeces(fit_count + stellar_rotation_angle_loop) = Number + stellar_rotation_angle_loop
      
      !A counting array that helps to place values of the Parameter_array in numerical order and find 
      !values at a certain position in the Parameter_array. 
      Number_array_total(Number + stellar_rotation_angle_loop) = Number + stellar_rotation_angle_loop
      !PRINT *, 'Made it to line 1345'
      !      PRINT *, "Got here 1"

!      Allocate(Lightcurve(total_interval), RMcurve(total_interval), RV_theory(total_interval, 2), RVcurve(total_interval), & 
!               Timestep(total_interval), Transit_LC_array(total_interval, 2), Time_of_loop(total_interval), &
!               data_points_compare(total_interval), RM_effect_array(total_interval, 2), Planet_star_distance_array(total_interval, 2), &
!               Phase_array(total_interval), Phase_act_array(total_interval), Phase_array_n(total_interval), &
!               True_anomaly_array(total_interval), ecc_anomaly_array(total_interval), Xpos_array(total_interval), &
!               Ypos_array(total_interval), Zpos_array(total_interval))
               
      Allocate(RV_theory(total_interval, 2), data_points_compare(total_interval))
                              

      Time_transit_start = 0                         !Set the transit start time to zero.
      Transit_start_no_inc = 0                       !Set the transit start time for 90 inc to zero.
      Time_occultation_start = 0                     !Set the occultation start time to zero.
      Time_transit_end = 0                           !Set the transit end time to zero.
      Transit_end_no_inc = 0                         !Set the transit end time for 90 inc to zero.
      Time_occultation_end = 0                       !Set the occultation end time to zero.
      Time_mid_transit = 0                           !Set the mid transit time to zero.
      Time_mid_occultation = 0                       !Set the mid occultation time to zero.
      num_compare = 0
      Planet_star_distance = 0                       !Set the distance between center of the star to the center of the planet (orbital radius) to zero.      
      Transit_start_no_inc_position = 1
      transit_start_position = 1
      Occultation_start_position = 1
      transit_End_position = 1
      Occultation_End_position = 1
      transit_mid_position = 1
      occultation_mid_position = 1
      transit_End_position_no_inc = 1
      occultation_End_position = 1

      DO k1 = 1, total_interval
         data_points_compare(k1) = 'N'
      END DO
      
!      PRINT *, "*****************************************************************************************************"
      
!      PRINT *, 'fit_flag0: ', fit_flag0, '   fit_flag1: ', fit_flag1, '   fit_flag2: ', fit_flag2, '   fit_flag3: ', fit_flag3
      !PRINT *, 'fit_flag1: ', fit_flag1, '   fit_flag2: ', fit_flag2, '   fit_flag3: ', fit_flag3, '   fit_flag4: ', fit_flag4
      !PRINT *, 'fit_flag5: ', fit_flag5, '   fit_flag6: ', fit_flag6, '   fit_flag7: ', fit_flag7, '   fit_flag8: ', fit_flag8
!      PRINT *, 'fit_flag4: ', fit_flag4, '   fit_flag5: ', fit_flag5, '   fit_flag6: ', fit_flag6, '   fit_flag7: ', fit_flag7, '   fit_flag8: ', fit_flag8
      
!      PRINT *, "array_resorted? ", array_resorted
      
!      PRINT *, "Orbital_period_1: ", Orbital_period_1
!      PRINT *, "Orbital period Number: ", Number_Orbital_period
      
!      PRINT *, "JD_time_mid_transit: ", JD_time_mid_transit
!      PRINT *, "JD time mid transit Number: ", Number_JD_time_mid_transit
      
!      PRINT *, "Mass of planet: ", Mpp
!      PRINT *, "Mp Number: ", Number_Mp
!      PRINT *, "Semi-major axis: ", Rorb/AU
      
!      PRINT *, "Radius of planet: ", Rpp
!      PRINT *, "Rp Number: ", Number_Rp
      
!      PRINT *, "Eccentricity of planet: ", Ecc
!      PRINT *, "Ecc Number: ", Number_Ecc
      
!      PRINT *, "Inclination of planet: ", Inc
!      PRINT *, "Inc Number: ", Number_Inc
      
!      PRINT *, "omega_arg_periastron: ", omega_arg_periastron
!      PRINT *, "omega_arg_periastron Number: ", Number_omega_arg_periastron
      
!      PRINT *, "RV_zero_offset: ", RV_zero_offset
!      PRINT *, "RV_zero_offset Number: ", Number_RV_zero_offset
      
!      PRINT *, "RV_offset_datasets: ", RV_offset_datasets
!      PRINT *, "RV_offset_datasets Number: ", Number_RV_offset_datasets

!******************************************************************************************************************************************* 

      DO Time_loop = 1, total_interval

         !time = Time_start + ((Time_loop - 1.0D0)*Time_plot_interval)
         !The time will come from the first data point which must be normalized as a fraction of a day before mid transit.
         time = time_data_array(Time_loop)

         !Set transit flag to N when outside of transit event.
         transit_flag = 'N'
         !Set occultation flag to N when outside of secondary transit event.
         occultation_flag = 'N'
         time_compare_flag = 'N'                   

         IF (Ecc == 0) THEN
            ecc_anomaly = ((2.0D0*pi)/Orbital_period)*(time - time_peri_passage)
         ELSE
            !Calculate the eccentric anomaly using the functions found near the end of the program.
            sum_ecc_anomaly = 0
            ecc_anomaly_before = 1000000
            DO order = 1, 20
               !Calculate the value of the Bessel function which is used to find the eccentric anomaly.
               !Bessel_value = BESSJ(order, order*Ecc)
               !Fortran 2008 Bessel function.
               Bessel_value = bessel_jn(order, order*Ecc)
               !PRINT *, "Bessel_value = ", Bessel_value
               IF (order > 1) THEN
                  ecc_anomaly_before = sum_ecc_anomaly
               END IF
               sum_ecc_anomaly = sum_ecc_anomaly + ((2.0D0/order)*Bessel_value*sin(order*(((2.0D0*pi) &
                                 /Orbital_period)*(time - time_peri_passage))))
               !PRINT *, "sum_ecc_anomaly = ", sum_ecc_anomaly
               ecc_anomaly_after = sum_ecc_anomaly
               IF ((order > 1) .AND. (ABS(ecc_anomaly_after - ecc_anomaly_before) <= Bessel_function_exit)) EXIT 
            END DO
            ecc_anomaly = (((2.0D0*pi)/Orbital_period)*(time - time_peri_passage)) + sum_ecc_anomaly
         END IF

         !Now use the ecc_anomaly to determine the True_anomaly and Phase_angle for a given data point.
         True_anomaly = 2.0D0*(atan(tan(ecc_anomaly/2.0D0)*(sqrt((1.0D0 + Ecc)/(1.0D0 - Ecc)))))

         IF (Ecc == 0) THEN
            !The time of the similation for a specific true_anomaly with the mid transit time equal to 0.
            Time_check = ((ecc_anomaly*Orbital_period)/(2.0D0*pi)) + time_peri_passage
            !The distance between the center of the planet to the center of the star in a circular orbit.
            !Planet_star_distance = Rorb
            Planet_star_distance = Rorb + Rorb_star
         ELSE
            !The time of the similation for a specific true_anomaly.
            Time_check = (((ecc_anomaly - (Ecc*sin(ecc_anomaly)))*Orbital_period)/(2.0D0*pi)) + time_peri_passage         
            !The distance between the center of the planet to the center of the star in an eccentric orbit.
            !Planet_star_distance = (Rorb*(1.0D0 - Ecc**2.0D0))/(1.0D0 + (Ecc*cos(True_anomaly)))
            Planet_star_distance = ((Rorb + Rorb_star)*(1.0D0 - Ecc**2.0D0))/(1.0D0 + (Ecc*cos(True_anomaly)))
         END IF
         !Time in reference to mid transit time
         Time_ref = Time - Time_transit

         !The position of the planet on the x-axis.
         Xpos = Planet_star_distance*((-sin(True_anomaly)*sin((omega_arg_periastron)*(pi/180.0D0))) &
                + (cos(True_anomaly)*cos((omega_arg_periastron)*(pi/180.0D0))))
         !The position of the planet on the y-axis.
         Ypos = Planet_star_distance*((cos(True_anomaly + pi)*cos((Inc)*(pi/180.0D0))*sin((omega_arg_periastron)*(pi/180.0D0))) &
                + (sin(True_anomaly + pi)*cos((Inc)*(pi/180.0D0))*cos((omega_arg_periastron)*(pi/180.0D0))))
         !The position of the planet on the z-axis.
         Zpos = Planet_star_distance*((-cos((omega_arg_periastron)*(pi/180.0D0))*sin((Inc)*(pi/180.0D0))*sin(True_anomaly)) &
                - (cos(True_anomaly)*sin((Inc)*(pi/180.0D0))*sin((omega_arg_periastron)*(pi/180.0D0))))
         Dist2 = Xpos**2.0D0 + Ypos**2.0D0                   !Square of the planet-star apparent seperation.
         Distance_center = sqrt(Dist2)               !Apparent seperation between the planet and the star.
         Lblocked = 0.0D0
         Lblocked2 = 0.0D0                           !A variable for the Anomalous velocity equation.
         v_rm = 0.0D0                                !Anomalous velocity of each pixel set to zero.
         Total_RM = 0.0D0

         IF (Xpos <= 0 .AND. Zpos >= 0) THEN
            !Planet is currently in quadrant three so add pi.
            phase_angle = atan(Xpos/Zpos) + pi
            !Taking into account orbital inclination.
            Phase_angle_observed = atan(-(sqrt(Xpos**2.0D0 + Ypos**2.0D0))/Zpos) + pi
         ELSE IF (Xpos >= 0 .AND. Zpos >= 0) THEN
            !Planet is currently in quadrant four so add pi.
            phase_angle = atan(Xpos/Zpos) + pi
            !Taking into account orbital inclination.
            Phase_angle_observed = atan((sqrt(Xpos**2.0D0 + Ypos**2.0D0))/Zpos) + pi
         ELSE IF (Xpos >= 0 .AND. Zpos <= 0) THEN
            !Planet is currently in quadrant one so add 2pi.
            phase_angle = 2.0D0*pi + atan(Xpos/Zpos)
            !Taking into account orbital inclination.
            Phase_angle_observed = 2.0D0*pi + atan((sqrt(Xpos**2.0D0 + Ypos**2.0D0))/Zpos)
         ELSE IF (Xpos <= 0 .AND. Zpos <= 0) THEN
            !Planet is currently in quadrant two so add 2pi.
            phase_angle = atan(Xpos/Zpos)
            !Taking into account orbital inclination.
            Phase_angle_observed = atan(-(sqrt(Xpos**2.0D0 + Ypos**2.0D0))/Zpos)
         END IF 

         True_phase = acos(sin(True_anomaly + (omega_arg_periastron*(pi/180.0D0)))*sin(Inc*(pi/180.0D0)))
         Phase_orbit_n = acos(sin(True_anomaly + (omega_arg_periastron*(pi/180.0D0))))
      
         !If the planet is neither infront of the star (transit) or behind the star (occultation), then calculate the flux being reflected 
         !off the surface of the exoplanet based on its bond albedo, radius, phase angle, etc.
         IF (Distance_center > (Rs + Rp)) THEN
            Aplan = pi*Rp2                                               !Surface area of the planet.
            Total_L = 1.0D0 + ((Aplan/(4.0D0*pi*Planet_star_distance**2.0D0))*Albedo*(0.5D0*(1.0D0+cos(True_phase))))
         END IF

         !Hopefully this will ensure that the same data points are compared despite the length of the transit changing due to verying the parameters.
         IF ((Time_ref >= -length_time_compare) .AND. (Time_ref <= length_time_compare)) THEN
            time_compare_flag = 'Y'
            !PRINT *, "Time_ref", Time_ref
            !PRINT *, "length_time_compare", length_time_compare
         END IF

         IF ((Distance_center <= (Rs + Rp)) .AND. (Zpos > 0)) THEN            
            !If the seperation between the disk of the star and planet are less then Rs + Rp and Zpos is positive then the transit begins.
            transit_flag = 'Y'                  !Planet is transiting. Will be used to determine when to do the chi squared analysis.
            
      
            IF ((Rp2/Rs2) >= 0.030) THEN
               Radius_planet_array = M_pixels / 2.0D0         !Radius of planet in the pixel array. 
               Center_of_planet_x = M_pixels / 2.0D0          !The center of the planet on the X-axis.
               Center_of_planet_y = M_pixels / 2.0D0          !The center of the planet on the Y-axis.
               Pixel = Rp/Radius_planet_array               !The number of meters per pixel.
               Area_pixel = Pixel**2.0D0                        !The area of each pixel.
               Io_Pixel = Io * Area_pixel                   !Variable to speed up calculations (also physically represents the
                                                      !luminosity of the brightest pixel).    

               DO i = 1, M_pixels
               !DO CONCURRENT (i = 1: M_pixels, j = 1: M_pixels)
                  X_pixel = (pixel * i) + (Xpos - Rp)           !Calculates the location of the pixel on the x-axis.
                  X_pixel2 = X_pixel**2.0D0
                  XpXp = X_pixel * Xpos                         !temporary var for speed calculation
                  DO j = 1, M_pixels
                     Y_pixel = (pixel * j) + abs(Ypos) - Rp        !Calculates the location of the pixel on the y-axis.
                     !Calculates the location of the pixel along the x-axis of the rotation axis of the star.
                     X_pixel_prime = (X_pixel*cos((stellar_rotation_angle*pi)/180.0D0)) + &
                                     (Y_pixel*sin((stellar_rotation_angle*pi)/180.0D0))     
                     Dist_center_pixel = X_pixel2 + Y_pixel**2.0D0         !squared distance of pixel from star
                     !Calculates the values of the pixels according to how far away they are from the center of the star and the plane.
                     !squared distance of pixel from planet using a limb darkening equation.
                     Dist_planet_pixel = Dist_center_pixel  - (2.0D0*(XpXp + (Y_pixel*Ypos))) + Dist2
                     Sub_planet_velocity = vsini*(X_pixel_prime/Rs)
                     IF ((Dist_center_pixel <= Rs2) .AND. (Dist_planet_pixel <= Rp2)) THEN
                        IF (linear_quadratic == 'q') THEN
                           Lblocked2 = Io_Pixel*(1.0D0-q_1*(1.0D0-sqrt(abs(1.0D0-(Dist_center_pixel/Rs2)))) - &
                                       q_2*(1.0D0 - sqrt(abs(1.0D0-(Dist_center_pixel/Rs2))))**2.0D0)          !Quadratic limb darkening equation.
                        ELSE
                           Lblocked2 = Io_Pixel*(1.0D0-u*(1.0D0-sqrt(abs(1.0D0-(Dist_center_pixel/Rs2)))))     !First order limb darkening equation.          
                        END IF
                        !Lblocked2 = Io_Pixel*(1.0D0-u*(1.0D0-sqrt(1.0D0-(Dist_center_pixel/Rs2))))          !First order limb darkening equation.
                        Lblocked = Lblocked + Lblocked2                                         !First order limb darkening equation.
                        v_rm = v_rm - ((Lblocked2*Sub_planet_velocity)*((((2.0D0*vturb**2.0D0)+(2.0D0*vsini**2.0D0))/((2.0D0*vturb**2.0D0) &
                               + vsini**2.0D0))**(3.0D0/2.0D0))*(1.0D0-((Sub_planet_velocity**2.0D0)/((2.0D0*vturb**2.0D0)+(vsini**2.0D0)))))   
                        !Anomalous velocity of each pixel.
                     END IF
                  END DO
               END DO
               Total_L = 1.0D0 - Lblocked
               Total_RM = 0.0D0 + v_rm                       !Total anomalous velocity for all the pixels.

            ELSE IF ((Rp2/Rs2) <= 0.030) THEN
  
               !Calculates the location of the center of the planet along the x-axis of the rotation axis of the star.
               X_prime = (Xpos*cos((stellar_rotation_angle*pi)/180.0D0)) &
                         + (Ypos*sin((stellar_rotation_angle*pi)/180.0D0))     
               set_distance_center = Distance_center         !The limb darkening equation will use this distance as long as the center of the 
                                                             !planet is inside the radius of the star.
               Sub_planet_velocity = vsini*(X_prime/Rs)      !Calculate the subplanetary velocity (the stellar velocity blocked by the 
                              !planetary disc) from vsini times the distance from the center of the planet to 
                             !the stellar rotation axis divided by the radius of the star.
               Io_planet = 0.0D0
               
               !Start the planet on the x-axis so that the planet is just far enough away not to touch the disk of the star. 
               IF ((Distance_center <= (Rs + Rp)) .AND. (Distance_center >= (Rs - Rp))) THEN    
                 
                  dist_cent1_int = (Distance_center**2.0D0 + Rs2 - Rp2)/(2.0D0*Distance_center)
                  !Location on the x-axis for the first intersection point.
                  X_int_1 = ((Xpos*dist_cent1_int)/Distance_center) &
                            + ((Ypos/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0)) 
                  !Location on the y-axis for the first intersection point.
                  Y_int_1 = ((Ypos*dist_cent1_int)/Distance_center) &
                            - ((Xpos/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0))            
                  !Location on the x-axis for the second intersection point.
                  X_int_2 = ((Xpos*dist_cent1_int)/Distance_center) &
                            - ((Ypos/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0))             
                  !Location on the y-axis for the second intersection point.
                  Y_int_2 = ((Ypos*dist_cent1_int)/Distance_center) &
                            + ((Xpos/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0)) 
                  !The limb darkening equation will use this distance if any part of the disc of the planet is outside the radius of the 
                  !star.
                  !This is the distance between the center of the star to the center of the area inside the star blocked by the planet.           
                  set_distance_center = ((abs(Distance_center) - Rp) + Rs)/2.0D0
                  !Next the program calculates the length between the two intersection point.
                  Length = sqrt((X_int_1 - X_int_2)**2.0D0 + (Y_int_1 - Y_int_2)**2.0D0)
                  !Calculate the angle between the Y position of the center of planet and the x position of the center of planet.
                  !This is used to determine the center of the area inside the star blocked by the planet in the stellar rotational
                  !axis coordinate system.
                  Theta_inside = atan(Ypos/Xpos)
                  !Calculates the distance from the center of the star (along x axis) to the center of the area inside the star blocked by 
                  !the planet.
                  !the pi factor added in this equation guarantees that Xpos_in has the correct sign.
                  Xpos_in = set_distance_center*cos(Theta_inside + pi)
                  Xpos_inside = Xpos_in

                  !This makes sure Xpos_inside has the correct sign.
                  IF (Xpos >= 0) THEN
                     Xpos_inside = -Xpos_in
                  END IF

                  !Calculates the distance from the center of the star (along y axis) to the center of the area inside the star blocked by 
                  !the planet.
                  Ypos_inside = abs(set_distance_center*sin(Theta_inside))
                  !Changes the x-coordinate to the stellar rotation axis by an angle formed between the orbital plane of the planet and the 
                  !stellar roatation plane of the star.
                  x_prime_distance = (Xpos_inside*cos((stellar_rotation_angle*pi)/180.0D0)) &
                                     + (Ypos_inside*sin((stellar_rotation_angle*pi)/180.0D0))
                  Sub_planet_velocity = vsini*(x_prime_distance/Rs)    !Calculate the subplanetary velocity (the stellar velocity blocked by 
                  !the planetary disc) from vsini times the distance from the center of the planet to the stellar rotation axis divided by 
                  !the radius of the star.
                  Beta1 = 2.0D0*asin((0.5D0*Length)/Rp)         !Angle used to calculate the partial area of the planet in of the star.
                  Alpha1 = 2.0D0*asin((0.5D0*Length)/Rs)        !Angle used to calculate the partial area of the planet in of the star.
                    
                  IF (Distance_center >= sqrt(Rs2 - Rp2)) THEN   
                     !The surface area of the planet when the center is outside the disk of the star.        
                     Aplan = (0.5D0 * Rs2 * (Alpha1 - sin(Alpha1))) + (0.5D0 * Rp2 * (Beta1 - sin(Beta1)))
                     !Normalized area of the planet with the given limb darkening function
                     Io_planet = Io * Aplan
                  END IF
     
                  IF (Distance_center < sqrt(Rs2 - Rp2)) THEN
                     Aplan = ((pi*Rp2 + (0.5D0 * Rs2 * (Alpha1 - sin(Alpha1)))) - (0.5D0 * Rp2 * (Beta1 - sin(Beta1))))
                     !The surface area of the planet when the center is inside the disk of the star.
                     Io_planet = Io * Aplan
                  END IF
               END IF

               IF (Distance_center < (Rs - Rp)) THEN
                  Aplan = (pi*Rp2)                        !The surface area of the planet.
                  Io_planet = Io * Aplan
               END IF
  
               !The ratio of the area of the planet blocking the star to the area of the star utilizing the first order limb darkening 
               !equation.
               IF (linear_quadratic == 'q') THEN
                  Lblocked = Io_planet*(1.0D0-q_1*(1.0D0-sqrt(abs(1.0D0-(set_distance_center**2.0D0/Rs2)))) - &
                              q_2*(1.0D0 - sqrt(abs(1.0D0-(set_distance_center**2.0D0/Rs2))))**2.0D0)          !Quadratic limb darkening equation.
               ELSE
                  Lblocked = (Io_planet)*(1.0D0-(u*(1.0D0-sqrt(abs(1.0D0-(set_distance_center**2.0D0/Rs2))))))  
               END IF
               !Lblocked = (Io_planet)*(1.0D0-(u*(1.0D0-sqrt(abs(1.0D0-(set_distance_center**2.0D0/Rs2))))))                                  
               v_rm = - ((Lblocked*Sub_planet_velocity)*((((2.0D0*vturb**2.0D0)+(2.0D0*vsini**2.0D0))/((2.0D0*vturb**2.0D0) &
                      + vsini**2.0D0))**(3.0D0/2.0D0))*(1.0D0-((Sub_planet_velocity**2.0D0)/((2.0D0*vturb**2.0D0)+(vsini**2.0D0)))))
               Total_L = 1.0D0 - Lblocked                                             !Total amount of light blocked by the planet.
               Total_RM = 0.0D0 + v_rm                                                  !Total anomalous velocity.
            END IF

         END IF

         !If the seperation between the disk of the star and planet are less then Rs + Rp and Zpos is negative then the secondary transit 
         !(occulation) begins.
         IF ((Distance_center <= (Rs + Rp)) .AND. (Zpos < 0)) THEN     

            occultation_flag = 'Y'                  !Planet is occulting.       

            IF (Distance_center >= (Rs - Rp)) THEN           !Partial secondary transit.             
               dist_cent1_int = (Distance_center**2.0D0 + Rs2 - Rp2)/(2.0D0*Distance_center)
               Ypos2 = abs(Ypos)   !Set Ypos to a positive value to avoid discontinuity.  
               !Location on the x-axis for the first intersection point.
               X_int_1 = ((Xpos*dist_cent1_int)/Distance_center) + ((Ypos2/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0)) 
               !Location on the y-axis for the first intersection point.
               Y_int_1 = ((Ypos2*dist_cent1_int)/Distance_center) - ((Xpos/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0))            
               !Location on the x-axis for the second intersection point.
               X_int_2 = ((Xpos*dist_cent1_int)/Distance_center) - ((Ypos2/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0))             
               !Location on the y-axis for the second intersection point.
               Y_int_2 = ((Ypos2*dist_cent1_int)/Distance_center) + ((Xpos/Distance_center)*sqrt(Rs2 - dist_cent1_int**2.0D0)) 
               !The limb darkening equation will use this distance if any part of the disc of the planet is outside the radius of the 
               !star. This is the distance between the center of the star to the center of the area inside the planet blocked by the star.           
                 
               set_distance_center = ((abs(Distance_center) - Rp) + Rs)/2.0D0
               !Next the program calculates the length between the two intersection point.
               Length = sqrt((X_int_1 - X_int_2)**2.0D0 + (Y_int_1 - Y_int_2)**2.0D0)
               !Calculate the angle between the Y position of the center of planet and the x position of the center of planet.
               Theta_inside = atan(Ypos2/Xpos)
               !Calculates the distance from the center of the star (along x axis) to the center of the area inside the star blocked by 
               !the planet.
               !the pi factor added in this equation guarneetes that Xpos_in has the correct sign.
               Xpos_in = set_distance_center*cos(Theta_inside + pi)
               Xpos_inside = Xpos_in
         
               !This makes sure Xpos_inside has the correct sign.
               IF (Xpos >= 0) THEN
                  Xpos_inside = -Xpos_in
               END IF

               !Calculates the distance from the center of the star (along y axis) to the center of the area inside the star blocked by 
               !the planet.
               Ypos_inside = set_distance_center*sin(Theta_inside)
               Beta1 = 2.0D0*asin((0.5D0*Length)/Rp)         !Angle used to calculate the partial area of the planet in the star.
               Alpha1 = 2.0D0*asin((0.5D0*Length)/Rs)        !Angle used to calculate the partial area of the planet in the star.
                    
               IF (Distance_center >= sqrt(Rs2 - Rp2)) THEN  
                  !The surface area of the planet that is not behind the star when the center of the planets disk is visible.        
                  Aplan = (pi*Rp2) - ((0.5D0 * Rp2 * (Beta1 - sin(Beta1))) + (0.5D0 * Rs2 * (Alpha1 - sin(Alpha1))))
               END IF
     
               IF (Distance_center < sqrt(Rs2 - Rp2)) THEN
                  !The surface area of the planet that is not behind the star when the center of the planets disk is behind the star.
                  Aplan = (0.5D0 * Rs2 * (Beta1 - sin(Beta1)))
               END IF
            END IF
  
            IF (Distance_center < (Rs - Rp)) THEN
               Aplan = 0.0D0                      !No flux coming from planet since it's behind the star.
            END IF
  
            !Total amount of light blocked by the planet.
            Total_L = 1.0D0 + ((Aplan/(4.0D0*pi*Planet_star_distance**2.0D0))*Albedo*(0.5D0*(1.0D0 + cos(True_phase))))
         END IF

         IF (Ecc == 0) THEN
            !The radial velocity of the star which includes adding the RM effect for a circular orbit.
            RV = RVamplitude*(cos(pi + True_anomaly + ((omega_arg_periastron)*(pi/180.0D0)))) + Total_RM
         ELSE
            !The radial velocity of the star which includes adding the RM effect for an eccentric orbit.
            RV = RVamplitude*(cos(pi + True_anomaly + ((omega_arg_periastron)*(pi/180.0D0))) &
                 + (Ecc*cos((omega_arg_periastron)*(pi/180.0D0) + pi))) + Total_RM
         END IF

         !If transit flag is Y, then indicate the data point to compare.
         IF ((transit_flag == 'Y') .OR. (time_compare_flag == 'Y')) THEN
            num_compare = num_compare + 1
            data_points_compare(Time_loop) = 'Y'
         END IF                   

         RV_theory(Time_loop,1) = Time_ref
         RV_theory(Time_loop,2) = RV

      END DO 
         
      !*******************************************************************************************************************************

      Number_points1 = 1
      chi_2 = 0          !Set chi squared equal to zero.
      model_data_difference = 0        

      DO s = 1, total_interval
         IF (data_points_compare(s) == 'Y') THEN
            !Now do chi squared analysis but only for data points during transit.
            !Now calculate chi squared from the adjusted RV model.
            !First calculate sigma squared using N - 1 (sample variance rather than theoretical variance).

            !chi_2 = chi_2 + ((RV_offset_data_array(s,2) - RV_Theory(s,2))/RV_offset_data_array(s,3))**2.0D0
            chi_2 = chi_2 + ((RV_synthetic_data_array(s,2) - RV_Theory(s,2))/RV_synthetic_data_array(s,3))**2.0D0
            !Model_theory_difference_array(Number_points1) = ABS(RV_offset_data_array(s,2) - RV_Theory(s,2))
            !model_data_difference = model_data_difference + ABS(RV_offset_data_array(s,2) - RV_Theory(s,2))
            model_data_difference = model_data_difference + ABS(RV_synthetic_data_array(s,2) - RV_Theory(s,2))

            Number_points1 = Number_points1 + 1              !Counter to determine the number of elements used to calculate sigma squared.
         END IF
      END DO

      !PRINT *,"Chi squared: ", chi_2
      r_Chi_2 = chi_2/(Number_points1 - Number_fit - 1)
      !PRINT *,"Reduced Chi squared: ", r_Chi_2
      
      !PRINT *, 'RV_offset_datasets: ', RV_offset_datasets
      !PRINT *, '(Orbital_period_1 - Orbital_period_prior)*10000.0D0: ', (Orbital_period_1 - Orbital_period_prior)*10000.0D0
      !PRINT *, 'Orbital_period_1: ', Orbital_period_1
      !PRINT *, '(JD_time_mid_transit - JD_time_mid_transit_prior)*100.0D0: ', (JD_time_mid_transit - JD_time_mid_transit_prior)*100.0D0
      !PRINT *, 'JD_time_mid_transit: ', JD_time_mid_transit
      !PRINT *, 'Mpp: ', Mpp
      !PRINT *, 'Rpp: ', Rpp
      !PRINT *, 'Ecc: ', Ecc
      !PRINT *, 'Inc: ', Inc
      !PRINT *, 'omega_arg_periastron: ', omega_arg_periastron
      !PRINT *, 'RV_zero_offset: ', RV_zero_offset
      !PRINT *, 'vsini: ', vsini
      !PRINT *, 'stellar_rotation_angle: ', stellar_rotation_angle

      Chi_squared_array_total(Number + stellar_rotation_angle_loop) = chi_2                           !Puts the chi squared values into an array.
      Reduced_Chi_Squared_array_total(Number + stellar_rotation_angle_loop) = r_Chi_2                 !Puts the reduced chi squared values into an array.
              
      Chi_squared_array_fit(fit_count  + stellar_rotation_angle_loop) = chi_2                           !Puts the chi squared values into an array.
      Reduced_Chi_Squared_array_fit(fit_count  + stellar_rotation_angle_loop) = r_Chi_2                 !Puts the reduced chi squared values into an array.
      
      DEALLOCATE(RV_theory)
      DEALLOCATE(data_points_compare)

!      STOP
!      PRINT *, 'Made it to line 1805'

   END DO
   !$OMP END PARALLEL DO
   
   !Moved counter outside of parallel loop. Each time the above do loop runs through all its iterations, 
   !number will increase by number of loops that have been iterated which is k.
   number = number + Number_stellar_rotation_angle_mc
   Number_stellar_rotation_angle = Number_stellar_rotation_angle + Number_stellar_rotation_angle_mc
   fit_count = fit_count + Number_stellar_rotation_angle_mc
   
   stellar_rotation_angle_flag = 'T'
   
   Number_iterations_left_total = Number_iterations_left_total - Number_stellar_rotation_angle_mc         !Let's the user know the progress of the program.
   
   !Iteration_number = number + Num_stellar_rotation_angle                           !Let's the user know the progress of the program.  
   Iteration_number = number                           !Let's the user know the progress of the program.  
   
   !PRINT *, 'Made it to line 1823'
   !PRINT *, 'Iteration_number: ', Iteration_number
   !PRINT *, 'fit_count: ', fit_count
   
END DO
     
!PRINT *, 'RV_offset_datasets: ', RV_offset_datasets
!PRINT *, '(Orbital_period_1 - Orbital_period_prior)*10000.0D0: ', (Orbital_period_1 - Orbital_period_prior)*10000.0D0
!PRINT *, 'Orbital_period_1: ', Orbital_period_1
!PRINT *, '(JD_time_mid_transit - JD_time_mid_transit_prior)*100.0D0: ', (JD_time_mid_transit - JD_time_mid_transit_prior)*100.0D0
!PRINT *, 'JD_time_mid_transit: ', JD_time_mid_transit
!PRINT *, 'Mpp: ', Mpp
!PRINT *, 'Rpp: ', Rpp
!PRINT *, 'Ecc: ', Ecc
!PRINT *, 'Inc: ', Inc
!PRINT *, 'omega_arg_periastron: ', omega_arg_periastron
!PRINT *, 'RV_zero_offset: ', RV_zero_offset
!PRINT *, 'vsini: ', vsini
!PRINT *, 'stellar_rotation_angle: ', stellar_rotation_angle     
     

!PRINT *, 'fit_flag9: ', fit_flag9
!PRINT *, 'fit_flag10: ', fit_flag10

PRINT *, "array_resorted? ", array_resorted

PRINT *, "Orbital_period_1: ", Orbital_period_1
!PRINT *, "Orbital period Number: ", Number_Orbital_period - 1

PRINT *, "JD_time_mid_transit: ", JD_time_mid_transit
!PRINT *, "JD time mid transit Number: ", Number_JD_time_mid_transit - 1

PRINT *, "Mass of planet: ", Mpp
!PRINT *, "Mp Number: ", Number_Mp - 1
PRINT *, "Semi-major axis: ", Rorb/AU
PRINT *, "Stellar Semi-major axis: ", Rorb_star/AU

PRINT *, "Radius of planet: ", Rpp
!PRINT *, "Rp Number: ", Number_Rp - 1

PRINT *, "Eccentricity of planet: ", Ecc
!PRINT *, "Ecc Number: ", Number_Ecc - 1

PRINT *, "Inclination of planet: ", Inc
!PRINT *, "Inc Number: ", Number_Inc - 1

PRINT *, "omega_arg_periastron: ", omega_arg_periastron
!PRINT *, "omega_arg_periastron Number: ", Number_omega_arg_periastron - 1

PRINT *, "RV_zero_offset: ", RV_zero_offset
!PRINT *, "RV_zero_offset Number: ", Number_RV_zero_offset - 1

PRINT *, "RV_offset_datasets: ", RV_offset_datasets
!PRINT *, "RV_offset_datasets Number: ", Number_RV_offset_datasets - 1    

PRINT *,"There are ", Number_iterations_left_total," interations left out of ", Number_iterations_total 

PRINT *,"Iteration number is ", Iteration_number
PRINT *,"Iteration fit number is ", fit_count
     
percent_processed = (100.0 * Iteration_number)/Number_iterations_total
string_bar = ''
string_bar_num = 1
IF ((percent_processed < 4.0) .AND. (percent_processed > 0.1)) THEN
   IF ((percent_processed - FLOOR(percent_processed) < 0.2499999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.0000001D0)) THEN
      string_bar(string_bar_num:) = ':'
   END IF
   IF ((percent_processed - FLOOR(percent_processed) < 0.4999999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.2499999D0)) THEN
      string_bar(string_bar_num:) = '/'
   END IF
   IF ((percent_processed - FLOOR(percent_processed) < 0.7499999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.4999999D0)) THEN
      string_bar(string_bar_num:) = '-'
   END IF
   IF ((percent_processed - FLOOR(percent_processed) < 0.9999999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.7499999D0)) THEN
      string_bar(string_bar_num:) = '\'
   END IF
END IF

IF (percent_processed >= 4) THEN
   DO percent_processed_loop = 1, FLOOR(percent_processed/4.0)
      string_bar(string_bar_num:) = '|'
      string_bar_num = string_bar_num + 1
      IF ((percent_processed < 99.9) .AND. (percent_processed_loop == FLOOR(percent_processed/4.0))) THEN
         IF ((percent_processed - FLOOR(percent_processed) < 0.2499999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.0000001D0)) THEN
            string_bar(string_bar_num:) = ':'
         END IF
         IF ((percent_processed - FLOOR(percent_processed) < 0.4999999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.2499999D0)) THEN
            string_bar(string_bar_num:) = '/'
         END IF
         IF ((percent_processed - FLOOR(percent_processed) < 0.7499999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.4999999D0)) THEN
            string_bar(string_bar_num:) = '-'
         END IF
         IF ((percent_processed - FLOOR(percent_processed) < 0.9999999D0) .AND. (percent_processed - FLOOR(percent_processed) >= 0.7499999D0)) THEN
            string_bar(string_bar_num:) = '\'
         END IF
      END IF
   END DO
END IF    
PRINT 1000, string_bar, FLOOR(percent_processed)
1000 FORMAT("Processing data: [", A25, "]", 1X, I3, "% complete")
IF (percent_processed < 4.0) THEN
   PRINT *, 'Estimated time remaining: ', 'calculating....'
ENDIF
IF (percent_processed >= 4.0) THEN
!   IF ((modulo(percent_processed,4.0) <= 0.1D0) .AND. (modulo(percent_processed,4.0) >= -0.1D0) .OR. (FLOOR(percent_processed) - FLOOR(percent_count*4) >= 0)) THEN
   IF (FLOOR(percent_processed) - (percent_count*4) >= 0) THEN
      call cpu_time(cpu_time_it2)
      time_remaining = ((cpu_time_it2 - cpu_time_it1)/percent_processed)*(100.0D0 - percent_processed)
      percent_count = percent_count + 1
   END IF
   WRITE (*,1001) 'Estimated time remaining: ', time_remaining/60.0D0, ' minutes'
   1001 FORMAT(A30,F8.2,A8)
END IF

!Let's find the minimum chi-square value for each Monte Carlo iteration. Then record the best vsini and lambda
!corresponding to the minimum chi-square as well as the 1 sigma uncertainty. The minimum chi-square value and
!1 sigma uncertainty will be used to weight the vsini and lambda values and derive the overall best vsini and
!lambda from the whole Monte Carlo simulation as well as there uncertainties.

PRINT *," "
PRINT *,"***********************************************************************************************"
PRINT *," "

!STOP

!Determine the lowest value for Chi square, reduced Chi square, and where these occur
!for the full parameter space and for just the parameters that will be used for the fit. 
!Also find the parameters which match these lowest values.

!Best fit vsini and lambda with other parameters at their prior values.
!DO l=1, Number_iterations_fit
!   print *, "Chi_squared_array_fit: ", Chi_squared_array_fit(l)
!END DO

IF (Monte_carlo_loop == mc_sample_size) THEN
   min_chi_squared_fit_prior = minval(Chi_squared_array_fit)
   print *, "mininum chi squared value for other parameters fixed to their prior values: ", min_chi_squared_fit_prior
   loc_min_chi_squared_fit_prior = minloc(Chi_squared_array_fit)
   !print *, "mininum chi squared location for other parameters fixed to their prior values: ", loc_min_chi_squared_fit_prior
   min_reduced_chi_squared_fit_prior = minval(Reduced_Chi_Squared_array_fit)
   !print *, "mininum reduced chi squared value for other parameters fixed to their prior values: ", min_reduced_chi_squared_fit_prior
   loc_min_reduced_chi_squared_fit_prior = minloc(Reduced_Chi_Squared_array_fit)
   !print *, "mininum reduced chi squared location for other parameters fixed to their prior values: ", loc_min_reduced_chi_squared_fit_prior

   IF (loc_min_chi_squared_fit_prior(1) /= loc_min_reduced_chi_squared_fit_prior(1)) THEN
      PRINT *, "The location of minimum value of chi squared and reduced chi squared not the same. Something went wrong!"
      PRINT *, "Minimum value location of chi squared prior: ", loc_min_chi_squared_fit_prior(1)
      PRINT *, "Minimum value location of reduced chi squared prior: ", loc_min_reduced_chi_squared_fit_prior(1)
   END IF

   min_chi_squared_mc(Monte_carlo_loop) = min_chi_squared_fit_prior
   min_reduced_chi_squared_mc(Monte_carlo_loop) = min_reduced_chi_squared_fit_prior

   best_vsini_fit_prior = Parameter_array_fit(loc_min_chi_squared_fit_prior(1),10)
   PRINT *, "The best vsini value corresponding to the lowest value of chi squared for other parameters fixed to their prior values: ", best_vsini_fit_prior, "m/s"
   best_spin_orbit_fit_prior = Parameter_array_fit(loc_min_chi_squared_fit_prior(1),11)
   PRINT *, "The best spin-orbit angle corresponding to the lowest value of chi squared for other parameters fixed to their prior values: ", best_spin_orbit_fit_prior, "deg"
   
ELSE
   min_chi_squared_fit = minval(Chi_squared_array_fit)
   print *, "mininum chi squared value for one MC iteration: ", min_chi_squared_fit
   loc_min_chi_squared_fit = minloc(Chi_squared_array_fit)
   print *, "mininum chi squared location for one MC iteration: ", loc_min_chi_squared_fit
   min_reduced_chi_squared_fit = minval(Reduced_Chi_Squared_array_fit)
   print *, "mininum reduced chi squared value for one MC iteration: ", min_reduced_chi_squared_fit
   loc_min_reduced_chi_squared_fit = minloc(Reduced_Chi_Squared_array_fit)
   print *, "mininum reduced chi squared location for one MC iteration: ", loc_min_reduced_chi_squared_fit

   IF (loc_min_chi_squared_fit(1) /= loc_min_reduced_chi_squared_fit(1)) THEN
      PRINT *, "The location of minimum value of chi squared and reduced chi squared not the same. Something went wrong!"
      PRINT *, "Minimum value location of chi squared: ", loc_min_chi_squared_fit(1)
      PRINT *, "Minimum value location of reduced chi squared: ", loc_min_reduced_chi_squared_fit(1)
   END IF

   min_chi_squared_mc(Monte_carlo_loop) = min_chi_squared_fit
   min_reduced_chi_squared_mc(Monte_carlo_loop) = min_reduced_chi_squared_fit

   !Best fit parameters for one MC iteration.
   best_RV_offset_datasets_fit = Parameter_array_fit(loc_min_chi_squared_fit(1),1)
   PRINT *, "The best RV offset between datasets corresponding to the lowest value of chi squared for one MC iteration: ", best_RV_offset_datasets_fit, "m/s"
   best_orbital_period_fit = (Parameter_array_fit(loc_min_chi_squared_fit(1),2)/10000.0D0) + Orbital_period_prior
   PRINT *, "The best orbital period corresponding to the lowest value of chi squared for one MC iteration: ", best_orbital_period_fit, "days"
   best_JD_time_mid_transit_fit = (Parameter_array_fit(loc_min_chi_squared_fit(1),3)/100.0D0) + JD_time_mid_transit_prior
   PRINT *, "The best mid transit time corresponding to the lowest value of chi squared for one MC iteration: ", best_JD_time_mid_transit_fit, "BJD"
   best_Mp_fit = Parameter_array_fit(loc_min_chi_squared_fit(1),4)
   PRINT *, "The best planetary mass corresponding to the lowest value of chi squared for one MC iteration: ", best_Mp_fit, "M_J"
   best_Rp_fit = Parameter_array_fit(loc_min_chi_squared_fit(1),5)
   PRINT *, "The best planetary radius corresponding to the lowest value of chi squared for one MC iteration: ", best_Rp_fit, "R_J"
   best_Ecc_fit = Parameter_array_fit(loc_min_chi_squared_fit(1),6)
   PRINT *, "The best eccentricity value corresponding to the lowest value of chi squared for one MC iteration: ", best_Ecc_fit
   best_Inc_fit = Parameter_array_fit(loc_min_chi_squared_fit(1),7)
   PRINT *, "The best inclination angle corresponding to the lowest value of chi squared for one MC iteration: ", best_Inc_fit, "deg"
   best_omega_arg_periastron_fit = Parameter_array_fit(loc_min_chi_squared_fit(1),8)
   PRINT *, "The best argument of periastron angle corresponding to the lowest value of chi squared for one MC iteration: ", best_omega_arg_periastron_fit, "deg"
   best_RV_zero_offset_fit = Parameter_array_fit(loc_min_chi_squared_fit(1),9)
   PRINT *, "The best RV zero offset velocity corresponding to the lowest value of chi squared for one MC iteration: ", best_RV_zero_offset_fit, "m/s"
   best_vsini_fit = Parameter_array_fit(loc_min_chi_squared_fit(1),10)
   PRINT *, "The best vsini value corresponding to the lowest value of chi squared for one MC iteration: ", best_vsini_fit, "m/s"
   best_spin_orbit_fit = Parameter_array_fit(loc_min_chi_squared_fit(1),11)
   PRINT *, "The best spin-orbit angle corresponding to the lowest value of chi squared for one MC iteration: ", best_spin_orbit_fit, "deg"
END IF

!STOP

!PRINT *," "
!PRINT *,"***********************************************************************************************"
!PRINT *," "

count_chi_fit = 0
DO i = 1, Number_iterations_fit
   IF (Chi_squared_array_fit(i) <= min_chi_squared_fit + chi_squared_change_fit) THEN
      count_chi_fit = count_chi_fit + 1
   END IF
END DO
Allocate(Parameter_chi_fit_array(count_chi_fit,12))

chi_counter = 0
DO i = 1, Number_iterations_fit
   IF (Chi_squared_array_fit(i) <= min_chi_squared_fit + chi_squared_change_fit) THEN
      chi_counter = chi_counter + 1
      Parameter_chi_fit_array(chi_counter,1) = Parameter_array_fit(i,1)
      Parameter_chi_fit_array(chi_counter,2) = Parameter_array_fit(i,2)
      Parameter_chi_fit_array(chi_counter,3) = Parameter_array_fit(i,3)
      Parameter_chi_fit_array(chi_counter,4) = Parameter_array_fit(i,4)
      Parameter_chi_fit_array(chi_counter,5) = Parameter_array_fit(i,5)
      Parameter_chi_fit_array(chi_counter,6) = Parameter_array_fit(i,6)
      Parameter_chi_fit_array(chi_counter,7) = Parameter_array_fit(i,7)
      Parameter_chi_fit_array(chi_counter,8) = Parameter_array_fit(i,8)
      Parameter_chi_fit_array(chi_counter,9) = Parameter_array_fit(i,9)
      Parameter_chi_fit_array(chi_counter,10) = Parameter_array_fit(i,10)
      Parameter_chi_fit_array(chi_counter,11) = Parameter_array_fit(i,11)
      Parameter_chi_fit_array(chi_counter,12) = Chi_squared_array_fit(i)
   END IF
END DO

!Positive uncertainty vsini
vsini_plus_error_fit = maxval(Parameter_chi_fit_array(:,10))
!PRINT *, "vsini plus error for one MC iteration: ", vsini_plus_error_fit - best_vsini_fit
loc_vsini_plus_error_fit = maxloc(Parameter_chi_fit_array(:,10))
!PRINT *, "vsini plus error location for one MC iteration: ", loc_vsini_plus_error_fit(1)
chi_square_vsini_plus_fit = Parameter_chi_fit_array(loc_vsini_plus_error_fit(1),12)
!PRINT *, "Chi squared value at vsini plus error for one MC iteration: ", chi_square_vsini_plus_fit

!Negative uncertainty vsini
vsini_minus_error_fit = minval(Parameter_chi_fit_array(:,10))
!PRINT *, "vsini minus error for one MC iteration: ", vsini_minus_error_fit - best_vsini_fit
loc_vsini_minus_error_fit = minloc(Parameter_chi_fit_array(:,10))
!PRINT *, "vsini minus error location for one MC iteration: ", loc_vsini_minus_error_fit(1)
chi_square_vsini_minus_fit = Parameter_chi_fit_array(loc_vsini_minus_error_fit(1),12)
!PRINT *, "Chi squared value at vsini minus error for one MC iteration: ", chi_square_vsini_minus_fit

mc_best_vsini_array(Monte_carlo_loop, 1) = best_vsini_fit
mc_best_vsini_array(Monte_carlo_loop, 2) = min_chi_squared_fit
mc_best_vsini_array(Monte_carlo_loop, 3) = min_reduced_chi_squared_fit
mc_best_vsini_array(Monte_carlo_loop, 4) = vsini_plus_error_fit - best_vsini_fit
mc_best_vsini_array(Monte_carlo_loop, 5) = vsini_minus_error_fit - best_vsini_fit
mc_best_vsini_array(Monte_carlo_loop, 6) = (abs(vsini_plus_error_fit - best_vsini_fit) + abs(vsini_minus_error_fit - best_vsini_fit))/2.0D0

!PRINT *," "
!PRINT *,"***********************************************************************************************"
!PRINT *," "

!Positive uncertainty spin_orbit
spin_orbit_plus_error_fit = maxval(Parameter_chi_fit_array(:,11))
!PRINT *, "spin_orbit plus error for one MC iteration: ", spin_orbit_plus_error_fit - best_spin_orbit_fit
loc_spin_orbit_plus_error_fit = maxloc(Parameter_chi_fit_array(:,11))
!PRINT *, "spin_orbit plus error location for one MC iteration: ", loc_spin_orbit_plus_error_fit(1)
chi_square_spin_orbit_plus_fit = Parameter_chi_fit_array(loc_spin_orbit_plus_error_fit(1),12)
!PRINT *, "Chi squared value at spin_orbit plus error for one MC iteration: ", chi_square_spin_orbit_plus_fit

!Negative uncertainty spin_orbit
spin_orbit_minus_error_fit = minval(Parameter_chi_fit_array(:,11))
!PRINT *, "spin_orbit minus error for one MC iteration: ", spin_orbit_minus_error_fit - best_spin_orbit_fit
loc_spin_orbit_minus_error_fit = minloc(Parameter_chi_fit_array(:,11))
!PRINT *, "spin_orbit minus error location for one MC iteration: ", loc_spin_orbit_minus_error_fit(1)
chi_square_spin_orbit_minus_fit = Parameter_chi_fit_array(loc_spin_orbit_minus_error_fit(1),12)
!PRINT *, "Chi squared value at spin_orbit minus error for one MC iteration: ", chi_square_spin_orbit_minus_fit

mc_best_lambda_array(Monte_carlo_loop, 1) = best_spin_orbit_fit
mc_best_lambda_array(Monte_carlo_loop, 2) = min_chi_squared_fit
mc_best_lambda_array(Monte_carlo_loop, 3) = min_reduced_chi_squared_fit
mc_best_lambda_array(Monte_carlo_loop, 4) = spin_orbit_plus_error_fit - best_spin_orbit_fit
mc_best_lambda_array(Monte_carlo_loop, 5) = spin_orbit_minus_error_fit - best_spin_orbit_fit
mc_best_lambda_array(Monte_carlo_loop, 6) = (abs(spin_orbit_plus_error_fit - best_spin_orbit_fit) + abs(spin_orbit_minus_error_fit - best_spin_orbit_fit))/2.0D0

PRINT *," "
PRINT *,"***********************************************************************************************"
PRINT *," "

!STOP

!Let's put these best fit values into the MC arrays

mc_best_RV_offset_datasets_array(Monte_carlo_loop, 1) = best_RV_offset_datasets_fit
mc_best_RV_offset_datasets_array(Monte_carlo_loop, 2) = min_chi_squared_fit
mc_best_RV_offset_datasets_array(Monte_carlo_loop, 3) = min_reduced_chi_squared_fit

mc_best_orbital_period_array(Monte_carlo_loop, 1) = best_orbital_period_fit
mc_best_orbital_period_array(Monte_carlo_loop, 2) = min_chi_squared_fit
mc_best_orbital_period_array(Monte_carlo_loop, 3) = min_reduced_chi_squared_fit

mc_best_JD_time_mid_transit_array(Monte_carlo_loop, 1) = best_JD_time_mid_transit_fit
mc_best_JD_time_mid_transit_array(Monte_carlo_loop, 2) = min_chi_squared_fit
mc_best_JD_time_mid_transit_array(Monte_carlo_loop, 3) = min_reduced_chi_squared_fit

mc_best_Mp_array(Monte_carlo_loop, 1) = best_Mp_fit
mc_best_Mp_array(Monte_carlo_loop, 2) = min_chi_squared_fit
mc_best_Mp_array(Monte_carlo_loop, 3) = min_reduced_chi_squared_fit

mc_best_Rp_array(Monte_carlo_loop, 1) = best_Rp_fit
mc_best_Rp_array(Monte_carlo_loop, 2) = min_chi_squared_fit
mc_best_Rp_array(Monte_carlo_loop, 3) = min_reduced_chi_squared_fit

mc_best_Ecc_array(Monte_carlo_loop, 1) = best_Ecc_fit
mc_best_Ecc_array(Monte_carlo_loop, 2) = min_chi_squared_fit
mc_best_Ecc_array(Monte_carlo_loop, 3) = min_reduced_chi_squared_fit

mc_best_Inc_array(Monte_carlo_loop, 1) = best_Inc_fit
mc_best_Inc_array(Monte_carlo_loop, 2) = min_chi_squared_fit
mc_best_Inc_array(Monte_carlo_loop, 3) = min_reduced_chi_squared_fit

mc_best_omega_arg_periastron_array(Monte_carlo_loop, 1) = best_omega_arg_periastron_fit
mc_best_omega_arg_periastron_array(Monte_carlo_loop, 2) = min_chi_squared_fit
mc_best_omega_arg_periastron_array(Monte_carlo_loop, 3) = min_reduced_chi_squared_fit

mc_best_RV_zero_offset_array(Monte_carlo_loop, 1) = best_RV_zero_offset_fit
mc_best_RV_zero_offset_array(Monte_carlo_loop, 2) = min_chi_squared_fit
mc_best_RV_zero_offset_array(Monte_carlo_loop, 3) = min_reduced_chi_squared_fit

!Deallocate all the arrays created in the MC loop.

DEALLOCATE(RV_offset_data_array, RV_synthetic_data_array)
DEALLOCATE(Data_my_rv_offset, data)
DEALLOCATE(Data_adjusted, sorted_data_array)
DEALLOCATE(Parameter_array_fit, Parameter_chi_fit_array)
DEALLOCATE(Chi_squared_array_fit)
DEALLOCATE(Reduced_Chi_Squared_array_fit)
DEALLOCATE(vsini_array, stellar_rotation_angle_array)
DEALLOCATE(time_data_array)
DEALLOCATE(fit_array_indeces)

!END DO Monte_carlo_loop
END DO 

!Let's use the Monte Carlo iterations to determine the overall best vsini and lambda including their uncertainties.
!First purge inactive memory to prevent swap memory being used.
!Purge memory? Works in OS X 10.8 or earlier. 10.9 and later require sudo purge and password.
!cmd1 = "purge"
!cmd1 = TRIM(cmd1)
!call execute_command_line(cmd1, wait=.true., exitstat=exitstat, cmdstat=cmdstat, cmdmsg=cmdmsg)

PRINT *," "
PRINT *,"***********************************************************************************************"
PRINT *," "

cmd1 = "say 'Finished Monte Carlo iterations.'"
cmd1 = TRIM(cmd1)
call execute_command_line(cmd1, wait=.false., exitstat=exitstat, cmdstat=cmdstat, cmdmsg=cmdmsg)

PRINT *, "number of total iterations: ", Number
PRINT *, "number of iterations there are suppose to be in total: ", Number_iterations_total

PRINT *," "
PRINT *,"***********************************************************************************************"
PRINT *," "

!Determine the lowest value for Chi square, reduced Chi square, and where these occur
!for the full MC parameter space.
min_chi_squared_total = minval(Chi_squared_array_total)
print *, "mininum chi squared value full MC parameter space: ", min_chi_squared_total
loc_min_chi_squared_total = minloc(Chi_squared_array_total)
print *, "mininum chi squared location full MC parameter space: ", loc_min_chi_squared_total
min_reduced_chi_squared_total = minval(Reduced_Chi_Squared_array_total)
print *, "mininum reduced chi squared value of full MC parameter space: ", min_reduced_chi_squared_total
loc_min_reduced_chi_squared_total = minloc(Reduced_Chi_Squared_array_total)
print *, "mininum reduced chi squared location full MC parameter space: ", loc_min_reduced_chi_squared_total

IF (loc_min_chi_squared_total(1) /= loc_min_reduced_chi_squared_total(1)) THEN
   PRINT *, "The location of minimum value of chi squared and reduced chi squared not the same. Something went wrong!"
   PRINT *, "Minimum value location of chi squared: ", loc_min_chi_squared_total(1)
   PRINT *, "Minimum value location of reduced chi squared: ", loc_min_reduced_chi_squared_total(1)
END IF

PRINT *," "
PRINT *,"***********************************************************************************************"
PRINT *," "

!Now determine the best parameters corresponding to the minimum chi square value in the MC parameter space.
best_RV_offset_datasets_total = Parameter_array_total(loc_min_chi_squared_total(1),1)
PRINT *, "The best RV offset between datasets corresponding to the lowest value of chi squared in full MC parameter space: ", best_RV_offset_datasets_total, "m/s"
best_orbital_period_total = (Parameter_array_total(loc_min_chi_squared_total(1),2)/10000.0D0) + Orbital_period_prior
PRINT *, "The best orbital period corresponding to the lowest value of chi squared in full MC parameter space: ", best_orbital_period_total, "days"
best_JD_time_mid_transit_total = (Parameter_array_total(loc_min_chi_squared_total(1),3)/100.0D0) + JD_time_mid_transit_prior
PRINT *, "The best mid transit time corresponding to the lowest value of chi squared in full MC parameter space: ", best_JD_time_mid_transit_total, "BJD"
best_Mp_total = Parameter_array_total(loc_min_chi_squared_total(1),4)
PRINT *, "The best planetary mass corresponding to the lowest value of chi squared in full MC parameter space: ", best_Mp_total, "M_J"
best_Rp_total = Parameter_array_total(loc_min_chi_squared_total(1),5)
PRINT *, "The best planetary radius corresponding to the lowest value of chi squared in full MC parameter space: ", best_Rp_total, "R_J"
best_Ecc_total = Parameter_array_total(loc_min_chi_squared_total(1),6)
PRINT *, "The best eccentricity value corresponding to the lowest value of chi squared in full MC parameter space: ", best_Ecc_total
best_Inc_total = Parameter_array_total(loc_min_chi_squared_total(1),7)
PRINT *, "The best inclination angle corresponding to the lowest value of chi squared in full MC parameter space: ", best_Inc_total, "deg"
best_omega_arg_periastron_total = Parameter_array_total(loc_min_chi_squared_total(1),8)
PRINT *, "The best argument of periastron angle corresponding to the lowest value of chi squared in full MC parameter space: ", best_omega_arg_periastron_total, "deg"
best_RV_zero_offset_total = Parameter_array_total(loc_min_chi_squared_total(1),9)
PRINT *, "The best RV zero offset velocity corresponding to the lowest value of chi squared in full MC parameter space: ", best_RV_zero_offset_total, "m/s"
best_vsini_total = Parameter_array_total(loc_min_chi_squared_total(1),10)
PRINT *, "The best vsini value corresponding to the lowest value of chi squared in full MC parameter space: ", best_vsini_total, "m/s"
best_spin_orbit_total = Parameter_array_total(loc_min_chi_squared_total(1),11)
PRINT *, "The best spin-orbit angle corresponding to the lowest value of chi squared in full MC parameter space: ", best_spin_orbit_total, "deg"

PRINT *, ''
call cpu_time(cpu_time_it3)
time_remaining = (cpu_time_it3 - cpu_time_it1)/1.5D0
WRITE (*,1002) 'Estimated time remaining: ', time_remaining/60.0D0, ' minutes'
1002 FORMAT(A30,F8.2,A8)

PRINT *," "
PRINT *,"***********************************************************************************************"
PRINT *," "

!STOP

!Probably the most robust way of getting the 1-sigma uncertainty on vsini and lambda is to compute standard deviation from the variance.
DO l = 1, mc_sample_size
   variance_vsini = variance_vsini + (mc_best_vsini_array(l, 1) - best_vsini_total)**2.0D0 + (mc_best_vsini_array(l, 6))**2.0D0
   variance_vsini_chi_squared = variance_vsini_chi_squared + (mc_best_vsini_array(l, 2) - min_chi_squared_total)**2.0D0
   !variance_vsini_chi_squared = variance_vsini_chi_squared + (chi_squared_change_fit)**2.0D0
   
   variance_spin_orbit = variance_spin_orbit + (mc_best_lambda_array(l, 1) - best_spin_orbit_total)**2.0D0 + (mc_best_lambda_array(l, 6))**2.0D0
   variance_spin_orbit_chi_squared = variance_spin_orbit_chi_squared + (mc_best_lambda_array(l, 2) - min_chi_squared_total)**2.0D0
   !variance_spin_orbit_chi_squared = variance_spin_orbit_chi_squared + (chi_squared_change_fit)**2.0D0
END DO

PRINT *, 'variance_vsini: ', variance_vsini
PRINT *, 'variance_vsini_chi_squared: ', variance_vsini_chi_squared
PRINT *, 'variance_spin_orbit: ', variance_spin_orbit
PRINT *, 'variance_spin_orbit_chi_squared: ', variance_spin_orbit_chi_squared

vsini_error_mc = (SQRT((1.0D0/(FLOAT(mc_sample_size) - 1.0D0))*variance_vsini))
vsini_delta_chi = (SQRT((1.0D0/(FLOAT(mc_sample_size) - 1.0D0))*variance_vsini_chi_squared))

spin_orbit_error_mc = (SQRT((1.0D0/(FLOAT(mc_sample_size) - 1.0D0))*variance_spin_orbit))
spin_orbit_delta_chi = (SQRT((1.0D0/(FLOAT(mc_sample_size) - 1.0D0))*variance_spin_orbit_chi_squared))

average_delta_chi = (vsini_delta_chi + vsini_delta_chi) / 2.0D0
PRINT *, 'average_delta_chi: ', average_delta_chi

!Write the delta chi square for IDL.
OPEN(unit=99, FILE=output_delta_chi_square_mc_filename, status='replace', action='write')
WRITE (99,64) vsini_delta_chi
WRITE (99,65) spin_orbit_delta_chi
WRITE (99,66) average_delta_chi
64 FORMAT(F10.4) 
65 FORMAT(F10.4) 
66 FORMAT(F10.4) 
CLOSE(99)

!Positive uncertainty vsini
vsini_plus_error_mc1 = best_vsini_total + vsini_error_mc
!Negative uncertainty vsini
vsini_minus_error_mc1 = best_vsini_total - vsini_error_mc
PRINT *, "vsini error for MC parameter space from standard deviation: ", vsini_error_mc

!Positive uncertainty lambda
spin_orbit_plus_error_mc1 = best_spin_orbit_total + spin_orbit_error_mc
!Negative uncertainty lambda
spin_orbit_minus_error_mc1 = best_spin_orbit_total - spin_orbit_error_mc
PRINT *, "spin_orbit error for MC parameter space from standard deviation: ", spin_orbit_error_mc

!STOP

!This part determines the 1 sigma uncertainty in vsini and lambda. This is done by first reordering the best fit chi square array
!in order of lowest minimum chi square along with the corresponding parameter array. Then we determine the chi square cutoff such
!that ~68.3% of the minimum chi square values and associated parameters lie within the cutoff.
CALL Selection_sort(min_chi_squared_mc,index_array)

Allocate(sorted_min_chi_squared_mc(mc_sample_size))
Allocate(sorted_min_reduced_chi_squared_mc(mc_sample_size))
Allocate(sorted_mc_best_vsini_array(mc_sample_size,6))
Allocate(sorted_mc_best_lambda_array(mc_sample_size,6))
Allocate(sorted_mc_best_RV_offset_datasets_array(mc_sample_size,3))
Allocate(sorted_mc_best_orbital_period_array(mc_sample_size,3))
Allocate(sorted_mc_best_JD_time_mid_transit_array(mc_sample_size,3))
Allocate(sorted_mc_best_Mp_array(mc_sample_size,3))
Allocate(sorted_mc_best_Rp_array(mc_sample_size,3))
Allocate(sorted_mc_best_Ecc_array(mc_sample_size,3))
Allocate(sorted_mc_best_Inc_array(mc_sample_size,3))
Allocate(sorted_mc_best_omega_arg_periastron_array(mc_sample_size,3))
Allocate(sorted_mc_best_RV_zero_offset_array(mc_sample_size,3))

DO l = 1, mc_sample_size
   sorted_min_chi_squared_mc(l) = min_chi_squared_mc(index_array(l))
   
   sorted_min_reduced_chi_squared_mc(l) = min_reduced_chi_squared_mc(index_array(l))
   
   sorted_mc_best_vsini_array(l, 1) = mc_best_vsini_array(index_array(l), 1)
   sorted_mc_best_vsini_array(l, 2) = mc_best_vsini_array(index_array(l), 2)
   sorted_mc_best_vsini_array(l, 3) = mc_best_vsini_array(index_array(l), 3)
   sorted_mc_best_vsini_array(l, 4) = mc_best_vsini_array(index_array(l), 4)
   sorted_mc_best_vsini_array(l, 5) = mc_best_vsini_array(index_array(l), 5)
   sorted_mc_best_vsini_array(l, 6) = mc_best_vsini_array(index_array(l), 6)
   
   sorted_mc_best_lambda_array(l, 1) = mc_best_lambda_array(index_array(l), 1)
   sorted_mc_best_lambda_array(l, 2) = mc_best_lambda_array(index_array(l), 2)
   sorted_mc_best_lambda_array(l, 3) = mc_best_lambda_array(index_array(l), 3)
   sorted_mc_best_lambda_array(l, 4) = mc_best_lambda_array(index_array(l), 4)
   sorted_mc_best_lambda_array(l, 5) = mc_best_lambda_array(index_array(l), 5)
   sorted_mc_best_lambda_array(l, 6) = mc_best_lambda_array(index_array(l), 6)
   
   sorted_mc_best_RV_offset_datasets_array(l, 1) = mc_best_RV_offset_datasets_array(index_array(l), 1)
   sorted_mc_best_RV_offset_datasets_array(l, 2) = mc_best_RV_offset_datasets_array(index_array(l), 2)
   sorted_mc_best_RV_offset_datasets_array(l, 3) = mc_best_RV_offset_datasets_array(index_array(l), 3)
   
   sorted_mc_best_orbital_period_array(l, 1) = mc_best_orbital_period_array(index_array(l), 1)
   sorted_mc_best_orbital_period_array(l, 2) = mc_best_orbital_period_array(index_array(l), 2)
   sorted_mc_best_orbital_period_array(l, 3) = mc_best_orbital_period_array(index_array(l), 3)
   
   sorted_mc_best_JD_time_mid_transit_array(l, 1) = mc_best_JD_time_mid_transit_array(index_array(l), 1)
   sorted_mc_best_JD_time_mid_transit_array(l, 2) = mc_best_JD_time_mid_transit_array(index_array(l), 2)
   sorted_mc_best_JD_time_mid_transit_array(l, 3) = mc_best_JD_time_mid_transit_array(index_array(l), 3)
   
   sorted_mc_best_Mp_array(l, 1) = mc_best_Mp_array(index_array(l), 1)
   sorted_mc_best_Mp_array(l, 2) = mc_best_Mp_array(index_array(l), 2)
   sorted_mc_best_Mp_array(l, 3) = mc_best_Mp_array(index_array(l), 3)
   
   sorted_mc_best_Rp_array(l, 1) = mc_best_Rp_array(index_array(l), 1)
   sorted_mc_best_Rp_array(l, 2) = mc_best_Rp_array(index_array(l), 2)
   sorted_mc_best_Rp_array(l, 3) = mc_best_Rp_array(index_array(l), 3)
   
   sorted_mc_best_Ecc_array(l, 1) = mc_best_Ecc_array(index_array(l), 1)
   sorted_mc_best_Ecc_array(l, 2) = mc_best_Ecc_array(index_array(l), 2)
   sorted_mc_best_Ecc_array(l, 3) = mc_best_Ecc_array(index_array(l), 3)
   
   sorted_mc_best_Inc_array(l, 1) = mc_best_Inc_array(index_array(l), 1)
   sorted_mc_best_Inc_array(l, 2) = mc_best_Inc_array(index_array(l), 2)
   sorted_mc_best_Inc_array(l, 3) = mc_best_Inc_array(index_array(l), 3)
   
   sorted_mc_best_omega_arg_periastron_array(l, 1) = mc_best_omega_arg_periastron_array(index_array(l), 1)
   sorted_mc_best_omega_arg_periastron_array(l, 2) = mc_best_omega_arg_periastron_array(index_array(l), 2)
   sorted_mc_best_omega_arg_periastron_array(l, 3) = mc_best_omega_arg_periastron_array(index_array(l), 3)
   
   sorted_mc_best_RV_zero_offset_array(l, 1) = mc_best_RV_zero_offset_array(index_array(l), 1)
   sorted_mc_best_RV_zero_offset_array(l, 2) = mc_best_RV_zero_offset_array(index_array(l), 2)
   sorted_mc_best_RV_zero_offset_array(l, 3) = mc_best_RV_zero_offset_array(index_array(l), 3)
   IF (index_array(l) /= l) THEN
      !PRINT *, "index_array(l): ", index_array(l)
      !PRINT *, "l: ", l
      array_resorted = 'Y'
   ELSE
      array_resorted = 'N'
   END IF
END DO

!to get the delta chi square, iterate over the sorted array until either the vsini or lambda value is outside the error range.
DO l = 1, mc_sample_size
   delta_chi_square = sorted_mc_best_vsini_array(l, 2) - min_chi_squared_total
   IF ((sorted_mc_best_vsini_array(l, 1) > best_vsini_total + vsini_error_mc) .OR. (sorted_mc_best_vsini_array(l, 1) < best_vsini_total - vsini_error_mc) &
      .OR. (sorted_mc_best_lambda_array(l, 1) > best_spin_orbit_total + spin_orbit_error_mc) &
      .OR. (sorted_mc_best_lambda_array(l, 1) < best_spin_orbit_total - spin_orbit_error_mc)) EXIT
END DO

PRINT *, 'delta_chi_square: ', delta_chi_square

!Write the delta chi square for IDL.
OPEN(unit=99, FILE=output_delta_chi_square_mc_filename2, status='replace', action='write')
WRITE (99,67) delta_chi_square 
67 FORMAT(F10.4) 
CLOSE(99)

cutoff_position = CEILING(mc_sample_size*0.683D0)

Allocate(cut_min_chi_squared_mc(cutoff_position))
Allocate(cut_min_reduced_chi_squared_mc(cutoff_position))
Allocate(cut_mc_best_vsini_array(cutoff_position,6))
Allocate(cut_mc_best_lambda_array(cutoff_position,6))
Allocate(cut_mc_best_RV_offset_datasets_array(cutoff_position,3))
Allocate(cut_mc_best_orbital_period_array(cutoff_position,3))
Allocate(cut_mc_best_JD_time_mid_transit_array(cutoff_position,3))
Allocate(cut_mc_best_Mp_array(cutoff_position,3))
Allocate(cut_mc_best_Rp_array(cutoff_position,3))
Allocate(cut_mc_best_Ecc_array(cutoff_position,3))
Allocate(cut_mc_best_Inc_array(cutoff_position,3))
Allocate(cut_mc_best_omega_arg_periastron_array(cutoff_position,3))
Allocate(cut_mc_best_RV_zero_offset_array(cutoff_position,3))

DO l = 1, cutoff_position
   cut_min_chi_squared_mc(l) = sorted_min_chi_squared_mc(l)
   
   cut_min_reduced_chi_squared_mc(l) = sorted_min_reduced_chi_squared_mc(l)
   
   cut_mc_best_vsini_array(l, 1) = sorted_mc_best_vsini_array(l, 1)
   cut_mc_best_vsini_array(l, 2) = sorted_mc_best_vsini_array(l, 2)
   cut_mc_best_vsini_array(l, 3) = sorted_mc_best_vsini_array(l, 3)
   cut_mc_best_vsini_array(l, 4) = sorted_mc_best_vsini_array(l, 4)
   cut_mc_best_vsini_array(l, 5) = sorted_mc_best_vsini_array(l, 5)
   cut_mc_best_vsini_array(l, 6) = sorted_mc_best_vsini_array(l, 6)
   
   cut_mc_best_lambda_array(l, 1) = sorted_mc_best_lambda_array(l, 1)
   cut_mc_best_lambda_array(l, 2) = sorted_mc_best_lambda_array(l, 2)
   cut_mc_best_lambda_array(l, 3) = sorted_mc_best_lambda_array(l, 3)
   cut_mc_best_lambda_array(l, 4) = sorted_mc_best_lambda_array(l, 4)
   cut_mc_best_lambda_array(l, 5) = sorted_mc_best_lambda_array(l, 5)
   cut_mc_best_lambda_array(l, 6) = sorted_mc_best_lambda_array(l, 6)
   
   cut_mc_best_RV_offset_datasets_array(l, 1) = sorted_mc_best_RV_offset_datasets_array(l, 1)
   cut_mc_best_RV_offset_datasets_array(l, 2) = sorted_mc_best_RV_offset_datasets_array(l, 2)
   cut_mc_best_RV_offset_datasets_array(l, 3) = sorted_mc_best_RV_offset_datasets_array(l, 3)
   
   cut_mc_best_orbital_period_array(l, 1) = sorted_mc_best_orbital_period_array(l, 1)
   cut_mc_best_orbital_period_array(l, 2) = sorted_mc_best_orbital_period_array(l, 2)
   cut_mc_best_orbital_period_array(l, 3) = sorted_mc_best_orbital_period_array(l, 3)
   
   cut_mc_best_JD_time_mid_transit_array(l, 1) = sorted_mc_best_JD_time_mid_transit_array(l, 1)
   cut_mc_best_JD_time_mid_transit_array(l, 2) = sorted_mc_best_JD_time_mid_transit_array(l, 2)
   cut_mc_best_JD_time_mid_transit_array(l, 3) = sorted_mc_best_JD_time_mid_transit_array(l, 3)
   
   cut_mc_best_Mp_array(l, 1) = sorted_mc_best_Mp_array(l, 1)
   cut_mc_best_Mp_array(l, 2) = sorted_mc_best_Mp_array(l, 2)
   cut_mc_best_Mp_array(l, 3) = sorted_mc_best_Mp_array(l, 3)
   
   cut_mc_best_Rp_array(l, 1) = sorted_mc_best_Rp_array(l, 1)
   cut_mc_best_Rp_array(l, 2) = sorted_mc_best_Rp_array(l, 2)
   cut_mc_best_Rp_array(l, 3) = sorted_mc_best_Rp_array(l, 3)
   
   cut_mc_best_Ecc_array(l, 1) = sorted_mc_best_Ecc_array(l, 1)
   cut_mc_best_Ecc_array(l, 2) = sorted_mc_best_Ecc_array(l, 2)
   cut_mc_best_Ecc_array(l, 3) = sorted_mc_best_Ecc_array(l, 3)
   
   cut_mc_best_Inc_array(l, 1) = sorted_mc_best_Inc_array(l, 1)
   cut_mc_best_Inc_array(l, 2) = sorted_mc_best_Inc_array(l, 2)
   cut_mc_best_Inc_array(l, 3) = sorted_mc_best_Inc_array(l, 3)
   
   cut_mc_best_omega_arg_periastron_array(l, 1) = sorted_mc_best_omega_arg_periastron_array(l, 1)
   cut_mc_best_omega_arg_periastron_array(l, 2) = sorted_mc_best_omega_arg_periastron_array(l, 2)
   cut_mc_best_omega_arg_periastron_array(l, 3) = sorted_mc_best_omega_arg_periastron_array(l, 3)
   
   cut_mc_best_RV_zero_offset_array(l, 1) = sorted_mc_best_RV_zero_offset_array(l, 1)
   cut_mc_best_RV_zero_offset_array(l, 2) = sorted_mc_best_RV_zero_offset_array(l, 2)
   cut_mc_best_RV_zero_offset_array(l, 3) = sorted_mc_best_RV_zero_offset_array(l, 3)
END DO

chi_squared_change_mc = maxval(cut_min_chi_squared_mc) - minval(cut_min_chi_squared_mc)

!Find the 1 sigma uncertainties in these parameters.
!Positive uncertainty vsini
vsini_plus_error_mc = maxval(cut_mc_best_vsini_array(:,1))
PRINT *, "vsini plus error MC parameter space: ", vsini_plus_error_mc - best_vsini_total
loc_vsini_plus_error_mc = maxloc(cut_mc_best_vsini_array(:,1))
PRINT *, "vsini plus error location MC parameter space: ", loc_vsini_plus_error_mc(1)
chi_square_vsini_plus_mc = cut_min_chi_squared_mc(loc_vsini_plus_error_mc(1))
PRINT *, "Chi squared value at vsini plus error MC parameter space: ", chi_square_vsini_plus_mc

!Negative uncertainty vsini
vsini_minus_error_mc = minval(cut_mc_best_vsini_array(:,1))
PRINT *, "vsini minus error MC parameter space: ", vsini_minus_error_mc - best_vsini_total
loc_vsini_minus_error_mc = minloc(cut_mc_best_vsini_array(:,1))
PRINT *, "vsini minus error location MC parameter space: ", loc_vsini_minus_error_mc(1)
chi_square_vsini_minus_mc = cut_min_chi_squared_mc(loc_vsini_minus_error_mc(1))
PRINT *, "Chi squared value at vsini minus error MC parameter space: ", chi_square_vsini_minus_mc

PRINT *," "
PRINT *,"***********************************************************************************************"
PRINT *," "

!Positive uncertainty spin_orbit
spin_orbit_plus_error_mc = maxval(cut_mc_best_lambda_array(:,1))
PRINT *, "spin_orbit plus error MC parameter space: ", spin_orbit_plus_error_mc - best_spin_orbit_total
loc_spin_orbit_plus_error_mc = maxloc(cut_mc_best_lambda_array(:,1))
PRINT *, "spin_orbit plus error location MC parameter space: ", loc_spin_orbit_plus_error_mc(1)
chi_square_spin_orbit_plus_mc = cut_min_chi_squared_mc(loc_spin_orbit_plus_error_mc(1))
PRINT *, "Chi squared value at spin_orbit plus error MC parameter space: ", chi_square_spin_orbit_plus_mc

!Negative uncertainty spin_orbit
spin_orbit_minus_error_mc = minval(cut_mc_best_lambda_array(:,1))
PRINT *, "spin_orbit minus error MC parameter space: ", spin_orbit_minus_error_mc - best_spin_orbit_total
loc_spin_orbit_minus_error_mc = minloc(cut_mc_best_lambda_array(:,1))
PRINT *, "spin_orbit minus error location MC parameter space: ", loc_spin_orbit_minus_error_mc(1)
chi_square_spin_orbit_minus_mc = cut_min_chi_squared_mc(loc_spin_orbit_minus_error_mc(1))
PRINT *, "Chi squared value at spin_orbit minus error MC parameter space: ", chi_square_spin_orbit_minus_mc

PRINT *," "
PRINT *,"***********************************************************************************************"
PRINT *," "

!STOP

!*****************************************************************************************************************************************

!Write output chi squared array so that IDL can read the array in.
OPEN(unit=99, FILE=output_chi_squared_array_total_filename_mc, status='replace', action='write')
!Now write array line by line into chi squared array file.
DO l = 1, Number_iterations_total
   WRITE (99,53) Chi_squared_array_total(l)
   53 FORMAT(F30.15)
END DO 
CLOSE(99)

!Write output chi squared array so that IDL can read the array in.
OPEN(unit=99, FILE=output_chi_squared_array_1sigcut_filename_mc, status='replace', action='write')
!Now write array line by line into chi squared array file.
DO l = 1, cutoff_position
   WRITE (99,54) cut_min_chi_squared_mc(l)
   54 FORMAT(F30.15)
END DO 
CLOSE(99)

!Write output reduced chi squared array so that IDL can read the array in.
OPEN(unit=99, FILE=output_reduced_chi_array_total_filename_mc, status='replace', action='write')
!Now write array line by line into reduced chi squared array file.
DO l = 1, Number_iterations_total
   WRITE (99,55) Reduced_Chi_Squared_array_total(l)
   55 FORMAT(F20.10)
END DO 
CLOSE(99)

!Write output reduced chi squared array so that IDL can read the array in.
OPEN(unit=99, FILE=output_reduced_chi_array_1sigcut_filename_mc, status='replace', action='write')
!Now write array line by line into reduced chi squared array file.
DO l = 1, cutoff_position
   WRITE (99,56) cut_min_reduced_chi_squared_mc(l)
   56 FORMAT(F20.10)
END DO 
CLOSE(99)

!Write output Mp array so that IDL can read the array in.
OPEN(unit=99, FILE=output_best_vsini_mc_array_filename, status='replace', action='write')
!Now write best vsini array line by line into a file.
DO l = 1, mc_sample_size
   WRITE (99,57) mc_best_vsini_array(l, 1), mc_best_vsini_array(l, 2), mc_best_vsini_array(l, 3), mc_best_vsini_array(l, 4), &
   mc_best_vsini_array(l, 5), mc_best_vsini_array(l, 6)
   57 FORMAT(F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6)
END DO
CLOSE(99)

!Write output Mp array so that IDL can read the array in.
OPEN(unit=99, FILE=output_best_spin_orbit_mc_array_filename, status='replace', action='write')
!Now write best vsini array line by line into a file.
DO l = 1, mc_sample_size
   WRITE (99,58) mc_best_lambda_array(l, 1), mc_best_lambda_array(l, 2), mc_best_lambda_array(l, 3), mc_best_lambda_array(l, 4), &
   mc_best_lambda_array(l, 5), mc_best_lambda_array(l, 6)
   58 FORMAT(F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6)
END DO
CLOSE(99)

!Write output Mp array so that IDL can read the array in.
OPEN(unit=99, FILE=output_1sigma_vsini_mc_array_filename, status='replace', action='write')
!Now write best vsini array line by line into a file.
DO l = 1, cutoff_position
   WRITE (99,59) cut_mc_best_vsini_array(l, 1), cut_mc_best_vsini_array(l, 2), cut_mc_best_vsini_array(l, 3), cut_mc_best_vsini_array(l, 4), &
   cut_mc_best_vsini_array(l, 5), cut_mc_best_vsini_array(l, 6)
   59 FORMAT(F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6)
END DO
CLOSE(99)
   
!Write output Mp array so that IDL can read the array in.
OPEN(unit=99, FILE=output_1sigma_spin_orbit_mc_array_filename, status='replace', action='write')
!Now write best vsini array line by line into a file.
DO l = 1, cutoff_position
   WRITE (99,60) cut_mc_best_lambda_array(l, 1), cut_mc_best_lambda_array(l, 2), cut_mc_best_lambda_array(l, 3), &
   cut_mc_best_lambda_array(l, 4), cut_mc_best_lambda_array(l, 5), cut_mc_best_lambda_array(l, 6)
   60 FORMAT(F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6)
END DO
CLOSE(99)

!Write parameter output array so that IDL can read the array in.
OPEN(unit=99, FILE=output_parameter_array_total_mc_filename, status='replace', action='write')
!Now write array line by line into parameter array file.
DO l = 1, Number_iterations_total
   WRITE (99,61) Parameter_array_total(l,1), (Parameter_array_total(l,2)/10000.0D0) + Orbital_period_prior, &
   (Parameter_array_total(l,3)/100.0D0) + JD_time_mid_transit_prior, Parameter_array_total(l,4), Parameter_array_total(l,5), Parameter_array_total(l,6), &
   Parameter_array_total(l,7), Parameter_array_total(l,8), Parameter_array_total(l,9), Parameter_array_total(l,10), Parameter_array_total(l,11)
   61 FORMAT(F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6, 1X, F15.6, &
   1X, F15.6, 1X, F15.6, 1X, F15.6)
END DO 
CLOSE(99)

!Write additional output parameters to IDL.
OPEN(unit=99, FILE=output_parameter_total_mc_idl_filename, status='replace', action='write')
!Now write out the parameters for IDL to read in.
WRITE (99,62) Number
62 FORMAT(I10)
CLOSE(99)

!Write additional output parameters to IDL.
OPEN(unit=99, FILE=output_parameter_cutoff_mc_idl_filename, status='replace', action='write')
!Now write out the parameters for IDL to read in.
WRITE (99,63) cutoff_position
63  FORMAT(I10)
CLOSE(99)

!Write out best vsini, best spin_orbit, minimum chi squared, minimum reduced chi squared, 
!location of minimum chi squared, uncertainty for vsini and spin_orbit, and locations for 
!these uncertainties.
OPEN(unit=99, FILE=output_best_parameters_mc_filename, status='replace', action='write')
WRITE(99,133) 'Minimum Chi Squared:', min_chi_squared_total
WRITE(99,134) 'Location of Minimum Chi Squared:', loc_min_chi_squared_total(1)
WRITE(99,135) 'Minimum Reduced Chi Squared:', min_reduced_chi_squared_total
WRITE(99,136) 'Location of Minimum Reduced Chi Squared:', loc_min_reduced_chi_squared_total(1)
WRITE(99,137) 'Best vsin(i) value:', best_vsini_total, 'm/s'
WRITE(99,138) 'Best vsin(i) value priors fixed:', best_vsini_fit_prior, 'm/s'
WRITE(99,139) 'Positive uncertainty on vsin(i) value:', vsini_error_mc, 'm/s'
WRITE(99,140) 'Negative uncertainty on vsin(i) value:', -vsini_error_mc, 'm/s'
WRITE(99,141) 'Best spin-orbit alignment angle:', best_spin_orbit_total, 'deg'
WRITE(99,142) 'Best spin-orbit alignment angle priors fixed:', best_spin_orbit_fit_prior, 'deg'
WRITE(99,143) 'Positive uncertainty on spin-orbit angle:', spin_orbit_error_mc, 'deg'
WRITE(99,144) 'Negative uncertainty on spin-orbit angle:', -spin_orbit_error_mc, 'deg'
WRITE(99,145) 'Best orbital period:', best_orbital_period_total, 'days'
WRITE(99,146) 'Best JD time mid transit:', best_JD_time_mid_transit_total, 'days'
WRITE(99,147) 'Best Mp:', best_Mp_total, 'Mj'
WRITE(99,148) 'Best Rp:', best_Rp_total, 'Rj'
WRITE(99,149) 'Best Ecc:', best_Ecc_total
WRITE(99,150) 'Best Inc:', best_Inc_total, 'deg'
WRITE(99,151) 'Best omega arg periastron:', best_omega_arg_periastron_total, 'deg'
WRITE(99,152) 'Best RV zero offset:', best_RV_zero_offset_total, 'm/s'
WRITE(99,153) 'Best RV offset between datasets:', best_RV_offset_datasets_total, 'm/s'
133  FORMAT(A50, 2X, F15.6)
134  FORMAT(A50, 2X, I10)
135  FORMAT(A50, 2X, F15.6)
136  FORMAT(A50, 2X, I10)
137  FORMAT(A50, 2X, F15.6, 5X, A3)
138  FORMAT(A50, 2X, F15.6, 5X, A3)
139  FORMAT(A50, 2X, F15.6, 5X, A3)
140  FORMAT(A50, 2X, F15.6, 5X, A3)
141  FORMAT(A50, 2X, F15.6, 5X, A3)
142  FORMAT(A50, 2X, F15.6, 5X, A3)
143  FORMAT(A50, 2X, F15.6, 5X, A3)
144  FORMAT(A50, 2X, F15.6, 5X, A3)
145  FORMAT(A50, 2X, F15.6, 5X, A4)
146  FORMAT(A50, 2X, F15.6, 5X, A4)
147  FORMAT(A50, 2X, F15.6, 5X, A2)
148  FORMAT(A50, 2X, F15.6, 5X, A2)
149  FORMAT(A50, 2X, F15.6)
150  FORMAT(A50, 2X, F15.6, 5X, A3)
151  FORMAT(A50, 2X, F15.6, 5X, A3)
152  FORMAT(A50, 2X, F15.6, 5X, A3)
153  FORMAT(A50, 2X, F15.6, 5X, A3)
CLOSE(99)

OPEN(unit=99, FILE=output_mc_uncert_vsini_lambda_filename, status='replace', action='write')
WRITE(99,154) 'vsini uncertainty from mc variance:', vsini_error_mc, 'm/s'
WRITE(99,155) 'spin orbit uncertainty from mc variance:', spin_orbit_error_mc, 'deg'
154  FORMAT(A50, 2X, F15.6, 5X, A3)
155  FORMAT(A50, 2X, F15.6, 5X, A3)
CLOSE(99)

OPEN(unit=99, FILE=output_mc_parm_uncert_vsini_lambda_filename, status='replace', action='write')
WRITE(99,156) 'Positive uncertainty on vsin(i) value:', vsini_plus_error_mc - best_vsini_total, 'm/s'
WRITE(99,157) 'Negative uncertainty on vsin(i) value:', vsini_minus_error_mc - best_vsini_total, 'm/s'
WRITE(99,158) 'Positive uncertainty on spin-orbit angle:', spin_orbit_plus_error_mc - best_spin_orbit_total, 'deg'
WRITE(99,159) 'Negative uncertainty on spin-orbit angle:', spin_orbit_minus_error_mc - best_spin_orbit_total, 'deg'
156  FORMAT(A50, 2X, F15.6, 5X, A3)
157  FORMAT(A50, 2X, F15.6, 5X, A3)
158  FORMAT(A50, 2X, F15.6, 5X, A3)
159  FORMAT(A50, 2X, F15.6, 5X, A3)
CLOSE(99)

OPEN(unit=99, FILE=output_mc_cut_array_indeces_filename, status='replace', action='write')
!Now write array line by line into stellar_rotation_angle array file.
DO l = 1, cutoff_position
   WRITE (99,229) index_array(l)
   229 FORMAT(I7)
END DO 
CLOSE(99)

PRINT *, 'Number of iterations ',  Number_iterations_total
PRINT *, 'Number_vsini_mc ',  Number_vsini_mc
PRINT *, 'Number_stellar_rotation_angle_mc ',  Number_stellar_rotation_angle_mc

call cpu_time(cpu_time_4)
PRINT *, 'The time taken to run is: ', cpu_time_4 - cpu_time_1, ' seconds'

DEALLOCATE(mc_parameter_array)
DEALLOCATE(Data_my_rv)
DEALLOCATE(Data_other_rv)
DEALLOCATE(min_chi_squared_mc)
DEALLOCATE(min_reduced_chi_squared_mc)
DEALLOCATE(mc_best_vsini_array)
DEALLOCATE(mc_best_lambda_array)
DEALLOCATE(mc_best_RV_offset_datasets_array)
DEALLOCATE(mc_best_orbital_period_array)
DEALLOCATE(mc_best_JD_time_mid_transit_array)
DEALLOCATE(mc_best_Mp_array)
DEALLOCATE(mc_best_Rp_array)
DEALLOCATE(mc_best_Ecc_array)
DEALLOCATE(mc_best_Inc_array)
DEALLOCATE(mc_best_omega_arg_periastron_array)
DEALLOCATE(mc_best_RV_zero_offset_array)

!Will determine the size of these other arrays soon. Placed them here for the time being.
DEALLOCATE(Number_array_total)
!Allocate(check_array_nettotal(Number_iterations_total))

!The size of the following array are based on the mc_sample_size.
DEALLOCATE(Inc_array)
DEALLOCATE(omega_arg_periastron_array)
DEALLOCATE(Ecc_array)
DEALLOCATE(Rp_array)
DEALLOCATE(Mp_array)
DEALLOCATE(Orbital_period_array)
DEALLOCATE(JD_time_mid_transit_array)
DEALLOCATE(RV_zero_offset_array)
DEALLOCATE(RV_offset_datasets_array)
DEALLOCATE(Chi_squared_array_total)
DEALLOCATE(Reduced_Chi_Squared_array_total)
DEALLOCATE(Parameter_array_total)
DEALLOCATE(index_array)

DEALLOCATE(sorted_min_chi_squared_mc)
DEALLOCATE(sorted_min_reduced_chi_squared_mc)
DEALLOCATE(sorted_mc_best_vsini_array)
DEALLOCATE(sorted_mc_best_lambda_array)
DEALLOCATE(sorted_mc_best_RV_offset_datasets_array)
DEALLOCATE(sorted_mc_best_orbital_period_array)
DEALLOCATE(sorted_mc_best_JD_time_mid_transit_array)
DEALLOCATE(sorted_mc_best_Mp_array)
DEALLOCATE(sorted_mc_best_Rp_array)
DEALLOCATE(sorted_mc_best_Ecc_array)
DEALLOCATE(sorted_mc_best_Inc_array)
DEALLOCATE(sorted_mc_best_omega_arg_periastron_array)
DEALLOCATE(sorted_mc_best_RV_zero_offset_array)

cmd1 = "say 'Finished processing data.'"
cmd1 = TRIM(cmd1)
call execute_command_line(cmd1, wait=.false., exitstat=exitstat, cmdstat=cmdstat, cmdmsg=cmdmsg)

!Purge memory? Works in OS X 10.8 or earlier. 10.9 and later require sudo purge and password.
!cmd1 = "purge"
!cmd1 = TRIM(cmd1)
!call execute_command_line(cmd1, wait=.false., exitstat=exitstat, cmdstat=cmdstat, cmdmsg=cmdmsg)

CONTAINS

!Subroutine for ranking an array in ascending order. 
SUBROUTINE Selection_sort(a,index_array)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: a(:)
INTEGER(4), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: index_array
DOUBLE PRECISION :: temp_array(size(a))
INTEGER(4) :: i, minIndex, temp2
DOUBLE PRECISION :: temp
Allocate(index_array(SIZE(a)))
 
DO i = 1, SIZE(a)
   index_array(i) = i
   temp_array(i) = a(i)
END DO
 
DO i = 1, SIZE(a) - 1
   minIndex = MINLOC(temp_array(i:), 1) + i - 1
   IF (temp_array(i) > temp_array(minIndex)) THEN
      temp = temp_array(i)
      temp2 = index_array(i)
      temp_array(i) = temp_array(minIndex)
      index_array(i) = index_array(minIndex)
      temp_array(minIndex) = temp
      index_array(minIndex) = temp2
   END IF
END DO

!PRINT *, 'Selection_sort routine invoked'

END SUBROUTINE Selection_sort

END PROGRAM RM_effect_mc_V05

module ran_mod 
! module contains three functions 
! ran1 returns a uniform random number between 0-1 
! spread returns random number between min - max 
! normal returns a normal distribution

contains 
    function ran1()  !returns random number between 0 - 1 
        use numz 
        implicit none 
        real(b8) ran1,x 
        CALL init_random_seed()
        call random_number(x) ! built in fortran 90 random number function 
        ran1=x 
    end function ran1

    function spread_rand(min,max)  !returns random number between min - max 
        use numz 
        implicit none 
        real(b8) spread_rand 
        real(b8) min,max 
        spread_rand=(max - min) * ran1() + min 
    end function spread_rand

    function normal(mean,sigma) !returns a normal distribution 
        use numz 
        implicit none 
        real(b8) normal,tmp 
        real(b8) mean,sigma 
        integer flag 
        real(b8) fac,gsave,rsq,r1,r2 
        save flag,gsave 
        data flag /0/ 
        if (flag.eq.0) then 
        rsq=2.0_b8 
            do while(rsq.ge.1.0_b8.or.rsq.eq.0.0_b8) ! new from for do 
                r1=2.0_b8*ran1()-1.0_b8 
                r2=2.0_b8*ran1()-1.0_b8 
                rsq=r1*r1+r2*r2 
            enddo 
            fac=sqrt(-2.0_b8*log(rsq)/rsq) 
            gsave=r1*fac 
            tmp=r2*fac 
            flag=1 
        else 
            tmp=gsave 
            flag=0 
        endif 
        normal=tmp*sigma+mean 
        return 
    end function normal

end module ran_mod
!end module

module numz 
  integer, parameter:: b8 = selected_real_kind(14)
end module

subroutine init_random_seed()
   !include ’system.inc’
   use iso_fortran_env, only: int64
   implicit none
   integer, allocatable :: seed(:)
   integer :: i, n, un, istat, dt(8), pid
   integer(int64) :: t
   integer(4) :: getpid
          
   call random_seed(size = n)
   !print *, 'n: ', n
   allocate(seed(n))
   ! First try if the OS provides a random number generator
   open(newunit=un, file="/dev/urandom", access="stream", &
        form="unformatted", action="read", status="old", iostat=istat)
   if (istat == 0) then
      read(un) seed
      close(un)
   else
      ! Fallback to XOR:ing the current time and pid. The PID is
      ! useful in case one launches multiple instances of the same
      ! program in parallel.
      call system_clock(t)
      if (t == 0) then
         call date_and_time(values=dt)
         t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
              + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
              + dt(3) * 24_int64 * 60 * 60 * 1000 &
              + dt(5) * 60 * 60 * 1000 &
              + dt(6) * 60 * 1000 + dt(7) * 1000 &
              + dt(8)
      end if
      !pid = getpid()
      pid = un
      t = ieor(t, int(pid, kind(t)))
      do i = 1, n
         seed(i) = lcg(t)
      end do
   end if
   call random_seed(put=seed)
          contains
   ! This simple PRNG might not be good enough for real work, but is
   ! sufficient for seeding a better PRNG.
   function lcg(s)
     integer :: lcg
     integer(int64) :: s
     if (s == 0) then
        s = 104729
     else
        s = mod(s, 4294967296_int64)
     end if
     s = mod(s * 279470273_int64, 4294967291_int64)
     lcg = int(mod(s, int(huge(0), int64)), kind(0))
   end function lcg
end subroutine init_random_seed

   FUNCTION BESSJ(N3,X)

!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows. 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.

      PARAMETER (IACC = 40,BIGNO = 1.0D10, BIGNI = 1.0D-10)
      DOUBLE PRECISION :: X,BESSJ,BESSJ0,BESSJ1,TOX,BJM,BJ,BJP,SUM
      IF (N3 .EQ. 0) THEN
      BESSJ = BESSJ0(X)
      RETURN
      ENDIF
      IF (N3 .EQ. 1) THEN
      BESSJ = BESSJ1(X)
      RETURN
      ENDIF
      IF (X == 0.) THEN
      BESSJ = 0.
      RETURN
      ENDIF
      TOX = 2./X
      IF (X > FLOAT(N3)) THEN
      BJM = BESSJ0(X)
      BJ  = BESSJ1(X)
      DO 11 J = 1,N3-1
      BJP = J*TOX*BJ-BJM
      BJM = BJ
      BJ  = BJP
   11 CONTINUE
      BESSJ = BJ
      ELSE
      M = 2*((N3+INT(SQRT(FLOAT(IACC*N3))))/2)
      BESSJ = 0.
      JSUM = 0
      SUM = 0.
      BJP = 0.
      BJ  = 1.
      DO 12 J = M, 1, -1
      BJM = J*TOX*BJ-BJP
      BJP = BJ
      BJ  = BJM
      IF (ABS(BJ) > BIGNO) THEN
      BJ  = BJ*BIGNI
      BJP = BJP*BIGNI
      BESSJ = BESSJ*BIGNI
      SUM = SUM*BIGNI
      ENDIF
      IF (JSUM /= 0) SUM = SUM+BJ
      JSUM = 1-JSUM
      IF (J == N3) BESSJ = BJP
   12 CONTINUE
      SUM = 2.*SUM-BJ
      BESSJ = BESSJ/SUM
      ENDIF
      RETURN
      END

   FUNCTION BESSJ0(X)
      DOUBLE PRECISION :: X,BESSJ0,AX,FR,FS,Z,FP,FQ,XX

!     This subroutine calculates the First Kind Bessel Function of
!     order 0, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.

      DOUBLE PRECISION :: Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6  &
                 ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4, &
      -.2073370639D-5,.2093887211D-6 /
      DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3, &
      -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
      DATA R1,R2,R3,R4,R5,R6 /57568490574.D0,-13362590354.D0, &
      651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0 /
      DATA S1,S2,S3,S4,S5,S6 /57568490411.D0,1029532985.D0, &
      9494680.718D0,59272.64853D0,267.8532712D0,1.D0 /
      IF(X == 0.D0) GO TO 1
      AX = ABS (X)
      IF (AX < 8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ0 = FR/FS
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-.785398164
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ0 = SQRT(.636619772/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
      ENDIF
      RETURN
    1 BESSJ0 = 1.D0
      RETURN
      END

! ---------------------------------------------------------------------------
   FUNCTION BESSJ1 (X)
      DOUBLE PRECISION X,BESSJ1,AX,FR,FS,Z,FP,FQ,XX
!     This subroutine calculates the First Kind Bessel Function of
!     order 1, for any real number X. The polynomial approximation by
!     series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
!     REFERENCES:
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.
      DOUBLE PRECISION :: Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6  &
               ,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4,  &
      .2457520174D-5,-.240337019D-6 /,P6 /.636619772D0 /
      DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,   &
      .8449199096D-5,-.88228987D-6,.105787412D-6 /
      DATA R1,R2,R3,R4,R5,R6 /72362614232.D0,-7895059235.D0, & 
      242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0 /
      DATA S1,S2,S3,S4,S5,S6 /144725228442.D0,2300535178.D0, &
      18583304.74D0,99447.43394D0,376.9991397D0,1.D0 /

      AX = ABS(X)
      IF (AX < 8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ1 = X*(FR/FS)
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-2.35619491
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ1 = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*SIGN(S6,X)
      ENDIF
      RETURN
      END