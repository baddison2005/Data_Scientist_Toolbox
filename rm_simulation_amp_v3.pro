PRO rm_simulation_amp_v3

;This program simulates the RV velocity curve, Rossiter_McLaughlin effect, transit light curve, and the secondary eclipse or occultation event.
;This is done with input parameters the user manually enters into this code. No parameter fitting or chi squared analysis is done here.

;Define constants used in the program.

Rss = 6.96D8                        ;Radius of our Sun (in meters).
Rj = 7.15D7                         ;Radius of Jupiter (in meters).
RE = 6.38D6                         ;Radius of Earth (in meters if planet is given in Earth radii).
AU = 1.50D11                        ;One astronomical unit in meters.
pi = !PI                            ;Constant value of pi.
G = 6.67D-11                        ;Gravitation constant.
Mss = 1.99D30                       ;Mass of the sun in kilograms.
Me = 5.97e24                        ;Mass of the Earth in kg.
Mj = 2.00e27                        ;Mass of Jupiter in kg.

;*****************************************Change Parameters Here***************************************************************
save_model_directory = '/Users/z3345211/Dropbox/RV_RM_transit_program/exofast/WASP-47/exosam_simulator/'
Planet_name = 'WASP-47'
save_data_directory = save_model_directory + Planet_name + '/run01/'
vel_overplot_name = save_model_directory + Planet_name + '/' + 'WASP-47_with_my_data_idl.csv'

Jupiter_units = 'Y'
Bessel_function_exit = 1.0D-8                              ;Exit condition for when the change in the eccentric anomaly after each iteration reaches an insignificant value.
;Model_plot_interval = 0.10                                 ;The interval change in eccentric anomaly.
;Phase_angle_start = 90                                    ;The starting location for modeling the planets orbit.
JD_time_mid_transit = 2457007.932132D0
;JD_time_mid_transit = 2455545.23479D0                      ;The mid transit time in JD. Will be normalized to 0 in model at mid transit point or when eccentric anomaly angle equals 0.
;JD_time_peri = 2455882D0                                   ;Estimate the amount of time (in seconds) after the passage of periastron has occurred when the planet
                                                           ;is at mid transit. Usually determined through RV modeling of data.
;Mean_anomaly_transit =  303.0D0                            ;Alternatively estimate the mean anomaly at transit if the time since passage of periastron at mid
                                                           ;transit is not known.
;Use_mean_anomaly = 'Y'                                     ;Y or N. If set to yes then the time of passage of periastron at mid transit is calculated.                                                                                                                     
Time_plot_interval = 10D0                                  ;Time interval plot seconds.
Ms_solar = 1.04D0                                          ;Mass of star in solar masses.
Long_ascend_node = 90                                      ;Assume 90 since we are only measuring the relative or projected inclination of orbit.
arg_periastron = 18.6d0                                       ;Argument of the periastron.
omega_arg_periastron = arg_periastron + Long_ascend_node   ;Omega argument of the periastron or longitude of pericentre.
M_pixels = 25.00                                           ;M number of pixels. Used in the rare case Rp is a significant fraction of Rs.
Inc_interval = 1.0d0                                         ;Orbital inclination angle interval.
Inc_begin = 89.41d0                                          ;Orbital inclination angle begin.
Inc_end = 89.41d0                                            ;Orbital inclination angle final.
Ecc = 0.00D0                                                ;Eccentricity of the planet.
;vsini = 19100.00D0                                        ;Stellar rotation rate in meters per second.
vsini = 3000.0D0                                         ;Stellar rotation rate in meters per second.
vmacro = 0.00D0                                            ;Stellar macroturbulence velocity in meters per second. Set to zero if unknown.
vmicro = 0.00D0                                         ;Stellar microturbulence velocity in meters per second. Set to zero if unknown.
Albedo = 0.30                                              ;The bond albedo of the planet.
stellar_rotation_angle_interval = 10                       ;Spin orbit alignment angle interval.
stellar_rotation_angle_begin = 0                         ;Spin orbit alignment angle begin.
stellar_rotation_angle_end = 350                           ;Spin orbit alignment angle end.
Rs_solar = 1.15D0                                         ;Radius of star (solar radii).
u = 0.466D0                                                ;Linear limb darkening coefficient.
Time_bef_aft = 5000                                        ;The amount of time in seconds to plot before and after the mid transit point,
                                                           ;RM, occultation events in zoomed in mode.
Mplan = 1.06D0                                            ;Mass of the planet in Jupiter or Earth masses.
Rplan = 1.158D0                                            ;Radius of the planet in Jupiter or Earth radii.
Orbital_period = 4.1591282D0                               ;Orbital period in days.
RV_zero_offset = -0.0D0                                 ;The velocity offset to be applied to the data.
title_on = 'T'                                             ;Display title on top of graph? T for true and F or any other character besides T for false.
overplot_velocities = 'F'                                  ;Overplot measured velocities on RM plots? T for true and F or any other character besides T for false.

;***********************************************************************************************************************************

;Create output directory if it doesn't exist.
   
FILE_MKDIR, save_data_directory

Rs = Rss * Rs_solar

IF Jupiter_units EQ 'Y' THEN BEGIN
   Mp = Mplan * Mj
   Rp = Rplan * Rj
ENDIF

IF Jupiter_units EQ 'N' THEN BEGIN
   Mp = Mplan * Me
   Rp = Rplan * RE
ENDIF

Ms = Ms_solar * Mss
;total_interval = ROUND(360/Model_plot_interval) + 1
Orbital_period = Orbital_period*(3600.0D0*24.0D0)
PRINT, "Orbital_period", Orbital_period
total_interval = ROUND(Orbital_period/Time_plot_interval) + 1

array_size = total_interval

RM_effect_array = dblarr(2,array_size)
RV_array = dblarr(2,array_size)
Transit_LC_array = dblarr(2,array_size)
Planet_star_distance_array  = dblarr(2,array_size)
Timestep = dblarr(array_size)
Phase_array = dblarr(array_size)
Phase_array_n = dblarr(array_size)
Phase_act_array = dblarr(array_size)
Xpos_array = dblarr(array_size)
Ypos_array = dblarr(array_size)
Zpos_array = dblarr(array_size)
True_anomaly_array = dblarr(array_size)
ecc_anomaly_array = dblarr(array_size)

Number_iterations = ROUND((((stellar_rotation_angle_end - stellar_rotation_angle_begin)/stellar_rotation_angle_interval) + 1)*(((Inc_end - Inc_begin)/Inc_interval) + 1))
Number_fit = 2

Num_Inc = ROUND(((Inc_end - Inc_begin)/Inc_interval) + 1)
Num_stellar_rotation_angle = ROUND(((stellar_rotation_angle_end - stellar_rotation_angle_begin)/stellar_rotation_angle_interval) + 1)

Total_L = 1.0                                ;Set the total amount of flux from the star to 1.
Total_RM = 0.0                               ;Set the Rossiter-McLaughlin anomaly to zero.
RV = 0.0                                     ;Set the RV to zero.
Rorb = ((Orbital_period^2*G*(Ms + Mp))/(4.0d0*pi^2))^(1.0d0/3.0d0)     ;Semi-major axis.
Rs2 = Rs^2                                   ;Variable to speed up calculations.
Rp2 = Rp^2                                   ;Variable to speed up calculations.
L_total_star = 1.0                           ;Normalize total luminosity of star to 1.
Io = 1.0D / (pi*Rs^2*(1.0D0-(u/3.0D0)))      ;The initial light intensity equation with limb darkening
                                             ;(normalize Io such that total star luminosity is 1).
Aplan = pi*Rp2                               ;Surface area of the planet.  
vturb = SQRT(vmicro^2 + vmacro^2)            ;Velocity width of spectral line due to mechanisms other than rotation (i.e. micro and macro turbulence).

True_anomaly_start = pi - omega_arg_periastron*(pi/180D0)
PRINT, "True_anomaly_start", True_anomaly_start*(180D0/pi)

IF True_anomaly_start GE pi THEN BEGIN
   ecc_anomaly_start = 2D0*atan(tan(True_anomaly_start/2D0)*sqrt((1-Ecc)/(1+Ecc))) + 2D0*pi 
ENDIF ELSE IF True_anomaly_start LE -pi THEN BEGIN
   ecc_anomaly_start = 2D0*atan(tan(True_anomaly_start/2D0)*sqrt((1-Ecc)/(1+Ecc))) - 2D0*pi
ENDIF ELSE BEGIN
   ecc_anomaly_start = 2D0*atan(tan(True_anomaly_start/2D0)*sqrt((1-Ecc)/(1+Ecc)))
ENDELSE

;Calculate the amount of time (in seconds) after the passage of periastron has occurred when the planet is at mid transit.
;This is determined from the transit mid time and the argument of periastron.
IF omega_arg_periastron GT 180.0D0 THEN BEGIN
   Mean_anomaly_transit = (2D0*pi) - (omega_arg_periastron*(pi/180D0)) + pi
ENDIF ELSE BEGIN
   Mean_anomaly_transit = pi - (omega_arg_periastron*(pi/180D0))
ENDELSE

JD_time_peri = JD_time_mid_transit - ((Mean_anomaly_transit*Orbital_period)/(2*pi))

time_peri_passage = (JD_time_mid_transit - JD_time_peri)*24.0D0*3600.0D0

IF Ecc EQ 0 THEN BEGIN
   Time_start = ((ecc_anomaly_start*Orbital_period)/(2D0*pi)) + DOUBLE(time_peri_passage) 
ENDIF ELSE BEGIN
   Time_start = (((ecc_anomaly_start - (Ecc*sin(ecc_anomaly_start)))*Orbital_period)/(2D0*pi)) + DOUBLE(time_peri_passage)
ENDELSE 
PRINT, "ecc_anomaly_start", ecc_anomaly_start*(180D0/pi)
PRINT, "Time_start", Time_start

True_anomaly_transit = ((3D0*pi)/2D0) - (omega_arg_periastron*(pi/180D0))
PRINT, "True_anomaly_transit", True_anomaly_transit*(180D0/pi)

IF True_anomaly_transit GE pi THEN BEGIN
   ecc_anomaly_transit = 2D0*atan(tan(True_anomaly_transit/2D0)*sqrt((1-Ecc)/(1+Ecc))) + 2D0*pi 
ENDIF ELSE IF True_anomaly_transit LE -pi THEN BEGIN
   ecc_anomaly_transit = 2D0*atan(tan(True_anomaly_transit/2D0)*sqrt((1-Ecc)/(1+Ecc))) - 2D0*pi
ENDIF ELSE BEGIN
   ecc_anomaly_transit = 2D0*atan(tan(True_anomaly_transit/2D0)*sqrt((1-Ecc)/(1+Ecc)))
ENDELSE

IF Ecc EQ 0 THEN BEGIN
   Time_transit = ((ecc_anomaly_transit*Orbital_period)/(2D0*pi)) + DOUBLE(time_peri_passage) 
ENDIF ELSE BEGIN
   Time_transit = (((ecc_anomaly_transit - (Ecc*sin(ecc_anomaly_transit)))*Orbital_period)/(2D0*pi)) + DOUBLE(time_peri_passage)
ENDELSE 
PRINT, "ecc_anomaly_transit", ecc_anomaly_transit*(180D0/pi)
PRINT, "Time_transit", Time_transit

;Inclination array used when Inc is iterated between Inc_end and Inc_begin.  Inc values are stored in this array.
Inc_array = dblarr(Num_Inc)
;stellar_rotation_angle array used when stellar_rotation_angle is iterated between stellar_rotation_angle_end and stellar_rotation_angle_begin.  stellar_rotation_angle values are stored in this array.
stellar_rotation_angle_array = dblarr(Num_stellar_rotation_angle)

Parameter_array = dblarr(8,Number_iterations) 
Number_array = Intarr(Number_iterations)

IF overplot_velocities EQ 'T' THEN BEGIN
   vel_overplot_data = read_csv(vel_overplot_name, COUNT=num_rv)
ENDIF

;Apply velocity offset to RV data.
IF overplot_velocities EQ 'T' THEN BEGIN
   vel_overplot_data.field2 = vel_overplot_data.field2 + RV_zero_offset
ENDIF

;****************************************************************************************************************************** 

e = 0
k = 0  

Number_Inc = 0
Number_stellar_rotation_angle = 0
Number = 0

FOR Inc_loop = e, Num_Inc - 1 DO BEGIN
   f = 1
   Inc = Inc_end - (Inc_interval * Inc_loop)
   Inc_array(Number_Inc) = Inc         ;Put Inc values into the vsini array at the position Number_Inc.
   k = 0
   
   IF Ecc EQ 0 THEN BEGIN
      ;Maximum amplitude caused by the exoplanet in a circular orbit.
      RVamplitude = (Mp/Ms)*sqrt((G*(Ms + Mp))/Rorb)*sin(Inc*(pi/180D0))
   ENDIF ELSE BEGIN
      ;Maximum amplitude caused by the exoplanet in an eccentric orbit.
      RVamplitude = (1.0D0/sqrt(1.0D0 - Ecc^2))*(Mp/Ms)*sqrt((G*(Ms + Mp))/Rorb)*sin(Inc*(pi/180D0))
   ENDELSE 
   
   ;The second do loop iterates values over the JD_mid_time of the planet's orbit.
   FOR stellar_rotation_angle_loop = k, Num_stellar_rotation_angle - 1 DO BEGIN
      stellar_rotation_angle = stellar_rotation_angle_end - (stellar_rotation_angle_interval * stellar_rotation_angle_loop)
      ;If Number_stellar_rotation_angle + 1 is less then or equal to ((stellar_rotation_angle_end - stellar_rotation_angle_begin) 
      ;/stellar_rotation_angle_interval) + 2 then the value for stellar_rotation_angle is put into the stellar_rotation_angle_array at the 
      ;location Number_stellar_rotation_angle. 

      IF (Number_stellar_rotation_angle LT ROUND(((stellar_rotation_angle_end - stellar_rotation_angle_begin) $
                                              /stellar_rotation_angle_interval) + 1)) THEN BEGIN
         stellar_rotation_angle_array(Number_stellar_rotation_angle) = stellar_rotation_angle
      ENDIF
      Number_stellar_rotation_angle = Number_stellar_rotation_angle + 1           ;Increase Number_stellar_rotation_angle by one every time
      ;stellar_rotation_angle_loop runs.

      ;Inc values at each iteration are placed in the Parameter_array at (Number,1).
      Parameter_array(0,Number) = Inc      
      ;stellar_rotation_angle values at each iteration are placed in the Parameter_array at (number,2).
      Parameter_array(1,Number) = stellar_rotation_angle
      Number_array(Number) = Number
      
      Time_transit_start = 0                ;Set the transit start time to zero.
      Transit_start_no_inc = 0              ;Set the transit start time for 90 inc to zero.
      Time_occultation_start = 0            ;Set the occultation start time to zero.
      Time_transit_end = 0                  ;Set the transit end time to zero.
      Transit_end_no_inc = 0                ;Set the transit end time for 90 inc to zero.
      Time_occultation_end = 0              ;Set the occultation end time to zero.
      Time_mid_transit = 0                  ;Set the mid transit time to zero.
      Time_mid_occultation = 0              ;Set the mid occultation time to zero.
      Planet_star_distance = 0              ;Set the distance between center of the star to the center of the planet (orbital radius) to zero.      
      Transit_start_no_inc_position = 1
      transit_start_position = 1
      Occultation_start_position = 1
      transit_End_position = 1
      Occultation_End_position = 1
      transit_mid_position = 1
      occultation_mid_position = 1
      transit_End_position_no_inc = 1
      occultation_End_position = 1
      
      FOR Time_loop = 1, total_interval DO BEGIN          ;A loop that steps through the planetary orbit
         Time = Time_start + ((Time_loop - 1D0)*Time_plot_interval)
   
         IF Ecc EQ 0 THEN BEGIN
            ecc_anomaly = ((2D0*pi)/Orbital_period)*(Time - time_peri_passage)
         ENDIF ELSE BEGIN
            ;Calculate the eccentric anomaly using the functions found near the end of the program.
            sum_ecc_anomaly = 0
            ecc_anomaly_before = 100000
            FOR order = 1, 20 DO BEGIN
               ;Calculate the value of the Bessel function which is used to find the eccentric anomaly.
               Bessel_value = BESELJ(order*Ecc, order, /DOUBLE)
               IF order GT 1 THEN BEGIN
                  ecc_anomaly_before = sum_ecc_anomaly
               ENDIF
               sum_ecc_anomaly = sum_ecc_anomaly + ((2D0/order)*Bessel_value*sin(order*(((2D0*pi)/Orbital_period)*(Time - time_peri_passage))))
               ecc_anomaly_after = sum_ecc_anomaly
               IF order GT 1 AND ABS(ecc_anomaly_after - ecc_anomaly_before) LE Bessel_function_exit THEN BREAK 
            ENDFOR
            ecc_anomaly = (((2D0*pi)/Orbital_period)*(Time - time_peri_passage)) + sum_ecc_anomaly
         ENDELSE 
         True_anomaly = 2D0*(atan(tan(ecc_anomaly/2D0)*(sqrt((1D0 + Ecc)/(1D0 - Ecc)))))
      
         IF Ecc EQ 0 THEN BEGIN
            ;The time of the similation for a specific true_anomaly with the mid transit time equal to 0.
            Time_check = ((ecc_anomaly*Orbital_period)/(2.0D0*pi)) + DOUBLE(time_peri_passage)
            ;The distance between the center of the planet to the center of the star in a circular orbit.
            Planet_star_distance = Rorb
         ENDIF ELSE BEGIN
            ;The time of the similation for a specific true_anomaly.
            Time_check = (((ecc_anomaly - (Ecc*sin(ecc_anomaly)))*Orbital_period)/(2.0D0*pi)) + DOUBLE(time_peri_passage)
            ;The distance between the center of the planet to the center of the star in an eccentric orbit.
            Planet_star_distance = (Rorb*(1 - Ecc^2))/(1.0D0 + (Ecc*cos(True_anomaly)))
         ENDELSE
         Time_ref = Time - DOUBLE(Time_transit)
      
         ;The position of the planet on the x-axis.
         Xpos = Planet_star_distance*((-sin(True_anomaly)*sin((omega_arg_periastron)*(pi/180.0D0))) + (cos(True_anomaly)*cos((omega_arg_periastron)*(pi/180.0D0))))
         ;The position of the planet on the y-axis.
         Ypos = Planet_star_distance*((cos(True_anomaly + pi)*cos((Inc)*(pi/180D0))*sin((omega_arg_periastron)*(pi/180D0))) + (sin(True_anomaly + pi)*cos((Inc)*(pi/180D0))*cos((omega_arg_periastron)*(pi/180D0))))
         ;The position of the planet on the z-axis.
         Zpos = Planet_star_distance*((-cos((omega_arg_periastron)*(pi/180.0D0))*sin((Inc)*(pi/180.0D0))*sin(True_anomaly)) - (cos(True_anomaly)*sin((Inc)*(pi/180.0D0))*sin((omega_arg_periastron)*(pi/180.0D0))))
         Dist2 = Xpos^2 + Ypos^2                   ;Square of the planet-star apparent seperation.
         Distance_center = sqrt(Dist2)             ;Apparent seperation between the planet and the star.
         Lblocked = 0.0D0
         Lblocked2 = 0.0D0                         ;A variable for the Anomalous velocity equation.
         v_rm = 0.0D0                              ;Anomalous velocity of each pixel set to zero.
         Total_RM = 0.0D
         IF Xpos LE 0 AND Zpos GE 0 THEN BEGIN
            ;Planet is currently in quadrant three so add pi.
            Phase_orbit = atan(Xpos/Zpos) + pi
            ;Taking into account orbital inclination.
            Phase_angle_observed = atan(-(sqrt(Xpos^2 + Ypos^2))/Zpos) + pi
         ENDIF ELSE IF Xpos GE 0 AND Zpos GE 0 THEN BEGIN
            ;Planet is currently in quadrant four so add pi.
            Phase_orbit = atan(Xpos/Zpos) + pi
            ;Taking into account orbital inclination.
            Phase_angle_observed = atan((sqrt(Xpos^2 + Ypos^2))/Zpos) + pi
         ENDIF ELSE IF Xpos GE 0 AND Zpos LE 0 THEN BEGIN
            ;Planet is currently in quadrant one so add 2pi.
            Phase_orbit = 2*pi + atan(Xpos/Zpos)
            ;Taking into account orbital inclination.
            Phase_angle_observed = 2*pi + atan((sqrt(Xpos^2 + Ypos^2))/Zpos)
         ENDIF ELSE IF Xpos LE 0 AND Zpos LE 0 THEN BEGIN
            ;Planet is currently in quadrant two so add 2pi.
            Phase_orbit = atan(Xpos/Zpos)
            ;Taking into account orbital inclination.
            Phase_angle_observed = atan(-(sqrt(Xpos^2 + Ypos^2))/Zpos)
         ENDIF 

         True_phase = acos(sin(True_anomaly + (omega_arg_periastron*(pi/180.0D0)))*sin(Inc*(pi/180.0D0)))
         Phase_orbit_n = acos(sin(True_anomaly + (omega_arg_periastron*(pi/180.0D0))))
               
         ;If the planet is neither in front of the star (transit) or behind the star (occultation), then calculate the flux being reflected off the surface of the exoplanet
         ;based on its bond albedo, radius, phase angle, etc.
         IF Distance_center GT (Rs + Rp) THEN BEGIN
            Aplan = pi*Rp2                                               ;Surface area of the planet.
            Total_L = 1.0D0 + ((Aplan/(4.0D0*pi*Planet_star_distance^2))*Albedo*(0.5D0*(1.0D0+cos(True_phase))))
         ENDIF
   
         IF ABS(Xpos) LE (Rs + Rp) AND Zpos GT 0 AND Transit_start_no_inc EQ 0 THEN BEGIN             ;If the seperation between the disk of the star and planet @ an inclination of 90 are less then 
                                                                                                      ;Rs + Rp and Zpos is positive then the transit begins.
            Transit_start_no_inc = Time_ref                ;Indicates when the transit starts at 90 inc. If statement prevents this variable from being overwritten every time the loop runs.
            Transit_start_no_inc_position = Time_loop      ;Location in array.
         ENDIF
               
         IF Distance_center LE (Rs + Rp) AND Zpos GT 0 THEN BEGIN       ;If the seperation between the disk of the star and planet are less then Rs + Rp and Zpos is positive
                                                                        ;then the transit begins.
       
            IF Time_transit_start EQ 0 THEN BEGIN                                                            
               Time_transit_start = Time_ref               ;Indicates when the transit starts. If statement prevents this variable from being overwritten every time the loop runs.  
               transit_start_position = Time_loop
            ENDIF       
                                                                  
            IF (Rp2/Rs2) GE 0.030 THEN BEGIN
               SArr = fltarr(M_pixels, M_pixels)            ;Creates the M_pixels by M_pixels array.
               Radius_planet_array = M_pixels / 2           ;Radius of planet in the pixel array.
               Center_of_planet_x = M_pixels / 2            ;The center of the planet on the X-axis.
               Center_of_planet_y = M_pixels / 2            ;The center of the planet on the Y-axis.
               Pixel = Rp/Radius_planet_array               ;The number of meters per pixel.
               Area_pixel = Pixel^2                         ;The area of each pixel. 
               Io_Pixel = Io * Area_pixel                   ;Variable to speed up calculations (also physically represents the
                                                            ;luminosity of the brightest pixel).                                                             
                  
               FOR i=0,M_pixels-1 DO BEGIN
                  X_pixel = (pixel * i) + (Xpos - Rp)     ;Calculates the location of the pixel on the x-axis.
                  X_pixel2 = X_pixel^2
                  XpXp = X_pixel * Xpos                   ;temporary var for speed calculation
                  FOR j=0,M_pixels-1 DO BEGIN
                     Y_pixel = (pixel * j) + abs(Ypos) - Rp            ;Calculates the location of the pixel on the y-axis.
                     X_pixel_prime = (X_pixel*cos((stellar_rotation_angle*pi)/180D0)) + (Y_pixel*sin((stellar_rotation_angle*pi)/180D0))     ;Calculates the location of the pixel along the x-axis of the rotation axis of the star.
                     Dist_center_pixel = X_pixel2 + Y_pixel^2         ;squared distance of pixel from star
                     ;Calculates the values of the pixels according to how far away they are from the center of the star and the plane.
                     ;squared distance of pixel from planet using a limb darkening equation.
                     Dist_planet_pixel = Dist_center_pixel  - 2.0D0*(XpXp + Y_pixel*Ypos) + Dist2
                     Sub_planet_velocity = vsini*(X_pixel_prime/Rs)
                     ;Sub_vturb = vturb*(X_pixel_prime/Rs)                  ;The turbulence velocity blocked by planet.
                     Sub_vturb = vturb                                    ;The turbulence velocity blocked by planet. 
                     IF Dist_center_pixel LE Rs2 AND Dist_planet_pixel LE Rp2 THEN BEGIN
                        Lblocked2 = Io_Pixel*(1.0D0-u*(1.0D0-sqrt(1.0D0-(Dist_center_pixel/Rs2))))          ;First order limb darkening equation.
                        IF ~ FINITE(Lblocked2) THEN PRINT, 'Lblocked2 Overflow occurred'
                        Lblocked = Lblocked + Lblocked2                                         ;First order limb darkening equation.
                        v_rm = v_rm - ((Lblocked2*Sub_planet_velocity)*((((2.0D0*Sub_vturb^2)+(2.0D0*vsini^2))/((2.0D0*Sub_vturb^2) + $
                               vsini^2))^(3.0D0/2.0D0))*(1.0D0-((Sub_planet_velocity^2)/((2.0D0*Sub_vturb^2)+(vsini^2)))))
                        ;v_rm = v_rm - ((Lblocked2*Sub_planet_velocity)*(1.33D0-(0.483D0*(Sub_planet_velocity/vsini)^2)))      ;Anomalous velocity of each pixel.
                        IF ~ FINITE(v_rm) THEN PRINT, 'V_rm Overflow occurred'
                     ENDIF
                  ENDFOR
               ENDFOR
               Total_L = 1.0D0 - Lblocked
               Total_RM = 0.0 + v_rm                       ;Total anomalous velocity for all the pixels.
         
            ENDIF ELSE IF (Rp2/Rs2) LE 0.030 THEN BEGIN

               ;Calculates the location of the center of the planet along the x-axis of the rotation axis of the star.
               X_prime = (Xpos*cos((stellar_rotation_angle*pi)/180.0D0)) + (Ypos*sin((stellar_rotation_angle*pi)/180.0D0))     
               set_distance_center = Distance_center         ;The limb darkening equation will use this distance as long as the center of the planet is inside the radius of the star.
               Sub_planet_velocity = vsini*(X_prime/Rs)      ;Calculate the subplanetary velocity (the stellar velocity blocked by the planetary disc) from vsini times
                                                             ;the distance from the center of the planet to the stellar rotation axis divided by the radius of the star.
               ;Sub_vturb = vturb*(X_prime/Rs)                ;The turbulence velocity blocked by planet.
               Sub_vturb = vturb                             ;The turbulence velocity blocked by planet. 
               Io_planet = 0.0D0
               
               IF Distance_center LE (Rs + Rp) AND Distance_center GE (Rs - Rp) THEN BEGIN    ;Start the planet on the x-axis so that the planet is
                                                                                              ;just far enough away not to touch the disk of the star.                
                  dist_cent1_int = (Distance_center^2 + Rs2 - Rp2)/(2.0D0*Distance_center)
                  ;Location on the x-axis for the first intersection point.
                  X_int_1 = ((Xpos*dist_cent1_int)/Distance_center) + ((Ypos/Distance_center)*sqrt(Rs2 - dist_cent1_int^2)) 
                  ;Location on the y-axis for the first intersection point.
                  Y_int_1 = ((Ypos*dist_cent1_int)/Distance_center) - ((Xpos/Distance_center)*sqrt(Rs2 - dist_cent1_int^2))            
                  ;Location on the x-axis for the second intersection point.
                  X_int_2 = ((Xpos*dist_cent1_int)/Distance_center) - ((Ypos/Distance_center)*sqrt(Rs2 - dist_cent1_int^2))             
                  ;Location on the y-axis for the second intersection point.
                  Y_int_2 = ((Ypos*dist_cent1_int)/Distance_center) + ((Xpos/Distance_center)*sqrt(Rs2 - dist_cent1_int^2)) 
                  ;The limb darkening equation will use this distance if any part of the disc of the planet is outside the radius of the star.
                  ;This is the distance between the center of the star to the center of the area inside the star blocked by the planet.           
                  set_distance_center = ((abs(Distance_center) - Rp) + Rs)/2.0D0
                  ;Next the program calculates the length between the two intersection point.
                  Length = sqrt((X_int_1 - X_int_2)^2 + (Y_int_1 - Y_int_2)^2)
                  ;Calculate the angle between the Y position of the center of planet and the x position of the center of planet.
                  ;This is used to determine the center of the area inside the star blocked by the planet in the stellar rotational
                  ;axis coordinate system.
                  Theta_inside = atan(Ypos/Xpos)
                  ;Calculates the distance from the center of the star (along x axis) to the center of the area inside the star blocked by the planet.
                  ;the pi factor added in this equation guarneetes that Xpos_in has the correct sign.
                  Xpos_in = set_distance_center*cos(Theta_inside + pi)
                  Xpos_inside = Xpos_in
         
                  ;This makes sure Xpos_inside has the correct sign.
                  IF Xpos GE 0 THEN BEGIN
                     Xpos_inside = -Xpos_in
                  ENDIF
         
                  ;Calculates the distance from the center of the star (along y axis) to the center of the area inside the star blocked by the planet.
                  Ypos_inside = abs(set_distance_center*sin(Theta_inside))
            
                  ;Changes the x-coordinate to the stellar rotation axis by an angle formed between the orbital plane of the planet and the 
                  ;stellar roatation plane of the star.
                  x_prime_distance = (Xpos_inside*cos((stellar_rotation_angle*pi)/180.0D0)) + (Ypos_inside*sin((stellar_rotation_angle*pi)/180.0D0))
                  Sub_planet_velocity = vsini*(x_prime_distance/Rs)    ;Calculate the subplanetary velocity (the stellar velocity blocked by the planetary disc) from vsini times
                                                                       ;the distance from the center of the planet to the stellar rotation axis divided by the radius of the star.
                  ;Sub_vturb = vturb*(x_prime_distance/Rs)              ;The turbulence velocity blocked by planet.  
                  Sub_vturb = vturb                                    ;The turbulence velocity blocked by planet.               
                  Beta1 = 2.0D0*asin((0.5D0*Length)/Rp)                ;Angle used to calculate the partial area of the planet in of the star.
                  Alpha1 = 2.0D0*asin((0.5D0*Length)/Rs)               ;Angle used to calculate the partial area of the planet in of the star.
                    
                  IF Distance_center GE sqrt(Rs2 - Rp2) THEN BEGIN   
                     ;The surface area of the planet when the center is outside the disk of the star.        
                     Aplan = ((0.5D0 * Rs2 * (Alpha1 - sin(Alpha1))) + (0.5D0 * Rp2 * (Beta1 - sin(Beta1))))
                     ;Normalized area of the planet with the given limb darkening function
                     Io_planet = Io * Aplan
                  ENDIF
     
                  IF Distance_center LT sqrt(Rs2 - Rp2) THEN BEGIN
                     Aplan = ((pi*Rp2 + (0.5D0 * Rs2 * (Alpha1 - sin(Alpha1)))) - (0.5D0 * Rp2 * (Beta1 - sin(Beta1))))
                     ;The surface area of the planet when the center is inside the disk of the star.
                     Io_planet = Io * Aplan
                  ENDIF
               ENDIF

  
               IF Distance_center LT (Rs - Rp) THEN BEGIN
                  Aplan = (pi*Rp2)                        ;The surface area of the planet.
                  Io_planet = Io * Aplan
               ENDIF
  
               Lblocked = (Io_planet)*(1.0D0-(u*(1.0D0-sqrt(abs(1.0D0-(set_distance_center^2/Rs2))))))         ;The ratio of the area of the planet blocking the star to   
                                                                                                               ;the area of the star utilizing the first order limb darkening equation.  
               IF ~ FINITE(Lblocked) THEN PRINT, 'Lblocked Overflow occurred'    
               v_rm = - ((Lblocked*Sub_planet_velocity)*((((2.0D0*Sub_vturb^2.0D0)+(2.0D0*vsini^2.0D0))/((2.0D0*Sub_vturb^2.0D0) + $
                        vsini^2.0D0))^(3.0D0/2.0D0))*(1.0D0-((Sub_planet_velocity^2.0D0)/((2.0D0*Sub_vturb^2.0D0)+(vsini^2.0D0)))))  
                                                                                                 
               ;v_rm = - ((Lblocked*Sub_planet_velocity)*(1.33D0-(0.483D0*(Sub_planet_velocity/vsini)^2)))      ;Anomalous velocity of each pixel.
               IF ~ FINITE(v_rm) THEN PRINT, 'V_rm Overflow occurred'
               Total_L = 1.0D0 - Lblocked                                             ;Total amount of light blocked by the planet.
               Total_RM = 0.0D0 + v_rm                                                ;Total anomalous velocity.
            ENDIF
         
         ENDIF      

         IF Distance_center LE (Rs + Rp) AND Zpos LT 0 THEN BEGIN       ;If the seperation between the disk of the star and planet are less then Rs + Rp and Zpos is negative
                                                                        ;then the secondary transit (occulation) begins.
 
            IF Time_occultation_start EQ 0 THEN BEGIN                                                            
               Time_occultation_start = Time_ref               ;Indicates when the secondary transit or occultation starts. If statement prevents this variable from being overwritten every time the loop runs.  
               Occultation_start_position = Time_loop
            ENDIF 
               
            IF Distance_center GE (Rs - Rp) THEN BEGIN    ;Partial secondary transit.                
               dist_cent1_int = (Distance_center^2 + Rs2 - Rp2)/(2.0D0*Distance_center)
               Ypos2 = abs(Ypos)   ;Set Ypos to a positive value to avoid discontinuity.  
               ;Location on the x-axis for the first intersection point.
               X_int_1 = ((Xpos*dist_cent1_int)/Distance_center) + ((Ypos2/Distance_center)*sqrt(Rs2 - dist_cent1_int^2)) 
               ;Location on the y-axis for the first intersection point.
               Y_int_1 = ((Ypos2*dist_cent1_int)/Distance_center) - ((Xpos/Distance_center)*sqrt(Rs2 - dist_cent1_int^2))            
               ;Location on the x-axis for the second intersection point.
               X_int_2 = ((Xpos*dist_cent1_int)/Distance_center) - ((Ypos2/Distance_center)*sqrt(Rs2 - dist_cent1_int^2))             
               ;Location on the y-axis for the second intersection point.
               Y_int_2 = ((Ypos2*dist_cent1_int)/Distance_center) + ((Xpos/Distance_center)*sqrt(Rs2 - dist_cent1_int^2)) 
               ;The limb darkening equation will use this distance if any part of the disc of the planet is outside the radius of the star.
               ;This is the distance between the center of the star to the center of the area inside the planet blocked by the star.           
               set_distance_center = ((abs(Distance_center) - Rp) + Rs)/2.0D0
               ;Next the program calculates the length between the two intersection point.
               Length = sqrt((X_int_1 - X_int_2)^2 + (Y_int_1 - Y_int_2)^2)
               ;Calculate the angle between the Y position of the center of planet and the x position of the center of planet.
               Theta_inside = atan(Ypos2/Xpos)
               ;Calculates the distance from the center of the star (along x axis) to the center of the area inside the star blocked by the planet.
               ;the pi factor added in this equation guarneetes that Xpos_in has the correct sign.
               Xpos_in = set_distance_center*cos(Theta_inside + pi)
               Xpos_inside = Xpos_in
         
               ;This makes sure Xpos_inside has the correct sign.
               IF Xpos GE 0 THEN BEGIN
                  Xpos_inside = -Xpos_in
               ENDIF
         
               ;Calculates the distance from the center of the star (along y axis) to the center of the area inside the star blocked by the planet.
               Ypos_inside = set_distance_center*sin(Theta_inside)
               Beta1 = 2.0D0*asin((0.5D0*Length)/Rp)         ;Angle used to calculate the partial area of the planet in the star.
               Alpha1 = 2.0D0*asin((0.5D0*Length)/Rs)        ;Angle used to calculate the partial area of the planet in the star.
                    
               IF Distance_center GE sqrt(Rs2 - Rp2) THEN BEGIN   
                  ;The surface area of the planet that is not behind the star when the center of the planets disk is visible.        
                  Aplan = (pi*Rp2) - ((0.5D0 * Rp2 * (Beta1 - sin(Beta1))) + (0.5D0 * Rs2 * (Alpha1 - sin(Alpha1))))
               ENDIF
     
               IF Distance_center LT sqrt(Rs2 - Rp2) THEN BEGIN
                  ;The surface area of the planet that is not behind the star when the center of the planets disk is behind the star.
                  Aplan = (0.5D0 * Rp2 * (Beta1 - sin(Beta1)))
               ENDIF
            ENDIF
  
            IF Distance_center LT (Rs - Rp) THEN BEGIN
               Aplan = 0                        ;No flux coming from planet since it's behind the star.
            ENDIF
  
            Total_L = 1.0D0 + ((Aplan/(4.0D0*pi*Planet_star_distance^2))*Albedo*(0.5D0*(1.0D0+cos(True_phase))))                 ;Total amount of light reflected by the planet.
         ENDIF 
    
         IF Time_transit_start NE 0 AND Time_transit_end EQ 0 AND Distance_center GE (Rs + Rp) THEN BEGIN                                                            
            Time_transit_end = Time_ref               ;Indicates when the transit Ends. If statement prevents this variable from being overwritten every time the loop runs.
            transit_End_position = Time_loop
         ENDIF  
               
         IF Time_mid_transit EQ 0 AND Xpos GE 0 THEN BEGIN
            Time_mid_transit = Time_ref               ;Indicates when the mid transit point occurs. If statement prevents this variable from being overwritten every time the loop runs.
            transit_mid_position = Time_loop
         ENDIF
               
         IF Time_mid_occultation EQ 0 AND Xpos LE 0 AND Time_mid_transit NE 0 THEN BEGIN
            Time_mid_occultation = Time_ref           ;Indicates when the mid occultation point occurs. If statement prevents this variable from being overwritten every time the loop runs.
            occultation_mid_position = Time_loop
         ENDIF
               
         IF Transit_start_no_inc NE 0 and Transit_end_no_inc EQ 0 AND ABS(Xpos) GE (Rs + Rp) THEN BEGIN                                                            
            Transit_end_no_inc = Time_ref               ;Indicates when the transit Ends @ 90 inc. If statement prevents this variable from being overwritten every time the loop runs. 
            transit_End_position_no_inc = Time_loop
         ENDIF
               
         IF Time_occultation_start NE 0 AND Time_occultation_end EQ 0 AND Distance_center GE (Rs + Rp) THEN BEGIN                                                            
            Time_occultation_end = Time_ref               ;Indicates when the transit Ends. If statement prevents this variable from being overwritten every time the loop runs.  
            occultation_End_position = Time_loop
         ENDIF
       
         IF Ecc EQ 0 THEN BEGIN
            ;The radial velocity of the star which includes adding the RM effect for a circular orbit.
            RV = RVamplitude*(cos(pi + True_anomaly + ((omega_arg_periastron)*(pi/180.0D0)))) + Total_RM
         ENDIF ELSE BEGIN
            ;The radial velocity of the star which includes adding the RM effect for an eccentric orbit.
            RV = RVamplitude*(cos(pi + True_anomaly + ((omega_arg_periastron)*(pi/180.0D0))) + (Ecc*cos((omega_arg_periastron)*(pi/180.0D0) + pi))) +  Total_RM
         ENDELSE 
                                                    
         ;Put data into the Lightcurve array from the total luminosity.
         Transit_LC_array(1, Time_loop - 1) = Total_L 
         Transit_LC_array(0, Time_loop - 1) = Time_ref
         ;Put all the anomalous velocity data into the anomalous velocity array.
         RM_effect_array(1, Time_loop - 1) = Total_RM
         RM_effect_array(0, Time_loop - 1) = Time_ref
         ;Put all the RV data into the RV array.
         RV_array(1, Time_loop - 1) = RV
         RV_array(0, Time_loop - 1) = Time_ref
         ;Put all the orbital radius values into the orbital radius array.
         Planet_star_distance_array(1, Time_loop - 1) = Planet_star_distance
         Planet_star_distance_array(0, Time_loop - 1) = Time_ref
         ;Set timestep to start at negative time.  The mid point of transit is at zero and everything after that is positive time in seconds.
         Timestep(Time_loop - 1) = Time_ref
         Phase_array(Time_loop - 1) = Phase_orbit
         Phase_act_array(Time_loop - 1) = True_phase
         Phase_array_n(Time_loop - 1) = Phase_orbit_n
         True_anomaly_array(Time_loop - 1) = True_anomaly
         ecc_anomaly_array(Time_loop - 1) = ecc_anomaly
         Xpos_array(Time_loop - 1) = Xpos
         Ypos_array(Time_loop - 1) = Ypos
         Zpos_array(Time_loop - 1) = Zpos
      ENDFOR
      
      Time_transit_actual = Time_transit_end - Time_transit_start
      ;The length of the occultation which cannot be calculated ahead of time due to the planet being in an eccentric orbit.
      Time_occultation_actual = Time_occultation_end - Time_occultation_start
      Time_transit_90inc = Transit_end_no_inc - Transit_start_no_inc
           
      Parameter_array(2,Number) = MIN(Transit_LC_array(1, *))
      Parameter_array(3,Number) = MAX(RM_effect_array(1, *)) + ABS(MIN(RM_effect_array(1, *)))
      Parameter_array(4,Number) = 1.0D0 - MIN(Transit_LC_array(1, *))
      Parameter_array(5,Number) = Time_transit_actual  
      Parameter_array(6,Number) = MAX(RV_array(1, *)) + ABS(MIN(RV_array(1, *))) 
      Parameter_array(7,Number) = MAX(Transit_LC_array(1, *)) - 1.0D0
            
      Number_string = 0
      Number_string = STRTRIM(number, 2)
      
      number_circ = 17.0 ; the circle will be "created" with 17 data points (vertices)
      theta_circ = findgen(number_circ)/(number_circ-1.0)*360.0*!DtoR ;
      x_circ = 1.0*sin(theta_circ)
      y_circ = 1.0*cos(theta_circ)
      usersym, x_circ, y_circ
      
      ;Adjust RV data time if it's outside model boundary.
      IF overplot_velocities EQ 'T' THEN BEGIN
         overplot_vel_time = dblarr(num_rv)
         FOR t=0, num_rv - 1 DO BEGIN
            IF vel_overplot_data.FIELD1(t)*(24.0D0*3600.0D0) LT RV_array(0,0) THEN BEGIN
               overplot_vel_time(t) = RV_array(0,total_interval - 1) - abs(vel_overplot_data.FIELD1(t)*(24.0D0*3600.0D0) - RV_array(0,0))
            ENDIF ELSE BEGIN
               overplot_vel_time(t) = vel_overplot_data.FIELD1(t)*(24.0D0*3600.0D0)
            ENDELSE   
         ENDFOR
      ENDIF
      
      IF title_on EQ 'T' THEN BEGIN
         cgplot, Transit_LC_array(0,*)/60, Transit_LC_array(1,*), xrange=[Transit_LC_array(0,0)/60,Transit_LC_array(0,array_size - 1)/60], $
         yrange=[MIN(Transit_LC_array(1,*)),2.0 - MIN(Transit_LC_array(1,*))], title = 'Predicted lightcurve', xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Relative magnitude', $
         output = 'PS', outfilename = save_data_directory + 'Predictedlightcurve_whole_orbit'  + Number_string + '.ps', axiscolor = 'black', background = 'white', color = 'black'
      ENDIF ELSE BEGIN
         cgplot, Transit_LC_array(0,*)/60, Transit_LC_array(1,*), xrange=[Transit_LC_array(0,0)/60,Transit_LC_array(0,array_size - 1)/60], $
         yrange=[MIN(Transit_LC_array(1,*)),2.0 - MIN(Transit_LC_array(1,*))], xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Relative magnitude', $
         output = 'PS', outfilename = save_data_directory + 'Predictedlightcurve_whole_orbit'  + Number_string + '.ps', axiscolor = 'black', background = 'white', color = 'black'
      ENDELSE      
      
;      thisDevice = !D.Name
;      Set_Plot, 'Z', /COPY
;      Device, Set_Resolution=[1200,800], Z_Buffer=0
;      Erase
;      plot, Transit_LC_array(0,*)/60, Transit_LC_array(1,*), xrange=[Transit_LC_array(0,0)/60,Transit_LC_array(0,array_size - 1)/60], yrange=[MIN(Transit_LC_array(1,*)),2.0 - MIN(Transit_LC_array(1,*))], title = 'Predicted lightcurve', $
;      color=plotcolor1, THICK='2', CHARSIZE = 1, background=bkcolor, xtitle = 'Time (Minutes from JD 0)', ytitle = 'Relative magnitude'
;      snapshot = TVRD()
;      TVLCT, red, green, blue, /Get
;      Device, Z_Buffer=1
;      Set_Plot, thisDevice
;      image24 = BytArr(3, 1200, 800)
;      image24[0,*,*] = red[snapshot]
;      image24[1,*,*] = green[snapshot]
;      image24[2,*,*] = blue[snapshot]
;      Write_jpeg, save_data_directory + 'Predictedlightcurve_whole_orbit'  + Number_string + '.jpeg', image24, True=1, QUALITY=100

      IF title_on EQ 'T' THEN BEGIN
         cgplot, Transit_LC_array(0,*)/60, Transit_LC_array(1,*), xrange=[Transit_LC_array(0,0)/60,Transit_LC_array(0,array_size - 1)/60], $
         yrange=[0.999 - (MAX(Transit_LC_array(1,*) - 1.0)),1.0001 + (MAX(Transit_LC_array(1,*) - 1.0))], title = 'Predicted lightcurve', xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Relative magnitude', $
         output = 'PS', outfilename = save_data_directory + 'Predictedlightcurve_whole_orbitz' + Number_string + '.ps', axiscolor = 'black', background = 'white', color = 'black'
      ENDIF ELSE BEGIN
         cgplot, Transit_LC_array(0,*)/60, Transit_LC_array(1,*), xrange=[Transit_LC_array(0,0)/60,Transit_LC_array(0,array_size - 1)/60], $
         yrange=[0.999 - (MAX(Transit_LC_array(1,*) - 1.0)),1.0001 + (MAX(Transit_LC_array(1,*) - 1.0))], xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Relative magnitude', $
         output = 'PS', outfilename = save_data_directory + 'Predictedlightcurve_whole_orbitz' + Number_string + '.ps', axiscolor = 'black', background = 'white', color = 'black'
      ENDELSE

;      thisDevice = !D.Name
;      Set_Plot, 'Z', /COPY
;      Device, Set_Resolution=[1200,800], Z_Buffer=0
;      Erase
;      plot, Transit_LC_array(0,*)/60, Transit_LC_array(1,*), xrange=[Transit_LC_array(0,0)/60,Transit_LC_array(0,array_size - 1)/60], yrange=[0.999 - (MAX(Transit_LC_array(1,*) - 1.0)),1.0001 + (MAX(Transit_LC_array(1,*) - 1.0))], title = 'Predicted lightcurve', $
;      color=plotcolor1, THICK='2', CHARSIZE = 1, background=bkcolor, xtitle = 'Time (Minutes from HJD 0)', ytitle = 'Relative magnitude'
;      snapshot = TVRD()
;      TVLCT, red, green, blue, /Get
;      Set_Plot, thisDevice
;      Device, Z_Buffer=1
;      image24 = BytArr(3, 1200, 800)
;      image24[0,*,*] = red[snapshot]
;      image24[1,*,*] = green[snapshot]
;      image24[2,*,*] = blue[snapshot]
;      Write_jpeg, save_data_directory + 'Predictedlightcurve_whole_orbit1'  + Number_string + '.jpeg', image24, True=1, QUALITY=100

      IF title_on EQ 'T' THEN BEGIN
         cgplot, Transit_LC_array(0,*)/60, Transit_LC_array(1,*), xrange=[(Transit_start_no_inc - Time_bef_aft - Time_mid_transit)/60,(Transit_end_no_inc + Time_bef_aft - Time_mid_transit)/60], $
         yrange=[MIN(Transit_LC_array(1,*)),2.0 - MIN(Transit_LC_array(1,*))], title = 'Predicted lightcurve', xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Relative magnitude', $
         output = 'PS', outfilename = save_data_directory + 'Predictedlightcurveplot'  + Number_string + '.ps', axiscolor = 'black', background = 'white', color = 'black'
      ENDIF ELSE BEGIN
         cgplot, Transit_LC_array(0,*)/60, Transit_LC_array(1,*), xrange=[(Transit_start_no_inc - Time_bef_aft - Time_mid_transit)/60,(Transit_end_no_inc + Time_bef_aft - Time_mid_transit)/60], $
         yrange=[MIN(Transit_LC_array(1,*)),2.0 - MIN(Transit_LC_array(1,*))], xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Relative magnitude', $
         output = 'PS', outfilename = save_data_directory + 'Predictedlightcurveplot'  + Number_string + '.ps', axiscolor = 'black', background = 'white', color = 'black'
      ENDELSE

;      thisDevice = !D.Name
;      Set_Plot, 'Z', /COPY
;      Device, Set_Resolution=[1200,800], Z_Buffer=0
;      Erase
;      plot, Transit_LC_array(0,*)/60, Transit_LC_array(1,*), xrange=[(Transit_start_no_inc - Time_bef_aft - Time_mid_transit)/60,(Transit_end_no_inc + Time_bef_aft - Time_mid_transit)/60], yrange=[MIN(Transit_LC_array(1,*)),2.0 - MIN(Transit_LC_array(1,*))], title = 'Predicted lightcurve', $
;      color=plotcolor1, THICK='2', CHARSIZE = 1, background=bkcolor, xtitle = 'Time (Minutes from HJD 0)', ytitle = 'Relative magnitude'
;      snapshot = TVRD()
;      TVLCT, red, green, blue, /Get
;      Device, Z_Buffer=1
;      Set_Plot, thisDevice
;      image24 = BytArr(3, 1200, 800)
;      image24[0,*,*] = red[snapshot]
;      image24[1,*,*] = green[snapshot]
;      image24[2,*,*] = blue[snapshot]
;      Write_jpeg, save_data_directory + 'Predictedlightcurveplot'  + Number_string + '.jpeg', image24, True=1, QUALITY=100

      IF title_on EQ 'T' THEN BEGIN
         cgplot, Transit_LC_array(0,*)/60, Transit_LC_array(1,*), xrange=[(Time_occultation_start - Time_bef_aft - Time_mid_transit)/60,(Time_occultation_end + Time_bef_aft - Time_mid_transit)/60], $
         yrange=[0.9999,MAX(Transit_LC_array(1,*))], title = 'Predicted lightcurve', xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Relative magnitude', output = 'PS', $
         outfilename = save_data_directory + 'Occultationlightcurveplot' + Number_string + '.ps', axiscolor = 'black', background = 'white', color = 'black'
      ENDIF ELSE BEGIN
         cgplot, Transit_LC_array(0,*)/60, Transit_LC_array(1,*), xrange=[(Time_occultation_start - Time_bef_aft - Time_mid_transit)/60,(Time_occultation_end + Time_bef_aft - Time_mid_transit)/60], $
         yrange=[0.9999,MAX(Transit_LC_array(1,*))], xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Relative magnitude', output = 'PS', $
         outfilename = save_data_directory + 'Occultationlightcurveplot' + Number_string + '.ps', axiscolor = 'black', background = 'white', color = 'black'
      ENDELSE

;      thisDevice = !D.Name
;      Set_Plot, 'Z', /COPY
;      Device, Set_Resolution=[1200,800], Z_Buffer=0
;      Erase
;      plot, Transit_LC_array(0,*)/60, Transit_LC_array(1,*), xrange=[(Time_occultation_start - Time_bef_aft - Time_mid_transit)/60,(Time_occultation_end + Time_bef_aft - Time_mid_transit)/60], yrange=[0.9999,MAX(Transit_LC_array(1,*))], title = 'Predicted lightcurve', $
;      color=plotcolor1, THICK='2', CHARSIZE = 1, background=bkcolor, xtitle = 'Time (Minutes from HJD 0)', ytitle = 'Relative magnitude'
;      snapshot = TVRD()
;      TVLCT, red, green, blue, /Get
;      Device, Z_Buffer=1
;      Set_Plot, thisDevice
;      image24 = BytArr(3, 1200, 800)
;      image24[0,*,*] = red[snapshot]
;      image24[1,*,*] = green[snapshot]
;      image24[2,*,*] = blue[snapshot]
;      Write_jpeg, save_data_directory + 'Occultationlightcurveplot'  + Number_string + '.jpeg', image24, True=1, QUALITY=100

      IF title_on EQ 'T' THEN BEGIN   
         cgplot, RM_effect_array(0,*)/60, RM_effect_array(1,*), xrange=[(Transit_start_no_inc - Time_bef_aft - Time_mid_transit)/60,(Transit_end_no_inc + Time_bef_aft - Time_mid_transit)/60], $
         yrange=[MIN(RM_effect_array(1,*)),MAX(RM_effect_array(1,*))], title = 'RM-effect anomaly', xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Radial velocity anomaly m/s', output = 'PS', $
         outfilename = save_data_directory + 'RVanomaly' + Number_string + '.ps', axiscolor = 'black', background = 'white', color = 'black'   
      ENDIF ELSE BEGIN
         cgplot, RM_effect_array(0,*)/60, RM_effect_array(1,*), xrange=[(Transit_start_no_inc - Time_bef_aft - Time_mid_transit)/60,(Transit_end_no_inc + Time_bef_aft - Time_mid_transit)/60], $
         yrange=[MIN(RM_effect_array(1,*)),MAX(RM_effect_array(1,*))], xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Radial velocity anomaly m/s', output = 'PS', $
         outfilename = save_data_directory + 'RVanomaly' + Number_string + '.ps', axiscolor = 'black', background = 'white', color = 'black' 
      ENDELSE
   
;      thisDevice = !D.Name
;      Set_Plot, 'Z', /COPY
;      Device, Set_Resolution=[1200,800], Z_Buffer=0
;      Erase
;      plot, RM_effect_array(0,*)/60, RM_effect_array(1,*), xrange=[(Transit_start_no_inc - Time_bef_aft - Time_mid_transit)/60,(Transit_end_no_inc + Time_bef_aft - Time_mid_transit)/60], yrange=[MIN(RM_effect_array(1,*)),MAX(RM_effect_array(1,*))], $
;      color=plotcolor1, THICK='2', CHARSIZE = 1, background=bkcolor, xtitle = 'Time (Minutes from HJD 0)', ytitle = 'Radial velocity anomaly m/s'
;      snapshot = TVRD()
;      TVLCT, red, green, blue, /Get
;      Device, Z_Buffer=1
;      Set_Plot, thisDevice
;      image24 = BytArr(3, 1200, 800)
;      image24[0,*,*] = red[snapshot]
;      image24[1,*,*] = green[snapshot]
;      image24[2,*,*] = blue[snapshot]
;      Write_jpeg, save_data_directory + 'RVanomaly'  + Number_string + '.jpeg', image24, True=1, QUALITY=100

      IF title_on EQ 'T' THEN BEGIN
         cgplot, RV_array(0,*)/60, RV_array(1,*), xrange=[RV_array(0,0)/60,RV_array(0,array_size - 1)/60], yrange=[MIN(RV_array(1,*)),MAX(RV_array(1,*))], $
         title = 'Radial velocity curve', xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Radial velocity m/s', output = 'PS', $
         outfilename = save_data_directory + 'RVcurve' + Number_string + '.ps', axiscolor = 'black', background = 'white', color = 'black'
      ENDIF ELSE BEGIN
         cgplot, RV_array(0,*)/60, RV_array(1,*), xrange=[RV_array(0,0)/60,RV_array(0,array_size - 1)/60], yrange=[MIN(RV_array(1,*)),MAX(RV_array(1,*))], $
         xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Radial velocity m/s', output = 'PS', $
         outfilename = save_data_directory + 'RVcurve' + Number_string + '.ps', axiscolor = 'black', background = 'white', color = 'black'
      ENDELSE

      IF overplot_velocities EQ 'T' THEN BEGIN
         IF title_on EQ 'T' THEN BEGIN   
;            RM_overplot = Obj_New('cgOverPlot', RV_array(0,*)/60, RV_array(1,*), color='black')
;            cgplot, vel_overplot_data.FIELD1*(24.0D0*60.0D0), vel_overplot_data.FIELD2, xrange=[RV_array(0,0)/60,RV_array(0,array_size - 1)/60], yrange=[MIN(RV_array(1,*)),MAX(RV_array(1,*))], $
;            title = 'Radial velocity curve', xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Radial velocity m/s', output = 'PS', $
;            outfilename = save_data_directory + 'RVcurve_op' + Number_string + '.ps', axiscolor = 'black', background = 'white', psym=2, symcolor = 'blue', OPLOTS=[RM_overplot]
;            cgErrPlot, vel_overplot_data.FIELD2+vel_overplot_data.FIELD3, vel_overplot_data.FIELD2-vel_overplot_data.FIELD3, Color='red'   
;            Obj_Destroy, [RM_overplot]            
            set_plot, 'PS'
            Device, Filename=save_data_directory + 'RVcurve_op' + Number_string + '.ps', /landscape
            RM_overplot = Obj_New('cgOverPlot', RV_array(0,*)/60.0D0, RV_array(1,*), color='black', Thick=4)
            PLOTERROR, overplot_vel_time/(60.0D0), vel_overplot_data.FIELD2, vel_overplot_data.FIELD3, Errcolor='red', axiscolor='black', background='white', psym=2, $
            symcolor='blue', xrange=[RV_array(0,0)/60.0D0,RV_array(0,array_size - 1)/60.0D0], yrange=[MIN(RV_array(1,*))-100.0D0,MAX(RV_array(1,*))+100.0D0], title = 'Radial velocity curve', $
            xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Radial velocity m/s', ERRThick=4, charsize=2, symsize = 2.5, font='TIMES BOLD', OPLOTS=[RM_overplot]            
            device, /Close_File
            set_plot, 'X'
            Obj_Destroy, [RM_overplot]
            
            set_plot, 'PS'
            Device, Filename=save_data_directory + 'RVcurvezoom_op' + Number_string + '.ps', /landscape
            RM_overplot = Obj_New('cgOverPlot', RV_array(0,*)/60.0D0, RV_array(1,*), color='black', Thick=4)
            PLOTERROR, overplot_vel_time/(60.0D0), vel_overplot_data.FIELD2, vel_overplot_data.FIELD3, Errcolor='red', axiscolor='black', background='white', psym=2, $
            symcolor='blue', xrange=[(Transit_start_no_inc - Time_bef_aft - Time_mid_transit)/60,(Transit_end_no_inc + Time_bef_aft - Time_mid_transit)/60], $
            yrange=[MIN(RV_array(1,*))-100.0D0,MAX(RV_array(1,*))+100.0D0], title = 'Radial velocity curve', $
            xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Radial velocity m/s', ERRThick=4, charsize=2, symsize = 2.5, font='TIMES BOLD', OPLOTS=[RM_overplot]            
            device, /Close_File
            set_plot, 'X'
            Obj_Destroy, [RM_overplot]
         ENDIF ELSE BEGIN
;            RM_overplot = Obj_New('cgOverPlot', RV_array(0,*)/60, RV_array(1,*), color='black')
;            cgplot, vel_overplot_data.FIELD1*(24.0D0*60.0D0), vel_overplot_data.FIELD2, xrange=[RV_array(0,0)/60,RV_array(0,array_size - 1)/60], yrange=[MIN(RV_array(1,*)),MAX(RV_array(1,*))], $
;            xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Radial velocity m/s', output = 'PS', $
;            outfilename = save_data_directory + 'RVcurve_op' + Number_string + '.ps', axiscolor = 'black', background = 'white', psym=2, symcolor = 'blue', OPLOTS=[RM_overplot]
;            cgErrPlot, vel_overplot_data.FIELD2+vel_overplot_data.FIELD3, vel_overplot_data.FIELD2-vel_overplot_data.FIELD3, Color='red'   
;            Obj_Destroy, [RM_overplot]
            set_plot, 'PS'
            Device, Filename=save_data_directory + 'RVcurve_op' + Number_string + '.ps', /landscape
            RM_overplot = Obj_New('cgOverPlot', RV_array(0,*)/60, RV_array(1,*), color='black', Thick=4)
            PLOTERROR, overplot_vel_time/(60.0D0), vel_overplot_data.FIELD2, vel_overplot_data.FIELD3, Errcolor='red', axiscolor='black', background='white', psym=2, $
            symcolor='blue', xrange=[RV_array(0,0)/60,RV_array(0,array_size - 1)/60], yrange=[MIN(RV_array(1,*))-100.0D0,MAX(RV_array(1,*))+100.0D0], $
            xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Radial velocity m/s', ERRThick=4, charsize=2, symsize = 2.5, font='TIMES BOLD', OPLOTS=[RM_overplot]            
            device, /Close_File
            set_plot, 'X'
            Obj_Destroy, [RM_overplot]
            
            set_plot, 'PS'
            Device, Filename=save_data_directory + 'RVcurvezoom_op' + Number_string + '.ps', /landscape
            RM_overplot = Obj_New('cgOverPlot', RV_array(0,*)/60, RV_array(1,*), color='black', Thick=4)
            PLOTERROR, overplot_vel_time/(60.0D0), vel_overplot_data.FIELD2, vel_overplot_data.FIELD3, Errcolor='red', axiscolor='black', background='white', psym=2, $
            symcolor='blue', xrange=[(Transit_start_no_inc - Time_bef_aft - Time_mid_transit)/60,(Transit_end_no_inc + Time_bef_aft - Time_mid_transit)/60], $
            yrange=[MIN(RV_array(1,*))-100.0D0,MAX(RV_array(1,*))+100.0D0], xtitle = 'Time (Minutes from MJD 0)', ytitle = 'Radial velocity m/s', $
            ERRThick=4, charsize=2, symsize = 2.5, font='TIMES BOLD', OPLOTS=[RM_overplot]            
            device, /Close_File
            set_plot, 'X'
            Obj_Destroy, [RM_overplot]
         ENDELSE
      ENDIF

;      thisDevice = !D.Name
;      Set_Plot, 'Z', /COPY
;      Device, Set_Resolution=[1200,800], Z_Buffer=0
;      Erase
;      plot, RV_array(0,*)/60, RV_array(1,*), xrange=[RV_array(0,0)/60,RV_array(0,array_size - 1)/60], yrange=[MIN(RV_array(1,*)),MAX(RV_array(1,*))], $
;      color=plotcolor1, THICK='2.5', CHARSIZE = 2.5, background=bkcolor, xtitle = 'Time (Min from HJD 0)', ytitle = 'Radial velocity anomaly m/s'
;      snapshot = TVRD()
;      TVLCT, red, green, blue, /Get
;      Device, Z_Buffer=1
;      Set_Plot, thisDevice
;      image24 = BytArr(3, 1200, 800)
;      image24[0,*,*] = red[snapshot]
;      image24[1,*,*] = green[snapshot]
;      image24[2,*,*] = blue[snapshot]
;      Write_jpeg, save_data_directory + 'RVcurve'  + Number_string + '.jpeg', image24, True=1, QUALITY=100
      
      Number = Number + 1 
      
      PRINT, 'Iteration number is: ', Number
      PRINT, 'Iterations left: ', Number_iterations - Number
      
   ENDFOR
   
   Number_Inc = Number_Inc + 1         ;Increase Number_vsini by 1 each time this do loop runs.
   
ENDFOR

l=0

;This part outputs all the values found in the model in a text file.
openw,1,save_data_directory + 'outputresults.txt'
PrintF, 1, 'run #', '|', 'Vsin(i)', '|', 'Spin-orbit (deg)', '|', 'Inclination (deg)', '|', 'Stellar Radius', '|', $
'Planet radius', '|', 'semi-major axis', '|', 'Ecc', '|', 'Orbital Period', '|', 'RM-amplitude', '|', 'RV-amplitude', $
'|', 'Transit length (sec)', '|', 'Transit depth', '|', 'Percent flux drop', '|', 'Occultation depth', $
FORMAT='(A6, x, A1, x, A10, 3x, A1, x, A17, 2x, A1, x, A19, 2x, A1, x, A15, x, A1, x, A14, 2x, A1, x, A16, x, A1, x, A9, 4x, A1, x, A15, x, A1, x, A14, 2x, A1, x, A15, 2x, A1, x, A22, 3x, A1, x, A14, x, A1, x, A18, x, A1, x, A18)'
FOR l = 0, Number - 1 DO BEGIN
   PrintF, 1, l, Vsini, Parameter_array(1,l), Parameter_array(0,l), Rs/Rss, Rp/Rj, Rorb/Au, Ecc, Orbital_period/(3600.0D0*24.0D0), Parameter_array(3,l), Parameter_array(6,l), Parameter_array(5,l), Parameter_array(2,l), $
   Parameter_array(4,l)*100D0, Parameter_array(7,l), FORMAT='(I3, 8x, F-10.2, 10x, F-10.2, 7x, F-15.8, 8x, F-15.8, 4x, F-14.8, 6x, F-10.4, 4x, F-15.8, x, F-16.8, x, F-16.8, 3x, F-16.8, 3x, F-20.8, 4x, F-16.8, 4x, F-16.8, 5x, F-12.8)'
ENDFOR
Close, 1
       
END