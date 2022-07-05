! ***********************************************************************
!
!   Copyright (C) 2010-2019  Bill Paxton & The MESA Team
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful,
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************

      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use math_lib

      implicit none

      real(dp) :: original_diffusion_dt_limit
      real(dp) :: postAGB_check = 0.0
      real(dp) :: rot_set_check = 0.0
      logical :: wd_diffusion = .false.
      real(dp) :: X_C_init, X_N_init

      integer :: time0, time1, clock_rate

      ! Parameters necessary for axion energy-loss
      real(dp) :: axion_g10
      real(dp) :: g_x_array(6)
      real(dp) :: g_y_array(10)
      real(dp) :: g_data(6,10)
      real(dp) :: wk_array(3,6,10)
      integer :: initiate_load, data_loaded
    
      ! these routines are called by the standard run_star check_model
      contains



      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)


         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).
         s% other_neu => other_neu


         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls
      
      subroutine f_function(y0, y1, log_y0, log_y1, g_data, g_x_list, g_y_list, wk_array, fval, ierr)
          use interp_2d_lib_db, only: interp_rgbi3p_db
          
          real(dp), intent(in) :: y0, y1, g_data(6,10), g_x_list(6), g_y_list(10), &
                                  log_y0(1), log_y1(1)
          real(dp), intent(inout) :: fval, wk_array(3, 6, 10)
          integer, intent(inout) :: ierr
          real(dp) :: g_value(1)
          
          call interp_rgbi3p_db(2, 6, 10, g_x_list, g_y_list, g_data, &
        1, log_y0, log_y1, g_value, ierr, wk_array)
        
        fval = 100*(1+y0**2)*g_value(1)/((1+y1**2)*(1+exp(y0)))
        
        !write(*,*) log_y0, log_y1, g_value, fval
        
      end subroutine
      
      
      !Friedland-Giannotti-Wise
      subroutine other_neu(  &
            id, k, T, log10_T, Rho, log10_Rho, abar, zbar, z2bar, log10_Tlim, flags, &
            loss, sources, ierr)
         use neu_lib, only: neu_get
         use neu_def
         use interp_2d_lib_db, only: interp_rgbi3p_db
         
         integer, intent(in) :: id ! id for star         
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell         
         real(dp), intent(in) :: T ! temperature
         real(dp), intent(in) :: log10_T ! log10 of temperature
         real(dp), intent(in) :: Rho ! density
         real(dp), intent(in) :: log10_Rho ! log10 of density
         real(dp), intent(in) :: abar ! mean atomic weight
         real(dp), intent(in) :: zbar ! mean charge
         real(dp), intent(in) :: z2bar ! mean charge squared
         real(dp), intent(in) :: log10_Tlim 
         logical, intent(inout) :: flags(num_neu_types) ! true if should include the type of loss
         real(dp), intent(inout) :: loss(num_neu_rvs) ! total from all sources
         real(dp), intent(inout) :: sources(num_neu_types, num_neu_rvs)
         integer, intent(out) :: ierr
            
         ! Stellar variables
         real(dp) :: ye ! Mean number of free electrons per nucleon
         real(dp) :: w ! Parameter for smoothly interpolating between deg and non-deg energy-loss
         real(dp) :: R_deg ! degeneracy factor - what's this actually called?
         real(dp) :: R_deg_max
         real(dp) :: zeta ! degeneracy parameter
         
         ! y factors
         real(dp) :: y_0, log_y_0(1) ! Temperature normalised plasma frequency
         real(dp) :: y_1, log_y_1(1) ! Temperature normalised Debye-Huckel wavenumber
         real(dp) :: y_ions, log_y_ions(1) ! y_1, but only including ions in the sum
         real(dp) :: y_TF, log_y_TF(1) ! kappa_TF/T - Temperature normalised Thomas-Fermi wavenumber
            
         ! Emission rates (erg/g/s)
         real(dp) :: f_nd ! f function for non-degenerate case
         real(dp) :: eps_nd ! non-degenerate Primakoff energy-loss rate
         real(dp) :: f_ions ! f function for case of ions only
         real(dp) :: eps_ions ! Primakoff energy-loss rates to ions only
         real(dp) :: f_el ! f function for (degenerate) electrons
         real(dp) :: eps_el ! Energy-loss rate to (degenerate) electrons
         real(dp) :: eps_deg ! Degenerate energy-loss rate (eps_ions + eps_el)
         real(dp) :: eps_total ! Total energy-loss rate, defined below
         
         ! Others
         integer :: i
         real(dp) :: y0_test(1), y1_test(1), ans(1), test_T, test_Rho
         
         type (star_info), pointer :: s
         
         include 'formats'

         ierr = 0         
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         call neu_get(  &
            T, log10_T, Rho, log10_Rho, abar, zbar, z2bar, log10_Tlim, flags, &
            loss, sources, ierr)
         if (ierr /= 0) return

         !..   Add axion production and energy loss through the  Primakoff process.
         
         if (initiate_load /= 1) then
             initiate_load = 1
         
             axion_g10 = s% job% extras_rpar(4)
             
             ! Load x and y arrays
             do i = 1,6
                 g_x_array(i) = -1.0 + (i-1)*0.5
             end do
             
             do i = 1,10
                 g_y_array(i) = -2.0 + (i-1)*0.5
             end do
             
             ! Load table from file
             open(1, file = "/data/gpfs/projects/punim0011/fhiskens/HB_stars/axions/bases/g_fun_data.txt", status = 'old', &
                  access = 'sequential', form = 'formatted', action = 'read')
             read(1,*) g_data
             close(1)
             
             ! Set up interpolation
             y0_test = 0.123
             y1_test = -0.43
             
             call interp_rgbi3p_db(1, 6, 10, g_x_array, g_y_array, g_data, & 
                             1, y0_test, y1_test, ans, ierr, wk_array)
             
             data_loaded = 1
         end if
         

         ! First calculate the degeneracy parameters
         ye = zbar/abar
         zeta = (3.01d5)*(Rho*ye)**(2./3.)/(T)
         w = atanpi(zeta-3)+0.5
         R_deg_max = max(1.50, zeta)
         R_deg = 1.50/R_deg_max


         ! Second: Calculate y factors
         y_0 = (3.33d5/T)*sqrt(Rho*ye)/(1+(1.019d-6*Rho*ye)**(2./3.))**(1./4.)
         log_y_0 = log10(y_0)
         y_1 = sqrt(6.632d20*((z2bar/abar)+ye)*Rho/(T*T*T))
         log_y_1 = log10(y_1)
         y_ions = sqrt(6.632d20*(z2bar/abar)*Rho/(T*T*T))
         log_y_ions = log10(y_ions)
         y_TF = (5.74d7/T)*(Rho*ye)**(0.16666666667)
         log_y_TF = log10(y_TF)
         
         ! Energy-loss rates - f functions
         if (data_loaded == 1) then
             call f_function(y_0, y_1, log_y_0, log_y_1, g_data, g_x_array, g_y_array, wk_array, f_nd, ierr)
             call f_function(y_0, y_ions, log_y_0, log_y_ions, g_data, g_x_array, g_y_array, wk_array, f_ions, ierr)
             call f_function(y_0, y_TF, log_y_0, log_y_TF, g_data, g_x_array, g_y_array, wk_array, f_el, ierr)
         else
             f_nd = 0.
             f_ions = 0.
             f_el = 0.
         end if
         
         ! Energy-loss rates
         eps_nd = (7.1d-52)*axion_g10**2*T**7*y_1**2*f_nd/Rho
         eps_ions = (7.1d-52)*axion_g10**2*T**7*y_ions**2*f_ions/Rho
         eps_el = 4.7d-31*axion_g10**2*R_deg*T**4*f_el*ye
         
         eps_deg = eps_ions + eps_el
         
         !Total energy-loss
         eps_total = (1-w)*eps_nd + w*eps_deg
         
         loss(ineu) = loss(ineu) + eps_total
         !loss(idneu_dT) = loss(idneu_dT) + d_sprimakoff_dT
         !loss(idneu_dRho) = loss(idneu_dRho) + d_sprimakoff_dRho
         
         
         
         !if (k == s% nz) write(*,3) 'axioncsi', s% model_number, k, axioncsi
         !if (k == s% nz) write(*,3) 'faxioncsi', s% model_number, k, faxioncsi
         !if (k == s% nz) write(*,3) 'sprimakoff', s% model_number, k, sprimakoff
         
         
            
            
      end subroutine other_neu


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
	 character(len=256) :: photosphere_summary, tau100_summary


         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

! set the correct summary file for the BC tables depending on [a/Fe]
        !photosphere_summary = 'table_' // trim(s% job% extras_cpar(2)) // 'summary.txt'
        !tau100_summary = 'table100_' // trim(s% job% extras_cpar(2)) // 'summary.txt'
        !call table_atm_init(.true., tau100_summary, photosphere_summary, ierr)





      end subroutine extras_startup


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if


         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 7
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n), c, hbar1, me, mH, wpl(50000), alpha, ne
         integer, intent(out) :: ierr
         integer :: i, nz, k(1), j
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         c = 3.d10
         me = 0.511d6
         hbar1 = 6.582d-16
         mH = 1.66054d-24
         alpha = 1.0/137.0

         nz = s% nz

         do i=1, nz
                  ne = (s% rho(i))*(s% ye(i))/mH
                  wpl(i) = sqrt(hbar1**3*c**3*4*3.14159265359*alpha*ne/me)
         end do

         names(1) = 'max_he_ye'

         if (s%   max_eps_he_k > 0) then
                vals(1) = s% ye(s% max_eps_he_k)
         else
                vals(1) = 0
         end if


         names(2) = 'central_wpl'
         vals(2) = wpl(nz)

         k = minloc(abs(wpl-10**(s% job% extras_rpar(4))))
         j = k(1)

         names(3) = 'nearest_wpl'
         vals(3) = wpl(j)

         names(4) = 'nearest_wpl_m'
         vals(4) = s% m(j)

         names(5) = 'nearest_wpl_T'
         vals(5) = s% T(j)

         names(6) = 'nearest_wpl_rho'
         vals(6) = s% rho(j)

         names(7) = 'nearest_wpl_z2bar'
         vals(7) = s% z2bar(j)

         ! note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in history column options.
         ! it must not include the new column names you are adding here.


      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 2
         
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n), c, me, hbar1, mH, alpha, wpl(50000), kb, const1, gam(50000), ne
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, i
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         
         c = 3.d10
         me = 0.511d6
         hbar1 = 6.582d-16
         mH = 1.66054d-24
         alpha = 1.0/137.0
         kb = 8.617d-5
         const1 = hbar1*hbar1*hbar1*c*c*c
         
         
         
         
         do i=1, nz
             ne = (s% rho(i))*(s% ye(i))/mH
             wpl(i) = sqrt(hbar1**3*c**3*4*3.14159265359*alpha*ne/me)
             gam(i) = const1*(4*3.14159265359*alpha**2/(3*me))*sqrt((2*3.14159265359*me)/(3*kb*(s% T(i))))* &
             (wpl(i))**2*(s% z2bar(i))*(s% Rho(i))/mH
         end do
         
         
         

         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do
         
         names(1) = 'plasma_frequency'
         vals(:, 1) = wpl(1:nz)
         
         names(2) = 'gamma_1'
         vals(:,2) = gam(1:nz)
         
         

      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
	 real(dp) :: envelope_mass_fraction, L_He, L_tot, min_center_h1_for_diff, &
            critmass, feh, rot_full_off, rot_full_on, frac2, mass_difference
         real(dp), parameter :: huge_dt_limit = 3.15d16 ! ~1 Gyr
         real(dp), parameter :: new_varcontrol_target = 1d-3
         real(dp), parameter :: Zsol = 0.0142
         type (star_info), pointer :: s
	    logical :: diff_test1, diff_test2, diff_test3


         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going


	 ! set BC: change to tables after running on simple photosphere
        !if (s% model_number == 100) then
        !   write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        !   write(*,*) 'switching from simple photosphere to ', s% job% extras_cpar(1)
        !   write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
        !   s% atm_option = s% job% extras_cpar(1)
        !endif


! set ROTATION: extra param are set in inlist: star_job
        rot_full_off = s% job% extras_rpar(1) !1.2
        rot_full_on = s% job% extras_rpar(2) !1.8

        if (rot_set_check == 0) then
            if ((s% job% extras_rpar(3) > 0.0) .and. (s% initial_mass > rot_full_off)) then
                !check if ZAMS is achieved, then set rotation
                if ((abs(log10(s% power_h_burn * Lsun / s% L(1))) < 0.01) .and. (s% star_age > 1d2)) then
                    if (s% initial_mass <= rot_full_on) then
                        frac2 = (s% initial_mass - rot_full_off) / (rot_full_on - rot_full_off)
                        frac2 = 0.5d0*(1 - cos(pi*frac2))
                    else
                        frac2 = 1.0
                    end if
                    write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                    write(*,*) 'new omega_div_omega_crit, fraction', s% job% extras_rpar(3) * frac2, frac2
                    write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                    s% job% new_omega_div_omega_crit = s% job% extras_rpar(3) * frac2
                    s% job% set_near_zams_omega_div_omega_crit_steps = 10
                    rot_set_check = 1
                end if
            else
                write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                write(*,*) 'no rotation'
                write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                rot_set_check = 1
            end if

        end if

! End run at TP AGB
        mass_difference = s% he_core_mass - s% c_core_mass
        !write(*,*) s% phase_of_evolution
        
        if (s% center_he4 < 1d-4) then
            if (s% phase_of_evolution == 8) then
                write(*,*) 'Star at TP-AGB phase - stopping run'
                call star_write_model(id, s% job% save_model_filename, ierr)
                extras_finish_step = terminate
            !else if (mass_difference < 0.5*(s% he_core_mass)) then
            !    if  (s% luminosity_by_category(2, 1) > s% luminosity_by_category(3, 1)) then
            !        write(*,*) 'Star at TP-AGB phase - stopping run'
            !        call star_write_model(id, s% job% save_model_filename, ierr)
            !        extras_finish_step = terminate
            !    end if
            end if
        end if

! check DIFFUSION: to determine whether or not diffusion should happen
! no diffusion for fully convective, post-MS, and mega-old models
! do diffusion during the WD phase
	    min_center_h1_for_diff = 1d-10
	    diff_test1 = abs(s% mass_conv_core - s% star_mass) < 1d-2 !fully convective
	    diff_test2 = s% star_age > 5d10 !really old
	    diff_test3 = s% center_h1 < min_center_h1_for_diff !past the main sequence
	    if( diff_test1 .or. diff_test2 .or. diff_test3 )then
            s% diffusion_dt_limit = huge_dt_limit
        else
            s% diffusion_dt_limit = original_diffusion_dt_limit
	    end if

        if (wd_diffusion) then
            s% diffusion_dt_limit = original_diffusion_dt_limit
        end if



         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve



      end module run_star_extras
