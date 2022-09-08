! ***********************************************************************
!
!   Copyright (C) 2012-2019  Bill Paxton & The MESA Team
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
!   gnu library general public license for more details.o;
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! *********************************************************************** 
      module run_binary_extras 

      use star_lib
      use star_def
      use const_def
      use const_def
      use chem_def
      use num_lib
      use binary_def
      use math_lib
      
      implicit none
      
      contains
      
     
      
      subroutine extras_binary_controls(binary_id, ierr)
         integer :: binary_id
         integer, intent(out) :: ierr
         type (binary_info), pointer :: b
         ierr = 0

         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         write(*,*) 'hello from extra_binary_controls'

         ! Set these function pointers to point to the functions you wish to use in
         ! your run_binary_extras. Any which are not set, default to a null_ version
         ! which does nothing.
         b% how_many_extra_binary_history_header_items => how_many_extra_binary_history_header_items
         b% data_for_extra_binary_history_header_items => data_for_extra_binary_history_header_items
         b% how_many_extra_binary_history_columns => how_many_extra_binary_history_columns
         b% data_for_extra_binary_history_columns => data_for_extra_binary_history_columns

         b% other_mdot_edd => eval_mdot_edd_routine
         !b% other_jdot_ml  => default_jdot_ml  
         b% other_jdot_ml  => rcl_jdot_ml        ! take the SAM of CL when radio turn on (disk instable) Jia
         !b% other_jdot_mb  => cz_jdot_mb          ! considering convective envelop 
         !b% other_mdot_edd => my_mdot_edd

         b% extras_binary_startup=> extras_binary_startup
         b% extras_binary_start_step=> extras_binary_start_step
         b% extras_binary_check_model=> extras_binary_check_model
         b% extras_binary_finish_step => extras_binary_finish_step
         b% extras_binary_after_evolve=> extras_binary_after_evolve

         ! Once you have set the function pointers you want, then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
         ! b% warn_binary_extra =.false.
         
      end subroutine extras_binary_controls

!      subroutine my_mdot_edd(binary_id, mdot_edd, mdot_edd_eta, ierr)
!         use const_def, only: dp
!         integer, intent(in) :: binary_id
!         real(dp), intent(out) :: mdot_edd
!         real(dp), intent(out) :: mdot_edd_eta
!         integer, intent(out) :: ierr
!         type (binary_info), pointer :: b
!         ierr = 0
!         call binary_ptr(binary_id, b, ierr)
!         if (ierr /= 0) then
!            write(*,*) 'failed in binary_ptr'
!            return
!         end if
!         mdot_edd = 0d0 
!         mdot_edd_eta = 0d0 
!      end subroutine my_mdot_edd
         

      subroutine eval_mdot_edd_routine(binary_id , mdot_edd , mdot_edd_eta, ierr)
         integer , intent (in) :: binary_id
         real ( dp ) , intent (out) :: mdot_edd
         real(dp), intent(out) :: mdot_edd_eta
         integer , intent ( out ) :: ierr
         real :: mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
                 cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
                 ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
                 rlo_start_age,rlo_age,radio_start_age,radio_age,vesc_sur_d,eva_wind,&
                 mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,hoc_accretor,hoc_donor
         !lp_chen: Chen,H-L et al.,2013 ApJ

         common /binary/   mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
                           cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
                           ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
                           rlo_start_age,rlo_age,radio_start_age,radio_age,vesc_sur_d,eva_wind,&
                           mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,&
                           hoc_accretor,hoc_donor
         type(binary_info), pointer :: b
         include 'formats.inc'
         ierr = 0
         
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
                 
         mdot_edd = 3.6d-8*(b% m(2)/(1.4*Msun))*(1.0d5*(clight**2.0)/(6.67428d-8*b% m(2)))*&
                   (1.7/(1.0+b% s1% X(1)))*Msun/secyer
!WHY The following line is there?????????????????????????????????
         if (mdot_act .lt. mdot_critical_dubus) then
            mdot_edd = 0.01 * mdot_edd
         end if
         mdot_edd_eta = 0.5
      end subroutine eval_mdot_edd_routine

      integer function how_many_extra_binary_history_header_items(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_header_items = 0
      end function how_many_extra_binary_history_header_items


      subroutine data_for_extra_binary_history_header_items( &
           binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id, n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
      end subroutine data_for_extra_binary_history_header_items


      integer function how_many_extra_binary_history_columns(binary_id)
         use binary_def, only: binary_info
         integer, intent(in) :: binary_id
         how_many_extra_binary_history_columns = 0
      end function how_many_extra_binary_history_columns


      subroutine data_for_extra_binary_history_columns(binary_id, n, names, vals, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(in) :: n
         character (len=maxlen_binary_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if
         
      end subroutine data_for_extra_binary_history_columns
      
      
      integer function extras_binary_startup(binary_id,restart,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         logical, intent(in) :: restart
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
         
!          b% s1% job% warn_run_star_extras = .false.
          extras_binary_startup = keep_going
      end function  extras_binary_startup
      
      integer function extras_binary_start_step(binary_id,ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr

         extras_binary_start_step = keep_going
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if
      
      end function  extras_binary_start_step
      
      !Return either keep_going, retry or terminate
      integer function extras_binary_check_model(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_check_model = keep_going
        
      end function extras_binary_check_model
      
      
      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_binary_finish_step(binary_id)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if  
         extras_binary_finish_step = keep_going
         
      end function extras_binary_finish_step
      
      subroutine extras_binary_after_evolve(binary_id, ierr)
         type (binary_info), pointer :: b
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then ! failure in  binary_ptr
            return
         end if      
         
 
      end subroutine extras_binary_after_evolve   
      
      
      subroutine rcl_jdot_ml(binary_id, ierr)
         !rlo_start_age: mdot_act > 1.0e-50
         !rlo_age : age - rlo_start_age
         !mdot_act : rlo + accreted wind
         !radio_start_age : when accreted mass > 0.05 (Antoniadis 2012) and rlo > mdot_critical (Dubus 1999)
         !radio_age : age - radio_start_age
         !
         INTEGER*4 getcwd, status
         character(LEN=100) dirname,dirname_postfix,dirname_prefix
         character n_dirname
         integer n_string
         integer, intent(in) :: binary_id
         integer, intent(out) :: ierr
         integer model_number_last,i,store_read
         integer recyelled
         real :: L1_to_acc, L1_to_cm, w_Kopal


         real :: mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
                 cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
                 ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
                 rlo_start_age,rlo_age,radio_start_age,radio_age,vesc_sur_d,eva_wind,&
                 mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,hoc_accretor,hoc_donor
         !lp_chen: Chen,H-L et al.,2013 ApJ
         common /dirname/ dirname
         common /binary/   mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
                           cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
                           ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
                           rlo_start_age,rlo_age,radio_start_age,radio_age,vesc_sur_d,eva_wind,&
                           mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,&
                           hoc_accretor,hoc_donor,model_number_last,recyelled

         type (binary_info), pointer :: b
         ierr = 0
         call binary_ptr(binary_id, b, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in binary_ptr'
            return
         end if

         if (b% s_donor% model_number .eq. 1)   then
		rlo_start_age = 0.0
		rlo_age = 0.0
		radio_start_age = 0.0
		radio_age = 0.0
		maccretor_i = (b% m(2)/Msun)
		eva_wind = 0.0
            mtransfer_rate_dot = 0.0
            mtransfer_rate_dot_div_mtransfer = 0.0
            recyelled = 0
         end if

         accreted_mass = (b% m(2)/Msun) - maccretor_i
         if (accreted_mass .ge. 0.1) recyelled = 1
	   mdot_critical_dubus = 3.2D-9*(((b% m(2)/Msun)/1.4)**0.5)*&
			(((b% m(1)/Msun)/1.0)**(-0.2))*((b% period/86400)**1.4)    !disk instability (Dubus et al. 1999)
         mdot_act = abs((b% mtransfer_rate) + (b% mdot_wind_transfer(b% d_i)))*secyer/Msun              
         ! rlo + wind_acc
         rlo = abs(b% mtransfer_rate)*secyer/Msun
         mdonor_wind = (b% mdot_system_transfer(b% d_i) + b% mdot_system_wind(b% d_i))*secyer/Msun
         maccretor_wind = (b% mdot_system_transfer(b% a_i) + b% mdot_system_wind(b% a_i))*secyer/Msun
         cbd_wind = b% mdot_system_cct*secyer/Msun
         mdot_edd_here = 3.6d-8*(b% m(2)/(1.4*Msun))*(1.0d5*(clight**2.0)/(6.67428d-8*b% m(2)))*&
                   (1.7/(1.0+b% s1% X(1)))
         vesc_sur_d = (2*6.672D-8*b% m(1)/b% r(1))**0.5

         if (rlo .ge. mdot_critical_dubus) radio_start_age = 0.0
         if ((rlo .gt. 1.0e-50) .and. (rlo_start_age .eq. 0.0)) rlo_start_age = b% s_donor% time/secyer
         if (rlo_start_age .gt. 0.0) rlo_age = b% s_donor% time/secyer - rlo_start_age
         if (radio_start_age .gt. 0.0) then
            radio_age = ( b% s_donor% time/secyer - radio_start_age)
            eva_wind = (b% s_donor% x_ctrl(1)*(1.0/(2.0*(vesc_sur_d**2)))*lp_chen*(b% r(1)/b% separation)**2)*secyer/Msun
         else
            radio_age = 0.0
            eva_wind = 0.0
         end if


         if (((radio_start_age .eq. 0) .and. (recyelled .eq. 1)   &
         .and. (rlo .le. mdot_critical_dubus)) ) then
         !(mtransfer_rate_dot_div_mtransfer .LE. 1E-7 .or. rlo .eq. 0)
            radio_start_age = b% s_donor% time/secyer

         end if



         !!!!!!!!!!!!!!!!!!!!!!!!
         ps_chen = 2.0*pi/(8.11156D11/(1.5D17+(radio_age)*3.15E7)**0.5)    !Chen et. al 2013 APJ 775:27 P0=3ms B~10^8G
         psdot_chen = (6.9813D-15/(2.0*3.14))*ps_chen**2.0
         lp_chen = 4.0*(pi**2)*(1.0D45)*psdot_chen/(ps_chen**3.0)
         rcl_chen = 2.99792458D10*ps_chen/(2.0*pi)

         ! wind in units of Msun/year  (value is >= 0)

         if (radio_age .eq. 0 .or. b% s_donor% x_ctrl(1) .eq. 0) then !no evaporation : disk stable + accreted_mass > 0.05

         	!mass lost from vicinity of donor
         	b% jdot_ml = (b% mdot_system_transfer(b% d_i) + b% mdot_system_wind(b% d_i))*&
             		(b% m(b% a_i)/(b% m(b% a_i)+b% m(b% d_i))*b% separation)**2*2*pi/b% period
         	!mass lost from vicinity of accretor
         	b% jdot_ml = b% jdot_ml + (b% mdot_system_transfer(b% a_i) + b% mdot_system_wind(b% a_i))*&
             		(b% m(b% d_i)/(b% m(b% a_i)+b% m(b% d_i))*b% separation)**2*2*pi/b% period
         	!mass lost from circumbinary coplanar toroid
         	b% jdot_ml = b% jdot_ml + b% mdot_system_cct * b% mass_transfer_gamma * &
             		sqrt(b% s_donor% cgrav(1) * (b% m(1) + b% m(2)) * b% separation)

         else   !evaporation start : no accretion b% m(2) = conste, take CL SAM
                b% m(2) = b% m_old(2)
         	                                       !mass lost from vicinity of donor
         	b% jdot_ml = (b% mdot_system_transfer(b% d_i) + b% mdot_system_wind(b% d_i))*&
             		(b% m(b% a_i)/(b% m(b% a_i)+b% m(b% d_i))*b% separation)**2*2*pi/b% period

            if (b% m(b% a_i)/b% m(b% d_i) .le. 10.0) then
            L1_to_acc = 0.5 + 0.227*log10(b% m(b% a_i)/b% m(b% d_i))               !distance from L1 to accretor in unit of separation
            else
            w_Kopal = (1.0/(3.0*(1+(b% m(b% a_i)/b% m(b% d_i)))))**(1.0/3.0)       !Beer 2007;Kopal1959
            L1_to_acc = 1.0 - w_Kopal + ((w_Kopal**2.0)/3.0) + ((w_Kopal**3.0)/9.0)           !distance from L1 to accretor in unit of separation
            endif
            L1_to_cm = abs(b% m(b% d_i)/(b% m(b% a_i)+b% m(b% d_i)) - L1_to_acc)   !distance from L1 to center of mass in unit of separation

         	!mass lost from L1 point
         	!b% jdot_ml = (b% mdot_system_transfer(b% d_i) + b% mdot_system_wind(b% d_i))*&
            !		(L1_to_cm*b% separation)**2*2*pi/b% period


         	!evaporation wind takes no AML
         	!b% jdot_ml = 0.0
         	!mass lost from vicinity of CL of accretor
         	b% jdot_ml = b% jdot_ml + (b% mdot_system_transfer(b% a_i) + b% mdot_system_wind(b% a_i))*&
             		((b% m(b% d_i)/(b% m(b% a_i)+b% m(b% d_i))*b% separation) - rcl_chen)**2*2*pi/b% period
         	!mass lost from circumbinary coplanar toroid
         	b% jdot_ml = b% jdot_ml + b% mdot_system_cct * b% mass_transfer_gamma * &
             		sqrt(b% s_donor% cgrav(1) * (b% m(1) + b% m(2)) * b% separation)
         end if

         mdonor_dot = ((b% m(1) - b% m_old(1))/(b% s_donor% time_step))/Msun
         maccretor_dot = ((b% m(2) - b% m_old(2))/(b% s_donor% time_step))/Msun
         mtransfer_rate_dot = abs(b% mtransfer_rate - b% mtransfer_rate_old)/b% s_donor% time_step

         if (b% mtransfer_rate .ne. 0) then
         	mtransfer_rate_dot_div_mtransfer = mtransfer_rate_dot/abs(b% mtransfer_rate)
         end if

         !if (b% s_donor% model_number .eq. 1) then
         !open (unit=10, file='fort.10')	
         !write (10,"(A7,100(A20))") 'No.','dt(yr)','age(yr)','Md(Msun)','Porb(d)','Sep','omega',&
         ! 	'Ma(Msun)','Md_dot','Ma_dot','Mdot_edd','Md_wind','Ma_wind','CBD_wind','Rlo(Msun/yr)',&
         ! 	'rlo_fraction','wind_fra(d_to_a)','wind_fra(a_to_d)','R_donor','Rl_donor',&
         ! 	'Ts_donor','L_donor','logP_sur_d','logrho_sur_d','logg_sur_d','vesc_sur_d','eva_wind',&
         !	'Jdot','Jdot_mb','Jdot_gr','Jdot_ml', 'jdot_ls', 'jdot_missing_wind', 'extra_jdot','orbital_am','rlo_dot','rlo_dot/Mdot'
         !close (10)
         !end if

         !open (unit=10, file='fort.10', access='append')	
         !write (10,"(I7,100(E20.4E3))") b% s_donor% model_number,&
	   !		b% s_donor% time_step/secyer,b% s_donor% time/secyer,b% m(1)/Msun,&
		!	b% period/86400,b% separation,b% s_donor% omega(1), b% m(2)/Msun,&
		!	mdonor_dot,maccretor_dot,mdot_edd_here,mdonor_wind,maccretor_wind,&
		!	cbd_wind,rlo,b% xfer_fraction,&
		!	b% wind_xfer_fraction(b% d_i),b% wind_xfer_fraction(b% a_i),&
		!	b% r(1),b% rl(1),b% s_donor% Teff,b% s_donor% L_phot,&
		!	b% s_donor% log_surface_pressure,b% s_donor% log_surface_density,&
		!	b% s_donor% log_surface_gravity,vesc_sur_d,eva_wind,&
		!	b% jdot, b% jdot_mb, b% jdot_gr, b% jdot_ml, b% jdot_ls, &
		!	b% jdot_missing_wind, b% extra_jdot, b% angular_momentum_j,&
		!	mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer
			
         !close (10)


         if (b% s_donor% model_number .eq. 1) then
         	model_number_last = 0
            status = getcwd(dirname)
         	do  n_string=1,len(dirname)-3
                  n_dirname = dirname(len(dirname)-n_string-2:len(dirname)-n_string-1)
         		if (n_dirname .eq. '/' ) dirname=dirname(len(dirname)-1-n_string:len(dirname))
            end do
            dirname_postfix = '.59'
            dirname = trim(dirname)//trim(dirname_postfix)
         	open (unit=59, file=dirname)
         	close (59,status='delete')
         end if





         if (b% s_donor% model_number .gt. model_number_last) then

         open (unit=59, file=dirname, access='append')	
         write (59,"(I7,2I5,100(E20.4E3))") &
            b% s_donor% model_number, b% d_i, b% accretion_mode,           &
            b% s_donor% star_age,                                          &
            radio_start_age, b% m(1)/Msun, b% m(2)/Msun, mdonor_dot,       &
            maccretor_dot, mdot_edd_here, rlo, b% wind_xfer_fraction,           &
            mdonor_wind, eva_wind, maccretor_wind, mdot_act,               &
            mdot_critical_dubus, b% period/86400, b% separation/Rsun,      &
            radio_age, ps_chen, lp_chen, b% eccentricity, b% r(1)/Rsun,    &
            b% r(2)/Rsun, b% rl(1)/Rsun, b% rl(2)/Rsun,                    &
            b% angular_momentum_j, b% s1% total_angular_momentum,          &
            b% jdot, b% jdot_mb, b% jdot_gr, b% jdot_ml, b% jdot_ls,       &
            b% jdot_missing_wind, b% extra_jdot,                           &
            b% s_donor% L_phot, b% s_donor% Teff,                          &
            b% s_donor% log_surface_density, b% s_donor% he_core_mass,     &
            b% s_donor% c_core_mass, b% s_donor% o_core_mass,              &
            b% s_donor% si_core_mass, b% s_donor% fe_core_mass,            &
            b% s_donor% power_h_burn, b% s_donor% power_he_burn,           &
            b% s_donor% power_c_burn

         else

         open (unit=59, file=dirname, position='REWIND')	
         do i=1,model_number_last-1
         read (59,*) store_read
         end do

         write (59,"(I7,2I5,100(E20.4E3))") &
            b% s_donor% model_number, b% d_i, b% accretion_mode,           &
            b% s_donor% star_age,                                          &
            radio_start_age, b% m(1)/Msun, b% m(2)/Msun, mdonor_dot,       &
            maccretor_dot, mdot_edd_here, rlo, b% wind_xfer_fraction,           &
            mdonor_wind, eva_wind, maccretor_wind, mdot_act,               &
            mdot_critical_dubus, b% period/86400, b% separation/Rsun,      &
            radio_age, ps_chen, lp_chen, b% eccentricity, b% r(1)/Rsun,    &
            b% r(2)/Rsun, b% rl(1)/Rsun, b% rl(2)/Rsun,                    &
            b% angular_momentum_j, b% s1% total_angular_momentum,          &
            b% jdot, b% jdot_mb, b% jdot_gr, b% jdot_ml, b% jdot_ls,       &
            b% jdot_missing_wind, b% extra_jdot,                           &
            b% s_donor% L_phot, b% s_donor% Teff,                          &
            b% s_donor% log_surface_density, b% s_donor% he_core_mass,     &
            b% s_donor% c_core_mass, b% s_donor% o_core_mass,              &
            b% s_donor% si_core_mass, b% s_donor% fe_core_mass,            &
            b% s_donor% power_h_burn, b% s_donor% power_he_burn,           &
            b% s_donor% power_c_burn

         end if

         model_number_last = b% s_donor% model_number

         close (59)

      end subroutine rcl_jdot_ml 
      
!       subroutine cz_jdot_mb(binary_id, ierr)      !consider convective envelop hoc factor
!         use mlt_def
!         integer model_number_last
!         real :: mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
!                 cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
!                 ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
!                 rlo_start_age,rlo_age,radio_start_age,radio_age,vesc_sur_d,eva_wind,&
!                 mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,hoc_accretor,hoc_donor
!         !lp_chen: Chen,H-L et al.,2013 ApJ
!         integer recyelled
!         common /binary/   mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
!                           cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
!                           ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
!                           rlo_start_age,rlo_age,radio_start_age,radio_age,vesc_sur_d,eva_wind,&
!                           mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,&
!                           hoc_accretor,hoc_donor,model_number_last,recyelled
!
!         integer, intent(in) :: binary_id
!         integer, intent(out) :: ierr
!         integer :: i, k, id
!         type (binary_info), pointer :: b
!         type (star_info), pointer :: s
!         real(dp) :: rsun4,two_pi_div_p3
!         ierr = 0
!         call binary_ptr(binary_id, b, ierr)
!         if (ierr /= 0) then
!            write(*,*) 'failed in binary_ptr'
!            return
!         end if
!         b% jdot_mb = 0
!         rsun4 = rsun*rsun*rsun*rsun
!         call check_radiative_core(b)
!         call check_convective_envelop(b)     !hoc_factor
!
!         two_pi_div_p3 = (2.0*pi/b% period)*(2.0*pi/b% period)*(2.0*pi/b% period)
!         ! use the formula from rappaport, verbunt, and joss.  apj, 275, 713-731. 1983.
!         if (b% have_radiative_core(b% d_i) .or. b% keep_mb_on) &
!            b% jdot_mb = -3.8d-30*b% m(b% d_i)*rsun4* &         
!                           pow_cr(min(b% r(b% d_i),b% rl(b% d_i))/rsun,b% magnetic_braking_gamma)* &
!                           two_pi_div_p3
!            b% jdot_mb = b% jdot_mb * hoc_donor
!         if (b% evolve_both_stars .and. b% include_accretor_mb .and. &
!             (b% have_radiative_core(b% a_i) .or. b% keep_mb_on)) then
!             b% jdot_mb = b% jdot_mb - &
!                           3.8d-30*b% m(b% a_i)*rsun4* &
!                           pow_cr(min(b% r(b% a_i),b% rl(b% a_i))/rsun,b% magnetic_braking_gamma)* &
!                           two_pi_div_p3
!             b% jdot_mb = b% jdot_mb * hoc_accretor
!         end if
!
!
!      end subroutine cz_jdot_mb
!
!      subroutine cz_jdot_mb(binary_id, ierr)      !consider convective envelop hoc factor
!         use mlt_def
!         integer model_number_last
!         real :: mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
!                 cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
!                 ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
!                 rlo_start_age,rlo_age,radio_start_age,radio_age,vesc_sur_d,eva_wind,&
!                 mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,hoc_accretor,hoc_donor
!         !lp_chen: Chen,H-L et al.,2013 ApJ
!         integer recyelled
!         common /binary/   mdonor_dot,maccretor_dot,mdonor_wind,maccretor_wind,&
!                           cbd_wind,rlo,mdot_act,mdot_critical_dubus,mdot_edd_here,&
!                           ps_chen,psdot_chen,lp_chen,rcl_chen,maccretor_i,accreted_mass,&
!                           mtransfer_rate_dot,mtransfer_rate_dot_div_mtransfer,&
!                           hoc_accretor,hoc_donor,model_number_last,recyelled
!
!         integer, intent(in) :: binary_id
!         integer, intent(out) :: ierr
!         integer :: i, k, id
!         type (binary_info), pointer :: b
!         type (star_info), pointer :: s
!         real(dp) :: rsun4,two_pi_div_p3
!         ierr = 0
!         call binary_ptr(binary_id, b, ierr)
!         if (ierr /= 0) then
!            write(*,*) 'failed in binary_ptr'
!            return
!         end if
!         b% jdot_mb = 0
!         rsun4 = rsun*rsun*rsun*rsun
!         call check_radiative_core(b)
!         call check_convective_envelop(b)     !hoc_factor
!
!         two_pi_div_p3 = (2.0*pi/b% period)*(2.0*pi/b% period)*(2.0*pi/b% period)
!         ! use the formula from rappaport, verbunt, and joss.  apj, 275, 713-731. 1983.
!         if (b% have_radiative_core(b% d_i) .or. b% keep_mb_on) &
!            b% jdot_mb = -3.8d-30*b% m(b% d_i)*rsun4* &         
!                           pow_cr(min(b% r(b% d_i),b% rl(b% d_i))/rsun,b% magnetic_braking_gamma)* &
!                           two_pi_div_p3
!            b% jdot_mb = b% jdot_mb * hoc_donor
!         if (b% evolve_both_stars .and. b% include_accretor_mb .and. &
!             (b% have_radiative_core(b% a_i) .or. b% keep_mb_on)) then
!             b% jdot_mb = b% jdot_mb - &
!                           3.8d-30*b% m(b% a_i)*rsun4* &
!                           pow_cr(min(b% r(b% a_i),b% rl(b% a_i))/rsun,b% magnetic_braking_gamma)* &
!                           two_pi_div_p3
!             b% jdot_mb = b% jdot_mb * hoc_accretor
!         end if
!
!
!      end subroutine cz_jdot_mb
! 
      
      end module run_binary_extras
