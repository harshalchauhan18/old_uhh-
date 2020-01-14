#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_uhh_clc - CLC cyanobacteria lifestage
!
! !INTERFACE:
   module fabm_uhh_clc
!
! !DESCRIPTION:
!
! Model of a cyanobacteria lifestage within a cyanobacteria lifecycle.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_uhh_clc
!     Variable identifiers
      type (type_state_variable_id)        :: id_c,id_s,id_g
      type (type_state_variable_id)        :: id_nitrate,id_ammonium
      type (type_state_variable_id)        :: id_phosphate,id_detritus,id_oxygen
      type (type_dependency_id)            :: id_par,id_temp
      type (type_state_variable_id)        :: id_nextc,id_nexts,id_nextg
      type (type_horizontal_dependency_id) :: id_I_0

!     Model parameters
      real(rk) :: sr
      real(rk) :: cya0
      real(rk) :: alpha
      real(rk) :: kc
      real(rk) :: s2
      real(rk) :: s3
      real(rk) :: mort
      real(rk) :: w
      real(rk) :: Qc
      real(rk) :: E_sc  ! energy storage capacity Emax in Hense&Beckmann 2006
      real(rk) :: kN
      real(rk) :: tscale
      real(rk) :: trange
      real(rk) :: theta_e
      real(rk) :: theta_q
      real(rk) :: rmatscale
      real(rk) :: omega0
      real(rk) :: scale=8.0_rk
      logical  :: n_fixation, lifecycling
      real(rk) :: m
      real(rk) :: uptake_factor
      real(rk) :: growth_factor
      real(rk) :: lightcapture_factor
      real(rk) :: e_max,e_min,q_max,q_min
      real(rk) :: Sflux_per_Cflux, Gflux_per_Cflux


      contains

      procedure :: initialize
      procedure :: do
      procedure :: get_light_extinction
   end type
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the CLC model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the uhh_clc namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_uhh_clc), intent(inout), target :: self
   integer,                        intent(in)          :: configunit
!
! !LOCAL VARIABLES:
   real(rk)                  :: c_initial,g_initial,s_initial
   real(rk)                  :: background_concentration
   real(rk)                  :: w
   real(rk)                  :: kc
   real(rk)                  :: mortality_rate,excretion_rate
   real(rk)                  :: i_min
   real(rk)                  :: rmax
   real(rk)                  :: alpha
   real(rk)                  :: tbg
   real(rk)                  :: beta_bg
   real(rk)                  :: sr,s2,s3
   real(rk)                  :: Qc
   real(rk)                  :: kN
   real(rk)                  :: e_max,e_min,q_max,q_min
   real(rk)                  :: Sflux_per_Cflux, Gflux_per_Cflux
   real(rk)                  :: e_storage_capacity
   real(rk)                  :: tscale
   real(rk)                  :: trange
   real(rk)                  :: theta_e
   real(rk)                  :: theta_q
   real(rk)                  :: omega0
   real(rk)                  :: rmatscale
   real(rk)                  :: scale=8.0_rk
   logical                   :: n_fixation
   real(rk)                  :: m
   real(rk)                  :: uptake_factor
   real(rk)                  :: growth_factor
   real(rk)                  :: lightcapture_factor
   character(len=64)         :: din_variable
   character(len=64)         :: next_in_cycle
   character(len=64)         :: phosphate_variable
   character(len=64)         :: ammonium_variable
   character(len=64)         :: nitrate_variable
   character(len=64)         :: detritus_variable
   character(len=64)         :: oxygen_variable
   character(len=64)         :: lifestage_name

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   namelist /uhh_clc/ &
             c_initial,s_initial,g_initial,background_concentration, &
             w, kc, sr, s2, s3, rmatscale, omega0, alpha, &
             Qc,kN,tscale,trange,theta_e,theta_q,scale, e_storage_capacity, &
             e_max,e_min,q_max,q_min,Sflux_per_Cflux, Gflux_per_Cflux, &
             mortality_rate, next_in_cycle, n_fixation, m, &
             uptake_factor, growth_factor, lightcapture_factor, &
             ammonium_variable,nitrate_variable, &
             phosphate_variable, detritus_variable, &
             oxygen_variable, lifestage_name
!EOP
!-----------------------------------------------------------------------
!BOC
   c_initial = 0.0_rk
   g_initial = 0.0_rk
   s_initial = 0.0_rk

!   background_concentration = 0.0225_rk
   w         = 0.1_rk
   kc        = 0.03_rk
   sr        = 0.0625_rk
   s2        = 6.625_rk
   s3        = 8.625_rk 
   Qc        = 0.0_rk
   kN        = 0.0_rk
   e_max     = 1.0_rk
   e_min     = 0.0_rk
   q_max     = 1.0_rk
   q_min     = 0.0_rk
   Sflux_per_Cflux = 1.0_rk
   Gflux_per_Cflux = 1.0_rk
   e_storage_capacity = 5.5_rk
   tscale    = 0.0_rk
   trange    = 0.0_rk
   theta_e   = 0.0_rk
   theta_q   = 0.0_rk
   scale     = 8.0_rk
   rmatscale = 45.0_rk
   omega0    = 2.8e-6_rk
   n_fixation=.false.
   m         = 3.0_rk
   uptake_factor = 0.0_rk
   growth_factor = 0.0_rk
   lightcapture_factor = 0.0_rk
   mortality_rate = 0.02_rk ! mort [1/d]
   next_in_cycle=''
   nitrate_variable = 'uhh_ergom_split_base_nit'
   ammonium_variable = 'uhh_ergom_split_base_amm'
   phosphate_variable = 'uhh_ergom_split_base_pho'
   detritus_variable = 'uhh_ergom_split_base_det'
   oxygen_variable = 'uhh_ergom_split_base_oxy'
   
   ! Read the namelist
   if (configunit>=0) read(configunit,nml=uhh_clc,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%cya0, 'background_concentration', default=background_concentration)
   call self%get_parameter(self%kc,     'kc',     default=kc)
   call self%get_parameter(self%sr,     'sr',     default=sr)
   call self%get_parameter(self%s2,     's2',     default=s2)
   call self%get_parameter(self%s3,     's3',     default=s3)
   call self%get_parameter(self%e_max,  'e_max',  default=e_max)
   call self%get_parameter(self%e_min,  'e_min',  default=e_min)
   call self%get_parameter(self%q_max,  'q_max',  default=q_max)
   call self%get_parameter(self%q_min,  'q_min',  default=q_min)

   ! set default S,G flux scales depending on E,Q range:
   if (self%q_min > 0.0_rk) Sflux_per_Cflux=0.5_rk ! e.g. V -> H
   if (self%e_min > 0.0_rk) Gflux_per_Cflux=0.5_rk ! e.g. H -> A
   if (self%e_max < 1.0_rk) Gflux_per_Cflux=1.5_rk ! e.g. R -> V
   if (self%q_max < 1.0_rk) Sflux_per_Cflux=1.5_rk ! e.g. A -> R
   call self%get_parameter(self%Sflux_per_Cflux, 'Sflux_per_Cflux', default=Sflux_per_Cflux)
   call self%get_parameter(self%Gflux_per_Cflux, 'Gflux_per_Cflux', default=Gflux_per_Cflux)

   call self%get_parameter(self%e_sc, 'e_storage_capacity', default=e_storage_capacity)
   call self%get_parameter(self%tscale, 'tscale', default=tscale)
   call self%get_parameter(self%trange, 'trange', default=trange)
   call self%get_parameter(self%theta_e,'theta_e',default=theta_e)
   call self%get_parameter(self%theta_q,'theta_q',default=theta_q)
   call self%get_parameter(self%kN,     'kN',     default=kN)
   call self%get_parameter(self%Qc,     'Qc',     default=Qc)
   call self%get_parameter(self%m,      'm',      default=m)
   call self%get_parameter(self%n_fixation, 'n_fixation', default=n_fixation)
   call self%get_parameter(self%uptake_factor, 'uptake_factor', default=uptake_factor)
   call self%get_parameter(self%growth_factor, 'growth_factor', default=growth_factor)
   call self%get_parameter(self%lightcapture_factor, 'lightcapture_factor', default=lightcapture_factor)
   call self%get_parameter(self%omega0, 'omega0', default=omega0)
   call self%get_parameter(self%rmatscale, 'rmatscale', default=rmatscale)
   call self%get_parameter(self%w,      'w',      default=w,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%mort,   'mortality_rate', default=mortality_rate,  scale_factor=1.0_rk/secs_pr_day)
   
   ! Register state variables
   call self%register_state_variable(self%id_c,'C'//lifestage_name(1:3), &
         'mmol n/m**3',trim(lifestage_name)//' biomass', c_initial, &
         standard_variable=type_bulk_standard_variable( &
             name='cyanobacteria_'//trim(lifestage_name)//'_biomass'), &
         minimum=0.0_rk,vertical_movement=self%w)

   call self%register_state_variable(self%id_g,'G'//lifestage_name(1:3), &
         'mmol n/m**3',trim(lifestage_name)//' energy', g_initial, &
         standard_variable=type_bulk_standard_variable( &
             name='cyanobacteria_'//trim(lifestage_name)//'_energy'), &
         minimum=0.0_rk,vertical_movement=self%w)

   call self%register_state_variable(self%id_s,'S'//lifestage_name(1:3), &
         'mmol n/m**3',trim(lifestage_name)//' quota', s_initial, &
         standard_variable=type_bulk_standard_variable( &
             name='cyanobacteria_'//trim(lifestage_name)//'_quota'), &
         minimum=0.0_rk,vertical_movement=self%w)


   ! Register dependencies on external standard variables
   call self%register_state_dependency(self%id_ammonium, 'ammonium_target', 'mmol/m**3','ammonium source')
   call self%register_state_dependency(self%id_nitrate, 'nitrate_target', 'mmol/m**3','nitrate source')
   call self%register_state_dependency(self%id_phosphate, 'phosphate_target',  'mmol/m**3','phosphate source')

   
   ! Register external state dependencies
   call self%register_state_dependency(self%id_detritus, 'mortality_target','mmol/m**3','sink for dead matter')
   call self%register_state_dependency(self%id_oxygen,   'oxygen_target'   ,'mmol-O2/m**3','dissolved oxygen pool')

   self%lifecycling = next_in_cycle /= ''
   if (self%lifecycling) then
      call self%register_state_variable(self%id_nextc,'next_C','mmol-N/m**3','next clc biomass', &
            standard_variable=type_bulk_standard_variable('cyanobacteria_'//trim(next_in_cycle)//'_biomass'))
      call self%register_state_variable(self%id_nextg,'next_G','mmol-N/m**3','next clc energy', &
            standard_variable=type_bulk_standard_variable('cyanobacteria_'//trim(next_in_cycle)//'_energy'))
      call self%register_state_variable(self%id_nexts,'next_S','mmol-N/m**3','next clc quota', &
            standard_variable=type_bulk_standard_variable('cyanobacteria_'//trim(next_in_cycle)//'_quota'))
!      call self%register_state_dependency(self%id_nextg,'next_G','mmol-N/m**3','next clc energy')
!      call self%register_state_dependency(self%id_nexts,'next_S','mmol-N/m**3','next clc quota')
!      call self%request_coupling(self%id_nextc,'cyanobacteria_'//trim(next_in_cycle)//'_biomass')
!      call self%request_coupling(self%id_nextg,'uhh_clc_G'//next_in_cycle(1:3))
!      call self%request_coupling(self%id_nexts,'uhh_clc_S'//next_in_cycle(1:3))
   end if

   call self%request_coupling(self%id_ammonium, ammonium_variable)
   call self%request_coupling(self%id_nitrate, nitrate_variable)
   call self%request_coupling(self%id_phosphate, phosphate_variable)
   call self%request_coupling(self%id_detritus,detritus_variable)
   call self%request_coupling(self%id_oxygen, oxygen_variable)
   
   ! Register environmental dependencies
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_I_0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)

   return

99 call self%fatal_error('fabm_uhh_clc','Error reading namelist uhh_clc')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of NPZD model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_uhh_clc), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: ni,am,c,g,s,po,par,I_0,temp,e,q,ep,er,fqc,fq,nut
   real(rk)                   :: uptake,light_capture,s_flux,g_flux,c_flux
   real(rk)                   :: growth,fixation,mortality
   real(rk)                   :: sigma_e,sigma_n,sigma_l,sigma_q,omega
   real(rk), parameter        :: secs_pr_day = 86400.0_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_c,c)                ! biomass
   _GET_(self%id_g,g)                ! energy
   _GET_(self%id_s,s)                ! quota
   _GET_(self%id_nitrate,ni)         ! nitrate
   _GET_(self%id_ammonium,am)        ! ammonium
   _GET_(self%id_phosphate,po)       ! phosphate
   nut = ni + am

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)             ! local photosynthetically active radiation
   _GET_(self%id_temp,temp)           ! local temperature
   _GET_HORIZONTAL_(self%id_I_0,I_0)  ! surface short wave radiation

   omega = (1.0_rk/self%rmatscale + &
     1.0_rk/(exp(3.0_rk/(temp-12.0_rk)-0.5_rk)+exp(-500.0_rk/(temp-12.0_rk)-25.0_rk)))
   sigma_l = omega * self%alpha * par / (omega**2 + self%alpha**2 * par**2)
   sigma_n = nut / (nut + self%kN)

   e = 0.0_rk
   q = 0.0_rk
   if (c > 0.0_rk) then
     e = g/c
     q = s/c
   endif
   ep       = tanh(self%scale-self%scale*e/self%e_sc)
   sigma_e  = tanh(           self%scale*e/self%e_sc)
   fqc      = tanh(self%scale-self%scale*q/self%Qc)
   sigma_q  = tanh(           self%scale*q/self%Qc)

   ! specific rates 
   light_capture = self%lightcapture_factor * sigma_l * ep
   uptake = self%uptake_factor * omega * sigma_n * sigma_e * fqc
   growth = self%growth_factor * omega * sigma_e * sigma_q
   ! fixation_factor is the energy consumption factor for fixation = 3 for heterocysts
   fixation = 0.0_rk
   if (self%n_fixation) fixation = omega * sigma_e

   if (self%lifecycling) then
     ! fluxes to next lifestage
     c_flux = (max(e-self%e_max,0.0_rk) + max(self%e_min-e,0.0_rk) + &
           max(q-self%q_max,0.0_rk) + max(self%q_min-q,0.0_rk))/self%tscale
     s_flux=self%Sflux_per_Cflux * c_flux
     g_flux=self%Gflux_per_Cflux * c_flux
   else
     c_flux=0.0_rk
     g_flux=0.0_rk
     s_flux=0.0_rk
   end if

     ! Set temporal derivatives in next lifecycle
   _SET_ODE_(self%id_nextc, c*c_flux)
   _SET_ODE_(self%id_nexts, s*s_flux)
   _SET_ODE_(self%id_nextg, g*g_flux)

   
   ! Set temporal derivatives
   _SET_ODE_(self%id_c,c*(growth + fixation - self%mort - c_flux))
   _SET_ODE_(self%id_s,c*(uptake - growth) - s*(self%mort + s_flux))
   _SET_ODE_(self%id_g,c*(light_capture - uptake - self%m*fixation - growth) - g*(self%mort + g_flux))

   ! take up only nitrate
   _SET_ODE_(self%id_nitrate,-c*uptake * ni/(ni+am))
   _SET_ODE_(self%id_ammonium,-c*uptake * am/(ni+am))
   _SET_ODE_(self%id_detritus,c*self%mort)
   _SET_ODE_(self%id_phosphate, -self%sr *c*uptake)
   ! add oxygen dynamics

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to biogeochemical
! variables
!
! !INTERFACE:
   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   class (type_uhh_clc), intent(in)     :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !LOCAL VARIABLES:
   real(rk)                     :: c
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_c,c) ! biomass concentration

   ! Self-shading
   _SET_EXTINCTION_(self%kc*c)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC


   end module fabm_uhh_clc

