#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_uhh_pnd - PND model
!
! !INTERFACE:
   module fabm_uhh_pnd
!
! !DESCRIPTION:
!
! PND model.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_uhh_pnd
!     Variable identifiers
      type (type_state_variable_id)        :: id_n,id_p,id_d
      type (type_state_variable_id)        :: id_dic
      type (type_dependency_id)            :: id_par,id_temp
      type (type_diagnostic_variable_id)   :: id_dpar
      type (type_diagnostic_variable_id)   :: id_primprod

!     Model parameters
      real(rk) :: kc
      real(rk) :: ws_phy
      real(rk) :: ws_det
      real(rk) :: alpha
      real(rk) :: mort,rem,graz
      real(rk) :: kN
      real(rk) :: tfc
      real(rk) :: omega0
      real(rk) :: depo
      real(rk) :: sfl_ni
      real(rk) :: exudation
      logical  :: do_bottom_flux=.false.
      real(rk) :: N_pool
      real(rk) :: N_pool_timescale

      contains

      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      procedure :: do_ppdd
      procedure :: do_bottom
      procedure :: do_bottom_ppdd
      procedure :: get_light_extinction
   end type

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk
   real(rk), parameter :: c_to_n = 106.0_rk/16.0_rk
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the PND model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the uhh_pnd namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_uhh_pnd), intent(inout), target :: self
   integer,                        intent(in)          :: configunit
!
! !LOCAL VARIABLES:
   real(rk)                  :: n_initial,p_initial,d_initial
   real(rk)                  :: ws_phy,ws_det
   real(rk)                  :: kc
   real(rk)                  :: mortality_rate
   real(rk)                  :: remineralization_rate
   real(rk)                  :: grazing_rate
   real(rk)                  :: alpha
   real(rk)                  :: kN
   real(rk)                  :: omega0
   real(rk)                  :: depo
   real(rk)                  :: tfc
   real(rk)                  :: sfl_ni
   real(rk)                  :: N_pool
   real(rk)                  :: N_pool_timescale
   real(rk)                  :: exudation
   character(len=64)         :: dic_variable


   namelist /uhh_pnd/ &
             n_initial,p_initial,d_initial, sfl_ni, N_pool, &
             ws_phy, ws_det, kc, omega0, alpha, kN, tfc, &
             dic_variable, exudation, N_pool_timescale, &
             mortality_rate, remineralization_rate, grazing_rate, depo

!EOP
!-----------------------------------------------------------------------
!BOC
   n_initial = 0.0_rk
   p_initial = 0.0_rk
   d_initial = 0.0_rk
   N_pool = 20.0_rk
   N_pool_timescale = 3600.0_rk
   exudation = 0.1_rk
   dic_variable = ''

   ! Read the namelist
   if (configunit>=0) read(configunit,nml=uhh_pnd,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%kc,     'kc',     default=kc)
   call self%get_parameter(self%tfc,    'tfc',    default=tfc)
   call self%get_parameter(self%alpha,  'alpha',  default=alpha)
   call self%get_parameter(self%omega0, 'omega0', default=omega0)
   call self%get_parameter(self%kN,     'kN',     default=kN)
   call self%get_parameter(self%rem,    'rem',    default=remineralization_rate)
   call self%get_parameter(self%mort,   'mort',   default=mortality_rate)
   call self%get_parameter(self%graz,   'graz',   default=grazing_rate)
   call self%get_parameter(self%depo,   'depo',   default=depo)
   call self%get_parameter(self%sfl_ni, 'sfl_ni', default=sfl_ni)
   call self%get_parameter(self%N_pool, 'N_pool', default=N_pool)
   call self%get_parameter(self%N_pool_timescale, 'N_pool_timescale', default=N_pool_timescale)
   call self%get_parameter(self%exudation, 'exudation', default=exudation,scale_factor=one_pr_day)
   self%do_bottom_flux = N_pool > 0.0_rk
   
   ! Register state variables
   call self%register_state_variable(self%id_n,'nit', &
         'mmol n/m**3','nutrient', n_initial, &
         minimum=1.0e-7_rk)

   call self%register_state_variable(self%id_p,'phy', &
         'mmol n/m**3','phytoplankton', p_initial, &
         minimum=1.0e-7_rk,vertical_movement=ws_phy)

   call self%register_state_variable(self%id_d,'det', &
         'mmol n/m**3','detritus', d_initial, &
         minimum=1.0e-7_rk,vertical_movement=ws_det)
   
   ! Register environmental dependencies
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_temp, standard_variables%temperature)

   ! Register optional link to external DIC pool
   call self%register_state_dependency(self%id_dic,'dic','mmol/m**3','total dissolved inorganic carbon',required=.false.)
   if (dic_variable/='') call self%request_coupling(self%id_dic,dic_variable)

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_dpar,'PAR','W/m**2', &
         'PAR', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_primprod,'primprod','mmol-N/m**3/s', &
         'primary production rate', output=output_instantaneous)
   return

99 call self%fatal_error('fabm_uhh_pnd','Error reading namelist uhh_pnd')

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Right hand sides of PND model
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
! !INPUT PARAMETERS:
   class (type_uhh_pnd), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: n,p,d,dic
   real(rk)                   :: par,temp
   real(rk)                   :: sigma_n,sigma_l,omega

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_n,n)                ! nutrient
   _GET_(self%id_p,p)                ! phytoplankton
   _GET_(self%id_d,d)                ! detritus

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)             ! local photosynthetically active radiation
   _GET_(self%id_temp,temp)           ! local temperature

   omega   = self%omega0*(self%tfc**temp)
   sigma_l = omega * self%alpha*par / sqrt(omega**2 + self%alpha**2 * par**2)
   sigma_n  = n/(self%kn + n)

   ! Set temporal derivatives
   _SET_ODE_(self%id_n,-sigma_n*sigma_l*p + self%rem*d)
   _SET_ODE_(self%id_p,sigma_n*sigma_l*p -(self%mort + self%graz*p)*p)
   _SET_ODE_(self%id_d,(self%mort + self%graz*p)*p - self%rem*d)

   ! set DIC sink, if available
   if (_AVAILABLE_(self%id_dic)) &
     _SET_ODE_(self%id_dic,c_to_n*(self%rem*d - sigma_n*sigma_l*p - self%exudation*p))

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dpar,par)
   _SET_DIAGNOSTIC_(self%id_primprod,sigma_n*sigma_l*p)
   
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
   class (type_uhh_pnd), intent(in)     :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !LOCAL VARIABLES:
   real(rk)                     :: p,d
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_p,p)                ! phytoplankton
   _GET_(self%id_d,d)                ! detritus

   ! Self-shading
   _SET_EXTINCTION_(self%kc*(p+d))

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC

   subroutine do_ppdd(self,_ARGUMENTS_DO_PPDD_)
   class (type_uhh_pnd),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_

   real(rk)                   :: n,p,d
   real(rk)                   :: par,temp
   real(rk)                   :: sigma_n,sigma_l,omega


   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_n,n)                ! nutrient
   _GET_(self%id_p,p)                ! phytoplankton
   _GET_(self%id_d,d)                ! detritus

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)             ! local photosynthetically active radiation
   _GET_(self%id_temp,temp)           ! local temperature

   omega   = self%omega0*(self%tfc**temp)
   sigma_l = omega * self%alpha*par / sqrt(omega**2 + self%alpha**2 * par**2)
   sigma_n  = n/(self%kn + n)

   ! Set temporal derivatives
   _SET_DD_SYM_(self%id_n,self%id_p,sigma_n*sigma_l*p)
   _SET_DD_SYM_(self%id_d,self%id_n,self%rem*d)
   _SET_DD_SYM_(self%id_p,self%id_d,(self%mort + self%graz*p)*p)

   if (_AVAILABLE_(self%id_dic)) &
     _SET_DD_(self%id_dic,self%id_dic,c_to_n*(sigma_n*sigma_l*p - self%rem*d - self%exudation*p))

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dpar,par)
   _SET_DIAGNOSTIC_(self%id_primprod,sigma_n*sigma_l*p)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do_ppdd


   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_uhh_pnd), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
   real(rk) :: d,n

   if (.not.(self%do_bottom_flux)) return
   _FABM_HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_d,d)                ! detritus
   _GET_(self%id_n,n)                ! nutrients

   ! in bio_clc code, self%depo was given in units 1/s per mmol-N/m3 biomass,
   ! here, self%depo is of units m/s per mmol-N/m3 biomass. The results should
   ! be comparable, since the bottom layer height in the calibrated setup
   ! is appr. 1.03 m.
   _SET_BOTTOM_EXCHANGE_(self%id_d,-self%depo*d*d)
   _SET_BOTTOM_EXCHANGE_(self%id_n,(self%N_pool-n)/self%N_pool_timescale)

   _FABM_HORIZONTAL_LOOP_END_

   end subroutine do_bottom


   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class(type_uhh_pnd),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
   real(rk),parameter :: secs_pr_day=86400._rk

   _FABM_HORIZONTAL_LOOP_BEGIN_
   _SET_SURFACE_EXCHANGE_(self%id_n,self%sfl_ni/secs_pr_day)
   _FABM_HORIZONTAL_LOOP_END_

   end subroutine do_surface


   subroutine do_bottom_ppdd(self,_ARGUMENTS_DO_BOTTOM_PPDD_)
   class (type_uhh_pnd), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_PPDD_
   real(rk) :: d,n

   if (.not.(self%do_bottom_flux)) return
   _FABM_HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_d,d)                ! detritus
   _GET_(self%id_n,n)                ! nutrients

   _SET_BOTTOM_EXCHANGE_(self%id_d,self%depo*d*d)
   _SET_BOTTOM_EXCHANGE_(self%id_n,(20.0-n)/3600.)

   _FABM_HORIZONTAL_LOOP_END_

   end subroutine do_bottom_ppdd

#if 0
   subroutine do_surface_ppdd(self,_ARGUMENTS_DO_SURFACE)
   class(type_uhh_pnd) :: self
  _DECLARE_ARGUMENTS_DO_SURFACE_
   real(rk),parameter :: secs_pr_day=86400._rk

   _FABM_HORIZONTAL_LOOP_BEGIN_
   _SET_SURFACE_EXCHANGE_(self%id_n,self%sfl_ni/secs_pr_day)
   _FABM_HORIZONTAL_LOOP_END_

   end subroutine do_surface_ppdd
#endif

   end module fabm_uhh_pnd

