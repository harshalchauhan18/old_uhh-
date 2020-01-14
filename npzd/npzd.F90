#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_uhh_npzd - NPZD model by Oschlies&Schartau(2005)
!
! !INTERFACE:
   module fabm_uhh_npzd
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
   type,extends(type_base_model),public :: type_uhh_npzd
!     Variable identifiers
      type (type_state_variable_id)        :: id_n,id_p,id_d,id_z
      type (type_state_variable_id)        :: id_dic
      type (type_dependency_id)            :: id_par,id_temp
      type (type_diagnostic_variable_id)   :: id_dpar
      type (type_diagnostic_variable_id)   :: id_primprod

!     Model parameters
      real(rk) :: kc
      real(rk) :: ws_phy
      real(rk) :: ws_det
      real(rk) :: alpha
      real(rk) :: mort,rem,agg
      real(rk) :: kN
      real(rk) :: tfc
      real(rk) :: mumax
      real(rk) :: gmax
      real(rk) :: eps
      real(rk) :: mort_z
      real(rk) :: slopf
      real(rk) :: depo
      real(rk) :: sfl_ni
      real(rk) :: exudation
      real(rk) :: excretion
      logical  :: do_bottom_flux=.false.
      real(rk) :: N_pool
      real(rk) :: D_pool
      real(rk) :: N_pool_timescale

      contains

      procedure :: initialize
      procedure :: do
      procedure :: do_surface
      procedure :: do_bottom
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
!  Here, the uhh_npzd namelist is read and variables exported
!  by the model are registered with FABM.
!
! !INPUT PARAMETERS:
   class (type_uhh_npzd), intent(inout), target :: self
   integer,                        intent(in)          :: configunit
!
! !LOCAL VARIABLES:
   real(rk)                  :: n_initial,p_initial,d_initial
   real(rk)                  :: z_initial
   real(rk)                  :: ws_phy,ws_det
   real(rk)                  :: kc
   real(rk)                  :: mortality_rate
   real(rk)                  :: remineralization_rate
   real(rk)                  :: aggregation_rate
   real(rk)                  :: gmax
   real(rk)                  :: eps
   real(rk)                  :: alpha
   real(rk)                  :: kN
   real(rk)                  :: mumax
   real(rk)                  :: depo
   real(rk)                  :: tfc
   real(rk)                  :: mort_z
   real(rk)                  :: slopf
   real(rk)                  :: sfl_ni
   real(rk)                  :: N_pool
   real(rk)                  :: D_pool
   real(rk)                  :: N_pool_timescale
   real(rk)                  :: exudation
   real(rk)                  :: excretion
   character(len=64)         :: dic_variable


   namelist /uhh_npzd/ &
             n_initial,p_initial,d_initial, sfl_ni, N_pool, &
             ws_phy, ws_det, kc, mumax, alpha, kN, tfc, &
             dic_variable, exudation, D_pool, N_pool_timescale, &
             z_initial, gmax, eps, mort_z, excretion, slopf, &
             mortality_rate, remineralization_rate, aggregation_rate, depo

!EOP
!-----------------------------------------------------------------------
!BOC
   n_initial = 0.0_rk
   p_initial = 0.0_rk
   d_initial = 0.0_rk
   z_initial = 0.0_rk
   N_pool = 20.0_rk
   D_pool = 0.03_rk
   N_pool_timescale = 3600.0_rk
   slopf = 1.0_rk
   exudation = 0.1_rk
   dic_variable = ''
   sfl_ni = 0.0_rk
   depo = 0.0_rk

   ! Read the namelist
   if (configunit>=0) read(configunit,nml=uhh_npzd,err=99)

   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%kc,     'kc',     default=kc)
   call self%get_parameter(self%tfc,    'tfc',    default=tfc)
   call self%get_parameter(self%alpha,  'alpha',  default=alpha,scale_factor=one_pr_day)
   call self%get_parameter(self%mumax,  'mumax',  default=mumax,scale_factor=one_pr_day)
   call self%get_parameter(self%kN,     'kN',     default=kN)
   call self%get_parameter(self%rem,    'rem',    default=remineralization_rate,scale_factor=one_pr_day)
   call self%get_parameter(self%mort,   'mort',   default=mortality_rate,scale_factor=one_pr_day)
   call self%get_parameter(self%slopf,  'slopf',  default=slopf)
   call self%get_parameter(self%mort_z, 'mort_z', default=mort_z,scale_factor=one_pr_day)
   call self%get_parameter(self%agg,    'agg',    default=aggregation_rate,scale_factor=one_pr_day)
   call self%get_parameter(self%depo,   'depo',   default=depo,scale_factor=one_pr_day)
   call self%get_parameter(self%sfl_ni, 'sfl_ni', default=sfl_ni,scale_factor=one_pr_day)
   call self%get_parameter(self%N_pool, 'N_pool', default=N_pool)
   call self%get_parameter(self%D_pool, 'D_pool', default=D_pool)
   call self%get_parameter(self%N_pool_timescale, 'N_pool_timescale', default=N_pool_timescale)
   call self%get_parameter(self%exudation, 'exudation', default=exudation,scale_factor=one_pr_day)
   call self%get_parameter(self%excretion, 'excretion', default=excretion,scale_factor=one_pr_day)
   call self%get_parameter(self%eps,    'eps',    default=eps,scale_factor=one_pr_day)
   call self%get_parameter(self%gmax,   'gmax',   default=gmax,scale_factor=one_pr_day)
   self%do_bottom_flux = N_pool > 0.0_rk
   
   ! Register state variables
   call self%register_state_variable(self%id_n,'nit', &
         'mmol n/m**3','nutrient', n_initial, &
         minimum=0.0e-7_rk)

   call self%register_state_variable(self%id_p,'phy', &
         'mmol n/m**3','phytoplankton', p_initial, &
         minimum=0.0e-7_rk,vertical_movement=ws_phy*one_pr_day)

   call self%register_state_variable(self%id_z,'zoo', &
         'mmol n/m**3','zooplankton', z_initial, &
         minimum=0.0e-7_rk)

   call self%register_state_variable(self%id_d,'det', &
         'mmol n/m**3','detritus', d_initial, &
         minimum=0.0e-7_rk,vertical_movement=ws_det*one_pr_day)
   
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

99 call self%fatal_error('fabm_uhh_npzd','Error reading namelist uhh_npzd')

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
   class (type_uhh_npzd), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_
!
! !LOCAL VARIABLES:
   real(rk)                   :: n,p,d,z,dic
   real(rk)                   :: par,temp
   real(rk)                   :: mumax,mu
   real(rk)                   :: sigma_n,sigma_l
   real(rk)                   :: graz,excretion,ta,zmort

!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_n,n)                ! nutrient
   _GET_(self%id_p,p)                ! phytoplankton
   _GET_(self%id_z,z)                ! zooplankton
   _GET_(self%id_d,d)                ! detritus

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)             ! local photosynthetically active radiation
   _GET_(self%id_temp,temp)           ! local temperature

   ta = self%tfc**temp
   mumax = self%mumax*ta
   sigma_l = mumax * self%alpha*par / sqrt(mumax**2 + (self%alpha*par)**2)
   sigma_n  = mumax * n/(self%kn + n)
   mu = min(sigma_l, sigma_n)
   graz = z * self%gmax*self%eps*p**2 / (self%gmax + self%eps*p**2)
   excretion = self%excretion * ta * z
   zmort = self%mort_z*z*z

   ! Set temporal derivatives
   _SET_ODE_(self%id_n,(self%mort*ta - mu)*p + self%rem*ta*d + excretion)
   _SET_ODE_(self%id_p,mu*p -(self%mort*ta + self%agg*p)*p - graz)
   _SET_ODE_(self%id_d,self%agg*p*p - self%rem*ta*d + zmort + (1.0_rk-self%slopf)*graz)
   _SET_ODE_(self%id_z,-zmort - excretion + self%slopf*graz)

   ! set DIC sink, if available
   if (_AVAILABLE_(self%id_dic)) &
     _SET_ODE_(self%id_dic,c_to_n*(self%rem*d - sigma_n*sigma_l*p - self%exudation*p))

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dpar,par)
   _SET_DIAGNOSTIC_(self%id_primprod,mu*p)
   
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
   class (type_uhh_npzd), intent(in)     :: self
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
   class (type_uhh_npzd),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_PPDD_

   real(rk)                   :: n,p,d,z
   real(rk)                   :: par,temp,ta
   real(rk)                   :: sigma_n,sigma_l,mumax,mu
   real(rk)                   :: graz, excretion, zmort

   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_n,n)                ! nutrient
   _GET_(self%id_p,p)                ! phytoplankton
   _GET_(self%id_z,z)                ! zooplankton
   _GET_(self%id_d,d)                ! detritus

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)             ! local photosynthetically active radiation
   _GET_(self%id_temp,temp)           ! local temperature

   ta = self%tfc**temp
   mumax   = self%mumax*ta
   sigma_l = mumax * self%alpha*par / sqrt(mumax**2 + self%alpha**2 * par**2)
   sigma_n  = mumax * n/(self%kn + n)
   mu = min(sigma_l, sigma_n)
   graz = z * self%gmax*self%eps*p**2 / (self%gmax + self%eps*p**2)
   excretion = self%excretion * ta * z
   zmort = self%mort_z*z*z

   ! Set temporal derivatives
   _SET_DD_SYM_(self%id_n,self%id_p,mu*p)
   _SET_DD_SYM_(self%id_d,self%id_n,self%rem*ta*d)
   _SET_DD_SYM_(self%id_z,self%id_n,excretion)
   _SET_DD_SYM_(self%id_p,self%id_z,graz)
   _SET_DD_SYM_(self%id_p,self%id_n,self%mort*ta*p)
   _SET_DD_SYM_(self%id_z,self%id_d,(1.0_rk-self%slopf)*graz + zmort)
   _SET_DD_SYM_(self%id_p,self%id_d,self%agg*p*p)

   if (_AVAILABLE_(self%id_dic)) &
     _SET_DD_(self%id_dic,self%id_dic,c_to_n*(sigma_n*sigma_l*p - self%rem*d - self%exudation*p))

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_dpar,par)
   _SET_DIAGNOSTIC_(self%id_primprod,mu*p)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do_ppdd


   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_uhh_npzd), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
   real(rk) :: d,n

   if (.not.(self%do_bottom_flux)) return
   _FABM_HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_d,d)                ! detritus
   _GET_(self%id_n,n)                ! nutrients

   _SET_BOTTOM_EXCHANGE_(self%id_d,(self%D_pool-d)/self%N_pool_timescale)
   _SET_BOTTOM_EXCHANGE_(self%id_n,(self%N_pool-n)/self%N_pool_timescale)

   _FABM_HORIZONTAL_LOOP_END_

   end subroutine do_bottom


   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class(type_uhh_npzd),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   _FABM_HORIZONTAL_LOOP_BEGIN_
   _SET_SURFACE_EXCHANGE_(self%id_n,self%sfl_ni)
   _FABM_HORIZONTAL_LOOP_END_

   end subroutine do_surface



   end module fabm_uhh_npzd

