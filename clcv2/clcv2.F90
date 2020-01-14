#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_uhh_clcv2 - CLC-v2 cyanobacteria lifestage base model
!
! !INTERFACE:
   module fabm_uhh_clcv2
!
! !DESCRIPTION:
!
! reduced version of a cyanobacteria lifestage model
! after Hense and Beckmann (2010)
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_uhh_clcv2
!     Variable identifiers
      type (type_state_variable_id)        :: id_ch,id_cn,id_aki,id_gh
      type (type_state_variable_id)        :: id_nitrate,id_ammonium
      type (type_state_variable_id)        :: id_phosphate,id_detritus,id_oxygen
      type (type_dependency_id)            :: id_par,id_temp
      type (type_diagnostic_variable_id)   :: id_sl,id_eh
      type (type_horizontal_diagnostic_variable_id) :: id_tmat0d,id_cmax
      type (type_bottom_state_variable_id) :: id_cmaxint,bsv_cmax
      type (type_global_dependency_id)     :: id_doy

!     Model parameters
      real(rk) :: sr
      real(rk) :: alpha
      real(rk) :: kc
      real(rk) :: s2
      real(rk) :: s3
      real(rk) :: mort
      real(rk) :: w
      real(rk) :: w_rec
      real(rk) :: w_aki
      real(rk) :: Qc
      real(rk) :: depo
      real(rk) :: omega0
      real(rk) :: rmatscale
      real(rk) :: tmat_crit
      real(rk) :: kN
      real(rk) :: tfc_c=25.0_rk
      real(rk) :: h2a_factor
      logical  :: use_phosphate
      logical  :: use_ammonium
      logical  :: use_oxygen

      contains

      procedure :: initialize
      procedure :: do
      procedure :: do_bottom
      procedure :: check_bottom_state
      procedure :: check_surface_state
      procedure :: get_light_extinction
      procedure :: sigma_from_forcing
   end type

   type :: type_sigma
     real(rk) :: l,n,t
   end type

   real(rk), parameter :: secs_pr_day = 86400.0_rk
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk
!EOP
!-----------------------------------------------------------------------

   contains


   !> Initialise the CLC model
   !
   !>  Here, the uhh_clc namelist is read and variables exported
   !!  by the model are registered with FABM.
   subroutine initialize(self,configunit)
   class (type_uhh_clcv2), intent(inout), target :: self
   integer,                        intent(in)          :: configunit

   real(rk)                  :: cn_initial,tmat_initial,ch_initial
   real(rk)                  :: aki_initial,gh_initial
   real(rk)                  :: w,w_aki,w_rec
   real(rk)                  :: h2a_factor
   real(rk)                  :: kc
   real(rk)                  :: mortality_rate
   real(rk)                  :: alpha
   real(rk)                  :: sr,s2,s3
   real(rk)                  :: Qc
   real(rk)                  :: kN
   real(rk)                  :: omega0
   real(rk)                  :: rmatscale
   real(rk)                  :: depo
   real(rk)                  :: tmat_crit
   character(len=64)         :: phosphate_variable
   character(len=64)         :: ammonium_variable
   character(len=64)         :: nitrate_variable
   character(len=64)         :: detritus_variable
   character(len=64)         :: oxygen_variable

   namelist /uhh_clcv2/ &
             cn_initial,tmat_initial,ch_initial, &
             aki_initial, depo, h2a_factor, &
             w, kc, sr, s2, s3, rmatscale, omega0, alpha, &
             Qc, kN, mortality_rate, tmat_crit, &
             ammonium_variable,nitrate_variable, w_rec, w_aki, &
             phosphate_variable, detritus_variable, &
             oxygen_variable

   cn_initial = 0.0_rk
   ch_initial = 0.0_rk

   w         = 0.1_rk     ! m/d
   w_rec     = 2.0_rk     ! m/d
   w_aki     = -10.0_rk   ! m/d
   tmat_crit = 100.0_rk   ! d
   kc        = 0.03_rk    ! 1/m/(mmolN/m3)
   sr        = 0.0625_rk  ! molP/molN
   s2        = 6.625_rk
   s3        = 8.625_rk 
   Qc        = 0.0_rk
   kN        = 0.0_rk     ! mmolN/m3
   rmatscale = 45.0_rk
   omega0    = 0.24_rk    ! 1/d
   depo      = 0.0_rk     ! 1/m2/s per (mmolN/m3)^2
   mortality_rate = 0.02_rk ! mort [1/d]
   h2a_factor = 4.0_rk    ! [1/1]

   nitrate_variable = 'uhh_ergom_split_base_nit'
   ammonium_variable = 'uhh_ergom_split_base_amm'
   phosphate_variable = 'uhh_ergom_split_base_pho'
   detritus_variable = 'uhh_ergom_split_base_det'
   oxygen_variable = 'uhh_ergom_split_base_oxy'

   ! Read the namelist
   if (configunit>=0) read(configunit,nml=uhh_clcv2,err=99)

   ! set dependency switches
   self%use_phosphate = phosphate_variable /= ''
   self%use_ammonium  = ammonium_variable /= ''
   self%use_oxygen    = oxygen_variable /= ''


   ! Store parameter values in our own derived type
   ! NB: all rates must be provided in values per day,
   ! and are converted here to values per second.
   call self%get_parameter(self%kc,     'kc',     default=kc)
   call self%get_parameter(self%alpha,  'alpha',  default=alpha)
   call self%get_parameter(self%sr,     'sr',     default=sr)
   call self%get_parameter(self%s2,     's2',     default=s2)
   call self%get_parameter(self%s3,     's3',     default=s3)
   call self%get_parameter(self%kN,     'kN',     default=kN)
   call self%get_parameter(self%Qc,     'Qc',     default=Qc)
   call self%get_parameter(self%depo,   'depo',   default=depo)
   call self%get_parameter(self%h2a_factor, 'h2a_factor', default=h2a_factor)
   call self%get_parameter(self%omega0, 'omega0', default=omega0, scale_factor=one_pr_day)
   call self%get_parameter(self%rmatscale, 'rmatscale', default=rmatscale)
   call self%get_parameter(self%tmat_crit, 'tmat_crit', default=tmat_crit)
   call self%get_parameter(self%w,      'w',      default=w,  scale_factor=one_pr_day)
   call self%get_parameter(self%w_rec,  'w_rec',  default=w_rec,  scale_factor=one_pr_day)
   call self%get_parameter(self%w_aki,  'w_aki',  default=w_aki,  scale_factor=one_pr_day)
   call self%get_parameter(self%mort,   'mortality_rate', default=mortality_rate,  scale_factor=one_pr_day)

   ! Register state variables
   call self%register_state_variable(self%id_cn,'Cn', &
         'mmol n/m**3','non-diazotrophic biomass', cn_initial, &
         minimum=0.0e-7_rk,vertical_movement=self%w_rec)

   call self%register_state_variable(self%id_ch,'Ch', &
         'mmol n/m**3','heterocysts biomass', ch_initial, &
         minimum=0.0e-7_rk,vertical_movement=self%w)

   call self%register_state_variable(self%id_aki,'Aki', &
         'mmol n/m**3','akinetes biomass', aki_initial, &
         minimum=0.0e-7_rk,vertical_movement=self%w_aki)

   call self%register_state_variable(self%id_cmaxint,'cmaxint', &
         'mmol n/m**3','time-integrated maximum bottom biomass', 1.0e3_rk*tmat_initial)

   call self%register_state_variable(self%bsv_cmax,'cmax', &
         'mmol n/m**3','maximum bottom biomass', 1.0e3_rk)

   ! Register dependencies on external standard variables
   if (self%use_ammonium) &
     call self%register_state_dependency(self%id_ammonium, 'ammonium_target', 'mmol/m**3','ammonium source')
   call self%register_state_dependency(self%id_nitrate, 'nitrate_target', 'mmol/m**3','nitrate source')
   if (self%use_phosphate) &
     call self%register_state_dependency(self%id_phosphate, 'phosphate_target',  'mmol/m**3','phosphate source')

   
   ! Register external state dependencies
   call self%register_state_dependency(self%id_detritus, 'mortality_target','mmol/m**3','sink for dead matter')
   if (self%use_oxygen) &
     call self%register_state_dependency(self%id_oxygen,   'oxygen_target'   ,'mmol-O2/m**3','dissolved oxygen pool')

   if (self%use_ammonium) &
     call self%request_coupling(self%id_ammonium, ammonium_variable)
   call self%request_coupling(self%id_nitrate, nitrate_variable)
   if (self%use_phosphate) &
     call self%request_coupling(self%id_phosphate, phosphate_variable)
   call self%request_coupling(self%id_detritus,detritus_variable)
   if (self%use_oxygen) &
     call self%request_coupling(self%id_oxygen, oxygen_variable)
   
   ! Register environmental dependencies
   call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_doy,standard_variables%number_of_days_since_start_of_the_year)

   ! Register diagnostic variables
   call self%register_diagnostic_variable(self%id_sl,'sigma_l','', &
      'light limitation factor', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tmat0d,'tmat0d','d', &
      '0d maturation time', output=output_instantaneous)

   return

99 call self%fatal_error('fabm_uhh_clcv2','Error reading namelist uhh_clcv2')

   end subroutine initialize



   !> Right hand sides of CLCv2 model
   subroutine do(self,_ARGUMENTS_DO_)
   class (type_uhh_clcv2), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_

   real(rk)                   :: ni,am,cn,ch,tmat,d,po,par,temp,nut
   real(rk)                   :: aki,light_capture,n_fixation
   real(rk)                   :: cmax,cmaxint
   type(type_sigma)           :: sigma
   real(rk)                   :: growth_h,growth_n
   real(rk)                   :: h2a_flux,a2n_flux,n2h_flux

   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_ch,ch)                ! heterocysts biomass
   _GET_(self%id_cn,cn)                ! non-diazotrophic biomass
   _GET_(self%id_aki,aki)              ! akinetes biomass
   _GET_(self%id_nitrate,ni)           ! nitrate
   _GET_(self%id_detritus,d)           ! detritus
   if (self%use_ammonium) then
     _GET_(self%id_ammonium,am)        ! ammonium
   else
     am=0.0_rk
   end if
   nut = ni + am

   if (self%use_phosphate) then
     _GET_(self%id_phosphate,po)       ! phosphate
   else
     po = nut/16.0_rk
   end if

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)              ! local photosynthetically active radiation
   _GET_(self%id_temp,temp)            ! local temperature

   sigma = self%sigma_from_forcing(nut,par,temp)

   ! specific rates 
   growth_h = sigma%l
   growth_n = 3.0_rk * growth_h * sigma%n
   n_fixation = self%omega0 * sigma%t

   ! get and calculate maturation time
   _GET_HORIZONTAL_(self%bsv_cmax,cmax)      ! maximum bottom biomass
   _GET_HORIZONTAL_(self%id_cmaxint,cmaxint) ! integrated maximum biomass
   tmat = cmaxint/(cmax+1.0d-7)

   ! lifestage fluxes
   a2n_flux=0.0_rk
   n2h_flux=0.0_rk
   h2a_flux=0.0_rk

   ! non-diazotrophic -> heretocysts
   if (growth_h > growth_n) n2h_flux = one_pr_day
   ! heterocysts -> akinetes
   if ((3.5_rk*growth_n > growth_h).and.(par>0.0_rk)) &
      h2a_flux = 4.0_rk * one_pr_day * sigma%t
   ! akinates -> non-diazotrophic
   if (tmat > self%tmat_crit) a2n_flux = one_pr_day

   ! Set temporal derivatives
   _SET_ODE_(self%id_ch,ch*(growth_h - self%mort) - h2a_flux*ch + n2h_flux*cn)
   _SET_ODE_(self%id_cn,cn*(growth_n - self%mort) - n2h_flux*cn + a2n_flux*aki)
   _SET_ODE_(self%id_aki,-aki*self%mort - a2n_flux*aki + h2a_flux*ch)

   ! external nutrients
   _SET_ODE_(self%id_nitrate,-cn*growth_n * ni/(ni+am))
   if (self%use_ammonium) then
     _SET_ODE_(self%id_ammonium,-cn*growth_n * am/(ni+am))
   end if
   _SET_ODE_(self%id_detritus,(cn+ch)*self%mort + aki*self%mort)
   if (self%use_phosphate) then
     _SET_ODE_(self%id_phosphate, -self%sr *cn*growth_n)
   end if
   ! add oxygen dynamics
   ! add dic dynamics

   ! Export diagnostic variables
   _SET_DIAGNOSTIC_(self%id_sl,sigma%l)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine do

   

   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
   class (type_uhh_clcv2), intent(in)     :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
   
   real(rk)                     :: cn,ch,aki

   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_cn,cn)
   _GET_(self%id_ch,ch)
   _GET_(self%id_aki,aki)

   ! Self-shading
   _SET_EXTINCTION_(self%kc*(cn+ch+aki))

   _LOOP_END_

   end subroutine get_light_extinction



   subroutine check_bottom_state(self,_ARGUMENTS_CHECK_BOTTOM_STATE_)
   class (type_uhh_clcv2), intent(in) :: self
   _DECLARE_ARGUMENTS_CHECK_BOTTOM_STATE_
   real(rk) :: cmax,aki

   _HORIZONTAL_LOOP_BEGIN_
   ! set cmax
   _GET_(self%id_aki,aki)                ! bottom biomass
   _GET_HORIZONTAL_(self%bsv_cmax,cmax)  ! previous maximum biomass
   _SET_HORIZONTAL_(self%bsv_cmax, max(aki,cmax))
   
   _HORIZONTAL_LOOP_END_
     
   end subroutine check_bottom_state



   subroutine check_surface_state(self,_ARGUMENTS_CHECK_SURFACE_STATE_)
   class (type_uhh_clcv2), intent(in) :: self
   _DECLARE_ARGUMENTS_CHECK_SURFACE_STATE_
   real(rk)         :: par,temp,nut,ch,cn,cmax,growth_n,growth_h
   type(type_sigma) :: sigma

   _HORIZONTAL_LOOP_BEGIN_
   _GET_(self%id_ch,ch)                ! heterocysts biomass
   _GET_(self%id_cn,cn)                ! non-diazotrophic biomass
   _GET_(self%id_nitrate,nut)          ! nitrate
   _GET_HORIZONTAL_(self%bsv_cmax, cmax) ! maximum bottom biomass

   ! Retrieve current environmental conditions.
   _GET_(self%id_par,par)              ! local photosynthetically active radiation
   _GET_(self%id_temp,temp)            ! local temperature

   ! calculate growth rates at surface
   sigma = self%sigma_from_forcing(nut,par,temp)
   growth_h = sigma%l
   growth_n = 3.0_rk * growth_h * sigma%n

   if ((3.5_rk*growth_n > growth_h).and.(par>0.0_rk).and.(cmax>1.0d-7).and.(ch>1.0d-1)) then
      _SET_HORIZONTAL_(self%bsv_cmax, 0.0_rk)
      _SET_HORIZONTAL_(self%id_cmaxint, 0.0_rk)
   end if

   _HORIZONTAL_LOOP_END_

   end subroutine check_surface_state



   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
   class (type_uhh_clcv2), intent(in) :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
   real(rk),save :: aki,cmax,cmaxint

   _FABM_HORIZONTAL_LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_aki,aki)                    ! biomass
   _GET_HORIZONTAL_(self%bsv_cmax,cmax)      ! maximum bottom biomass
   _GET_HORIZONTAL_(self%id_cmaxint,cmaxint) ! integrated maximum biomass
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tmat0d,cmaxint/(cmax+1.e-7_rk))

   _SET_BOTTOM_EXCHANGE_(self%id_aki,-self%depo*aki*aki)
   _SET_ODE_BEN_(self%id_cmaxint,cmax*one_pr_day)

   _FABM_HORIZONTAL_LOOP_END_

   end subroutine do_bottom


   function sigma_from_forcing(self,nut,par,temp) result(sigma)
   class (type_uhh_clcv2), intent(in) :: self
   real(rk), intent(in) :: nut,par,temp
   type (type_sigma)    :: sigma
   real(rk)             :: omega

   !RH: nan for t==12.0_rk !
   sigma%t = (1.0_rk/self%rmatscale + &
     1.0_rk/(0.25_rk+exp(3.0_rk/(temp-12.0_rk)-0.5_rk)+exp(-500.0_rk/(temp-12.0_rk)+self%tfc_c)))

   omega = self%omega0 * sigma%t
   sigma%l = omega * self%alpha * par / sqrt(omega**2 + self%alpha**2 * par**2)
   sigma%n = nut / (nut + self%kN)

   end function sigma_from_forcing


   end module fabm_uhh_clcv2

