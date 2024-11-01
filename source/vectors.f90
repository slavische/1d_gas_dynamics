module vectors
   use iso_fortran_env, only: real64
   implicit none
   private

   public::conservative_vars_t, nonconservative_vars_t, calculate_fluxes
   !> вектор консервативных переменных
   !> TODO: Добавить размерность к полям
   type::conservative_vars_t
      private
      real(kind=real64), public::mass     !< масса
      real(kind=real64), public::momentum !< импульс
      real(kind=real64), public::energy   !< энергия
      real(kind=real64), public::gamma    !< показатель адиабаты
   contains
      private
      procedure, public::to_nonconservative => conservative_t_to_nonconservative_t
   end type

   !> вектор примитивных переменных
   !> TODO: Добавить размерность к полям
   type::nonconservative_vars_t
      private
      real(kind=real64), public::density    ! плотность
      real(kind=real64), public::velocity   ! скорость
      real(kind=real64), public::pressure   ! давление
      real(kind=real64), public::gamma      ! показатель адиабаты
   contains
      private
      procedure, public::to_conservative => nonconservative_t_to_conservative_t
   end type

contains

   function conservative_t_to_nonconservative_t(self) result(nonconservative_vars)
      ! Входные данные
      class(conservative_vars_t), intent(in)::self
      ! Выходные данные
      type(nonconservative_vars_t)::nonconservative_vars
      ! Промежуточные данные
      real(kind=real64)::internal_energy

      nonconservative_vars%density = self%mass
      nonconservative_vars%velocity = self%momentum/self%mass
      internal_energy = self%energy/self%mass - 0.5*nonconservative_vars%velocity**2.
      nonconservative_vars%pressure = nonconservative_vars%density*internal_energy*(self%gamma - 1.)
      nonconservative_vars%gamma = self%gamma

   end function conservative_t_to_nonconservative_t

   ! =======================================================================================================
   ! =======================================================================================================
   ! =======================================================================================================

   function nonconservative_t_to_conservative_t(self) result(conservative_vars)
      ! Входные данные
      class(nonconservative_vars_t), intent(in)::self
      ! Выходные данные
      type(conservative_vars_t) :: conservative_vars
      ! Промежуточные данные
      real(kind=real64)::internal_energy

      conservative_vars%mass = self%density
      conservative_vars%momentum = self%density*self%velocity
      internal_energy = self%pressure/self%density/(self%gamma - 1.)
      conservative_vars%energy = self%density*(internal_energy + 0.5*self%velocity**2.)
      conservative_vars%gamma = self%gamma

   end function nonconservative_t_to_conservative_t

   function calculate_fluxes(self) result(fluxes)
      ! Входные данные
      class(nonconservative_vars_t)::self
      ! Выходные данные
      type(conservative_vars_t)::fluxes
      ! Промежуточные данные
      real(kind=real64)::internal_energy

      fluxes%mass = self%density*self%velocity
      fluxes%momentum = fluxes%mass*self%velocity + self%pressure
      internal_energy = self%pressure/self%density/(self%gamma - 1.)
      fluxes%energy = fluxes%mass*(internal_energy + 0.5*self%velocity**2.) + self%pressure*self%velocity
      fluxes%gamma = self%gamma*self%velocity

   end function calculate_fluxes

end module vectors
