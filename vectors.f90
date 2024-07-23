module vectors
   implicit none
   ! вектор консервативных переменных
   type :: vector_conservative_vars
      double precision :: mass         ! масса
      double precision :: momentum     ! импульс
      double precision :: energy       ! энергия
      double precision :: gamma        ! показатель адиабаты
   end type

   ! вектор примитивных переменных
   type :: vector_nonconservative_vars
      double precision :: density      ! плотность
      double precision :: velocity     ! скорость
      double precision :: pressure     ! давление
      double precision :: gamma        ! показатель адиабаты
   end type


contains

   function convert_cons_to_noncons(cons_vars) result(noncons_vars)
      type (vector_conservative_vars), intent(in) :: cons_vars
      double precision :: internal_energy
      type (vector_nonconservative_vars) :: noncons_vars

      noncons_vars%density = cons_vars%mass
      noncons_vars%velocity = cons_vars%momentum / cons_vars%mass
      internal_energy = cons_vars%energy / cons_vars%mass - 0.5 * noncons_vars%velocity ** 2.
      noncons_vars%pressure = noncons_vars%density * internal_energy * ( cons_vars%gamma - 1.)
      noncons_vars%gamma = cons_vars%gamma

   end function convert_cons_to_noncons

   function convert_noncons_to_cons(noncons_vars) result(cons_vars)
      type (vector_nonconservative_vars), intent(in) ::  noncons_vars
      double precision :: internal_energy
      type (vector_conservative_vars) :: cons_vars

      cons_vars%mass = noncons_vars%density
      cons_vars%momentum = noncons_vars%density * noncons_vars%velocity
      internal_energy = noncons_vars%pressure / noncons_vars%density / (noncons_vars%gamma - 1.)
      cons_vars%energy = noncons_vars%density * ( internal_energy + 0.5 * noncons_vars%velocity ** 2. )
      cons_vars%gamma = noncons_vars%gamma

   end function convert_noncons_to_cons

   function calc_vector_fluxes( noncons_vars ) result(fluxes)
      type (vector_nonconservative_vars) :: noncons_vars
      type (vector_conservative_vars) :: fluxes
      double precision :: internal_energy

      fluxes%mass = noncons_vars%density * noncons_vars%velocity
      fluxes%momentum = fluxes%mass * noncons_vars%velocity + noncons_vars%pressure
      internal_energy = noncons_vars%pressure / noncons_vars%density / (noncons_vars%gamma - 1.)
      fluxes%energy = fluxes%mass * ( internal_energy + 0.5 * noncons_vars%velocity ** 2.) + noncons_vars%pressure*noncons_vars%velocity
      fluxes%gamma = noncons_vars%gamma * noncons_vars%velocity
      
   end function calc_vector_fluxes

end module vectors
