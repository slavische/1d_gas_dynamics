module gas_dynamics
   use vectors
   implicit none

contains

   function calc_intercell_fluxes( W ) result(intercell_fluxes)
      use utils
      use riemann_problem
      type (vector_nonconservative_vars), intent(in) :: W(:)
      type (vector_conservative_vars), allocatable :: intercell_fluxes(:)
      type (vector_nonconservative_vars) :: WL, WR
      double precision :: cL, cR
      double precision :: pressure_cont, velocity_cont
      integer :: i, n_cells

      n_cells = size(W, dim = 1)
      allocate( intercell_fluxes(n_cells + 1) )

      do i = 1, n_cells + 1

         if ( i == 1) then
            ! реализация граничного условия типа "стенка" на левой границе расчетной области
            WL%density = W(i)%density
            WL%velocity = W(i)%velocity
            WL%pressure =  W(i)%pressure
            WL%gamma = W(i)%gamma
         else
            WL = W(i - 1)
         end if

         if ( i == n_cells + 1 ) then
            ! реализация граничного условия типа "стенка" на правой границе расчетной области
            WR%density = W(i - 1)%density
            WR%velocity = W(i - 1)%velocity
            WR%pressure =  W(i - 1)%pressure
            WR%gamma = W(i - 1)%gamma
         else
            WR = W(i)
         end if

         cL = calc_sound_speed( WL%gamma, WL%density, WL%pressure )
         cR = calc_sound_speed( WR%gamma, WR%density, WR%pressure )
         call calc_contact_pressure_and_velocity( cL, WL, cR, WR, pressure_cont, velocity_cont )
         intercell_fluxes(i) = calc_vector_fluxes( sample_solution( 0.d0, cL, WL, cR, WR, pressure_cont, velocity_cont ) )

      end do

   end function calc_intercell_fluxes

   subroutine godunov_solver(dt, x, W)
      double precision, intent(in) :: dt
      double precision, intent(in) :: x(:)
      type (vector_nonconservative_vars), intent(inout) :: W(:)
      type (vector_conservative_vars), allocatable :: U(:), F(:)
      integer :: i, n_cells
      double precision :: dx

      n_cells = size( W, dim = 1 )
      allocate( U(n_cells), F(n_cells + 1) )
      F = calc_intercell_fluxes( W )

      do i = 1, n_cells
         U(i) = convert_noncons_to_cons( W(i) )
         dx = (x(i + 1) - x(i))
         U(i)%mass = U(i)%mass - dt/dx * ( F(i + 1)%mass - F(i)%mass )
         U(i)%momentum = U(i)%momentum - dt/dx * ( F(i + 1)%momentum - F(i)%momentum )
         U(i)%energy = U(i)%energy - dt/dx * ( F(i + 1)%energy - F(i)%energy )
         W(i) = convert_cons_to_noncons( U(i) )
      end do

   end subroutine godunov_solver

end module gas_dynamics
