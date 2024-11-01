module utils
   implicit none
   private

   public::generate_grid, init_solution, calc_sound_speed, calc_time_step

contains

   subroutine generate_grid(n, x_min, x_max, x, xc)
      integer, intent(in) :: n
      double precision, intent(in) :: x_min, x_max
      double precision, allocatable, intent(out) :: x(:), xc(:)
      double precision :: h
      integer :: i

      allocate (x(n + 1), xc(n))

      h = (x_max - x_min)/n                       ! шаг по пространственной коорлдинате
      x = [(x_min + (i - 1)*h, i=1, n + 1)]   ! координаты границ расчетных ячеек
      xc = [(0.5*(x(i) + x(i + 1)), i=1, n)]  ! координаты центров расчетных ячеек

   end subroutine generate_grid

   subroutine init_solution(xc, gamma_left, density_left, velocity_left, pressure_left, &
                            gamma_right, density_right, velocity_right, pressure_right, solution)
      use vectors
      double precision, intent(in) :: xc(:)
      double precision, intent(in) :: gamma_left, density_left, velocity_left, pressure_left
      double precision, intent(in) :: gamma_right, density_right, velocity_right, pressure_right
      type(nonconservative_vars_t), allocatable, intent(out) :: solution(:)
      integer :: n_cells
      integer :: i

      n_cells = size(xc, dim=1)
      allocate (solution(n_cells))
      do i = 1, n_cells
         if (xc(i) < 0.) then
            solution(i)%density = density_left
            solution(i)%velocity = velocity_left
            solution(i)%pressure = pressure_left
            solution(i)%gamma = gamma_left
         else
            solution(i)%density = density_right
            solution(i)%velocity = velocity_right
            solution(i)%pressure = pressure_right
            solution(i)%gamma = gamma_right
         end if
      end do

   end subroutine init_solution

   function calc_sound_speed(gamma, density, pressure) result(sound_speed)
      double precision, intent(in) :: gamma
      double precision, intent(in) :: density
      double precision, intent(in) :: pressure
      double precision :: sound_speed

      sound_speed = dsqrt(gamma*pressure/density)
   end function calc_sound_speed

   function calc_time_step(x, W, cfl) result(time_step)
      use vectors
      double precision, intent(in) :: x(:)
      type(nonconservative_vars_t), intent(in) :: W(:)
      double precision, intent(in) :: cfl
      double precision :: time_step
      integer :: i
      double precision :: current_step = 1.

      do i = 1, size(W, dim=1)
         time_step = (x(i + 1) - x(i))/dabs(W(i)%velocity + calc_sound_speed(W(i)%gamma, W(i)%density, W(i)%pressure))
         if (time_step < current_step) then
            current_step = time_step
         end if
      end do
      time_step = cfl*current_step
   end function calc_time_step

end module utils
