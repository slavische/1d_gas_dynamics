module io
   implicit none
   double precision :: X_MIN            ! координата левой границы расчетной области
   double precision :: X_MAX            ! координата правой границы расчетной области
   integer :: N_CELLS                   ! число ячеек в расчетной области
   double precision :: DENSITY_LEFT     ! плотность слева от разрыва
   double precision :: VELOCITY_LEFT    ! скорость слева от разрыва
   double precision :: PRESSURE_LEFT    ! давление слева от разрыва
   double precision :: GAMMA_LEFT       ! показатель адиабаты слева от разрыва
   double precision :: DENSITY_RIGHT    ! плотность справа от разрыва
   double precision :: VELOCITY_RIGHT   ! скорость справа от разрыва
   double precision :: PRESSURE_RIGHT   ! давление справа от разрыва
   double precision :: GAMMA_RIGHT      ! показатель адиабаты справа от разрыва
   double precision :: END_TIME         ! конечный момент времени
   double precision :: TIME_STEP        ! шаг по времени
   double precision :: CFL              ! число Куранта
   integer :: CFL0_NUM_STEP             ! число шагов по времени, когда используется маеньшее число Куранта CFL0
   double precision :: CFL0             ! число Куранта

contains

   subroutine read_input_file(file_name)
      implicit none
      character(len = 20), intent(in) :: file_name

      open(unit = 10, file = file_name, status = 'old', action = 'read')

      read(10,*)
      read(10, '(22x, ES12.5)') X_MIN
      read(10, '(22x, ES12.5)') X_MAX
      read(10, '(22x, I12)') N_CELLS
      read(10,*)
      read(10,*)
      read(10, '(22x, ES12.5)') DENSITY_LEFT
      read(10, '(22x, ES12.5)') VELOCITY_LEFT
      read(10, '(22x, ES12.5)') PRESSURE_LEFT
      read(10, '(22x, ES12.5)') GAMMA_LEFT
      read(10,*)
      read(10, '(22x, ES12.5)') DENSITY_RIGHT
      read(10, '(22x, ES12.5)') VELOCITY_RIGHT
      read(10, '(22x, ES12.5)') PRESSURE_RIGHT
      read(10, '(22x, ES12.5)') GAMMA_RIGHT
      read(10,*)
      read(10, '(22x, ES12.5)') END_TIME
      read(10, '(22x, ES12.5)') TIME_STEP
      read(10, '(22x, ES12.5)') CFL
      read(10, '(22x, I12)') CFL0_NUM_STEP
      read(10, '(22x, ES12.5)') CFL0

      close(10)

   end subroutine read_input_file

   subroutine write_solution(xc, solution, file_name)
      use vectors
      double precision, intent(in) :: xc(:)
      type (vector_nonconservative_vars), intent(in) :: solution(:)
      character(len = 20), intent(in) :: file_name
      integer :: i

      open(unit = 10, file = file_name, status = 'REPLACE', action = 'WRITE')
100   format(5(2x, ES12.5))
      do i = 1, size(xc, dim = 1)
         write(10, 100) xc(i), solution(i)%density, solution(i)%velocity, solution(i)%pressure, solution(i)%gamma
      end do


   end subroutine write_solution
end module io
