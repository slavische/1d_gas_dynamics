! Программа предназначеная для решения простейшей газодинамической задачи типа "ударная труба"
! Задаются параметры слева и справа от разрыва, разрыв расположен в точке x = 0 
! Для сравнения результатов предусмотренно построение точного решения.
! Для построения графиков можно использовать gnuplot, готовый скрипт находится в plot.gp
program sample
   use io
   use utils
   use vectors
   use gas_dynamics
   use riemann_problem
   implicit none
   character(len = 40) :: input_file_name
   character(len = 40) :: output_file_name
   double precision, allocatable :: x(:)
   double precision, allocatable :: xc(:)
   type (vector_nonconservative_vars), allocatable :: W(:)
   type (vector_conservative_vars), allocatable :: U(:)
   double precision :: t, dt
   type (vector_nonconservative_vars) :: WL, WR
   double precision :: cL, cR
   integer :: i

   ! Чтение файла с исходными данными
   print *, 'Введите имя файла исходными данными (пример test1.dat):'
   read(*, '(A)') input_file_name
   call read_input_file(trim(adjustl(input_file_name)))

   ! Создание расчетной сетки
   call generate_grid(N_CELLS, X_MIN, X_MAX, x, xc)

   ! Иницилизация решения
   t = 0.
   call init_solution(xc, GAMMA_LEFT, DENSITY_LEFT, VELOCITY_LEFT, PRESSURE_LEFT, &
      GAMMA_RIGHT, DENSITY_RIGHT, VELOCITY_RIGHT, PRESSURE_RIGHT, W)

   do
      dt = calc_time_step(x, W, CFL)
      if ( t + dt > END_TIME ) then
         dt = END_TIME - t
         t = END_TIME
         call godunov_solver(dt, x, W)
         exit
      end if
      t = t + dt
      call godunov_solver(dt, x, W)
   end do

   ! Запись решения в файл
   output_file_name = 'solution.dat'
   call write_solution( xc, W, output_file_name )

   ! Построение точного решения для сравнения результатов
   cL = calc_sound_speed(GAMMA_LEFT, DENSITY_LEFT, PRESSURE_LEFT)
   cR = calc_sound_speed(GAMMA_RIGHT, DENSITY_RIGHT, PRESSURE_RIGHT)

   WL%density = density_left
   WL%velocity = velocity_left
   WL%pressure = pressure_left
   WL%gamma = gamma_left

   WR%density = density_right
   WR%velocity = velocity_right
   WR%pressure = pressure_right
   WR%gamma = gamma_right

   W = exact_solution( xc, END_TIME, cL, WL, cR, WR )

   ! Запись решения в файл
   output_file_name = 'exact_solution.dat'
   call write_solution(xc, W, output_file_name)

end program sample
