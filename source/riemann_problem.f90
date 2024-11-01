module riemann_problem
   use vectors
   use iso_fortran_env, only: real64
   implicit none
   private

   public::sample_solution, exact_solution, calc_contact_pressure_and_velocity

contains

   pure function calculate_pressure(c_l, w_l, c_r, w_r, gamma) result(p)
      ! Входные данные
      type(nonconservative_vars_t), intent(in)::w_l !< Вектор потоков "слева"
      type(nonconservative_vars_t), intent(in)::w_r !< Вектор потоков "справа"
      real(kind=real64), intent(in)::c_l
      real(kind=real64), intent(in)::c_r
      real(kind=real64), intent(in)::gamma
      ! Выходные данные
      real(kind=real64)::p
      ! Промежуточные данные
      real(kind=real64), parameter::q_user = 2.0
      real(kind=real64)::g1
      real(kind=real64)::g2
      real(kind=real64)::g3
      real(kind=real64)::g4
      real(kind=real64)::g5
      real(kind=real64)::g6
      real(kind=real64)::g7
      real(kind=real64)::g8
      real(kind=real64)::cup
      real(kind=real64)::ppv
      real(kind=real64)::p_min
      real(kind=real64)::p_max
      real(kind=real64)::q_max
      real(kind=real64)::pq
      real(kind=real64)::um
      real(kind=real64)::ptl
      real(kind=real64)::ptr
      real(kind=real64)::gel
      real(kind=real64)::ger

      ! Eleuterio F. Toro Riemann Solvers and Numerical Methods for Fluid Dynamics
      g1 = (gamma - 1.0)/(2.0*gamma)
      g2 = (gamma + 1.0)/(2.0*gamma)
      g3 = 2.0*gamma/(gamma - 1.0)
      g4 = 2.0/(gamma - 1.0)
      g5 = 2.0/(gamma + 1.0)
      g6 = (gamma - 1.0)/(gamma + 1.0)
      g7 = (gamma - 1.0)/2.0
      g8 = gamma - 1.0

      ! Compute guess pressure from PVRS Riemann solver
      ! Compute gamma related constants
      cup = 0.25*(w_l%density + w_r%density)*(c_l + c_r)
      ppv = 0.5*(w_l%pressure + w_r%pressure) + 0.5*(w_l%velocity - w_r%velocity)*cup
      ppv = max(0.0, ppv)
      p_min = min(w_l%pressure, w_r%pressure)
      p_max = max(w_l%pressure, w_r%pressure)
      q_max = p_max/p_min

      IF (q_max .LE. q_user .AND. (p_min .LE. ppv .AND. ppv .LE. p_max)) THEN
         ! Select PVRS Riemann solver
         p = ppv
      ELSE
         IF (ppv .LT. p_min) THEN
            ! Select Two-Rarefaction Riemann solver
            pq = (w_l%pressure/w_r%pressure)**G1
            um = (pq*w_l%velocity/c_l + w_r%velocity/c_r + G4*(pq - 1.0))/(pq/c_l + 1.0/c_r)
            ptl = 1.0 + G7*(w_l%velocity - um)/c_l
            ptr = 1.0 + G7*(um - w_r%velocity)/c_r
            p = 0.5*(w_l%pressure*ptl**G3 + w_r%pressure*ptr**G3)
         ELSE
            ! Select Two-Shock Riemann solver with PVRS as estimate
            gel = DSQRT((G5/w_l%density)/(G6*w_l%pressure + ppv))
            ger = DSQRT((G5/w_r%density)/(G6*w_r%pressure + ppv))
            p = (gel*w_l%pressure + ger*w_r%pressure - (w_r%velocity - w_l%velocity))/(gel + ger)
         END IF
      END IF
   end function calculate_pressure

   !> расчет давления на контактном разрыве в первом приближении
   !> TODO: нет размерностей
   function initial_contact_pressure(c_l, w_l, c_r, w_r) result(pressure_cont)
      ! Входные данные
      type(nonconservative_vars_t), intent(in)::w_l !< Вектор потоков "слева"
      type(nonconservative_vars_t), intent(in)::w_r !< Вектор потоков "справа"
      real(kind=real64), intent(in)::c_l
      real(kind=real64), intent(in)::c_r
      ! Выходные данные
      real(kind=real64)::pressure_cont

      pressure_cont = 0.5*(calculate_pressure(c_l, w_l, c_r, w_r, w_l%gamma) + calculate_pressure(c_l, w_l, c_r, w_r, w_r%gamma))

      ! "Звуковой распад"
      ! для произхвольных гамма
      ! pressure_cont = ( WL%pressure * WR%density * cR + WR%pressure * WL%density * cL + ( WL%velocity - WR%velocity ) * WL%density * cL * WR%density * cR )/( WL%density * cL + WR%density * cR )

   end function initial_contact_pressure

   ! расчет функции и производной функции давления
   subroutine calc_pressure_func(pressure_cont, c, w, f, df)
      real(kind=real64), intent(in) :: pressure_cont         ! давление на контактном разрыве
      real(kind=real64), intent(in) :: c                     ! скорость звука
      type(nonconservative_vars_t), intent(in) :: W   ! вектор примитивных переменных
      real(kind=real64), intent(out) :: f                    ! функция давления
      real(kind=real64), intent(out) :: df                   ! производная функции давления
      real(kind=real64) :: p_ratio, G1, G2, G3, G4, A

      p_ratio = (pressure_cont)/(W%pressure)
      G1 = W%gamma + 1.0
      G2 = W%gamma - 1.0
      G3 = 0.5*(W%gamma + 1.0)/W%gamma
      G4 = 0.5*(W%gamma - 1.0)/W%gamma
      if (pressure_cont >= W%pressure) then
         ! ударная волна
         A = G3*p_ratio + G4
         F = (pressure_cont - W%pressure)/(W%density*c*dsqrt(A))
         DF = 0.25*(G1*p_ratio + 3.*W%gamma - 1.)/(W%gamma*W%density*c*A**1.5)
      else
         ! волна разрежения
         A = p_ratio**G4
         F = 2.0*c/G2*(A - 1.0)
         DF = c*A/(W%gamma*(pressure_cont))
      end if

   end subroutine calc_pressure_func

   ! расчет давления и скорости на контактном разрыве с помощью итерационного метода
   subroutine calc_contact_pressure_and_velocity(cL, WL, cR, WR, pressure_cont, velocity_cont)
      implicit none
      real(kind=real64), intent(in) :: cL, cR                         ! скорость звука слева и справа от контактного разрыва
      type(nonconservative_vars_t), intent(in) :: WL, WR       ! вектор примитивных переменных слева и справа от контактного разрыва
      real(kind=real64), intent(out) :: pressure_cont                 ! давление на контактном разрыве
      real(kind=real64), intent(out) :: velocity_cont                 ! скорость на контактном разрыве
      real(kind=real64) :: p_prev
      real(kind=real64) :: fL, dfL, fR, dfR
      integer :: iter_num
      real(kind=real64) :: criteria

      real(kind=real64), parameter :: EPS = 1.d-6
      integer, parameter :: MAX_ITER_NUM = 500

      if (2.*(cL/(WL%gamma - 1.) + cR/(WR%gamma - 1.)) <= (WR%velocity - WL%velocity)) then
         WRite (*, *) " calc_contact_pressure_velocity -> vacuum is generated "
         return
      end if

      iter_num = 0
      p_prev = initial_contact_pressure(cL, WL, cR, WR)

      if (p_prev < 0.) then
         WRite (*, *) " calc_contact_pressure_velocity -> initial pressure guess is negative "
         return
      end if

      do
         call calc_pressure_func(p_prev, cL, WL, fL, dfL)
         call calc_pressure_func(p_prev, cR, WR, fR, dfR)
         pressure_cont = p_prev - (fL + fR + WR%velocity - WL%velocity)/(dfL + dfR)
         criteria = 2.*dabs((pressure_cont - p_prev))/(pressure_cont + p_prev)
         iter_num = iter_num + 1
         if (iter_num > MAX_ITER_NUM) then
            WRite (*, *) " ncalc_contact_pressure_velocity -> number of iterations exceeds the maximum value "
            return
         end if
         p_prev = pressure_cont
         if (criteria <= EPS) then
            exit
         end if
      end do

      velocity_cont = 0.5*(WL%velocity + WR%velocity - fL + fR)

   end subroutine calc_contact_pressure_and_velocity

   ! отбор решения задачи о распаде произвольного разрыва
   function sample_solution(ksi, cL, WL, cR, WR, pressure_cont, velocity_cont) result(W)
      use vectors
      real(kind=real64), intent(in) :: ksi                         ! автомодельная переменная ( x/t )
      real(kind=real64), intent(in) :: pressure_cont               ! давление на контактном разрыве
      real(kind=real64), intent(in) :: velocity_cont               ! скорость  на контактном разрыве
      real(kind=real64), intent(in) :: cL, cR                      ! скорости звука слева и справа от разрыва
      type(nonconservative_vars_t), intent(in) :: WL, WR    ! векторы примитивных переменных слева и справа от разрыва
      type(nonconservative_vars_t) :: W                     ! вектор примитивных переменных
      real(kind=real64) :: G1, G2
      real(kind=real64) :: shL, cmL, stL, aL, sL
      real(kind=real64) :: shR, cmR, stR, aR, sR
      real(kind=real64) :: c

      if (ksi <= velocity_cont) then

         G1 = 0.5*(WL%gamma + 1.)
         G2 = 0.5*(WL%gamma - 1.)

         if (pressure_cont < WL%pressure) then
            shL = WL%velocity - cL
            if (ksi <= shL) then
               W = WL
            else
               cmL = cL + G2*(WL%velocity - velocity_cont)
               stL = velocity_cont - cmL
               if (ksi > stL) then
                  W%density = WL%gamma*pressure_cont/cmL**2.
                  W%velocity = velocity_cont
                  W%pressure = pressure_cont
                  W%gamma = WL%gamma
               else
                  c = (2.*cL + (WL%gamma - 1.)*(WL%velocity - ksi))/(WL%gamma + 1.)
                  w%velocity = ksi + c
                  w%density = WL%density*(c/cL)**(2./(WL%gamma - 1.))
                  w%pressure = w%density*c**2./WL%gamma
                  W%gamma = WL%gamma
               end if
            end if
         else
            aL = dsqrt(WL%density*(G1*pressure_cont + G2*WL%pressure))
            sL = WL%velocity - aL/WL%density
            if (ksi <= sL) then
               W = WL
            else
               w%density = WL%density*aL/(aL - WL%density*(WL%velocity - velocity_cont))
               w%velocity = velocity_cont
               w%pressure = pressure_cont
               W%gamma = WL%gamma
            end if
         end if
      else
         G1 = 0.5*(WR%gamma + 1.)
         G2 = 0.5*(WR%gamma - 1.)

         if (pressure_cont > WR%pressure) then
            aR = dsqrt(WR%density*(G1*pressure_cont + G2*WR%pressure))
            sR = WR%velocity + aR/WR%density
            if (ksi >= sR) then
               W = WR
            else
               w%density = WR%density*aR/(aR + WR%density*(WR%velocity - velocity_cont))
               w%velocity = velocity_cont
               w%pressure = pressure_cont
               W%gamma = WR%gamma
            end if
         else
            shR = WR%velocity + cR
            if (ksi >= shR) then
               w = WR
            else
               cmR = cR - G2*(WR%velocity - velocity_cont)
               stR = velocity_cont + cmR
               if (ksi <= stR) then
                  w%density = WR%gamma*pressure_cont/cmR**2.
                  w%velocity = velocity_cont
                  w%pressure = pressure_cont
                  W%gamma = WR%gamma
               else
                  c = (2.*cR + (WR%gamma - 1.)*(ksi - WR%velocity))/(WR%gamma + 1.)
                  w%velocity = ksi - c
                  w%density = WR%density*(c/cR)**(2./(WR%gamma - 1.))
                  W%pressure = w%density*c**2./WR%gamma
                  W%gamma = WR%gamma
               end if
            end if
         end if
      end if

   end function sample_solution

   ! построение точного решения задачи о распаде произвольного разрыва
   function exact_solution(coordinates, time, cL, WL, cR, WR) result(solution)
      real(kind=real64), intent(in) :: coordinates(:)                    ! массив координат
      real(kind=real64), intent(in) :: time                              ! момент времени, для которого строится точное решение
      real(kind=real64), intent(in) :: cL, cR                            ! скорости звука слева и справа от разрыва
      type(nonconservative_vars_t), intent(in) :: WL, WR          ! векторы примитивных переменных слева и справа от разрыва
      type(nonconservative_vars_t), allocatable :: solution(:)    ! точное решение задачи о распаде произвольного разрыва
      integer :: i
      real(kind=real64) :: pressure_cont, velocity_cont

      allocate (solution(size(coordinates, dim=1)))
      call calc_contact_pressure_and_velocity(cL, WL, cR, WR, pressure_cont, velocity_cont)
      solution = [(sample_solution(coordinates(i)/time, cL, WL, cR, WR, pressure_cont, velocity_cont), &
                   i=1, size(coordinates, dim=1))]

   end function exact_solution

end module riemann_problem
