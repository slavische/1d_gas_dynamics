module riemann_problem
   use vectors
   implicit none

contains

   ! расчет давления на контактном разрыве в первом приближении
   function initial_contact_pressure(cL, WL, cR, WR) result(pressure_cont)
      type (vector_nonconservative_vars), intent(in) :: WL, WR
      double precision, intent(in) :: cL, cR
      double precision :: pressure_cont
      double precision, parameter :: QUSER = 2.0
      double precision :: GAMMA, G1, G2, G3, G4, G5, G6, G7, G8
      double precision :: CUP, PPV, PMIN, PMAX, QMAX, PQ, UM, PTL, PTR, GEL,GER

      ! Eleuterio F. Toro Riemann Solvers and Numerical Methods for Fluid Dynamics
      ! текущая реализация только для постоянной гаммы
      ! Compute gamma related constants
      ! TODO: переработать для произвольных гамма
      GAMMA = WL%gamma
      G1 = (GAMMA - 1.0)/(2.0*GAMMA)
      G2 = (GAMMA + 1.0)/(2.0*GAMMA)
      G3 = 2.0*GAMMA/(GAMMA - 1.0)
      G4 = 2.0/(GAMMA - 1.0)
      G5 = 2.0/(GAMMA + 1.0)
      G6 = (GAMMA - 1.0)/(GAMMA + 1.0)
      G7 = (GAMMA - 1.0)/2.0
      G8 = GAMMA - 1.0
      ! Compute guess pressure from PVRS Riemann solver
      CUP  = 0.25 * ( WL%density + WR%density ) * ( cL + cR )
      PPV  = 0.5 * ( WL%pressure + WR%pressure ) + 0.5 * ( WL%velocity - WR%velocity ) * CUP
      PPV  = max( 0.0, PPV )
      PMIN = min( WL%pressure,  WR%pressure )
      PMAX = max( WL%pressure,  WR%pressure )
      QMAX = PMAX / PMIN

      IF(QMAX.LE.QUSER.AND.(PMIN.LE.PPV.AND.PPV.LE.PMAX))THEN
         ! Select PVRS Riemann solver
         pressure_cont = PPV
      ELSE
         IF(PPV.LT.PMIN)THEN
            ! Select Two-Rarefaction Riemann solver
            PQ  = ( WL%pressure / WR%pressure ) ** G1
            UM  = ( PQ * WL%velocity / cL + WR%velocity / cR + G4 * (PQ - 1.0) ) / ( PQ / cL + 1.0 / cR )
            PTL = 1.0 + G7 * ( WL%velocity - UM ) / cL
            PTR = 1.0 + G7 * ( UM - WR%velocity ) / cR
            pressure_cont = 0.5 * ( WL%pressure * PTL ** G3 + WR%pressure * PTR ** G3)
         ELSE
            ! Select Two-Shock Riemann solver with PVRS as estimate
            GEL = DSQRT( ( G5 / WL%density ) / ( G6 * WL%pressure + PPV ) )
            GER = DSQRT( ( G5 / WR%density ) / ( G6 * WR%pressure + PPV ) )
            pressure_cont = ( GEL * WL%pressure + GER * WR%pressure - ( WR%velocity - WL%velocity ) ) / ( GEL + GER )
         ENDIF
      ENDIF

      ! "Звуковой распад" 
      ! для произхвольных гамма
      ! pressure_cont = ( WL%pressure * WR%density * cR + WR%pressure * WL%density * cL + ( WL%velocity - WR%velocity ) * WL%density * cL * WR%density * cR )/( WL%density * cL + WR%density * cR )

   end function initial_contact_pressure

   ! расчет функции и производной функции давления
   subroutine calc_pressure_func(pressure_cont, c, w, f, df)
      double precision, intent(in) :: pressure_cont         ! давление на контактном разрыве
      double precision, intent(in) :: c                     ! скорость звука
      type (vector_nonconservative_vars), intent(in) :: W   ! вектор примитивных переменных
      double precision, intent(out) :: f                    ! функция давления
      double precision, intent(out) :: df                   ! производная функции давления
      double precision :: p_ratio, G1, G2, G3, G4, A

      p_ratio = ( pressure_cont )/( W%pressure  )
      G1 = W%gamma + 1.0
      G2 = W%gamma - 1.0
      G3 = 0.5*(W%gamma + 1.0)/W%gamma
      G4 = 0.5*(W%gamma - 1.0)/W%gamma
      if ( pressure_cont >= W%pressure ) then
         ! ударная волна
         A = G3*p_ratio + G4
         F = (pressure_cont - W%pressure)/(W%density*c*dsqrt(A))
         DF = 0.25*( G1*p_ratio + 3.*W%gamma - 1. )/( W%gamma*W%density*c*A**1.5)
      else
         ! волна разрежения
         A = p_ratio**G4
         F = 2.0*c/G2*(A - 1.0)
         DF = c*A/(W%gamma*(pressure_cont ))
      end if

   end subroutine calc_pressure_func

   ! расчет давления и скорости на контактном разрыве с помощью итерационного метода
   subroutine calc_contact_pressure_and_velocity( cL, WL, cR, WR, pressure_cont, velocity_cont )
      implicit none
      double precision, intent(in) :: cL, cR                         ! скорость звука слева и справа от контактного разрыва
      type (vector_nonconservative_vars), intent(in) :: WL, WR       ! вектор примитивных переменных слева и справа от контактного разрыва
      double precision, intent(out) :: pressure_cont                 ! давление на контактном разрыве
      double precision, intent(out) :: velocity_cont                 ! скорость на контактном разрыве
      double precision :: p_prev
      double precision :: fL, dfL, fR, dfR
      integer :: iter_num
      double precision :: criteria

      double precision, parameter :: EPS = 1.d-6
      integer, parameter :: MAX_ITER_NUM = 500

      if ( 2. * ( cL / ( WL%gamma - 1.) + cR / ( WR%gamma - 1. ) ) <= (WR%velocity - WL%velocity) ) then
         WRite(*,*) " calc_contact_pressure_velocity -> vacuum is generated "
         return
      end if

      iter_num = 0
      p_prev = initial_contact_pressure(cL, WL, cR, WR)

      if ( p_prev < 0. ) then
         WRite(*,*) " calc_contact_pressure_velocity -> initial pressure guess is negative "
         return
      end if

      do
         call calc_pressure_func( p_prev, cL, WL, fL, dfL )
         call calc_pressure_func( p_prev, cR, WR, fR, dfR )
         pressure_cont = p_prev - ( fL + fR + WR%velocity - WL%velocity ) / ( dfL + dfR )
         criteria = 2. * dabs( ( pressure_cont - p_prev ) )/( pressure_cont + p_prev )
         iter_num = iter_num + 1
         if ( iter_num > MAX_ITER_NUM )then
            WRite(*,*) " ncalc_contact_pressure_velocity -> number of iterations exceeds the maximum value "
            return
         end if
         p_prev = pressure_cont
         if ( criteria <= EPS ) then
            exit
         end if
      end do

      velocity_cont = 0.5 * ( WL%velocity + WR%velocity - fL + fR)

   end subroutine calc_contact_pressure_and_velocity

   ! отбор решения задачи о распаде произвольного разрыва
   function sample_solution(ksi, cL, WL, cR, WR, pressure_cont, velocity_cont) result(W)
      use vectors
      double precision, intent(in) :: ksi                         ! автомодельная переменная ( x/t )
      double precision, intent(in) :: pressure_cont               ! давление на контактном разрыве
      double precision, intent(in) :: velocity_cont               ! скорость  на контактном разрыве
      double precision, intent(in) :: cL, cR                      ! скорости звука слева и справа от разрыва
      type (vector_nonconservative_vars), intent(in) :: WL, WR    ! векторы примитивных переменных слева и справа от разрыва
      type (vector_nonconservative_vars) :: W                     ! вектор примитивных переменных
      double precision :: G1, G2
      double precision :: shL, cmL, stL, aL, sL
      double precision :: shR, cmR, stR, aR, sR
      double precision :: c


      if ( ksi <= velocity_cont ) then

         G1 = 0.5*(WL%gamma + 1.)
         G2 = 0.5*(WL%gamma - 1.)

         if ( pressure_cont < WL%pressure) then
            shL = WL%velocity - cL
            if ( ksi <= shL ) then
               W = WL
            else
               cmL = cL + G2 * ( WL%velocity - velocity_cont )
               stL = velocity_cont - cmL
               if ( ksi > stL ) then
                  W%density = WL%gamma * pressure_cont / cmL ** 2.
                  W%velocity = velocity_cont
                  W%pressure = pressure_cont
                  W%gamma = WL%gamma
               else
                  c = ( 2. * cL + ( WL%gamma - 1. ) * ( WL%velocity - ksi ) ) / ( WL%gamma + 1. )
                  w%velocity = ksi + c
                  w%density = WL%density * ( c / cL ) ** ( 2. / ( WL%gamma - 1. ) )
                  w%pressure = w%density *c ** 2. / WL%gamma
                  W%gamma = WL%gamma
               end if
            end if
         else
            aL = dsqrt( WL%density * ( G1 * pressure_cont + G2 * WL%pressure ))
            sL = WL%velocity - aL / WL%density
            if ( ksi <= sL ) then
               W = WL
            else
               w%density = WL%density * aL / ( aL - WL%density * ( WL%velocity - velocity_cont ) )
               w%velocity = velocity_cont
               w%pressure = pressure_cont
               W%gamma = WL%gamma
            end if
         end if
      else
         G1 = 0.5*(WR%gamma + 1.)
         G2 = 0.5*(WR%gamma - 1.)

         if ( pressure_cont > WR%pressure ) then
            aR = dsqrt( WR%density * ( G1 * pressure_cont + G2 * WR%pressure ))
            sR = WR%velocity + aR / WR%density
            if ( ksi >= sR ) then
               W = WR
            else
               w%density = WR%density * aR / ( aR + WR%density * ( WR%velocity - velocity_cont ) )
               w%velocity = velocity_cont
               w%pressure = pressure_cont
               W%gamma = WR%gamma
            end if
         else
            shR = WR%velocity + cR
            if ( ksi >= shR ) then
               w = WR
            else
               cmR = cR - G2 * ( WR%velocity - velocity_cont )
               stR = velocity_cont + cmR
               if ( ksi <= stR ) then
                  w%density = WR%gamma * pressure_cont / cmR ** 2.
                  w%velocity = velocity_cont
                  w%pressure = pressure_cont
                  W%gamma = WR%gamma
               else
                  c = ( 2. * cR + ( WR%gamma - 1. )*( ksi - WR%velocity )) / ( WR%gamma + 1. )
                  w%velocity = ksi - c
                  w%density = WR%density*( c / cR ) ** ( 2. / ( WR%gamma - 1. ) )
                  W%pressure = w%density * c ** 2. / WR%gamma
                  W%gamma = WR%gamma
               end if
            end if
         end if
      end if

   end function sample_solution

   ! построение точного решения задачи о распаде произвольного разрыва
   function exact_solution(coordinates, time, cL, WL, cR, WR) result(solution)
      double precision, intent(in) :: coordinates(:)                    ! массив координат
      double precision, intent(in) :: time                              ! момент времени, для которого строится точное решение
      double precision, intent(in) :: cL, cR                            ! скорости звука слева и справа от разрыва
      type (vector_nonconservative_vars), intent(in) :: WL, WR          ! векторы примитивных переменных слева и справа от разрыва
      type (vector_nonconservative_vars), allocatable :: solution(:)    ! точное решение задачи о распаде произвольного разрыва
      integer :: i
      double precision :: pressure_cont, velocity_cont

      allocate( solution( size( coordinates, dim = 1 ) ))
      call calc_contact_pressure_and_velocity(cL, WL, cR, WR, pressure_cont, velocity_cont)
      solution = [ ( sample_solution( coordinates(i)/time, cL, WL, cR, WR, pressure_cont, velocity_cont ), &
         i = 1, size( coordinates, dim = 1 ) ) ]

   end function exact_solution

end module riemann_problem
