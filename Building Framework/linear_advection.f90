program linear_advection
  use iso_fortran_env, only: real32, real64, real128
! Define working precision based on compilation flag
#ifdef FP16
  ! Note: REAL16 here refers to 16-bit (half precision) if supported by the compiler
  integer, parameter :: WP = real16
#elif defined(FP32)
  integer, parameter :: WP = real32      ! 32-bit single precision
#elif defined(FP64)
  integer, parameter :: WP = real64      ! 64-bit double precision
#elif defined(FP128)
  integer, parameter :: WP = real128     ! 128-bit quadruple precision
#else
  integer, parameter :: WP = real64      ! default to double precision if not specified
#endif

  implicit none
  ! Simulation parameters
  real(WP), parameter :: a       = 2.0_wp        ! Advection speed
  real(WP), parameter :: L       = 1.0_wp        ! Domain length [0, L)
  integer,  parameter :: N       = 100           ! Number of spatial grid points
  real(WP), parameter :: CFL     = 0.5_wp        ! CFL number for stability
  real(WP), parameter :: dx      = L / real(N, WP)  ! Spatial resolution
  real(WP), parameter :: T_final = 0.5_wp        ! Final simulation time
  real(WP) :: dt                                   
  integer :: nsteps_full, nsteps_total, output_interval
  logical :: has_last_step

  ! Arrays for solution
  real(WP), dimension(:), allocatable :: u, u_new
  integer :: i, step, left, right

  ! Compute time step based on CFL criterion and determine number of steps
  dt = CFL * dx / abs(a)
  nsteps_full = int(T_final / dt)           ! integer number of full dt steps
  if (abs(real(nsteps_full,WP)*dt - T_final) < 1.0e-12_wp) then
    has_last_step = .false.
    nsteps_total = nsteps_full
  else
    has_last_step = .true.
    nsteps_total = nsteps_full + 1         ! one extra step for the last partial step
  end if

  ! Determine output interval (approximately 50 snapshots for animation)
  output_interval = max(1, nsteps_total / 50)
  
  allocate(u(1:N))
  allocate(u_new(1:N))

  ! Initialize the solution: u(x,0) = sin(2πx)
  do i = 1, N
    ! Compute position x_i in [0, L)
    real(WP) :: x
    x = (real(i-1, WP) * dx)
    u(i) = sin(2.0_wp * acos(-1.0_wp) * x)   ! acos(-1) provides π in Fortran
  end do

  ! Output initial condition to CSV (step 0)
  call output_to_csv(u, 0)

  ! Time integration loop
  do step = 1, nsteps_full
     ! Upwind scheme update (first-order explicit upwind)
     real(WP) :: Cplus, Cminus
     Cplus  = max(a, 0.0_wp) * dt / dx      ! coefficient for upwind difference if a > 0
     Cminus = min(a, 0.0_wp) * dt / dx      ! coefficient for upwind difference if a < 0
     
     do i = 1, N
        left  = i - 1; if (i == 1)  left  = N    ! index to the left neighbor (periodic)
        right = i + 1; if (i == N)  right = 1    ! index to the right neighbor (periodic)
        ! Update formula:
        ! If a > 0: u_new(i) = u(i) - a*dt/dx * [u(i) - u(left)]
        ! If a < 0: u_new(i) = u(i) - a*dt/dx * [u(right) - u(i)]
        u_new(i) = u(i) - Cplus * (u(i) - u(left)) - Cminus * (u(right) - u(i))
     end do

     u(:) = u_new(:)  ! advance solution to next time level

     ! Output snapshot at this step if it aligns with the output interval
     if (mod(step, output_interval) == 0 .or. step == nsteps_full) then
        call output_to_csv(u, step)
     end if
  end do

  ! Handle final fractional time step (if T_final is not an integer multiple of dt)
  if (has_last_step) then
     real(WP) :: last_dt, Cplus, Cminus
     last_dt = T_final - real(nsteps_full, WP) * dt   ! remaining time to reach T_final
     Cplus  = max(a, 0.0_wp) * last_dt / dx
     Cminus = min(a, 0.0_wp) * last_dt / dx

     do i = 1, N
        left  = i - 1; if (i == 1)  left = N
        right = i + 1; if (i == N)  right = 1
        u_new(i) = u(i) - Cplus * (u(i) - u(left)) - Cminus * (u(right) - u(i))
     end do

     u(:) = u_new(:)
     step = nsteps_total
     call output_to_csv(u, step)   ! output final state at T_final
  else
     ! If no fractional step, ensure final state is output if not already done
     if (mod(nsteps_total, output_interval) /= 0) then
        call output_to_csv(u, nsteps_total)
     end if
  end if

contains

  !---------------------------------------------------------------
  ! Subroutine to output the solution array `u_array` to a CSV file
  ! The file is named as "advection_<precision>_XXXX.csv", where XXXX is the time-step index.
  !---------------------------------------------------------------
  subroutine output_to_csv(u_array, step_no)
    implicit none
    real(WP), intent(in) :: u_array(:)
    integer,  intent(in) :: step_no
    integer :: j, ios, unit
    character(len=10) :: tag
    character(len=50) :: filename

    ! Determine precision tag for filename
  #ifdef FP16
    tag = "fp16"
  #elif defined(FP32)
    tag = "fp32"
  #elif defined(FP64)
    tag = "fp64"
  #elif defined(FP128)
    tag = "fp128"
  #endif

    write(filename, '(A, A, "_", I4.4, ".csv")') "advection_", tag, step_no
    open(newunit=unit, file=filename, status="replace", action="write", iostat=ios)
    if (ios /= 0) stop "Error opening output file"

    ! Write each grid point's coordinate and solution value as "x,u"
    do j = 1, size(u_array)
       write(unit, '(F8.5, ",", F8.5)') real(j-1, WP)*dx, u_array(j)
    end do

    close(unit)
  end subroutine

end program linear_advection
