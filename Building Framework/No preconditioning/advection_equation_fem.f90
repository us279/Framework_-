! 1D Advection Equation Solver using stabilized FEM
! Equation: du/dt + C * du/dx = 0
! Solves with linear elements and SUPG stabilization

! Compile with:
! gfortran -cpp -DPRECISION_FP16 advection_equation_fem.f90 -o advection_fem_fp16
! gfortran -cpp -DPRECISION_FP32 advection_equation_fem.f90 -o advection_fem_fp32
! gfortran -cpp advection_equation_fem.f90 -o advection_fem_fp64
! gfortran -cpp -DPRECISION_FP128 advection_equation_fem.f90 -o advection_fem_fp128

! Compile with ROOFLINE support:
! gfortran -cpp -DROOFLINE -DPRECISION_FP16 roofline_counter.o advection_equation_fem.f90 -o advection_roof_fp16
! gfortran -cpp -DROOFLINE -DPRECISION_FP32 roofline_counter.o advection_equation_fem.f90 -o advection_roof_fp32
! gfortran -cpp -DROOFLINE roofline_counter.o advection_equation_fem.f90 -o advection_roof_fp64
! gfortran -cpp -DROOFLINE -DPRECISION_FP128 roofline_counter.o advection_equation_fem.f90 -o advection_roof_fp128

module precision_mod
  implicit none
  
#if defined(PRECISION_FP16)
  integer, parameter :: wp = selected_real_kind(3, 4)    ! Half precision (FP16)
  character(len=4), parameter :: precision_label = "FP16"
#elif defined(PRECISION_FP32)
  integer, parameter :: wp = selected_real_kind(6, 37)   ! Single precision (FP32)
  character(len=4), parameter :: precision_label = "FP32"
#elif defined(PRECISION_FP128)
  integer, parameter :: wp = selected_real_kind(33, 4931)! Quad precision (FP128)
  character(len=5), parameter :: precision_label = "FP128"
#else
  integer, parameter :: wp = selected_real_kind(15, 307) ! Double precision (FP64)
  character(len=4), parameter :: precision_label = "FP64"
#endif

end module precision_mod

module entropy_mod
  use precision_mod, only: wp
  implicit none
  
contains
  
  ! Calculate quadratic entropy (energy) of solution
  function entropy_energy(u, dx) result(e)
    real(wp), intent(in) :: u(:)    ! Solution vector
    real(wp), intent(in) :: dx      ! Element size
    real(wp)             :: e       ! Resulting entropy
    
    ! Quadratic entropy: E(t) = 0.5 * ∫ u(x,t)² dx
    e = 0.5_wp * sum(u**2) * dx
  end function entropy_energy
  
  ! Trapezoidal‐rule version of quadratic entropy:
  function entropy_energy_trap(u, dx) result(e)
    use precision_mod, only: wp
    implicit none
    real(wp), intent(in) :: u(:)   ! nodal values
    real(wp), intent(in) :: dx     ! mesh spacing
    real(wp)             :: e
    integer              :: n

    n = size(u)
    if (n < 2) then
      e = 0.0_wp
    else
      e = 0.5_wp*dx * ( 0.5_wp*(u(1)**2 + u(n)**2) + sum(u(2:n-1)**2) )
    end if
  end function entropy_energy_trap
  
end module entropy_mod

module fem_advection_mod
  use precision_mod
#ifdef ROOFLINE
  use roofline_counter, only: add_flops, add_bytes
#endif
  implicit none

  !— where to put all CSV / error / summary files —
  character(len=*), parameter :: output_dir = 'csv outputs'


  ! Problem parameters
  real(wp), parameter :: C = 2.0_wp        ! Advection velocity
  real(wp), parameter :: x_min = 0.0_wp    ! Domain start
  real(wp), parameter :: x_max = 1.0_wp    ! Domain end
  real(wp), parameter :: t_end = 0.5_wp    ! End time
  real(wp), parameter :: pi = 3.14159265358979323846_wp
  
  ! Discretization parameters
  integer, parameter :: n_elements = 100    ! Number of elements
  integer, parameter :: n_nodes = n_elements + 1 ! Number of nodes
  real(wp), parameter :: dx = (x_max - x_min) / n_elements
  
  ! Arrays for solution and mesh
  real(wp), allocatable :: x(:)            ! Node coordinates
  real(wp), allocatable :: u(:), u_old(:)  ! Solution vectors
  real(wp), allocatable :: M(:,:)          ! Mass matrix
  real(wp), allocatable :: K(:,:)          ! Advection matrix (includes SUPG)
  real(wp), allocatable :: rhs(:)          ! Right-hand side vector
  real(wp), allocatable :: mlump(:)
  
  ! Time stepping parameters
  real(wp), parameter :: cfl = 0.1_wp      ! CFL number (kept small for stability)
  real(wp) :: dt                          ! Time step
  integer :: n_time_steps                  ! Number of time steps
  
contains
  
  ! Initialize mesh, matrices, and initial condition
  subroutine initialize()
    integer :: i, j
    real(wp) :: tau_supg                  ! SUPG stabilization parameter
    
    ! Allocate arrays
    allocate(x(n_nodes))
    allocate(u(n_nodes), u_old(n_nodes))
    allocate(M(n_nodes, n_nodes))
    allocate(K(n_nodes, n_nodes))
    allocate(rhs(n_nodes))
    allocate(mlump(n_nodes))
    
    ! Generate mesh
    do i = 1, n_nodes
      x(i) = x_min + (i-1) * dx
    end do
    
    ! Set initial condition: u(x,0) = sin(2*pi*x)
    do i = 1, n_nodes
      u(i) = sin(2.0_wp * pi * x(i))
#ifdef ROOFLINE
      call add_flops(20)  ! sin and multiply are complex operations
      call add_bytes(wp)  ! Write to u(i)
#endif
    end do

#if defined(PRECISION_FP16)
  call quantize_to_half(u)
#endif

    
    ! Calculate time step based on CFL
    dt = cfl * dx / abs(C)
    n_time_steps = ceiling(t_end / dt)
    
    print *, "Time step size:", dt
    print *, "Number of time steps:", n_time_steps
    
    ! Assemble mass matrix (consistent mass matrix)
    M = 0.0_wp
    do i = 1, n_elements
      ! Local mass matrix for element i (linear elements)
      M(i,i) = M(i,i) + dx / 3.0_wp
      M(i,i+1) = M(i,i+1) + dx / 6.0_wp
      M(i+1,i) = M(i+1,i) + dx / 6.0_wp
      M(i+1,i+1) = M(i+1,i+1) + dx / 3.0_wp
#ifdef ROOFLINE
      call add_flops(4)  ! 4 additions
      call add_bytes(4*wp)  ! 4 writes to M
#endif
    end do

    ! Lumped‑mass vector (never changes, so compute once) -- FIX
    do i = 1, n_nodes
      mlump(i) = sum(M(i,:))
#ifdef ROOFLINE
      call add_flops(n_nodes)  ! n_nodes additions for sum
      call add_bytes((n_nodes+1)*wp)  ! n_nodes reads, 1 write
#endif
    end do
    
    ! Assemble advection matrix with SUPG stabilization
    K = 0.0_wp
    do i = 1, n_elements
      ! Standard Galerkin contribution for advection (linear elements)
      K(i,i) = K(i,i) + C * 0.5_wp
      K(i,i+1) = K(i,i+1) - C * 0.5_wp
      K(i+1,i) = K(i+1,i) + C * 0.5_wp
      K(i+1,i+1) = K(i+1,i+1) - C * 0.5_wp
      
      ! Add SUPG stabilization
      tau_supg = cfl * dx / (2.0_wp * abs(C))  ! Stabilization parameter
      
      K(i,i) = K(i,i) + tau_supg * C * C / dx
      K(i,i+1) = K(i,i+1) - tau_supg * C * C / dx
      K(i+1,i) = K(i+1,i) - tau_supg * C * C / dx
      K(i+1,i+1) = K(i+1,i+1) + tau_supg * C * C / dx
#ifdef ROOFLINE
      call add_flops(12)  ! 12 operations for SUPG assembly
      call add_bytes(8*wp)  ! 4 reads, 4 writes to K
#endif
    end do
  end subroutine initialize
  
  ! Apply periodic boundary conditions
  subroutine apply_periodic_bc()
    ! For periodic BC, the solution at the first and last nodes should be equal
    ! We modify the system to enforce u(1) = u(n_nodes)
    
    ! Add the last row to the first row
    K(1,:) = K(1,:) + K(n_nodes,:)
    M(1,:) = M(1,:) + M(n_nodes,:)
    
    ! Set the last row to enforce u(n_nodes) = u(1)
    K(n_nodes,:) = 0.0_wp
    M(n_nodes,:) = 0.0_wp
    K(n_nodes,1) = 1.0_wp
    K(n_nodes,n_nodes) = -1.0_wp
    
    ! Update the RHS for the last row
    rhs(n_nodes) = 0.0_wp
  end subroutine apply_periodic_bc
  
  ! Solve one time step using explicit Euler method
  subroutine solve_time_step()
    integer :: i
    real(wp), allocatable :: k1(:), k2(:), k3(:), k4(:)
    allocate(k1(n_nodes), k2(n_nodes), k3(n_nodes), k4(n_nodes))

    ! rhs = − M⁻¹ K u   (element‑wise division by mlump)
    k1 = -( matmul(K, u) ) / mlump
#ifdef ROOFLINE
    call add_flops(n_nodes * n_nodes + n_nodes) ! matmul + division by mlump
    call add_bytes(3 * n_nodes * wp) ! read K, u, write k1
#endif

    k2 = -( matmul(K, u + 0.5_wp*dt*k1) ) / mlump
#ifdef ROOFLINE
    call add_flops(2 * n_nodes + n_nodes * n_nodes + n_nodes) ! u + k1*dt/2, matmul, division
    call add_bytes(4 * n_nodes * wp) ! read K, u, k1, write k2
#endif

    k3 = -( matmul(K, u + 0.5_wp*dt*k2) ) / mlump
#ifdef ROOFLINE
    call add_flops(2 * n_nodes + n_nodes * n_nodes + n_nodes) ! u + k2*dt/2, matmul, division
    call add_bytes(4 * n_nodes * wp) ! read K, u, k2, write k3
#endif

    k4 = -( matmul(K, u + dt*k3) ) / mlump
#ifdef ROOFLINE
    call add_flops(2 * n_nodes + n_nodes * n_nodes + n_nodes) ! u + k3*dt, matmul, division
    call add_bytes(4 * n_nodes * wp) ! read K, u, k3, write k4
#endif

    u = u + dt*(k1 + 2.0_wp*k2 + 2.0_wp*k3 + k4)/6.0_wp
#ifdef ROOFLINE
    call add_flops(5 * n_nodes) ! weighted sum of k terms
    call add_bytes(5 * n_nodes * wp) ! read k1,k2,k3,k4,u, write u
#endif

    u(n_nodes) = u(1)
  end subroutine solve_time_step
  
  ! Run the simulation with explicit time stepping
  subroutine time_stepping()
    integer :: t
    real(wp) :: current_time
    
    ! Output initial condition
    call output_csv(0, 0.0_wp)
    

    ! Emulate half precision by rounding the new solution
#if defined(PRECISION_FP16)
    call quantize_to_half(u)
#endif

    ! Time stepping loop
    do t = 1, n_time_steps
      current_time = t * dt
      
      ! Solve one time step
      call solve_time_step()
      
      ! Output solution at specified intervals
      if (mod(t, max(1, n_time_steps/10)) == 0 .or. t == n_time_steps) then
        call output_csv(t, current_time)
      end if
    end do
    
    ! Calculate and output errors
    call calculate_errors(t_end)
  end subroutine time_stepping
  
  ! Calculate errors relative to exact solution
  subroutine calculate_errors(end_time)
    real(wp), intent(in) :: end_time
    integer :: i
    real(wp) :: x_shifted, u_exact, error_l1, error_l2, error_linf
    character(len=100) :: filename
    
    ! Calculate errors
    error_l1 = 0.0_wp
    error_l2 = 0.0_wp
    error_linf = 0.0_wp
    
    do i = 1, n_nodes
      ! For advection problem, the exact solution at time t is u(x-C*t, 0) with periodic BC
      x_shifted = x(i) - C * end_time
      
      ! Apply periodic boundary conditions
      do while (x_shifted < x_min)
        x_shifted = x_shifted + (x_max - x_min)
      end do
      do while (x_shifted >= x_max)
        x_shifted = x_shifted - (x_max - x_min)
      end do
      
      ! Exact solution at this point (using our initial condition)
      u_exact = sin(2.0_wp * pi * x_shifted)
      
      ! Accumulate errors
      error_l1 = error_l1 + abs(u(i) - u_exact)
      error_l2 = error_l2 + (u(i) - u_exact)**2
      error_linf = max(error_linf, abs(u(i) - u_exact))
    end do
    
    ! Normalize L1 and L2 errors
    error_l1 = error_l1 / n_nodes
    error_l2 = sqrt(error_l2 / n_nodes)
    
    ! Write errors to file
    write(filename, '(a,a,a,a)') &
     trim(output_dir)//'/',     &
     'advection_',              &
     trim(precision_label),     &
     '_errors.txt'


    open(unit=12, file=filename, status='replace')
    write(12, '(a,G0)')   'L1 Error: ',    error_l1
    write(12, '(a,G0)')   'L2 Error: ',    error_l2
    write(12, '(a,G0)')   'Linf Error: ',  error_linf
    close(12)
    
    ! Display errors to console
    print *, 'Error metrics for precision: ', precision_label
    print '(a,es16.8)', 'L1 Error: ', error_l1
    print '(a,es16.8)', 'L2 Error: ', error_l2
    print '(a,es16.8)', 'Linf Error: ', error_linf
  end subroutine calculate_errors
  
  ! Output solution to CSV file
  subroutine output_csv(step, time)
    integer, intent(in) :: step
    real(wp), intent(in) :: time
    character(len=100) :: filename
    integer :: i
    
    ! Create filename with precision label
    write(filename, '(a,a,a,a,i6.6,a)') &
     trim(output_dir)//'/',     & ! <- new prefix
     'advection_',              &
     trim(precision_label),     &
     '_',                       &
     step,                     &
     '.csv'

    
    ! Open file
    open(unit=10, file=filename, status='replace')
    
    ! Write header
    write(10, '(a,G0)')    'x,u,time=', time
    
    ! Write data
    do i = 1, n_nodes
      write(10, '(G0,a,G0)') x(i), ',', u(i)
    end do
    
    ! Close file
    close(10)
    
    ! Create a summary file for the current time step
    if (step == n_time_steps) then
      write(filename, '(a,a,a,a)') &
      trim(output_dir)//'/',     &
      'advection_',              &
      trim(precision_label),     &
      '_summary.csv'

      open(unit=11, file=filename, status='replace')
      write(11, '(a)') 'x,u'
      do i = 1, n_nodes
        write(11, '(G0,a,G0)') x(i), ',', u(i)
      end do
      close(11)
    end if
  end subroutine output_csv

  !---------------------------------------------------------------
  !> Emulate IEEE-754 binary16 by rounding every entry to 10-bit mantissa
  subroutine quantize_to_half(v)
    use precision_mod, only: wp
    implicit none
    real(wp), intent(inout) :: v(:)
    integer :: i
    real(wp), parameter :: scale = 2.0_wp**10  ! 10 bits of mantissa
    do i = 1, size(v)
      v(i) = nint(v(i) * scale) / scale
    end do
  end subroutine quantize_to_half
  !---------------------------------------------------------------


  
end module fem_advection_mod

program advection_equation_fem
  use precision_mod
  use entropy_mod, only: entropy_energy, entropy_energy_trap
  use fem_advection_mod
#ifdef ROOFLINE
  use roofline_counter, only : flop_count, byte_count
#endif
  implicit none


  integer :: istat
  integer :: t, unit_ent
  real(wp) :: E0, E, current_time   ! Initial and running entropy
  character(len = 100) :: filename
#ifdef ROOFLINE
  real     :: roof_t_start, roof_t_end
#endif

  !— Create output directory if it doesn't exist —
  call execute_command_line('mkdir -p "'//trim(output_dir)//'"', exitstat=istat)
  if (istat /= 0) then
    print *, 'Warning: could not create directory '//trim(output_dir)
  end if
  
#ifdef ROOFLINE
  call cpu_time(roof_t_start)
#endif
  
  ! Display precision being used
  print *, '=================================='
  print *, 'Running Advection Equation FEM solver with precision: ', precision_label
  print *, 'Working precision kind: ', wp
  print *, 'Advection velocity C: ', C
  print *, 'CFL number: ', cfl
  print *, '=================================='

  ! Initialize the problem
  call initialize()
  
  ! Initialize entropy tracking
  E0 = entropy_energy_trap(u, dx)
  unit_ent = 42
  write(filename, '(a,a,a,a)') &
      trim(output_dir)//'/', &
      'advection_', &
      trim(precision_label), &
      '_entropy.csv'
  open(unit=unit_ent, file=filename, status='replace', action='write')
  write(unit_ent, '(a)') 'step,time,E'
  write(unit_ent, '(i0,1x,G0,1x,G0)')         0, 0.0_wp, abs(E0)
  
  ! Custom time stepping loop with entropy monitoring
  ! Output initial condition
  call output_csv(0, 0.0_wp)
  
  ! Emulate half precision by rounding the new solution
#if defined(PRECISION_FP16)
  call quantize_to_half(u)
#endif

  ! Time stepping loop
  do t = 1, n_time_steps
    current_time = t * dt
    
    ! Solve one time step
    call solve_time_step()
    
    ! Output solution and entropy at specified intervals
    if (mod(t, max(1, n_time_steps/100)) == 0 .or. t == n_time_steps) then
      call output_csv(t, current_time)
      
      ! Calculate and log entropy
      E = entropy_energy_trap(u, dx)
      
      
      write(unit_ent, '(i0,1x,G0,1x,G0)') &
         t, current_time, E
      
      
    end if
  end do
  
  ! Calculate and output errors
  call calculate_errors(t_end)
  
  ! Final entropy calculation
  E = entropy_energy_trap(u, dx)
  
  ! Output entropy conservation summary
  print *, '----------------------------------------'
  print *, 'Entropy check: initial E =', E0
  print *, '               final   E =', E
  print *, 'Detailed trace in ', trim(output_dir)//'/'//&
           'advection_'//trim(precision_label)//'_entropy.csv'
  print *, '----------------------------------------'
  
  ! Close entropy log file
  close(unit_ent)
  
#ifdef ROOFLINE
  call cpu_time(roof_t_end)
  open(unit=99, file=trim(output_dir)//'/advection_'//trim(precision_label)//'_roofline.csv', &
       status='replace', action='write')
  write(99,'(a)') 'runtime_s,flops,bytes'
  write(99,'(G0,1x,i0,1x,i0)') roof_t_end - roof_t_start, flop_count, byte_count
  close(99)
  print *, 'Roofline data written to ', trim(output_dir)//'/advection_'//trim(precision_label)//'_roofline.csv'
#endif
  
  print *, 'Simulation completed. CSV files have been written.'
  
end program advection_equation_fem



! gfortran -cpp -DPRECISION_FP16 advection_equation_fem.f90 -o advection_fem_fp16
! gfortran -cpp -DPRECISION_FP32 advection_equation_fem.f90 -o advection_fem_fp32
! gfortran -cpp advection_equation_fem.f90 -o advection_fem_fp64
! gfortran -cpp -DPRECISION_FP128 advection_equation_fem.f90 -o advection_fem_fp128

! With roofline instrumentation:
! gfortran -cpp -DROOFLINE -DPRECISION_FP16 roofline_counter.o advection_equation_fem.f90 -o advection_roof_fp16
! gfortran -cpp -DROOFLINE -DPRECISION_FP32 roofline_counter.o advection_equation_fem.f90 -o advection_roof_fp32
! gfortran -cpp -DROOFLINE roofline_counter.o advection_equation_fem.f90 -o advection_roof_fp64
! gfortran -cpp -DROOFLINE -DPRECISION_FP128 roofline_counter.o advection_equation_fem.f90 -o advection_roof_fp128

! ./advection_roof_fp16
! ./advection_roof_fp32
! ./advection_roof_fp64
! ./advection_roof_fp128