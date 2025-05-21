! 1D Heat Equation Solver using Continuous Galerkin FEM
! Equation: du/dt = alpha * d^2u/dx^2
! Solves with linear elements and Crank-Nicolson time stepping

! Compile with different precision:
! gfortran -cpp -DPRECISION_FP16 heat_equation_fem.f90 -o heat_fem_fp16
! gfortran -cpp -DPRECISION_FP32 heat_equation_fem.f90 -o heat_fem_fp32
! gfortran -cpp heat_equation_fem.f90 -o heat_fem_fp64
! gfortran -cpp -DPRECISION_FP128 heat_equation_fem.f90 -o heat_fem_fp128

! Compile with ROOFLINE support:
! gfortran -cpp -DROOFLINE -DPRECISION_FP16 roofline_counter.o heat_equation_fem.f90 -o heat_roof_fp16
! gfortran -cpp -DROOFLINE -DPRECISION_FP32 roofline_counter.o heat_equation_fem.f90 -o heat_roof_fp32
! gfortran -cpp -DROOFLINE roofline_counter.o heat_equation_fem.f90 -o heat_roof_fp64
! gfortran -cpp -DROOFLINE -DPRECISION_FP128 roofline_counter.o heat_equation_fem.f90 -o heat_roof_fp128

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

module energy_mod
  use precision_mod, only: wp
  implicit none
  
contains
  
  
  !-------------------------------------------------------------------------------
  !> Compute M = ∫ u dx via trapezoidal rule
  function mass_trap(u, dx) result(m)
    real(wp), intent(in) :: u(:)   ! Solution vector
    real(wp), intent(in) :: dx     ! Element size
    real(wp)             :: m      ! Resulting mass
    integer              :: n
    
    n = size(u)
    if (n > 1) then
      m = dx * ( 0.5_wp*(u(1) + u(n)) + sum(u(2:n-1)) )
    else
      m = u(1) * dx
    end if
  end function mass_trap
  !-------------------------------------------------------------------------------
  
end module energy_mod

module fem_heat_mod
  use precision_mod
#ifdef ROOFLINE
  use roofline_counter, only: add_flops, add_bytes
#endif
  implicit none
  !— where to put all CSV / error / summary files —
  character(len=*), parameter :: output_dir = 'csv outputs'

  ! Problem parameters
  real(wp), parameter :: alpha = 0.01_wp      ! Diffusion coefficient
  real(wp), parameter :: x_min = 0.0_wp       ! Domain start
  real(wp), parameter :: x_max = 1.0_wp       ! Domain end
  real(wp), parameter :: t_end = 0.5_wp       ! End time
  real(wp), parameter :: pi = 3.14159265358979323846_wp
  
  ! Discretization parameters
  integer, parameter :: n_elements = 100      ! Number of elements
  integer, parameter :: n_nodes = n_elements+1 ! Number of nodes
  real(wp), parameter :: dx = (x_max - x_min) / n_elements
  real(wp), parameter :: dt = 0.5_wp * dx**2 / alpha ! Stable time step based on CFL
  integer, parameter :: n_time_steps = ceiling(t_end / dt)
  
  ! Arrays for solution and mesh
  real(wp), allocatable :: x(:)               ! Node coordinates
  real(wp), allocatable :: u(:), u_new(:)     ! Solution
  real(wp), allocatable :: M(:,:), K(:,:), A(:,:) ! System matrices
  real(wp), allocatable :: b(:), rhs(:)       ! Right-hand side vector
  
contains
  
  ! Initialize mesh and matrices
  subroutine initialize()
    integer :: i, j
    
    ! Allocate arrays
    allocate(x(n_nodes))
    allocate(u(n_nodes), u_new(n_nodes))
    allocate(M(n_nodes,n_nodes), K(n_nodes,n_nodes), A(n_nodes,n_nodes))
    allocate(b(n_nodes), rhs(n_nodes))
    
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
    
    ! Apply quantization for half precision if needed
#if defined(PRECISION_FP16)
    call quantize_to_half(u)
#endif
    
    ! Assemble mass matrix (lumped for simplicity)
    M = 0.0_wp
    do i = 1, n_nodes
      M(i,i) = dx
      if (i > 1) M(i,i) = M(i,i) + dx/2.0_wp
      if (i < n_nodes) M(i,i) = M(i,i) + dx/2.0_wp
#ifdef ROOFLINE
      call add_flops(2)  ! 2 conditional additions
      call add_bytes(wp)  ! Update M(i,i)
#endif
    end do
    
    ! Assemble stiffness matrix
    K = 0.0_wp
    do i = 1, n_elements
      K(i,i) = K(i,i) + 1.0_wp/dx
      K(i,i+1) = K(i,i+1) - 1.0_wp/dx
      K(i+1,i) = K(i+1,i) - 1.0_wp/dx
      K(i+1,i+1) = K(i+1,i+1) + 1.0_wp/dx
#ifdef ROOFLINE
      call add_flops(8)  ! 4 additions/subtractions, 4 divisions
      call add_bytes(4*wp)  ! Update 4 K matrix elements
#endif
    end do
    
    ! Assemble Crank-Nicolson matrix: A = M + 0.5*dt*alpha*K
    A = M + 0.5_wp * dt * alpha * K
#ifdef ROOFLINE
    call add_flops(2*n_nodes*n_nodes + n_nodes*n_nodes)  ! Matrix scaling and addition
    call add_bytes(3*n_nodes*n_nodes*wp)  ! Read M and K, write A
#endif
    
    ! Compute the right-hand side matrix: M - 0.5*dt*alpha*K
    b = matmul(M - 0.5_wp * dt * alpha * K, u)
#ifdef ROOFLINE
    call add_flops(2*n_nodes*n_nodes + n_nodes*n_nodes*n_nodes)  ! Matrix subtraction and matrix-vector multiplication
    call add_bytes((2*n_nodes*n_nodes + n_nodes + n_nodes)*wp)  ! Read matrices and vector, write b
#endif
  end subroutine initialize
  
  ! Solve the system using direct solution (for simplicity)
  subroutine solve_system()
    integer :: i, info
    real(wp), allocatable :: A_copy(:,:), b_copy(:)
    integer, allocatable :: ipiv(:)
    
    allocate(A_copy(n_nodes,n_nodes), b_copy(n_nodes), ipiv(n_nodes))
    
    ! Copy matrices for LAPACK solver
    A_copy = A
    b_copy = b
#ifdef ROOFLINE
    call add_bytes(n_nodes*n_nodes*wp + n_nodes*wp)  ! Copy A and b
#endif
    
    ! Simple Gaussian elimination (could use LAPACK for more robust solution)
    do i = 1, n_nodes-1
      A_copy(i+1:n_nodes,i) = A_copy(i+1:n_nodes,i) / A_copy(i,i)
#ifdef ROOFLINE
      call add_flops(n_nodes-i)  ! Division operations
      call add_bytes((n_nodes-i+1)*wp)  ! Read and write operation results
#endif
      
      A_copy(i+1:n_nodes,i+1:n_nodes) = A_copy(i+1:n_nodes,i+1:n_nodes) - &
                                        matmul(A_copy(i+1:n_nodes,i:i), A_copy(i:i,i+1:n_nodes))
#ifdef ROOFLINE
      call add_flops((n_nodes-i)*(n_nodes-i) + (n_nodes-i)*(n_nodes-i))  ! Matrix multiplication and subtraction
      call add_bytes(3*(n_nodes-i)*(n_nodes-i)*wp)  ! Read two matrices, write result
#endif
      
      b_copy(i+1:n_nodes) = b_copy(i+1:n_nodes) - A_copy(i+1:n_nodes,i) * b_copy(i)
#ifdef ROOFLINE
      call add_flops(2*(n_nodes-i))  ! Multiplication and subtraction
      call add_bytes(3*(n_nodes-i)*wp)  ! Read two vectors, write result
#endif
    end do
    
    ! Back substitution
    u_new(n_nodes) = b_copy(n_nodes) / A_copy(n_nodes,n_nodes)
#ifdef ROOFLINE
    call add_flops(1)  ! Division
    call add_bytes(2*wp + wp)  ! Read two values, write result
#endif
    
    do i = n_nodes-1, 1, -1
      u_new(i) = (b_copy(i) - dot_product(A_copy(i,i+1:n_nodes), u_new(i+1:n_nodes))) / A_copy(i,i)
#ifdef ROOFLINE
      call add_flops(2*(n_nodes-i) + 2)  ! Dot product, subtraction, division
      call add_bytes((2*(n_nodes-i) + 3)*wp)  ! Read vectors for dot product, read/write other values
#endif
    end do
    
    ! Apply quantization for half precision if needed
#if defined(PRECISION_FP16)
    call quantize_to_half(u_new)
#endif
    
    deallocate(A_copy, b_copy, ipiv)
  end subroutine solve_system
  
  ! Apply Dirichlet boundary conditions
  subroutine apply_bc()
    ! Zero Dirichlet BCs at both ends
    u(1) = 0.0_wp
    u(n_nodes) = 0.0_wp
    u_new(1) = 0.0_wp
    u_new(n_nodes) = 0.0_wp
    
    ! Adjust the RHS vector for BCs
    b(1) = 0.0_wp
    b(n_nodes) = 0.0_wp
  end subroutine apply_bc
  
  ! Time-stepping routine
  subroutine time_stepping()
    integer :: i, t
    real(wp) :: current_time
    
    ! Output initial condition
    call output_csv(0, 0.0_wp)
    
    ! Time stepping loop
    do t = 1, n_time_steps
      current_time = t * dt
      
      ! Compute the right-hand side
      b = matmul(M - 0.5_wp * dt * alpha * K, u)
#ifdef ROOFLINE
      call add_flops(2*n_nodes*n_nodes + n_nodes*n_nodes*n_nodes)  ! Matrix subtraction and matrix-vector multiplication
      call add_bytes((2*n_nodes*n_nodes + n_nodes + n_nodes)*wp)  ! Read matrices and vector, write b
#endif
      
      ! Apply boundary conditions
      call apply_bc()
      
      ! Solve the system
      call solve_system()
      
      ! Update solution
      u = u_new
#ifdef ROOFLINE
      call add_bytes(n_nodes*wp)  ! Copy u_new to u
#endif
      
      ! Output solution at specific time steps
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
    real(wp) :: u_exact, error_l1, error_l2, error_linf, diff
    character(len=100) :: filename
    
    ! Calculate errors
    error_l1 = 0.0_wp
    error_l2 = 0.0_wp
    error_linf = 0.0_wp
    
    do i = 1, n_nodes
      ! For heat equation, we can compare against the analytical solution
      ! Using the initial condition with exponential decay
      u_exact = sin(2.0_wp * pi * x(i)) * exp(-4.0_wp * pi * pi * alpha * end_time)
      
      ! Accumulate errors
      diff = abs(u(i) - u_exact)
      error_l1 = error_l1 + diff
      error_l2 = error_l2 + diff**2
      error_linf = max(error_linf, diff)
    end do
    
    ! Normalize L1 and L2 errors
    error_l1 = error_l1 / n_nodes
    error_l2 = sqrt(error_l2 / n_nodes)
    
    ! Write errors to file
    write(filename, '(a,a,a,a)') &
      trim(output_dir)//'/', &
      'heat_', &
      trim(precision_label), &
      '_errors.txt'
    
    open(unit=12, file=filename, status='replace')
    write(12, '(a,G0)') 'L1 Error: ', error_l1
    write(12, '(a,G0)') 'L2 Error: ', error_l2
    write(12, '(a,G0)') 'Linf Error: ', error_linf
    close(12)
    
    ! Display errors to console
    print *, 'Error metrics for precision: ', precision_label
    print '(a,es16.8)', 'L1 Error: ', error_l1
    print '(a,es16.8)', 'L2 Error: ', error_l2
    print '(a,es16.8)', 'Linf Error: ', error_linf
  end subroutine calculate_errors
  
  ! Output solution to CSV
  subroutine output_csv(step, time)
    integer, intent(in) :: step
    real(wp), intent(in) :: time
    character(len=100) :: filename
    integer :: i
    integer :: istat
    
    ! Create filename with precision label
    write(filename,'(a,a,i6.6,a)') &
      trim(output_dir)//'/', &
      'heat_'//trim(precision_label)//'_', &
       step, &
      '.csv'

  
    ! Open file
    open(unit=10, &
       file=filename, &
       status='replace', &
       action='write', &
       iostat=istat)
    if (istat /= 0) then
      print *, 'Error opening file ', filename, ' iostat=', istat
    end if

    
    ! Write header
    write(10, '(a,G0)') 'x,u,time=', time
    
    ! Write data
    do i = 1, n_nodes
      write(10, '(G0,a,G0)') x(i), ',', u(i)
    end do
    
    ! Close file
    close(10)
    
    ! Create a summary file for the current time step
    if (step == n_time_steps) then
      write(filename,'(a)') trim(output_dir)//'/'//'heat_'//trim(precision_label)//'_summary.csv'
      
      open(unit=11, file=filename, status='replace', action='write', iostat=istat)
      if (istat /= 0) print *, 'Error opening summary file ', filename
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
  
end module fem_heat_mod

program heat_equation_fem
  use precision_mod
  use energy_mod, only : mass_trap
  use fem_heat_mod
#ifdef ROOFLINE
  use roofline_counter, only : flop_count, byte_count
#endif
  implicit none

  integer :: istat
  integer :: t, unit_mass
  real(wp) :: M_init, M_now      ! Initial and current mass
  character(len=100) :: filename
#ifdef ROOFLINE
  real :: roof_t_start, roof_t_end
#endif

  !— ensure the CSV outputs folder exists —
  call execute_command_line('mkdir -p "'//trim(output_dir)//'"', exitstat=istat)
  if (istat /= 0) then
    print *, 'Warning: could not create directory '//trim(output_dir)
  end if

#ifdef ROOFLINE
  call cpu_time(roof_t_start)
#endif
  
  ! Display precision being used
  print *, '=================================='
  print *, 'Running Heat Equation FEM solver with precision: ', precision_label
  print *, 'Working precision kind: ', wp
  print *, 'Diffusion coefficient alpha: ', alpha
  print *, 'Time step dt: ', dt
  print *, '=================================='

  ! Initialize the problem
  call initialize()
  
  print *, "Time step size:", dt
  print *, "Number of time steps:", n_time_steps
  
  ! Initialize mass tracking
  M_init = mass_trap(u, dx)
  unit_mass = 44
  write(filename,'(a,a,a,a)') trim(output_dir)//'/', &
                            'heat_', trim(precision_label), '_mass.csv'
  open(unit=unit_mass, file=filename, status='replace', action='write')
  write(unit_mass,'(a)') 'step,time,mass'
  write(unit_mass,'(i0,1x,G0,1x,G0,1x,G0)') 0, 0.0_wp, abs(M_init)
  
  ! Run time stepping with mass monitoring
  ! Output initial condition
  call output_csv(0, 0.0_wp)
  
  ! Time stepping loop with mass monitoring
  do t = 1, n_time_steps
    ! Compute the right-hand side
    b = matmul(M - 0.5_wp * dt * alpha * K, u)
    
    ! Apply boundary conditions
    call apply_bc()
    
    ! Solve the system
    call solve_system()
    
    ! Update solution
    u = u_new
    
    ! Output solution and mass at specific time steps
    if (mod(t, max(1, n_time_steps/100)) == 0 .or. t == n_time_steps) then
      call output_csv(t, t*dt)
      
      ! Calculate and log mass
      M_now = mass_trap(u, dx)
      write(unit_mass,'(i0,1x,es16.8,1x,es16.8,1x,es16.8)') t, t*dt, M_now
      
      ! Print mass conservation to screen
      print '(a,i6,2x,es16.8)', 'Mass at step ', t, M_now
    end if
  end do
  
  ! Final mass calculation
  M_now = mass_trap(u, dx)
  
  ! Calculate and output errors
  call calculate_errors(t_end)
  
  ! Output mass conservation summary
  print *, '----------------------------------------'
  print *, 'Mass check: initial M =', M_init
  print *, '            final   M =', M_now
  print *, 'Detailed trace in ', trim(filename)
  print *, '----------------------------------------'
  
  ! Close mass log file
  close(unit_mass)
  
#ifdef ROOFLINE
  call cpu_time(roof_t_end)
  open(unit=99, file=trim(output_dir)//'/heat_'//trim(precision_label)//'_roofline.csv', &
       status='replace', action='write')
  write(99,'(a)') 'runtime_s,flops,bytes'
  write(99,'(G0,1x,i0,1x,i0)') roof_t_end-roof_t_start, flop_count, byte_count
  close(99)
  print *, 'Roofline data written to ', trim(output_dir)//'/heat_'//trim(precision_label)//'_roofline.csv'
#endif
  
  print *, 'Simulation completed. CSV files have been written.'
  
end program heat_equation_fem

! gfortran -cpp -DPRECISION_FP16 heat_equation_fem.f90 -o heat_fem_fp16
! gfortran -cpp -DPRECISION_FP32 heat_equation_fem.f90 -o heat_fem_fp32
! gfortran -cpp heat_equation_fem.f90 -o heat_fem_fp64
! gfortran -cpp -DPRECISION_FP128 heat_equation_fem.f90 -o heat_fem_fp128

! With roofline instrumentation:
! gfortran -cpp -DROOFLINE -DPRECISION_FP16 roofline_counter.o heat_equation_fem.f90 -o heat_roof_fp16
! gfortran -cpp -DROOFLINE -DPRECISION_FP32 roofline_counter.o heat_equation_fem.f90 -o heat_roof_fp32
! gfortran -cpp -DROOFLINE roofline_counter.o heat_equation_fem.f90 -o heat_roof_fp64
! gfortran -cpp -DROOFLINE -DPRECISION_FP128 roofline_counter.o heat_equation_fem.f90 -o heat_roof_fp128

! ./heat_roof_fp16
! ./heat_roof_fp32
! ./heat_roof_fp64
! ./heat_roof_fp128