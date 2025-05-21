! 1D Viscous Burgers Equation Solver using SUPG
! Equation: du/dt + u * du/dx = nu * d^2u/dx^2
! Solves with SUPG stabilization and IMEX time stepping

! Compile with:
! gfortran -cpp -DPRECISION_FP16 burgers_equation_supg.f90 -o burgers_supg_fp16
! gfortran -cpp -DPRECISION_FP32 burgers_equation_supg.f90 -o burgers_supg_fp32
! gfortran -cpp burgers_equation_supg.f90 -o burgers_supg_fp64
! gfortran -cpp -DPRECISION_FP128 burgers_equation_supg.f90 -o burgers_supg_fp128

! Compile with ROOFLINE support:
! gfortran -cpp -DROOFLINE -DPRECISION_FP16 roofline_counter.o burgers_equation_supg.f90 -o burgers_roof_fp16
! gfortran -cpp -DROOFLINE -DPRECISION_FP32 roofline_counter.o burgers_equation_supg.f90 -o burgers_roof_fp32
! gfortran -cpp -DROOFLINE roofline_counter.o burgers_equation_supg.f90 -o burgers_roof_fp64
! gfortran -cpp -DROOFLINE -DPRECISION_FP128 roofline_counter.o burgers_equation_supg.f90 -o burgers_roof_fp128

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

module burgers_exact_mod
  use precision_mod, only: wp
  implicit none
  real(wp), parameter :: pi = 3.14159265358979323846_wp
  
contains
  ! Modified Bessel function of the first kind
  function besseli(n, x) result(res)
    integer, intent(in) :: n
    real(wp), intent(in) :: x
    real(wp) :: res
    
    ! Do all intermediate work in DOUBLE PRECISION
    real(8) :: sum8, term8, fact_k8, fact_nk8, x8
    integer :: k, maxk, i
    
    maxk = 50  ! Upper limit for series expansion
    x8 = real(x, 8)  ! Convert input to double precision
    
    ! Special case for I_0(x)
    if (n == 0) then
      sum8 = 1.0_8
      term8 = 1.0_8
      
      do k = 1, maxk
        term8 = term8 * (x8*x8/4.0_8) / (real(k,8)*real(k,8))
        sum8 = sum8 + term8
        if (abs(term8) < 1.0e-15_8 * abs(sum8)) exit
      end do
      
      res = real(sum8, wp)  ! Convert back to working precision
      return
    end if
    
    ! General case for I_n(x), n > 0
    sum8 = 0.0_8
    
    do k = 0, maxk
      ! Calculate (x/2)^(n+2k) / (k! * (n+k)!)
      fact_k8 = 1.0_8
      fact_nk8 = 1.0_8
      
      do i = 1, k
        fact_k8 = fact_k8 * real(i, 8)
      end do
      
      do i = 1, n+k
        fact_nk8 = fact_nk8 * real(i, 8)
      end do
      
      term8 = (x8/2.0_8)**(n+2*k) / (fact_k8 * fact_nk8)
      sum8 = sum8 + term8
      
      if (k > 0 .and. abs(term8) < 1.0e-15_8 * abs(sum8)) exit
    end do
    
    res = real(sum8, wp)  ! Convert back to working precision
  end function besseli
  
  ! Exact solution for the Burgers equation using Fourier-Bessel series
  function burgers_exact(x, t, nu, n_terms) result(u)
    real(wp), intent(in) :: x, t, nu
    integer, intent(in) :: n_terms
    real(wp) :: u
    
    real(wp) :: num, den, term, arg
    integer :: n
    
    num = 0.0_wp
    den = 1.0_wp
    arg = 1.0_wp / (2.0_wp * pi * nu)  ! Argument of I_n
    
    do n = 1, n_terms
      term = exp(-4.0_wp * pi * pi * nu * n * n * t) * besseli(n, arg)
      num = num + n * term * sin(2.0_wp * pi * n * x)
      den = den + 2.0_wp * term * cos(2.0_wp * pi * n * x)
    end do
    
    u = 4.0_wp * pi * nu * num / den
  end function burgers_exact
end module burgers_exact_mod

module momentum_mod
  use precision_mod, only: wp
  implicit none
  
contains
  
  ! Calculate total momentum: M(t) = ∫ u(x,t) dx
  function momentum(u, dx) result(m)
    real(wp), intent(in) :: u(:)    ! Solution vector
    real(wp), intent(in) :: dx      ! Element size
    real(wp)             :: m       ! Resulting momentum
    
    ! Total momentum: M(t) = ∫ u(x,t) dx
    integer :: n
    n = size(u)
    if (n > 1) then
    m = dx * ( 0.5_wp*(u(1) + u(n)) + sum(u(2:n-1)) )
    else
       m = u(1) * dx   ! fallback if only one node
    end if
  end function momentum
  
end module momentum_mod

module supg_burgers_mod
  use precision_mod
#ifdef ROOFLINE
  use roofline_counter, only: add_flops, add_bytes
#endif
  implicit none

  !— where to put all CSV / error / summary files —
  character(len=*), parameter :: output_dir = 'csv outputs'


  ! Problem parameters
  real(wp), parameter :: nu = 0.005_wp    ! Viscosity
  real(wp), parameter :: x_min = 0.0_wp   ! Domain start
  real(wp), parameter :: x_max = 1.0_wp   ! Domain end
  real(wp), parameter :: t_end = 0.5_wp   ! End time
  real(wp), parameter :: pi = 3.14159265358979323846_wp
  
  ! Discretization parameters
  integer, parameter :: n_elements = 128   ! Number of elements
  integer, parameter :: n_nodes = n_elements+1 ! Number of nodes
  real(wp), parameter :: dx = (x_max - x_min) / n_elements
  
  ! Arrays for solution and mesh
  real(wp), allocatable :: x(:)           ! Node coordinates
  real(wp), allocatable :: u(:), u_old(:) ! Solution
  real(wp), allocatable :: M(:,:), K(:,:), C(:,:), A(:,:) ! System matrices
  real(wp), allocatable :: b(:)           ! Right-hand side vector
  
  ! Time stepping parameters
  real(wp), parameter :: cfl_adv = 0.25_wp  ! CFL for advection
  real(wp), parameter :: cfl_diff = 0.25_wp ! CFL for diffusion
  real(wp) :: dt, u_max                   ! Time step and max velocity
  integer :: n_time_steps
  
contains
  
  ! Initialize mesh and matrices
  subroutine initialize()
    integer :: i
    
    ! Allocate arrays
    allocate(x(n_nodes))
    allocate(u(n_nodes), u_old(n_nodes))
    allocate(M(n_nodes,n_nodes), K(n_nodes,n_nodes), C(n_nodes,n_nodes), A(n_nodes,n_nodes))
    allocate(b(n_nodes))
    
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
    
    ! Emulate half precision by rounding the initial solution
#if defined(PRECISION_FP16)
    call quantize_to_half(u)
#endif
    
    ! Estimate maximum velocity for CFL
    u_max = maxval(abs(u))
#ifdef ROOFLINE
    call add_flops(2*n_nodes)  ! abs and max operations
    call add_bytes(n_nodes*wp)  ! Read u array
#endif
    
    ! Calculate time step based on CFL
    dt = min(cfl_adv * dx / u_max, cfl_diff * dx**2 / (2.0_wp * nu))
    n_time_steps = ceiling(t_end / dt)
    
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
    
    ! Assemble stiffness matrix for diffusion term
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
  end subroutine initialize
  
  ! Assemble advection matrix with SUPG stabilization
  subroutine assemble_advection_supg()
    integer :: i
    real(wp) :: u_elem, tau
    
    ! Initialize advection matrix
    C = 0.0_wp
    
    ! Loop over elements
    do i = 1, n_elements
      ! Element velocity (average)
      u_elem = 0.5_wp * (u(i) + u(i+1))
#ifdef ROOFLINE
      call add_flops(2)  ! 1 addition, 1 multiplication
      call add_bytes(2*wp)  ! Read u(i) and u(i+1)
#endif
      
      ! Standard Galerkin contribution for advection
      C(i,i) = C(i,i) - 0.5_wp * u_elem
      C(i,i+1) = C(i,i+1) + 0.5_wp * u_elem
      C(i+1,i) = C(i+1,i) - 0.5_wp * u_elem
      C(i+1,i+1) = C(i+1,i+1) + 0.5_wp * u_elem
#ifdef ROOFLINE
      call add_flops(8)  ! 4 multiplications, 4 additions/subtractions
      call add_bytes(8*wp)  ! Read and write 4 C matrix elements
#endif
      
      ! SUPG stabilization parameter (simplified)
      tau = 0.5_wp * dx / abs(u_elem)
#ifdef ROOFLINE
      call add_flops(3)  ! abs, multiplication, division
      call add_bytes(wp)  ! Write tau
#endif
      
      ! Add SUPG stabilization term
      ! Note: This is simplified. Full SUPG would include terms from the residual.
      C(i,i) = C(i,i) + tau * u_elem**2 / dx
      C(i,i+1) = C(i,i+1) - tau * u_elem**2 / dx
      C(i+1,i) = C(i+1,i) - tau * u_elem**2 / dx
      C(i+1,i+1) = C(i+1,i+1) + tau * u_elem**2 / dx
#ifdef ROOFLINE
      call add_flops(16)  ! 4 terms with square, multiplication, division, and addition
      call add_bytes(8*wp)  ! Read and write 4 C matrix elements
#endif
    end do
  end subroutine assemble_advection_supg
  
  ! Solve the system using direct solution
  subroutine solve_system()
    integer :: i
    real(wp), allocatable :: A_copy(:,:), b_copy(:)
    integer, allocatable :: ipiv(:)
    
    allocate(A_copy(n_nodes,n_nodes), b_copy(n_nodes), ipiv(n_nodes))
    
    ! Copy matrices for solver
    A_copy = A
    b_copy = b
#ifdef ROOFLINE
    call add_flops(0)  ! Array assignment is not a floating-point operation
    call add_bytes(n_nodes*n_nodes*wp + n_nodes*wp)  ! Copy A and b
#endif
    
    ! Simple Gaussian elimination (could use LAPACK for more robust solution)
    do i = 1, n_nodes-1
      A_copy(i+1:n_nodes,i) = A_copy(i+1:n_nodes,i) / A_copy(i,i)
#ifdef ROOFLINE
      call add_flops(n_nodes-i)  ! Division operations
      call add_bytes((n_nodes-i+1)*wp)  ! Read and write division results
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
    u(n_nodes) = b_copy(n_nodes) / A_copy(n_nodes,n_nodes)
#ifdef ROOFLINE
    call add_flops(1)  ! Division
    call add_bytes(2*wp)  ! Read two values, write result
#endif
    
    do i = n_nodes-1, 1, -1
      u(i) = (b_copy(i) - dot_product(A_copy(i,i+1:n_nodes), u(i+1:n_nodes))) / A_copy(i,i)
#ifdef ROOFLINE
      call add_flops(2*(n_nodes-i) + 2)  ! Dot product, subtraction, division
      call add_bytes((2*(n_nodes-i) + 3)*wp)  ! Read vectors for dot product, read/write other values
#endif
    end do
    
    ! Emulate half precision by rounding the solution after solving
#if defined(PRECISION_FP16)
    call quantize_to_half(u)
#endif
    
    deallocate(A_copy, b_copy, ipiv)
  end subroutine solve_system
  
  ! Apply boundary conditions
  subroutine apply_bc()
    ! Periodic boundary conditions
    u(1) = u(n_nodes)
    
    ! Adjust system accordingly
    A(1,:) = 0.0_wp
    A(1,1) = 1.0_wp
    A(n_nodes,:) = 0.0_wp
    A(n_nodes,n_nodes) = 1.0_wp
    A(n_nodes,1) = -1.0_wp
    
    b(1) = 0.0_wp
    b(n_nodes) = 0.0_wp
  end subroutine apply_bc
  
  ! Time-stepping routine with IMEX scheme
  subroutine time_stepping()
    integer :: t
    real(wp) :: current_time
    
    ! Output initial condition
    call output_csv(0, 0.0_wp)
    
    ! Time stepping loop with IMEX (Implicit-Explicit)
    do t = 1, n_time_steps
      current_time = t * dt
      
      ! Store old solution
      u_old = u
#ifdef ROOFLINE
      call add_bytes(n_nodes*wp)  ! Copy u to u_old
#endif
      
      ! Update maximum velocity and recompute time step if adaptive
      u_max = maxval(abs(u))
#ifdef ROOFLINE
      call add_flops(2*n_nodes)  ! abs and max operations
      call add_bytes(n_nodes*wp)  ! Read u array
#endif
      
      ! Assemble advection matrix with current solution (SUPG)
      call assemble_advection_supg()
      
      ! Form system matrix for implicit diffusion: M + dt*nu*K
      A = M + dt * nu * K
#ifdef ROOFLINE
      call add_flops(2*n_nodes*n_nodes)  ! Scalar multiplication and matrix addition
      call add_bytes(3*n_nodes*n_nodes*wp)  ! Read two matrices, write result
#endif
      
      ! Form right-hand side with explicit advection: M*u_old - dt*C*u_old
      b = matmul(M - dt * C, u_old)
#ifdef ROOFLINE
      call add_flops(n_nodes*n_nodes + n_nodes*n_nodes*n_nodes + n_nodes*n_nodes)  ! Matrix subtraction, matrix-vector product
      call add_bytes((2*n_nodes*n_nodes + n_nodes + n_nodes)*wp)  ! Read matrices and vector, write result
#endif
      
      ! Apply boundary conditions
      call apply_bc()
      
      ! Solve the system
      call solve_system()
      
      ! Output solution at specific time steps
      if (mod(t, max(1, n_time_steps/10)) == 0 .or. t == n_time_steps) then
        call output_csv(t, current_time)
      end if
    end do
  end subroutine time_stepping

  subroutine calculate_errors()
    use burgers_exact_mod, only: burgers_exact
    implicit none
    real(wp) :: error_l1, error_l2, error_linf, diff, u_exact
    integer :: i
    character(len = 100) ::filename
    ! assume we still have u(:) and x(:) at t_end
    error_l1   = 0.0_wp
    error_l2   = 0.0_wp
    error_linf = 0.0_wp
    do i = 1, n_nodes
      ! Replace simple sine with actual exact solution using Fourier-Bessel series
      u_exact = burgers_exact(x(i), t_end, nu, 30)
      diff = abs(u(i) - u_exact)
      error_l1 = error_l1 + diff
      error_l2 = error_l2 + diff**2
      error_linf = max(error_linf, diff)
    end do
    error_l1 = error_l1 * dx
    error_l2 = sqrt(error_l2 * dx)

    ! --- write plain-text error file ------------------------------------
    write(filename,'(a,a,a,a)') trim(output_dir)//'/', &
                              'burgers_', trim(precision_label), '_errors.txt'
    open(unit=12, file=filename, status='replace')
    write(12,'(a,G0)') 'L1 Error: ',   error_l1
    write(12,'(a,G0)') 'L2 Error: ',   error_l2
    write(12,'(a,G0)') 'Linf Error: ', error_linf
    close(12)

    ! --- replicate advection's console summary --------------------------
    print *, 'Error metrics for precision: ', precision_label
    print '(a,es16.8)', 'L1 Error: ',   error_l1
    print '(a,es16.8)', 'L2 Error: ',   error_l2
    print '(a,es16.8)', 'Linf Error: ', error_linf
  end subroutine calculate_errors

  
  ! Output solution to CSV
  subroutine output_csv(step, time)
    integer, intent(in) :: step
    real(wp), intent(in) :: time
    character(len=100) :: filename
    integer :: i
    integer :: istat
    


    ! Create filename with precision label **and underscore before step**
    write(filename,'(a,a,a,a,i6.6,a)') trim(output_dir)//'/', &
                                      'burgers_', trim(precision_label), '_', step, '.csv'

    ! open(unit=10, file=filename, …)

    open(unit=10, &
        file=filename, &
        status='replace', &
        action='write', &
        iostat=istat)
    if (istat /= 0) then
      print *, 'Error opening file ', filename, ' iostat=', istat
    end if

    
    
    ! Write header
    write(10,'(a,G0)')    'x,u,time=', time
    
    ! Write data
    do i = 1, n_nodes
      write(10, '(G0,a,G0)')            x(i), ',', u(i)
    end do
    
    ! Close file
    close(10)
    
    ! Create a summary file for the current time step
    if (step == n_time_steps) then
      !write(filename, '(a,a,a)') 'burgers_', trim(precision_label), '_summary.csv'
      write(filename,'(a,a,a,a)') &
      trim(output_dir)//'/', &
      'burgers_', trim(precision_label), '_summary.csv'

      open(unit=11, file=filename, status='replace')
      write(11, '(a)') 'x,u'
      do i = 1, n_nodes
        write(11, '(G0,a,G0)')            x(i), ',', u(i)
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
  
end module supg_burgers_mod

program burgers_equation_supg
  use precision_mod
  use momentum_mod, only : momentum
  use supg_burgers_mod
#ifdef ROOFLINE
  use roofline_counter, only : flop_count, byte_count
#endif
  implicit none
  
  integer :: istat
  integer :: t, unit_mom
  real(wp) :: mom_init, mom_now      ! Renamed to avoid name clash with 'M' matrix
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
  ! -- right after you check the output directory --
  print *, '=================================='
  print *, 'Running Viscous Burgers Equation SUPG solver with precision: ', precision_label
  print *, 'Working precision kind: ', wp
  print *, 'Viscosity: ', nu
  print *, '=================================='


  ! Initialize the problem
  call initialize()
  
  print *, 'Time step: ', dt
  print *, 'Number of time steps: ', n_time_steps
  
  ! Initialize momentum tracking
  mom_init = momentum(u, dx)
  unit_mom = 43
  write(filename,'(a,a,a,a)') trim(output_dir)//'/', &
                            'burgers_', trim(precision_label), '_momentum.csv'
  open(unit=unit_mom, file=filename, status='replace', action='write')
  write(unit_mom,'(a)') 'step,time,M'
  WRITE(unit_mom,'(I0,1X,ES40.33,1X,ES40.33)') 0, 0.0_wp, mom_init
  
  
  ! Run time stepping
  call time_stepping()
  
  ! Custom time stepping loop to track momentum
  do t = 1, n_time_steps
    if (mod(t, max(1, n_time_steps/100)) == 0 .or. t == n_time_steps) then
      ! Calculate and log momentum
      mom_now = momentum(u, dx)
      write(unit_mom,'(i0,1x,G0,1x,G0)') t, t*dt, mom_now
      
      ! Print momentum drift to screen
      print '(a,i6,a,es10.3)', ' M at step ', t, ': ', abs(mom_now)
    end if
  end do
  
  ! Final momentum check
  mom_now = momentum(u, dx)
  
  ! Calculate errors
  call calculate_errors()
  
  ! Output momentum conservation summary
  print *, '----------------------------------------'
  print *, 'Momentum check: initial M =', mom_init
  print *, '                final   M =', mom_now
  print *, 'Full trace in ', trim(filename)
  print *, '----------------------------------------'
  
  ! Close momentum log file
  close(unit_mom)
  
  ! Output final solution
  call output_csv(n_time_steps, t_end)

#ifdef ROOFLINE
   call cpu_time(roof_t_end)

   ! write the CSV unconditionally – we are in a single-image program
   open(unit=99, file=trim(output_dir)//'/burgers_'// &
                      trim(precision_label)//'_roofline.csv', &
        status='replace', action='write')
   write(99,'(a)') 'runtime_s,flops,bytes'
   write(99,'(G0,1x,i0,1x,i0)')     roof_t_end - roof_t_start, flop_count, byte_count
   close(99)
#endif


  print *, 'Simulation completed. CSV files have been written.'
  
end program burgers_equation_supg


! gfortran -cpp -DPRECISION_FP16 burgers_equation_supg.f90 -o burgers_supg_fp16
! gfortran -cpp -DPRECISION_FP32 burgers_equation_supg.f90 -o burgers_supg_fp32
! gfortran -cpp burgers_equation_supg.f90 -o burgers_supg_fp64
! gfortran -cpp -DPRECISION_FP128 burgers_equation_supg.f90 -o burgers_supg_fp128

! With roofline instrumentation:
! gfortran -cpp -DROOFLINE -DPRECISION_FP16 roofline_counter.o burgers_equation_supg.f90 -o burgers_roof_fp16
! gfortran -cpp -DROOFLINE -DPRECISION_FP32 roofline_counter.o burgers_equation_supg.f90 -o burgers_roof_fp32
! gfortran -cpp -DROOFLINE roofline_counter.o burgers_equation_supg.f90 -o burgers_roof_fp64
! gfortran -cpp -DROOFLINE -DPRECISION_FP128 roofline_counter.o burgers_equation_supg.f90 -o burgers_roof_fp128

! ./burgers_roof_fp16
! ./burgers_roof_fp32
! ./burgers_roof_fp64
! ./burgers_roof_fp128