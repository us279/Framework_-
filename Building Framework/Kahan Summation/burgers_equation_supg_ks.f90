! 1D Viscous Burgers Equation Solver using SUPG
! Equation: du/dt + u * du/dx = nu * d^2u/dx^2
! Solves with SUPG stabilization and IMEX time stepping

! Compile with:
! gfortran -cpp -DPRECISION_FP16 burgers_equation_supg_ks.f90 -o burgers_supg_ks_fp16
! gfortran -cpp -DPRECISION_FP32 burgers_equation_supg_ks.f90 -o burgers_supg_ks_fp32
! gfortran -cpp burgers_equation_supg_ks.f90 -o burgers_supg_ks_fp64
! gfortran -cpp -DPRECISION_FP128 burgers_equation_supg_ks.f90 -o burgers_supg_ks_fp128

! Compile with ROOFLINE support:
! gfortran -cpp -DROOFLINE -DPRECISION_FP16 roofline_counter.o burgers_equation_supg_ks.f90 -o burgers_roof_ks_fp16
! gfortran -cpp -DROOFLINE -DPRECISION_FP32 roofline_counter.o burgers_equation_supg_ks.f90 -o burgers_roof_ks_fp32
! gfortran -cpp -DROOFLINE roofline_counter.o burgers_equation_supg_ks.f90 -o burgers_roof_ks_fp64
! gfortran -cpp -DROOFLINE -DPRECISION_FP128 roofline_counter.o burgers_equation_supg_ks.f90 -o burgers_roof_ks_fp128

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

! Module for exact solution of Burgers equation
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

! Kahan summation module for compensated arithmetic
module kahan_mod
  use precision_mod, only: wp
  implicit none
contains
!---------------------------------------------------------------
function kahan_sum(a) result(s)
  real(wp), intent(in) :: a(:)
  real(wp)             :: s
  real(wp)             :: c, y, t
  integer              :: i
  s = 0.0_wp
  c = 0.0_wp
  do i = 1, size(a)
     y = a(i) - c
     t = s + y
     c = (t - s) - y
     s = t
  end do
end function kahan_sum
!---------------------------------------------------------------
pure subroutine kahan_axpy(x, a, y, c)   ! x := x + a*y  (in-place)
  real(wp), intent(inout) :: x(:)        ! target vector
  real(wp), intent(in)    :: a           ! scalar
  real(wp), intent(in)    :: y(:)        ! increment
  real(wp), intent(inout) :: c(:)        ! same size as x – running compensation
  real(wp) :: yi, tmp
  integer  :: i
  do i = 1, size(x)
     yi  = a*y(i) - c(i)
     tmp = x(i) + yi
     c(i) = (tmp - x(i)) - yi
     x(i) = tmp
  end do
end subroutine kahan_axpy
!---------------------------------------------------------------
pure subroutine kahan_add_scalar(x, dx, c) ! x := x + dx  (scalar version)
  real(wp), intent(inout) :: x
  real(wp), intent(in)    :: dx
  real(wp), intent(inout) :: c
  real(wp) :: y, t
  y = dx - c
  t = x + y
  c = (t - x) - y
  x = t
end subroutine kahan_add_scalar
!---------------------------------------------------------------
end module kahan_mod

module momentum_mod
  use precision_mod, only: wp
  use kahan_mod, only: kahan_sum  ! Add kahan module
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
      m = dx * ( 0.5_wp*(u(1) + u(n)) + kahan_sum(u(2:n-1)) )  ! Use kahan_sum
    else
      m = u(1) * dx   ! fallback if only one node
    end if
  end function momentum
  
end module momentum_mod

module supg_burgers_mod
  use precision_mod
  use kahan_mod, only: kahan_add_scalar, kahan_axpy, kahan_sum  ! Add Kahan module
  use momentum_mod, only: momentum
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
  
  ! Kahan compensation variables
  real(wp) :: time_c = 0.0_wp            ! Scalar time compensation
  real(wp), allocatable :: kahan_u(:)     ! Per-node compensation for u-updates
  
  ! New: Persistent compensation arrays for matrix operations
  real(wp), allocatable :: compA(:,:)     ! Compensation for A matrix 
  real(wp), allocatable :: compB(:)       ! Compensation for b vector
  real(wp), allocatable :: compTemp(:,:)  ! Compensation for temporary matrix
  
  ! Time stepping parameters
  real(wp), parameter :: cfl_adv = 0.25_wp  ! CFL for advection
  real(wp), parameter :: cfl_diff = 0.25_wp ! CFL for diffusion
  real(wp) :: dt, u_max                   ! Time step and max velocity
  integer :: n_time_steps
  real(wp) :: current_time = 0.0_wp       ! Current simulation time
  
contains
  
  ! Initialize mesh and matrices
  subroutine initialize()
    integer :: i
    
    ! Allocate arrays
    allocate(x(n_nodes))
    allocate(u(n_nodes), u_old(n_nodes))
    allocate(M(n_nodes,n_nodes), K(n_nodes,n_nodes), C(n_nodes,n_nodes), A(n_nodes,n_nodes))
    allocate(b(n_nodes))
    
    ! Allocate and initialize Kahan compensation arrays
    allocate(kahan_u(n_nodes))
    kahan_u = 0.0_wp
    time_c = 0.0_wp
    
    ! Allocate persistent compensation arrays for matrix operations
    allocate(compA(n_nodes,n_nodes), compB(n_nodes), compTemp(n_nodes,n_nodes))
    compA = 0.0_wp
    compB = 0.0_wp
    compTemp = 0.0_wp
    
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
  
  ! Solve the system using direct solution with Kahan summation
  subroutine solve_system()
    integer :: i, j, k
    real(wp), allocatable :: A_copy(:,:), b_copy(:)
    integer, allocatable :: ipiv(:)
    real(wp) :: dot_result
    real(wp), allocatable :: dot_comp_array(:)
    
    allocate(A_copy(n_nodes,n_nodes), b_copy(n_nodes), ipiv(n_nodes))
    allocate(dot_comp_array(n_nodes))
    
    ! Copy matrices for solver
    A_copy = A
    b_copy = b
#ifdef ROOFLINE
    call add_flops(0)  ! Array assignment is not a floating-point operation
    call add_bytes(n_nodes*n_nodes*wp + n_nodes*wp)  ! Copy A and b
#endif
    
    ! Reset compensation arrays for this solve
    compA = 0.0_wp
    compB = 0.0_wp
    
    ! Simple Gaussian elimination (could use LAPACK for more robust solution)
    do i = 1, n_nodes-1
      A_copy(i+1:n_nodes,i) = A_copy(i+1:n_nodes,i) / A_copy(i,i)
#ifdef ROOFLINE
      call add_flops(n_nodes-i)  ! Division operations
      call add_bytes((n_nodes-i+1)*wp)  ! Read and write division results
#endif
      
      ! Matrix update with Kahan summation for each element
      do j = i+1, n_nodes
        do k = i+1, n_nodes
          ! Use persistent compensation array instead of temporary
          call kahan_add_scalar(A_copy(j,k), -A_copy(j,i) * A_copy(i,k), compA(j,k))
        end do
      end do
      
#ifdef ROOFLINE
      call add_flops(3*(n_nodes-i)*(n_nodes-i))  ! More operations due to Kahan
      call add_bytes(4*(n_nodes-i)*(n_nodes-i)*wp)  ! More memory accesses due to Kahan
#endif
      
      ! Update RHS with Kahan summation
      do j = i+1, n_nodes
        ! Use persistent compensation array
        call kahan_add_scalar(b_copy(j), -A_copy(j,i) * b_copy(i), compB(j))
      end do
      
#ifdef ROOFLINE
      call add_flops(3*(n_nodes-i))  ! More operations for Kahan
      call add_bytes(4*(n_nodes-i)*wp)  ! More memory accesses for Kahan
#endif
    end do
    
    ! Back substitution with Kahan summation
    u(n_nodes) = b_copy(n_nodes) / A_copy(n_nodes,n_nodes)
#ifdef ROOFLINE
    call add_flops(1)  ! Division
    call add_bytes(3*wp)  ! Read two values, write result
#endif
    
    ! Initialize dot product compensation array
    dot_comp_array = 0.0_wp
    
    do i = n_nodes-1, 1, -1
      dot_result = 0.0_wp
      
      ! Compute dot product with Kahan summation using persistent compensation
      do j = i+1, n_nodes
        call kahan_add_scalar(dot_result, A_copy(i,j) * u(j), dot_comp_array(i))
      end do
      
      u(i) = (b_copy(i) - dot_result) / A_copy(i,i)
#ifdef ROOFLINE
      call add_flops(3*(n_nodes-i) + 2)  ! More ops for Kahan plus div/sub
      call add_bytes((4*(n_nodes-i) + 3)*wp)  ! More memory access for Kahan
#endif
    end do
    
    ! Emulate half precision by rounding the solution after solving
#if defined(PRECISION_FP16)
    call quantize_to_half(u)
#endif
    
    deallocate(A_copy, b_copy, ipiv, dot_comp_array)
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
  
  ! Time-stepping routine with IMEX scheme and Kahan summation
  subroutine time_stepping(unit_mom)
    integer, intent(in) :: unit_mom  ! File unit for momentum output
    integer :: t, i, j
    real(wp) :: mom_now  ! Current momentum
    real(wp), allocatable :: temp_matrix(:,:)
    
    ! Output initial condition
    call output_csv(0, 0.0_wp)
    
    ! Initialize time
    current_time = 0.0_wp
    
    ! Allocate temporary storage for matrix operations
    allocate(temp_matrix(n_nodes, n_nodes))
    
    ! Initialize compensation arrays for matrix operations in time loop
    compTemp = 0.0_wp
    
    ! Time stepping loop with IMEX (Implicit-Explicit)
    do t = 1, n_time_steps
      ! Update time with Kahan compensation
      call kahan_add_scalar(current_time, dt, time_c)
      
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
      ! Calculate with Kahan summation
      do i = 1, n_nodes
        do j = 1, n_nodes
          ! Use persistent compensation array
          A(i,j) = M(i,j)
          call kahan_add_scalar(A(i,j), dt * nu * K(i,j), compA(i,j))
        end do
      end do
      
#ifdef ROOFLINE
      call add_flops(3*n_nodes*n_nodes)  ! More operations for Kahan
      call add_bytes(4*n_nodes*n_nodes*wp)  ! More memory accesses for Kahan
#endif
      
      ! Form right-hand side with explicit advection: M*u_old - dt*C*u_old
      ! First calculate M - dt * C with Kahan
      do i = 1, n_nodes
        do j = 1, n_nodes
          temp_matrix(i,j) = M(i,j)
          call kahan_add_scalar(temp_matrix(i,j), -dt * C(i,j), compTemp(i,j))
        end do
      end do
      
      ! Now calculate matrix-vector product with Kahan
      b = 0.0_wp
      
      do i = 1, n_nodes
        do j = 1, n_nodes
          call kahan_add_scalar(b(i), temp_matrix(i,j) * u_old(j), compB(i))
        end do
      end do
      
#ifdef ROOFLINE
      call add_flops(3*n_nodes*n_nodes + 3*n_nodes*n_nodes)  ! Matrix-matrix and matrix-vector with Kahan
      call add_bytes(4*n_nodes*n_nodes*wp + 4*n_nodes*n_nodes*wp)  ! Memory access for both operations
#endif
      
      ! Apply boundary conditions
      call apply_bc()
      
      ! Solve the system (solution stored directly in u)
      call solve_system()
      
      ! Output solution and momentum at specific time steps
      if (mod(t, max(1, n_time_steps/100)) == 0 .or. t == n_time_steps) then
        call output_csv(t, current_time)
        
        ! Calculate and log momentum
        mom_now = momentum(u, dx)
        write(unit_mom,'(I0,1X,G0,1X,G0)') t, current_time, mom_now
        
        ! Print momentum
        print '(a,i6,a,es10.3)', ' M at step ', t, ': ', mom_now
      end if
    end do
    
    ! Clean up
    deallocate(temp_matrix)
    
    ! Calculate and output errors
    call calculate_errors()
  end subroutine time_stepping

  subroutine calculate_errors()
    use burgers_exact_mod, only: burgers_exact
    implicit none
    real(wp) :: error_l1, error_l2, error_linf, diff, u_exact
    integer :: i
    character(len = 100) ::filename
    real(wp), allocatable :: diff_array(:), diff_squared(:)
    
    ! Allocate arrays for Kahan summation
    allocate(diff_array(n_nodes), diff_squared(n_nodes))
    
    ! Initialize errors
    error_l1   = 0.0_wp
    error_l2   = 0.0_wp
    error_linf = 0.0_wp
    
    ! Calculate differences and prepare arrays for Kahan summation
    do i = 1, n_nodes
      ! Use the exact solution from the Fourier-Bessel series
      u_exact = burgers_exact(x(i), t_end, nu, 30)
      diff = abs(u(i) - u_exact)
      diff_array(i) = diff
      diff_squared(i) = diff**2
      error_linf = max(error_linf, diff)
    end do
    
    ! Use Kahan summation for L1 and L2 errors
    error_l1 = kahan_sum(diff_array) * dx
    error_l2 = sqrt(kahan_sum(diff_squared) * dx)
    
    ! Clean up
    deallocate(diff_array, diff_squared)

    ! --- write plain-text error file ------------------------------------
    write(filename,'(a,a,a,a)') trim(output_dir)//'/', &
                              'burgers_', trim(precision_label), '_errors.txt'
    open(unit=12, file=filename, status='replace')
    write(12,'(A,G0)') 'L1 Error: ',   error_l1
    write(12,'(A,G0)') 'L2 Error: ',   error_l2
    write(12,'(A,G0)') 'Linf Error: ', error_linf
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

    open(unit=10, &
        file=filename, &
        status='replace', &
        action='write', &
        iostat=istat)
    if (istat /= 0) then
      print *, 'Error opening file ', filename, ' iostat=', istat
    end if
    
    ! Write header
    write(10, '(A,G0)') 'x,u,time=', time
    
    ! Write data
    do i = 1, n_nodes
      write(10, '(G0,A,G0)') x(i), ',', u(i)
    end do
    
    ! Close file
    close(10)
    
    ! Create a summary file for the current time step
    if (step == n_time_steps) then
      write(filename,'(a,a,a,a)') &
      trim(output_dir)//'/', &
      'burgers_', trim(precision_label), '_summary.csv'

      open(unit=11, file=filename, status='replace')
      write(11, '(a)') 'x,u'
      do i = 1, n_nodes
        write(11, '(f12.6,a,es20.10)') x(i), ',', u(i)
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
  integer :: unit_mom
  real(wp) :: mom_init, mom_final      ! Initial and final momentum
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
  write(unit_mom,'(I0,1X,G0,1X,G0)') 0, 0.0_wp, mom_init
  
  ! Run time stepping with Kahan summation - pass the momentum file unit
  call time_stepping(unit_mom)
  
  ! Final momentum calculation
  mom_final = momentum(u, dx)
  
  ! Output momentum conservation summary
  print *, '----------------------------------------'
  print *, 'Momentum check: initial M =', mom_init
  print *, '                final   M =', mom_final
  print *, 'Full trace in ', trim(filename)
  print *, '----------------------------------------'
  
  ! Close momentum log file
  close(unit_mom)

#ifdef ROOFLINE
   call cpu_time(roof_t_end)

   ! write the CSV unconditionally – we are in a single-image program
   open(unit=99, file=trim(output_dir)//'/burgers_'// &
                      trim(precision_label)//'_roofline.csv', &
        status='replace', action='write')
   write(99,'(a)') 'runtime_s,flops,bytes'
   write(99,'(G0,1X,I0,1X,I0)')  &
        roof_t_end - roof_t_start, flop_count, byte_count
   close(99)
#endif

  print *, 'Simulation completed. CSV files have been written.'
  
end program burgers_equation_supg


! gfortran -cpp -DPRECISION_FP16 burgers_equation_supg_ks.f90 -o burgers_supg_ks_fp16
! gfortran -cpp -DPRECISION_FP32 burgers_equation_supg_ks.f90 -o burgers_supg_ks_fp32
! gfortran -cpp burgers_equation_supg_ks.f90 -o burgers_supg_ks_fp64
! gfortran -cpp -DPRECISION_FP128 burgers_equation_supg_ks.f90 -o burgers_supg_ks_fp128

! With roofline instrumentation:
! gfortran -cpp -DROOFLINE -DPRECISION_FP16 roofline_counter.o burgers_equation_supg_ks.f90 -o burgers_roof_ks_fp16
! gfortran -cpp -DROOFLINE -DPRECISION_FP32 roofline_counter.o burgers_equation_supg_ks.f90 -o burgers_roof_ks_fp32
! gfortran -cpp -DROOFLINE roofline_counter.o burgers_equation_supg_ks.f90 -o burgers_roof_ks_fp64
! gfortran -cpp -DROOFLINE -DPRECISION_FP128 roofline_counter.o burgers_equation_supg_ks.f90 -o burgers_roof_ks_fp128

! ./burgers_roof_ks_fp16
! ./burgers_roof_ks_fp32
! ./burgers_roof_ks_fp64
! ./burgers_roof_ks_fp128