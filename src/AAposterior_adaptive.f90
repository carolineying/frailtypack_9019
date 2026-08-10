module posterior_adaptive_settings
    implicit none

    integer,parameter :: adaptive_nidiag = 20
    integer,parameter :: adaptive_nddiag = 10
    integer,save :: adaptive_method = 0
    integer,save :: adaptive_mode_maxit = 80
    integer,save :: adaptive_fallback = 1
    integer,save :: adaptive_diagnostics_enabled = 0
    integer,save :: adaptive_warm_start = 1
    integer,save :: adaptive_nonpositive_action = 0
    double precision,save :: adaptive_mode_tol = 1.d-8
    double precision,save :: adaptive_hessian_eps = 1.d-10

    integer,save :: posterior_adaptive_active = 0
    integer,save :: posterior_adaptive_fallback_active = 0
    integer,save :: posterior_callback_returns_log = 0
    integer,save :: posterior_mode_status = 0
    integer,save :: posterior_mode_iterations = 0
    integer,save :: posterior_hessian_regularized = 0
    double precision,save :: posterior_log_mode = 0.d0
    double precision,save :: posterior_log_jacobian = 0.d0
    double precision,save :: posterior_proposal_scale = 1.d0
    double precision,save :: posterior_gradient_norm = 0.d0
    double precision,save :: posterior_min_eigenvalue = 0.d0
    double precision,save :: posterior_min_stabilized = 0.d0
    double precision,save :: posterior_condition_number = 0.d0

    double precision,dimension(:),allocatable,save :: posterior_mode_z
    double precision,dimension(:,:),allocatable,save :: posterior_chol_h
    double precision,dimension(:),allocatable,save :: posterior_mode_cache_b
    integer,dimension(:),allocatable,save :: posterior_mode_cache_valid
    integer,dimension(:),allocatable,save :: adaptive_floor_count
    integer,dimension(:,:),allocatable,save :: adaptive_floor_reason_count
    integer,dimension(:,:),allocatable,save :: adaptive_diag_i
    double precision,dimension(:,:),allocatable,save :: adaptive_diag_d

    ! The Genz--Keister rules above dimension four have signed rule weights.
    ! In adaptive mode HRMTRL accumulates positive and negative contributions
    ! separately on the log scale.  These variables persist across successive
    ! embedded orders within one HRMSYM call.
    double precision,dimension(:),allocatable,save :: gk_log_positive
    double precision,dimension(:),allocatable,save :: gk_log_negative

contains

subroutine initialize_posterior_adaptive(nsujet,ng,adaptive_i,adaptive_d)
    implicit none
    integer,intent(in)::nsujet,ng
    integer,dimension(6),intent(in)::adaptive_i
    double precision,dimension(2),intent(in)::adaptive_d

    adaptive_method = adaptive_i(1)
    adaptive_mode_maxit = max(1,adaptive_i(2))
    adaptive_fallback = adaptive_i(3)
    adaptive_diagnostics_enabled = adaptive_i(4)
    adaptive_warm_start = adaptive_i(5)
    adaptive_nonpositive_action = adaptive_i(6)
    adaptive_mode_tol = adaptive_d(1)
    adaptive_hessian_eps = adaptive_d(2)

    posterior_adaptive_active = 0
    posterior_adaptive_fallback_active = 0
    posterior_callback_returns_log = 0
    if(allocated(posterior_mode_z))deallocate(posterior_mode_z)
    if(allocated(posterior_chol_h))deallocate(posterior_chol_h)
    if(allocated(posterior_mode_cache_b))deallocate(posterior_mode_cache_b)
    if(allocated(posterior_mode_cache_valid))deallocate(posterior_mode_cache_valid)
    if(allocated(adaptive_floor_count))deallocate(adaptive_floor_count)
    if(allocated(adaptive_floor_reason_count))deallocate(adaptive_floor_reason_count)
    if(allocated(adaptive_diag_i))deallocate(adaptive_diag_i)
    if(allocated(adaptive_diag_d))deallocate(adaptive_diag_d)
    if(allocated(gk_log_positive))deallocate(gk_log_positive)
    if(allocated(gk_log_negative))deallocate(gk_log_negative)

    allocate(posterior_mode_cache_b(nsujet),posterior_mode_cache_valid(ng))
    allocate(adaptive_floor_count(ng))
    allocate(adaptive_floor_reason_count(ng,4))
    allocate(adaptive_diag_i(ng,adaptive_nidiag))
    allocate(adaptive_diag_d(ng,adaptive_nddiag))
    posterior_mode_cache_b = 0.d0
    posterior_mode_cache_valid = 0
    adaptive_floor_count = 0
    adaptive_floor_reason_count = 0
    adaptive_diag_i = 0
    adaptive_diag_d = 0.d0
end subroutine initialize_posterior_adaptive


subroutine finalize_posterior_adaptive()
    implicit none
    posterior_adaptive_active = 0
    posterior_adaptive_fallback_active = 0
    posterior_callback_returns_log = 0
    if(allocated(posterior_mode_z))deallocate(posterior_mode_z)
    if(allocated(posterior_chol_h))deallocate(posterior_chol_h)
    if(allocated(posterior_mode_cache_b))deallocate(posterior_mode_cache_b)
    if(allocated(posterior_mode_cache_valid))deallocate(posterior_mode_cache_valid)
    if(allocated(adaptive_floor_count))deallocate(adaptive_floor_count)
    if(allocated(adaptive_floor_reason_count))deallocate(adaptive_floor_reason_count)
    if(allocated(adaptive_diag_i))deallocate(adaptive_diag_i)
    if(allocated(adaptive_diag_d))deallocate(adaptive_diag_d)
    if(allocated(gk_log_positive))deallocate(gk_log_positive)
    if(allocated(gk_log_negative))deallocate(gk_log_negative)
end subroutine finalize_posterior_adaptive


subroutine reset_posterior_family()
    implicit none
    posterior_adaptive_active = 0
    posterior_adaptive_fallback_active = 0
    posterior_callback_returns_log = 0
    posterior_mode_status = 0
    posterior_mode_iterations = 0
    posterior_hessian_regularized = 0
    posterior_log_mode = 0.d0
    posterior_log_jacobian = 0.d0
    posterior_proposal_scale = 1.d0
    posterior_gradient_norm = 0.d0
    posterior_min_eigenvalue = 0.d0
    posterior_min_stabilized = 0.d0
    posterior_condition_number = 0.d0
    if(allocated(posterior_mode_z))deallocate(posterior_mode_z)
    if(allocated(posterior_chol_h))deallocate(posterior_chol_h)
end subroutine reset_posterior_family


subroutine setup_posterior_adaptive(family,offset,n,prior_chol,delta,hazard,status)
    implicit none
    integer,intent(in)::family,offset,n
    double precision,dimension(n,n),intent(in)::prior_chol
    double precision,dimension(n),intent(in)::delta,hazard
    integer,intent(out)::status

    integer::attempt,local_status,iterations
    double precision::grad_norm,log_value,min_eig,min_stab,condition
    integer::regularized
    double precision,dimension(:),allocatable::start,z,gradient,eigenvalues
    double precision,dimension(:,:),allocatable::hessian,chol,eigenvectors

    call reset_posterior_family()
    status = 0
    if(adaptive_method.ne.1.or.n.le.0)return

    allocate(start(n),z(n),gradient(n),eigenvalues(n))
    allocate(hessian(n,n),chol(n,n),eigenvectors(n,n))

    start = 0.d0
    if(adaptive_warm_start.eq.1.and.allocated(posterior_mode_cache_valid))then
        if(posterior_mode_cache_valid(family).eq.1)then
            start = posterior_mode_cache_b(offset:offset+n-1)
            call solve_lower(prior_chol,start,local_status)
            if(local_status.ne.0)start = 0.d0
        end if
    end if

    local_status = -1
    do attempt=1,3
        if(attempt.eq.2)start = 0.d0
        if(attempt.eq.3)then
            ! Explicit one-step Newton fallback from the prior mean.
            start = 0.d0
            call posterior_state(n,prior_chol,delta,hazard,start,log_value,&
                gradient,hessian,local_status)
            if(local_status.ne.0)cycle
            call cholesky_lower(n,hessian,chol,local_status)
            if(local_status.ne.0)cycle
            call solve_cholesky(n,chol,gradient,start)
            start = 0.5d0*start
        end if
        call safeguarded_newton(n,prior_chol,delta,hazard,start,z,iterations,&
            grad_norm,local_status)
        posterior_mode_iterations = iterations
        posterior_gradient_norm = grad_norm
        if(local_status.eq.1)exit
    end do

    if(local_status.ne.1)then
        ! A proposal need not be centred at the exact conditional mode for
        ! adaptive quadrature to remain exact.  Preserve the last finite
        ! safeguarded-Newton iterate and its curvature as an approximate
        ! proposal before resorting to the much poorer prior-centred rule.
        status = local_status
        call posterior_state(n,prior_chol,delta,hazard,z,log_value,gradient,&
            hessian,attempt)
        if(attempt.eq.0)then
            grad_norm = maxval(abs(gradient))
            call stabilize_hessian(n,hessian,adaptive_hessian_eps,chol,&
                eigenvectors,eigenvalues,min_eig,min_stab,condition,&
                regularized,attempt)
            if(attempt.eq.0)then
                allocate(posterior_mode_z(n),posterior_chol_h(n,n))
                posterior_mode_z = z
                posterior_chol_h = chol
                posterior_log_mode = log_value
                posterior_gradient_norm = grad_norm
                posterior_min_eigenvalue = min_eig
                posterior_min_stabilized = min_stab
                posterior_condition_number = condition
                posterior_hessian_regularized = regularized
                posterior_mode_status = status
                posterior_adaptive_active = 1
                call set_posterior_proposal_scale(n,1.d0)
                call store_mode_diagnostics(family,n)
                deallocate(start,z,gradient,eigenvalues,hessian,chol,eigenvectors)
                return
            end if
        end if
        posterior_mode_status = status
        call store_mode_diagnostics(family,n)
        deallocate(start,z,gradient,eigenvalues,hessian,chol,eigenvectors)
        return
    end if

    call posterior_state(n,prior_chol,delta,hazard,z,log_value,gradient,&
        hessian,local_status)
    if(local_status.ne.0)then
        posterior_mode_status = -4
        status = -4
        call store_mode_diagnostics(family,n)
        deallocate(start,z,gradient,eigenvalues,hessian,chol,eigenvectors)
        return
    end if

    grad_norm = maxval(abs(gradient))
    call stabilize_hessian(n,hessian,adaptive_hessian_eps,chol,eigenvectors,&
        eigenvalues,min_eig,min_stab,condition,regularized,local_status)
    if(local_status.ne.0)then
        posterior_mode_status = -5
        status = -5
        call store_mode_diagnostics(family,n)
        deallocate(start,z,gradient,eigenvalues,hessian,chol,eigenvectors)
        return
    end if

    allocate(posterior_mode_z(n),posterior_chol_h(n,n))
    posterior_mode_z = z
    posterior_chol_h = chol
    posterior_log_mode = log_value
    posterior_gradient_norm = grad_norm
    posterior_min_eigenvalue = min_eig
    posterior_min_stabilized = min_stab
    posterior_condition_number = condition
    posterior_hessian_regularized = regularized
    posterior_mode_iterations = iterations
    posterior_mode_status = 1
    posterior_adaptive_active = 1
    call set_posterior_proposal_scale(n,1.d0)
    status = 1

    if(adaptive_warm_start.eq.1)then
        posterior_mode_cache_b(offset:offset+n-1) = matmul(prior_chol,z)
        posterior_mode_cache_valid(family) = 1
    end if
    call store_mode_diagnostics(family,n)

    deallocate(start,z,gradient,eigenvalues,hessian,chol,eigenvectors)
end subroutine setup_posterior_adaptive


subroutine safeguarded_newton(n,prior_chol,delta,hazard,start,z,iterations,&
    grad_norm,status)
    implicit none
    integer,intent(in)::n
    double precision,dimension(n,n),intent(in)::prior_chol
    double precision,dimension(n),intent(in)::delta,hazard,start
    double precision,dimension(n),intent(out)::z
    integer,intent(out)::iterations,status
    double precision,intent(out)::grad_norm

    integer::iter,line_iter,local_status
    double precision::old_log,new_log,alpha,directional,relative_step,max_du
    double precision,dimension(n)::gradient,step,z_trial,trial_gradient,du
    double precision,dimension(n,n)::hessian,chol,trial_hessian

    z = start
    iterations = 0
    grad_norm = huge(1.d0)
    status = -1
    relative_step = huge(1.d0)

    do iter=1,adaptive_mode_maxit
        call posterior_state(n,prior_chol,delta,hazard,z,old_log,gradient,&
            hessian,local_status)
        if(local_status.ne.0)then
            status = -2
            return
        end if
        grad_norm = maxval(abs(gradient))
        if(grad_norm.le.adaptive_mode_tol*(1.d0+maxval(abs(z))).and.&
            (iter.eq.1.or.relative_step.le.dsqrt(adaptive_mode_tol)))then
            status = 1
            iterations = iter-1
            return
        end if

        call cholesky_lower(n,hessian,chol,local_status)
        if(local_status.ne.0)then
            status = -3
            return
        end if
        call solve_cholesky(n,chol,gradient,step)
        directional = dot_product(gradient,step)
        if(directional.le.0.d0)then
            status = -3
            return
        end if

        ! Limit the change in the frailty linear predictor before Armijo
        ! backtracking.  This prevents otherwise valid Newton directions from
        ! proposing exp(Uz) overflow and repeatedly collapsing the step.
        du = matmul(prior_chol,step)
        max_du = maxval(abs(du))
        alpha = 1.d0
        if(max_du.gt.8.d0)alpha = 8.d0/max_du
        new_log = -huge(1.d0)
        do line_iter=1,40
            z_trial = z+alpha*step
            call posterior_state(n,prior_chol,delta,hazard,z_trial,new_log,&
                trial_gradient,trial_hessian,local_status)
            if(local_status.eq.0)then
                if(new_log.ge.old_log+1.d-4*alpha*directional)exit
            end if
            alpha = 0.5d0*alpha
        end do
        if(line_iter.gt.40)then
            status = -2
            return
        end if

        relative_step = maxval(abs(alpha*step))/(1.d0+maxval(abs(z_trial)))
        z = z_trial
        iterations = iter
        grad_norm = maxval(abs(trial_gradient))
        if(grad_norm.le.adaptive_mode_tol*(1.d0+maxval(abs(z))).and.&
            relative_step.le.dsqrt(adaptive_mode_tol))then
            status = 1
            return
        end if
    end do
    status = -1
end subroutine safeguarded_newton


subroutine posterior_state(n,prior_chol,delta,hazard,z,log_value,gradient,&
    hessian,status)
    implicit none
    integer,intent(in)::n
    double precision,dimension(n,n),intent(in)::prior_chol
    double precision,dimension(n),intent(in)::delta,hazard,z
    double precision,intent(out)::log_value
    double precision,dimension(n),intent(out)::gradient
    double precision,dimension(n,n),intent(out)::hessian
    integer,intent(out)::status

    integer::i,j,k
    double precision,dimension(n)::u,risk

    status = 0
    do i=1,n
        if(hazard(i).lt.-1.d-10.or.hazard(i).ne.hazard(i))then
            status = -1
            return
        end if
    end do

    u = matmul(prior_chol,z)
    do i=1,n
        if(hazard(i).le.0.d0)then
            risk(i) = 0.d0
        else if(u(i).gt.dlog(huge(1.d0))-dlog(hazard(i))-2.d0)then
            status = -1
            return
        else
            risk(i) = hazard(i)*dexp(u(i))
        end if
    end do

    log_value = dot_product(delta,u)-sum(risk)-0.5d0*dot_product(z,z)
    if(log_value.ne.log_value)then
        status = -1
        return
    end if
    gradient = matmul(transpose(prior_chol),delta-risk)-z

    hessian = 0.d0
    do j=1,n
        hessian(j,j) = 1.d0
    end do
    do i=1,n
        do j=1,n
            do k=1,n
                hessian(j,k) = hessian(j,k)+prior_chol(i,j)*risk(i)*&
                    prior_chol(i,k)
            end do
        end do
    end do
    hessian = 0.5d0*(hessian+transpose(hessian))
end subroutine posterior_state


subroutine stabilize_hessian(n,hessian,floor_eps,chol,eigenvectors,&
    eigenvalues,min_original,min_stabilized,condition,regularized,status)
    implicit none
    integer,intent(in)::n
    double precision,dimension(n,n),intent(inout)::hessian
    double precision,intent(in)::floor_eps
    double precision,dimension(n,n),intent(out)::chol,eigenvectors
    double precision,dimension(n),intent(out)::eigenvalues
    double precision,intent(out)::min_original,min_stabilized,condition
    integer,intent(out)::regularized,status

    integer::info,lwork,i,j,k
    double precision::lambda_max,lambda_floor
    double precision,dimension(:),allocatable::work

    hessian = 0.5d0*(hessian+transpose(hessian))
    eigenvectors = hessian
    lwork = max(1,3*n-1)
    allocate(work(lwork))
    call dsyev('V','U',n,eigenvectors,n,eigenvalues,work,lwork,info)
    deallocate(work)
    if(info.ne.0)then
        status = -1
        return
    end if

    min_original = eigenvalues(1)
    lambda_max = eigenvalues(n)
    call cholesky_lower(n,hessian,chol,info)
    regularized = 0
    if(info.ne.0)then
        regularized = 1
        lambda_floor = max(floor_eps*max(lambda_max,1.d-300),1.d-300)
        do i=1,n
            eigenvalues(i) = max(eigenvalues(i),lambda_floor)
        end do
        hessian = 0.d0
        do i=1,n
            do j=1,n
                do k=1,n
                    hessian(i,j) = hessian(i,j)+eigenvectors(i,k)*&
                        eigenvalues(k)*eigenvectors(j,k)
                end do
            end do
        end do
        hessian = 0.5d0*(hessian+transpose(hessian))
        call cholesky_lower(n,hessian,chol,info)
        if(info.ne.0)then
            status = -2
            return
        end if
    end if
    min_stabilized = minval(eigenvalues)
    condition = maxval(eigenvalues)/max(min_stabilized,1.d-300)
    status = 0
end subroutine stabilize_hessian


subroutine cholesky_lower(n,a,l,status)
    implicit none
    integer,intent(in)::n
    double precision,dimension(n,n),intent(in)::a
    double precision,dimension(n,n),intent(out)::l
    integer,intent(out)::status
    integer::i,j,k
    double precision::value

    l = 0.d0
    status = 0
    do i=1,n
        do j=1,i
            value = a(i,j)
            do k=1,j-1
                value = value-l(i,k)*l(j,k)
            end do
            if(i.eq.j)then
                if(value.le.0.d0.or.value.ne.value)then
                    status = i
                    return
                end if
                l(i,j) = dsqrt(value)
            else
                l(i,j) = value/l(j,j)
            end if
        end do
    end do
end subroutine cholesky_lower


subroutine solve_lower(l,x,status)
    implicit none
    double precision,dimension(:,:),intent(in)::l
    double precision,dimension(:),intent(inout)::x
    integer,intent(out)::status
    integer::i,j,n
    n = size(x)
    status = 0
    do i=1,n
        if(l(i,i).le.0.d0)then
            status = i
            return
        end if
        do j=1,i-1
            x(i) = x(i)-l(i,j)*x(j)
        end do
        x(i) = x(i)/l(i,i)
    end do
end subroutine solve_lower


subroutine solve_cholesky(n,l,rhs,solution)
    implicit none
    integer,intent(in)::n
    double precision,dimension(n,n),intent(in)::l
    double precision,dimension(n),intent(in)::rhs
    double precision,dimension(n),intent(out)::solution
    integer::i,j

    solution = rhs
    do i=1,n
        do j=1,i-1
            solution(i) = solution(i)-l(i,j)*solution(j)
        end do
        solution(i) = solution(i)/l(i,i)
    end do
    do i=n,1,-1
        do j=i+1,n
            solution(i) = solution(i)-l(j,i)*solution(j)
        end do
        solution(i) = solution(i)/l(i,i)
    end do
end subroutine solve_cholesky


subroutine adaptive_node_transform(n,x,z,status)
    implicit none
    integer,intent(in)::n
    double precision,dimension(n),intent(in)::x
    double precision,dimension(n),intent(out)::z
    integer,intent(out)::status
    integer::i,j

    status = 0
    z = posterior_proposal_scale*x
    ! Solve transpose(C) y = x, with C the lower Cholesky factor of H.
    do i=n,1,-1
        if(posterior_chol_h(i,i).le.0.d0)then
            status = i
            return
        end if
        do j=i+1,n
            z(i) = z(i)-posterior_chol_h(j,i)*z(j)
        end do
        z(i) = z(i)/posterior_chol_h(i,i)
    end do
    z = posterior_mode_z+z
end subroutine adaptive_node_transform


subroutine set_posterior_proposal_scale(n,scale)
    implicit none
    integer,intent(in)::n
    double precision,intent(in)::scale
    integer::i
    posterior_proposal_scale = max(scale,1.d-6)
    posterior_log_jacobian = dble(n)*dlog(posterior_proposal_scale)
    if(allocated(posterior_chol_h))then
        do i=1,n
            posterior_log_jacobian = posterior_log_jacobian-&
                dlog(posterior_chol_h(i,i))
        end do
    end if
end subroutine set_posterior_proposal_scale


subroutine store_mode_diagnostics(family,dimension)
    implicit none
    integer,intent(in)::family,dimension
    if(.not.allocated(adaptive_diag_i))return
    adaptive_diag_i(family,1) = family
    adaptive_diag_i(family,2) = dimension
    adaptive_diag_i(family,3) = adaptive_method
    adaptive_diag_i(family,4) = posterior_mode_status
    adaptive_diag_i(family,5) = posterior_mode_iterations
    adaptive_diag_i(family,6) = posterior_hessian_regularized
    adaptive_diag_d(family,1) = posterior_gradient_norm
    adaptive_diag_d(family,2) = posterior_condition_number
    adaptive_diag_d(family,3) = posterior_min_eigenvalue
    adaptive_diag_d(family,4) = posterior_min_stabilized
end subroutine store_mode_diagnostics


subroutine store_quadrature_diagnostics(family,method_used,rule,neval,ifail,&
    fallback_used,quadrature_failed,reached_cap,floor_used,abserr,relerr,&
    log_integral)
    use cubature_settings,only:cubature_signed_sum_status,&
        cubature_log_positive,cubature_log_negative,cubature_cancellation_ratio
    implicit none
    integer,intent(in)::family,method_used,rule,neval,ifail,fallback_used
    integer,intent(in)::quadrature_failed,reached_cap,floor_used
    double precision,intent(in)::abserr,relerr,log_integral
    if(.not.allocated(adaptive_diag_i))return
    adaptive_diag_i(family,3) = method_used
    adaptive_diag_i(family,7) = rule
    adaptive_diag_i(family,8) = neval
    adaptive_diag_i(family,9) = ifail
    adaptive_diag_i(family,10) = fallback_used
    adaptive_diag_i(family,11) = reached_cap
    adaptive_diag_i(family,12) = quadrature_failed
    adaptive_diag_i(family,14) = cubature_signed_sum_status
    adaptive_diag_i(family,15) = floor_used
    if(allocated(adaptive_floor_count))then
        adaptive_diag_i(family,16) = adaptive_floor_count(family)
    end if
    if(allocated(adaptive_floor_reason_count))then
        adaptive_diag_i(family,17:20) = adaptive_floor_reason_count(family,1:4)
    end if
    adaptive_diag_d(family,5) = abserr
    adaptive_diag_d(family,6) = relerr
    adaptive_diag_d(family,7) = log_integral
    adaptive_diag_d(family,8) = cubature_log_positive
    adaptive_diag_d(family,9) = cubature_log_negative
    adaptive_diag_d(family,10) = cubature_cancellation_ratio
end subroutine store_quadrature_diagnostics


subroutine record_integration_floor(family,reason)
    implicit none
    integer,intent(in)::family,reason
    if(allocated(adaptive_floor_count))then
        adaptive_floor_count(family) = adaptive_floor_count(family)+1
    end if
    if(allocated(adaptive_floor_reason_count))then
        if(reason.ge.1.and.reason.le.4)then
            adaptive_floor_reason_count(family,reason) = &
                adaptive_floor_reason_count(family,reason)+1
        end if
    end if
end subroutine record_integration_floor


subroutine store_covariance_diagnostics(family,status)
    implicit none
    integer,intent(in)::family,status
    if(.not.allocated(adaptive_diag_i))return
    adaptive_diag_i(family,13) = status
end subroutine store_covariance_diagnostics


subroutine initialize_signed_gk(n)
    implicit none
    integer,intent(in)::n
    if(allocated(gk_log_positive))deallocate(gk_log_positive)
    if(allocated(gk_log_negative))deallocate(gk_log_negative)
    allocate(gk_log_positive(n),gk_log_negative(n))
    gk_log_positive = -huge(1.d0)
    gk_log_negative = -huge(1.d0)
end subroutine initialize_signed_gk


double precision function logspace_add(a,b)
    implicit none
    double precision,intent(in)::a,b
    double precision::m
    if(a.le.-0.5d0*huge(1.d0))then
        logspace_add = b
    else if(b.le.-0.5d0*huge(1.d0))then
        logspace_add = a
    else
        m = max(a,b)
        logspace_add = m+dlog(dexp(a-m)+dexp(b-m))
    end if
end function logspace_add


subroutine signed_log_value(logpos,logneg,value,status)
    implicit none
    double precision,intent(in)::logpos,logneg
    double precision,intent(out)::value
    integer,intent(out)::status
    double precision::logabs

    if(logpos.le.logneg)then
        value = 0.d0
        status = -1
        return
    end if
    if(logneg.le.-0.5d0*huge(1.d0))then
        logabs = logpos
    else
        logabs = logpos+dlog(1.d0-dexp(logneg-logpos))
    end if
    if(logabs.gt.dlog(huge(1.d0)))then
        value = huge(1.d0)
        status = -2
    else
        value = dexp(logabs)
        status = 0
    end if
end subroutine signed_log_value

end module posterior_adaptive_settings
