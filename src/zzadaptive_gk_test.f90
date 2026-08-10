module adaptive_gk_test_state
    implicit none
    integer::test_case=1,test_method=0,test_dimension=0
    double precision::test_ridge=0.d0
    double precision,dimension(:),allocatable::test_shift
    integer::matrix_test_case=1
    double precision::matrix_test_constant=1.d0
    double precision,dimension(:),allocatable::matrix_target_mean
    double precision,dimension(:,:),allocatable::matrix_target_precision_chol
    double precision,dimension(:),allocatable::survival_delta,survival_hazard
    double precision,dimension(:,:),allocatable::survival_prior_chol
end module adaptive_gk_test_state


subroutine adaptive_gk_test(ndim,testcase,method,minpts,maxpts,epsabs,epsrel,&
    shift,ridge,result,abserr,neval,ifail,analytic,status)
    use lois_normales,only:hrmsym
    use adaptive_gk_test_state
    use cubature_settings,only:cubature_signed_sum_status,&
        cubature_result_is_finite_positive
    use posterior_adaptive_settings
    implicit none
    integer,intent(in)::ndim,testcase,method,minpts,maxpts
    integer,intent(out)::neval,ifail,status
    double precision,intent(in)::epsabs,epsrel,ridge
    double precision,dimension(ndim),intent(in)::shift
    double precision,intent(out)::result,abserr,analytic
    double precision,dimension(1)::values,errors
    double precision,dimension(:),allocatable::work
    integer,dimension(6)::icontrol
    double precision,dimension(2)::dcontrol
    double precision::precision,log_scale
    integer::i,restar
    external::adaptive_gk_test_callback

    status = 0
    if(ndim.lt.1.or.ridge.le.-1.d0.or.&
        (method.lt.0.or.method.gt.4))then
        status = -1
        result = 0.d0
        abserr = huge(1.d0)
        analytic = 0.d0
        neval = 0
        ifail = 1
        return
    end if

    test_case = testcase
    test_method = method
    test_dimension = ndim
    test_ridge = ridge
    if(allocated(test_shift))deallocate(test_shift)
    allocate(test_shift(ndim))
    if(testcase.eq.1)then
        test_shift = 0.d0
        test_ridge = 0.d0
    else
        test_shift = shift
    end if

    precision = 1.d0+test_ridge
    analytic = precision**(-0.5d0*dble(ndim))*&
        dexp(0.5d0*sum(test_shift*test_shift)/precision)
    log_scale = 0.d0
    posterior_adaptive_active = 0
    posterior_adaptive_fallback_active = 0

    if(method.eq.1.or.method.eq.3.or.method.eq.4)then
        icontrol = (/1,1,0,0,0,0/)
        dcontrol = (/1.d-10,1.d-10/)
        call initialize_posterior_adaptive(ndim,1,icontrol,dcontrol)
        if(allocated(posterior_mode_z))deallocate(posterior_mode_z)
        if(allocated(posterior_chol_h))deallocate(posterior_chol_h)
        allocate(posterior_mode_z(ndim),posterior_chol_h(ndim,ndim))
        posterior_mode_z = test_shift/precision
        posterior_chol_h = 0.d0
        do i=1,ndim
            posterior_chol_h(i,i) = dsqrt(precision)
        end do
        posterior_log_mode = 0.5d0*sum(test_shift*posterior_mode_z)
        posterior_adaptive_active = 1
        if(method.eq.4)then
            call set_posterior_proposal_scale(ndim,1.5d0)
        else
            call set_posterior_proposal_scale(ndim,1.d0)
        end if
        log_scale = posterior_log_mode+posterior_log_jacobian
    else if(method.eq.2)then
        ! Exercise the production prior-centred fallback arithmetic without
        ! changing the directly requested legacy branch (method 0).
        posterior_adaptive_fallback_active = 1
    end if

    allocate(work(20000))
    work = 0.d0
    values = 0.d0
    errors = 0.d0
    restar = 0
    if(method.eq.3.and.ndim.gt.4)then
        call hrmsym(ndim,1,minpts,max(minpts,int(0.70d0*dble(maxpts))),&
            adaptive_gk_test_callback,epsabs,epsrel,restar,values,errors,&
            neval,ifail,work)
        i = neval
        restar = 1
        call hrmsym(ndim,1,min(minpts,maxpts-i),maxpts-i,&
            adaptive_gk_test_callback,epsabs,epsrel,restar,values,errors,&
            neval,ifail,work)
        neval = neval+i
    else
        call hrmsym(ndim,1,minpts,maxpts,adaptive_gk_test_callback,epsabs,epsrel,&
            restar,values,errors,neval,ifail,work)
    end if
    result = values(1)*dexp(log_scale)
    abserr = errors(1)*dexp(log_scale)
    if(method.ge.1)then
        status = cubature_signed_sum_status
    else if(.not.cubature_result_is_finite_positive(result,abserr))then
        ! This is the same validation used after a production legacy fallback.
        status = -3
    end if

    posterior_adaptive_active = 0
    posterior_adaptive_fallback_active = 0
    if(method.eq.1.or.method.eq.3.or.method.eq.4)&
        call finalize_posterior_adaptive()
    if(allocated(test_shift))deallocate(test_shift)
    deallocate(work)
end subroutine adaptive_gk_test


subroutine adaptive_gk_test_callback(ndim,x,nf,funvls)
    use adaptive_gk_test_state
    use posterior_adaptive_settings,only:posterior_adaptive_active,&
        posterior_callback_returns_log,posterior_log_mode,adaptive_node_transform
    implicit none
    integer,intent(in)::ndim,nf
    double precision,dimension(ndim),intent(in)::x
    double precision,dimension(nf),intent(out)::funvls
    double precision,dimension(ndim)::b
    double precision::log_value
    integer::transform_status

    if(posterior_adaptive_active.eq.1)then
        call adaptive_node_transform(ndim,x,b,transform_status)
        if(transform_status.ne.0)then
            funvls = -huge(1.d0)
            return
        end if
        log_value = sum(test_shift*b)-0.5d0*test_ridge*sum(b*b)-&
            0.5d0*sum(b*b)-posterior_log_mode+0.5d0*sum(x*x)
    else
        b = x
        log_value = sum(test_shift*b)-0.5d0*test_ridge*sum(b*b)
    end if

    if(posterior_callback_returns_log.eq.1)then
        funvls(1) = log_value
    else
        funvls(1) = dexp(log_value)
    end if
end subroutine adaptive_gk_test_callback


subroutine adaptive_gk_matrix_test(ndim,testcase,method,minpts,maxpts,epsabs,&
    epsrel,target_mean,target_precision_chol,proposal_mean,&
    proposal_precision_chol,proposal_scale,constant,result,abserr,neval,&
    ifail,status,logpositive,lognegative)
    use lois_normales,only:hrmsym
    use adaptive_gk_test_state
    use cubature_settings,only:cubature_signed_sum_status
    use posterior_adaptive_settings
    implicit none
    integer,intent(in)::ndim,testcase,method,minpts,maxpts
    integer,intent(out)::neval,ifail,status
    double precision,intent(in)::epsabs,epsrel,proposal_scale,constant
    double precision,dimension(ndim),intent(in)::target_mean,proposal_mean
    double precision,dimension(ndim,ndim),intent(in)::target_precision_chol,&
        proposal_precision_chol
    double precision,intent(out)::result,abserr,logpositive,lognegative
    integer,dimension(6)::icontrol
    double precision,dimension(2)::dcontrol
    double precision,dimension(1)::values,errors
    double precision,dimension(:),allocatable::work
    double precision::log_scale
    integer::restar
    external::adaptive_gk_matrix_callback

    status = 0
    result = 0.d0
    abserr = huge(1.d0)
    logpositive = -huge(1.d0)
    lognegative = -huge(1.d0)
    neval = 0
    ifail = 1
    if(ndim.lt.1.or.(testcase.ne.1.and.testcase.ne.2).or.&
        (method.ne.0.and.method.ne.1).or.proposal_scale.le.0.d0.or.&
        constant.le.0.d0)then
        status = -1
        return
    end if

    matrix_test_case = testcase
    matrix_test_constant = constant
    if(allocated(matrix_target_mean))deallocate(matrix_target_mean)
    if(allocated(matrix_target_precision_chol))&
        deallocate(matrix_target_precision_chol)
    allocate(matrix_target_mean(ndim),matrix_target_precision_chol(ndim,ndim))
    matrix_target_mean = target_mean
    matrix_target_precision_chol = target_precision_chol

    posterior_adaptive_active = 0
    posterior_adaptive_fallback_active = 0
    log_scale = 0.d0
    if(method.eq.1)then
        icontrol = (/1,1,0,0,0,0/)
        dcontrol = (/1.d-10,1.d-10/)
        call initialize_posterior_adaptive(ndim,1,icontrol,dcontrol)
        allocate(posterior_mode_z(ndim),posterior_chol_h(ndim,ndim))
        posterior_mode_z = proposal_mean
        posterior_chol_h = proposal_precision_chol
        posterior_log_mode = 0.d0
        posterior_adaptive_active = 1
        call set_posterior_proposal_scale(ndim,proposal_scale)
        log_scale = posterior_log_jacobian
    end if

    allocate(work(20050))
    work = 0.d0
    values = 0.d0
    errors = 0.d0
    restar = 0
    call hrmsym(ndim,1,minpts,maxpts,adaptive_gk_matrix_callback,epsabs,&
        epsrel,restar,values,errors,neval,ifail,work)
    result = values(1)*dexp(log_scale)
    abserr = errors(1)*dexp(log_scale)
    if(method.eq.1)then
        status = cubature_signed_sum_status
        if(allocated(gk_log_positive))logpositive = gk_log_positive(1)+log_scale
        if(allocated(gk_log_negative))lognegative = gk_log_negative(1)+log_scale
        call finalize_posterior_adaptive()
    end if

    posterior_adaptive_active = 0
    posterior_adaptive_fallback_active = 0
    if(allocated(matrix_target_mean))deallocate(matrix_target_mean)
    if(allocated(matrix_target_precision_chol))&
        deallocate(matrix_target_precision_chol)
    deallocate(work)
end subroutine adaptive_gk_matrix_test


subroutine adaptive_gk_matrix_callback(ndim,x,nf,funvls)
    use adaptive_gk_test_state
    use posterior_adaptive_settings,only:posterior_adaptive_active,&
        posterior_callback_returns_log,adaptive_node_transform
    implicit none
    integer,intent(in)::ndim,nf
    double precision,dimension(ndim),intent(in)::x
    double precision,dimension(nf),intent(out)::funvls
    double precision,dimension(ndim)::z,difference,rotated
    double precision::log_likelihood,log_complete,logdet_precision
    integer::i,transform_status

    if(posterior_adaptive_active.eq.1)then
        call adaptive_node_transform(ndim,x,z,transform_status)
        if(transform_status.ne.0)then
            funvls = -huge(1.d0)
            return
        end if
    else
        z = x
    end if

    if(matrix_test_case.eq.1)then
        log_likelihood = 0.d0
        log_complete = -0.5d0*dot_product(z,z)
    else
        difference = z-matrix_target_mean
        rotated = matmul(transpose(matrix_target_precision_chol),difference)
        logdet_precision = 0.d0
        do i=1,ndim
            logdet_precision = logdet_precision+&
                2.d0*dlog(matrix_target_precision_chol(i,i))
        end do
        log_complete = dlog(matrix_test_constant)+0.5d0*logdet_precision-&
            0.5d0*dot_product(rotated,rotated)
        log_likelihood = log_complete+0.5d0*dot_product(z,z)
    end if

    if(posterior_adaptive_active.eq.1)then
        funvls(1) = log_complete+0.5d0*dot_product(x,x)
    else if(posterior_callback_returns_log.eq.1)then
        funvls(1) = log_likelihood
    else
        funvls(1) = dexp(log_likelihood)
    end if
end subroutine adaptive_gk_matrix_callback


subroutine adaptive_mode_test(ndim,prior_chol,delta,hazard,mode_tol,maxit,&
    hessian_eps,mode_z,log_mode,gradient,hessian,chol_h,status,iterations,&
    gradnorm,mineig,minstab,condition,regularized)
    use posterior_adaptive_settings
    implicit none
    integer,intent(in)::ndim,maxit
    integer,intent(out)::status,iterations,regularized
    double precision,intent(in)::mode_tol,hessian_eps
    double precision,dimension(ndim,ndim),intent(in)::prior_chol
    double precision,dimension(ndim),intent(in)::delta,hazard
    double precision,dimension(ndim),intent(out)::mode_z,gradient
    double precision,dimension(ndim,ndim),intent(out)::hessian,chol_h
    double precision,intent(out)::log_mode,gradnorm,mineig,minstab,condition
    integer,dimension(6)::icontrol
    double precision,dimension(2)::dcontrol
    integer::state_status

    icontrol = (/1,maxit,0,0,0,0/)
    dcontrol = (/mode_tol,hessian_eps/)
    call initialize_posterior_adaptive(ndim,1,icontrol,dcontrol)
    call setup_posterior_adaptive(1,1,ndim,prior_chol,delta,hazard,status)
    iterations = posterior_mode_iterations
    gradnorm = posterior_gradient_norm
    mineig = posterior_min_eigenvalue
    minstab = posterior_min_stabilized
    condition = posterior_condition_number
    regularized = posterior_hessian_regularized
    if(posterior_adaptive_active.eq.1)then
        mode_z = posterior_mode_z
        log_mode = posterior_log_mode
        chol_h = posterior_chol_h
        call posterior_state(ndim,prior_chol,delta,hazard,mode_z,log_mode,&
            gradient,hessian,state_status)
        if(state_status.ne.0)status = -10
    else
        mode_z = 0.d0
        gradient = huge(1.d0)
        hessian = 0.d0
        chol_h = 0.d0
        log_mode = -huge(1.d0)
    end if
    call finalize_posterior_adaptive()
end subroutine adaptive_mode_test


subroutine adaptive_survival_test(ndim,method,minpts,maxpts,epsabs,epsrel,&
    prior_chol,delta,hazard,proposal_scale,proposal_shift,result,abserr,&
    neval,ifail,status,laplace,mode_z,gradient,hessian,chol_h,iterations,&
    gradnorm,mineig,minstab,condition,regularized)
    use lois_normales,only:hrmsym
    use adaptive_gk_test_state
    use posterior_adaptive_settings
    implicit none
    integer,intent(in)::ndim,method,minpts,maxpts
    integer,intent(out)::neval,ifail,status,iterations,regularized
    double precision,intent(in)::epsabs,epsrel,proposal_scale
    double precision,dimension(ndim,ndim),intent(in)::prior_chol
    double precision,dimension(ndim),intent(in)::delta,hazard,proposal_shift
    double precision,intent(out)::result,abserr,laplace,gradnorm,mineig,&
        minstab,condition
    double precision,dimension(ndim),intent(out)::mode_z,gradient
    double precision,dimension(ndim,ndim),intent(out)::hessian,chol_h
    integer,dimension(6)::icontrol
    double precision,dimension(2)::dcontrol
    double precision,dimension(1)::values,errors
    double precision,dimension(:),allocatable::work
    double precision::log_scale,log_mode_copy
    integer::restar,state_status,i
    external::adaptive_survival_callback

    result = 0.d0
    abserr = huge(1.d0)
    laplace = 0.d0
    neval = 0
    ifail = 1
    status = -1
    iterations = 0
    regularized = 0
    gradient = 0.d0
    hessian = 0.d0
    chol_h = 0.d0
    mode_z = 0.d0
    gradnorm = huge(1.d0)
    mineig = 0.d0
    minstab = 0.d0
    condition = 0.d0
    if(ndim.lt.1.or.(method.ne.0.and.method.ne.1).or.&
        proposal_scale.le.0.d0)return

    if(allocated(survival_delta))deallocate(survival_delta)
    if(allocated(survival_hazard))deallocate(survival_hazard)
    if(allocated(survival_prior_chol))deallocate(survival_prior_chol)
    allocate(survival_delta(ndim),survival_hazard(ndim),&
        survival_prior_chol(ndim,ndim))
    survival_delta = delta
    survival_hazard = hazard
    survival_prior_chol = prior_chol

    posterior_adaptive_active = 0
    posterior_adaptive_fallback_active = 0
    log_scale = 0.d0
    if(method.eq.1)then
        icontrol = (/1,100,0,0,0,0/)
        dcontrol = (/1.d-11,1.d-10/)
        call initialize_posterior_adaptive(ndim,1,icontrol,dcontrol)
        call setup_posterior_adaptive(1,1,ndim,prior_chol,delta,hazard,status)
        if(posterior_adaptive_active.ne.1)goto 900
        mode_z = posterior_mode_z
        log_mode_copy = posterior_log_mode
        chol_h = posterior_chol_h
        iterations = posterior_mode_iterations
        gradnorm = posterior_gradient_norm
        mineig = posterior_min_eigenvalue
        minstab = posterior_min_stabilized
        condition = posterior_condition_number
        regularized = posterior_hessian_regularized
        call posterior_state(ndim,prior_chol,delta,hazard,mode_z,log_mode_copy,&
            gradient,hessian,state_status)
        if(state_status.ne.0)then
            status = -10
            goto 900
        end if
        laplace = dexp(log_mode_copy)
        do i=1,ndim
            laplace = laplace/chol_h(i,i)
        end do
        posterior_mode_z = posterior_mode_z+proposal_shift
        call set_posterior_proposal_scale(ndim,proposal_scale)
        log_scale = posterior_log_mode+posterior_log_jacobian
    else
        status = 0
    end if

    allocate(work(20050))
    work = 0.d0
    values = 0.d0
    errors = 0.d0
    restar = 0
    call hrmsym(ndim,1,minpts,maxpts,adaptive_survival_callback,epsabs,&
        epsrel,restar,values,errors,neval,ifail,work)
    result = values(1)*dexp(log_scale)
    abserr = errors(1)*dexp(log_scale)
    deallocate(work)

900 continue
    posterior_adaptive_active = 0
    posterior_adaptive_fallback_active = 0
    if(method.eq.1)call finalize_posterior_adaptive()
    if(allocated(survival_delta))deallocate(survival_delta)
    if(allocated(survival_hazard))deallocate(survival_hazard)
    if(allocated(survival_prior_chol))deallocate(survival_prior_chol)
end subroutine adaptive_survival_test


subroutine adaptive_survival_callback(ndim,x,nf,funvls)
    use adaptive_gk_test_state
    use posterior_adaptive_settings,only:posterior_adaptive_active,&
        posterior_callback_returns_log,posterior_log_mode,adaptive_node_transform
    implicit none
    integer,intent(in)::ndim,nf
    double precision,dimension(ndim),intent(in)::x
    double precision,dimension(nf),intent(out)::funvls
    double precision,dimension(ndim)::z,u
    double precision::log_likelihood,log_correction
    integer::transform_status

    if(posterior_adaptive_active.eq.1)then
        call adaptive_node_transform(ndim,x,z,transform_status)
        if(transform_status.ne.0)then
            funvls = -huge(1.d0)
            return
        end if
    else
        z = x
    end if
    u = matmul(survival_prior_chol,z)
    log_likelihood = dot_product(survival_delta,u)-&
        sum(survival_hazard*dexp(u))
    if(posterior_adaptive_active.eq.1)then
        log_correction = log_likelihood-0.5d0*dot_product(z,z)-&
            posterior_log_mode+0.5d0*dot_product(x,x)
    else
        log_correction = log_likelihood
    end if
    if(posterior_callback_returns_log.eq.1)then
        funvls(1) = log_correction
    else
        funvls(1) = dexp(log_correction)
    end if
end subroutine adaptive_survival_callback


subroutine ascertainment_gaussian_test(cumulative_hazard,variance,result)
    use donnees,only:x2,w2
    implicit none
    double precision,intent(in)::cumulative_hazard,variance
    double precision,intent(out)::result
    double precision,parameter::pi=3.141592653589793d0
    integer::j

    result = 0.d0
    do j=1,20
        result = result+w2(j)*dexp(-cumulative_hazard*&
            dexp(dsqrt(variance)*x2(j))-0.5d0*x2(j)*x2(j))/dsqrt(2.d0*pi)
    end do
end subroutine ascertainment_gaussian_test
