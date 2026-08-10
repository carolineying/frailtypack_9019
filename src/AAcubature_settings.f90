module cubature_settings
    implicit none
    integer,save :: cubature_minpts,cubature_maxpts,cubature_calls
    integer,save :: cubature_failures,cubature_max_neval
    integer,save :: cubature_max_dimension,cubature_last_rule
    integer,save :: cubature_signed_sum_status
    integer,save :: cubature_cap_blocked
    double precision,save :: cubature_epsabs,cubature_epsrel
    double precision,save :: cubature_max_abserr,cubature_max_relerr
    double precision,save :: cubature_log_positive,cubature_log_negative
    double precision,save :: cubature_cancellation_ratio

contains

    logical function cubature_result_is_finite_positive(value,error)
        implicit none
        double precision,intent(in)::value,error
        cubature_result_is_finite_positive = value.eq.value.and.&
            error.eq.error.and.value.gt.0.d0.and.value.le.huge(1.d0).and.&
            error.ge.0.d0.and.error.le.huge(1.d0)
    end function cubature_result_is_finite_positive
end module cubature_settings
