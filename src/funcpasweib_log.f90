

!========================          FUNCPA_WEIB          ====================
    double precision function funcpasweib_log(b,np,id,thi,jd,thj,k0)
    use tailles
    !use comon,only:etaR,etaD,betaR,betaD
    use comon,only:t0,t1,c,nsujet,nva, &
    nst,stra,ve,effet,ng,g,nig,AG,kkapa,sig2,sig22, &
    indictronq,auxig,res3,res5,res6,sig_prob,sig_prob2,  &
    etaT,betaT,res2_ind,res3_ind,res1,k0T,Kmat1,Kmat2, corrRE, Ut,Ut_det,Ut22,&
    Ut_table,nmemb,currentage,proband,agemin,sig_fixed
    use residusM
    use lois_normales
    use donnees_indiv,only:it_rec
    use optim
    use cubature_settings
    use posterior_adaptive_settings
    
    implicit none

! *** NOUVELLLE DECLARATION F90 :

    integer::nb,np,id,jd,i,j,k,cptg,ig,choix,ii,jj,integration_status
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,dnb,res,vet,int,int_safe,finddet,eps,slogasc
    double precision,dimension(np)::b,bh
    double precision,dimension(ngmax)::res2!res1,
    double precision,dimension(2)::k0
    double precision,dimension(ngmax)::integrale1,integrale2,integrale3,integrale
    double precision,dimension(ngmax)::log_integrale
    double precision::family_abserr
    double precision,parameter::pi=3.141592653589793d0
    double precision,parameter::eps_int=1.d-12,tol_int=1.d-8

    ! for the numerical integral hrmsym
        integer :: restar,nf,ier
        double precision:: epsabs,epsrel
        double precision,dimension(2):: result, abserr2
        double precision,dimension(1000) :: work
        external :: vraistot
        integer ::neval,ifail
        
    kkapa=k0 ! inutile dans le calcul mais est quand meme argument de la fonction donc doit etre present
    j=0
    eps=1.d-12

    ! These module work arrays must be local to one likelihood evaluation.
    ! A hard family-integral status can jump to label 123 before the normal
    ! per-family deallocation.  Clean stale state defensively on entry and at
    ! the common exit so consecutive frailtyPenal() calls in one R process are
    ! safe.
    call cleanup_correlated_weib_work_arrays()


    bh=b

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    !bh(1) =   1.7769365120637861    
    !    bh(2) = 11.408506591304478  
            
        !        bh(3) = 0.45347120527691670


    ii=1 !en plus strates A.Lafourcade 05/2014
    do jj=1,nst
        betaT(jj)=bh(ii)**2
        etaT(jj)=bh(ii+1)**2
        ii=ii+2
    end do

    if(effet.eq.1.and.corrRE.le.1) then
        sig2 = bh(np-nva)*bh(np-nva)
    else if(corrRE.eq.2) then
        sig2 = sig_fixed!1.d0!
        sig22 = bh(np-nva)*bh(np-nva)!bh(np-nva-1)*bh(np-nva-1)!1.d0!
    else if(corrRE.eq.3) then 
        sig2 = sig_fixed
    endif

!-------------------------------------------------------
!--------- calcul de la vraisemblance ------------------
!-------------------------------------------------------

!--- avec ou sans variable explicative  ------cc

    res1 = 0.d0
    res2 = 0.d0
    res3 = 0.d0
    res2_ind = 0.d0
    res3_ind = 0.d0
    integrale = 1.d0
    integrale1 = 1.d0
    integrale2 = 1.d0
    integrale3 = 1.d0
    log_integrale = 0.d0
    res5 = 0.d0
    res6 = 0.d0
    sig_prob = 0.d0
    cpt = 0

!*******************************************
!---- sans effet aleatoire dans le modele
!*******************************************

    if (effet.eq.0) then
        do i=1,nsujet
            cpt(g(i))=cpt(g(i))+1

            if(nva.gt.0)then
                vet = 0.d0
                do j=1,nva
                    vet =vet + bh(np-nva+j)*dble(ve(i,j))
                end do
                vet = dexp(vet)
            else
                vet=1.d0
            endif

            if(c(i).eq.1)then !en plus strates A.Lafourcade 05/2014
                res2(g(i)) = res2(g(i))+(betaT(stra(i))-1.d0)*dlog(t1(i)-agemin)+ &
                dlog(betaT(stra(i)))-betaT(stra(i))*dlog(etaT(stra(i)))+dlog(vet)
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpasweib_log=-1.d9
                goto 123
            end if
           !en plus strates A.Lafourcade 05/2014
            res1(g(i)) = res1(g(i)) + (((t1(i)-agemin)/etaT(stra(i)))**betaT(stra(i)))*vet -&
            ((t0(i)/etaT(stra(i)))**betaT(stra(i)))*vet
            RisqCumul(i) = (((t1(i)-agemin)/etaT(stra(i)))**betaT(stra(i)))*vet

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpasweib_log=-1.d9
                goto 123
            end if
        end do

        res = 0.d0
        cptg = 0

! k indice les groupes
        do k=1,ng
            if(cpt(k).gt.0)then
                nb = nig(k)
                dnb = dble(nig(k))
                res = res-res1(k)+res2(k)
                cptg = cptg + 1
                if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    funcpasweib_log=-1.d9
                    goto 123
                end if
            endif
        end do

!*********************************************
!----avec un effet aleatoire dans le modele
!*********************************************

    else

! i indice les sujets
        do i=1,nsujet

            cpt(g(i))=cpt(g(i))+1

            if(nva.gt.0)then
                vet = 0.d0
                do j=1,nva
                    vet =vet + bh(np-nva+j)*dble(ve(i,j))
                end do
                vet = dexp(vet)
            else
                vet=1.d0
            endif

            if(c(i).eq.1)then !en plus strates A.Lafourcade 05/2014
                res2(g(i)) = res2(g(i))+(betaT(stra(i))-1.d0)*dlog((t1(i)-agemin))+ &
                dlog(betaT(stra(i)))-betaT(stra(i))*dlog(etaT(stra(i)))+dlog(vet)
            endif
!write(*,*)res2(g(i)),betaT(stra(i)),dlog(t1(i)),etaT(stra(i)),dlog(vet)
            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge.1.d30)) then
                 funcpasweib_log=-1.d9
                 goto 123
            end if
!if(res2(1).ge.1.33772.and.res2(1).le.1.33773)write(*,*)i,res2(g(i)),betaT(stra(i)),dlog(t1(i)),etaT(stra(i)),dlog(vet)
 
! modification pour nouvelle vraisemblance / troncature:
            !en plus strates A.Lafourcade 05/2014
                res3_ind(i) = ((t0(i)/etaT(stra(i)))**betaT(stra(i)))*vet
          res3(g(i)) = res3(g(i)) + ((t0(i)/etaT(stra(i)))**betaT(stra(i)))*vet ! en plus
           res5(i) = (((t1(i)-agemin)/etaT(stra(i)))**betaT(stra(i)))*vet
            res1(g(i)) = res1(g(i)) + res5(i) ! pour les résidus
			
        res6(i) = ((currentage(i)-agemin)/etaT(stra(i)))**betaT(stra(i))*vet !cumulative hazard for probands
            
        
            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge.1.d30)) then
                print*,"here6"
                funcpasweib_log=-1.d9
                goto 123
            end if

            if ((res5(i).ne.res5(i)).or.(abs(res5(i)).ge.1.d30)) then
                print*,"here7",res5(i),i,t1(i),agemin,vet,((t1(i)-agemin)/etaT(stra(i))),&
                (((t1(i)-agemin)/etaT(stra(i)))**betaT(stra(i)))
                funcpasweib_log=-1.d9
                goto 123
            end if
        end do

!**************INTEGRALES ****************************
       if(corrRE.eq.0) then
    
            do ig=1,ng
                auxig = ig
                choix = 1
                call gauherS(int,choix)
                integrale1(ig) = int
                if (AG.eq.1) then
                    choix = 3
                    call gauherS(int,choix)
                    integrale3(ig) = int
    !                auxig = ig
    !            epsabs = 1.d-100
    !            epsrel = 1.d-100
    !            restar = 0
    !            nf = 1
    !            allocate(Ut(1,1))
    !            Ut = dsqrt(sig2)
    !                call  hrmsym(1, nf, 30, 500,vraistot, epsabs, &
    !                epsrel, restar, result, abserr2, neval, ifail, work)
    !    integrale(ig) =result(1)
    !        deallocate(ut)
                endif
                if (indictronq.eq.1) then
                    choix = 2
                    call gauherS(int,choix)
                    integrale2(ig) = int
                endif
            end do
        else if(corrRE.eq.1.or.corrRE.eq.3) then
            it_rec = 1
        !    det_all = 0.d0
            do ig = 1,ng
                auxig = ig
                epsabs = cubature_epsabs
                epsrel = cubature_epsrel
                restar = 0
                nf = 1
                allocate(Ut(nmemb(ig),nmemb(ig)),Ut_det(nmemb(ig),nmemb(ig)),Ut_table(nmemb(ig)*(nmemb(ig)+1)/2))
                Ut = 0.d0
                Ut_det = 0.d0
                do j=1,nmemb(ig)
                    do k=1,nmemb(ig)
                        Ut(j,k) =sig2* Kmat1(it_rec+j-1,it_rec+k-1)!sqrt(sig2)* 
                        Ut_det(j,k)=sig2* Kmat1(it_rec+j-1,it_rec+k-1)
                    end do
                end do
        
        !        write(*,*)ut 
        !        write(*,*)ut_det
        !        do j=1,nmemb(ig)
        !            do k=1,j
        !                Ut_det(k,j)=Ut(j,k)
        !            write(*,*)k,j,Ut_det(k,j),Ut_det(j,k),Ut(j,k),Ut(k,j)
        !            end do
        !        end do
        !    write(*,*)'ut',ut_det
        !        det_all(ig) = finddet(Ut_det,nmemb(ig))
                
        
                    Ut = 0.d0
                        do j=1,nmemb(ig)
                    do k=1,nmemb(ig)
                        Ut(j,k) =sig2* Kmat1(it_rec+j-1,it_rec+k-1)!
                        
                     if(j.eq.k)sig_prob(it_rec+j-1) = Ut(j,k)
                    end do
                end do
                jj = 0
                do j=1,nmemb(ig)
                    do k=j,nmemb(ig)
                        jj=j+k*(k-1)/2
                        Ut_table(jj)=Ut(j,k)
                    end do
                end do
        
                call dmfsdj(Ut_table, nmemb(ig),eps,ier)
                call store_covariance_diagnostics(ig,ier)
                if(ier.ne.0)then
                    call record_covariance_factorization_failure(ig,nmemb(ig),ier)
                    funcpasweib_log=-1.d9
                    goto 123
                end if
                
            
                Ut=0.d0
                do j=1,nmemb(ig)
                  do k=1,j
                     Ut(j,k)=Ut_table(k+j*(j-1)/2)
                  end do
                end do

                call integrate_correlated_family(ig,it_rec,nmemb(ig),&
                    log_integrale(ig),integrale(ig),family_abserr,neval,&
                    ifail,integration_status)
                if(integration_status.ne.0)then
                    funcpasweib_log=-1.d9
                    goto 123
                end if
            !    write(*,*)ig,integrale(ig)
                it_rec = it_rec + nmemb(ig)
                deallocate(Ut,Ut_det,Ut_table)
            end do
        else if(corrRE.eq.2) then 
        it_rec = 1
        !    det_all = 0.d0
            do ig = 1,ng
                auxig = ig
                epsabs = cubature_epsabs
                epsrel = cubature_epsrel
                restar = 0
                nf = 1
                allocate(Ut(nmemb(ig),nmemb(ig)),Ut22(nmemb(ig),nmemb(ig)),Ut_table(nmemb(ig)*(nmemb(ig)+1)/2))
                Ut = 0.d0
                Ut22 = 0.d0
                do j=1,nmemb(ig)
                    do k=1,nmemb(ig)
            !        if(ig.eq.9)write(*,*)j,k,sig2,Kmat1(it_rec+j-1,it_rec+k-1),sig22, Kmat2(it_rec+j-1,it_rec+k-1)
                        Ut(j,k) =sig2*Kmat1(it_rec+j-1,it_rec+k-1)+sig22* Kmat2(it_rec+j-1,it_rec+k-1)
                    Ut22(j,k) =sig22* Kmat2(it_rec+j-1,it_rec+k-1)
                    if(j.eq.k) then 
                    sig_prob(it_rec+j-1) =  Ut(j,k)!sig2*Kmat1(it_rec+j-1,it_rec+k-1)!
                    sig_prob2(it_rec+j-1) = Ut22(j,k)
                    end if
                end do
                end do

                jj = 0
                do j=1,nmemb(ig)
                    do k=j,nmemb(ig)
                        jj=j+k*(k-1)/2
                        Ut_table(jj)=Ut(j,k)
                    end do
                end do

                call dmfsdj(Ut_table, nmemb(ig),eps,ier)
                call store_covariance_diagnostics(ig,ier)
                if(ier.ne.0)then
                    call record_covariance_factorization_failure(ig,nmemb(ig),ier)
                    funcpasweib_log=-1.d9
                    goto 123
                end if

            
                Ut=0.d0
                do j=1,nmemb(ig)
                  do k=1,j
                     Ut(j,k)=Ut_table(k+j*(j-1)/2)
                  end do
                end do
                
                
                jj = 0
                do j=1,nmemb(ig)
                    do k=j,nmemb(ig)
                        jj=j+k*(k-1)/2
                        Ut_table(jj)=Ut22(j,k)
                    end do
                end do
        
                call dmfsdj(Ut_table, nmemb(ig),eps,ier)
                if(ier.ne.0)then
                    call store_covariance_diagnostics(ig,ier)
                    call record_covariance_factorization_failure(ig,nmemb(ig),ier)
                    funcpasweib_log=-1.d9
                    goto 123
                end if
                
            
                Ut22=0.d0
                do j=1,nmemb(ig)
                  do k=1,j
                     Ut22(j,k)=Ut_table(k+j*(j-1)/2)
                  end do
                end do
               
                call integrate_correlated_family(ig,it_rec,nmemb(ig),&
                    log_integrale(ig),integrale(ig),family_abserr,neval,&
                    ifail,integration_status)
                if(integration_status.ne.0)then
                    funcpasweib_log=-1.d9
                    goto 123
                end if
            !write(*,*)ig,integrale(ig),nmemb(ig)
            it_rec = it_rec + nmemb(ig)
                deallocate(Ut,Ut22,Ut_table)
            end do
        end if
!************* FIN INTEGRALES ************************
!stop
        res = 0.d0
        cptg = 0

!     gam2 = gamma(inv)
! k indice les groupes

        do k=1,ng
                if(corrRE.eq.0) then
!ccccc ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
                    if(AG.EQ.1)then
                        res = res+res2(k)-dlog(dsqrt(sig2))-dlog(2.d0*pi)/2.d0+ &
                        dlog(integrale3(k))
!ccccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    else
                        res = res+res2(k)+dlog(integrale1(k)) &
                        -dlog(integrale2(k))
                    endif
     !           write(*,*)k,res,res2(k), -dlog(dsqrt(sig2))-dlog(2.d0*pi)/2.d0,dlog(integrale3(k))
                else
                res = res+res2(k)+log_integrale(k)
    !    write(*,*)k,res,res2(k),dlog(integrale(k))
    !            write(*,*)k,res,res2(k), -nmemb(k)*dlog(2.d0*pi)/2.d0,dlog(integrale(k)),0.5d0*dlog(abs(det_all(k))),&
    !            det_all(k)
             !       write(*,*)k,res,res2(k),nmemb(k)*dlog(2.d0*pi)/2.d0,dlog(integrale(k))        
                end if
                !write(*,*)k,res,res2(k),dlog(integrale(k))
                   if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        print*,"here8",k,res,res2(k),integrale(k),bh**2.d0
                          funcpasweib_log=-1.d9
                          goto 123
                    end if
!     developpement de taylor d ordre 3
!                   write(*,*)'************** TAYLOR *************'
!cccc ancienne vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    !if(AG.EQ.1)then
                    !    res = res-dnb*dlog(theta*(res1(k)-res3(k))+1.d0) &
                    !    -(res1(k)-res3(k))*(1.d0-theta*(res1(k)-res3(k))/2.d0 &
                    !    +theta*theta*(res1(k)-res3(k))*(res1(k)-res3(k))/3.d0)+res2(k)+sum

!cccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    !else
                    !    res = res-dnb*dlog(theta*res1(k)+1.d0)-res1(k)*(1.d0-theta*res1(k) &
                    !    /2.d0+theta*theta*res1(k)*res1(k)/3.d0) &
                    !    +res2(k)+sum &
                    !    +res3(k)*(1.d0-theta*res3(k)/2.d0 &
                    !        +theta*theta*res3(k)*res3(k)/3.d0)
                    !endif
                       !if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                        !  funcpasweib_log=-1.d9
                         ! goto 123
                       !end if
                !endif
            !endif
        end do
    endif !fin boucle effet=0
!    stop
!if(res.ge.0.d0)then 
!write(*,*)res 
!write(*,*)res2 
!write(*,*)integrale 
!write(*,*)-0.5d0*dlog(abs(det_all(k)))-nmemb(k)*dlog(2.d0*pi)/2.d0 
!stop 
!end if

!    Changed JRG 25 May 05
    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
      funcpasweib_log=-1.d9
      goto 123
    end if
        
    !!!! Ascertainment correction
    slogasc  = 0.d0
    do i = 1,nsujet 
            if(proband(i).eq.1) then 
                auxig = i
                if(corrRE.ge.1) then
                    ! All correlated configurations use the proband marginal
                    ! variance already stored in sig_prob.  In particular,
                    ! corrRE=3 is the one-matrix, fixed-variance case.
                    int = -huge(1.d0)
                    call gauherASC(int,1,sig_prob(i))
                    if (int.ne.int.or.int.lt.-tol_int.or.&
                        int.gt.1.d0+tol_int) then
                    write(*,*) "SERIOUS ASCERTAINMENT INTEGRAL FAILURE"
                    write(*,*) "int = ", int
                    write(*,*) "Valid range is [0, 1]"
                    write(*,*) "i = ", i, " proband = ", proband(i), " currentage = ", currentage(i)
                    write(*,*) "sig_prob = ", sig_prob(i), " corrRE = ", corrRE, " sig2 = ", sig2
                    if(corrRE.eq.2) write(*,*) "sig_prob2 = ", sig_prob2(i), " sig22 = ", sig22
                    funcpasweib_log=-1.d9
                    goto 123
                end if
                int_safe = min(max(int, eps_int), 1.d0-eps_int)
                if (int_safe.ne.int.and.adaptive_diagnostics_enabled.eq.1) then
                    write(*,*) "CLIPPED ASCERTAINMENT int from ", int, " to ", int_safe
                    write(*,*) "i = ", i, " proband = ", proband(i), " currentage = ", currentage(i)
                    write(*,*) "sig_prob = ", sig_prob(i), " corrRE = ", corrRE, " sig2 = ", sig2
                    if(corrRE.eq.2) write(*,*) "sig_prob2 = ", sig_prob2(i), " sig22 = ", sig22
                end if
                if(t1(i).lt.currentage(i)) then 
                    slogasc = slogasc + dlog(1.d0-int_safe)
                else 
                    slogasc = slogasc + dlog(int_safe)
                end if
            else 
                if(t1(i).lt.currentage(i)) then
                    slogasc = slogasc + dlog(1-dexp(-res6(i)))
                else 
                    slogasc  = slogasc - res6(i)
                end if
            end if 
        end if
    end do
 !   write(*,*)res, slogasc
  !  stop
!write(*,*)sig2,res,res2(1),dlog(integrale(1)),res5(1),res3_ind(1)
  funcpasweib_log = res - slogasc

     if ((funcpasweib_log.ne.funcpasweib_log).or.(abs(funcpasweib_log).ge. 1.d30)) then
      funcpasweib_log=-1.d9
      goto 123
    end if
    do k=1,ng
        cumulhaz(k)=res1(k)
    end do
        do k = 1,nsujet 
        cumulhaz_corr(k) = res5(k) 
    end do

123     continue

    call cleanup_correlated_weib_work_arrays()

    return

    end function funcpasweib_log


subroutine cleanup_correlated_weib_work_arrays()
    use comon,only:Ut,Ut_det,Ut22,Ut_table
    implicit none

    if(allocated(Ut))deallocate(Ut)
    if(allocated(Ut_det))deallocate(Ut_det)
    if(allocated(Ut22))deallocate(Ut22)
    if(allocated(Ut_table))deallocate(Ut_table)
end subroutine cleanup_correlated_weib_work_arrays


subroutine record_covariance_factorization_failure(family,dimension,status)
    use cubature_settings
    use posterior_adaptive_settings
    implicit none
    integer,intent(in)::family,dimension,status

    call reset_posterior_family()
    cubature_signed_sum_status = -9
    cubature_log_positive = -huge(1.d0)
    cubature_log_negative = -huge(1.d0)
    cubature_cancellation_ratio = huge(1.d0)
    if(adaptive_method.eq.1)posterior_mode_status = -6
    posterior_gradient_norm = huge(1.d0)
    call store_mode_diagnostics(family,dimension)
    call store_covariance_diagnostics(family,status)
    call store_quadrature_diagnostics(family,adaptive_method,0,0,5,0,1,0,0,&
        huge(1.d0),huge(1.d0),-huge(1.d0))
    cubature_failures = cubature_failures+1
    cubature_max_dimension = max(cubature_max_dimension,dimension)
end subroutine record_covariance_factorization_failure


subroutine correlated_weib_cleanup_test(status)
    use comon,only:Ut,Ut_det,Ut22,Ut_table
    implicit none
    integer,intent(out)::status

    status = 0
    call cleanup_correlated_weib_work_arrays()
    allocate(Ut(2,2),Ut_det(2,2),Ut22(2,2),Ut_table(3))
    call cleanup_correlated_weib_work_arrays()
    if(allocated(Ut))status = status+1
    if(allocated(Ut_det))status = status+2
    if(allocated(Ut22))status = status+4
    if(allocated(Ut_table))status = status+8

    ! A second allocation proves that no stale allocation survives cleanup.
    allocate(Ut(1,1),Ut_det(1,1),Ut22(1,1),Ut_table(1))
    call cleanup_correlated_weib_work_arrays()
    if(allocated(Ut))status = status+16
    if(allocated(Ut_det))status = status+32
    if(allocated(Ut22))status = status+64
    if(allocated(Ut_table))status = status+128
end subroutine correlated_weib_cleanup_test


subroutine integrate_correlated_family(family,offset,n,log_integral,&
    integral_value,final_abserr,final_neval,final_ifail,integration_status)
    use comon,only:Ut,c,res5,res3_ind
    use lois_normales,only:hrmsym
    use cubature_settings
    use posterior_adaptive_settings

    implicit none
    integer,intent(in)::family,offset,n
    double precision,intent(out)::log_integral,integral_value,final_abserr
    integer,intent(out)::final_neval,final_ifail,integration_status

    integer::i,mode_status,restar,nf,method_used,fallback_used
    integer::quadrature_failed,reached_cap,call_neval,call_ifail
    integer::evaluations_used,remaining_pts,call_minpts,any_cap_blocked
    integer::stage_budget,stage_rule,stage_signed_status
    integer::floor_used,floor_reason
    double precision::epsabs,epsrel,log_scale,relative_error
    double precision::stage_result,stage_error,stage_log_positive
    double precision::stage_log_negative,stage_cancellation_ratio
    double precision,parameter::integration_floor_value=1.d-300
    double precision,parameter::integration_log_floor=dlog(1.d-300)
    logical::adaptive_result_ready
    double precision,dimension(2)::result,abserr
    ! HRMSYM documents up to 10,000 stored symmetric sums and requires
    ! 2*NF+(NF+1)*NUMSMS work elements.  The historical length 1,000 could
    ! overwrite adjacent state at higher rules; use the documented bound.
    double precision,dimension(20050)::work
    double precision,dimension(n)::delta,hazard
    external::vraistot

    do i=1,n
        delta(i) = dble(c(offset+i-1))
        hazard(i) = res5(offset+i-1)-res3_ind(offset+i-1)
    end do

    call reset_posterior_family()
    mode_status = 0
    if(adaptive_method.eq.1)then
        call setup_posterior_adaptive(family,offset,n,Ut,delta,hazard,mode_status)
    else
        call store_mode_diagnostics(family,n)
    end if

    epsabs = cubature_epsabs
    epsrel = cubature_epsrel
    nf = 1
    fallback_used = 0
    method_used = adaptive_method
    integration_status = 0
    evaluations_used = 0
    any_cap_blocked = 0
    result = 0.d0
    abserr = 0.d0
    work = 0.d0
    adaptive_result_ready = .false.
    floor_used = 0
    floor_reason = 0

    if(adaptive_method.eq.1.and.mode_status.ne.1)then
        if(adaptive_fallback.ne.1)then
            log_integral = -huge(1.d0)
            integral_value = 0.d0
            final_abserr = huge(1.d0)
            final_neval = 0
            final_ifail = 2
            integration_status = -1
            call store_quadrature_diagnostics(family,adaptive_method,0,0,2,0,&
                1,0,0,huge(1.d0),huge(1.d0),log_integral)
            return
        else if(posterior_adaptive_active.eq.1)then
            ! The mode tolerance was not met, but setup retained a finite
            ! centre and positive-definite curvature.  Using that proposal is
            ! still an exact change of variables and is preferable to the
            ! prior-centred rule.
            fallback_used = 1
            method_used = 5
        else
            posterior_adaptive_active = 0
            posterior_adaptive_fallback_active = 1
            fallback_used = 1
            method_used = 2
        end if
    end if

    ! Above dimension four the fully symmetric rules have signed weights.
    ! Run adaptive quadrature in two stages.  The first stage preserves a
    ! positive completed rule before the final budget is spent; if the later
    ! rule becomes nonpositive we can return the earlier valid approximation
    ! instead of leaving only a few evaluations for prior-centred fallback.
    if(adaptive_method.eq.1.and.posterior_adaptive_active.eq.1.and.n.gt.4)then
        stage_budget = max(cubature_minpts,int(0.70d0*dble(cubature_maxpts)))
        stage_budget = min(stage_budget,cubature_maxpts)
        restar = 0
        call_minpts = min(cubature_minpts,stage_budget)
        call hrmsym(n,nf,call_minpts,stage_budget,vraistot,epsabs,&
            epsrel,restar,result,abserr,call_neval,call_ifail,work)
        cubature_calls = cubature_calls+1
        if(call_ifail.ne.0)cubature_failures = cubature_failures+1
        evaluations_used = call_neval
        if(cubature_cap_blocked.eq.1)any_cap_blocked = 1

        if(cubature_signed_sum_status.eq.0.and.&
            cubature_result_is_finite_positive(result(1),abserr(1)))then
            stage_result = result(1)
            stage_error = abserr(1)
            stage_rule = cubature_last_rule
            stage_signed_status = cubature_signed_sum_status
            stage_log_positive = cubature_log_positive
            stage_log_negative = cubature_log_negative
            stage_cancellation_ratio = cubature_cancellation_ratio
            remaining_pts = max(0,cubature_maxpts-evaluations_used)
            if(remaining_pts.gt.0)then
                restar = 1
                call_minpts = min(cubature_minpts,remaining_pts)
                call hrmsym(n,nf,call_minpts,remaining_pts,vraistot,epsabs,&
                    epsrel,restar,result,abserr,call_neval,call_ifail,work)
                cubature_calls = cubature_calls+1
                if(call_ifail.ne.0)cubature_failures = cubature_failures+1
                evaluations_used = evaluations_used+call_neval
                if(cubature_cap_blocked.eq.1)any_cap_blocked = 1
            end if
            if(cubature_signed_sum_status.ne.0.or.&
                .not.cubature_result_is_finite_positive(result(1),abserr(1)))then
                result(1) = stage_result
                abserr(1) = stage_error
                cubature_last_rule = stage_rule
                cubature_signed_sum_status = stage_signed_status
                cubature_log_positive = stage_log_positive
                cubature_log_negative = stage_log_negative
                cubature_cancellation_ratio = stage_cancellation_ratio
                call_ifail = 1
                fallback_used = 1
                method_used = 3
            end if
            final_neval = evaluations_used
            final_ifail = call_ifail
            adaptive_result_ready = .true.
        end if

        if(.not.adaptive_result_ready)then
            ! The initial Hessian proposal was already unstable.  Spend the
            ! reserved budget on a broader exact Gaussian proposal before trying
            ! the prior-centred fallback.
            remaining_pts = max(0,cubature_maxpts-evaluations_used)
            if(adaptive_fallback.eq.1.and.remaining_pts.gt.0)then
                call set_posterior_proposal_scale(n,1.5d0)
                fallback_used = 1
                method_used = 4
                result = 0.d0
                abserr = 0.d0
                work = 0.d0
                restar = 0
                call_minpts = min(cubature_minpts,remaining_pts)
                call hrmsym(n,nf,call_minpts,remaining_pts,vraistot,epsabs,&
                    epsrel,restar,result,abserr,call_neval,call_ifail,work)
                cubature_calls = cubature_calls+1
                if(call_ifail.ne.0)cubature_failures = cubature_failures+1
                evaluations_used = evaluations_used+call_neval
                if(cubature_cap_blocked.eq.1)any_cap_blocked = 1
                if(cubature_signed_sum_status.eq.0.and.&
                    cubature_result_is_finite_positive(result(1),abserr(1)))then
                    final_neval = evaluations_used
                    final_ifail = call_ifail
                    adaptive_result_ready = .true.
                end if
            end if

            if(.not.adaptive_result_ready)then
                ! Last resort: log-stabilized prior-centred quadrature using only the
                ! still-unspent part of the same strict family budget.
                posterior_adaptive_active = 0
                posterior_callback_returns_log = 0
                posterior_adaptive_fallback_active = 1
                fallback_used = 1
                method_used = 2
                result = 0.d0
                abserr = 0.d0
                work = 0.d0
            end if
        end if
    end if

100 continue
    if(.not.adaptive_result_ready)then
        restar = 0
        remaining_pts = max(0,cubature_maxpts-evaluations_used)
        if(remaining_pts.eq.0)then
            call_neval = 0
            call_ifail = 1
            abserr = huge(1.d0)
            cubature_cap_blocked = 1
        else
            call_minpts = min(cubature_minpts,remaining_pts)
            call hrmsym(n,nf,call_minpts,remaining_pts,vraistot,epsabs,&
                epsrel,restar,result,abserr,call_neval,call_ifail,work)
        end if
        evaluations_used = evaluations_used+call_neval
        final_neval = evaluations_used
        final_ifail = call_ifail
        if(cubature_cap_blocked.eq.1)any_cap_blocked = 1

        cubature_calls = cubature_calls+1
        if(final_ifail.ne.0)cubature_failures = cubature_failures+1

        if(posterior_adaptive_active.eq.1)then
            if(result(1).le.0.d0.or.result(1).ne.result(1).or.&
                cubature_signed_sum_status.ne.0)then
                if(adaptive_fallback.eq.1)then
                    posterior_adaptive_active = 0
                    posterior_callback_returns_log = 0
                    posterior_adaptive_fallback_active = 1
                    fallback_used = 1
                    method_used = 2
                    goto 100
                else
                    log_integral = -huge(1.d0)
                    integral_value = 0.d0
                    final_abserr = huge(1.d0)
                    final_ifail = 3
                    integration_status = -2
                    relative_error = huge(1.d0)
                    if(adaptive_nonpositive_action.eq.1.and.&
                        cubature_signed_sum_status.eq.-1.and.&
                        abserr(1).eq.abserr(1).and.abserr(1).ge.0.d0.and.&
                        abs(abserr(1)).lt.huge(1.d0))then
                        floor_reason = 1
                        goto 150
                    end if
                    goto 200
                end if
            end if
        end if
    end if

    if(posterior_adaptive_active.eq.1)then
        log_scale = posterior_log_mode+posterior_log_jacobian
        log_integral = dlog(result(1))+log_scale
        relative_error = abserr(1)/max(abs(result(1)),1.d-300)
        if(log_scale+dlog(max(abserr(1),tiny(1.d0))).gt.&
            dlog(huge(1.d0)))then
            final_abserr = huge(1.d0)
        else
            final_abserr = abserr(1)*dexp(log_scale)
        end if
    else
        if(cubature_signed_sum_status.ne.0.or.&
            .not.cubature_result_is_finite_positive(result(1),abserr(1)))then
            log_integral = -huge(1.d0)
            integral_value = 0.d0
            final_abserr = huge(1.d0)
            final_ifail = 4
            integration_status = -3
            relative_error = huge(1.d0)
            cubature_failures = cubature_failures+1
            if(adaptive_nonpositive_action.eq.1)then
                if(cubature_signed_sum_status.eq.-1.and.&
                    abserr(1).eq.abserr(1).and.abserr(1).ge.0.d0.and.&
                    abs(abserr(1)).lt.huge(1.d0))then
                    floor_reason = 1
                    goto 150
                else if(cubature_signed_sum_status.eq.0.and.&
                    result(1).eq.result(1).and.&
                    abs(result(1)).lt.huge(1.d0).and.result(1).le.0.d0.and.&
                    abserr(1).eq.abserr(1).and.abserr(1).ge.0.d0.and.&
                    abs(abserr(1)).lt.huge(1.d0))then
                    floor_reason = 2
                    goto 150
                end if
            end if
            goto 200
        end if
        ! Both direct and fallback paths preserve every finite positive value
        ! in strict mode.  The optional floor is applied on the original
        ! family-integral scale below, never silently at 1e-150.
        log_integral = dlog(result(1))
        final_abserr = abserr(1)
        relative_error = abserr(1)/max(abs(result(1)),1.d-300)
    end if

    if(adaptive_nonpositive_action.eq.1.and.&
        log_integral.lt.integration_log_floor)then
        if(posterior_adaptive_active.eq.1)then
            floor_reason = 4
        else
            floor_reason = 3
        end if
        goto 150
    end if

    if(log_integral.lt.dlog(tiny(1.d0)))then
        integral_value = tiny(1.d0)
    else if(log_integral.gt.dlog(huge(1.d0)))then
        integral_value = huge(1.d0)
    else
        integral_value = dexp(log_integral)
    end if
    goto 200

150 continue
    log_integral = integration_log_floor
    integral_value = integration_floor_value
    integration_status = 0
    floor_used = 1
    call record_integration_floor(family,floor_reason)

200 continue
    ! Restore the locally returned mode status after quadrature.  No valid
    ! native path returns a positive status other than one.
    if(adaptive_method.eq.1)then
        select case(mode_status)
        case(1,-1,-2,-3,-4,-5)
            posterior_mode_status = mode_status
        case default
            posterior_mode_status = -9
        end select
        call store_mode_diagnostics(family,n)
    end if
    quadrature_failed = 0
    if(final_ifail.ne.0.or.integration_status.ne.0)quadrature_failed = 1
    reached_cap = 0
    if(final_neval.ge.cubature_maxpts.or.any_cap_blocked.eq.1)reached_cap = 1
    cubature_max_abserr = max(cubature_max_abserr,final_abserr)
    cubature_max_relerr = max(cubature_max_relerr,relative_error)
    cubature_max_neval = max(cubature_max_neval,final_neval)
    cubature_max_dimension = max(cubature_max_dimension,n)
    call store_quadrature_diagnostics(family,method_used,cubature_last_rule,&
        final_neval,final_ifail,fallback_used,quadrature_failed,reached_cap,&
        floor_used,&
        final_abserr,relative_error,log_integral)
    posterior_adaptive_active = 0
    posterior_callback_returns_log = 0
    posterior_adaptive_fallback_active = 0
end subroutine integrate_correlated_family


!==========================  DISTANCE   =================================
! fonction supprimée car déjà définie dans funcpassplines.f90


  subroutine gausshermiteBIS2017(ss,npg)

        double precision::ss,ss1
       double precision::auxfunca
      double precision,external::func301
        INTEGER::ii,jj,npg
        double precision,dimension(npg)::X,W

    

        x(1)=-6.863345294d0
        x(2)=-6.13827922d0
        x(3)=-5.533147152d0
        x(4)=-4.988918969d0
        x(5)=-4.483055357d0
        x(6)=-4.003908604d0
        x(7)=-3.544443873d0
        x(8)=-3.09997053d0
        x(9)=-2.667132125d0
        x(10)=-2.243391468d0
        x(11)=-1.826741144d0
        x(12)=-1.4155278d0
        x(13)=-1.008338271d0
        x(14)=-0.603921059d0
        x(15)=-0.201128577d0
        x(16)=0.201128577d0
        x(17)=0.603921059d0
        x(18)=1.008338271d0
        x(19)=1.4155278d0
        x(20)=1.826741144d0
        x(21)=2.243391468d0
        x(22)=2.667132125d0
        x(23)=3.09997053d0
        x(24)=3.544443873d0
        x(25)=4.003908604d0
        x(26)=4.483055357d0
        x(27)=4.988918969d0
        x(28)=5.533147152d0
        x(29)=6.13827922d0
        x(30)=6.863345294d0

        w(1)=0.834247471d0
        w(2)=0.649097981d0
        w(3)=0.569402692d0
        w(4)=0.522525689d0
        w(5)=0.491057996d0
        w(6)=0.468374813d0
        w(7)=0.451321036d0
        w(8)=0.438177023d0
        w(9)=0.427918063d0
        w(10)=0.419895004d0
        w(11)=0.413679364d0
        w(12)=0.408981575d0
        w(13)=0.405605123d0
        w(14)=0.403419817d0
        w(15)=0.402346067d0
        w(16)=0.402346067d0
        w(17)=0.403419817d0
        w(18)=0.405605123d0
        w(19)=0.408981575d0
        w(20)=0.413679364d0
        w(21)=0.419895004d0
        w(22)=0.427918063d0
        w(23)=0.438177023d0
        w(24)=0.451321036d0
        w(25)=0.468374813d0
        w(26)=0.491057996d0
        w(27)=0.522525689d0
        w(28)=0.569402692d0
        w(29)=0.649097981d0
        w(30)=0.834247471d0

    
        auxfunca=0.d0
    ss=0.d0
    ss1=0.d0
!      call gaussher(gauss,npg)
    ii=0
!write(*,*)"okkkkk"
!stop

    do ii=1,npg
        jj=0
        do jj=1,npg
            auxfunca=func301(x(ii),x(jj))
            ss1 = ss1+w(jj)*(auxfunca)
        end do
        ss = ss+w(ii)*ss1
!   write(*,*)'',ss!func30(x(ii),x(jj))
!    pause
    end do

    end subroutine gausshermiteBIS2017
    
    
        double precision  function func301(frail1,frail2)
! calcul de l integrant, pour un effet aleatoire donne frail et un groupe donne auxig (cf funcpa)
 use optim
    use taillesmultiv
    use comonmultiv,only:nigmeta,alpha1,alpha2,res1meta,res3meta,&
    nig,auxig,alpha,theta,eta,aux1,res1,res3,cdc 
    use comon,only:c,res5,res3_ind,ut
    use donnees_indiv,only:it_rec 
    implicit none

    double precision::frail1,frail2
            double precision :: yscalar,eps,finddet,det,res22,alnorm,prod_cag
            integer :: j,i,jj,k,ier,ii
            double precision,dimension(2*(2+1)/2)::matv
            double precision,dimension(2,1)::  Xea2
            double precision,dimension(2):: uii, Xea22,Xea
            double precision,dimension(1)::uiiui
            double precision,dimension(2,2)::mat
    double precision,dimension(2,2)::Utt
            double precision,parameter::pi=3.141592653589793d0

    

    Xea(1) = frail1 
    Xea(2) = frail2 
    
            Xea2(1:2,1) = Xea(1:2)
    
            Xea22(1:2) = Xea(1:2)
    Utt = 0.d0
    
    do j=1,2
        do k=1,2
            Utt(j,k)=Ut(k,j)
        end do
    end do
            mat = matmul(ut,utt)
    
            jj=0
    do j=1,2
        do k=j,2
        jj=j+k*(k-1)/2
        matv(jj)=mat(j,k)
            end do
            end do
            ier = 0
            eps = 1.d-10
    
            call dsinvj(matv,2,eps,ier)
    
            mat=0.d0
        do j=1,2
                    do k=1,2
                                        if (k.ge.j) then
                mat(j,k)=matv(j+k*(k-1)/2)
                else
                mat(j,k)=matv(k+j*(j-1)/2)
                end if
            end do
                        end do
    
                    uii = matmul(Xea22,mat)
                    det = finddet(matmul(ut,utt),2)
    
    
                    uiiui=matmul(uii,Xea2)
    
    
    func301=0.d0
 !   func30 = frail1*(cdc(auxig)*alpha1+nig(auxig)) &
 !   +frail2*(cdc(auxig)*alpha2+nigmeta(auxig)) &
 !   -dexp(frail1)*(res1(auxig)-res3(auxig))-dexp(frail2)*(res1meta(auxig)-res3meta(auxig)) &
 !   -dexp(frail1*alpha1+frail2*alpha2)*aux1(auxig) &
 !   +(2.d0*((2.d0*dexp(alpha)/(dexp(alpha)+1.d0))-1.d0) &
 !   *frail1*frail2/sqrt(theta*eta)  &
 !   -(frail1**2.d0)/theta -(frail2**2.d0)/eta) &
 !   /(2.d0*(1.d0-((2.d0*dexp(alpha)/(dexp(alpha)+1.d0))-1.d0)**2.d0))
 !   func30 = dexp(func30)
 
 do j=1,2
        
    func301 = func301+dlog(dexp((xea(j))*c(it_rec+j-1)-dexp(xea(j))*(res5(it_rec+j-1)-&
            res3_ind(it_rec+j-1))))-uiiui(1)/2.d0-0.5d0*dlog(det)&
                       -dlog(2.d0*pi)!&

!    vrais= dexp(ui(1)*nig(auxig) &
  !  -dexp(ui(1))*(res1(auxig)-res3(auxig)))
    

        end do
        
 func301 = dexp(func301)
    return

    end function func301
