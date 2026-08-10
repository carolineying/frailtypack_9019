

!========================          FUNCPA_SPLINES          ====================
    double precision function funcpassplines_log(b,np,id,thi,jd,thj,k0)
    use tailles
    !use comon,only:nz1,nz2,t0,t1
    use comon,only:m3m3,m2m2,m1m1,mmm,m3m2,m3m1,m3m,m2m1,m2m,m1m, &
    mm3,mm2,mm1,mm,im3,im2,im1,im,date,zi,c,t1,nt0,nt1,nt2,nsujet,nva, &
    ndate,kkapa,nst,stra,ve,pe,effet,ng,g,nig,nmemb,AG,resnonpen,sig2,sig22, &
    indictronq,auxig,res3,res5,res6,sig_prob,res1,res2_ind,res3_ind,k0T,Kmat1,Kmat2, corrRE, Ut,&
    Ut_table,nmemb,cpt2,currentage,proband,agemin
    use residusM
    use lois_normales
    use donnees_indiv,only:it_rec
    use optim

    implicit none

! *** NOUVELLLE DECLARATION F90 :

    integer::nb,n,np,id,jd,i,j,k,vj,cptg,ig,choix,jj
    integer,dimension(ngmax)::cpt
    double precision::thi,thj,dnb,res,vet,h1,int,eps,finddet,slogasc
    double precision,dimension(nst)::peT,somT
    double precision,dimension(-2:npmax,nst)::theT
    double precision,dimension(np)::b,bh
    double precision,dimension(ngmax)::res2
!     double precision,dimension(2)::k0
    double precision,dimension(nst)::k0 !en plus
    double precision,dimension(ndatemax,nst)::dutT
    double precision,dimension(0:ndatemax,nst)::utT
    double precision,dimension(ngmax)::integrale1,integrale2,integrale3,integrale,det_all
    double precision,parameter::pi=3.141592653589793d0
    
    ! for the numerical integral hrmsym
        integer :: restar,nf,ier
        double precision:: epsabs,epsrel
        double precision,dimension(2):: result, abserr2
        double precision,dimension(1000) :: work
        external :: vraistot
        integer ::neval,ifail
        

    kkapa=k0
    j=0
    sig2=0.d0
    do i=1,np
        bh(i)=b(i)
    end do

    if (id.ne.0) bh(id)=bh(id)+thi
    if (jd.ne.0) bh(jd)=bh(jd)+thj

    n = (np-nva-effet)/nst

    do jj=1,nst !en plus strates A.Lafourcade 05/2014
        do i=1,n
            theT(i-3,jj)=(bh((jj-1)*n+i))*(bh((jj-1)*n+i))
        end do
    end do

    if(effet.eq.1) then
        sig2 = bh(np-nva)*bh(np-nva)
    endif

!---------  calcul de ut1(ti) et ut2(ti) ---------------------------
!    attention the(1)  sont en nz=1
!        donc en ti on a the(i)
allocate(cpt2(ngmax))
    vj = 0

    somT=0.d0 !en plus
    do jj=1,nst !en plus strates A.Lafourcade 05/2014
        dutT(1,jj) = (theT(-2,jj)*4.d0/(zi(2)-zi(1)))
        utT(0,jj) = 0.d0
        utT(1,jj) = theT(-2,jj)*dutT(1,jj)*0.25d0*(zi(1)-zi(-2))
    end do

    do i=2,ndate-1
        do k = 2,n-2
            if (((date(i)).ge.(zi(k-1))).and.(date(i).lt.zi(k)))then
                j = k-1
                if ((j.gt.1).and.(j.gt.vj))then
                    do jj=1,nst !en plus strates A.Lafourcade 05/2014
                        somT(jj) = somT(jj)+theT(j-4,jj)
                    end do
                    vj  = j
                endif
            endif
        end do
        do jj=1,nst !en plus strates A.Lafourcade 05/2014
            utT(i,jj) = somT(jj) +(theT(j-3,jj)*im3(i))+(theT(j-2,jj)*im2(i)) &
            +(theT(j-1,jj)*im1(i))+(theT(j,jj)*im(i))
            dutT(i,jj) = (theT(j-3,jj)*mm3(i))+(theT(j-2,jj)*mm2(i)) &
            +(theT(j-1,jj)*mm1(i))+(theT(j,jj)*mm(i))
        end do
    end do

    i = n-2
    h1 = (zi(i)-zi(i-1))
    do jj=1,nst !en plus strates A.Lafourcade 05/2014
        utT(ndate,jj)=somT(jj)+theT(i-4,jj)+theT(i-3,jj)+theT(i-2,jj)+theT(i-1,jj)
        dutT(ndate,jj) = (4.d0*theT(i-1,jj)/h1)
    end do

!-------------------------------------------------------
!--------- calcul de la vraisemblance ------------------
!-------------------------------------------------------

!--- avec ou sans variable explicative  ------cc

    res1 = 0.d0
    res2 = 0.d0
    res3 = 0.d0
    res2_ind = 0.d0
    res3_ind = 0.d0
    res6 = 0.d0
    integrale = 1.d0
    integrale1 = 1.d0
    integrale2 = 1.d0
    integrale3 = 1.d0
    res5 = 0.d0
    sig_prob = 0.d0
    cpt = 0
    cpt2 = 0
!*******************************************
!---- sans effet aleatoire dans le modele
!*******************************************

    if (effet.eq.0) then
    
        do i=1,nsujet
            cpt(g(i))=cpt(g(i))+1
            cpt2(g(i))=cpt2(g(i))+1

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
                
                res2(g(i)) = res2(g(i))+dlog(dutT(nt1(i),stra(i))*vet)!log(a*b)=log(a)+log(b)
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge. 1.d30)) then
                funcpassplines_log=-1.d9
                goto 123
            end if

            !en plus strates A.Lafourcade 05/2014
            res1(g(i)) = res1(g(i)) + utT(nt1(i),stra(i))*vet-utT(nt0(i),stra(i))*vet !en plus
            RisqCumul(i) = utT(nt1(i),stra(i))*vet

            if ((res1(g(i)).ne.res1(g(i))).or.(abs(res1(g(i))).ge. 1.d30)) then
                funcpassplines_log=-1.d9
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
                    funcpassplines_log=-1.d9
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
cpt2(g(i))=cpt2(g(i))+1

            if(nva.gt.0)then
                vet = 0.d0
                do j=1,nva
                    vet =vet + bh(np-nva+j)*dble(ve(i,j))
                end do
                vet = dexp(vet)
            else
                vet=1.d0
            endif

          !en plus strates A.Lafourcade 05/2014
            if(c(i).eq.1)then !en plus
                res2_ind(i) = dlog(dutT(nt1(i),stra(i))*vet)
                res2(g(i)) = res2(g(i))+dlog(dutT(nt1(i),stra(i))*vet)
            endif

            if ((res2(g(i)).ne.res2(g(i))).or.(abs(res2(g(i))).ge.1.d30)) then
                funcpassplines_log=-1.d9
                goto 123
            end if

! modification pour nouvelle vraisemblance / troncature:
            res3_ind(i) = utT(nt0(i),stra(i))*vet
            res3(g(i)) = res3(g(i)) + utT(nt0(i),stra(i))*vet !en plus
            res5(i) = utT(nt1(i),stra(i))*vet
            res1(g(i)) = res1(g(i)) + res5(i) ! pour les residus
            res6(i) = utT(nt2(i),stra(i))*vet
        !    res6(i) = ((currentage(i)-agemin)/etaT(stra(i)))**betaT(stra(i))*vet !cumulative hazard for probands
        
!if(res5(i).ge.1000)write(*,*)i,res5(i) ,bh
!write(*,*)i,utT(nt1(i),stra(i)),utT(nt2(i),stra(i)),nt1(i),nt2(i),agemin
            if ((res3(g(i)).ne.res3(g(i))).or.(abs(res3(g(i))).ge.1.d30)) then
            write(*,*)"here1"
                funcpassplines_log=-1.d9
                goto 123
            end if

            if ((res5(i).ne.res5(i)).or.(abs(res5(i)).ge.1.d30)) then
             write(*,*)"here2"
              funcpassplines_log=-1.d9
                goto 123
            end if
        end do

!**************INTEGRALES ****************************

        if(corrRE.eq.0) then
    !    it_rec = 1
            do ig=1,ng
                auxig = ig
                choix = 1
                call gauherS(int,choix)
                integrale1(ig) = int
                if (AG.eq.1) then
                    choix = 3
                    call gauherS(int,choix)
                    integrale3(ig) = int
                    
    !            auxig = ig
    !            epsabs = 1.d-100
    !            epsrel = 1.d-100
    !            restar = 0
    !            nf = 1
    !            allocate(Ut(1,1))
    !            Ut = dsqrt(sig2)
    !                call  hrmsym(1, nf, 30, 500,vraistot, epsabs, &
    !                epsrel, restar, result, abserr2, neval, ifail, work)
    !        !        write(*,*)result(1),int,it_rec,cpt2(ig),ifail
    !        !        stop
    !    integrale3(ig) =result(1)
    !        deallocate(ut)
    !        
    !        it_rec = it_rec + cpt(ig)
                endif
                if (indictronq.eq.1) then
                    choix = 2
                    call gauherS(int,choix)
                    integrale2(ig) = int
                endif
            end do
        else if(corrRE.eq.1) then
            it_rec = 1
            det_all = 0.d0
            do ig = 1,ng
                auxig = ig
                epsabs = 1.d-100
                epsrel = 1.d-100
                restar = 0
                nf = 1
                allocate(Ut(nmemb(ig),nmemb(ig)),Ut_table(nmemb(ig)*(nmemb(ig)+1)/2))
                Ut = 0.d0
                do j=1,nmemb(ig)
                    do k=1,nmemb(ig)
                        Ut(j,k) =sig2* Kmat1(it_rec+j-1,it_rec+k-1)!sqrt(sig2)* 
                    end do
                end do
            
                det_all(ig) = finddet(Ut,nmemb(ig))
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
                
            
                Ut=0.d0
                do j=1,nmemb(ig)
                  do k=1,j
                     Ut(j,k)=Ut_table(k+j*(j-1)/2)
                  end do
                end do
               
                call  hrmsym(nmemb(ig), nf, 30, 750,vraistot, epsabs, &
                    epsrel, restar, result, abserr2, neval, ifail, work)
         
                if (result(1).le.1.d-150) then
                    result(1)=1.d-150
                end if
                integrale(ig) =result(1)
                it_rec = it_rec + nmemb(ig)
                deallocate(Ut,Ut_table)
            end do
        else if(corrRE.eq.2) then 
        it_rec = 1
            det_all = 0.d0
            do ig = 1,ng
                auxig = ig
                epsabs = 1.d-100
                epsrel = 1.d-100
                restar = 0
                nf = 1
                allocate(Ut(nmemb(ig),nmemb(ig)),Ut_table(nmemb(ig)*(nmemb(ig)+1)/2))
                Ut = 0.d0
                do j=1,nmemb(ig)
                    do k=1,nmemb(ig)
                        Ut(j,k) =sig2*sig22* Kmat1(it_rec+j-1,it_rec+k-1)*Kmat2(it_rec+j-1,it_rec+k-1)
                    end do
                end do
            
                det_all(ig) = finddet(Ut,nmemb(ig))
                    Ut = 0.d0
                        do j=1,nmemb(ig)
                    do k=1,nmemb(ig)
                        Ut(j,k) =sig2*sig22* Kmat1(it_rec+j-1,it_rec+k-1)*Kmat2(it_rec+j-1,it_rec+k-1)
                        
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
                
            
                Ut=0.d0
                do j=1,nmemb(ig)
                  do k=1,j
                     Ut(j,k)=Ut_table(k+j*(j-1)/2)
                  end do
                end do
               
                call  hrmsym(nmemb(ig), nf, 30, 750,vraistot, epsabs, &
                    epsrel, restar, result, abserr2, neval, ifail, work)
         
                if (result(1).le.1.d-150) then
                    result(1)=1.d-150
                end if
                integrale(ig) =result(1)
                it_rec = it_rec + nmemb(ig)
                deallocate(Ut,Ut_table)
            end do
        end if
!************* FIN INTEGRALES ************************
!stop
        res = 0.d0

!     gam2 = gamma(inv)
! k indice les groupes
        do k=1,ng
                !if(theta.gt.(1.d-5)) then
!ccccc ancienne vraisemblance : ANDERSEN-GILL ccccccccccccccccccccccccc
                if(corrRE.eq.0) then
                    if(AG.EQ.1)then
                        res = res+res2(k)-dlog(dsqrt(sig2))-dlog(2.d0*pi)/2.d0+ &!
                        dlog(integrale3(k))
!ccccc nouvelle vraisemblance :ccccccccccccccccccccccccccccccccccccccccccccccc
                    else
                        res = res+res2(k)+dlog(integrale1(k)) &
                        -dlog(integrale2(k))
                    endif
                else 
                    res = res+res2(k)-nmemb(k)*dlog(2.d0*pi)/2.d0+&!-0.5d0*dlog(abs(det_all(k)))+&
                        dlog(integrale(k))!+res2(k)
                end if
                    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
                    write(*,*)res,res2(k),0.5d0*dlog(abs(det_all(k))),nmemb(k)*dlog(2.d0*pi)/2.d0, &
                        integrale(k)
                          funcpassplines_log=-1.d9
                          goto 123
                    end if
                    
                !else
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
                        !  funcpassplines_log=-1.d9
                         ! goto 123
                       !end if
                !endif
            !endif
        end do
    endif !fin boucle effet=0

!--------- calcul de la penalisation -------------------

    peT=0.d0    !en plus strates A.Lafourcade 05/2014
    pe=0.d0
    do jj=1,nst !en plus strates A.Lafourcade 05/2014
        do i=1,n-3

            peT(jj) = peT(jj)+(theT(i-3,jj)*theT(i-3,jj)*m3m3(i))+(theT(i-2,jj) &
            *theT(i-2,jj)*m2m2(i))+(theT(i-1,jj)*theT(i-1,jj)*m1m1(i))+( &
            theT(i,jj)*theT(i,jj)*mmm(i))+(2.d0*theT(i-3,jj)*theT(i-2,jj)* &
            m3m2(i))+(2.d0*theT(i-3,jj)*theT(i-1,jj)*m3m1(i))+(2.d0* &
            theT(i-3,jj)*theT(i,jj)*m3m(i))+(2.d0*theT(i-2,jj)*theT(i-1,jj)* &
            m2m1(i))+(2.d0*theT(i-2,jj)*theT(i,jj)*m2m(i))+(2.d0*theT(i-1,jj) &
            *theT(i,jj)*m1m(i))  
        end do
    end do

    do jj=1,nst !en plus strates A.Lafourcade 05/2014
        pe=pe+pet(jj)*k0T(jj)
    end do

    resnonpen = res

    res = res - pe

    slogasc  = 0.d0
    if(effet.eq.1.and.corrRE.ge.1)then
        do i = 1,nsujet 
            if(proband(i).eq.1) then 
                auxig = i
                call gauherASC(int,sig_prob(i))
                if(t1(i).lt.currentage(i)) then
                    slogasc = slogasc + dlog(1-int)
                else 
                    slogasc = slogasc + dlog(int)
                end if
            end if 
        end do
        res = res - slogasc
    end if
    if ((res.ne.res).or.(abs(res).ge. 1.d30)) then
      funcpassplines_log=-1.d9
      goto 123
    end if

    funcpassplines_log = res

    do k=1,ng
        cumulhaz(k)=res1(k)
    end do
    do k = 1,nsujet 
        cumulhaz_corr(k) = res5(k) 
    end do

123     continue
deallocate(cpt2)
    return

    end function funcpassplines_log


!==========================  DISTANCE   =================================
! fonction supprimée car déjà définie dans funcpassplines.f90



   !****************** for GENZ algorithm ********************
    subroutine vraistot(nea,xea,nf,funvls)
    
    use lois_normales
    use comon,only:auxig,sig2,res1,res3,res3_ind,res2_ind,res5,nig, ut,ut22,c,ag,cpt2,&
                    nmemb,corrRE
    use optim
    use random_effect
    use comongroup,only:vet,vet2
    use donnees_indiv,only:it_rec 
    use posterior_adaptive_settings,only:posterior_adaptive_active,&
        posterior_callback_returns_log,posterior_log_mode,&
        adaptive_node_transform
    
    implicit none
    integer :: nea,nea2 ,j,nf,k,transform_status
    double precision,dimension(nea):: xea!,xea2
    double precision,dimension(nmemb(auxig)):: xea1,xea2
    double precision :: funvls,vrais,det,finddet,log_vrais
    double precision,dimension(:),allocatable::ui,ui2
    double precision,parameter::pi=3.141592653589793d0
    double precision,dimension(nmemb(auxig),nmemb(auxig))::Utt,Utt2,matrix,SUtt

    nea2 = nmemb(auxig)
  
    Xea1=0.d0
    Xea2=0.d0
    if(posterior_adaptive_active.eq.1)then
        call adaptive_node_transform(nea2,Xea,Xea1,transform_status)
        if(transform_status.ne.0)then
            if(posterior_callback_returns_log.eq.1)then
                funvls = -huge(1.d0)
            else
                funvls = 0.d0
            end if
            return
        end if
    else
        do j=1,nea2
            Xea1(j)=Xea(j)
        end do
    end if
    Xea2 = 0.d0
    
    


    allocate(ui(nea2),ui2(nea2))
    ui2 = 0.d0

     Utt = 0.d0
     Utt2 = 0.d0
    
    do j=1,nea2
        do k=1,nea2
            Utt(j,k)=Ut(k,j)!+Ut22(k,j)
    !        SUtt(j,k) = Ut(j,k) + Ut22(j,k)
        end do
    end do
    !det = finddet(MATMUL(SUtt,Utt),nea2)


    ui = MATMUL(Ut,Xea1)
    if(corrRE.eq.2)ui2 = 0!MATMUL(Ut22,Xea2)
    if(corrRE.eq.1)ui2 = 0


    
    
    vrais = 1.d0
    log_vrais = 0.d0
    if(AG.eq.1) then !calendar time
    
        do j=1,nea2

            if(ui(j)+ui2(j).gt.700.d0.and.&
                res5(it_rec+j-1)-res3_ind(it_rec+j-1).gt.0.d0)then
                log_vrais = -huge(1.d0)
                exit
            else
                log_vrais = log_vrais+(ui(j)+ui2(j))*c(it_rec+j-1)-&
                    dexp(ui(j)+ui2(j))*(res5(it_rec+j-1)-res3_ind(it_rec+j-1))
            end if

!    vrais= dexp(ui(1)*nig(auxig) &
  !  -dexp(ui(1))*(res1(auxig)-res3(auxig)))
    

        end do
        if(posterior_adaptive_active.eq.1)then
            log_vrais = log_vrais-0.5d0*dot_product(Xea1,Xea1)-&
                posterior_log_mode+0.5d0*dot_product(Xea,Xea)
        end if
        if(posterior_callback_returns_log.eq.1)then
            funvls = log_vrais
        else if(log_vrais.lt.dlog(tiny(1.d0)))then
            funvls = 0.d0
        else if(log_vrais.gt.dlog(huge(1.d0)))then
            funvls = huge(1.d0)
        else
            funvls = dexp(log_vrais)
        end if
    else !gap time
        if(posterior_callback_returns_log.eq.1)then
            funvls = 0.d0
        else
            funvls = 1.d0
        end if
    end if
!    if(auxig.eq.22.or.auxig.eq.23)write(*,*)vrais,nea2,ui,ui2
!write(*,*)vrais,dsqrt(det)
!    if(auxig.eq.9)write(*,*)ut
    !write(*,*)"xea2",xea2
    !write(*,*)"ut",ut
   deallocate(ui,ui2)

    end subroutine vraistot
   
