    module tranpsort_module

        use Set_Precision, only: Prec



        implicit none


        contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This Part is used to calculate bed roughness                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine bedrough(mbc,mx,my,u,v,h)

                use sediment_module, only: pbbed,hcr,D,g,m0,rhos,rho
                use Set_Precision, only: Prec

                implicit none

                ! argument
                integer, intent(in) :: mbc,mx,my
                real(kind=Prec), intent(in) :: h(1-mbc:mx+mbc,1-mbc:my+mbc),u(1-mbc:mx+mbc,1-mbc:my+mbc),&
                        v(1-mbc:mx+mbc,1-mbc:my+mbc)

                !local

                integer, :: i,j
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) ::   vmag2,hloc,Dmm,zon,a2,ustarc,ustarcrit, &
                                taub,taucrit,Tstar,delb,zos
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: frc
                Real(kind=Prec) :: a1,gammaWs

                !output
                Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: z0

                gammaWs = 0.056
                frc = pbbed(:,:,1,:)
                vmag2 = u**2.0+v**2.0
                hloc = max(h,hcr)
                Dmm = meansize(D,frc,mx,my,mbc,gmax)
                zon = Dmm/30.0
                a1 = 0.68
                a2 = 0.0204*(log(Dmm*100))**2.0+0.0220*log(Dmm*100)+0.0709
                do i = 1-mbc, mx+mbc
                    do j = 1-mbc, my+mbc
                        ustarc(i,j) = sqrt(g*m0**2.0*vmag2(i,j)/((rhos-rho)*Dmm(i,j)*hloc(i,j)**(1.0/3.0)))
                        ustarcrit(i,j) = sqrt(g*m0**2.0*ub_cr(i,j,int(gmax/2.0+1.0))**2.0/((rhos-rho)*Dmm(i,j)&
                            *hloc(i,j)**(1.0/3.0)))
                    enddo
                enddo
                taub = rho*ustarc**2.0
                taucrit=rho*ustarcrit**2.0
                do i = 1-mbc,mx+mbc
                    do j =1-mbc, my+mbc
                        Tstar(i,j)=taub(i,j)/taucrit(i,j)
                        delb(i,j)=Dmm(i,j)*a1*Tstar(i,j)/(1.0+a2(i,j)*Tstar(i,j))
                    enddo
                enddo
                zos = gammaWs*delb
                z0 = zon+zos
            end subroutine bedrough
            !*******************************************************************
            !
            !
            !This function is used to caculate mean grain size
            !
            !
            !********************************************************************
            function meansize(D,fr,mx,my,mbc,gmax)



                implicit none

                Real(kind=Prec) :: mx,my,mbc,gmax
                Real(kind=Prec), Dimension(gmax) :: D
                Real(kind=Prec), Dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: fr
                Real(kind=Prec), Dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: meansize
                integer, :: i,j,k
                meansize = 0.0

                do i = 1-mbc, mx+mbc
                    do j = 1-mbc, my+mbc
                        do k = 1, gmax
                            meansize(i,j) = meansize(i,j)+D(k)*fr(i,j,k)
                        end do
                    end do
                end do
            end function meansize
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! This Part is used to calculate critical velocity for each grain size classes!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine crtical_velocity1(mbc,mx,my,h)

                use sediment_module, only: rhos,rho,gmax,g,D,hcr,k0,m0
                use Set_Precision, only: Prec

                implicit none

                ! argument
                integer, intent(in) :: mbc,mx,my
                real(kind=Prec), intent(in) :: h(1-mbc:mx+mbc,1-mbc:my+mbc)

                !local
                integer, :: i,j,k
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: hloc
                Real(kind=Prec) :: delta
                !output
                Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: ub_cr, us_cr1, us_cr2

                call settling_velocity(mbc,mx,my)

                delta  = (rhos-rho)/rho

                hloc = max(h,hcr)

                do i = 1-mbc, mx+mbc
                    do j = 1-mbc, my+mbc
                        do k = 1, gmax
                            if (D(k) <= 0.0005) then
                                ub_cr(i,j,k) = 0.19*D(k)**0.1*dlog10(4*hloc(i,j)/D(k))
                            elseif (D(k) > 0.0005 .AND. D(k) <= 0.002) then
                                ub_cr(i,j,k) = 8.5*D(k)**0.6*dlog10(4*hloc(i,j)/D(k))
                            else
                                print *, "The Grain size is out of range"
                            end if
                            us_cr1(i,j,k) = ws(i,j,k)*sqrt(delta)*hloc(i,j)**(1.0/6.0)/(2.5*k0*sqrt(g)*m0)
                            us_cr2(i,j,k) = ws(i,j,k)*sqrt(delta)*hloc(i,j)**(1.0/6.0)/(1.2*k0*sqrt(g)*m0)
                        enddo
                    enddo
                enddo
            end subroutine crtical_velocity1

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! This Part is used to calculate critical velocity for each grain size classes!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine crtical_velocity2(mbc,mx,my,h)

                use sediment_module, only: rhos,rho,gmax,g,D,hcr,Trep,k0,m0
                use Set_Precision, only: Prec

                implicit none

                ! argument
                integer, intent(in) :: mbc,mx,my
                real(kind=Prec), intent(in) :: h(1-mbc:mx+mbc,1-mbc:my+mbc)

                !local
                integer, :: i,j,k
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) :: hloc
                Real(kind=Prec) :: delta
                !output
                Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: ub_cr, us_cr1, us_cr2

                call settling_velocity(mbc,mx,my)

                delta  = (rhos-rho)/rho

                hloc = max(h,hcr)

                do i = 1-mbc, mx+mbc
                    do j = 1-mbc, my+mbc
                        do k = 1, gmax

                            if (D(k) <= 0.0005) then
                                ub_cr(i,j,k) = 0.24*(delta*g)**(2.0/3.0)*(D(k)*Trep)**(1.0/3.0)
                            else
                                ub_cr(i,j,k) = 0.95*(delta*g)**(0.57)*D(k)**0.43*Trep**0.14
                            end if
                            us_cr1(i,j,k) = ws(i,j,k)*sqrt(delta)*hloc(i,j)**(1.0/6.0)/(2.5*k0*sqrt(g)*m0)
                            us_cr2(i,j,k) = ws(i,j,k)*sqrt(delta)*hloc(i,j)**(1.0/6.0)/(1.2*k0*sqrt(g)*m0)
                        end do
                    end do
                end do
            end subroutine crtical_velocity2

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! This Part is used to calculate settling velocity for each grain size classes!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            subroutine settling_velocity(mbc,mx,my)

                use sediment_module, only: rhos,rho,gmax,g,D
                use Set_Precision, only: Prec

                !use params

                implicit none

                ! argument
                integer, intent(in) :: mbc,mx,my

                !local
                integer, :: i,j,k
                Real(kind=Prec) :: delta,Te,vis,Sster,c1,c2,wster
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: C
                Real(kind=Prec),dimension(gmax) :: w,R,alpha
                !output
                Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: ws

                delta  = (rhos-rho)/rho

                C = ccbg+cc

                do k = 1, gmax
                    Te    = 20.0
                    vis  = 4.0/(20.0+Te)*1e-5 ! Van rijn, 1993
                    Sster = D(k)/(4*vis)*sqrt((rhos/rho-1)*g*D(k))
                    c1    = 1.06*tanh(0.064*Sster*exp(-7.5/Sster**2.0))
                    c2    = 0.220*tanh(2.34*Sster**(-1.180)*exp(-0.0064*Sster**2.0))
                    wster = c1+c2*Sster
                    w(k) = wster*sqrt((rhos/rho-1.0)*g*D(k))
                    R(k) = w(k)*D(k)/vis
                    alpha(k) = 2.35*(2.0+0.175*R(k)**(3.0/4.0))/(1.0+0.175*R(k)**(3.0/4.0))
                    do i =1-mbc, mx+mbc
                        do j = 1-mbc, my+mbc
                            ws(i,j,k) = (1-C(i,j,k))**alpha(k)*w(k)
                        enddo
                    enddo
                enddo

            end subroutine settling_velocity

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Soulsby-VanRijn Method                                                             !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine sb_vr(mbc,mx,my,u,v,h)

                use sediment_module, only: eps,rhos,rho,gmax,g,D,hcr,tsfac,Tsmin,cmax,sws,pbbed
                use Set_Precision, only: Prec

                implicit none

                ! Arguments
                integer, intent(in) :: mbc,mx,my
                real(kind=Prec), intent(in) ::u(1-mbc:mx+mbc,1-mbc:my+mbc),v(1-mbc:mx+mbc,1-mbc:my+mbc),h(1-mbc:mx+mbc,1-mbc:my+mbc)

                !local
                integer, :: i,j,k
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) ::   wet,vmg,urms,urms2,hloc,Cd,Asb,term1,term2
                Real(kind=Prec),dimension(gmax) :: dster
                Real(kind=Prec) :: delta,Ass,perc
                Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: Ts,ceq,ceqs,ceqb

                !output
                Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: Tsg,ceqbg,ceqsg

                wet = 0.0

                do i = 1-mbc, mx+mbc
                    do j = 1-mb, my+mbc
                        if (h(i,j)>eps) then
                            wet(i,j) = 1.0
                        endif
                    enddo
                enddo

                delta = rhos - rho

                do k=1,gmax
                    dster(k)=(delta*g/1e-12)**(1.0/3.0)*D(k)
                enddo

                urms = 0.0

                vmg  = dsqrt(u**2+v**2)

                urms2  = urms**2.0

                hloc = max(h,hcr)

                ! calculate threshold velocity Ucr for bedload

                call settling_velocity(mbc,mx,my)

                call crtical_velocity1(mbc,mx,my,h)

                !print *,us_cr2
                ! bed roughness

                call bedrough(mbc,mx,my,u,v,h)

                do k = 1, gmax

                    Ts(:,:,k) = tsfac*hloc/ws(:,:,k)
                    Tsg(:,:,k) = max(Ts(:,:,k),Tsmin)

                    ! drag coefficient
                    do i = 1-mbc, mx+mbc
                        do j = 1-mbc, my+mbc
                            Cd(i,j)=(0.40/(log(max(hloc(i,j),10.0*z0(i,j))/z0(i,j))-1.0))**2.0
                        end do
                    end do

                    ! transport parameters
                    Asb=0.005*hloc*(D(k)/hloc/(delta*g*D(k)))**1.20         ! bed load coefficent
                    Ass=0.0120*D(k)*dster(k)**(-0.60)/(delta*g*D(k))**1.20  ! suspended load coeffient

                    term1 = (vmg**2.0+0.018/Cd*sws*urms2)
                    !term1 = min(term1,smax*g/cf*D(k)*delta)
                    term1 = sqrt(term1)
                    term2 = 0.0*term1
                    !ceqb = 0.0 !initialize ceqb
                    !ceqs = 0.0 !initialize ceqs

                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc
                            if(term1(i,j)>Ub_cr(i,j,k) .and. hloc(i,j)>eps) then
                                term2(i,j)=(term1(i,j)-Ub_cr(i,j,k))**2.40
                            end if
                            ceq(i,j,k) = (Asb(i,j)+Ass)*term2(i,j)/hloc(i,j)
                        enddo
                    enddo
                    do j = 1-mbc,my+mbc
                        do i = 1-mbc,mx+mbc
                            if(term1(i,j)<Us_cr1(i,j,k)) then
                                ceqb(i,j,k) = min(ceq(i,j,k),cmax/gmax/2.0)
                                ceqs(i,j,k) = 0.0
                            elseif(term1(i,j)>Us_cr2(i,j,k)) then
                                ceqb(i,j,k) = 0.0
                                ceqs(i,j,k) = min(ceq(i,j,k),cmax/gmax/2.0)
                            else
                                perc = term1(i,j)/Us_cr2(i,j,k)
                                ceqb(i,j,k) = (1-perc)*min(ceq(i,j,k),cmax/gmax/2.0)
                                ceqs(i,j,k) = perc*min(ceq(i,j,k),cmax/gmax/2.0)
                            end if
                            ceqbg(i,j,k) = ceqb(i,j,k)*wet(i,j)
                            ceqsg(i,j,k) = ceqs(i,j,k)*wet(i,j)
                        enddo
                    enddo
                enddo
            end subroutine sb_vr

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Van Thiel-Van Rijn Method                                                          !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            subroutine vt_vr(mbc,mx,my,u,v,h)

                use sediment_module, only: eps,rhos,rho,gmax,g,D,hcr,tsfac,Tsmin,cmax,sws,pbbed
                use Set_Precision, only: Prec

                implicit none

                ! Arguments
                integer, intent(in) :: mbc,mx,my
                real(kind=Prec), intent(in) ::u(1-mbc:mx+mbc,1-mbc:my+mbc),v(1-mbc:mx+mbc,1-mbc:my+mbc),h(1-mbc:mx+mbc,1-mbc:my+mbc)

                !local

                integer, :: i,j,k
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) ::   wet,vmg,urms,urms2,hloc,Cd,Asb,term1,term2,term3
                Real(kind=Prec),dimension(gmax) :: dster
                Real(kind=Prec) :: delta,Ass,perc
                Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: Ts,ceq,ceqs,ceqb
                !output
                Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: Tsg,ceqbg,ceqsg


                wet = 0.0

                do i = 1-mbc, mx+mbc
                    do j = 1-mb, my+mbc
                        if (h(i,j)>eps) then
                            wet(i,j) = 1.0
                        endif
                    enddo
                enddo

                delta = rhos - rho

                do k=1,gmax
                    dster(k)=(delta*g/1e-12)**(1.0/3.0)*D(k)
                enddo

                urms = 0.0

                vmg  = dsqrt(u**2+v**2)

                urms2  = urms**2.0

                hloc = max(h,hcr)

                ! calculate threshold velocity Ucr for bedload

                call settling_velocity(mbc,mx,my)

                call crtical_velocity2(mbc,mx,my,h)

                do k = 1,gmax

                    ! transport parameters

                    Asb=0.015*hloc*(D(k)/hloc)**1.20/(delta*g*D(k))**0.75        !bed load coefficent
                    Ass=0.012*D(k)*dster(k)**(-0.60)/(delta*g*D(k))**1.20        !suspended load coeffient

                    term1=vmg**2+0.640*sws*urms2

                    !term1=min(term1,smax*g/cf*D(k)*delta)
                    term1=sqrt(term1)

                    !ceqb = 0.0*term1                                                                     !initialize ceqb
                    !ceqs = 0.0*term1 
                    !initialize ceqs
                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc
                            if(term1(i,j)>Ub_cr(i,j,k) .and. h(i,j)>eps) then
                                term2(i,j)=(term1(i,j)-Ub_cr(i,j,k))**1.50
                                term3(i,j)=(term1(i,j)-Ub_cr(i,j,k))**2.40
                            endif
                        enddo
                    enddo
                    ceq(:,:,k) = (Asb*term2+Ass*term3)/hloc
                    do j = 1-mbc,my+mbc
                        do i = 1-mbc,mx+mbc
                            if(term1(i,j)<Us_cr1(i,j,k)) then
                                ceqb(i,j,k) = min(ceq(i,j,k),cmax/gmax/2.0)
                                ceqs(i,j,k) = 0.0
                            elseif(term1(i,j)>Us_cr2(i,j,k)) then
                                ceqb(i,j,k) = 0.0
                                ceqs(i,j,k) = min(ceq(i,j,k),cmax/gmax/2.0)
                            else
                                !perc = term1(i,j)/Us_cr2(i,j,k)
                                ceqb(i,j,k) = Asb(i,j)*term2(i,j)/hloc(i,j)!(1-perc)*min(ceq(i,j,k),cmax/gmax/2.0)
                                ceqs(i,j,k) = Ass*term3(i,j)/hloc(i,j)!perc*min(ceq(i,j,k),cmax/gmax/2.0)
                            endif
                            ceqbg(i,j,k) = ceqb(i,j,k)*wet(i,j)
                            ceqsg(i,j,k) = ceqs(i,j,k)*wet(i,j)
                        enddo
                    enddo
                enddo
            end subroutine vt_vr
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! This Part is used to calculate Source term for finite volume method                !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine transus(mbc,mx,my,dx,dy,time,u,v,h,dt)

                use flux, only: Flux_vector
                use sediment_module, only: trim,gmax,morfac,por,D,thetanum,cmax,lmax,pbbed,dzbed
                use Set_Precision, only: Prec

                implicit none

                ! Arguments
                integer, intent(in) :: mbc,mx,my
                real(kind=Prec), intent(in) :: dx,dy,dt,time
                real(kind=Prec), intent(in) ::u(1-mbc:mx+mbc,1-mbc:my+mbc),v(1-mbc:mx+mbc,1-mbc:my+mbc),h(1-mbc:mx+mbc,1-mbc:my+mbc)

                !local

                integer, :: i,j,k
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc) ::   vmag2,dcsdy,dcbdy,dcsdx,dcbdx,hold !hold?
                Real(kind=Prec),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: frc,fac,ero1,ero2,depo_ex1,depo_ex2, &
                                cc,ccb
                Real(kind=Prec) :: exp_ero

                !out
                Real(kind=Prec),intent(inout),dimension(1-mbc:mx+mbc,1-mbc:my+mbc,gmax) :: cu,cub,cv,cvb,ccg,ccbg,susg,Subg,Svbg,Svsg

                vmag2     = u**2+v**2
                ! calculate equibrium sediment concentration
                if (trim=='soulsby_vanrijn') then           ! Soulsby van Rijn
                    call sb_vr(mbc,mx,my,u,v,h)
                elseif (trim=='vanthiel_vanrijn') then       ! Van Thiel de Vries & Reniers 2008
                    call vt_vr(mbc,mx,my,u,v,h)
                end if
                ! compute reduction factor for sediment sources due to presence of hard layers
                frc = pbbed(:,:,1,:)
                do k = 1,gmax
                    do j= 1-mbc,my+mbc
                        do i= 1-mbc,mx+mbc
                            exp_ero = morfac*dt/(1.0-por)*h(i,j)*(ceqsg(i,j,k)*frc(i,j,k)/Tsg(i,j,k) &
                                    + ceqbg(i,j,k)*frc(i,j,k)/dt) !ceqsg, ceqbg from vt_vr or sb_vr
                            fac(i,j,k) =min(1.0,dzbed(i,j,1)*frc(i,j,k)/max(tiny(0.0),exp_ero))! what's laythick? TODO
                            !print *, exp_ero

                        enddo
                    enddo
                enddo
                ! compute diffusion coefficient
                cc = ccg
                ccb = ccbg
                do k = 1,gmax

                    if (D(k)>0.002) then
                        print *, "WARNING: Grain size is larger than 2 mm"
                    endif
                    ! x-direction
                    do j=1-mbc,my+mbc
                        do i=1-mbc,mx+mbc-1
                            if(u(i,j)>0.0) then
                                cu(i,j,k)=thetanum*cc(i,j,k)+(1.0-thetanum)*cc(min(i+1,mx+mbc),j,k)
                                cub(i,j,k)=thetanum*ccb(i,j,k)+(1.0-thetanum)*ccb(min(i+1,mx+mbc),j,k)
                            elseif(u(i,j)<0.0) then
                                cu(i,j,k)=thetanum*cc(i+1,j,k)+(1.0-thetanum)*cc(max(i,2-mbc),j,k)
                                cub(i,j,k)=thetanum*ccb(i+1,j,k)+(1.0-thetanum)*ccb(max(i,2-mbc),j,k)
                            else
                                cu(i,j,k)=0.50*(cc(i,j,k)+cc(i+1,j,k))
                                cub(i,j,k)=0.50*(ccb(i,j,k)+ccb(i+1,j,k))
                            endif
                            dcsdx(i,j)=(sum(cc(i+1,j,:))-sum(cc(i,j,:)))/dx
                            dcbdx(i,j)=(sum(ccb(i+1,j,:))-sum(ccb(i,j,:)))/dx
                        enddo
                    enddo
                    cu(mx+mbc,:,:) = cc(mx+mbc,:,:)
                    cub(mx+mbc,:,:) = ccb(mx+mbc,:,:)
                    !y-direction
                    if (jmax>0) then
                        do j=1-mbc,my+mbc
                            do i=1-mbc,mx+mbc
                                if(v(i,j)>0) then
                                    cv(i,j,k)=thetanum*cc(i,j,k)+(1.0-thetanum)*cc(i,min(j+1,my+mbc),k)
                                    cvb(i,j,k)=thetanum*ccb(i,j,k)+(1.0-thetanum)*ccb(i,min(j+1,my+mbc),k)
                                elseif(v(i,j)<0) then
                                    cv(i,j,k)=thetanum*cc(i,j+1,k)+(1.0-thetanum)*cc(i,max(j,2-mbc),k)
                                    cvb(i,j,k)=thetanum*ccb(i,j+1,k)+(1.0-thetanum)*ccb(i,max(j,2-mbc),k)
                                else
                                    cv(i,j,k)=0.50*(cc(i,j,k)+cc(i,j+1,k)) !Jaap: cc instead of cv
                                    cvb(i,j,k)=0.50*(ccb(i,j,k)+ccb(i,j+1,k))
                                end if
                                dcsdy(i,j)=(sum(cc(i,j+1,:))-sum(cc(i,j,:)))/dy !Jaap
                                dcbdy(i,j)=(sum(ccb(i,j+1,:))-sum(ccb(i,j,:)))/dy
                            end do
                        end do
                        cv(:,my+mbc,:) = cc(:,my+mbc,:)
                        cvb(:,my+mbc,:) = ccb(:,my+mbc,:)
                    else
                        cv = cc
                        cvb = ceqbg(:,:,:)
                    endif ! my>0

                    call Flux_vector(mbc,mx,my,u,v,h,cc,ccb)

                    if (my>0) then
                        do j=1-mbc,my+mbc
                            do i=1-mbc,mx+mbc
                                !suspended sediment transport
                                ero1(i,j,k) = fac(i,j,k)*h(i,j)*ceqsg(i,j,k)*pbbed(i,j,1,k)/Tsg(i,j,k)
                                cc(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                        (hold(i,j)*cc(i,j,k)/dt -((Sus(i,j,k)*dx-Sus(i-1,j,k)*dx+&
                                        Svs(i,j,k)*dy-Svs(i,j-1,k)*dy)/(dx*dy)-ero1(i,j,k)))
                                cc(i,j,k)=max(cc(i,j,k),0.00)
                                cc(i,j,k)=min(cc(i,j,k),cmax/2.0*h(i,j))
                                depo_ex1(i,j,k) = cc(i,j,k)/Tsg(i,j,k)
                                !bed sediment tranpsort
                                ero2(i,j,k) = fac(i,j,k)*h(i,j)*ceqbg(i,j,k)*pbbed(i,j,1,k)/Tsg(i,j,k)
                                ccb(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                        (hold(i,j)*cc(i,j,k)/dt -((Sub(i,j,k)*dx-Sub(i-1,j,k)*dx+&
                                        Svb(i,j,k)*dy-Svb(i,j-1,k)*dy)/(dx*dy)-ero2(i,j,k)))
                                ccb(i,j,k)=max(ccb(i,j,k),0.00)
                                ccb(i,j,k)=min(ccb(i,j,k),cmax/2.0*h(i,j))
                                depo_ex2(i,j,k) = ccb(i,j,k)/Tsg(i,j,k)
                            enddo
                        enddo
                    else
                        j=1
                        do i=1-mbc,mx+mbc
                            !suspended sediment transport
                            ero1(i,j,k) = fac(i,j,k)*h(i,j)*ceqsg(i,j,k)*pbbed(i,j,1,k)/Tsg(i,j,k)
                            cc(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                        (hold(i,j)*cc(i,j,k)/dt -((Sus(i,j,k)*dx-Sus(i-1,j,k)*dx)/(dx*dy)-&
                                        ero1(i,j,k)))
                            cc(i,j,k)=max(cc(i,j,k),0.00)
                            cc(i,j,k)=min(cc(i,j,k),cmax*h(i,j))
                            depo_ex1(i,j,k) = cc(i,j,k)/Tsg(i,j,k)
                            !bed sediment tranpsort
                            ero2(i,j,k) = fac(i,j,k)*h(i,j)*ceqbg(i,j,k)*pbbed(i,j,1,k)/Tsg(i,j,k)
                            ccb(i,j,k) = (dt*Tsg(i,j,k))/(dt+Tsg(i,j,k))* &
                                        (hold(i,j)*ccb(i,j,k)/dt -((Sub(i,j,k)*dx-Sub(i-1,j,k)*dx)/(dx*dy)-&
                                        ero2(i,j,k)))
                            ccb(i,j,k)=max(ccb(i,j,k),0.00)
                            ccb(i,j,k)=min(ccb(i,j,k),cmax*h(i,j))
                            depo_ex2(i,j,k) = ccb(i,j,k)/Tsg(i,j,k)
                        enddo
                    endif
                    cc(:,:,k) = cc(:,:,k)/h
                    ccg = cc
                    ccb(:,:,k) = ccb(:,:,k)/h
                    ccbg = ccb
                    !print *, k
                end do
                Svsg = Svs
                Susg = Sus
                Svbg = Svb
                Subg = Sub
            end subroutine transus
        end module tranpsort_module
