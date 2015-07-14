
!*******************************************************************
!
!           TEST
!
!*******************************************************************

    program Test

        use Set_Precision
        use params
        use sed
        use flux
        use morphevolution
        use update

        implicit none

        real(Kind= Prec) :: dzbs = 0.01

        ! transus
        allocate (C(imax,jmax,gmax),ws(imax,jmax,gmax),Ts(imax,jmax,gmax),Tsg(imax,jmax,gmax))
        allocate (cu(imax,jmax,gmax),cub(imax,jmax,gmax),cv(imax,jmax,gmax),cvb(imax,jmax,gmax))
        allocate (cc(imax,jmax,gmax),ccb(imax,jmax,gmax),ccg(imax,jmax,gmax),ccbg(imax,jmax,gmax))
        allocate(D(gmax),A(gmax),alpha1(gmax),alpha2(gmax),alpha(gmax),wst(gmax),R(gmax),w(gmax),dster(gmax))
        allocate(hloc(imax,jmax),h(imax,jmax),Dmm(imax,jmax),zon(imax,jmax),vmag2(imax,jmax),ustarc(imax,jmax) &
                ,ustarcrit(imax,jmax))
        allocate(taub(imax,jmax),taucrit(imax,jmax),Tstar(imax,jmax),delb(imax,jmax),zos(imax,jmax),z0(imax,jmax))
        allocate(ub_cr(imax,jmax,gmax),us_cr1(imax,jmax,gmax),us_cr2(imax,jmax,gmax),frc(imax,jmax,gmax))
        allocate(pbbed(imax,jmax,lmax,gmax))
        allocate(u(imax,jmax),v(imax,jmax),vmg(imax,jmax),urms(imax,jmax),urms2(imax,jmax),Cd(imax,jmax),dcsdx(imax,jmax))
        allocate(dcsdy(imax,jmax),dcbdy(imax,jmax),dcbdx(imax,jmax))
        allocate(Asb(imax,jmax),wet(imax,jmax),hold(imax,jmax))
        allocate(term1(imax,jmax),term2(imax,jmax),term3(imax,jmax))
        allocate(ceqb(imax,jmax,gmax),ceqs(imax,jmax,gmax),ceq(imax,jmax,gmax),ceqbg(imax,jmax,gmax),ceqsg(imax,jmax,gmax))
        allocate(ero1(imax,jmax,gmax),depo_ex1(imax,jmax,gmax),ero2(imax,jmax,gmax),depo_ex2(imax,jmax,gmax))
        allocate(sus(imax,jmax,gmax),sub(imax,jmax,gmax),svs(imax,jmax,gmax),svb(imax,jmax,gmax))
        allocate(susg(imax,jmax,gmax),subg(imax,jmax,gmax),svsg(imax,jmax,gmax),svbg(imax,jmax,gmax),fac(imax,jmax,gmax))
        allocate(laythick(imax,jmax,lmax))

        !flux
        allocate(V_new(0:imax+1,0:jmax+1,3))
        allocate(Vx_l(imax,jmax,3),Vx_r(imax,jmax,3),Vy_l(imax,jmax,3),Vy_r(imax,jmax,3))
        allocate(Cx_l(imax,jmax,gmax,2),Cx_r(imax,jmax,gmax,2),Cy_l(imax,jmax,gmax,2),Cy_r(imax,jmax,gmax,2))
        allocate(C_new_x(0:imax+1,0:jmax+1,gmax,2),C_new_y(0:imax+1,0:jmax+1,gmax,2))
        allocate(psix1(-1:imax+2,-1:jmax+2,gmax+3),psix2(-1:imax+2,-1:jmax+2,gmax+3),psiy1(-1:imax+2,-1:jmax+2,gmax+3) &
                ,psiy2(-1:imax+2,-1:jmax+2,gmax+3))
        !morph
        allocate(indx(jmax))
        allocate(totalnum(imax,jmax))
        allocate(dzbdt(imax,jmax))
        allocate(indSus(imax,jmax,gmax),indSvs(imax,jmax,gmax),indSub(imax,jmax,gmax),indSvb(imax,jmax,gmax))
        allocate(Sout(imax,jmax,gmax),fre(imax,jmax,gmax),dzg(imax,jmax,gmax))
        allocate(zb(imax,jmax),sedero(imax,jmax))
        allocate(dzbdx(imax,jmax),dzbdy(imax,jmax),totalthick(imax,jmax),zs(imax,jmax),dzav(imax,jmax))
        allocate(dz(lmax),sedcal(gmax))
        allocate(pb(gmax,lmax))
        allocate(dzbed(imax,jmax,lmax))
        allocate(z0bed(imax,jmax))
        allocate(edg(imax,jmax,gmax))

        cu = 0.0001
        cub = 0.0001
        cv = 0.0000
        cvb = 0.0000
        ccg = 0.0001
        ccbg = 0.0001
        cc = ccg
        ccb = ccbg
        pbbed(:,:,1:2,1) = 0.8
        pbbed(:,:,1:2,2) = 0.2
        !do j = 1, 2
        !    do i = 1, imax
        !        pbbed(i,j,1,:) = 1
        !    end do
        !end do
        !do j = 3, jmax
        !    do i = 1, imax
        !        pbbed(i,j,2,:) = 1
        !    end do
        !end do
        D(1) = 0.0003
        D(2) = 0.0001
        urms = 0.0 !long wave
        u = 0.0
        do i = 1, imax
            do j = 1, jmax
                !u(i,j) = 1.0*j
                v(i,j) = 1.0*j
                zb(i,j) = 0.30*j
                h(i,j) = 5*j
            enddo
        enddo
        !z0bed = zb - totalthick
        sedcal = 1.0
        !print *,zb
        !laythick = 0.01
        hold = h
        !dzbed(:,:,1) = 0.05
        !dzbed(:,:,2) = 0.05
        !print *, dzbed
        !call transus
        susg = 0.01
        svsg = 0.01
        subg = 0.01
        svbg = 0.01
        indx = 5
        do i = 1, imax
            do j = 1, jmax
                zb(i,j) =0.5+0.03*i
            end do
        end do

        call transus

        !susg = 0.05
        !svsg = 0.05
        !subg = 0.05
        !svbg = 0.05


        !zb = 1.0
        totalthick = 0.27
        zs = hold
        z0bed = zb - totalthick
        i = 1
        j = 1
        dzbed = 0.0
        pbbed = 0.0
        dzbed(:,:,6)= 0.02
        dzbed(:,:,1)= 0.05
        dzbed(:,:,2:5) = 0.05
        pbbed(:,:,1:6,1) = 0.8
        pbbed(:,:,1:6,2) = 0.2
        !dzg(1,1,1) = 0.08*(1-0.4)*0.5/dt
        !dzg(1,1,2) = 0.08*(1-0.4)*0.5/dt
        nd_var = 6
        !print *,totalthick(i,j)

        !dzbs = 0.08
        !call update_fractions(i,j,dzbed(1,1,:),pbbed(1,1,:,:),dzg(1,1,:),dzbs)
        call bed_update


        !print *, dzg


    end program Test






























