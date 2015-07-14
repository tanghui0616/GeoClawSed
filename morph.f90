        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Copyright (C) 2015 Virginia Tech, Sediment transport Processes Group    !
        ! Hui Tang, Wei Cheng and Robert Weiss                                    !
        !                                                                         !
        ! tanghui@vt.edu                                                          !
        !                                                                         !
        ! This library is free software; you can redistribute it and/or           !
        ! modify it under the terms of the GNU Lesser General Public              !
        ! License as published by the Free Software Foundation; either            !
        ! version 2.1 of the License, or (at your option) any later version.      !
        !                                                                         !
        ! This library is distributed in the hope that it will be useful,         !
        ! but WITHOUT ANY WARRANTY; without even the implied warranty of          !
        ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        !
        ! Lesser General Public License for more details.                         !
        !                                                                         !
        ! You should have received a copy of the GNU Lesser General Public        !
        ! License along with this library; if not, write to the Free Software     !
        ! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307     !
        ! USA                                                                     !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        module morphevolution

            use Set_precision
            use params
            use sed
            use update

            IMPLICIT NONE

            contains

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !This part is used for containting the bed updating part                                    !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            subroutine bed_update

                use params

                IMPLICIT NONE

                dzbdt  = 0.0
                !print *, pbbed

                if (t>=morstart .and. morfac > .9990) then

                    if (struct == 1) then
                        indSus = 0
                        indSub = 0
                        indSvs = 0
                        indSvb = 0
                        Sout   = 0.0
                        do k = 1,gmax
                            do j=1,jmax
                                do i=2,imax
                                    ! fluxes at i,j
                                    if (Subg(i,j,k) > 0.0) then      ! bed load x-direction
                                        indSub(i,j,k) = 1
                                        Sout(i,j,k) = Sout(i,j,k) + Subg(i,j,k)*dx
                                    endif
                                    ! fluxes at i-1,j
                                    if (Subg(i-1,j,k) < 0.0 ) then   ! bed load x-direction
                                        Sout(i,j,k) = Sout(i,j,k) - Subg(i-1,j,k)*dx
                                    endif
                                    if (sourcesink==0) then
                                    ! fluxes at i,j
                                        if (Susg(i,j,k) > 0.0 ) then     ! suspended load x-direction
                                            indSus(i,j,k) = 1
                                            Sout(i,j,k) = Sout(i,j,k) + Susg(i,j,k)*dx
                                        endif
                                        ! fluxes at i-1,j
                                        if (Susg(i-1,j,k) < 0.0 ) then   ! suspended load x-direction
                                            Sout(i,j,k) = Sout(i,j,k) - Susg(i-1,j,k)*dx
                                        endif
                                    endif !sourcesink==0
                                enddo! imax
                            enddo! jmax
                            if (jmax>0) then
                                do j=2,jmax
                                    do i=1,imax
                                        if (Svbg(i,j,k) > 0.0 ) then     ! bed load y-direction
                                            indSvb(i,j,k) = 1
                                            Sout(i,j,k) = Sout(i,j,k) + Svbg(i,j,k)*dy
                                        endif
                                        ! fluxes at i,j-1
                                        if (Svbg(i,j-1,k) < 0.0 ) then   ! bed load y-direction
                                            Sout(i,j,k) = Sout(i,j,k) - Svbg(i,j-1,k)*dy
                                        endif
                                        if (sourcesink==0) then
                                            if (Svsg(i,j,k) > 0.0 ) then     ! suspended load y-direction
                                                indSvs(i,j,k) = 1
                                                Sout(i,j,k) = Sout(i,j,k) + Svsg(i,j,k)*dy
                                            endif
                                            ! fluxes at i,j-1
                                            if (Svsg(i,j-1,k) < 0.0 ) then   ! suspended load y-direction
                                                Sout(i,j,k) = Sout(i,j,k) - Svsg(i,j-1,k)*dy
                                            endif
                                        endif ! sourcesink = 0
                                    enddo !imax
                                enddo !jmax
                            endif !jmax>0
                            ! reduce the sediment flux/ transport due to hard structure layers
                            do j=1,jmax
                                do i=1,imax
                                    ! reduction factor for cell outgoing sediment transports
                                    Savailable = laythick(i,j,1)*pbbed(i,j,1,k)/morfac/dt*(1.0-por)*dx*dy
                                    fre(i,j,k)  = min(1.0,Savailable/max(Sout(i,j,k),tiny(0.0)) )
                                enddo
                                do i=1,imax
                                    ! fix sediment transports for the presence of a hard layer; remind indSus etc are 1 in cases of cell outgoing transports
                                        ! updated S         all outgoing transports                  cell incoming transports
                                    if (fre(i,j,k) < 1.0)then
                                        Subg(i,j,k)   = fre(i,j,k)*indSub(i,j,k)*Subg(i,j,k) &
                                                        + (1-indSub(i,j,k))*Subg(i,j,k)
                                        Subg(i-1,j,k) = fre(i-1,j,k)*(1-indSub(i-1,j,k))*Subg(i-1,j,k) &
                                                        + indSub(i-1,j,k)*Subg(i-1,j,k)
                                        if (jmax>0) then
                                            Svbg(i,j,k)   = fre(i,j,k)*indSvb(i,j,k)*Svbg(i,j,k) &
                                                        + (1-indSvb(i,j,k))*Svbg(i,j,k)
                                            Svbg(i,j-1,k) = fre(i,j-1,k)*(1-indSvb(i,j-1,k))*Svbg(i,j-1,k) &
                                                        + indSvb(i,j-1,k)*Svbg(i,j-1,k)
                                        endif
                                        if (sourcesink==0) then
                                            Susg(i,j,k)   = fre(i,j,k)*indSus(i,j,k)*Susg(i,j,k) &
                                                        + (1-indSus(i,j,k))*Susg(i,j,k)
                                            Susg(i-1,j,k) = fre(i-1,j,k)*(1-indSus(i-1,j,k))*Susg(i-1,j,k) &
                                                        + indSus(i-1,j,k)*Susg(i-1,j,k)
                                            if (jmax>0) then
                                                Svsg(i,j,k)   = fre(i,j,k)*indSvs(i,j,k)*Svsg(i,j,k) &
                                                        + (1-indSvs(i,j,k))*Svsg(i,j,k)
                                                Svsg(i,j-1,k) = fre(i,j,k)*(1-indSvs(i,j-1,k))*Svsg(i,j-1,k) &
                                                        + indSvs(i,j-1,k)*Svsg(i,j-1,k)
                                            endif !jmax > 0
                                        endif ! sourcesink = 0
                                    endif !fac<1.0
                                enddo ! imax
                            enddo !jmax
                        enddo !gmax
                    endif !struct == 1
                    if (jmax>0) then
                        do j=1,jmax
                            do i=1,imax
                                !print *, i,j
                                !print *, dzbed(i,j,:)
                            ! bed level changes per fraction in this morphological time step in meters sand including pores
                            ! positive in case of erosion
                                if (sourcesink==0) then
                                    dzg(i,j,:)=morfac*dt/(1.0-por)*( &
                                        ! dz from sus transport gradients
                                        Susg(i,j,:)*dx*h(i,j)-Susg(i-1,j,:)*dx*h(i-1,j) +&
                                        Svsg(i,j,:)*dy*h(i,j)-Svsg(i,j-1,:)*dy*h(i,j-1) +&
                                        ! dz from bed load transport gradients
                                        Subg(i,j,:)*dx*h(i,j)-Subg(i-1,j,:)*dx*h(i-1,j) +&
                                        Svbg(i,j,:)*dy*h(i,j)-Svbg(i,j-1,:)*dy*h(i,j-1)  +&
                                        !source term
                                        (cu(i,j,:)+cv(i,j,:)+cub(i,j,:)+cvb(i,j,:)-ceqs(i,j,:)-ceqb(i,j,:))*h(i,j)/Tsg(i,j,:))
                                elseif (sourcesink==1) then
                                    dzg(i,j,:)=morfac*dt/(1.0-por)*( &
                                        ! dz from sus transport gradients
                                        Susg(i,j,:)*dx*h(i,j)-Susg(i-1,j,:)*dx*h(i-1,j) +&
                                        Svsg(i,j,:)*dy*h(i,j)-Svsg(i,j-1,:)*dy*h(i,j-1) +&
                                        !source term
                                        (cu(i,j,:)+cv(i,j,:)-ceqs(i,j,:))*h(i,j)/Tsg(i,j,:))
                                endif
                                if (gmax==1) then ! Simple bed update in case one fraction
                                    zb(i,j) = zb(i,j)-sum(dzg(i,j,:))
                                    dzbdt(i,j) = dzbdt(i,j)-sum(dzg(i,j,:))
                                    sedero(i,j) = sedero(i,j)-sum(dzg(i,j,:))
                                    totalthick(i,j) = max(0.0,totalthick(i,j)-sum(dzg(i,j,:)))
                                else
                                    edg(i,j,:) = dzg(i,j,:)*(1.0-por)/dt
                                    if (totalthick(i,j)>thick) then
                                        totalnum(i,j) = nint(totalthick(i,j)/thick)
                                        if (totalthick(i,j)>totalnum(i,j)*thick) then
                                            totalnum(i,j)=totalnum(i,j)+1
                                        endif
                                    endif
                                    if(dzbed(i,j,totalnum(i,j))<toler*thick) then
                                        totalnum(i,j) = totalnum(i,j) - 1
                                    endif
                                    nd_var = totalnum(i,j)
                                    !print *, i,j
                                    !print *, dzbed(i,j,:)
                                    !print *, pbbed(i,j,:,:)
                                    !print *, edg(i,j,:)
                                    !print *, sum(dzg(i,j,:))
                                    !print *, nd_var
                                    call update_fractions(i,j,dzbed(i,j,:),pbbed(i,j,:,:),edg(i,j,:),sum(dzg(i,j,:)))
                                    !nd_var = maxval(totalnum)
                                    !print *, "**************************************************"
                                    !print *, dzbed(i,j,:)
                                    !print *, pbbed(i,j,:,:)

                                endif
                            enddo ! imax+1
                        enddo ! jmax+1
                    else
                        j=1
                        do i=1,imax
                            ! bed level changes per fraction in this morphological time step in meters sand including pores
                            ! positive in case of erosion
                            if (sourcesink==0) then
                                dzg(i,j,:)=morfac*dt/(1.0-por)*( &
                                    ! dz from sus transport gradients
                                    Susg(i,j,:)*dx*h(i,j)-Susg(i-1,j,:)*dx*h(i-1,j) +&
                                    ! dz from bed load transport gradients
                                    Subg(i,j,:)*dx*h(i,j)-Subg(i-1,j,:)*dx*h(i-1,j) +&
                                    !source term
                                    (cu(i,j,:)+cub(i,j,:)-ceqs(i,j,:)-ceqb(i,j,:))*h(i,j)/Tsg(i,j,:))
                            elseif (sourcesink==1) then
                                dzg(i,j,:)=morfac*dt/(1.0-por)*( &
                                    ! dz from sus transport gradients
                                    Susg(i,j,:)*dx*h(i,j)-Susg(i-1,j,:)*dx*h(i-1,j) +&
                                    !source term
                                    (cu(i,j,:)-ceqs(i,j,:))*h(i,j)/Tsg(i,j,:))
                            endif
                            if (gmax==1) then ! Simple bed update in case one fraction
                                zb(i,j) = zb(i,j)-sum(dzg(i,j,:))
                                dzbdt(i,j) = dzbdt(i,j)-sum(dzg(i,j,:))
                                sedero(i,j) = sedero(i,j)-sum(dzg(i,j,:))
                                totalthick(i,j) = max(0.0,totalthick(i,j)-sum(dzg(i,j,:)))
                            else ! multiple fractions...
                                edg(i,j,:) = dzg(i,j,:)*(1.0-por)/dt
                                if (totalthick(i,j)>thick) then
                                    totalnum(i,j) = nint(totalthick(i,j)/thick)
                                    if (totalthick(i,j)>totalnum(i,j)*thick) then
                                        totalnum(i,j)=totalnum(i,j)+1
                                    endif
                                endif
                                if(dzbed(i,j,totalnum(i,j))<toler*thick) then
                                    totalnum(i,j) = totalnum(i,j) - 1
                                endif
                                nd_var = totalnum(i,j)
                                call update_fractions(i,j,dzbed(i,j,:),pbbed(i,j,:,:),edg(i,j,:),sum(dzg(i,j,:)))
                            endif
                        enddo ! imax
                    endif !jmax = 1
                endif !t
                !print *, zb
                !print *, pbbed(:,:,1,1)
                call avalanch
                !print *, zb
                !print *, pbbed(:,:,1,1)

            end subroutine bed_update

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! This part is designed for avanlanching scheme
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            subroutine avalanch

                use params

                IMPLICIT NONE

                Integer         :: ie,id,jdz,je,jd
                Real(kind=Prec) :: sdz,dzb,one=1.00
                Real(kind=Prec),dimension(gmax) :: edg1,edg2
                !Real(kind=Prec),dimension(gmax,lmax) :: pb

                if (avalanching==1) then

                    do ii=1,nint(morfac)
                        aval=.false.
                        dzbdx=0.0
                        dzbdy=0.0
                        do j=1,jmax
                            do i=1,imax-1
                                dzbdx(i,j)=(zb(i+1,j)-zb(i,j))/dx
                            enddo
                        enddo
                        do j=1,jmax-1
                            do i=1,imax
                                dzbdy(i,j)=(zb(i,j+1)-zb(i,j))/dy
                            enddo
                        enddo
                        !print *, dzbdx
                        !print *, dzbdy
                        indx = imax+1
                        do i=1,imax-1
                            do j=1,jmax
                                !decide Maximum bedlevel change due to avalanching
                                if(max(h(i,j),h(i+1,j))>hswitch+eps) then ! Jaap instead of hh
                                    dzmax=wetslp
                                    if (i>indx(j)) then ! tricks: seaward of indx (transition from sand to structure) wetslope is set to 0.03; TODO
                                        dzmax = max(wetslp,abs(dzbdx(i,j))*0.990)
                                    endif
                                else
                                    dzmax=dryslp
                                endif
                                if(abs(dzbdx(i,j))>dzmax .and. h(i+nint(max(0.0,sign(one,dzbdx(i,j)))),j)>eps) then
                                    aval=.true.
                                    ! the total mass need to be moved
                                    dzb=sign(one,dzbdx(i,j))*(abs(dzbdx(i,j))-dzmax)*dx
                                    if (dzb >= 0.0) then
                                        ie = i+1                                        ! index erosion point
                                        id = i                                          ! index deposition point
                                        dzb=min(dzb,dzmax*dt/dx)                        ! make sure dzb is not in conflict with maximum erosion rate dzmax
                                        dzb=min(dzb,totalthick(i+1,j))                 ! make sure dzb is not larger than sediment layer thickness
                                    else
                                        ie = i                                          ! index erosion point
                                        id = i+1                                        ! index deposition point
                                        dzb=max(dzb,-dzmax*dt/dx)
                                        dzb=max(dzb,-totalthick(i,j))
                                    endif
                                    if (gmax == 1) then ! Simple bed update in case one fraction
                                        dzleft = abs(dzb)
                                        zb(id,j) = zb(id,j)+dzleft
                                        zb(ie,j) = zb(ie,j)-dzleft
                                        dzbdt(id,j) = dzbdt(id,j)+dzleft
                                        dzbdt(ie,j) = dzbdt(ie,j)-dzleft
                                        sedero(id,j) = sedero(id,j)+dzleft
                                        sedero(ie,j) = sedero(ie,j)-dzleft
                                        totalthick(id,j) = max(0.0,totalthick(id,j)+dzleft)
                                        totalthick(ie,j) = max(0.0,totalthick(ie,j)-dzleft)
                                        zs(id,j)  = zs(id,j)+dzleft
                                        zs(ie,j)  = zs(ie,j)-dzleft
                                        dzav(id,j)= dzav(id,j)+dzleft
                                        dzav(ie,j)= dzav(ie,j)-dzleft
                                    else ! multiple fractions...

                                        ! now fix fractions....
                                        !dz = dzbed(ie,j,:)
                                        !pb = pbbed(ie,j,:,:)
                                        ! figure out how many sediment layers (ndz) are eroded in point iii
                                        sdz = 0.0
                                        ndz = 0
                                        do while (sdz<abs(dzb))
                                            ndz = ndz+1
                                            sdz = sdz+dzbed(ie,j,ndz)
                                        enddo
                                        !print *, ndz
                                        ! now update bed and fractions by stepping through each layer seperately
                                        dzleft = abs(dzb)
                                        dzavt  = 0.0
                                        do jdz=1,ndz
                                            dzt = min(dzbed(ie,j,jdz),dzleft)
                                            dzleft = dzleft-dzt
                                            ! erosion deposition per fraction (upwind or downwind); edg is positive in case of erosion
                                            do k=1,gmax
                                                edg2(k) =  sedcal(k)*dzt*pbbed(ie,j,jdz,k)*(1.0-por)/dt                ! erosion
                                                edg1(k) = -sedcal(k)*edg2(k)                                   ! deposition
                                            enddo
                                            dzavt = dzavt + sum(edg2)*dt/(1.0-por)
                                            nd_var = totalnum(ie,j)
                                            !print *,ie,j
                                            !print *,dzbed(ie,j,:)
                                            !print *,pbbed(ie,j,:,:)
                                            !print *,edg2
                                            !print *,dzavt
                                            call update_fractions(ie,j,dzbed(ie,j,:),pbbed(ie,j,:,:),edg2,dzavt)! update bed in eroding point
                                            nd_var = totalnum(id,j)
                                            call update_fractions(id,j,dzbed(id,j,:),pbbed(id,j,:,:),edg1,-dzavt) ! update bed in deposition point
                                            nd_var = maxval(totalnum)
                                        enddo
                                        ! update water levels and dzav
                                        zs(ie,j)  = zs(ie,j)-dzavt
                                        dzav(ie,j)= dzav(ie,j)-dzavt
                                        zs(id,j)  = zs(id,j)+dzavt
                                        dzav(id,j)= dzav(id,j)+dzavt
                                    end if ! yes/no multiple fractions
                                end if !dzmax
                            end do !jmax
                        end do !imax
                        !JJ: update y slopes after avalanching in X-direction seems more appropriat

                        do j=1,jmax-1 !
                            do i=1,imax
                                if(max(h(i,j),h(i,j+1))>hswitch+eps) then
                                    dzmax=wetslp
                                else
                                    dzmax=dryslp
                                end if

                                if(abs(dzbdy(i,j))>dzmax .and. h(i,j+nint(max(0.0,sign(one,dzbdy(i,j)))))>eps) then

                                    aval=.true.
                                    dzb=sign(one,dzbdy(i,j))*(abs(dzbdy(i,j))-dzmax)*dy
                                    if (dzb >= 0.0) then
                                            je = j+1                                        ! index erosion point
                                            jd = j                                          ! index deposition point
                                            dzb=min(dzb,dzmax*dt/dy)
                                            dzb=min(dzb,totalthick(i,j+1))
                                    else
                                            je = j                                          ! index erosion point
                                            jd = j+1                                        ! index deposition point
                                            dzb=max(dzb,-dzmax*dt/dy)
                                            dzb=max(dzb,-totalthick(i,j))
                                    endif
                                    if (gmax == 1) then ! Simple bed update in case one fraction
                                        dzleft = abs(dzb)
                                        zb(i,jd) = zb(i,jd)+dzleft
                                        zb(i,je) = zb(i,je)-dzleft
                                        dzbdt(i,jd) = dzbdt(i,jd)+dzleft
                                        dzbdt(i,je) = dzbdt(i,je)-dzleft
                                        sedero(i,jd) = sedero(i,jd)+dzleft
                                        sedero(i,je) = sedero(i,je)-dzleft
                                        totalthick(i,jd) = max(0.0,totalthick(i,jd)+dzleft)
                                        totalthick(i,je) = max(0.0,totalthick(i,je)-dzleft)
                                        zs(i,jd)  = zs(i,jd)+dzleft
                                        zs(i,je)  = zs(i,je)-dzleft
                                        dzav(i,jd)= dzav(i,jd)+dzleft
                                        dzav(i,je)= dzav(i,je)-dzleft
                                    else ! multiple fractions...\

                                        ! figure out how many depth layers (ndz) are affected
                                        sdz = 0
                                        ndz = 0

                                        do while (sdz<abs(dzb))
                                            ndz = ndz+1
                                            sdz = sdz+dzbed(i,je,ndz)
                                        enddo
                                        ! now update bed and fractions by stepping through each layer seperately
                                        dzleft = abs(dzb)
                                        dzavt  = 0.0
                                        do jdz=1,ndz
                                            dzt = min(dzbed(i,je,jdz),dzleft)
                                            dzleft = dzleft-dzt;
                                            !print *, dzbed(i,je,jdz)
                                            ! erosion deposition per fraction (upwind or downwind); edg is positive in case of erosion
                                            do k=1,gmax
                                                edg2(k) = sedcal(k)*dzt*pbbed(i,je,jdz,k)*(1.0-por)/dt        ! erosion
                                                edg1(k) = -sedcal(k)*edg2(k)                            ! deposition
                                            enddo
                                            dzavt = dzavt + sum(edg2)*dt/(1.0-por)
                                            nd_var = totalnum(i,je)
                                            !print *,i,je
                                            !print *,dzbed(i,je,:)
                                            !print *, pbbed(i,je,:,:)
                                            !print *,edg2
                                            !print *,dzt
                                            call update_fractions(i,je,dzbed(i,je,:),pbbed(i,je,:,:),edg2,dzavt)           ! upwind point
                                            nd_var = totalnum(i,jd)
                                            call update_fractions(i,jd,dzbed(i,jd,:),pbbed(i,jd,:,:),edg1,-dzavt)    ! downwind point
                                            nd_var = maxval(totalnum)
                                        enddo
                                        ! update water levels and dzav
                                        zs(i,je)  = zs(i,je)-dzavt
                                        dzav(i,je)= dzav(i,je)-dzavt
                                        zs(i,jd)  = zs(i,jd)+dzavt
                                        dzav(i,jd)= dzav(i,jd)+dzavt
                                    endif !yes/no multiple fractions
                                end if !dzmax
                            end do !imax
                        end do !jmax
                        if (.not.aval) exit
                    end do !morfac
                end if ! avalanching
            end subroutine avalanch


        end module morphevolution









