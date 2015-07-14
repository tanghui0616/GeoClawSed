!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This Part is used to calculate sediment Flux for finite volume method              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module Flux

        use params
        use Set_Precision

        implicit none

        contains
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !This part is used to caculate flux term                                             !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine Flux_vector

            use params
            use Set_Precision

            implicit none
            integer :: jg

            !flux limitor

            V_new(1:imax,1:jmax,1) = u(:,:)
            V_new(1:imax,1:jmax,2) = v(:,:)
            V_new(1:imax,1:jmax,3) = h(:,:)
            C_new_x(1:imax,1:jmax,:,1)= cc
            C_new_x(1:imax,1:jmax,:,2)= ccb
            C_new_y(1:imax,1:jmax,:,1)= cc
            C_new_y(1:imax,1:jmax,:,2)= ccb
            V_new(0,:,:) = V_new(1,:,:)
            V_new(imax+1,:,:) = V_new(imax,:,:)
            V_new(:,0,:) = V_new(:,1,:)
            V_new(:,jmax+1,:) = V_new(:,jmax,:)
            C_new_x(0,:,:,:) = C_new_x(1,:,:,:)
            C_new_x(imax+1,:,:,:) = C_new_x(imax,:,:,:)
            C_new_x(:,0,:,:) = C_new_x(:,1,:,:)
            C_new_x(:,jmax+1,:,:) =C_new_x(:,jmax,:,:)
            C_new_y(0,:,:,:) = C_new_y(1,:,:,:)
            C_new_y(imax+1,:,:,:) = C_new_y(imax,:,:,:)
            C_new_y(:,0,:,:) = C_new_y(:,1,:,:)
            C_new_y(:,jmax+1,:,:) =C_new_y(:,jmax,:,:)

            !flux limitor
            call flux_limitor

            do i=1, imax
                do j=1, jmax
                    do jg = 1, 3
                        Vx_l(i,j,jg) = V_new(i-1,j,jg) + 0.25*vareps* &
                                        ((1-k1)*psix1(i-1,j,jg)*(V_new(i,j,jg)-V_new(i-1,j,jg))+ &
                                        (1+k1)*psix2(i,j,jg)*(V_new(i+1,j,jg)-V_new(i,j,jg))      )
                        Vx_r(i,j,jg) = V_new(i,j,jg)  + 0.25*vareps* &
                                        ((1-k1)*psix2(i+1,j,jg)*(V_new(i+2,j,jg)-V_new(i+1,j,jg)) + &
                                        (1+k1)*psix1(i,j,jg)*(V_new(i+1,j,jg)-V_new(i,j,jg))      )
                        Vy_l(i,j,jg) = V_new(i,j-1,jg) + 0.25*vareps* &
                                        ((1-k1)*psiy1(i,j-1,jg)*(V_new(i,j,jg)-V_new(i,j-1,jg))+ &
                                        (1+k1)*psiy2(i,j,jg)*(V_new(i+1,j,jg)-V_new(i,j,jg))      )
                        Vy_r(i,j,jg) = V_new(i,j,jg)  + 0.25*vareps* &
                                        ((1-k1)*psiy2(i,j+1,jg)*(V_new(i,j+2,jg)-V_new(i,j+1,jg)) + &
                                        (1+k1)*psiy1(i,j,jg)*(V_new(i,j+1,jg)-V_new(i,j,jg))      )
                    end do
                    do jg=1,gmax
                            Cx_l(i,j,jg,:) = C_new_x(i-1,j,jg,:) + 0.25*vareps* &
                                        ((1-k1)*psix1(i-1,j,jg+3)*(C_new_x(i,j,jg,:)-C_new_x(i-1,j,jg,:))+ &
                                        (1+k1)*psix2(i,j,jg+3)*(C_new_x(i+1,j,jg,:)-C_new_x(i,j,jg,:))      )
                            Cx_r(i,j,jg,:) = C_new_x(i,j,jg,:)  + 0.25*vareps* &
                                        ((1-k1)*psix2(i+1,j,jg+3)*(C_new_x(i+2,j,jg,:)-C_new_x(i+1,j,jg,:)) + &
                                        (1+k1)*psix1(i,j,jg+3)*(C_new_x(i+1,j,jg,:)-C_new_x(i,j,jg,:))      )
                            Cy_l(i,j,jg,:) = C_new_y(i,j-1,jg,:) + 0.25*vareps* &
                                        ((1-k1)*psiy1(i,j-1,jg+3)*(C_new_y(i,j,jg,:)-C_new_y(i,j-1,jg,:))+ &
                                        (1+k1)*psiy2(i,j,jg+3)*(C_new_y(i+1,j,jg,:)-C_new_y(i,j,jg,:))      )
                            Cy_r(i,j,jg,:) = C_new_y(i,j,jg,:)  + 0.25*vareps* &
                                        ((1-k1)*psiy2(i,j+1,jg+3)*(C_new_y(i,j+2,jg,:)-C_new_y(i,j+1,jg,:)) + &
                                        (1+k1)*psiy1(i,j,jg+3)*(C_new_y(i,j+1,jg,:)-C_new_y(i,j,jg,:))      )
                    end do
                end do
            end do
            do i=1, imax
                do j=1, jmax
                    if (method == 'SVL') then

                        Sus(i,j,:)  = flux_calc_SVL(Vx_l(i,j,1),Vx_l(i,j,2),Vx_l(i,j,3),Cx_l(i,j,:,:),Vx_r(i,j,1), &
                                        Vx_r(i,j,2),Vx_r(i,j,3),Cx_r(i,j,:,:))
                        Svs(i,j,:)  = flux_calc_SVL(Vy_l(i,j,1),Vy_l(i,j,2),Vy_l(i,j,3),Cy_l(i,j,:,:),Vy_r(i,j,1), &
                                        Vy_r(i,j,2),Vy_r(i,j,3),Cy_r(i,j,:,:))
                        Sub(i,j,:)  = flux_calc_SVL(Vx_l(i,j,1),Vx_l(i,j,2),Vx_l(i,j,3),Cx_l(i,j,:,:),Vx_r(i,j,1), &
                                        Vx_r(i,j,2),Vx_r(i,j,3),Cx_r(i,j,:,:))
                        Svb(i,j,:)  = flux_calc_SVL(Vy_l(i,j,1),Vy_l(i,j,2),Vy_l(i,j,3),Cy_l(i,j,:,:),Vy_r(i,j,1), &
                                        Vy_r(i,j,2),Vy_r(i,j,3),Cy_r(i,j,:,:))
                    else
                        print *, 'WARNING: DECIDE TURN ON OR OFF FLUX LIMITER'
                    endif
                enddo
            enddo

        end subroutine Flux_vector
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !This part is used to Standard Van Leer method
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function flux_calc_SVL(u_l,v_l,h_l,c_l,u_r,v_r,h_r,c_r)

            implicit none

            Real(kind=Prec)  :: u_l,v_l,h_l,u_r,v_r,h_r
            Real(kind=Prec)  :: speed_l,speed_r,fr_l,fr_r,Fr1,Fr2,Beta_r,Beta_l,Alpha1,Alpha2,C1,C2,one
            Real(kind=Prec), Dimension(gmax) :: flux_calc_SVL
            Real(kind=Prec), Dimension(gmax) :: c_l, c_r
            integer :: jg
            one = 1.0
            speed_l = dsqrt(dabs(g*h_l))
            speed_r = dsqrt(dabs(g*h_r))
            fr_l     = u_l/speed_l
            fr_r     = u_r/speed_r
            Fr1      = 0.25*(fr_l+1.0)**2.0
            Fr2      = -0.25*(fr_r-1.0)**2.0
            Beta_l  = -max(0,1-INT(dabs(fr_l)))
            Beta_r  = -max(0,1-INT(dabs(fr_r)))
            Alpha1  = 0.5*(1.0+SIGN(one,fr_l))
            Alpha2  = 0.5*(1.0-SIGN(one,fr_r))
            C1      = Alpha1*(1.0+Beta_l)*fr_l-Beta_l*Fr1
            C2      = Alpha2*(1.0+Beta_r)*fr_r-Beta_r*Fr2
            do jg = 1, gmax
                flux_calc_SVL(jg) = c_l(jg)*speed_l*C1+c_r(jg)*speed_r*C2
            end do
        end function flux_calc_SVL
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !This part is used to put flux limiter
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine flux_limitor

            implicit none

            Real(kind=Prec),Dimension(1:imax-1,1:jmax-1,gmax+3) :: Rx1,Rx2 !
            Real(kind=Prec),Dimension(1:imax-1,1:jmax-1,gmax+3) :: Ry1,Ry2  !
            Real(kind=Prec),Dimension(1:imax-1,1:jmax-1,gmax+3) :: denx, deny
            Real(kind=Prec),Dimension(0:imax,0:jmax,gmax) :: C_new_x_t,C_new_y_t
            integer :: ii

            do i = 1, imax
                do j = 1, jmax
                    do ii = 1, gmax
                        C_new_x_t(i,j,ii) = sum(C_new_x(i,j,ii,:))
                        C_new_y_t(i,j,ii) = sum(C_new_y(i,j,ii,:))
                    enddo
                enddo
            enddo

            do i = 1, imax-1
                do j = 1, jmax-1
                    do ii = 1, 3
                        denx(i,j,ii) = V_new(i+1,j,ii)-V_new(i,j,ii)
                        denx(i,j,ii) = sign(max(abs(deny(i,j,ii)),toler),deny(i,j,ii))
                        deny(i,j,ii) = V_new(i,j+1,ii)-V_new(i,j,ii)
                        deny(i,j,ii) = sign(max(abs(deny(i,j,ii)),toler),deny(i,j,ii))
                        Rx1(i,j,ii) = (V_new(i+2,j,ii)-V_new(i+1,j,ii))/denx(i,j,ii)
                        Rx2(i,j,ii) = (V_new(i,j,ii)-V_new(i-1,j,ii))/denx(i,j,ii)
                        Ry1(i,j,ii) = (V_new(i,j+2,ii)-V_new(i,j+1,ii))/deny(i,j,ii)
                        Ry2(i,j,ii) = (V_new(i,j,ii)-V_new(i,j-1,ii))/deny(i,j,ii)
                    end do
                    do ii = 1, gmax
                        denx(i,j,ii+3) = C_new_x_t(i+1,j,ii)-C_new_x_t(i,j,ii)
                        denx(i,j,ii+3) = sign(max(abs(denx(i,j,ii)),toler),denx(i,j,ii))
                        deny(i,j,ii+3) = C_new_y_t(i,j+1,ii)-C_new_y_t(i,j,ii)
                        deny(i,j,ii+3) = sign(max(abs(deny(i,j,ii)),toler),deny(i,j,ii))
                        Rx1(i,j,ii+3) = (C_new_x_t(i+2,j,ii)-C_new_x_t(i+1,j,ii))/denx(i,j,ii)
                        Rx2(i,j,ii+3) = (C_new_x_t(i,j,ii)-C_new_x_t(i-1,j,ii))/denx(i,j,ii)
                        Ry1(i,j,ii+3) = (C_new_y_t(i,j+2,ii)-C_new_y_t(i,j+1,ii))/deny(i,j,ii)
                        Ry2(i,j,ii+3) = (C_new_y_t(i,j,ii)-C_new_y_t(i,j-1,ii))/deny(i,j,ii)
                    end do
                    if (limit_method == 'Vanleer') then
                        psix1(i,j,:) = (Rx1(i,j,:)+abs(Rx1(i,j,:)))/(1.0+Rx1(i,j,:))
                        psix2(i,j,:) = (Rx2(i,j,:)+abs(Rx2(i,j,:)))/(1.0+Rx2(i,j,:))
                        psiy1(i,j,:) = (Ry1(i,j,:)+abs(Ry1(i,j,:)))/(1.0+Ry1(i,j,:))
                        psiy2(i,j,:) = (Ry2(i,j,:)+abs(Ry2(i,j,:)))/(1.0+Ry2(i,j,:))
                    elseif (limit_method == 'VanAlbada') then
                        psix1(i,j,:) = (Rx1(i,j,:)+(Rx1(i,j,:))**2.0)/(1.0+(Rx1(i,j,:))**2.0)
                        psix2(i,j,:) = (Rx2(i,j,:)+(Rx2(i,j,:))**2.0)/(1.0+(Rx2(i,j,:))**2.0)
                        psiy1(i,j,:) = (Ry1(i,j,:)+(Ry1(i,j,:))**2.0)/(1.0+(Ry1(i,j,:))**2.0)
                        psiy2(i,j,:) = (Ry2(i,j,:)+(Ry2(i,j,:))**2.0)/(1.0+(Ry2(i,j,:))**2.0)
                    elseif (limit_method == 'minmod') then
                        do ii = 1, gmax+3
                            if (Rx1(i,j,ii)>0) then
                                psix1(i,j,ii) = min(1.0, Rx1(i,j,ii))
                            else
                                psix1(i,j,ii) = 0.0
                            endif
                            if (Rx2(i,j,ii)>0) then
                                psix2(i,j,ii) = min(1.0, Rx2(i,j,ii))
                            else
                                psix2(i,j,ii) = 0.0
                            endif
                            if (Ry1(i,j,ii)>0) then
                                psiy1(i,j,ii) = min(1.0, Ry1(i,j,ii))
                            else
                                psiy1(i,j,ii) = 0.0
                            endif
                            if (Ry2(i,j,ii)>0) then
                                psiy2(i,j,ii) = min(1.0, Ry2(i,j,ii))
                            else
                                psiy2(i,j,ii) = 0.0
                            endif
                        enddo
                    elseif(limit_method == 'beta_limiter') then
                        do ii = 1, 3
                            psix1(i,j,ii) = max(0.0,min(beta*Rx1(i,j,ii),1.0),min(Rx1(i,j,ii),beta))
                            psix2(i,j,ii) = max(0.0,min(beta*Rx2(i,j,ii),1.0),min(Rx2(i,j,ii),beta))
                            psiy1(i,j,ii) = max(0.0,min(beta*Ry1(i,j,ii),1.0),min(Ry1(i,j,ii),beta))
                            psiy2(i,j,ii) = max(0.0,min(beta*Ry2(i,j,ii),1.0),min(Ry2(i,j,ii),beta))
                        enddo
                    endif
                end do
            end do
            psix1(0,:,:)=psix1(1,:,:)
            psix1(-1,:,:)=psix1(0,:,:)
            psix1(imax,:,:)=psix1(imax-1,:,:)
            psix1(imax+1,:,:)=psix1(imax,:,:)
            psix1(imax+2,:,:)=psix1(imax+1,:,:)
            psix2(0,:,:)=psix2(1,:,:)
            psix2(-1,:,:)=psix2(0,:,:)
            psix2(imax,:,:)=psix2(imax-1,:,:)
            psix2(imax+1,:,:)=psix2(imax,:,:)
            psix2(imax+2,:,:)=psix1(imax+1,:,:)
            psix1(:,0,:)=psix1(:,1,:)
            psix1(:,-1,:)=psix1(:,0,:)
            psix1(:,jmax,:)=psix1(:,jmax-1,:)
            psix1(:,jmax+1,:)=psix1(:,jmax,:)
            psix1(:,jmax+2,:)=psix1(:,jmax+1,:)
            psix2(:,0,:)=psix2(:,1,:)
            psix2(:,-1,:)=psix2(:,0,:)
            psix2(:,jmax,:)=psix2(:,jmax-1,:)
            psix2(:,jmax+1,:)=psix2(:,jmax,:)
            psix2(:,jmax+2,:)=psix2(:,jmax+1,:)
            psiy1(0,:,:)=psiy1(1,:,:)
            psiy1(-1,:,:)=psiy1(0,:,:)
            psiy1(imax,:,:)=psiy1(imax-1,:,:)
            psiy1(imax+1,:,:)=psiy1(imax,:,:)
            psiy1(imax+2,:,:)=psiy1(imax+1,:,:)
            psiy2(0,:,:)=psiy2(1,:,:)
            psiy2(-1,:,:)=psiy2(0,:,:)
            psiy2(imax,:,:)=psiy2(imax-1,:,:)
            psiy2(imax+1,:,:)=psiy2(imax,:,:)
            psiy2(imax+2,:,:)=psiy2(imax+1,:,:)
            psiy1(:,0,:)=psiy1(:,1,:)
            psiy1(:,-1,:)=psiy1(:,0,:)
            psiy1(:,jmax,:)=psiy1(:,jmax-1,:)
            psiy1(:,jmax+1,:)=psiy1(:,jmax,:)
            psiy1(:,jmax+2,:)=psiy1(:,jmax+1,:)
            psiy2(:,0,:)=psiy2(:,1,:)
            psiy2(:,-1,:)=psiy2(:,0,:)
            psiy2(:,jmax,:)=psiy2(:,jmax-1,:)
            psiy2(:,jmax+1,:)=psiy2(:,jmax,:)
            psiy2(:,jmax+2,:)=psiy2(:,jmax+1,:)
            if (vareps ==0)  then
                psix1 = 0.0
                psix2 = 0.0
                psiy1 = 0.0
                psiy2 = 0.0
            end if
        end subroutine flux_limitor
    end module flux
