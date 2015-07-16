
!  -------------------------------------------------------------

    subroutine tick(nvar,cut,nstart,vtime,time,naux,start_time,&
                rest,dt_max)

        use geoclaw_module
        use refinement_module, only: varRefTime
        use amr_module
        use topo_module, only: dt_max_dtopo, num_dtopo, topo_finalized,&
                       aux_finalized, topo0work

        implicit double precision (a-h,o-z)

        logical vtime,dumpout/.false./,dumpchk/.false./,rest,dump_final
        dimension dtnew(maxlv), ntogo(maxlv), tlevel(maxlv)
        integer clock_start, clock_finish, clock_rate

        ! :::::::::::::::::::::::::::: TICK :::::::::::::::::::::::::::::
        !  main driver routine.  controls:
        !        integration  of all grids.
        !        error estimation / regridding
        !        output counting
        !        updating of fine to coarse grids

        !  parameters:
        !     nstop   = # of coarse grid time steps to be taken
        !     iout    = output interval every 'iout' coarse time steps
        !               (if 0, not used - set to inf.)
        !     vtime   = true for variable timestep, calculated each coarse step
        !
        !  integration strategy is to advance a fine grid until it catches
        !  up to the coarse grid. this strategy is applied recursively.
        !  coarse grid goes first.
        !
        !  nsteps: used to count how number steps left for a level to be
        !          integrated before it catches up with the next coarser level.
        !  ncycle: counts number of coarse grid steps = # cycles.
        !
        !  icheck: counts the number of steps (incrementing by 1
        !          each step) to keep track of when that level should
        !          have its error estimated and finer levels should be regridded.
        ! ::::::::::::::::::::::::::::::::::::;::::::::::::::::::::::::::

        ncycle         = nstart
        call setbestsrc()     ! need at very start of run, including restart
        if (iout .eq. 0) then
            !output_style 1 or 2
            iout  = iinfinity
            nextout = 0
            if (nout .gt. 0) then
                nextout = 1
                if (nstart .gt. 0) then
                    ! restart: make sure output times start after restart time
                    do ii = 1, nout
                        if (tout(ii) .gt. time) then
                            nextout = ii
                            exit
                        endif
                    end do
                endif
            endif
        endif
        nextchk = 1
        if ((nstart .gt. 0) .and. (checkpt_style.eq.2)) then
        !if this is a restart, make sure chkpt times start after restart time
            do ii = 1, nchkpt
                if (tchk(ii) .gt. time) then
                    nextchk = ii
                    exit
                endif
            enddo
            continue
        endif

        tlevel(1)      = time

        do  i       = 2, mxnest
            tlevel(i) = tlevel(1)
        enddo

!  ------ start of coarse grid integration loop. ------------------

        if (ncycle .le. nstop .and. time .le. tfinal) then

            if (nout .gt. 0) then
                if (nextout  .le. nout) then
                    outtime       = tout(nextout)
                else
                    outtime       = rinfinity
                endif
            else
                outtime = tfinal
            endif

            if (nextchk  .le. nchkpt) then
                chktime       = tchk(nextchk)
            else
                chktime       = rinfinity
            endif

            dumpout = .false.  !# may be reset below

            if (time.lt.outtime .and. time+1.001*possk(1) .ge. outtime) then
                ! adjust time step  to hit outtime exactly, and make output
                !  apr 2010 mjb: modified so allow slightly larger timestep to
                !  hit output time exactly, instead of taking minuscule timestep
                !  should still be stable since increase dt in only 3rd digit.
                oldposs = possk(1)
                possk(1) = outtime - time
                write(*,*)" old possk is ", possk(1)
                diffdt = oldposs - possk(1)  ! if positive new step is smaller

                if (.false.) then
                    write(SED_PARM_UNIT,*) &  ! notify of change
                        " Adjusting timestep by ",e10.3,&
                        " to hit output time of ",e12.6,diffdt,outtime
                    write(*,*)" new possk is ", possk(1)
                    if (diffdt .lt. 0.) then ! new step is slightly larger
                        pctIncrease = -100.*diffdt/oldposs   ! minus sign to make whole expr. positive
                            write(*,123)" New step is ",e8.2," % larger.",
                           "  Should still be stable", pctIncrease

                    endif
                endif
                do i = 2, mxnest
                    possk(i) = possk(i-1) / kratio(i-1)
                enddo
                if (nout .gt. 0) then
                    nextout = nextout + 1
                    dumpout = .true.
                endif
            endif
            if (time.lt.chktime .and. time + possk(1) .ge. chktime) then
                ! adjust time step  to hit chktime exactly, and do checkpointing
                possk(1) = chktime - time
                do  i = 2, mxnest
                    possk(i) = possk(i-1) / kratio(i-1)
                enddo
                nextchk = nextchk + 1
                dumpchk = .true.
            else
                dumpchk = .false.
            endif
            level        = 1
            ntogo(level) = 1
            do i = 1, maxlv
                dtnew(i)  = rinfinity
            enddo

        !     ------------- regridding  time?  ---------
        !
        ! check if either
        !   (i)  this level should have its error estimated before being advanced
        !   (ii) this level needs to provide boundary values for either of
        !        next 2 finer levels to have their error estimated.
        !        this only affects two grid levels higher, occurs because
        !        previous time step needs boundary vals for giant step.
        !       no error estimation on finest possible grid level
        !

































