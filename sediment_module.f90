! ============================================================================
!  File:        sediment_module.f90
!  Program:     sediment_module
!  Author:      Hui Tang tanghui@vt.edu
!  This module is developed based on previous work from Kyle Mandli and 
!  David George. This module is designed to read sediment data and setup sediment
!  condition to Geoclaw
! ============================================================================
!    Copyright (C) 2010-04-21 Clawpack Developers http://www.clawpack.org
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD)
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================
! ============================================================================
! Read sediment files as specified in sediment.data
! including two files:
! sediment thickness file and grainsize distribution file
! Each sediment thickness file has a type stored in sedtype(i).
!   sedtype = 1:  standard GIS format: 3 columns: lon,lat,thickness(m)
!   sedtype = 2:  Header as in DEM file, thickness(m) one value per line
!   sedtype = 3:  Header as in DEM file, thickness(m) one row per line
! For other formats modify readtopo routine.
! Grainsize distribution file is similar:
!   sedtype = 1:  standard GIS format: 2+gmax columns: lon,lat,P1...Pgmax
!   sedtype = 2:  Header as in DEM file, percentages for one location per line
!   sedtype = 3:  Header as in DEM file, percentages for one location per line
! advancing northwest to northeast then from north to south. Values should
! be uniformly spaced.
!
! Sediment data will update during computatio
! =============================================================================
    module sediment_module

        use Set_Precision
        use amr_module, only: tstart_thisrun
        use topo_mocule, only: rectintegral, intersection

        implicit none
        save

        ! Work array for Sediment for all t
        real(kind=8), allocatable :: sedtwork(:),sedpwork(:)

        ! Sediment file data
        integer :: test_sediment !don't need in this case actually, only for test
        character(len=150), allocatable :: sedtfname(:),sedpfname(:)
        integer :: mtsedfiles,mtsedsize,mpsedsize,mpsedfiles
        real(kind=Prec), allocatable :: xlowsed(:), ylowsed(:), tlowsed(:)
        real(kind=Prec), allocatable :: xhised(:), yhised(:), thised(:)
        real(kind=Prec), allocatable :: dxsed(:), dysed(:)
        real(kind=Prec), allocatable :: sedtime(:)
        integer, allocatable ::  mxsed(:), mysed(:)

        integer, allocatable :: i0sed(:), mtsed(:), msedorder(:)
        integer, allocatable :: minlevelsed(:), maxlevelsed(:), isedtype(:)
        integer, allocatable :: sedID(:),sed0save(:)
        logical :: sed_finalized ! don't need here 

        ! sediment transport support
        integer :: isedtran, aux_finalized

        !Analytic sediment add later

        ! Initial topography
        ! Work array for initial topography (only arrays where topo evolves)
        real(kind=Prec), allocatable :: sed0work(:)
        integer, allocatable :: i0sed0(:),sed0ID(:)
        integer :: msed0size,msed0files


        ! ========================================================================
        !  Constants and parameter
        ! ========================================================================
        Real(kind=Prec) ::      rhos,rho,por,cmax,facDc,hcr,thick
        Integer :: gmax,lmax
        Real(kind=Prec), dimension(:),allocatable :: D
        Real(kind=Prec) ::      g,vis,Te,Trep,eps,k0,m0,tsfac,Tsmin,smax,cf,facDc
        Real(kind=Prec) ::      nuh,nuhfac,gammaWs,a1,hswitch,wetslp,dryslp
        Integer :: sourcesink,avalanching,struct,morfac,sws
        Real(kind=Prec) ::      morstart
        Real(kind=Prec) ::       vareps,split,merge,beta
        CHARACTER(120)  ::      limit_method,trim,method
        contains

            subroutine read_sed_setting(file_name)

                use geoclaw_module
                use amr_module, only: xlower,xupper,xlower,yupper
                use Set_Precision
                use topo_module, only: topoarea

                implicit none

                !Input arguments
                character(len=*), intent(in), optional :: file_name

                !locals
                integer, parameter :: iunit = 7
                integer :: i,j,ised,finer_than,rank
                real(kind=Prec) :: area_i,area_j,x_junk,y_junk
                real(kind=Prec) :: area, area_domain

                open(unit=SED_PARM_UNIT,file='fort.geo',status="unknown",action="write")
                ! Open and begin parameter file output
                write(SED_PARM_UNIT,*) ' '
                write(SED_PARM_UNIT,*) '--------------------------------------------'
                write(SED_PARM_UNIT,*) 'SETSED:'
                write(SED_PARM_UNIT,*) '---------'

                if (present(file_name)) then
                    call opendatafile(iunit, file_name)
                else
                    call opendatafile(iunit, 'sed.data')
                endif

                ! Read in sediment specification type
                read(iunit,"(i1)") test_sediment

                ! Primary sedimenr type, read in sediment files specified
                if (test_sediment == 0) then
                    read(iunit,*) mtsedfiles
                    read(iunit,*) mpsedfiles

                    if (mtsedfiles == 0) then
                        write(SED_PARM_UNIT,*) '   mtsedfiles = 0'
                        write(SED_PARM_UNIT,*) '   No sediment thickness files specified, '
                        write(SED_PARM_UNIT,*) '          will set totalthick(x,y) = 0 in setaux'
                        return
                    endif
                    if (mpsedfiles == 0) then
                        write(SED_PARM_UNIT,*) '   mpsedfiles = 0'
                        write(SED_PARM_UNIT,*) '   No sediment grainsize distribution files specified, '
                        write(SED_PARM_UNIT,*) '          will set pebbed(x,y,n) = 1/gmax in setaux'
                        return
                    endif
                    if (mtsedfiles .ne. ptsedfiles) then
                        write(SED_PARM_UNIT,*) 'Warning:'
                        write(SED_PARM_UNIT,*) '   The number of thickness and grainsize distribution file is not same! '
                        write(SED_PARM_UNIT,*) '   Check your data! '
                        return
                    endif
                    write(SED_PARM_UNIT,*) '   msedfiles = ',msedfiles
                    write(SED_PARM_UNIT,*) '   psedfiles = ',psedfiles

                    ! Read and allocate data parameters for each file, as thickness file and grainsize file is same, we use one to put, two file share the same data may not good
                    allocate(mxsed(mtsedfiles),mysed(mtsedfiles))
                    allocate(xlowsed(mtsedfiles),ylowsed(mtsedfiles))
                    allocate(tlowsed(mtsedfiles),xhised(mtsedfiles),yhised(mtsedfiles))
                    allocate(thised(mtsedfiles),dxsed(mtsedfiles),dysed(mtsedfiles))
                    allocate(sedtfname(mtsedfiles),isedtype(mtsedfiles))
                    allocate(sedpfname(mpsedfiles),ipsedtype(mpsedfiles)) !add to file
                    allocate(minlevelsed(mtsedfiles),maxlevelsed(mtsedfiles))
                    allocate(i0sed(mtsedfiles),msed(mtsedfiles),msedorder(mtsedfiles))
                    allocate(sedID(mtsedfiles),sedtime(mtsedfiles),sed0save(mtsedfiless))
                    allocate(i0sed0(mtsedfiles),sed0ID(mtsedfiles))

                    do i=1,mtsedfiles
                        read(iunit,*) sedtfname(i)
                        read(iunit,*) isedtype(i),minlevelsed(i), maxlevelsed(i), &
                            tlowsed(i),thised(i)

                        write(SED_PARM_UNIT,*) '   '
                        write(SED_PARM_UNIT,*) '   ',sedtfname(i)
                        write(SED_PARM_UNIT,*) '  isedtype = ', isedtype(i)
                        write(SED_PARM_UNIT,*) '  minlevel, maxlevel = ', &
                            minlevelsed(i), maxlevelsed(i)
                        write(SED_PARM_UNIT,*) '  tlow, thi = ', tlowsed(i),thised(i)
                        call read_sed_header(sedtfname(i),isedtype(i),mxsed(i), &
                                mysed(i),xlowsed(i),ylowsed(i),xhised(i),yhised(i), &
                                dxsed(i),dysed(i))
                        sedID(i) = i
                        msed(i) = mxsed(i)*mysed(i)
                    enddo
                    ! Indexing into work array
                    i0topo(1)=1
                    if (mtsedfiles > 1) then
                        do i=2,mtsedfiles
                            i0topo(i)=i0topo(i-1) + mtopo(i-1)
                        enddo
                    endif

                    do i=1,mpsedfiles
                        read(iunit,*) sedpfname(i)
                        read(iunit,*) ipsedtype(i),minlevelsed(i), maxlevelsed(i), &
                                    tlowsed(i),thised(i)
                        write(SED_PARM_UNIT,*) '   '
                        write(SED_PARM_UNIT,*) '   ',sedpfname(i)
                        write(SED_PARM_UNIT,*) '  isedtype = ', isedtype(i)
                        write(SED_PARM_UNIT,*) '  minlevel, maxlevel = ', &
                                minlevelsed(i), maxlevelsed(i)
                        write(SED_PARM_UNIT,*) '  tlow, thi = ', tlowsed(i),thised(i)
                        call read_sed_header(sedpfname(i),isedtype(i),mxsed(i), &
                                mysed(i),xlowsed(i),ylowsed(i),xhised(i),yhised(i), &
                                dxsed(i),dysed(i))
                    enddo

                    ! Read sediment information and allocate space for each file
                    msedsize = sum(msed)
                    allocate(sedtwork(msedsize))
                    allocate(sedpwork(psedsize))
                    do i=1,msedfiles
                        sedID(i) = i
                        sedtime(i) = -huge(1.0)
                        call read_tsed_file(mxsed(i),mysed(i),isedtype(i),sedtfname(i), &
                            sedwork(i0sed(i):i0sed(i)+msed(i)-1))
                    enddo

                    do i=1,psedfiles
                        sedID(i) = i
                        sedtime(i) = -huge(1.0)
                        call read_psed_file(mxsed(i),mysed(i),isedtype(i),sedpfname(i), &
                            sedwork(i0sed(i):i0sed(i)+msed(i)-1))
                    enddo
                    ! Sediment order...This determines which order to process sediment data
                    !
                    ! The finest one will be given priority in any region
                    ! msedorder(rank) = i means that i'th sediment file has rank rank,
                    ! where the file with rank=1 is the finest and considered first.
                    do i=1,mtsedfiles
                        finer_than = 0
                        do j=1,mtsedfiles
                            if (j /= i) then
                                area_i=dxsed(i)*dysed(i)
                                area_j=dxsed(j)*dysed(j)
                                if (area_i < area_j) finer_than = finer_than + 1
                                ! if two files have the same resolution, order is
                                ! arbitrarily chosen
                                if ((area_i == area_j).and.(j < i)) then
                                    finer_than = finer_than + 1
                                endif
                            endif
                        enddo
                        ! ifinerthan tells how many other files, file i is finer than
                        rank = mtsedfiles - finer_than
                        msedorder(rank) = i
                    enddo

                    write(SED_PARM_UNIT,*) ' '
                    write(SED_PARM_UNIT,*) '  Ranking of sediment files', &
                                            '  finest to coarsest: ', &
                                (msedorder(rank),rank=1,mtsedfiles)
                    write(SED_PARM_UNIT,*) ' '

                    msed0size = dot_product(msed,sed0save)
                    allocate(sed0work(sed0size))
                    do i = 2,mtsedfiles
                        i0sed0(i)= i0sed0(i-1) + msed(i-1)*sed0save(i-1)
                    enddo

                    do i = 1,mtsedfiles
                        if (sed0save(i)>0) then
                            sed0work(i0sed0(i):i0sed0(i)+msed(i)-1) = &
                            sedwork(i0sed(i):i0sed(i)+msed(i)-1)
                        endif
                    enddo



                    ! Check that topo arrays cover full domain: need to check the area but not now!
                    !call topoarea(xlower,xupper,ylower,yupper,1,area)
                    !area_domain = (yupper-ylower)*(xupper-xlower)
                    !if (abs(area - area_domain) > 1e-12*area_domain) then
                    !    write(6,*) '**** sediment arrays do not cover domain'
                    !    write(6,*) '**** area of overlap = ', area
                    !    write(6,*) '**** area of domain  = ', area_domain
                    !    stop
                    !endif
                endif
            end subroutine read_sed_setting
            ! ========================================================================
            !  read_tsed_file(mx,my,sed_type,fname,totalthick)
            !
            !  Read sediment thicness file.
            ! ========================================================================

            subroutine read_tsed_file(mx,my,sed_type,fname,totalthick)

                use geoclaw_module
                use Set_Precision

                implicit none

                ! Arguments
                integer, intent(in) :: mx,my,sed_type
                character(len=150), intent(in) :: fname
                real(kind=Prec), intent(inout) :: totalthick(1:mx*my)

                ! Locals
                integer, parameter :: iunit = 19, miss_unit = 17
                real(kind=Prec), parameter :: sed_missing = -150.d0
                logical, parameter :: maketype2 = .false.
                integer :: i,j,num_points,missing,status,sed_start
                real(kind=Prec) :: no_data_value,x,y,z,th_temp

                print *, ' '
                print *, 'Reading sediment thickness file  ', fname

                open(unit=iunit, file=fname, status='unknown',form='formatted')

                select case(abs(sed_type))
                ! ASCII file with x,y,z values on each line.
                ! (progressing from upper left corner across rows, then down)
                ! Assumes a uniform rectangular grid of data values.
                    case(1)
                        i = 0
                        status = 0
                        do i=1,mx*my
                            read(iunit,fmt=*,iostat=status) x,y,th_temp
                            if ((i > mx * my) .and. (status == 0)) then
                                print *,'*** Error: i > mx*my = ',mx*my
                                print *,'*** i, mx, my: ',i,mx,my
                                print *,'*** status = ',status
                                stop
                            endif

                            if (status /= 0) then
                                print *,"Error reading sediment thickness file, reached EOF."
                                print *,"  File = ",fname
                                stop
                            else
                                totalthick(i) = th_temp
                            endif
                        enddo

                    ! ================================================================
                    ! ASCII file with header followed by thickness data
                    ! (progressing from upper left corner across rows, then down)
                    ! one value per line if topo_type=2 or
                    ! mx values per line if topo_type=3
                    ! ================================================================
                    case(2:3)
                    ! Read header
                        do i=1,5
                            read(iunit,*)
                        enddo
                        read(iunit,*) no_data_value

                        ! Read in data
                        missing = 0
                        select case(abs(sed_type))
                            case(2)
                                do i=1,mx*my
                                    read(iunit,*) totalthick(i)
                                    if (totalthick(i) == no_data_value) then
                                        missing = missing + 1
                                        totalthick(i) = sed_missing
                                    endif
                                enddo
                            case(3)
                                do j=1,my
                                    read(iunit,*) (totalthick((j-1)*mx + i),i=1,mx)
                                    do i=1,mx
                                        if (totalthick((j-1)*mx + i) == no_data_value) then
                                            missing = missing + 1
                                            totalthick((j-1)*mx + i) = sed_missing
                                        endif
                                    enddo
                                enddo
                        end select

                        ! Write a warning if we found and missing values
                        if (missing > 0)  then
                            print *, '   WARNING...some missing data values this file'
                            print *, '       ',missing,' missing data values'
                            print *, '              (see fort.missing)'
                            print *, '   These values have arbitrarily been set to ',&
                                    sed_missing
                        endif
                    end select

                    close(unit=iunit)

                endif

            end subroutine read_tsed_file


            ! ========================================================================
            !  read_psed_file(mx,my,sed_type,fname,totalthick)
            !
            !  Read sediment grain size distribution file.
            ! ========================================================================

            subroutine read_psed_file(mx,my,sed_type,fname,pbed)

                use geoclaw_module
                use Set_Precision

                implicit none

                ! Arguments
                integer, intent(in) :: mx,my,sed_type
                character(len=150), intent(in) :: fname
                real(kind=Prec), intent(inout) :: pbed(1:mx*my,gmax)

                ! Locals
                integer, parameter :: iunit = 19, miss_unit = 17
                real(kind=Prec), parameter :: sed_missing = -150.d0
                logical, parameter :: maketype2 = .false.
                integer :: i,j,num_points,missing,status,sed_start
                real(kind=Prec) :: no_data_value,x,y,z
                real(kind=Prec) :: p_temp(gmax)

                print *, ' '
                print *, 'Reading sediment grain size distribution file  ', fname

                open(unit=iunit, file=fname, status='unknown',form='formatted')

                select case(abs(sed_type))
                    ! ASCII file with x,y,z values on each line.
                    ! (progressing from upper left corner across rows, then down)
                    ! Assumes a uniform rectangular grid of data values.
                    case(1)
                        i = 0
                        status = 0
                        do i=1,mx*my
                            read(iunit,fmt=*,iostat=status) x,y,p_temp !may have problem
                            if ((i > mx * my) .and. (status == 0)) then
                                print *,'*** Error: i > mx*my = ',mx*my
                                print *,'*** i, mx, my: ',i,mx,my
                                print *,'*** status = ',status
                                stop
                            endif

                            if (status /= 0) then
                                print *,"Error reading sediment thickness file, reached EOF."
                                print *,"  File = ",fname
                                stop
                            else
                                pbed(i,j,:) = p_temp
                            endif
                        enddo

                        ! ================================================================
                        ! ASCII file with header followed by thickness data
                        ! (progressing from upper left corner across rows, then down)
                        ! one value per line if topo_type=2 or
                        ! mx values per line if topo_type=3
                        ! ================================================================
                        case(2:3)
                            ! Read header
                            do i=1,5
                                read(iunit,*)
                            enddo
                            read(iunit,*) no_data_value

                            ! Read in data
                            missing = 0
                            select case(abs(sed_type))
                                case(2)
                                    do i=1,mx*my
                                        read(iunit,*) pbed(i,:)
                                        do j = 1, gmax
                                            if (pbed(i,j) == no_data_value) then
                                                missing = missing + 1
                                                pbed(i,j) = sed_missing
                                            endif
                                        enddo
                                    enddo
                                case(3)
                                    do j=1,my
                                        read(iunit,*) (pbed((j-1)*mx + i,:),i=1,mx)
                                        do i=1,mx
                                            do j = 1, gmax
                                                if (pbed((j-1)*mx + i,:) == no_data_value) then
                                                    missing = missing + 1
                                                    pbed((j-1)*mx + i,j) = sed_missing
                                                endif
                                            enddo
                                    enddo
                            end select

                            ! Write a warning if we found and missing values
                        if (missing > 0)  then
                            print *, '   WARNING...some missing data values this file'
                            print *, '       ',missing,' missing data values'
                            print *, '              (see fort.missing)'
                            print *, '   These values have arbitrarily been set to ',&
                                            sed_missing
                        endif
                    end select
                    close(unit=iunit)
                endif

            end subroutine read_tsed_file


            ! ========================================================================
            ! subroutine read_sed_header(fname,sed_type,mx,my,xll,yll,xhi,yhi,dx,dy)
            ! ========================================================================
            !  Read sediment file header to determine space needed in allocatable array
            !
            !  :Input:
            !   - fname - (char) Name of file
            !   - sed_type - (int) Type of sediment file (-3 < topo_type < 3)
            !
            !  :Output:
            !   - mx,my - (int) Number of grid points
            !   - xll,yll,xhi,yhi - (float) Lower and upper coordinates for grid
            !   - dx,dy - (float) Spatial resolution of grid
            ! ========================================================================
            subroutine read_sed_header(fname,sed_type,mx,my,xll,yll,xhi,yhi,dx,dy)

                use geoclaw_module

                implicit none

                ! Input and Output
                character(len=150), intent(in) :: fname
                integer, intent(in) :: sed_type
                integer, intent(out) :: mx,my
                real(kind=8), intent(out) :: xll,yll,xhi,yhi,dx,dy

                ! Local
                integer, parameter :: iunit = 19
                integer :: sed_size, status
                real(kind=8) :: x,y,th,nodata_value
                logical :: found_file

                inquire(file=fname,exist=found_file)
                if (.not. found_file) then
                    print *, 'Missing sed file:'
                    print *, '   ', fname
                    stop
                endif

                open(unit=iunit, file=fname, status='unknown',form='formatted')

                select case(abs(sed_type))
                    ! ASCII file with 3 columns
                    ! determine data size
                    case(1)
                    ! Initial size variables
                    sed_size = 0
                    mx = 0

                    ! Read in first values, determines xlow and yhi
                    read(iunit,*) xll,yhi
                    sed_size = sed_size + 1
                    mx = mx + 1

                    ! Go through first row figuring out mx, continue to count
                    y = yhi
                    do while (yhi == y)
                        read(iunit,*) x,y,th
                        sed_size = sed_size + 1
                        mx = mx + 1
                    enddo
                    mx = mx - 1
                    ! Continue to count the rest of the lines
                    do
                        read(iunit,fmt=*,iostat=status) x,y,th
                        if (status /= 0) exit
                        sed_size = sed_size + 1
                    enddo
                    if (status > 0) then
                        print *,"ERROR:  Error reading header of sediment file ",fname
                        stop
                    endif

                    ! Calculate remaining values
                    my = sed_size / mx
                    xhi = x
                    yll = y
                    dx = (xhi-xll) / (mx-1)
                    dy = (yhi-yll) / (my-1)

                    ! ASCII file with header followed by z data
                    case(2:3)
                        read(iunit,*) mx
                        read(iunit,*) my
                        read(iunit,*) xll
                        read(iunit,*) yll
                        read(iunit,*) dx
                        read(iunit,*) nodata_value
                        dy = dx
                        xhi = xll + (mx-1)*dx
                        yhi = yll + (my-1)*dy

                    case default
                        print *, 'ERROR:  Unrecognized sed_type'
                        print *, '    sed_type = ',sed_type
                        print *, '  for sediment file:'
                        print *, '   ', fname
                        stop
                end select

                close(iunit)

                write(GEO_PARM_UNIT,*) '  mx = ',mx,'  x = (',xll,',',xhi,')'
                write(GEO_PARM_UNIT,*) '  my = ',my,'  y = (',yll,',',yhi,')'
                write(GEO_PARM_UNIT,*) '  dx, dy (meters/degrees) = ', dx,dy

            end subroutine read_sed_header
            ! ========================================================================
            !  set_geo(fname)
            !  Reads in user parameters from the given file name if provided
            subroutine set_sed(file_name)

                use geoclaw_module, only : pi

                implicit none

                ! Input
                character(len=*), intent(in), optional :: file_name

                open(unit=SED_PARM_UNIT,file='fort.geo',status="unknown",action="write")

                write(SED_PARM_UNIT,*) ' '
                write(SED_PARM_UNIT,*) '--------------------------------------------'
                write(SED_PARM_UNIT,*) 'Sediment Parameters:'
                write(SED_PARM_UNIT,*) '-------------------'

                ! Read user parameters from setsed.data or sediment.data, see detail in these documents
                if (present(file_name)) then
                    call opendatafile(unit, file_name)
                else
                    call opendatafile(unit, 'set_sed.data')
                endif
                ! Sediment parameters
                read(unit,"(1d16.8)") rhos
                read(unit,"(1d16.8)") rho
                read(unit,"(1d16.8)") por
                read(unit,"(1d16.8)") cmax
                read(unit,"(1d16.8)") facDc
                read(unit,"(1d16.8)") thick
                read(unit,*) gmax
                read(unit,*) lmax
                read(unit,*) hcr
                if (gmax) then
                    allocate(D(gmax))
                    read(unit,*) D(:)
                endif
                !physical parameter
                read(unit,*) g
                read(unit,*) vis
                read(unit,*) Te
                read(unit,*) Trep
                read(unit,*) eps
                read(unit,*) k0
                read(unit,*) m0
                read(unit,*) tsfac
                read(unit,*) Tsmin
                read(unit,*) smax
                read(unit,*) cf
                read(unit,*) nuh
                read(unit,*) nuhfac
                read(unit,*) hswitch
                read(unit,*) wetslp
                read(unit,*) dryslp
                read(unit,*) toler
                !processes control paramter 
                read(unit,*) sourcesink
                read(unit,*) avalanching
                read(unit,*) struct
                read(unit,*) morfac
                read(unit,*) thetanum
                read(unit,*) sws
                read(unit,*) morstart
                read(unit,*) split
                read(unit,*) merge
                read(unit,*) beta
                !algorithm parameter
                read(unit,*) vareps
                read(unit,*) k1
                !method
                read(unit,*) limit_method
                read(unit,*) trim
                read(unit,*) method
                close(unit)
                if (gmax<1) then
                    print *, 'ERROR in setsed:  must set up grain size for sediment'
                    stop
                endif
                write(SED_PARM_UNIT,*) '   Sediment density:',rhos
                write(SED_PARM_UNIT,*) '   Water density:',rho
                write(SED_PARM_UNIT,*) '   Porosity:',por
                write(SED_PARM_UNIT,*) '   Maximum allowed sediment concentration:',cmax
                write(SED_PARM_UNIT,*) '   Control sediment diffusion coefficient:',facDc
                write(SED_PARM_UNIT,*) '   Sediment thickness for each layer:' thick
                write(SED_PARM_UNIT,*) '   Number of Grain size classes:'gmax
                write(SED_PARM_UNIT,*) '   Number of sediment layers:'lmax
                write(SED_PARM_UNIT,*) '   Water depth consider sediment transport:'hrc
                if (friction_forcing) then
                    write(SEO_PARM_UNIT,*) '   Sediment grain size classes:', D(:)
                endif
                write(SED_PARM_UNIT,*) '   gravity:',g
                write(SEO_PARM_UNIT,*) '   kinematic viscosity:',vis
                write(SEO_PARM_UNIT,*) '   Water temperature:',Te
                write(SEO_PARM_UNIT,*) '   Representative wave period:',Trep
                write(SEO_PARM_UNIT,*) '   Threshold water depth above which cells are considered wet:', eps
                write(SEO_PARM_UNIT,*) '   von kaman coefficient:',k0
                write(SEO_PARM_UNIT,*) '   mining coeffcient:',m0
                write(SEO_PARM_UNIT,*) '   Coefficient determining Ts = tsfac * h/ws in sediment source term:',tsfac
                write(SED_PARM_UNIT,*) '   Minimum adaptation time scale in advection diffusion equation sediment:', Tsmin
                write(SED_PARM_UNIT,*) '   maximum Shields parameter for ceq Diane Foster:',smax
                write(SED_PARM_UNIT,*) '   Friction coefficient flow:',cf
                write(SED_PARM_UNIT,*) '   horizontal background viscosity:',nuh
                write(SED_PARM_UNIT,*) '   Viscosity switch for roller induced turbulent horizontal viscosity:',nuhfac
                write(SED_PARM_UNIT,*) '   Water depth at which is switched from wetslp to dryslp',hswitch
                write(SED_PARM_UNIT,*) '   Critical avalanching slope under water:',wetslp
                write(SED_PARM_UNIT,*) '   Critical avalanching slope above water:',dryslp
                write(SED_PARM_UNIT,*) '   toler for sediment flux limitor:',toler
                write(SED_PARM_UNIT,*) '    use source-sink terms to calculate bed level change:',sourcesink
                write(SED_PARM_UNIT,*) '    Include avalanching:',avalanching
                write(SED_PARM_UNIT,*) '    Switch for hard structures:',struct
                write(SED_PARM_UNIT,*) '    morphological acceleration factor:',morfac
                write(SED_PARM_UNIT,*) '    Coefficient determining scheme:',thetanum
                write(SED_PARM_UNIT,*) '    short wave & roller stirring and undertow:',sws
                write(SED_PARM_UNIT,*) '    Start time:',morstart
                write(SED_PARM_UNIT,*) '    Split threshold for variable sediment layer:',split
                write(SED_PARM_UNIT,*) '    Mergethreshold for variable sediment layer ',merge
                write(SED_PARM_UNIT,*) '    Coefficient determining order of accuracy:',vareps
                write(SED_PARM_UNIT,*) '    Coefficient determining whether fully upwind',k1
                write(SED_PARM_UNIT,*) '    method to caculate sediment flux limiter',limit_method
                write(SED_PARM_UNIT,*) '    method to caculate equibrium sediment concentration',trim
                write(SED_PARM_UNIT,*) '    method to caculate sediment flux',method
            end subroutine set_sed
    end module sediment_module



























