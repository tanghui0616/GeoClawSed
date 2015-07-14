
    !*****************************************************************
    !*
    !* This program  is used for set up precision
    !*
    !*****************************************************************
    module Set_Precision !Sets the default precision for floating point numbers

        implicit none

        public :: sngl, dbl, quad, Prec

        integer, parameter :: Sngl = Selected_Real_kind( 6, 37)
        integer, parameter :: Quad = Selected_Real_kind(33, 4931)
        integer, parameter :: Dbl  = Selected_Real_Kind(15, 307)!double precision
        integer, parameter :: Prec = Dbl ! Set default precision to double precision
        integer, parameter :: Charl = 30

    end module Set_Precision

    module params

        use Set_Precision
    ! These parameters are constants, variables read from params.txt, or are scalars derived directly from read input
    !
    !  Rules are:
    !  - [Section] markers indicate new sets of parameters, followed by the name of the set
    !  - Fortran declaration always as "kind  ::  name   = initial value"
    !  - After the declaration of a variable add exclamation mark, followed by [unit] and decription.
    !     If the parameter is essentially only for advanced users, follow the unit declation by "(advanced)"
    !    If the parameter is deprecated, but still used for backwards-compatibility, follow the unit declation by "(deprecated)"
    !  - Description of a variable may continue on a new line as long as the first character is "!" and the
    !    position of the first "!" is greater than 50 characters from the start of the line. Best practice
    !    is to keep in line with the start of the description on the line above.
    !  - Please keep the declaration of "globalvars", "meanvars" and "pointvars" on the line directly after their respective related size
    !    parameter declaration (i.e. declaration of "nglobalvars", "nmeanvar" and "npointvar"). This is needed for the autogeneration of
    !    parameters.inc and subsequent params.dat file.
    !  - To enable parsing of "all_input" subroutine please only use "if () then, ... , endif/else/elseif" blocks,
    !    rather than one line "if () ..." commands in the "all_input" subroutine.
    !
    ! constant
    !  Type                     name                   initialize    !  [unit] (advanced/deprecated) description
        Integer         ::      imax                   =10            !  [-]      Number of points in the x-direction
        Integer         ::      jmax                   =10            !  [-]      Number of points in the y-direction
        Integer         ::      lmax                   =10            !  [-]      Number of points in the y-direction
        Integer         ::      gmax                   =2            !  [-]      Number of Grain size classes
        Integer         ::      i                                    !  [-]      index
        Integer         ::      j                                    !  [-]      index
        Integer         ::      k                                    !  [-]      index
        Integer         ::      nl                                   !  [-]      index
        Integer         ::      ii                                   !  [-]      index
        Integer         ::      suc                    = 0           !  [-]     Calibration factor for suspensions transports
        Integer         ::      bed                    = 0           !  [-]     Calibration factor for bed transports
        Integer         ::      sourcesink             = 0           !  [-]     (advanced) In suspended transport use source-sink terms to calculate bed level change (1) or sus transport gradients (0)
        Integer         ::      avalanching            = 1           !  [-]     Include avalanching (1) or exclude (0)
        Integer         ::      struct                 = 1           !  [-]     Switch for hard structures (1) or close (0)
        Integer         ::      ndz                    = 0           !  [-]    the total number of sediment layer erode during avanlanching
        Integer         ::      nd_var                 = 5           !  [-]     (advanced) Index of layer with variable thickness

        Logical         ::      aval                   = .false.     !  [-]
        Real(kind=Prec) ::      rhos                   = 2650.0      !  [kgm^-3] Solid sediment density (no pores)
        Real(kind=Prec) ::      rho                    = 1025.0      !  [kgm^-3] Density of water
        Real(kind=Prec) ::      g                      = 9.81        !  [ms^-2] Gravitational acceleration
        Real(kind=Prec) ::      vis                    = 8.9e-4      !  [m^2/s] kinematic viscosity
        Real(kind=Prec) ::      delta                  = 0           !  [-]     Parameter for settling velocity
        Real(kind=Prec) ::      Sster                  = 0           !  [-]     Parameter for settling velocity
        Real(kind=Prec) ::      wster                  = 0           !  [-]     Parameter for settling velocity
        Real(kind=Prec) ::      c1                     = 0           !  [-]     Parameter for settling velocity
        Real(kind=Prec) ::      c2                     = 0           !  [-]     Parameter for settling velocity
        Real(kind=Prec) ::      Trep                   = 10.0           !  [s]     Representative wave period
        Real(kind=Prec) ::      Te                     = 20.0        !  [C]     Water temperature
        Real(kind=Prec) ::      px                     = 3.1415926   !  [-]     Pi
        Real(kind=Prec) ::      tsfac                  = 0.1         !  [-]     Coefficient determining Ts = tsfac * h/ws in sediment source term
        Real(kind=Prec) ::      Tsmin                  = 0.5         !  [s]     Minimum adaptation time scale in advection diffusion equation sediment
        Real(kind=Prec) ::      Ass                    = 0.0         !  [-]     suspended load coeffient
        Real(kind=Prec) ::      smax                   = -1.0        !  [-]     maximum Shields parameter for ceq Diane Foster
        Real(kind=Prec) ::      cf                     = 0.003         !  [-]     Friction coefficient flow
        Real(kind=Prec) ::      eps                    = 0.005       !  [m]     Threshold water depth above which cells are considered wet
        Real(kind=Prec) ::      cmax                   = 0.1         !  [m^3/m^3]  Maximum allowed sediment concentration
        Real(kind=Prec) ::      perc                   = 0.0         !  [-]     Percentage of bed sediment concentration
        Real(kind=Prec) ::      exp_ero                = 0.0         !  [-]     predict erosion rate per fraction 
        Real(kind=Prec) ::      morfac                 = 1.0         !  [-]     morphological acceleration factor
        Real(kind=Prec) ::      por                    = 0.4         !  [-]     porosity
        Real(kind=Prec) ::      facDc                  = 1.0         !  [-]     control sediment diffusion coefficient
        Real(kind=Prec) ::      nuh                    = 0.1         !  [m^2s^-1]     horizontal background viscosity
        Real(kind=Prec) ::      nuhfac                 = 1.0         !  [m^2s^-1]     Viscosity switch for roller induced turbulent horizontal viscosity
        Real(kind=Prec) ::      thetanum               = 1.0         !  [-]     Coefficient determining whether upwind (1) or central scheme (0.5) is used.
        Real(kind=Prec) ::      facsl                  = 0.0         !  [-]     Factor bedslope effect
        Real(kind=Prec) ::      vareps                 = 1.0         !  [-]     Coefficient determining order of accuracy
        Real(kind=Prec) ::      k0                     = 0.41         !  [-]     von kaman coefficient
        Real(kind=Prec) ::      k1                     = -1.0        !  [-]     Coefficient determining whether fully upwind (-1) or central scheme (0.5) is used.
        Real(kind=Prec) ::      hcr                    = 0.01        !  [m]     water depth consider sediment transport
        Real(kind=Prec) ::      gammaWs                = 0.056       !  [-]     constant to calculate bed roughness
        Real(kind=Prec) ::      a1                     = 0.0         !  [-]     constant for caculating bed roughness
        Real(kind=Prec) ::      m0                     = 0.015         !  [-]     mining coeffcient
        Real(kind=Prec) ::      sws                    = 0.0         !  [-]     (advanced) 1 = short wave & roller stirring and undertow, 0 = no short wave & roller stirring and undertow
        Real(kind=Prec) ::      morstart               = 120.0       !  [s]     Start time morphology, in morphological time
        Real(kind=Prec) ::      hswitch                = 0.1         !  [m]      Water depth at which is switched from wetslp to dryslp
        Real(kind=Prec) ::      dzmax                  = 0.0         !  [m/s/m] Maximum bedlevel change due to avalanching
        Real(kind=Prec) ::      wetslp                 = 0.03         !  [-] Critical avalanching slope under water (dz/dx and dz/dy)
        Real(kind=Prec) ::      dryslp                 = 0.01         !  [-] Critical avalanching slope above water (dz/dx and dz/dy)
        Real(kind=Prec) ::      dx                     = 1.0         !  [m]  cell size of x direction
        Real(kind=Prec) ::      dy                     = 1.0         !  [m]  cell size of y direction
        Real(kind=Prec) ::      split                  = 1.01         !  [-]  Split threshold for variable sediment layer (ratio to nominal thickness)
        Real(kind=Prec) ::      merge                  = 0.01        !  [-]  Merge threshold for variable sediment layer (ratio to nominal thickness)
        Real(kind=Prec) ::      toler                  = 1e-6        !  [-]  toler for sediment flux limitor
        Real(kind=Prec) ::      beta                   = 1           !  [-]  parameter for beta sediment flux limitor
        Real(kind=Prec) ::      dzleft                 = 0.0         !  [m]  thickness change due avanlanching
        Real(kind=Prec) ::      dzt                    = 0.0         !  [m]  sediment thickness move for each step during avanlanching
        Real(kind=Prec) ::      dzavt                  = 0.0         !  [m]  sediment thickness already move during avanlanching
        Real(kind=Prec) ::      frac_dz                = 0.7         !  [m]  (advanced) Relative thickness to split time step for bed
        !Real(kind=Prec) ::      fc                     = 0.0         !  [-]  factor over mass change in cells to limit erosion and deposition
        !Real(kind=Prec) ::      dzbt                   = 0.0         !  [-]  temperary for dzb
        Real(kind=Prec) ::      dt                     = 0.1         !  [m]  time step
        Real(kind=Prec) ::       t                     = 300.0         !  [m]  simulation time
        Real(kind=Prec) ::      Savailable             = 1.0         !  [m]  sediment available to erode or deposit
        Real(kind=Prec) ::      thick                  = 0.05         !  [m]  toal sediment thickness for each layer

        CHARACTER(120)       ::      limit_method           = 'VanAlbada'                  !  [-]     method to caculate sediment flux.
        CHARACTER(120)       ::      trim                   = 'soulsby_vanrijn'          !  [-]     method to caculate equibrium sediment concentration.
        CHARACTER(120)       ::      method                 = 'SVL'                      !  [-]     method to caculate sediment flux.






    !variable
    !  Type                                             name                   initialize    !  [unit] (advanced/deprecated) description        size
        Integer,dimension(:,:),allocatable ::    totalnum             != 0          !  [m] total sediment erodible                     imax:jmax
        Integer,dimension(:,:,:),allocatable ::  indSus                         != 0          !  [-] index of suspended load in x direction     imax:jmax:gmax
        Integer,dimension(:,:,:),allocatable ::  indSub                         != 0          !  [-] index of bed load in x directio            imax:jmax:gmax
        Integer,dimension(:,:,:),allocatable ::  indSvs                         != 0          !  [-] index of suspended load in y direction     imax:jmax:gmax
        Integer,dimension(:,:,:),allocatable ::  indSvb                         != 0          !  [-] index of bed load in y direction           imax:jmax:gmax

        Real(kind=Prec),dimension(:),allocatable ::      D                      != 0          !  [m]         Grain size classes                  gmax
        Real(kind=Prec),dimension(:),allocatable ::      A                      != 0          !  [-]     Parameter for settling velocity         gmax
        Real(kind=Prec),dimension(:),allocatable ::      alpha1                 != 0          !  [-]     Parameter for settling velocity         gmax
        Real(kind=Prec),dimension(:),allocatable ::      alpha2                 != 0          !  [-]     Parameter for settling velocity         gmax
        Real(kind=Prec),dimension(:),allocatable ::      wst                    != 0          !  [-]     Parameter for settling velocity         gmax
        Real(kind=Prec),dimension(:),allocatable ::      alpha                  != 0          !  [-]     Parameter for settling velocity         gmax
        Real(kind=Prec),dimension(:),allocatable ::      R                      != 0          !  [-]     Parameter for settling velocity         gmax
        Real(kind=Prec),dimension(:),allocatable ::      w                      != 0          !  [-]     Parameter for settling velocity         gmax
        Real(kind=Prec),dimension(:),allocatable ::      dster                  != 0          !  [-]     Parameter for settling velocity         gmax
        Real(kind=Prec),dimension(:),allocatable ::      sedcal                 != 0          !  [-] Sediment transport calibration coefficient per grain         gmax
        Real(kind=Prec),dimension(:),allocatable ::      indx                   != 0          !  [-]  indicate the water runup jmax
        Real(kind=Prec),dimension(:),allocatable ::      dz                     != 0          !  [-]     bed gradient for one location           lmax
        Real(kind=Prec),dimension(:,:),allocatable ::    h                      != 0          !  [m]          water depth                        imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::    zs                     != 0          !  [m]          water level                        imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::    z0bed                  != 0          !  [m]          hard structure location                        imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::    dzav                   != 0          !  [m]  total bed level change due to avalanching  imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::    Dmm                    != 0          !  [m]     mean grain size for top layer           imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::    hold                   != 0          !  [m]     water depth at previous time step       imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   hloc                    != 0          !  [m]         modified water depth                imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   Cd                      != 0          !  [-]          Drag Coefficient                   imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   Asb                     != 0          !  [-]          bed load coeffient                 imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   vmg                     != 0          !  [m/s]       velocity maginitude                 imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   vmag2                   != 0          !  [m/s]       velocity maginitude                 imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   term1                   != 0          !  [-]       term1 for sediment transport          imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   term2                   != 0          !  [-]       term2 for sediment transport          imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   term3                   != 0          !  [-]       term3 for sediment transport          imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   wet                     != 0          !  [-]       mark wet/dry points                   imax:jmax
        !Real(kind=Prec),dimension(:,:),allocatable ::   dt                      != 0          !  [s]       computational time step               imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   Dc                      != 0          !  [-]       Diffusion Coefficient                 imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   u                       != 0          !  [m/s]           u velocity                      imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   v                       != 0          !  [m/s]           v velocity                      imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   a2                      != 0          !  [-]constant for caculating bed roughness        imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   urms                    != 0          !  [m/s]   root mean square velocity           imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   urms2                   != 0          !  [m^2/s^2]   root mean square velocity           imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   ustarc                  != 0          !  [m/s]        current shear velocity             imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   taub                    != 0          !  [m/s]        current shear stress               imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   ustarcrit               != 0          !  [m/s]        critical shear velocity            imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   taucrit                 != 0          !  [m/s]        critical shear stress              imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   Tstar                   != 0          !  [-]        parameter for bed roughness          imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   delb                    != 0          !  [m]    average saltation height                 imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   zos                     != 0          !  [m]   The roughness from saltating sediment     imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::   z0                      != 0          !  [m]   The bed roughness                         imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::    structdepth            != 0          !  [m] sediment available to erode                 imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::    dzbdt                  != 0          !  [m/s] bed elevation change rate                 imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::    sedero                 != 0          !  [m] total deposition or erosion                 imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::    totalthick             != 0          !  [m] total sediment erodible                     imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::    dzbdx                  != 0          !  [-]bed level gradient in x direction            imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::    dzbdy                  != 0          !  [-]bed level gradient in y direction            imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::    zon                    != 0          !  [m]     Nikuradse bed roughness                 imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::    zb                     != 0          !  [m]     bed elevation                           imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::    pb                     != 0          !  [m]  grain size distribution for each layer     gmax:lmax
        Real(kind=Prec),dimension(:,:),allocatable ::  dcsdx                    != 0          !  [-]sediment concentration gradient in x direction             imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::  dcsdy                    != 0          !  [-]sediment concentration gradient in x direction             imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::  dcbdx                    != 0          !  [-]sediment concentration gradient in x direction             imax:jmax
        Real(kind=Prec),dimension(:,:),allocatable ::  dcbdy                    != 0          !  [-]sediment concentration gradient in x direction             imax:jmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  frc                    != 0          !  [-]fraction for each Grain size classes on top layer   imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  ws                     != 0          !  [m/s]         settling velocity                 imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  C                      != 0          !  [m^3/m^3]   Sediment concentration              imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  ub_cr                  != 0          !  [m/s]   crtical velocity for bedload            imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  us_cr1                 != 0          !  [m/s] crtical velocity for beginning suspended load            imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  us_cr2                 != 0          !  [m/s] crtical velocity for fully suspended bedload            imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Ts                     != 0          !  [s]   parameter for adaption time               imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Tsg                    != 0          !  [s]            adaption time                    imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Ceq                    != 0          !  [m^3/m^3]total eqiubrium sediment concentration imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Ceqb                   != 0          !  [m^3/m^3]bed eqiubrium sediment concentration   imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Ceqs                   != 0          !  [m^3/m^3]suspended eqiubrium sediment concentration   imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Ceqbg                  != 0          !  [m^3/m^3]final bed eqiubrium sediment concentration   imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Ceqsg                  != 0          !  [m^3/m^3]final suspended eqiubrium sediment concentration   imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  cc                     != 0          !  [m^3/m^3]temperary suspended sediment concentration   imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  ccg                    != 0          !  [m^3/m^3]suspended sediment concentration             imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  ccb                    != 0          !  [m^3/m^3]temperary bed sediment concentration   imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  ccbg                   != 0          !  [m^3/m^3]bed sediment concentration             imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  cu                     != 0          !  [m^3/m^3]suspended sediment concentration caused by u velocity            imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  cub                    != 0          !  [m^3/m^3]bed sediment concentration caused by u velocity            imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  cv                     != 0          !  [m^3/m^3]suspended sediment concentration caused by v velocity            imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  cvb                    != 0          !  [m^3/m^3]bed sediment concentration caused by v velocity            imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Sus                    != 0          !  [m^3/m^3]suspended sediment load in x direction                     imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Sub                    != 0          !  [m^3/m^3]bed sediment load in x direction                     imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Svs                    != 0          !  [m^3/m^3]suspended sediment load in y direction                      imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Svb                    != 0          !  [m^3/m^3]bed sediment load in y direction                     imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  fac                    != 0          !  [-]sediment fraction in water imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  fre                    != 0          !  [-]sediment erosion recution factor imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  dzbed                  != 0          !  [-]bed level gradient                                 imax:jmax:lmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  pbbedu                 != 0          !  [-]Sediment grain class fraction in x direction sediment transport              imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  pbbedv                 != 0          !  [-]Sediment grain class fraction in y direction sediment transport              imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  ero1                   != 0          !  [-]bed erosion rate from susepneded sediment                       imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  ero2                   != 0          !  [-]bed erosion rate from bed sediment                       imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  depo_ex1               != 0          !  [-] exmplicit bed deposition rate from suspended sediment                       imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  depo_ex2               != 0          !  [-] exmplicit bed deposition rate from bed sediment                       imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Susg                   != 0          !  [m^3/m^3]suspended sediment transport in x direction                      imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Subg                   != 0          !  [m^3/m^3]bed sediment transport in x direction                      imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Svsg                   != 0          !  [m^3/m^3]suspended sediment load transport in y direction                     imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Svbg                   != 0          !  [m^3/m^3]bed sediment load transport in y direction                     imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Sout                   != 0          !  [m^3/m^3]total sediment flux                     imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  V_new                  != 0          !  [-] temporary variable for right and left state including u,v,h                     imax:jmax:3
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Vx_l                   != 0          !  [-] temporary variable for left state including u,v,h in x direction                  imax:jmax:3
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Vx_r                   != 0          !  [-] temporary variable for right state including u,v,h in x direction                  imax:jmax:3
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Vy_l                   != 0          !  [-] temporary variable for left state including u,v,h in y direction                  imax:jmax:3
        Real(kind=Prec),dimension(:,:,:),allocatable ::  Vy_r                   != 0          !  [-] temporary variable for right state including u,v,h in y direction                  imax:jmax:3
        Real(kind=Prec),dimension(:,:,:),allocatable ::  psix1                  != 0          !  [-] sediment flux limitor in x direction                        imax:jmax:gmax+3
        Real(kind=Prec),dimension(:,:,:),allocatable ::  psix2                  != 0          !  [-] sediment flux limitor in x direction                        imax:jmax:gmax+3
        Real(kind=Prec),dimension(:,:,:),allocatable ::  psiy1                  != 0          !  [-] sediment flux limitor in y direction                        imax:jmax:gmax+3
        Real(kind=Prec),dimension(:,:,:),allocatable ::  psiy2                  != 0          !  [-] sediment flux limitor in y direction                        imax:jmax:gmax+3
        Real(kind=Prec),dimension(:,:,:),allocatable ::  dzg                    != 0          !  [m] bed elevation change for each grain size class                        imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  laythick               != 0          !  [m] sediment thickness for each layer                       imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:),allocatable ::  edg                    != 0          !  [m] bed elevation change for each grain size class                        imax:jmax:gmax
        Real(kind=Prec),dimension(:,:,:,:),allocatable ::  Cx_l                   != 0          !  [-] temporary variable for left state including sediment concentration in x direction                    imax:jmax:gmax:2
        Real(kind=Prec),dimension(:,:,:,:),allocatable ::  Cx_r                   != 0          !  [-] temporary variable for right state including sediment concentration in x direction                    imax:jmax:gmax:2
        Real(kind=Prec),dimension(:,:,:,:),allocatable ::  Cy_l                   != 0          !  [-] temporary variable for left state including sediment concentration in y direction                    imax:jmax:gmax:2
        Real(kind=Prec),dimension(:,:,:,:),allocatable ::  Cy_r                   != 0          !  [-] temporary variable for right state including sediment concentration in y direction                    imax:jmax:gmax:2
        Real(kind=Prec),dimension(:,:,:,:),allocatable ::  C_new_x                != 0          !  [-] temporary variable for right and left state including sediment concentration in x direction                    imax:jmax:gmax:2
        Real(kind=Prec),dimension(:,:,:,:),allocatable ::  C_new_y                != 0          !  [-] temporary variable for right and left state including sediment concentration in y direction                    imax:jmax:gmax:2
        Real(kind=Prec),dimension(:,:,:,:),allocatable ::  pbbed                  != 0          !  [-]Sediment grain class fraction in first sediment layer              imax:jmax:gmax:lmax

    end module params


















