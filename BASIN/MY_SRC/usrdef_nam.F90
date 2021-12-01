MODULE usrdef_nam
   !!======================================================================
   !!                       ***  MODULE  usrdef_nam  ***
   !!
   !!                      ===  BASIN configuration  ===
   !!
   !! User defined : set the domain characteristics of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2017-10  (J. Chanut)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_nam   : read user defined namelist and set global domain size
   !!   usr_def_hgr   : initialize the horizontal mesh 
   !!----------------------------------------------------------------------
   USE dom_oce, ONLY: nimpp , njmpp            ! i- & j-indices of the local domain
   USE par_oce, ONLY:        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_nam   ! called by nemogcm.F90

   !                              !!* namusr_def namelist *!!
   REAL(wp), PUBLIC ::   rn_e1_deg      =     1     ! Resolution in degrees of longitude (Mercator grid)
   REAL(wp), PUBLIC ::   rn_phi0        =    20.    ! Latitude of the south frontier (T point) [degrees]
   INTEGER , PUBLIC ::   nn_iglo        =    42     ! Number of grid points along i direction
   INTEGER , PUBLIC ::   nn_jglo        =    62     ! Number of grid points along j direction
   INTEGER , PUBLIC ::   nn_k           =    42     ! Number of grid points along k direction
   INTEGER , PUBLIC ::   nn_forcingtype =     0     ! Surface forcing type
   LOGICAL , PUBLIC ::   ln_diu_cyc     =  .TRUE.   ! Use diurnal cycle for qsr or not
   LOGICAL , PUBLIC ::   ln_ann_cyc     =  .TRUE.   ! Use an annual cycle or not
   REAL(wp), PUBLIC ::   rn_emp_prop    =     1.    ! Proportionality factor to apply on the EMP
   REAL(wp), PUBLIC ::   rn_trp         =   -40.    ! Retroaction term on T*, must be negative [W/m2/K]
   REAL(wp), PUBLIC ::   rn_srp         =     0.0   ! Retroaction term on S*, must be negative []
   LOGICAL , PUBLIC ::   ln_qsr         =  .TRUE.   ! Solar heat or not
   INTEGER , PUBLIC ::   nn_botcase     =     0     ! bottom definition (0:flat, 1:bowl in cosh, 2: bowl in 1-x**4)
   REAL(wp), PUBLIC ::   rn_H           =  4000.    ! Maximum depth of the basin or asymptotical depth of the basin (depending on the basin shape)
   REAL(wp), PUBLIC ::   rn_hborder     =  2000.    ! Depth of the borders of the basin
   REAL(wp), PUBLIC ::   rn_distLam     =     3.    ! Typical length scale of the slope in longitude
   INTEGER , PUBLIC ::   nn_initcase    =     0     ! initial condition case (0=rest+uniform, 1=rest+stratification)
   INTEGER , PUBLIC ::   nn_perio       =     0     ! periodicity of the domain
   LOGICAL , PUBLIC ::   ln_zco_nam     = .FALSE.   ! z               vertical coordinate
   LOGICAL , PUBLIC ::   ln_zps_nam     = .FALSE.   ! z-partial steps vertical coordinate
   LOGICAL , PUBLIC ::   ln_sco_nam     =  .TRUE.   ! s               vertical coordinate
   INTEGER , PUBLIC ::   nn_ztype       =     0     ! type of vertical grid spacing (0: uniform) (z-coordinate or s pure coordinate)
   LOGICAL , PUBLIC ::   ln_equ_flat    = .true.    ! Flattening the topography in the vicinity of the equator (s-coordinate) to limit the error on the pressure gradient 
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_nam.F90 10074 2018-08-28 16:15:49Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_nam( ldtxt, ldnam, cd_cfg, kk_cfg, kpi, kpj, kpk, kperio )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!                    
      !! ** Purpose :   read user defined namelist and define the domain size
      !!
      !! ** Method  :   read in namusr_def containing all the user specific namelist parameter
      !!
      !!                Here BASIN configuration
      !!
      !! ** input   : - namusr_def namelist found in namelist_cfg
      !!----------------------------------------------------------------------
      CHARACTER(len=*), DIMENSION(:), INTENT(out) ::   ldtxt, ldnam    ! stored print information
      CHARACTER(len=*)              , INTENT(out) ::   cd_cfg          ! configuration name
      INTEGER                       , INTENT(out) ::   kk_cfg          ! configuration resolution
      INTEGER                       , INTENT(out) ::   kpi, kpj, kpk   ! global domain sizes 
      INTEGER                       , INTENT(out) ::   kperio          ! lateral global domain b.c.
      !
      INTEGER ::   ios, ii               ! Local integer
      REAL(wp)::   zh                    ! Local scalars
      !!
      NAMELIST/namusr_def/ rn_e1_deg, rn_phi0, nn_iglo, nn_jglo, nn_k    &
         &                 , rn_emp_prop                                 &
         &                 , nn_botcase, nn_initcase, nn_forcingtype     &
         &                 , nn_perio, ln_zco_nam, ln_zps_nam            &
         &                 , ln_sco_nam, nn_ztype, rn_H, rn_hborder      &
         &                 , rn_distLam, ln_equ_flat, ln_ann_cyc         &
         &                 , rn_trp, rn_srp, ln_qsr, ln_diu_cyc
      !!----------------------------------------------------------------------
      !
      ii = 1
      !
      REWIND( numnam_cfg )          ! Namelist namusr_def (exist in namelist_cfg only)
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist', .TRUE. )
      !
      !
      WRITE( ldnam(:), namusr_def )
      !
      cd_cfg = 'BASIN'             ! name & resolution (not used)
      kk_cfg = rn_e1_deg
      !
      ! Global Domain size
      kpi = nn_iglo
      kpj = nn_jglo
      kpk = nn_k
      !
      kperio = nn_perio   ! 0: closed basin, 8: south symmetrical
      IF( kperio == 8 )   rn_phi0 = -rn_e1_deg   ! manually set the longitude of the global first (southern) T point
      !                             ! control print
      WRITE(ldtxt(ii),*) '   '                                                                                                  ;   ii = ii + 1
      WRITE(ldtxt(ii),*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'                           ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '~~~~~~~~~~~ '                                                                                         ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Namelist namusr_def : BASIN test case'                                                             ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Resolution in degrees of longitude (Mercator grid)      rn_e1_deg      = ', rn_e1_deg, 'degrees'   ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Latitude of the south frontier (T point) [degrees]      rn_phi0        = ', rn_phi0, 'degrees'     ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Number of grid points along i direction                 nn_piglo       = ', kpi                    ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Number of grid points along j direction                 nn_pjglo       = ', kpj                    ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Number of grid points along k direction                 nn_k           = ', kpk                    ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Use surface forcing like Gyre (1) or remade one (0)     nn_forcingtype = ', nn_forcingtype         ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Use annual cycle or not                                 ln_ann_cyc     = ', ln_ann_cyc             ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Proportionality factor applied on EMP                   rn_emp_prop    = ', rn_emp_prop            ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Bottom definition                                       nn_botcase     = ', nn_botcase             ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      (0:flat, 1:bowl in cosh, 2: bowl in 1-x**4)'                                                    ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Maximum / asymptotical depth of the basin               rn_H           = ', rn_H, 'm'              ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Depth of the bathymetry at the coast                    rn_hborder     = ', rn_hborder, 'm'        ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Typical length scale of the slope in longitude          rn_distLam     = ', rn_distLam, 'degrees'
      WRITE(ldtxt(ii),*) '   Initial condition case                                  nn_initcase    = ', nn_initcase            ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      (0:rest and constant T and S, '                                                                 ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '       1: rest and stratification)'                                                                   ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Periodicity of the domain                               nn_perio       = ', nn_perio               ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '      (0:closed, 8:south symmetrical)'                                                                ;   ii = ii + 1
      WRITE(ldtxt(ii),*) '   Flattening the topography at the equator                 ln_equ_flat   = ', ln_equ_flat            ;   ii = ii + 1
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
