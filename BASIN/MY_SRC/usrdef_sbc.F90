MODULE usrdef_sbc
   !!======================================================================
   !!                       ***  MODULE  usrdef_sbc  ***
   !! 
   !!                      ===  BASIN configuration  ===
   !!
   !! User defined :   surface forcing of a user configuration
   !!======================================================================
   !! History :  4.0   ! 2017-11  (J.Chanut)  user defined interface
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_sbc    : user defined surface bounday conditions
   !!----------------------------------------------------------------------
   USE oce, ONLY : tsb                                             ! ocean dynamics and tracers
   USE dom_oce, ONLY:                                              ! ocean space and time domain
   USE sbc_oce, ONLY: utau, vtau, taum, wndm, emp, sfx, qns, qsr   ! Surface boundary condition: ocean fields
   USE phycst                                                      ! physical constants
   USE sbcdcy, ONLY: sbc_dcy, nday_qsr                             ! surface boundary condition: diurnal cycle
   !
   USE usrdef_nam, ONLY : nn_forcingtype, rn_emp_prop, ln_ann_cyc, ln_diu_cyc, rn_trp, rn_srp, ln_qsr
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! 
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined) 

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usrdef_sbc_oce      ! routine called in sbcmod module
   PUBLIC   usrdef_sbc_ice_tau  ! routine called by icestp.F90 for ice dynamics
   PUBLIC   usrdef_sbc_ice_flx  ! routine called by icestp.F90 for ice thermo

   REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   ztstar, zqsr_dayMean, zsstar   !: ztstar used in the heat forcing, zqsr_dayMean is the dayly averaged solar heat flux, zsstar for sfx

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_sbc.F90 10074 2018-08-28 16:15:49Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

  SUBROUTINE usrdef_sbc_oce( kt )
     !!---------------------------------------------------------------------
     !!                    ***  ROUTINE usrdef_sbc_oce  ***
     !!              
     !! ** Purpose :   provide at each time-step the GYRE surface boundary
     !!              condition, i.e. the momentum, heat and freshwater fluxes.
     !!
     !! ** Method  :   analytical configuration.
     !!                CAUTION : never mask the surface stress field !
     !!
     !! ** Action  : - set the ocean surface boundary condition, i.e.   
     !!                   utau, vtau, taum, wndm, qns, qsr, emp, sfx
     !!---------------------------------------------------------------------
     INTEGER, INTENT(in) ::   kt      ! ocean time step
     !
     INTEGER  ::   ji, jj
     REAL(wp) ::   zL, z1_2L 
     REAL(wp) ::   zdeltasst, zts_eq, zts_n, zdeltaT
     REAL(wp) ::   zdeltaemp, zemp_mean, zF0, zconv, zaemp, zb, zdeltaemp2, zaemp2
     REAL(wp) ::   zdeltatau, ztau0, zatau
     REAL(wp) ::   zf1, zf2, zf3, zweight2, zweight3
     REAL(wp) ::   za1, za2, za3
     REAL(wp) ::   zphi01, zphi02, zphi03
     REAL(wp) ::   z1_d1, zd2, z1_d2, zd3, z1_d3
     REAL(wp) ::   z1_s
     REAL(wp), DIMENSION(jpi,jpj) ::   zSurf
     REAL(wp) ::   zcos_sais1, zcos_sais2
     !!---------------------------------------------------------------------
     !
     SELECT CASE( nn_forcingtype )
     CASE(0)
        CALL ctl_stop( 'usrdef_sbc_oce : option not available anymore' )
     CASE(1)
        CALL ctl_stop( 'usrdef_sbc_oce : option not available anymore' )
     CASE(2)
        ! just a zonal wind stress, going eastward
        IF( kt .EQ. nit000 .AND. lwp ) THEN
           WRITE(numout,*)'usrdef_sbc_oce : analytical surface fluxes'
           WRITE(numout,*) 'Using just a zonal wind stress, goind eastward'
        ENDIF
        utau(:,:) = -0.1_wp
        vtau(:,:) = 0._wp
        taum(:,:) = 0._wp
        wndm(:,:) = 0._wp
        !
        emp (:,:) = 0._wp
        sfx (:,:) = 0._wp
        qns (:,:) = 0._wp
        qsr (:,:) = 0._wp
     CASE(3)
        ! no forcing
        IF( kt .EQ. nit000 .AND. lwp ) THEN
           WRITE(numout,*)'usrdef_sbc_oce : analytical surface fluxes'
           WRITE(numout,*) 'no forcing'
        ENDIF
        utau(:,:) = 0._wp
        vtau(:,:) = 0._wp
        taum(:,:) = 0._wp
        wndm(:,:) = 0._wp
        !
        emp (:,:) = 0._wp
        sfx (:,:) = 0._wp
        qns (:,:) = 0._wp
        qsr (:,:) = 0._wp
     CASE(4)
        ! Forcing close to Wolfe and Cessi 2014, JPO
        IF( kt == nit000 .AND. lwp ) THEN
           WRITE(numout,*) 'usrdef_sbc_oce : analytical surface fluxes'
           WRITE(numout,*) '~~~~~~~~~~~~~~'
           WRITE(numout,*) '   Forcing close to Wolfe and Cessi 2014, JPO'
           WRITE(numout,*) '      Zonal annual wind'
           WRITE(numout,*) '      Zonal annual E-P'
           IF( ln_qsr )   THEN
              WRITE(numout,*) '      Solar heat flux'
              IF( ln_ann_cyc )   THEN
                 WRITE(numout,*) '         Annual cycle of solar heat flux'
              ELSE
                 WRITE(numout,*) '         Annual average of solar heat flux (December = June)'
              ENDIF
              IF( ln_diu_cyc )   THEN
                 WRITE(numout,*) '         Diurnal cycle of solar heat flux'
              ELSE
                 WRITE(numout,*) '         Daily average of solar heat flux (day = night)'
              ENDIF
           ELSE
              WRITE(numout,*) '      No solar heat flux'
           ENDIF
           WRITE(numout,*) '      Zonal annual T* (heat forcing proportional to (SST - T*)'
           IF( rn_srp /= 0._wp )   THEN
              WRITE(numout,*) '      Zonal annual S* (salt flux proportional to (SSS - S*)'
           ELSE
              WRITE(numout,*) '      No salt restoring'
           ENDIF
        ENDIF
        ! Initialization of parameters
        zL = 61                        ! [degrees] Approximative meridional extend of the basin
        !
        ztau0         =   0.1_wp       ! [Pa]
        zatau         =   0.8_wp       ! [no unit]
        zdeltatau     =   5.77_wp      ! [degrees]
        !ztrp          = -40._wp        ! [W/m2/K] retroaction term on heat fluxes 
        zts_eq        =  25._wp        ! [deg C]
        zts_n         =   0._wp        ! [deg C]
        zdeltasst     =  16.22_wp      ! [degrees]
        !
        zconv         =   3.16e-5_wp   ! convertion factor: 1 m/yr => 3.16e-5 mm/s
        zaemp      =   2._wp           ! [no unit]
        zdeltaemp  =   8.11_wp         ! [degrees]
        zF0        =   0.81_wp         ! [m/yr] before conversion
        zF0 = zF0 * zconv              ! [mm/s] after  conversion
        IF( kt == nit000 ) THEN
           ALLOCATE( ztstar(jpi,jpj) )   ! Allocation of ztstar
           ALLOCATE( zqsr_dayMean(jpi,jpj) )   ! Allocation of zqsr_dayMean
           IF( rn_srp /= 0._wp )   THEN
              ALLOCATE( zsstar(jpi,jpj) )   ! Allocation of zsstar
              zsstar(:,:) = 37.12_wp * EXP( - gphit(:,:)**2 / 260._wp**2 ) - 1.1_wp * EXP( - gphit(:,:)**2 / 7.5_wp**2 )
           ENDIF
        ENDIF
        ! necessary to compute at each time step because seasonnal variation of Ztstar and solar heat flux
        vtau(:,:) = 0._wp   ! no meridional wind stress
        wndm(:,:) = 0._wp   ! no use of 10 m wind speed
        !
        ! SALT FLUX
        IF( rn_srp /= 0._wp )   THEN
           sfx(:,:) = ( rn_srp * ( tsb(:,:,1,jp_sal) - zsstar(:,:) ) ) * tmask(:,:,1)   ! Restoring salt flux
        ELSE
           sfx (:,:) = 0._wp   ! no salt flux
        ENDIF
        !
        DO jj = 1, jpj
           DO ji = 1, jpi
              utau(ji,jj)    = ztau0 * (             - COS( ( 3 * rpi * gphit(ji,jj) )  / ( 2 * zL )    )  + zatau * EXP( -  gphit(ji,jj)**2        / zdeltatau**2 ) )
              taum(ji,jj)    = ABS( utau(ji,jj) )
              IF( utau(ji,jj) > 0 )   taum(ji,jj) = taum(ji,jj) * 1.3_wp   ! Boost in westerlies for TKE
              !
              ! EMP from Wolfe and Cessi 2014, JPO
              emp (ji,jj) =   zF0 * (               COS( (     rpi * gphit(ji,jj) )  / (     zL )    )  - zaemp * EXP( -  gphit(ji,jj)**2        / zdeltaemp**2 ) )
              ! Seasonnal cycle on T* coming from zcos_sais2
              ztstar(ji,jj)   = (zts_eq - zts_n - zcos_sais2 * zdeltaT ) * COS( ( rpi * gphit(ji,jj) ) * z1_2L )**2 + zts_n + zcos_sais2 * zdeltaT
           END DO
        END DO
        !
        emp(:,:) = rn_emp_prop * emp(:,:)   ! taking the proportionality factor into account
        CALL remove_emp_mean()
        !
        ! Q SOLAR (from Gyre)
        IF( ln_qsr )   THEN
           zqsr_dayMean(:,:) = 230._wp * COS( rpi * (gphit(:,:) - 23.5 * zcos_sais1 ) / ( 0.9_wp * 180._wp ) ) * tmask(:,:,1)
        ELSE
           zqsr_dayMean(:,:) = 0._wp
        ENDIF
        CALL compute_diurn_cycle( kt, zqsr_dayMean, ln_diu_cyc )   ! Adding diurnal cycle if needed
        !
        ! QNS
        ! take (SST - T*) into account, heat content of emp, remove qsr
        qns(:,:) = (  rn_trp * ( tsb(:,:,1,jp_tem) - ztstar(:,:) ) &
             &      - emp(:,:) * tsb(:,:,1,jp_tem) * rcp           &
             &      - zqsr_dayMean(:,:)                         ) * tmask(:,:,1)
     CASE(5)
        ! Forcing inspired from Wolfe and Cessi 2014, JPO
        IF( kt == nit000 .AND. lwp ) THEN
           WRITE(numout,*) 'usrdef_sbc_oce : analytical surface fluxes'
           WRITE(numout,*) '~~~~~~~~~~~~~~'
           WRITE(numout,*) '   Forcing inspired from *Wolfe and Cessi 2014, JPO*, and from the *GYRE configuration*'
           WRITE(numout,*) '      Zonal annual wind'
           WRITE(numout,*) '      Zonal annual E-P'
           IF( ln_qsr )   THEN
              WRITE(numout,*) '      Solar heat flux'
              IF( ln_ann_cyc )   THEN
                 WRITE(numout,*) '         Annual cycle of solar heat flux'
              ELSE
                 WRITE(numout,*) '         Annual average of solar heat flux (December = June)'
              ENDIF
              IF( ln_diu_cyc )   THEN
                 WRITE(numout,*) '         Diurnal cycle of solar heat flux'
              ELSE
                 WRITE(numout,*) '         Daily average of solar heat flux (day = night)'
              ENDIF
           ELSE
              WRITE(numout,*) '      No solar heat flux'
           ENDIF
           WRITE(numout,*) '      Zonal T* (heat forcing proportional to (SST - T*)'
           IF( ln_ann_cyc )   THEN
              WRITE(numout,*) '         Annual cycle for T*'
           ELSE
              WRITE(numout,*) '         Annual average for T* (December = June)'
           ENDIF
           IF( rn_srp /= 0._wp )   THEN
              WRITE(numout,*) '      Zonal annual S* (salt flux proportional to (SSS - S*)'
           ELSE
              WRITE(numout,*) '      No salt restoring'
           ENDIF
        ENDIF
        ! Initialization of parameters
        ! Computation of the day of the year (from Gyre)
        CALL compute_day_of_year( kt, zcos_sais1, zcos_sais2, ln_ann_cyc )
        ! 
        zL = 61                     ! [degrees] Approximative meridional extend of the basin
        ! Wind stress
        ztau0         =   0.1_wp    ! [Pa]
        zatau         =   0.8_wp    ! [no unit]
        zdeltatau     =   5.77_wp   ! [deg North]
        ! T star and qns
        !ztrp          = -40._wp     ! [W/m2/K] retroaction term on heat fluxes 
        zts_eq        =  25._wp     ! [deg C] Temperature at the equator
        zts_n         =   0._wp     ! [deg C] Temperature in the north
        zdeltasst     =  16.22_wp   ! [deg North]
        z1_2L         =   1._wp / (2._wp * zL)
        zdeltaT = 2                 ! [deg C] half difference of temperature during winter and summer in the north (magnitude of the cos) !!rc TODO set in namelist
        ! EMP
        zconv         =   1._wp / ( 86400._wp)   ! convertion factor: 1 mm/day => 1/(3600*24) mm/s
        !!rc TODO put a1, a2 and a3 in namelist
        za1 = -3.24_wp              ! [mm/day] Set the amplitude of EMP at the equator
        za2 = 4.15_wp               ! [mm/day] Set the amplitude of EMP at mid latitude
        za3 = -1.59_wp              ! [mm/day] Set the amplitude of EMP at the northern part
        zphi01 = 0._wp              ! [deg North]
        zphi02 = 20._wp             ! [deg North]
        zphi03 = 50._wp             ! [deg North]
        z1_d1 = 1._wp / 8._wp       ! [1 / deg North]
        zd2 = 30._wp                ! [deg North]
        z1_d2 = 1._wp / zd2         ! [1 / deg North]
        zd3 = 40._wp                ! [deg North]
        z1_d3 = 1._wp / zd3         ! [1 / deg North]
        z1_s = 1._wp / 10._wp       ! streching of the tanh function (i.e. smoothness of the filter)
        za1 = za1 * zconv           ! [mm/s] after  conversion
        za2 = za2 * zconv           ! [mm/s] after  conversion
        za3 = za3 * zconv           ! [mm/s] after  conversion
        !
        IF( kt == nit000 ) THEN
           ALLOCATE( ztstar(jpi,jpj) )   ! Allocation of ztstar
           ALLOCATE( zqsr_dayMean(jpi,jpj) )   ! Allocation of zqsr_dayMean
           IF( rn_srp /= 0._wp )   THEN
              ALLOCATE( zsstar(jpi,jpj) )   ! Allocation of zsstar
              ! See https://www.desmos.com/calculator/qrapqqrbfa
              zsstar(:,:) = 37.12_wp * EXP( - gphit(:,:)**2 / 260._wp**2 ) - 1.1_wp * EXP( - gphit(:,:)**2 / 7.5_wp**2 )
           ENDIF
        ENDIF
        ! necessary to compute at each time step because seasonnal variation of ztstar and solar heat flux
        vtau(:,:) = 0._wp   ! no meridional wind stress
        wndm(:,:) = 0._wp   ! no use of 10 m wind speed
        !
        ! SALT FLUX
        IF( rn_srp /= 0._wp )   THEN
           sfx(:,:) = ( rn_srp * ( tsb(:,:,1,jp_sal) - zsstar(:,:) ) ) * tmask(:,:,1)   ! Restoring salt flux
        ELSE
           sfx (:,:) = 0._wp   ! no salt flux
        ENDIF
        !
        DO jj = 1, jpj
           DO ji = 1, jpi
              utau(ji,jj) = ztau0 * ( -COS((3 * rpi * gphit(ji,jj))/(2 * zL)) + zatau * EXP(-gphit(ji,jj)**2/zdeltatau**2) )
              taum(ji,jj) = ABS( utau(ji,jj) )
              IF( utau(ji,jj) > 0 )   taum(ji,jj) = taum(ji,jj) * 1.3_wp   ! Boost in westerlies for TKE
              !
              ! EMP
              ! See https://www.desmos.com/calculator/v0vpbcc81h
              ! weights
              zweight2 = 0.5 * ( TANH((gphit(ji,jj) - zphi02 + zd2 * 0.5) * z1_s) - TANH((gphit(ji,jj) - zphi02 - zd2 * 0.5) * z1_s) )
              zweight3 = 0.5 * ( TANH((gphit(ji,jj) - zphi03 + zd3 * 0.5) * z1_s) - TANH((gphit(ji,jj) - zphi03 - zd3 * 0.5) * z1_s) )
              ! each component
              zf1 = za1 * EXP(-(gphit(ji,jj) - zphi01)**2 * z1_d1**2               )
              zf2 = za2 * SIN( (gphit(ji,jj) - zphi02)    * z1_d2 * rpi + 0.5 * rpi) * zweight2
              zf3 = za3 * SIN( (gphit(ji,jj) - zphi03)    * z1_d3 * rpi + 0.5 * rpi) * zweight3
              ! total
              emp (ji,jj) = zf1 + zf2 + zf3
              ! The mean is removed later on (no simple analytical function)
              !
              ! T*
              ! See https://www.desmos.com/calculator/zij8tgy5yr
              ! Seasonnal cycle on T* coming from zcos_sais2
              ztstar(ji,jj)   = (zts_eq - (zts_n + zcos_sais2 * zdeltaT) ) * COS( ( rpi * gphit(ji,jj) ) * z1_2L )**2 + (zts_n + zcos_sais2 * zdeltaT)
           END DO
        END DO
        !
        emp(:,:) = rn_emp_prop * emp(:,:)   ! taking the proportionality factor into account
        CALL remove_emp_mean()
        !
        ! Q SOLAR (flux from GYRE)
        ! see https://www.desmos.com/calculator/87duqiuxsf
        IF( ln_qsr )   THEN
           zqsr_dayMean(:,:) = 230._wp * COS( rpi * (gphit(:,:) - 23.5 * zcos_sais1 ) / ( 0.9_wp * 180._wp ) ) * tmask(:,:,1)
           CALL compute_diurn_cycle( kt, zqsr_dayMean, ln_diu_cyc )   ! Adding diurnal cycle if needed
        ELSE
           zqsr_dayMean(:,:) = 0._wp
           qsr(:,:) = 0._wp
        ENDIF
        !
        ! QNS
        ! take (SST - T*) into account, heat content of emp, remove zqsr_dayMean
        qns(:,:) = (  rn_trp * ( tsb(:,:,1,jp_tem) - ztstar(:,:) ) &
             &      - emp(:,:) * tsb(:,:,1,jp_tem) * rcp           &
             &      - zqsr_dayMean(:,:)                         ) * tmask(:,:,1)
     END SELECT
     ! We call lbc_lnk to take the boundaries into account (especially the equator symmetrical condition)
     CALL lbc_lnk( 'usrdef_sbc', taum, 'T',  1. )
     CALL lbc_lnk( 'usrdef_sbc', wndm, 'T',  1. )
     CALL lbc_lnk( 'usrdef_sbc', utau, 'U', -1. )
     CALL lbc_lnk( 'usrdef_sbc', vtau, 'V', -1. )
     !
     CALL lbc_lnk( 'usrdef_sbc', emp , 'T',  1. )
     CALL lbc_lnk( 'usrdef_sbc', sfx , 'T',  1. )
     CALL lbc_lnk( 'usrdef_sbc', qns , 'T',  1. )
     CALL lbc_lnk( 'usrdef_sbc', qsr , 'T',  1. )
   END SUBROUTINE usrdef_sbc_oce

   
   SUBROUTINE usrdef_sbc_ice_tau( kt )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
   END SUBROUTINE usrdef_sbc_ice_tau

   
   SUBROUTINE usrdef_sbc_ice_flx( kt )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
   END SUBROUTINE usrdef_sbc_ice_flx

    
   SUBROUTINE remove_emp_mean()     
     !!---------------------------------------------------------------------
     !!                    ***  ROUTINE remove_emp_mean  ***
     !!              
     !! ** Purpose :   Remove the average on emp (leading to a balanced flux)
     !!
     !! ** Method  :   emp = emp - MEAN(emp)
     !!                where MEAN is an average taking the surface of the cells into account
     !!
     !! ** Action  :   set emp
     !!---------------------------------------------------------------------
     REAL(wp), DIMENSION(jpi,jpj) ::   zSurf
     REAL(wp) ::   zemp_mean
     !!---------------------------------------------------------------------
     !
     ! Setting global E-P to 0
     SELECT CASE ( jperio )
     CASE( 8 )
        ! 8: equatorial symmetric
        !
        ! Computing the surface of the cells
        zSurf( :     , :     ) = 0._wp                                         ! masking the halo + boundary points
        zSurf(2:jpim1,2:jpjm1) = e1t(2:jpim1,2:jpjm1) * e2t(2:jpim1,2:jpjm1)   ! surface of the cells
        WHERE( gphit(:,:) == 0._wp ) zSurf(:,:) = zSurf(:,:) * 0.5_wp          ! At the equator, only half a cell must be counted
        zemp_mean = glob_sum( 'usrdef_sbc', emp(:,:) * zSurf(:,:) ) / glob_sum( 'usrdef_sbc', zSurf(:,:) )
        !
     CASE( 0, 1, 2, 7 )
        ! 0: closed
        ! 1: cyclic east-west
        ! 2:cyclic north-south
        ! 7: cyclic east-west and north-south
        !
        ! Computing the surface of the cells
        zSurf( :     , :     ) = 0._wp                                         ! masking the halo + boundary points
        zSurf(2:jpim1,2:jpjm1) = e1t(2:jpim1,2:jpjm1) * e2t(2:jpim1,2:jpjm1)   ! surface of the cells
        zemp_mean = glob_sum( 'usrdef_sbc', emp(:,:) * zSurf(:,:) ) / glob_sum( 'usrdef_sbc', zSurf(:,:) )
     CASE( 3, 4, 5, 6 )
        ! 3: north fold with T-point pivot
        ! 4: cyclic east-west and north fold with T-point pivot
        ! 5: north fold with F-point pivot
        ! 6: cyclic east-west and north fold with F-point pivot
        CALL ctl_stop( 'usrdef_sbc: The northfold domain boundary conditions (jperio = 3, 4, 5, 6) are not compatible with the implemented EMP' )
     CASE DEFAULT
        CALL ctl_stop( 'usrdef_sbc: jperio is out of accepted range [0, 1, 2, 7, 8]' )
     END SELECT
     emp(:,:) = emp(:,:) - zemp_mean                     ! freshwater flux (=0 in domain average)
     !!rc WRITE(numout,*) 'DEBUG GLOB_SUM(emp*zSurf) after    = ', GLOB_SUM( 'usrdef_sbc', emp(:,:) * zSurf(:,:) )
     !
   END SUBROUTINE remove_emp_mean

   
   SUBROUTINE compute_day_of_year( kt, pcos_sais1, pcos_sais2, ll_ann_cyc )
     !!---------------------------------------------------------------------
     !!                    ***  SUBROUTINE compute_day_of_year  ***
     !!              
     !! ** Purpose :   Computation of the day of the year (from Gyre) as a cosine
     !!                pcos_sais1 is for heat flux, and is min the 21th of December, max the 21th of June
     !!                pcos_sais2 is for T*       , and is min the 21th of January , max the 21th of July
     !!
     !! ** Method  :   
     !!
     !! ** Action  : 
     !!---------------------------------------------------------------------
     INTEGER , INTENT(in   ) ::   kt      ! ocean time step
     REAL(wp), INTENT(  out) ::   pcos_sais1, pcos_sais2   ! cosine of the day of year (1 is for solar heat flux, 2 is for T* cycle)
     LOGICAL , INTENT(in   ), optional ::   ll_ann_cyc    ! if .false., the cos are set to zero. Default behaviour is true
     !
     LOGICAL  ::   ld_compute   ! local variable to know if we need to compute the cosines or set them to 0
     ! Variables to get the day of the year
     INTEGER  ::   zyear0                 ! initial year 
     INTEGER  ::   zmonth0                ! initial month
     INTEGER  ::   zday0                  ! initial day
     INTEGER  ::   zday_year0             ! initial day since january 1st
     REAL(wp) ::   ztime                  ! time in hour
     REAL(wp) ::   ztimemax , ztimemin    ! 21th June, and 21th decem. if date0 = 1st january
     REAL(wp) ::   ztimemax1, ztimemin1   ! 21th June, and 21th decem. if date0 = 1st january
     REAL(wp) ::   ztimemax2, ztimemin2   ! 21th June, and 21th decem. if date0 = 1st january
     REAL(wp) ::   zyydd                 ! number of days in one year
     !!---------------------------------------------------------------------
     !
     IF( PRESENT(ll_ann_cyc) )   THEN
        ld_compute = ll_ann_cyc
     ELSE
        ld_compute = .TRUE.
     ENDIF
     !
     IF( ld_compute )   THEN
        zyydd = REAL(nyear_len(1),wp)
        zyear0     =   ndate0 / 10000._wp                                ! initial year
        zmonth0    = ( ndate0 - zyear0 * 10000._wp ) / 100._wp           ! initial month
        zday0      =   ndate0 - zyear0 * 10000._wp - zmonth0 * 100._wp   ! initial day betwen 1 and 30
        zday_year0 = ( zmonth0 - 1._wp ) * 30._wp + zday0                ! initial day betwen 1 and 360
        !
        ! current day (in hours) since january the 1st of the current year
        ztime = REAL( kt ) * rdt / (rmmss * rhhmm)   &       !  total incrementation (in hours)
             &      - (nyear  - 1) * rjjhh * zyydd              !  minus years since beginning of experiment (in hours)

        ztimemax1 = ((5.*30.)+21.)* 24.                      ! 21th june     at 24h in hours
        ztimemin1 = ztimemax1 + rjjhh * zyydd / 2            ! 21th december        in hours
        ztimemax2 = ((6.*30.)+21.)* 24.                      ! 21th july     at 24h in hours
        ztimemin2 = ztimemax2 - rjjhh * zyydd / 2            ! 21th january         in hours
        !                                                    ! NB: rjjhh * zyydd / 4 = one seasonal cycle in hours
        !
        ! 1/2 period between 21th June and 21th December and between 21th July and 21th January (1 for solar heat flux, 2 for T*)
        pcos_sais1 = COS( (ztime - ztimemax1) / (ztimemin1 - ztimemax1) * rpi )
        pcos_sais2 = COS( (ztime - ztimemax2) / (ztimemax2 - ztimemin2) * rpi )
     ELSE
        pcos_sais1 = 0._wp
        pcos_sais2 = 0._wp
     ENDIF
   END SUBROUTINE compute_day_of_year


   SUBROUTINE compute_diurn_cycle( kt, pqsr_dayMean, ll_diu_cyc )
     !!---------------------------------------------------------------------
     !!                    ***  SUBROUTINE compute_diurn_cycle  ***
     !!              
     !! ** Purpose :   Set the value of qsr.
     !!                If ll_diu_cyc is .true. or is not present, use the diurnal cycle.
     !!                If ll_diu_cyc is .false. use the daily mean.
     !!
     !! ** Method  :   
     !!
     !! ** Action  : 
     !!---------------------------------------------------------------------
     INTEGER ,                     INTENT(in   ) ::   kt      ! ocean time step
     REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pqsr_dayMean
     LOGICAL , optional          , INTENT(in   ) ::   ll_diu_cyc    ! if .false. use the mean value, if .true. use the diurnal cycle.
     !
     LOGICAL  ::   ld_compute   ! local variable to know if we need to compute the diurnal cycle
     !!---------------------------------------------------------------------
     !
     ! Adding diurnal cycle if needed
     IF( PRESENT(ll_diu_cyc) )   THEN
        ld_compute = ll_diu_cyc
     ELSE
        ld_compute = .TRUE.
     ENDIF
     IF( ld_compute )   THEN
        IF(  kt == nit000 )   nday_qsr = -1
        qsr(:,:) = sbc_dcy( pqsr_dayMean(:,:) ) * tmask(:,:,1)
     ELSE
        qsr(:,:) =  pqsr_dayMean(:,:)
     ENDIF
   END SUBROUTINE compute_diurn_cycle
   !!======================================================================
END MODULE usrdef_sbc
