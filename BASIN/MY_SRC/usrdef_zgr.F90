MODULE usrdef_zgr
   !!======================================================================
   !!                       ***  MODULE  usrdef_zgr  ***
   !!
   !!                      ===  BASIN configuration  ===
   !!
   !! User defined : vertical coordinate system of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2017-11  (J. Chanut)  Original code
   !!                 ! 2019-05  (R.Caneill and G.Madec) Adaptation
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_zgr   : user defined vertical coordinate system
   !!      zgr_z      : reference 1D z-coordinate 
   !!      zgr_top_bot: ocean top and bottom level indices
   !!      zgr_zco    : 3D vertical coordinate in pure z-coordinate case
   !!---------------------------------------------------------------------
   USE oce    , ONLY:            ! ocean variables
   USE dom_oce, ONLY:        ! ocean domain
   USE phycst         ! physical constants
   !
   USE usrdef_nam, ONLY: nn_botcase, rn_distlam, rn_H, rn_hborder, ln_sco_nam, ln_zco_nam, ln_zps_nam, nn_ztype, ln_equ_flat
   !!rc USE depth_e3       ! depth <=> e3
   USE zgr_lib    ! tools for vertical coordinate
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O library !!rc
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_zgr        ! called by domzgr.F90

  !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_zgr.F90 10425 2018-12-19 21:54:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS             

   SUBROUTINE usr_def_zgr( ld_zco  , ld_zps  , ld_sco  , ld_isfcav,    &   ! type of vertical coordinate
      &                    pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d  ,    &   ! 1D reference vertical coordinate
      &                    pdept , pdepw ,                             &   ! 3D t & w-points depth
      &                    pe3t  , pe3u  , pe3v   , pe3f ,             &   ! vertical scale factors
      &                    pe3w  , pe3uw , pe3vw         ,             &   !     -      -      -
      &                    k_top , k_bot    )                              ! top & bottom ocean level
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE usr_def_zgr  ***
      !!
      !! ** Purpose :   User defined the vertical coordinates
      !!
      !!----------------------------------------------------------------------
      LOGICAL                   , INTENT(out) ::   ld_zco, ld_zps, ld_sco      ! vertical coordinate flags
      LOGICAL                   , INTENT(out) ::   ld_isfcav                   ! under iceshelf cavity flag
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pe3t_1d , pe3w_1d           ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pdept, pdepw                ! grid-point depth        [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors  [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pe3w , pe3uw, pe3vw         ! i-scale factors 
      INTEGER , DIMENSION(:,:)  , INTENT(out) ::   k_top, k_bot                ! first & last ocean level
      !
      REAL(wp) ::   zdzmin    ! minimum value of e3 at the surface   [m]
      REAL(wp) ::   zkth      ! position of the inflexion point
      INTEGER  ::   ikconst   ! Number of levels with pure z-coordinate (i.e. constant depth)
      REAL(wp) ::   zacr      ! slope of the tanh
      REAL(wp) ::   zHmax
      LOGICAL  ::   ll_sco_pure, ll_sco_mi96
      REAL(wp), DIMENSION(jpi,jpj)           ::   zbathy, z2d   ! bathymetry, 2D variable
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr : BASIN configuration'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      ! Vertical coordinate type
      ld_zco      = ln_zco_nam
      ld_zps      = ln_zps_nam
      ld_sco      = ln_sco_nam
      ld_isfcav   = .FALSE.   ! no ice-shelves
      ll_sco_pure = .FALSE.   ! not pure s stretching
      ll_sco_mi96 = .TRUE.    ! using mi96 function for stretching, with layer of z-coordinates at the surface
      !
      zHmax = rn_H   ! Maximum depth of the basin
      CALL zgr_bat( zHmax, zbathy )   ! User creation of bathymetry
      !
      ! Build the vertical coordinate system
      ! ------------------------------------
      !
      !   === z-coordinate   ===   !
      IF( ld_zco ) THEN
         !
         CALL zgr_z1D( zHmax, nn_ztype, pdept_1d, pdepw_1d )   ! Reference z-coordinate system
         !
         !                                                       ! z-coordinate (3D arrays) from the 1D z-coord.
         CALL zgr_zco( pdept_1d, pdepw_1d,                   &   ! <==>  1D reference vertical coordinate
            &          pe3t_1d , pe3w_1d ,                   &   ! ==>>  1D t & w-points vertical scale factors  
            &          pdept   , pdepw   ,                   &   !       3D t & w-points depth
            &          pe3t    , pe3u    , pe3v    , pe3f,   &   !       vertical scale factors
            &          pe3w    , pe3uw   , pe3vw           )     !           -      -      -
         !
         CALL zgr_msk_top_bot( pdept_1d, zbathy, k_top, k_bot )                 ! masked top and bottom ocean t-level indices
         !
      ENDIF
      !
      !   ===   s-coordinate   ===   !
      IF( ld_sco ) THEN
         !
         IF( ll_sco_pure ) THEN   ! z-coordinate (3D arrays) from the 1D z-coord., pure s coordinate
            CALL zgr_z1D( zHmax, nn_ztype, pdept_1d, pdepw_1d )   ! Reference s-coordinate system
            !
            CALL zgr_sco_pure( pdept_1d, pdepw_1d,                   &   ! <==>  1D reference vertical coordinate
               &               zbathy  , zHmax   ,                   &   ! <<==  bathymetry and its maximum
               &               pe3t_1d , pe3w_1d ,                   &   ! ==>>  1D t & w-points vertical scale factors  
               &               pdept   , pdepw   ,                   &   !       3D t & w-points depth
               &               pe3t    , pe3u    , pe3v    , pe3f,   &   !       vertical scale factors
               &               pe3w    , pe3uw   , pe3vw           )     !           -      -      -
         ENDIF
         !
         IF( ll_sco_mi96 ) THEN   ! z-coordinate (3D arrays), mixed z and s coordinates
            zdzmin  = 10._wp
            zkth    = 26.8_wp
            ikconst = 5
            zacr    = 10.5_wp
            CALL zgr_sco_mi96( pdept_1d, pdepw_1d,                   &   ! ==>>  1D reference vertical coordinate
               &               zbathy  ,                             &   ! <<==  bathymetry
               &               zdzmin  , zkth    , ikconst , zacr,   &   ! <<==  parameters for the tanh stretching
               &               pe3t_1d , pe3w_1d ,                   &   ! ==>>  1D t & w-points vertical scale factors  
               &               pdept   , pdepw   ,                   &   !       3D t & w-points depth
               &               pe3t    , pe3u    , pe3v    , pe3f,   &   !       vertical scale factors
               &               pe3w    , pe3uw   , pe3vw           )     !           -      -      -
         ENDIF
         !
         ! Slope tests
         CALL zgr_test_slopes( zbathy, pdept, pe3u, pe3v )
         !
         !                                           ! masked top and bottom ocean t-level indices
         zbathy(:,:) = 1._wp                                     ! Use zbathy as land-sea mask to compute k_top and k_bot
         CALL lbc_lnk( 'usr_def_zgr', zbathy, 'T', 1._wp )
         k_top(:,:) =         INT( zbathy(:,:) )
         k_bot(:,:) = jpkm1 * INT( zbathy(:,:) )  
         !
      ENDIF
   END SUBROUTINE usr_def_zgr

   
   SUBROUTINE zgr_bat( pH, pbathy )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_bat  ***
      !!
      !! ** Purpose :  Create analyticaly the bathymetry at T points, according to the nn_botcase option
      !!
      !! ** Method  :  Depending on nn_botcase will create different types of bathymetry.
      !!               The lon/lat coordinates are normalized between -1 and 1, and an analytical function is applied
      !!                  If nn_botcase == 0: flat bottom at a depth pH
      !!                  If nn_botcase == 1: bowl shape with a cosh
      !!                  If nn_botcase == 2: bowl shape with a (1-x**4)
      !!----------------------------------------------------------------------
      REAL(wp)                , INTENT(inout) ::   pH       ! ocean depth maximum   [m]
      REAL(wp), DIMENSION(:,:), INTENT(  out) ::   pbathy   ! Ocean bathymetry      [m]
      !
      INTEGER  ::   ji, jj             ! dummy loop indices
      REAL(wp) ::   zmaxlam, zmaxphi, zminlam, zminphi
      REAL(wp) ::   zx, zy                                   ! local scalars
      REAL(wp) ::   z1_H, z1_dLam, z1_dPhi, zxnorm, zynorm   ! local scalars
      REAL(wp) ::   zdistPhi, ztaper                         ! local scalars
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d
      !!----------------------------------------------------------------------
      !
      IF(lwp)   THEN
         WRITE(numout,*)
         WRITE(numout,*)    '    zgr_bat : defines the bathymetry at T points.'
         WRITE(numout,*)    '    ~~~~~~~'
         WRITE(numout,*)    '       BASIN case : flat bottom or bowl shape'
         WRITE(numout,*)    '       nn_botcase = ', nn_botcase
      ENDIF
      !
      SELECT CASE(nn_botcase)
      !
      CASE(0)   ! flat bottom
         IF(lwp)   THEN
            WRITE(numout,*) '          i.e. flat bottom'
            WRITE(numout,*) '          Depth of the bathymetry   pH   = ', pH
         ENDIF
         pbathy(:,:) = pH
      CASE(1)   ! new bowl shape with exponentials
         IF(lwp)   THEN
            WRITE(numout,*) '          i.e. bowl shape in cosh(x)'
            WRITE(numout,*) '          see https://www.desmos.com/calculator/onczibqzb4'
            WRITE(numout,*) '          Asymptotical depth of the exponentials              pHasymp   = ', pH,      ' m'
            WRITE(numout,*) '          Depth of the bathymetry at the coast             rn_hborder   = ', rn_hborder, ' m'
            WRITE(numout,*) '          Typical length scale of the slope in longitude   rn_distLam   = ', rn_distLam,   ' degrees'
            WRITE(numout,*) '                                                                     '
            WRITE(numout,*) '             |    ^                                   ^              '
            WRITE(numout,*) '             | hborder                                |              '
            WRITE(numout,*) '             |    v                                   |              '
            WRITE(numout,*) '             |*                                       |              '
            WRITE(numout,*) '             | *                                     pH              '
            WRITE(numout,*) '             |  *                                     |              '
            WRITE(numout,*) '             |   *                                    |              '
            WRITE(numout,*) '             |    \**                                 |              '
            WRITE(numout,*) '             |     \ ***                              |              '
            WRITE(numout,*) '             |      \   *****                         v              '
            WRITE(numout,*) '             |       \       *********************************       '
            WRITE(numout,*) '              < zdist>                                               '
         ENDIF
         !
         CALL zgr_get_boundaries( zmaxlam, zminlam, zmaxphi, zminphi )
         !rn_distLam = 3._wp
         zdistPhi =  COS( rad * zmaxphi ) * rn_distLam
         z1_dLam = 1._wp / rn_distLam
         z1_dPhi = 1._wp / zdistPhi
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! Here we use zx and zy not as normalization of coordinates, but as normalyzed bathy on longitude and latitude
               zxnorm = 1._wp + EXP( ( zminlam - zmaxlam) * z1_dLam )
               zynorm = 1._wp + EXP( ( zminphi - zmaxphi) * z1_dPhi )
               !
               zx = 1._wp - 1._wp / zxnorm * (  EXP(   ( glamt(ji,jj) - zmaxlam ) * z1_dLam )   &
                  &                           + EXP( - ( glamt(ji,jj) - zminlam ) * z1_dLam ) )
               zy = 1._wp - 1._wp / zynorm * (  EXP(   ( gphit(ji,jj) - zmaxphi ) * z1_dPhi )   &
                  &                           + EXP( - ( gphit(ji,jj) - zminphi ) * z1_dPhi ) )
               !
               pbathy(ji,jj) = zx * zy * ( pH - rn_hborder ) + rn_hborder
            END DO
         END DO
      CASE(2)   ! bowl shape 
         IF(lwp)   THEN
            WRITE(numout,*) '          i.e. bowl shape in (1-x**4)'
         ENDIF
         CALL zgr_get_boundaries( zmaxlam, zminlam, zmaxphi, zminphi )
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! Normalization between -1 and 1
               zx = 2._wp*(glamt(ji,jj) - zminlam)/(zmaxlam - zminlam) - 1._wp
               zy = 2._wp*(gphit(ji,jj) - zminphi)/(zmaxphi - zminphi) - 1._wp
               pbathy(ji,jj) =  0.5_wp * pH * ( 1 + (1-zx**4)*(1-zy**4) )
            END DO
         END DO
      CASE DEFAULT
         CALL ctl_stop( 'zgr_bat: The chosen case for the bathymetry does not exist (nn_botcase).' )
      !
      END SELECT
      !
      ! Check that zHmax is really the maximum on the inner domain
      z2d(:,:) = pbathy(:,:)
      CALL lbc_lnk( 'usr_def_zgr', z2d, 'T', 1._wp )   ! put the land boundaries to 0
      pH = MAXVAL(z2d(:,:))
      IF( lk_mpp )   CALL mpp_max( 'zgr_bat', pH )
      IF(lwp)   WRITE(numout,*) '      Effective Hmax   = ', pH, ' m      while rn_H   = ', rn_H, ' m'
      !
      ! Flattening the bathymetry at the equator
      IF( ln_equ_flat )   THEN
         IF(lwp)   WRITE(numout,*) '      s-coordinates are tapered in vicinity of the Equator'
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztaper = EXP( -(gphit(ji,jj)/8._wp)**2 )
               pbathy(ji,jj) = pH * ztaper + pbathy(ji,jj) * ( 1._wp - ztaper )
            END DO
         END DO
      ENDIF
      !
   END SUBROUTINE zgr_bat


   SUBROUTINE zgr_get_boundaries( pmaxlam, pminlam, pmaxphi, pminphi , cdpoint)
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_get_boundaries  ***
      !!
      !! ** Purpose :   Compute the maximum and minimum values of latitude and longitude in the global inner domain
      !!
      !! ** Method  :   We take only the inner domain for the min/max of phi and lambda
      !!                   max phi is taken on gphiv( 2:jpim1 , 1:jpjm1 ) so the last unmasked V point represents the boundary of the bathymetry
      !!                   max lam -    -   -  glamu( 1:jpim1 , 2:jpjm1 ) -  -   last   -      U   -       -       -   -       -   -      -
      !!                   min lam -    -   -  glamu( 1:jpim1 , 2:jpjm1 ) -  -  first   -      U   -       -       -   -       -   -      -
      !!                   min phi -    -   -  gphiv( 2:jpim1 , 1:jpjm1 ) -  -  first   -      V   -       -       -   -       -   -      -
      !!                If jperio==8 (symmetrical at the equator), returns pminphi=-pmaxphi
      !!
      !!                If cdpoint is present (can be any character), we do the calculations on the inner T points of the domain
      !!----------------------------------------------------------------------
      CHARACTER, OPTIONAL, INTENT(in   ) ::   cdpoint
      REAL(wp)           , INTENT(  out) ::   pmaxlam, pminlam, pmaxphi, pminphi
      !
      REAL(wp), DIMENSION(4) :: zmaxmin
      !!----------------------------------------------------------------------
      IF( PRESENT(cdpoint) ) THEN   ! We take the min and max values at T points
         pmaxlam    = MAXVAL( glamt( 2:jpim1 , 2:jpjm1 ) )
         pminlam    = MINVAL( glamt( 2:jpim1 , 2:jpjm1 ) )
         pmaxphi    = MAXVAL( gphit( 2:jpim1 , 2:jpjm1 ) )
         IF( jperio /= 8 ) THEN   ! not symmetrical at the equator
            pminphi = MINVAL( gphit( 2:jpim1 , 2:jpjm1 ) )
         ELSE
            pminphi = -pmaxphi
         ENDIF
      ELSE                        ! We take min and max values at the coast (i.e. U points for longitude, V points for latitude)
         pmaxlam    = MAXVAL( glamu( 1:jpim1 , 2:jpjm1 ) )
         pminlam    = MINVAL( glamu( 1:jpim1 , 2:jpjm1 ) )
         pmaxphi    = MAXVAL( gphiv( 2:jpim1 , 1:jpjm1 ) )
         IF( jperio /= 8 ) THEN   ! not symmetrical at the equator
            pminphi = MINVAL( gphiv( 2:jpim1 , 1:jpjm1 ) )
         ELSE
            pminphi = -pmaxphi
         ENDIF
      ENDIF
      IF( lk_mpp ) THEN
         ! max and min over the global domain
         ! the array is just used for the mpp communications (to reduce the time cost)
         !
         ! We use here the propertie for the min: MIN(a,b) = -MAX(-a,-b) if a,b >= 0
         ! we know that  90+minphi > 0,
         !              360+minlam > 0
         zmaxmin(1) = pmaxphi
         zmaxmin(2) = pmaxlam
         zmaxmin(3) = - (  90 + pminphi )
         zmaxmin(4) = - ( 360 + pminlam )
         CALL mpp_max( 'usrdef_zgr', zmaxmin )
         pmaxphi =   zmaxmin(1)
         pmaxlam =   zmaxmin(2)
         pminphi = - zmaxmin(3) -  90
         pminlam = - zmaxmin(4) - 360
      ENDIF
      IF( lwp ) WRITE(numout,*)
      IF( lwp ) WRITE(numout,*) '      zgr_get_boundaries : '
      IF( lwp ) WRITE(numout,*) '      ~~~~~~~~~~~~~~~~~~ '
      IF( lwp ) WRITE(numout,*) '         pminlam = ', pminlam
      IF( lwp ) WRITE(numout,*) '         pmaxlam = ', pmaxlam
      IF( lwp ) WRITE(numout,*) '         pminphi = ', pminphi
      IF( lwp ) WRITE(numout,*) '         pmaxphi = ', pmaxphi
   END SUBROUTINE zgr_get_boundaries

       
   SUBROUTINE zgr_z1D( pHmax, pztype, pdept_1d, pdepw_1d )   ! 1D reference vertical coordinate
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_z  ***
      !!
      !! ** Purpose :   set the 1D depth of model levels and the resulting 
      !!              vertical scale factors.
      !!
      !! ** Method  :   1D z-coordinate system (use in all type of coordinate)
      !!       The depth of model levels is set from dep(k), an analytical function:
      !!                   w-level: depw_1d  = dep(k)
      !!                   t-level: dept_1d  = dep(k+0.5)
      !!       The scale factors are the discrete derivative of the depth:
      !!                   e3w_1d(jk) = dk[ dept_1d ] 
      !!                   e3t_1d(jk) = dk[ depw_1d ]
      !!           with at top and bottom :
      !!                   e3w_1d( 1 ) = 2 * ( dept_1d( 1 ) - depw_1d( 1 ) )
      !!                   e3t_1d(jpk) = 2 * ( dept_1d(jpk) - depw_1d(jpk) )
      !!       The depth are then re-computed from the sum of e3. This ensures 
      !!    that depths are identical when reading domain configuration file. 
      !!    Indeed, only e3. are saved in this file, depth are compute by a call
      !!    to the e3_to_depth subroutine.
      !!
      !!       Here the Madec & Imbard (1996) function is used.
      !!
      !! ** Action  : - pdept_1d, pdepw_1d : depth of T- and W-point (m)
      !!              - pe3t_1d , pe3w_1d  : scale factors at T- and W-levels (m)
      !!
      !! Reference : Marti, Madec & Delecluse, 1992, JGR, 97, No8, 12,763-12,766.
      !!             Madec and Imbard, 1996, Clim. Dyn.
      !!----------------------------------------------------------------------
      REAL(wp)               , INTENT(in   ) ::   pHmax                ! ocean depth maximum   [m]
      INTEGER                , INTENT(in   ) ::   pztype               ! type of z grid (0: constant)
      REAL(wp), DIMENSION(:) , INTENT(  out) ::   pdept_1d, pdepw_1d   ! 1D grid-point depth   [m]
      !
      INTEGER  ::   jk       ! dummy loop indices
      REAL(wp) ::   zd       ! local scalar
      REAL(wp) ::   zt, zw   ! local scalars
      REAL(wp) ::   zsur, za0, za1, zkth, zacr   ! Values for the Madec & Imbard (1996) function  
      !!----------------------------------------------------------------------
      !
      zd = pHmax / REAL(jpkm1,wp)
      !
      IF(lwp) THEN            ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) '    zgr_z   : Reference vertical z-coordinates '
         WRITE(numout,*) '    ~~~~~~~'
      ENDIF
      !
      ! 1D Reference z-coordinate    (using Madec & Imbard 1996 function)   !!rc TODO set non uniform spacing cf ORCA
      ! -------------------------
      !
      SELECT CASE( pztype )
      CASE( 0 )   ! Uniform vertical grid
         IF(lwp) THEN
            WRITE(numout,*) '       BASIN case : uniform vertical grid :'
            WRITE(numout,*) '                     with thickness = ', zd
         ENDIF
         pdepw_1d(1) = 0._wp
         pdept_1d(1) = 0.5_wp * zd
         ! 
         DO jk = 2, jpk          ! depth at T and W-points
            pdepw_1d(jk) = pdepw_1d(jk-1) + zd 
            pdept_1d(jk) = pdept_1d(jk-1) + zd 
         END DO
      CASE( 1 )   ! Non uniform spacing
         IF(lwp)   WRITE(numout,*) '       BASIN case : non uniform vertical grid'
         CALL ctl_stop( 'zgr_z1D: The chosen case (pztype = 1) has not been implemeted yet' )
      CASE DEFAULT
         CALL ctl_stop( 'zgr_z1D: The chosen case for the vertical 1D grid does not exist' )
      END SELECT
      !
   END SUBROUTINE zgr_z1D


   SUBROUTINE zgr_msk_top_bot( pdept_1d , pbathy, k_top , k_bot )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE zgr_msk_top_bot  ***
      !!
      !! ** Purpose :   set the masked top and bottom ocean t-levels for z-coordinates
      !!
      !! ** Method  :   BASIN case = closed or south symmetrical box ocean without ocean cavities
      !!                      k_top = 1                                       except along boundaries according to jperio
      !!                      k_bot = {first point under the ocean bottom}    except along boundaries according to jperio
      !!
      !! ** Action  : - k_top : first wet ocean level index
      !!              - k_bot : last  wet ocean level index
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)  , INTENT(in   ) ::   pdept_1d
      REAL(wp), DIMENSION(:,:), INTENT(in   ) ::   pbathy          ! bathymetry   [m]
      INTEGER , DIMENSION(:,:), INTENT(  out) ::   k_top , k_bot   ! first & last wet ocean level
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d      ! local 2D variable
      INTEGER :: jk   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '    zgr_top_bot : defines the top and bottom wet ocean levels for z-coordinates.'
      IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) '       BASIN case : = closed or south symmetrical box ocean without ocean cavities'
      IF(lwp) WRITE(numout,*) '          k_top = 1                                       except along boundaries according to jperio'
      IF(lwp) WRITE(numout,*) '          k_bot = {first point under the ocean bottom}    except along boundaries according to jperio'
      !
      !
      ! bottom ocean mask computed from the depth of grid-points
      k_bot(:,:) = jpkm1   ! initialisation of k_bot to its maximum authorized value
      DO jk = 1, jpkm1
         WHERE( pdept_1d(jk) < pbathy(:,:) .AND. pbathy(:,:) <= pdept_1d(jk+1) )   z2d(:,:) = REAL(jk,wp)
      END DO
      !
      CALL lbc_lnk( 'usrdef_zgr', z2d, 'T', 1. )           ! set surrounding land point to zero (depending on jperio value)
      !
      k_bot(:,:) = INT( z2d(:,:) )
      k_top(:,:) = MIN( 1 , k_bot(:,:) )     ! = 1 over the ocean point, =0 elsewhere
      !
   END SUBROUTINE zgr_msk_top_bot


   SUBROUTINE zgr_test_slopes( pbathy, pdept, pe3u, pe3v )
!!rc TODO implement multi processors
!!rc TODO remove masked U and V points
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE zgr_test_slopes  ***
      !!
      !! ** Purpose :   TODO
      !!
      !! ** Method  :   TODO
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:)  , INTENT(in) ::   pbathy
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   pdept       ! grid-point depth        [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(in) ::   pe3u , pe3v  ! vertical scale factors  [m]
      !
      INTEGER  ::   inum, ji, jj, jk, i_criteria_not_respected_i, i_criteria_not_respected_j
      REAL(WP) ::   zsi_bat, zsj_bat, zsmax_bat   ! vars for checking the slopes (    bathymetry       relative to the    horizontal          )
      REAL(WP) ::   zsi_dia, zsj_dia, zsmin_dia   ! vars for checking the slopes (diagonal of the grid relative to the slope of the grid      )
      REAL(WP) ::   zsi_gri, zsj_gri, zsmax_gri   ! vars for checking the slopes (      grid           relative to a horizontal stratification)
      REAL(wp), DIMENSION(jpkm1)             ::   zsi1d_dia, zsj1d_dia, zsi1d_gri, zsj1d_gri
      REAL(wp), DIMENSION(jpi-2,jpj-2,jpkm1) ::   zglo_ri, zglo_rj    ! global ratio of slopes of the grid
      REAL(wp), DIMENSION(jpi,jpj)           ::   z2d   ! 2D variable
      !!----------------------------------------------------------------------
      ! maximum slope of bathymetry
      ! We take here only the inner values
      zsi_bat   = MAXVAL(        ABS( (pbathy(3:jpi   , 2:jpjm1)      - pbathy(2:jpim1 , 2:jpjm1))      * r1_e1u(2:jpim1 , 2:jpjm1) ) )
      zsj_bat   = MAXVAL(        ABS( (pbathy(2:jpim1 , 3:jpj  )      - pbathy(2:jpim1 , 2:jpjm1))      * r1_e2v(2:jpim1 , 2:jpjm1) ) )
      !
      DO jk = 1, jpkm1
         ! Maximum slope of the grid relative to a horizontal stratification 
         zsi1d_gri(jk) = MAXVAL( ABS( (pdept( 3:jpi   , 2:jpjm1 , jk) - pdept( 2:jpim1 , 2:jpjm1 , jk)) * r1_e1u(2:jpim1 , 2:jpjm1) ) )
         zsj1d_gri(jk) = MAXVAL( ABS( (pdept( 2:jpim1 , 3:jpj   , jk) - pdept( 2:jpim1 , 2:jpjm1 , jk)) * r1_e2v(2:jpim1 , 2:jpjm1) ) )
         ! minimum slopes in the diagonals of the grid cells
         zsi1d_dia(jk) = MINVAL( pe3u(2:jpim1 , 2:jpjm1 , jk) * r1_e1u(2:jpim1 , 2:jpjm1) )
         zsj1d_dia(jk) = MINVAL( pe3v(2:jpim1 , 2:jpjm1 , jk) * r1_e2v(2:jpim1 , 2:jpjm1) )
         ! global ratio of slope (must be <= 1)
         !!rc WRITE(numout,*) 'usrdef_zgr DEBUG   jpi, jpj', jpi, jpj
         
         zglo_ri(:,:,jk) = ABS( (pdept( 3:jpi   , 2:jpjm1 , jk) - pdept( 2:jpim1 , 2:jpjm1 , jk)) * r1_e1u(2:jpim1 , 2:jpjm1) ) * ( e1u(2:jpim1 , 2:jpjm1) / pe3u(2:jpim1 , 2:jpjm1 , jk))
         zglo_rj(:,:,jk) = ABS( (pdept( 2:jpim1 , 3:jpj   , jk) - pdept( 2:jpim1 , 2:jpjm1 , jk)) * r1_e2v(2:jpim1 , 2:jpjm1) ) * ( e2v(2:jpim1 , 2:jpjm1) / pe3v(2:jpim1 , 2:jpjm1 , jk))
         !
         !!rc zglo_ri(:,:,jk) = ABS( (pdept( 2:jpi   , 1:jpjm1 , jk) - pdept( 1:jpim1 , 1:jpjm1 , jk)) * r1_e1u(1:jpim1 , 1:jpjm1) ) * ( e1u(1:jpim1 , 1:jpjm1) / pe3u(1:jpim1 , 1:jpjm1 , jk))
         !!rc zglo_rj(:,:,jk) = ABS( (pdept( 1:jpim1 , 2:jpj   , jk) - pdept( 1:jpim1 , 1:jpjm1 , jk)) * r1_e2v(1:jpim1 , 1:jpjm1) ) * ( e2v(1:jpim1 , 1:jpjm1) / pe3v(1:jpim1 , 1:jpjm1 , jk))
      END DO
      !
      zsmax_bat = MAX(          zsi_bat ,          zsj_bat  )
      zsmax_gri = MAX( MAXVAL(zsi1d_gri), MAXVAL(zsj1d_gri) )
      zsmin_dia = MIN( MINVAL(zsi1d_dia), MINVAL(zsj1d_dia) )
      !
      i_criteria_not_respected_i = SUM( MERGE( 1,0,SUM(MERGE(1,0,zglo_ri > 3), 3) >= 1 ) ) ! get 1 if zglo_ri > 3 somewhere in the column, 0 else
      i_criteria_not_respected_j = SUM( MERGE( 1,0,SUM(MERGE(1,0,zglo_rj > 3), 3) >= 1 ) )
      !
      ! Open and create a netcdf file
      CALL iom_open( 'slope_ratio', inum, ldwrt = .TRUE.)
      z2d(:,:) = 0._wp
      z2d(2:jpim1,2:jpjm1) = SUM(MERGE(1,0,zglo_ri > 3), 3)
      CALL iom_rstput( 0, 0, inum, 'ratio_slopei', z2d)
      z2d(2:jpim1,2:jpjm1) = SUM(MERGE(1,0,zglo_rj > 3), 3)
      CALL iom_rstput( 0, 0, inum, 'ratio_slopej', z2d)
      CALL iom_close( inum )
      !
      ! If mpp, we compute the min / max over the whole domain
      IF( lk_mpp ) THEN
         CALL mpp_max( 'usr_def_zgr', zsmax_bat )
         CALL mpp_max( 'usr_def_zgr', zsmax_gri )
         CALL mpp_min( 'usr_def_zgr', zsmin_dia )
      ENDIF
      !
      IF( lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' usr_def_zgr : Diag for the slopes (used in lateral diffusion)'
         WRITE(numout,*) ' ~~~~~~~~~~~'
         WRITE(numout,*) '    Maximum slope of the bathymetry                        = ', zsmax_bat * 100, ' per cent'
         WRITE(numout,*) '    Maximum slope of the grid                              = ', zsmax_gri * 100, ' per cent'
         WRITE(numout,*) '    Minimum value of e3/e1 (i.e. slope of the diagonals)   = ', zsmin_dia * 100    , ' per cent'
         WRITE(numout,*) '    Corresponding maximum accepted slope                   = ', zsmin_dia * 100 * 3, ' per cent, with a factor 3'
         WRITE(numout,*) '    Corresponding maximum accepted slope                   = ', zsmin_dia * 100 * 5, ' per cent, with a factor 5'
         WRITE(numout,*) '    Maximum value of ratio (slope grid relative to horiz) / (slope diagonal relative to slope of the grid)'
         WRITE(numout,*) '                                                           = ', MAXVAL(zglo_ri), ' must be <= to 3 or 5 depending on the factor TODO use a var + mpp_max'
         WRITE(numout,*) '    Number of i point where the criteria is not respected  = ', i_criteria_not_respected_i
         WRITE(numout,*) '    Number of j point where the criteria is not respected  = ', i_criteria_not_respected_j
         WRITE(numout,*) '    Number tot of none masked horizontal points            = ', (jpj-2)*(jpi-2)
      ENDIF
      !
   END SUBROUTINE zgr_test_slopes
       
   !!======================================================================
END MODULE usrdef_zgr
