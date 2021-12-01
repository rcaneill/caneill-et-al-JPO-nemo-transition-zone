MODULE zgr_lib
   !!==============================================================================
   !!                       ***  MODULE usrlib_zgr   ***
   !! Ocean domain : library routines used to build a vertical coordinate 
   !!==============================================================================
   !! History :  4.0  ! 2019-05  ( R. Caneil, G. Madec)  Original code 
   !....
   !!             -   ! 2005-10  (A. Beckmann)  modifications for hybrid s-ccordinates & new stretching function
   !!            3.0  ! 2008-06  (G. Madec)  insertion of domzgr_zps.h90 & coding style
   !!            3.3  ! 2010-11  (G. Madec) add mbk. arrays associated to the deepest ocean level
   !!            3.4  ! 2012-08  (J. Siddorn) added Siddorn and Furner stretching function
   !!            3.6  ! 2014-11  (P. Mathiot and C. Harris) add ice shelf capabilitye  
   !!            3.?  ! 2015-11  (H. Liu) Modifications for Wetting/Drying
   !....
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zgr_bot_level: deepest ocean level for t-, u, and v-points
   !!   zgr_z        : reference z-coordinate 
   !!   zgr_zps      : z-coordinate with partial steps
   !!   zgr_sco      : s-coordinate
   !!   fssig        : tanh stretch function
   !!   fssig1       : Song and Haidvogel 1994 stretch function
   !!   fgamma       : Siddorn and Furner 2012 stretching function
   !!---------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zgr_zco      : z-coordinate
   !!
   !!   depth_to_e3    : use the depth of t- and w-points to calculate e3t & e3w
   !!                    (generic interface for 1D and 3D fields)
   !!   e3_to_depth    : use e3t & e3w to calculate the depth of t- and w-points
   !!                    (generic interface for 1D and 3D fields)
   !!  e3tw_to_other_e3: set e3u, e3v, e3f from e3t and e3uw, e3vw, e3w from e3w
   !!---------------------------------------------------------------------
   USE oce               ! ocean variables
   USE dom_oce           ! ocean domain
   !
   USE in_out_manager    ! I/O manager
   USE lbclnk            ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp           ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE
  
   INTERFACE depth_to_e3
      MODULE PROCEDURE depth_to_e3_1d, depth_to_e3_3d
   END INTERFACE

   INTERFACE e3_to_depth
      MODULE PROCEDURE e3_to_depth_1d, e3_to_depth_3d
   END INTERFACE

   PUBLIC   zgr_zco            ! called by usrdef_zgr
   PUBLIC   zgr_sco_pure       ! called by usrdef_zgr
   PUBLIC   zgr_sco_mi96       ! called by usrdef_zgr
   PUBLIC   depth_to_e3        ! called by usrdef_zgr
   PUBLIC   e3_to_depth        ! called by usrdef_zgr and domzgr.F90
      
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: depth_e3.F90 10069 2018-08-28 14:12:24Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS             

   SUBROUTINE zgr_zco( pdept_1d, pdepw_1d,                   &   ! <==>  1D reference vertical coordinate
      &                pe3t_1d , pe3w_1d ,                   &   ! ==>>  1D t & w-points vertical scale factors  
      &                pdept   , pdepw   ,                   &   ! ==>>  3D t & w-points depth
      &                pe3t    , pe3u    , pe3v    , pe3f,   &   ! ==>>  vertical scale factors
      &                pe3w    , pe3uw   , pe3vw           )     ! ==>>      -      -      -
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE zgr_zco  ***
      !!
      !! ** Purpose :   User defined the full step vertical coordinate  (zco)
      !!              from a 1D t- and w-points depth
      !!                recompute the depth as the sum of scale factor in 
      !!              order to preserve the truncation between gdepw and ht_0
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)    , INTENT(inout) ::   pdept_1d, pdepw_1d          ! 1D t- and w-depth          [m]
      REAL(wp), DIMENSION(:)    , INTENT(  out) ::   pe3t_1d, pe3w_1d            ! 1D t and w-scale factors   [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! 3D t and w-depth           [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors     [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         ! at t, u, v, f, w, uw, vw points
      !
      INTEGER  ::   jk   ! do-loop index
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'zgr_lib: zgr_zco : define full step vertical coord. system from 1D z-coordinate'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~'
      !
      DO jk = 1, jpk          !==  3D depth = 1D depth  ==!
         pdept(:,:,jk) = pdept_1d(jk) 
         pdepw(:,:,jk) = pdepw_1d(jk)
      END DO
      !                       !==  t- and w- scale factors from depth  ==!
      CALL depth_to_e3( pdept   , pdepw   , pe3t   , pe3w    )
      CALL depth_to_e3( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d )
      !
      !                       !==  update t- and w-depths from SUM(e3)  <== needed to get the same last digit in ht_0 calculation
      CALL e3_to_depth( pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d ) 
      CALL e3_to_depth( pe3t   ,  pe3w  , pdept   , pdepw    ) 
      !
      !                       !==  other scale factors  ==!
      pe3u (:,:,:) = pe3t(:,:,:)          ! t-level 
      pe3v (:,:,:) = pe3t(:,:,:)
      pe3f (:,:,:) = pe3t(:,:,:)
      pe3uw(:,:,:) = pe3w(:,:,:)          ! w-level 
      pe3vw(:,:,:) = pe3w(:,:,:)
      !
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) '              Reference 1D z-coordinate depth and scale factors:'
         WRITE(numout, "(9x,' level  gdept_1d  gdepw_1d  e3t_1d   e3w_1d  ')" )
         WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, pdept_1d(jk), pdepw_1d(jk), pe3t_1d(jk), pe3w_1d(jk), jk = 1, jpk )
      ENDIF
      !
   END SUBROUTINE zgr_zco


   SUBROUTINE zgr_sco_pure( pdept_1d, pdepw_1d,                   &   ! <==>  1D reference vertical coordinate
      &                     pbathy  , pHmax   ,                   &   ! <<==  bathymetry and its maximum
      &                     pe3t_1d , pe3w_1d ,                   &   ! ==>>  1D t & w-points vertical scale factors  
      &                     pdept   , pdepw   ,                   &   !       3D t & w-points depth
      &                     pe3t    , pe3u    , pe3v    , pe3f,   &   !       vertical scale factors
      &                     pe3w    , pe3uw   , pe3vw           )     !           -      -      -
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zgr_sco  ***
      !!                     
      !! ** Purpose :   
      !!
      !! ** Method  :   
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)    , INTENT(inout) ::   pdept_1d, pdepw_1d          ! 1D t- and w-depth          [m]
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   pbathy                      ! bathymetry                 [m]
      REAL(wp)                  , INTENT(in   ) ::   pHmax                       ! ocean depth maximum        [m]
      REAL(wp), DIMENSION(:)    , INTENT(  out) ::   pe3t_1d, pe3w_1d            ! 1D t and w-scale factors   [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! 3D t and w-depth           [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors     [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         ! at t, u, v, f, w, uw, vw points
      !
      INTEGER  ::   jk        ! do-loop index
      REAL(wp) ::   z1_Hmax   ! local scalar
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'zgr_lib: zgr_zco : define full step vertical coord. system from 1D z-coordinate'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~'
      !
      !                       !==  3D depth = 1D depth  ==!
      z1_Hmax = 1._wp / pHmax
      !
      DO jk = 1, jpk
         pdept(:,:,jk) = pdept_1d(jk) * pbathy(:,:) * z1_Hmax
         pdepw(:,:,jk) = pdepw_1d(jk) * pbathy(:,:) * z1_Hmax
      END DO
      !                       !==  t- and w- scale factors from depth  ==!
      CALL depth_to_e3( pdept   , pdepw   , pe3t   , pe3w    )
      CALL depth_to_e3( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d )
      !
      !                       !==  update t- and w-depths from SUM(e3)  <== needed to get the same last digit in ht_0 calculation
      CALL e3_to_depth( pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d ) 
      CALL e3_to_depth( pe3t   ,  pe3w  , pdept   , pdepw    ) 
      !
      !                       !==  other scale factors  ==!
      CALL e3tw_to_other_e3( pe3t , pe3w ,         &   ! <<==  e3 at t- and w-levels
         &                   pe3u , pe3v , pe3f,   &   ! ==>>  e3 at u- ,v- ,f-points
         &                   pe3uw, pe3vw        )     !         and uw-,vw-points
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) '              Reference 1D z-coordinate depth and scale factors:'
         WRITE(numout, "(9x,' level  gdept_1d  gdepw_1d  e3t_1d   e3w_1d  ')" )
         WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, pdept_1d(jk), pdepw_1d(jk), pe3t_1d(jk), pe3w_1d(jk), jk = 1, jpk )
      ENDIF
      !
      !
   END SUBROUTINE zgr_sco_pure


   SUBROUTINE zgr_sco_mi96( pdept_1d, pdepw_1d,                   &   ! ==>>  1D reference vertical coordinate
      &                     pbathy  ,                             &   ! <<==  bathymetry
      &                     pdzmin  , pkth    , kkconst , pacr,   &   ! <<==  parameters for the tanh stretching
      &                     pe3t_1d , pe3w_1d ,                   &   ! ==>>  1D t & w-points vertical scale factors  
      &                     pdept   , pdepw   ,                   &   !       3D t & w-points depth
      &                     pe3t    , pe3u    , pe3v    , pe3f,   &   !       vertical scale factors
      &                     pe3w    , pe3uw   , pe3vw           )     !           -      -      -
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zgr_sco  ***
      !!                     
      !! ** Purpose :   TODO
      !!
      !! ** Method  :   TODO
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)    , INTENT(  out) ::   pdept_1d, pdepw_1d          ! 1D t- and w-depth                    [m]
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   pbathy                      ! bathymetry                           [m]
      REAL(wp)                  , INTENT(in   ) ::   pdzmin                      ! minimum value of e3 at the surface   [m]
      REAL(wp)                  , INTENT(in   ) ::   pkth                        ! position of the inflexion point
      INTEGER                   , INTENT(in   ) ::   kkconst                     ! Number of levels with pure z-coordinate (i.e. constant depth)
      REAL(wp)                  , INTENT(in   ) ::   pacr                        ! slope of the tanh
      REAL(wp), DIMENSION(:)    , INTENT(  out) ::   pe3t_1d, pe3w_1d            ! 1D t and w-scale factors (computed at 4000m)  [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! 3D t and w-depth           [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors     [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         ! at t, u, v, f, w, uw, vw points
      !
      INTEGER  ::   ji, jj, jk        ! do-loop index
      REAL(wp) ::   z1_Hmax   ! local scalar
      REAL(wp), DIMENSION(2,jpk)   ::   zlev_dep   ! depth of the levels
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'zgr_lib: zgr_sco : define full vertical s-coord. system using the 2d bathymetry field'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~'
      !
      DO jj = 1, jpj
         DO ji = 1, jpi
            zlev_dep = mi96_1d( pbathy(ji,jj), pdzmin, pkth, kkconst, pacr )
            pdepw(ji,jj,:) = zlev_dep(1,:)
            pdept(ji,jj,:) = zlev_dep(2,:)
         END DO
      END DO
      !
      ! 1D profile
      zlev_dep = mi96_1d( 4000._wp, pdzmin, pkth, kkconst, pacr, ldprint=.TRUE. )
      pdepw_1d(:) = zlev_dep(1,:)
      pdept_1d(:) = zlev_dep(2,:)
      !
      !                       !==  t- and w- scale factors from depth  ==!
      CALL depth_to_e3( pdept   , pdepw   , pe3t   , pe3w    )
      CALL depth_to_e3( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d )
      !
      !                       !==  update t- and w-depths from SUM(e3)  <== needed to get the same last digit in ht_0 calculation
      CALL e3_to_depth( pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d ) 
      CALL e3_to_depth( pe3t   ,  pe3w  , pdept   , pdepw    ) 
      !
      !                       !==  other scale factors  ==!
      CALL e3tw_to_other_e3( pe3t , pe3w ,         &   ! <<==  e3 at t- and w-levels
         &                   pe3u , pe3v , pe3f,   &   ! ==>>  e3 at u- ,v- ,f-points
         &                   pe3uw, pe3vw        )     !         and uw-,vw-points
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) '              Reference 1D z-coordinate depth and scale factors:'
         WRITE(numout, "(9x,' level  gdept_1d  gdepw_1d  e3t_1d   e3w_1d  ')" )
         WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, pdept_1d(jk), pdepw_1d(jk), pe3t_1d(jk), pe3w_1d(jk), jk = 1, jpk )
      ENDIF
      !
      !
   END SUBROUTINE zgr_sco_mi96 


   SUBROUTINE sco_sh94( pdept_1d, pdepw_1d,                   &   ! <==>  1D reference vertical coordinate
      &                 pbathy  ,                             &   ! <<==  bathymetry
      &                 pn_theta, pn_bb   , pn_hc   ,         &   !       and sh94 stretching parameters
      &                 pe3t_1d , pe3w_1d ,                   &   ! ==>>  1D t & w-points vertical scale factors  
      &                 pdept   , pdepw   ,                   &   !       3D t & w-points depth
      &                 pe3t    , pe3u    , pe3v    , pe3f,   &   !       vertical scale factors
      &                 pe3w    , pe3uw   , pe3vw           )     !           -      -      -
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sco_sh94  ***
      !!                     
      !! ** Purpose :   stretch the s-coordinate system
      !!
      !! ** Method  :   s-coordinate stretch using the Song and Haidvogel 1994
      !!                mixed S/sigma coordinate
      !!
      !! Reference : Song and Haidvogel 1994. 
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)    , INTENT(inout) ::   pdept_1d, pdepw_1d          ! 1D t- and w-depth         [m]
      REAL(wp), DIMENSION(:,:)  , INTENT(in   ) ::   pbathy                      ! ocean bathymetry          [m]
      REAL(wp)                  , INTENT(in   ) ::   pn_theta, pn_bb, pn_hc      ! sh94 strecting parameters [-]
      REAL(wp), DIMENSION(:)    , INTENT(  out) ::   pe3t_1d, pe3w_1d            ! 1D t and w-scale factors  [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! 3D t and w-depth          [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors    [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         ! at t, u, v, f, w, uw, vw points
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop argument
      REAL(wp) ::   zct, zcw, zcoeft, zcoefw     ! local scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zs_w, zs_t   ! 3D workspace
      REAL(wp), DIMENSION(jpi,jpj)     ::   hbatt, scosrf   !!rc TODO put in argument
      !!----------------------------------------------------------------------
      WRITE(numout,*) '   s_sh94:   Song and Haidvogel 1994 s-coordinate'
      WRITE(numout,*) '   ~~~~~~~'
      WRITE(numout,*) '      Critical depth (h<hc : pure sigma)           rn_hc         = ', pn_hc
      WRITE(numout,*) '      surface control parameter (0<=rn_theta<=20)  rn_theta      = ', pn_theta
      WRITE(numout,*) '      stretching parameter (song and haidvogel)    rn_bb         = ', pn_bb
      !
      zs_w(:,:,:) = 0._wp
      zs_t(:,:,:) = 0._wp
      !
      DO ji = 1, jpi
         DO jj = 1, jpj
            IF( pbathy(ji,jj) <= pn_hc ) THEN   !-  shallow water, pure sigma-coordinate (stretching in z/h)
               DO jk = 1, jpk
                  zs_w(ji,jj,jk) =   REAL(jk-1,wp)            / REAL(jpkm1,wp)
                  zs_t(ji,jj,jk) = ( REAL(jk-1,wp) + 0.5_wp ) / REAL(jpkm1,wp)
               END DO
            ELSE                                !-  deep water, st
               DO jk = 1, jpk
                  zct = ( REAL(jk,wp) - 0.5_wp ) / REAL(jpkm1,wp)
                  zcw = ( REAL(jk,wp) - 1.0_wp ) / REAL(jpkm1,wp)
                  pdept(ji,jj,jk) = ( scosrf(ji,jj) + (hbatt(ji,jj)-pn_hc)*zs_t(ji,jj,jk)+pn_hc*zcoeft )
                  pdepw(ji,jj,jk) = ( scosrf(ji,jj) + (hbatt(ji,jj)-pn_hc)*zs_w(ji,jj,jk)+pn_hc*zcoefw )
               END DO
            ENDIF
            !
         END DO
      END DO
      !
      !                       !==  t- and w- scale factors from depth  ==!
      CALL depth_to_e3( pdept   , pdepw   , pe3t   , pe3w    )
      CALL depth_to_e3( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d )
      !
      !                       !==  update t- and w-depths from SUM(e3)  <== needed to get the same last digit in ht_0 calculation
      CALL e3_to_depth( pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d ) 
      CALL e3_to_depth( pe3t   , pe3w   , pdept   , pdepw    ) 
      !
      !                       !==  other scale factors  ==!
      CALL e3tw_to_other_e3( pe3t , pe3t ,           &   ! <<==  in
         &                   pe3u , pe3v , pe3f,     &   ! ==>>  out
         &                   pe3uw, pe3vw          )
      !
   END SUBROUTINE sco_sh94


   SUBROUTINE e3tw_to_other_e3( pe3t , pe3w ,         &   ! <<==  e3 at t- and w-levels
      &                         pe3u , pe3v , pe3f,   &   ! ==>>  e3 at u- ,v- ,f-points
      &                         pe3uw, pe3vw        )     !         and uw-,vw-points
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE e3tw_to_other_e3  ***
      !!                     
      !! ** Purpose :   set the vertical scale factors used in the computation
      !!
      !! ** Method  :   s-coordinate : scale factors are simple t-level (w-level) 
      !!              averaging from e3t (e3w) neighbours.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pe3t, pe3w           ! t- and w-scale factors   [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3u , pe3v , pe3f   ! vertical scale factors    [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3uw, pe3vw         ! at t, u, v, f, w, uw, vw points
      !
      INTEGER  ::   ji, jj, jk           ! dummy loop argument
      !!----------------------------------------------------------------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '   e3tw_to_other_e3:   set e3u ,e3v , e3f from e3t'
      IF(lwp) WRITE(numout,*) '   ~~~~~~~~~~~~~~~~~   and e3uw,e3vw, e3w from e3w'
      !
      !                       !* set values except on last row or/and column
      DO ji = 1, jpim1
         pe3u (ji,:,:) = 0.50_wp * ( pe3t(ji,:,:) + pe3t(ji+1,:,:) )   
         pe3uw(ji,:,:) = 0.50_wp * ( pe3w(ji,:,:) + pe3w(ji+1,:,:) )
      END DO      
      DO jj = 1, jpjm1
         pe3v (:,jj,:) = 0.50_wp * ( pe3t(:,jj,:) + pe3t(:,jj+1,:) )
         pe3vw(:,jj,:) = 0.50_wp * ( pe3w(:,jj,:) + pe3w(:,jj+1,:) )
      END DO
      DO jk = 1, jpk
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               pe3f (ji,jj,jk) = 0.25_wp * (  pe3t(ji  ,jj,jk) + pe3t(ji  ,jj+1,jk)  &
                  &                         + pe3t(ji+1,jj,jk) + pe3t(ji+1,jj+1,jk)  )
            END DO
         END DO
      END DO
      ! 
      !                       !*  Apply lateral boundary condition   (use of optional argument cd_mpp)
      !                           !   MPI-halos exchange + applied lbc with specific treatment of closed boundaries :
      !                           !     - North and East land boundaries set to inner values
      !                           !     - South and West land boundaries not updated (they have already a good value
!!rc      CALL lbc_lnk_multi( 'e3tw_to_other_e3', pe3u , 'U', 1., pe3v , 'V', 1., pe3f , 'F', 1.,     &
!!rc         &                                    pe3uw, 'U', 1., pe3uw, 'V', 1., cd_mpp = 'land_halo_set_to_inner_values_NE')
      CALL lbc_lnk_multi( 'e3tw_to_other_e3', pe3u , 'U', 1., pe3v , 'V', 1., pe3f , 'F', 1.,     &
          &                                    pe3uw, 'U', 1., pe3uw, 'V', 1., cd_mpp = 'only_communications' )
      IF( mig(jpi) == jpiglo ) THEN   ! last column of the local domain is the global domain one (which is closed here)
         !                            ! Extend inner domain value on the last column
         pe3u (jpi, :     ,:) = pe3u (jpim1, :     ,:)
         pe3uw(jpi, :     ,:) = pe3uw(jpim1, :     ,:)
         pe3f (jpi,1:jpjm1,:) = pe3f (jpim1,1:jpjm1,:)
      ENDIF
      !
      IF( mjg(jpj) == jpjglo ) THEN   ! last row of the local domain is the global domain one (which is closed here)
         !                            ! Extend inner domain value on the last row
         pe3v (:,jpj,:) = pe3v (:,jpjm1,:)
         pe3vw(:,jpj,:) = pe3vw(:,jpjm1,:)
         pe3f (:,jpj,:) = pe3f (:,jpjm1,:)
      ENDIF
      
      !
   END SUBROUTINE e3tw_to_other_e3


   FUNCTION fssig1( pk1, pn_theta, pbb ) RESULT( pf1 )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE fssig1 ***
      !!
      !! ** Purpose :   provide the Song and Haidvogel version of the analytical function in s-coordinate
      !!
      !! ** Method  :   the function provides the non-dimensional position of
      !!                T and W (i.e. between 0 and 1)
      !!                T-points at integer values (between 1 and jpk)
      !!                W-points at integer values - 1/2 (between 0.5 and jpk-0.5)
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::   pk1        ! continuous "k" coordinate
      REAL(wp), INTENT(in) ::   pn_theta   ! Stretching coefficient
      REAL(wp), INTENT(in) ::   pbb        ! Stretching coefficient
      REAL(wp)             ::   pf1        ! sigma value
      !!----------------------------------------------------------------------
      !
      IF( pn_theta == 0 ) THEN     !==  uniform sigma  ==!
         pf1 = - ( pk1 - 0.5_wp ) / REAL( jpkm1 )
      ELSEIF( pbb == 0  ) THEN     !==  surface increased stretching  ==! 
         pf1 =                       SINH( -pn_theta*(pk1-0.5_wp)/REAL(jpkm1)          ) / SINH( pn_theta )
      ELSE                         !==  surface and bottom increased stretching
         pf1 =   ( 1._wp - pbb ) *   SINH( -pn_theta*(pk1-0.5_wp)/REAL(jpkm1)          ) / SINH( pn_theta )               &
            &  +           pbb   * ( TANH(  pn_theta*(pk1-0.5_wp)/REAL(jpkm1) + 0.5_wp ) - TANH( 0.5_wp * pn_theta )  )   &
            &                    / ( 2._wp * TANH( 0.5_wp * pn_theta ) )
      ENDIF
      !
   END FUNCTION fssig1

   FUNCTION mi96_1d( pHmax, pdzmin, pkth, kkconst, pacr, ldprint ) RESULT( zlev_dep )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE  ***
      !!
      !! ** Purpose :   
      !!
      !! ** Method  :   
      !!                
      !!                
      !!                
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in   )             ::   pHmax   ! depth of the bottom   [meters]
      REAL(wp), INTENT(in   )             ::   pdzmin   ! minimum value of e3 at the surface   [meters]
      REAL(wp), INTENT(in   )             ::   pkth   ! position of the inflexion point
      INTEGER , INTENT(in   )             ::   kkconst   ! Number of levels with pure z-coordinate (i.e. constant depth)
      REAL(wp), INTENT(in   )             ::   pacr   ! slope of the tanh
      LOGICAL , INTENT(in   ), OPTIONAL   ::   ldprint   ! Whether to print infos about the arguments (default .false.)
      REAL(wp), DIMENSION(2,jpk)          ::   zlev_dep   ! depth of the levels
      !!----------------------------------------------------------------------
      INTEGER  ::   jk
      REAL(wp) ::   za0, za1, zsur   ! Computed parameters for the tanh
      REAL(wp) ::   zh_co   ! depth of the connection between pure z-coordinate and s-coordinates   [meters]
      REAL(wp) ::   zw, zt   ! local scalar
      !
      zh_co = REAL( kkconst , wp ) * pdzmin
      !
      za1  = (  pdzmin - (pHmax - zh_co) / REAL(jpkm1-kkconst,wp)  )                                                      &
         & / ( TANH((1-pkth)/pacr) - pacr / REAL(jpkm1-kkconst) * (  LOG( COSH( (jpk - kkconst - pkth) / pacr) )      &
         &                                                         - LOG( COSH( ( 1            - pkth) / pacr) )  )  )
      za0  = pdzmin  - za1 *             TANH( (1-pkth) / pacr )
      zsur =   - za0 - za1 * pacr * LOG( COSH( (1-pkth) / pacr )  )
      !
      IF( lwp .AND. PRESENT(ldprint) ) THEN                         ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) '    mi96_id : Vertical coordinate ... TODO'
         WRITE(numout,*) '    ~~~~~~~'
         WRITE(numout,*) '       Depth of the bottom                                   pHmax = ', pHmax
         WRITE(numout,*) '       Value of e3 in the pure z-coordinate'
         WRITE(numout,*) '          (= value of e3 at the first s-coordinate level)   pdzmin = ', pdzmin
         WRITE(numout,*) '       Position of the inflexion point                        pkth = ', pkth
         WRITE(numout,*) '       Number of pure z levels                             kkconst = ', kkconst
         WRITE(numout,*) '       Slope of the tanh                                      pacr = ', pacr
      ENDIF
      !
      DO jk = 1, kkconst      ! z-coordinate part: depth at T and W-points
         zw = REAL( jk , wp )
         zt = REAL( jk , wp ) + 0.5_wp
         !
         zlev_dep(1,jk) = (zw - 1._wp) * pdzmin
         zlev_dep(2,jk) = (zt - 1._wp) * pdzmin
      END DO
      DO jk = kkconst+1, jpk          ! s-coordinate part: depth at T and W-points
         zw = REAL( jk - kkconst, wp )
         zt = REAL( jk - kkconst, wp ) + 0.5_wp
         !
         zlev_dep(1,jk) = zsur + za0 * zw + za1 * pacr *  LOG( COSH( (zw-pkth) / pacr ) ) + zh_co
         zlev_dep(2,jk) = zsur + za0 * zt + za1 * pacr *  LOG( COSH( (zt-pkth) / pacr ) ) + zh_co
      END DO
      !
   END FUNCTION mi96_1d





   SUBROUTINE zgr_bat_env    !    ( pbat, pbat_env )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zgr_sco  ***
      !!                     
      !! ** Purpose :   define an envelop 
      !!
      !! ** Method  :   s-coordinate
      !!----------------------------------------------------------------------
      
      !!
      REAL(wp) ::   rn_rmax, rn_sbot_max, rn_sbot_min    !!rc TODO argument
      !
      INTEGER  ::   ji, jj, jk, jl           ! dummy loop argument
      INTEGER  ::   iip1, ijp1, iim1, ijm1   ! temporary integers
      INTEGER  ::   ios                      ! Local integer output status for namelist read
      REAL(wp) ::   zrmax, ztaper   ! temporary scalars
      REAL(wp) ::   zrfact
      !
      REAL(wp), DIMENSION(jpi,jpj) :: ztmpi1, ztmpi2, ztmpj1, ztmpj2
      REAL(wp), DIMENSION(jpi,jpj) :: zenv, ztmp, zmsk, zri, zrj, zhbat
      REAL(wp), DIMENSION(jpi,jpj) :: bathy, hbatt, hbatu, hbatv, hbatf
      REAL(wp), DIMENSION(jpi,jpj) :: hift, hifu, hifv, hiff, scosrf, scobot
      !!
      !!----------------------------------------------------------------------
      !

      hift(:,:) = rn_sbot_min                     ! set the minimum depth for the s-coordinate
      hifu(:,:) = rn_sbot_min
      hifv(:,:) = rn_sbot_min
      hiff(:,:) = rn_sbot_min

      !                                        ! set maximum ocean depth
      bathy(:,:) = MIN( rn_sbot_max, bathy(:,:) )

         DO jj = 1, jpj
            DO ji = 1, jpi
              IF( bathy(ji,jj) > 0._wp )   bathy(ji,jj) = MAX( rn_sbot_min, bathy(ji,jj) )
            END DO
         END DO
      !                                        ! =============================
      !                                        ! Define the envelop bathymetry   (hbatt)
      !                                        ! =============================
      ! use r-value to create hybrid coordinates
      zenv(:,:) = bathy(:,:)
      !
      ! set first land point adjacent to a wet cell to sbot_min as this needs to be included in smoothing
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( bathy(ji,jj) == 0._wp ) THEN
                  iip1 = MIN( ji+1, jpi )
                  ijp1 = MIN( jj+1, jpj )
                  iim1 = MAX( ji-1, 1 )
                  ijm1 = MAX( jj-1, 1 )
!!gm BUG fix see ticket #1617
                  IF( ( + bathy(iim1,ijm1) + bathy(ji,ijp1) + bathy(iip1,ijp1)              &
                     &  + bathy(iim1,jj  )                  + bathy(iip1,jj  )              &
                     &  + bathy(iim1,ijm1) + bathy(ji,ijm1) + bathy(iip1,ijp1)  ) > 0._wp ) &
                     &    zenv(ji,jj) = rn_sbot_min
!!gm
!!gm               IF( ( bathy(iip1,jj  ) + bathy(iim1,jj  ) + bathy(ji,ijp1  ) + bathy(ji,ijm1) +         &
!!gm                  &  bathy(iip1,ijp1) + bathy(iim1,ijm1) + bathy(iip1,ijp1) + bathy(iim1,ijm1)) > 0._wp ) THEN
!!gm               zenv(ji,jj) = rn_sbot_min
!!gm             ENDIF
!!gm end
              ENDIF
            END DO
         END DO

      ! apply lateral boundary condition   CAUTION: keep the value when the lbc field is zero
      CALL lbc_lnk( 'zgr_bat_env', zenv, 'T', 1._wp, 'no0' )
      ! 
      ! smooth the bathymetry (if required)
      scosrf(:,:) = 0._wp             ! ocean surface depth (here zero: no under ice-shelf sea)
      scobot(:,:) = bathy(:,:)        ! ocean bottom  depth
      !
      jl = 0
      zrmax = 1._wp
      !   
      !     
      ! set scaling factor used in reducing vertical gradients
      zrfact = ( 1._wp - rn_rmax ) / ( 1._wp + rn_rmax )
      !
      ! initialise temporary evelope depth arrays
      ztmpi1(:,:) = zenv(:,:)
      ztmpi2(:,:) = zenv(:,:)
      ztmpj1(:,:) = zenv(:,:)
      ztmpj2(:,:) = zenv(:,:)
      !
      ! initialise temporary r-value arrays
      zri(:,:) = 1._wp
      zrj(:,:) = 1._wp
      !                                                            ! ================ !
      DO WHILE( jl <= 10000 .AND. ( zrmax - rn_rmax ) > 1.e-8_wp ) !  Iterative loop  !
         !                                                         ! ================ !
         jl = jl + 1
         zrmax = 0._wp
         ! we set zrmax from previous r-values (zri and zrj) first
         ! if set after current r-value calculation (as previously)
         ! we could exit DO WHILE prematurely before checking r-value
         ! of current zenv
         DO jj = 1, nlcj
            DO ji = 1, nlci
               zrmax = MAX( zrmax, ABS(zri(ji,jj)), ABS(zrj(ji,jj)) )
            END DO
         END DO
         zri(:,:) = 0._wp
         zrj(:,:) = 0._wp
         DO jj = 1, nlcj
            DO ji = 1, nlci
               iip1 = MIN( ji+1, nlci )      ! force zri = 0 on last line (ji=ncli+1 to jpi)
               ijp1 = MIN( jj+1, nlcj )      ! force zrj = 0 on last raw  (jj=nclj+1 to jpj)
               IF( (zenv(ji,jj) > 0._wp) .AND. (zenv(iip1,jj) > 0._wp)) THEN
                  zri(ji,jj) = ( zenv(iip1,jj  ) - zenv(ji,jj) ) / ( zenv(iip1,jj  ) + zenv(ji,jj) )
               END IF
               IF( (zenv(ji,jj) > 0._wp) .AND. (zenv(ji,ijp1) > 0._wp)) THEN
                  zrj(ji,jj) = ( zenv(ji  ,ijp1) - zenv(ji,jj) ) / ( zenv(ji  ,ijp1) + zenv(ji,jj) )
               END IF
               IF( zri(ji,jj) >  rn_rmax )   ztmpi1(ji  ,jj  ) = zenv(iip1,jj  ) * zrfact
               IF( zri(ji,jj) < -rn_rmax )   ztmpi2(iip1,jj  ) = zenv(ji  ,jj  ) * zrfact
               IF( zrj(ji,jj) >  rn_rmax )   ztmpj1(ji  ,jj  ) = zenv(ji  ,ijp1) * zrfact
               IF( zrj(ji,jj) < -rn_rmax )   ztmpj2(ji  ,ijp1) = zenv(ji  ,jj  ) * zrfact
            END DO
         END DO
         IF( lk_mpp )   CALL mpp_max( 'zgr_bat_env', zrmax )   ! max over the global domain
         !
         IF(lwp)WRITE(numout,*) 'zgr_sco :   iter= ',jl, ' rmax= ', zrmax
         !
         DO jj = 1, nlcj
            DO ji = 1, nlci
               zenv(ji,jj) = MAX(zenv(ji,jj), ztmpi1(ji,jj), ztmpi2(ji,jj), ztmpj1(ji,jj), ztmpj2(ji,jj) )
            END DO
         END DO
         ! apply lateral boundary condition   CAUTION: keep the value when the lbc field is zero
         CALL lbc_lnk( 'zgr_bat_env',  zenv, 'T', 1._wp, 'no0' )
         !                                                  ! ================ !
      END DO                                                !     End loop     !
      !                                                     ! ================ !
      DO jj = 1, jpj
         DO ji = 1, jpi
            zenv(ji,jj) = MAX( zenv(ji,jj), rn_sbot_min ) ! set all points to avoid undefined scale value warnings
         END DO
      END DO
      !
      ! Envelope bathymetry saved in hbatt
      hbatt(:,:) = zenv(:,:) 
      IF( MINVAL( gphit(:,:) ) * MAXVAL( gphit(:,:) ) <= 0._wp ) THEN
         CALL ctl_warn( ' s-coordinates are tapered in vicinity of the Equator' )
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztaper = EXP( -(gphit(ji,jj)/8._wp)**2._wp )
               hbatt(ji,jj) = rn_sbot_max * ztaper + hbatt(ji,jj) * ( 1._wp - ztaper )
            END DO
         END DO
      ENDIF
      !
      
            !
      !
      !                                        ! ==============================
      !                                        !   hbatu, hbatv, hbatf fields
      !                                        ! ==============================
      IF(lwp) THEN
         WRITE(numout,*)
           WRITE(numout,*) ' zgr_sco: minimum depth of the envelop topography set to : ', rn_sbot_min
      ENDIF
      hbatu(:,:) = rn_sbot_min
      hbatv(:,:) = rn_sbot_min
      hbatf(:,:) = rn_sbot_min
      DO jj = 1, jpjm1
        DO ji = 1, jpim1   ! NO vector opt.
           hbatu(ji,jj) = 0.50_wp * ( hbatt(ji  ,jj) + hbatt(ji+1,jj  ) )
           hbatv(ji,jj) = 0.50_wp * ( hbatt(ji  ,jj) + hbatt(ji  ,jj+1) )
           hbatf(ji,jj) = 0.25_wp * ( hbatt(ji  ,jj) + hbatt(ji  ,jj+1)   &
              &                     + hbatt(ji+1,jj) + hbatt(ji+1,jj+1) )
        END DO
      END DO

      ! 
      ! Apply lateral boundary condition
!!gm  ! CAUTION: retain non zero value in the initial file this should be OK for orca cfg, not for EEL
      zhbat(:,:) = hbatu(:,:)   ;   CALL lbc_lnk( 'zgr_bat_env',  hbatu, 'U', 1._wp )
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( hbatu(ji,jj) == 0._wp ) THEN
               !No worries about the following line when ln_wd == .true.
               IF( zhbat(ji,jj) == 0._wp )   hbatu(ji,jj) = rn_sbot_min
               IF( zhbat(ji,jj) /= 0._wp )   hbatu(ji,jj) = zhbat(ji,jj)
            ENDIF
         END DO
      END DO
      zhbat(:,:) = hbatv(:,:)   ;   CALL lbc_lnk( 'zgr_bat_env',  hbatv, 'V', 1._wp )
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( hbatv(ji,jj) == 0._wp ) THEN
               IF( zhbat(ji,jj) == 0._wp )   hbatv(ji,jj) = rn_sbot_min
               IF( zhbat(ji,jj) /= 0._wp )   hbatv(ji,jj) = zhbat(ji,jj)
            ENDIF
         END DO
      END DO
      zhbat(:,:) = hbatf(:,:)   ;   CALL lbc_lnk( 'zgr_bat_env',  hbatf, 'F', 1._wp )
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( hbatf(ji,jj) == 0._wp ) THEN
               IF( zhbat(ji,jj) == 0._wp )   hbatf(ji,jj) = rn_sbot_min
               IF( zhbat(ji,jj) /= 0._wp )   hbatf(ji,jj) = zhbat(ji,jj)
            ENDIF
         END DO
      END DO

!!bug:  key_helsinki a verifer
        hift(:,:) = MIN( hift(:,:), hbatt(:,:) )
        hifu(:,:) = MIN( hifu(:,:), hbatu(:,:) )
        hifv(:,:) = MIN( hifv(:,:), hbatv(:,:) )
        hiff(:,:) = MIN( hiff(:,:), hbatf(:,:) )

      IF( nprint == 1 .AND. lwp )   THEN
         WRITE(numout,*) ' MAX val hif   t ', MAXVAL( hift (:,:) ), ' f ', MAXVAL( hiff (:,:) ),  &
            &                        ' u ',   MAXVAL( hifu (:,:) ), ' v ', MAXVAL( hifv (:,:) )
         WRITE(numout,*) ' MIN val hif   t ', MINVAL( hift (:,:) ), ' f ', MINVAL( hiff (:,:) ),  &
            &                        ' u ',   MINVAL( hifu (:,:) ), ' v ', MINVAL( hifv (:,:) )
         WRITE(numout,*) ' MAX val hbat  t ', MAXVAL( hbatt(:,:) ), ' f ', MAXVAL( hbatf(:,:) ),  &
            &                        ' u ',   MAXVAL( hbatu(:,:) ), ' v ', MAXVAL( hbatv(:,:) )
         WRITE(numout,*) ' MIN val hbat  t ', MINVAL( hbatt(:,:) ), ' f ', MINVAL( hbatf(:,:) ),  &
            &                        ' u ',   MINVAL( hbatu(:,:) ), ' v ', MINVAL( hbatv(:,:) )
      ENDIF
!! helsinki
      !
   END SUBROUTINE zgr_bat_env
   
   
   SUBROUTINE depth_to_e3_1d( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE depth_to_e3_1d  ***
      !!
      !! ** Purpose :   compute e3t & e3w scale factors from t- & w-depths of model levels
      !!
      !! ** Method  :   The scale factors are given by the discrete derivative 
      !!              of the depth:
      !!                               e3w(jk) = dk[ dept_1d ] 
      !!                               e3t(jk) = dk[ depw_1d ]
      !!              with, at top and bottom :
      !!                      e3w( 1 ) = 2 * ( dept( 1 ) - depw( 1 ) )
      !!                      e3t(jpk) = 2 * ( dept(jpk) - depw(jpk) )   
      !!
      !! ** Action  : - pe3t_1d , pe3w_1d  : scale factors at T- and W-levels (m)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:), INTENT(in   ) ::   pdept_1d, pdepw_1d   ! depths          [m]
      REAL(wp), DIMENSION(:), INTENT(  out) ::   pe3t_1d , pe3w_1d    ! e3.=dk[depth]   [m]
      !
      INTEGER  ::   jk           ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      ! use pdep. at w- and t-points to compute e3. (e3. = dk[depth])
      !
      pe3w_1d( 1 ) = 2._wp * ( pdept_1d(1) - pdepw_1d(1) ) 
      DO jk = 1, jpkm1
         pe3w_1d(jk+1) = pdept_1d(jk+1) - pdept_1d(jk) 
         pe3t_1d(jk  ) = pdepw_1d(jk+1) - pdepw_1d(jk) 
      END DO
      pe3t_1d(jpk) = 2._wp * ( pdept_1d(jpk) - pdepw_1d(jpk) )
      !
   END SUBROUTINE depth_to_e3_1d
   
      
   SUBROUTINE depth_to_e3_3d( pdept_3d, pdepw_3d, pe3t_3d, pe3w_3d )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE depth_to_e3_3d  ***
      !!
      !! ** Purpose :   compute e3t & e3w scale factors from t- & w-depths of model levels
      !!
      !! ** Method  :   The scale factors are given by the discrete derivative 
      !!              of the depth:
      !!                               e3w(jk) = dk[ dept_1d ] 
      !!                               e3t(jk) = dk[ depw_1d ]
      !!              with, at top and bottom :
      !!                      e3w( 1 ) = 2 * ( dept( 1 ) - depw( 1 ) )
      !!                      e3t(jpk) = 2 * ( dept(jpk) - depw(jpk) )   
      !!
      !! ** Action  : - pe3t_1d , pe3w_1d  : scale factors at T- and W-levels (m)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pdept_3d, pdepw_3d   ! depth           [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t_3d , pe3w_3d    ! e3.=dk[depth]   [m]
      !
      INTEGER  ::   jk           ! dummy loop indices
      !!----------------------------------------------------------------------      
      pe3w_3d(:,:, 1 ) = 2._wp * ( pdept_3d(:,:,1) - pdepw_3d(:,:,1) ) 
      DO jk = 1, jpkm1
         pe3w_3d(:,:,jk+1) = pdept_3d(:,:,jk+1) - pdept_3d(:,:,jk) 
         pe3t_3d(:,:,jk  ) = pdepw_3d(:,:,jk+1) - pdepw_3d(:,:,jk) 
      END DO
      pe3t_3d(:,:,jpk) = 2._wp * ( pdept_3d(:,:,jpk) - pdepw_3d(:,:,jpk) )   
      !
   END SUBROUTINE depth_to_e3_3d


   SUBROUTINE e3_to_depth_1d( pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE e3_to_depth_1d  ***
      !!
      !! ** Purpose :   compute t- & w-depths of model levels from e3t & e3w scale factors
      !!
      !! ** Method  :   The t- & w-depth are given by the summation of e3w & e3t, resp. 
      !!
      !! ** Action  : - pe3t_1d, pe3w_1d : scale factor of t- and w-point (m)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:), INTENT(in   ) ::   pe3t_1d , pe3w_1d    ! vert. scale factors   [m]
      REAL(wp), DIMENSION(:), INTENT(  out) ::   pdept_1d, pdepw_1d   ! depth = SUM( e3 )     [m]
      !
      INTEGER  ::   jk           ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      pdepw_1d(1) = 0.0_wp
      pdept_1d(1) = 0.5_wp * pe3w_1d(1)
      DO jk = 2, jpk
         pdepw_1d(jk) = pdepw_1d(jk-1) + pe3t_1d(jk-1) 
         pdept_1d(jk) = pdept_1d(jk-1) + pe3w_1d(jk  ) 
      END DO
      !
   END SUBROUTINE e3_to_depth_1d
   
      
   SUBROUTINE e3_to_depth_3d( pe3t_3d, pe3w_3d, pdept_3d, pdepw_3d )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE e3_to_depth_3d  ***
      !!
      !! ** Purpose :   compute t- & w-depths of model levels from e3t & e3w scale factors
      !!
      !! ** Method  :   The t- & w-depth are given by the summation of e3w & e3t, resp. 
      !!
      !! ** Action  : - pe3t_1d, pe3w_1d : scale factor of t- and w-point (m)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(in   ) ::   pe3t_3d , pe3w_3d    ! vert. scale factors   [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept_3d, pdepw_3d   ! depth = SUM( e3 )     [m]
      !
      INTEGER  ::   jk           ! dummy loop indices
      !!----------------------------------------------------------------------      
      !
      pdepw_3d(:,:,1) = 0.0_wp
      pdept_3d(:,:,1) = 0.5_wp * pe3w_3d(:,:,1)
      DO jk = 2, jpk
         pdepw_3d(:,:,jk) = pdepw_3d(:,:,jk-1) + pe3t_3d(:,:,jk-1) 
         pdept_3d(:,:,jk) = pdept_3d(:,:,jk-1) + pe3w_3d(:,:,jk  ) 
      END DO
      !
   END SUBROUTINE e3_to_depth_3d

   !!======================================================================
END MODULE zgr_lib
