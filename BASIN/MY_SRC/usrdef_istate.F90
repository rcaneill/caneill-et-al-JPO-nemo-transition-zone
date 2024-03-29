MODULE usrdef_istate
   !!======================================================================
   !!                     ***  MODULE usrdef_istate   ***
   !!
   !!                      ===  CANAL configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2017-11  (J. Chanut) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce, ONLY:         ! ocean space and time domain
   USE dom_oce, ONLY: gphit        
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   !   
   USE usrdef_nam, ONLY: nn_initcase
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate   ! called by istate.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_istate.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here CANAL configuration 
      !!
      !! ** Method  :   Set a gaussian anomaly of pressure and associated
      !!                geostrophic velocities
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height
      !
      INTEGER  ::   jk, jj, ji   ! dummy loop
      REAL(wp) ::   zTtop, zTbot, zphiMAX, z1_phiMAX   ! local scalar
      REAL(wp) ::   z1_700, z1_150, z1_1500, z7_1500, z1_100, z1_460, z1_5000, z1_650
      REAL(wp), DIMENSION(2) ::   zmpparr
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : CANAL configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   '
      !
      SELECT CASE(nn_initcase)
      CASE(0)    ! rest
         ! sea level:
         pssh(:,:) = 0._wp
         ! temperature:
         pts(:,:,:,jp_tem) = 10._wp * ptmask(:,:,:)
         ! salinity:  
         pts(:,:,:,jp_sal) = 35._wp * ptmask(:,:,:)
         ! velocities:
         pu(:,:,:) = 0._wp
         pv(:,:,:) = 0._wp
      CASE(1)   ! strati
         pu  (:,:,:) = 0._wp        ! ocean at rest
         pv  (:,:,:) = 0._wp
         pssh(:,:)   = 0._wp
         !
         DO jk = 1, jpk             ! horizontally uniform T & S profiles
            DO jj = 1, jpj
               DO ji = 1, jpi
                  pts(ji,jj,jk,jp_tem) =  (  (  16. - 12. * TANH( (pdept(ji,jj,jk) - 400) / 700 ) )   &
                       &           * (-TANH( (500. - pdept(ji,jj,jk)) / 150. ) + 1.) / 2.             &
                       &           + ( 15. * ( 1. - TANH( (pdept(ji,jj,jk)-50.) / 1500.) )            &
                       &           - 1.4 * TANH((pdept(ji,jj,jk)-100.) / 100.)                        &
                       &           + 7.  * (1500. - pdept(ji,jj,jk) ) / 1500.)                        &
                       &           * (-TANH( (pdept(ji,jj,jk) - 500.) / 150.) + 1.) / 2.  ) * ptmask(ji,jj,jk)
                  
                  pts(ji,jj,jk,jp_sal) =  (  (  36.25 - 1.13 * TANH( (pdept(ji,jj,jk) - 305) / 460 ) )  &
                       &         * (-TANH((500. - pdept(ji,jj,jk)) / 150.) + 1.) / 2                  &
                       &         + ( 35.55 + 1.25 * (5000. - pdept(ji,jj,jk)) / 5000.                 &
                       &         - 1.62 * TANH( (pdept(ji,jj,jk) - 60.  ) / 650. )                    &
                       &         + 0.2  * TANH( (pdept(ji,jj,jk) - 35.  ) / 100. )                    &
                       &         + 0.2  * TANH( (pdept(ji,jj,jk) - 1000.) / 5000.) )                  &
                       &         * (-TANH( (pdept(ji,jj,jk) - 500.) / 150.) + 1.) / 2  ) * ptmask(ji,jj,jk)
               END DO
            END DO
         END DO
      CASE(2)   ! linear stratification in T
         pu  (:,:,:) = 0._wp        ! ocean at rest
         pv  (:,:,:) = 0._wp
         pssh(:,:)   = 0._wp
         !
         zTtop = 7._wp
         zTbot =  3._wp
         pts(:,:,:,jp_sal) = ptmask(:,:,:) * 35._wp
         pts(:,:,:,jp_tem) = ptmask(:,:,:) * (zTtop + pdept(:,:,:) * (zTbot - zTtop) / 4000._wp)
      CASE(3)   ! equatorial stratification on T, S = cst
         pu  (:,:,:) = 0._wp        ! ocean at rest
         pv  (:,:,:) = 0._wp
         pssh(:,:)   = 0._wp
         !
         DO jk = 1, jpk             ! horizontally uniform T & S profiles
            DO jj = 1, jpj
               DO ji = 1, jpi
                  pts(ji,jj,jk,jp_tem) =  (  (  16. - 12. * TANH( (pdept(ji,jj,jk) - 400) / 700 ) )   &
                       &           * (-TANH( (500. - pdept(ji,jj,jk)) / 150. ) + 1.) / 2.             &
                       &           + ( 15. * ( 1. - TANH( (pdept(ji,jj,jk)-50.) / 1500.) )            &
                       &           - 1.4 * TANH((pdept(ji,jj,jk)-100.) / 100.)                        &
                       &           + 7.  * (1500. - pdept(ji,jj,jk) ) / 1500.)                        &
                       &           * (-TANH( (pdept(ji,jj,jk) - 500.) / 150.) + 1.) / 2.  ) * ptmask(ji,jj,jk)
                  
                  pts(ji,jj,jk,jp_sal) =  35._wp * ptmask(ji,jj,jk)
               END DO
            END DO
         END DO
      CASE(4)   ! strati, with surface T going from about 23 degrees at the equator to 4 degrees at 60N !!rc TODO
         pu  (:,:,:) = 0._wp        ! ocean at rest
         pv  (:,:,:) = 0._wp
         pssh(:,:)   = 0._wp
         !
         ! vars for the profiles
         z1_700  = 1._wp /  700._wp
         z1_150  = 1._wp /  150._wp
         z1_1500 = 1._wp / 1500._wp
         z7_1500 = 7._wp / 1500._wp
         z1_100  = 1._wp /  100._wp
         z1_460  = 1._wp /  460._wp
         z1_5000 = 1._wp / 5000._wp
         z1_650  = 1._wp /  650._wp
         !
         ! horizontally uniform T & S profiles
         pts(:,:,:,jp_tem) =  (  (  16._wp - 12._wp * TANH( (pdept(:,:,:) - 400._wp) * z1_700 ) )            &
            &                  * (   1._wp -          TANH( (500._wp - pdept(:,:,:)) * z1_150 ) ) * 0.5_wp   &
            &                  + (  15._wp    * ( 1._wp - TANH( (  pdept(:,:,:) -    50._wp) * z1_1500) )    &
            &                      - 1.4_wp   * (         TANH( (  pdept(:,:,:) -   100._wp) *  z1_100) )    &
            &                      + z7_1500  * (             ( ( -pdept(:,:,:) +  1500._wp)          ) )    &
            &                    ) * (-TANH( (pdept(:,:,:) - 500._wp) * z1_150) + 1._wp) * 0.5_wp            &
            &                 ) * ptmask(:,:,:)
         !                   
         pts(:,:,:,jp_sal) =  (  (  36.25_wp - 1.13_wp * TANH( (pdept(:,:,:) - 305._wp) * z1_460 ) )  &
            &         * (-TANH((500._wp - pdept(:,:,:)) / 150._wp) + 1._wp) * 0.5_wp                   &
            &         + ( 35.55_wp + 1.25_wp * (5000._wp - pdept(:,:,:)) * z1_5000                  &
            &         - 1.62_wp * TANH( (pdept(:,:,:) - 60._wp  ) * z1_650 )                     &
            &         + 0.2_wp  * TANH( (pdept(:,:,:) - 35._wp  ) * z1_100 )                     &
            &         + 0.2_wp  * TANH( (pdept(:,:,:) - 1000._wp) * z1_5000) )                   &
            &         * (-TANH( (pdept(:,:,:) - 500._wp) * z1_150) + 1._wp) * 0.5_wp  ) * ptmask(:,:,:)
         ! Create horizontal gradient of T (going linearly from equatorial profile to uniform T=4 degC profile)
         zphiMAX = MAXVAL( gphit(:,:) )
         zTbot = MINVAL( pts(2:jpim1 , 2:jpjm1 , 1:jpkm1 , jp_tem) )   ! Must be a positive temperature
         IF( lk_mpp )   THEN
            zmpparr(1) = zphiMAX
            zmpparr(2) = - zTbot
            CALL mpp_max( 'usr_def_istate', zmpparr )
            zphiMAX = zmpparr(1)
            zTbot = - zmpparr(2)
         ENDIF
         z1_phiMAX = 1._wp / zphiMAX
         DO jk = 1, jpkm1
            pts(:,:,jk,jp_tem) = ( ( pts(:,:,jk,jp_tem) - zTbot ) * ( zphiMAX - gphit(:,:) ) * z1_phiMAX + zTbot ) * ptmask(:,:,jk)
         END DO
      CASE(5)    ! Test of geostrophic adjustment
         ! sea level:
         pssh(:,:) = -8._wp
         WHERE( gphit(:,:) < 30._wp )   pssh(:,:) = 8._wp
         ! temperature:
         pts(:,:,:,jp_tem) = 10._wp * ptmask(:,:,:)
         ! salinity:  
         pts(:,:,:,jp_sal) = 35._wp * ptmask(:,:,:)
         ! velocities:
         pu(:,:,:) = 0._wp
         pv(:,:,:) = 0._wp
      END SELECT

      CALL lbc_lnk( 'usrdef_istate', pssh, 'T',  1. )
      CALL lbc_lnk(  'usrdef_istate', pts, 'T',  1. )
      CALL lbc_lnk(   'usrdef_istate', pu, 'U', -1. )
      CALL lbc_lnk(   'usrdef_istate', pv, 'V', -1. )

   END SUBROUTINE usr_def_istate

   !!======================================================================
END MODULE usrdef_istate
