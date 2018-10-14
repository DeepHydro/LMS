!==============================================================================|
   SUBROUTINE VDIF_WASP(F)
!==============================================================================!
!									       !
!   This subroutine is used to calculate the eight variables of water 	       !
!   quality model in the Satilla River. They are: 			       !
!     (1) Dissolved Oxygen (DO)						       !
!     (2) Carbonaceous Biochemical Oxygen Demand (CBOD)			       !
!     (3) Phytoplankton (PHYT)						       !
!     (4) Ammonia Nitrogen (NH4)					       !
!     (5) Nitrate and Nitrite Nitrogen (NO3+NO2)			       !
!     (6) Organic Nitrogen (ON) 					       !
!     (7) Orthophosphorus or Inorganic Phosphorus (OPO4) 		       !
!     (8) Organic Phosphorus (OP) 					       !
!									       !
!   This subroutine is used to calculate the true water quality variables      !
!   by including vertical diffusion implicitly.		                       !
! 									       !
!==============================================================================!
   USE MOD_WQM
   USE ALL_VARS
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE
   REAL(SP), DIMENSION(0:MT,KB,NB) :: F
   REAL(DP), DIMENSION(M,KB,NB)    :: FF, VHF, VHPF 
   REAL(DP), DIMENSION(M,KB)       :: AF, CF, RAD
   REAL(SP), DIMENSION(M,NB)       :: BENFLUX,WFSURF
   REAL(SP), DIMENSION(M)          :: SOURCE1,SOURCE2,SOURCE3
   REAL(SP), DIMENSION(M)          :: TBOT
   REAL(SP) :: FKH,UMOLPR  
   REAL(SP) :: TEMPWUVBOT,TMP
   INTEGER  :: I,K,J,KI,N1

   UMOLPR = UMOL*1.E0_SP

!--- CALCULATE BOTTOM FLUX TERM --------------------------

   BENFLUX = 0.0_SP

   IF(BENWQM_KEY)THEN
     DO I = 1, M
#  if !defined (WET_DRY)
       IF(D(I) > 0.0_SP) THEN
#  else
       IF(ISWETN(I) == 1)THEN
#  endif
         SOURCE1(I) = 0.0_SP
         SOURCE2(I) = 0.0_SP
         SOURCE3(I) = 0.0_SP

         TEMPWUVBOT = 0.0_SP
         DO J = 1, NTVE(I)
           TMP=SQRT(WUBOT(NBVE(I,J))**2+WVBOT(NBVE(I,J))**2)
           TEMPWUVBOT = TEMPWUVBOT + TMP
         END DO
         TEMPWUVBOT = TEMPWUVBOT/FLOAT(NTVE(I))
         TBOT(I) = TEMPWUVBOT*1.0E+3_SP

         IF(H(I) >= 1.9_SP .AND. TBOT(I) > TCE2) THEN
           SOURCE1(I) = (TBOT(I) - TCE2) * 1.E-3_SP * RSED1(I)
           SOURCE2(I) = (TBOT(I) - TCE2) * 1.E-3_SP * RSED2(I)
           SOURCE3(I) = (TBOT(I) - TCE2) * 1.E-3_SP * RSED3(I)
         ELSE IF(H(I) >= 1.9_SP .AND. TBOT(I) < TCS2) THEN
           SOURCE1(I) = (TBOT(I) - TCS2) * 1.E-3_SP * RSED1(I)
           SOURCE2(I) = (TBOT(I) - TCS2) * 1.E-3_SP * RSED2(I)
           SOURCE3(I) = (TBOT(I) - TCS2) * 1.E-3_SP * RSED3(I)
           SOURCE1(I) = SOURCE1(I) * 0.05_SP
           SOURCE2(I) = SOURCE2(I) * 0.05_SP
           SOURCE3(I) = SOURCE3(I) * 0.05_SP
         END IF

         BENFLUX(I,4) = SOURCE1(I)
         BENFLUX(I,5) = SOURCE2(I)
         BENFLUX(I,7) = SOURCE3(I)

         BENFLUX(I,1) = (-DIFF_Z*(F(I,KB,1)-SEDWQM(I,1))/            &
                        DEP_BEN)/DAY_SEC + BENFLUX(I,1)
         BENFLUX(I,2) = (-DIFF_Z*(F(I,KB,2)*FD2-SEDWQM(I,2)*FDB2)/   &
                        DEP_BEN)/DAY_SEC + BENFLUX(I,2)
         BENFLUX(I,3) = (-DIFF_Z*F(I,KB,3)/DEP_BEN)/DAY_SEC          &
                        + BENFLUX(I,3)
         BENFLUX(I,4) = (-DIFF_Z*(F(I,KB,4)-SEDWQM(I,4)*FDB4)/       &
                        DEP_BEN)/DAY_SEC + BENFLUX(I,4)
         BENFLUX(I,5) = (-DIFF_Z*(F(I,KB,5)-SEDWQM(I,5)*FDB5)/       & 
                        DEP_BEN)/DAY_SEC + BENFLUX(I,5)
         BENFLUX(I,6) = (-DIFF_Z*(F(I,KB,6)*FD6-SEDWQM(I,6)*FDB5)/   &
                        DEP_BEN)/DAY_SEC + BENFLUX(I,6)
         BENFLUX(I,7) = (-DIFF_Z*(F(I,KB,7)-SEDWQM(I,7)*FDB7)/       & 
                        DEP_BEN)/DAY_SEC + BENFLUX(I,7)
         BENFLUX(I,8) = (-DIFF_Z*(F(I,KB,8)*FD8-SEDWQM(I,8)*FDB8)/   &
                        DEP_BEN)/DAY_SEC + BENFLUX(I,8)
       END IF
     END DO
   END IF
!----------------------------------------------------------------
!                                                                
!  the following section solves the equation               
!  dti*(kh*f')' -f=-fb
!                                                                
!----------------------------------------------------------------

   DO K = 2, KBM1
     DO I = 1, M
         FKH=KH(I,K)
         AF(I,K-1)=-DTI*(FKH+UMOLPR)/(DZ(I,K-1)*DZZ(I,K-1)*D(I)*D(I))
         CF(I,K)=-DTI*(FKH+UMOLPR)/(DZ(I,K)*DZZ(I,K-1)*D(I)*D(I))
     END DO
   END DO

   WFSURF = 0.0_SP

!------------------------------------------------
!  Surface BCs; WFSURF
!----------------------------------------------- 

   DO N1=1,NB
     DO I = 1, M
         VHF(I,1,N1) = AF(I,1) / (AF(I,1)-1.)
         VHPF(I,1,N1) = -DTI * WFSURF(I,N1) / (-DZ(I,1)*D(I)) - F(I,1,N1)
         VHPF(I,1,N1) = VHPF(I,1,N1) / (AF(I,1)-1.)
     END DO
   END DO
       
   DO N1=1,NB
     DO K = 2, KBM2
       DO I = 1, M
           VHPF(I,K,N1)=1./ (AF(I,K)+CF(I,K)*(1.-VHF(I,K-1,N1))-1.)
           VHF(I,K,N1) = AF(I,K) * VHPF(I,K,N1)
           VHPF(I,K,N1) = (CF(I,K)*VHPF(I,K-1,N1)- DBLE(F(I,K,N1)))*VHPF(I,K,N1)
       END DO
     END DO
   END DO

      DO N1=1,NB
     DO K = 1, KBM1 
       DO I = 1, M
           FF(I,K,N1) = F(I,K,N1)
       END DO
     END DO
      END DO
      
   !  FF = F
     
   DO N1=1,NB
     DO I = 1, M
       IF (D(I) > 0.0_SP .AND. ISONB(I) /= 2) THEN
         FF(I,KBM1,N1) = (CF(I,KBM1)*VHPF(I,KBM2,N1)-FF(I,KBM1,N1)   &
                         -DTI*BENFLUX(I,N1)/(D(I)*DZ(I,KBM1)))/  &
                         (CF(I,KBM1)*(1.-VHF(I,KBM2,N1))-1.)
       END IF
     END DO
   END DO

   DO N1=1,NB
     DO K = 2, KBM1
       KI = KB - K
       DO I = 1, M
           IF (D(I) > 0.0_SP .AND. ISONB(I) /= 2) THEN
           FF(I,KI,N1) = (VHF(I,KI,N1)*FF(I,KI+1,N1)+VHPF(I,KI,N1))
         END IF
       END DO
     END DO
   END DO

 DO N1 = 1, NB
         DO K = 1, KBM1
             DO I = 1, M      
                     F(I,K,N1) = FF(I,K,N1)
                     IF(F(I,K,N1)<=0)THEN
                         F(I,K,N1)=WQ_DEFAULT(N1)
                     ENDIF
                     IF(F(I,K,N1) > WQ_MAX(N1))THEN
                         F(I,K,N1) = WQ_MAX(N1)
                     ENDIF
             END DO
         ENDDO
     ENDDO
!
!----------------- CALCULATE BOTTOM CONCENTRATION ----------------------
!
   IF (BENWQM_KEY) THEN
     DO N1 = 1, NB
       DO K = 1, KBM1
         DO I = 1, M
#  if !defined (WET_DRY)
           IF (D(I) > 0.0_SP) THEN
#  else
           IF(ISWETN(I)*ISWETNT(I) == 1 )then
#  endif
             F(I,KB,N1) = F(I,KBM1,N1)-DZZ(I,KBM1)*D(I)*      &
                          BENFLUX(I,N1)/(KH(I,KBM1)+UMOLPR)
           END IF
         END DO
       END DO
     END DO
     END IF

   RETURN
   END SUBROUTINE VDIF_WASP