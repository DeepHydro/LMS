!==============================================================================|
!     this subroutine is used to calculate the true temperature                !
!     and salinity by including vertical diffusion implicitly.                 !
!==============================================================================|
! NOTE:  This subroutine has a lot of room for optimization.
!        Suggestions: 
!        Remove extra if/then statements inside double do loops
!        Switch loop order to only check if(wetdry) once per node            
!        Remove extra do loops- use array assignment!

   SUBROUTINE VDIF_TS(NCON1,F)                

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE BCS
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE
   INTEGER :: I,K,NCON1,J,KI,hutemp
!   REAL(SP) :: TMP,TMP1,TMP2,TMP3,QTMP,GW,ZDEP,FKH,UMOLPR,WETFAC
!   REAL(SP), DIMENSION(0:MT,KB)  :: F
!   REAL(SP), DIMENSION(M,KB)     :: FF,AF,CF,VHF,VHPF,RAD
!   REAL(SP), DIMENSION(M)        :: KHBOTTOM,WFSURF,SWRADF
   REAL(DP) :: TMP,TMP1,TMP2,TMP3,QTMP,GW,ZDEP,FKH,UMOLPR,WETFAC
   REAL(SP), DIMENSION(0:MT,KB)  :: F
   REAL(DP), DIMENSION(M,KB)     :: FF,AF,CF,VHF,VHPF,RAD
   REAL(DP), DIMENSION(M)        :: KHBOTTOM,WFSURF,SWRADF,SASURF
   REAL(DP), DIMENSION(M)        :: COSGAMA1,COSGAMA2

   UMOLPR = UMOL*1.E0_SP
   SASURF = 0.0_SP
!
!------------------------------------------------------------------------------!
!                                                                              !
!        the following section solves the equation                             !
!         dti*(kh*f')'-f=-fb                                                   !
!                                                                              !
!------------------------------------------------------------------------------!


   DO K = 2, KBM1
     DO I = 1, M
#  if !defined (WET_DRY)
       IF (D(I) > 0.0_SP) THEN
#  else
       IF(ISWETN(I) == 1)THEN
#  endif
         FKH = KH(I,K)

         IF(K == KBM1) THEN
           KHBOTTOM(I)=FKH
         END IF

         AF(I,K-1)=-DTI*(FKH+UMOLPR)/(DZ(I,K-1)*DZZ(I,K-1)*D(I)*D(I))
         CF(I,K)=-DTI*(FKH+UMOLPR)/(DZ(I,K)*DZZ(I,K-1)*D(I)*D(I))
       END IF
     END DO
   END DO


!------------------------------------------------------------------------------!
!     the net heat flux input.                                                 !
!     the method shown below can be used when we turn off the                  !
!     body force in subroutine advt. be sure this method could                 !
!     cause the surface overheated if the vertical resolution                  !
!     is not high enough.                                                      !
!------------------------------------------------------------------------------!

   IF(NCON1 == 1) THEN
       
   IF(H_TYPE /= 'body_h')THEN
!     DO I=1,M
!       WFSURF(I)=WTSURF(I)
!       SWRADF(I)=SWRAD(I)
!       SASURF(I)=0.0_SP
!     END DO
      WFSURF(1:M)=WTSURF(1:M)
      SWRADF(1:M)=SWRAD(1:M)

      SASURF=0.0_SP

     DO K = 1, KB
       DO I = 1, M
#    if !defined (WET_DRY)
         IF(D(I) > 0.0_SP) THEN
#    else
         IF(ISWETN(I) == 1)THEN
#    endif
           ZDEP = Z(I,K)*D(I)
           RAD(I,K)=SWRADF(I)*(RHEAT*EXP(ZDEP/ZETA1)+(1-RHEAT)*EXP(ZDEP/ZETA2))
         END IF
       END DO
     END DO
   ELSE 
     RAD    = 0.0_SP
     WFSURF = 0.0_SP
     SWRADF = 0.0_SP
   END IF 
        
   ELSE IF(NCON1 == 2) THEN

     DO I = 1, M
       SWRADF(I)= 0.0_SP
       WFSURF(I)=0.0_SP
       COSGAMA1(I)=1.0_SP
!       SASURF(I)=-(QEVAP3(I)-QPREC3(I))*F(I,1)*ROFVROS*COSGAMA1(I)
!---Considering the salinity conservation, the sea surface salinity flux-----!
!---can be set as zero, that is----------------------------------------------!
       SASURF(I) = 0.0_SP

!----------------------------------------------------------------------------!
       DO K=1,KB
         RAD(I,K)=0.0_SP
       END DO
     END DO

   ELSE
        
     PRINT*,'NCON NOT CORRECT IN VDIF_TS',NCON1
     CALL PSTOP

   END IF

!------------------------------------------------------------------------------!
!   surface bcs; wfsurf                                                        !
!------------------------------------------------------------------------------!

   DO I = 1, M
#  if !defined (WET_DRY)
     IF (D(I) > 0.0_SP) THEN
#  else
     IF(ISWETN(I) == 1)THEN
#  endif
       VHF(I,1) = AF(I,1) / (AF(I,1)-1.)
       VHPF(I,1) = -DTI *(SASURF(I)+WFSURF(I)-SWRADF(I) &
                   +RAD(I,1)-RAD(I,2)) / (-DZ(I,1)*D(I)) - F(I,1)
        VHPF(I,1) = VHPF(I,1) / (AF(I,1)-1.)
        END IF
   END DO

   DO K = 2, KBM2
     DO I = 1, M
#  if !defined (WET_DRY)
     IF (D(I) > 0.0_SP) THEN
#  else
     IF(ISWETN(I) == 1)THEN
#  endif
         VHPF(I,K)=1./ (AF(I,K)+CF(I,K)*(1.-VHF(I,K-1))-1.)
         VHF(I,K) = AF(I,K) * VHPF(I,K)
         VHPF(I,K) = (CF(I,K)*VHPF(I,K-1)-DBLE(F(I,K)) &
                     +DTI*(RAD(I,K)-RAD(I,K+1))/(D(I)*DZ(I,K)))*VHPF(I,K)
        END IF
     END DO
   END DO


   DO  K = 1, KBM1
     DO  I = 1, M
#  if !defined (WET_DRY)
     IF (D(I) > 0.0_SP) THEN
#  else
     IF(ISWETN(I) == 1)THEN
#  endif
         FF(I,K) = F(I,K)
       END IF
     END DO
   END DO

   DO I = 1, M
#  if !defined (WET_DRY)
     IF (D(I) > 0.0_SP .AND.ISONB(I) /= 2) THEN
#  else
     IF(ISWETN(I) == 1 .AND.ISONB(I) /= 2)THEN
#  endif
       TMP1=PFPXB(I)*COS(SITA_GD(I))+PFPYB(I)*SIN(SITA_GD(I))
       TMP2=AH_BOTTOM(I)*PHPN(I)
       TMP3=KHBOTTOM(I)+UMOLPR+AH_BOTTOM(I)*PHPN(I)*PHPN(I)
       TMP=TMP1*TMP2/TMP3*(KHBOTTOM(I)+UMOLPR)

! -------------------------------------------------------------------
       TMP = 0.0_SP
       !also try IF (TMP1 > 0.0_SP) TMP=0.0_SP
! -------------------------------------------------------------------


       GW=0.0_SP
!       IF(IBFW > 0) THEN
!         DO J=1,IBFW
!           IF(I == NODE_BFW(J).AND.(NCON1 == 2)) THEN
!!             QTMP=-(F(I,KBM1)*D1(I)*DZ(KBM1)*BFWDIS(J))/ &
!!                   (D1(I)*DZ(KBM1)*ART1(I)+BFWDIS(J))
!!             GW=DTI/D1(I)/DZ(KBM1)*QTMP
!!             TMP = 0.0_SP
!             COSGAMA2(I)=1.0_SP
!             TMP = BFWDIS(J)*F(I,KBM1)*ROFVROS*COSGAMA2(I)
!           END IF
!         END DO
!       END IF

       FF(I,KBM1) = ((CF(I,KBM1)*VHPF(I,KBM2)-FF(I,KBM1)-GW &
               +DTI*(RAD(I,KBM1)-RAD(I,KB)-TMP)/(D(I)*DZ(I,KBM1))) &
                /(CF(I,KBM1)*(1.-VHF(I,KBM2))-1._SP))
     END IF
   END DO

   DO  K = 2, KBM1
     KI = KB - K
     DO  I = 1, M
#  if !defined (WET_DRY)
     IF (D(I) > 0.0_SP .AND.ISONB(I) /= 2) THEN
#  else
     IF(ISWETN(I) == 1 .AND.ISONB(I) /= 2)THEN
#  endif
         FF(I,KI) = (VHF(I,KI)*FF(I,KI+1)+VHPF(I,KI))
       END IF
     END DO
   END DO

   DO I = 1, M
#  if defined (WET_DRY)
     IF(ISWETN(I)*ISWETNT(I) == 1 )then
#  endif
       DO K = 1, KBM1
         F(I,K) = FF(I,K)
       END DO
#  if defined (WET_DRY)
     END IF
#  endif
   END DO


   RETURN
   END SUBROUTINE VDIF_TS
!==============================================================================|
