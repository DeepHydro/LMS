
    SUBROUTINE ADV_WASP
    USE MOD_WQM
    USE ALL_VARS
    USE BCS
    USE MOD_OBCS

#  if defined (SPHERICAL)
    USE MOD_SPHERICAL
#  if defined (NORTHPOLE)
    USE MOD_NORTHPOLE
#  endif   
#  endif
#  if defined (SEMI_IMPLICIT)
    USE MOD_SEMI_IMPLICIT
#  endif

#  if defined (WET_DRY)
    USE MOD_WD
#  endif

    IMPLICIT NONE
    REAL(SP), DIMENSION(0:MT,KB,NB)  :: RF,XFLUX,XFLUX_ADV
    REAL(SP), DIMENSION(0:MT)           :: PUPX,PUPY,PVPX,PVPY  
    REAL(SP), DIMENSION(0:MT)           :: PSPX,PSPY,PSPXD,PSPYD,VISCOFF
    REAL(SP), DIMENSION(3*(NT),KBM1) :: DTIJ 
    REAL(SP), DIMENSION(3*(NT),KBM1) :: UVN
    REAL(SP) :: UTMP,VTMP,SITAI,FFD,FF1,X11,Y11,X22,Y22,X33,Y33,TMP1,TMP2,XI,YI
    REAL(SP) :: DXA,DYA,DXB,DYB,FIJ1,FIJ2,UN
    REAL(SP) :: TXX,TYY,FXX,FYY,VISCOF,EXFLUX,TEMP,STPOINT
    REAL(SP) :: FACT,FM1
    REAL(SP) ::TT,TTIME
    INTEGER  :: I,I1,I2,IA,IB,J,J1,J2,K,JTMP,JJ,II,N1
    REAL(SP) :: S1MIN, S1MAX, S2MIN, S2MAX
    REAL(SP) :: WQM1MIN, WQM1MAX, WQM2MIN, WQM2MAX
#  if defined (SPHERICAL)
    REAL(DP) :: TY,TXPI,TYPI
    REAL(DP) :: XTMP1,XTMP
    REAL(DP) :: X1_DP,Y1_DP,X2_DP,Y2_DP,XII,YII
    REAL(DP) :: X11_TMP,Y11_TMP,X33_TMP,Y33_TMP
#  if defined (NORTHPOLE)
    REAL(DP) :: VX1_TMP,VY1_TMP,VX2_TMP,VY2_TMP
    REAL(DP) :: TXPI_TMP,TYPI_TMP
#  endif
#  endif
#  if defined (SEMI_IMPLICIT)
    REAL(SP) :: UN1
    REAL(SP), DIMENSION(3*(NT),KBM1) :: UVN1
    REAL(SP), DIMENSION(3*(NT),KBM1) :: DTIJ1
#  endif
#  if defined (MPDATA)
    REAL(SP) :: WQMMIN,WQMMAX,XXXX
    REAL(SP), DIMENSION(0:MT,KB)     :: WQM_S    
    REAL(SP), DIMENSION(0:MT,KB)     :: WQM_SF   
    REAL(SP), DIMENSION(0:MT,KB)     :: WWWS     
    REAL(SP), DIMENSION(0:MT,KB)     :: WWWSF   
    REAL(SP), DIMENSION(0:MT)        :: DTWWWS  
    REAL(SP), DIMENSION(0:MT,KB)     :: ZZZFLUX !! temporary total flux in corrected part
    REAL(SP), DIMENSION(0:MT,KB)     :: BETA    !! temporary beta coefficient in corrected part
    REAL(DP), DIMENSION(0:MT,KB)     :: BETAIN  !! temporary beta coefficient in corrected part
    REAL(SP), DIMENSION(0:MT,KB)     :: BETAOUT !! temporary beta coefficient in corrected part
    REAL(SP), DIMENSION(0:MT,KB)     :: WQM_FRESH    
    REAL(SP):: MY_TEMP,KR1
    INTEGER ITERA, NTERA
#  endif

    FACT = 0.0_SP
    FM1  = 1.0_SP
    IF(HORZMIX == 'closure') THEN
        FACT = 1.0_SP
        FM1  = 0.0_SP
    END IF

    !
    !--Initialize Fluxes-----------------------------------------------------------!
    !
    XFLUX     = 0.0_SP
    XFLUX_ADV = 0.0_SP

    !
    !--Loop Over Control Volume Sub-Edges And Calculate Normal Velocity------------!
    !

    DO I=1,NCV
        I1=NTRG(I)
        !     DTIJ(I)=DT1(I1)
        DO K=1,KBM1
            DTIJ(I,K) = DT1(I1)*DZ1(I1,K)
            UVN(I,K)  = V(I1,K)*DLTXE(I) - U(I1,K)*DLTYE(I)
#      if defined (SEMI_IMPLICIT)
            DTIJ1(I,K) = D1(I1)*DZ1(I1,K)
            UVN1(I,K) = VF(I1,K)*DLTXE(I) - UF(I1,K)*DLTYE(I)
#      endif
        END DO
    END DO

    TTIME=THOUR

    RF = 0.0_SP
    !
    !------- CALCULATE SOURCE AND SINK TERMS FOR EVERY VARIABLE -----------
    !
    DO I = 1,M
        IF(D(I) > 0.0_SP) Then
            DO K = 1, KBM1
                TT = 0.0_SP
                TT = T1(I,K)-20.0_SP !!JQI comment for test  !! JQI, In Zhengs code, TT=0

                 SODD(I,K) = SOD 
                 K_REAE(I,K) = K_DO_REA
                 !KR1 =   K_REAE(I,K)*TEMP_REAE**TT
                 KR1 =K_DO_REA
                 CS(I,K) = DO_Saturation
                 K_NITRR(I,K) = K_NITR
                DPP(I,K) = K_death
                F_ONN(I) = FON
                F_OPP(I) = FOP
                GPP(I,K) = K_GROW
                DPP(I,K) = K_mort
                   
                !-------------------- For dissolved oxygen ----------------------------
                !IF(H(I) <= 0.5_SP) THEN
                !    RF(I,K,1) =                                                         &
                !    K_REAE(I,K)*TEMP_REAE**TT*(CS(I,K)-WQM(I,K,1))-                     &
                !    K_DEOX*TEMP_DEOX**TT*WQM(I,K,1)*WQM(I,K,2)/(KBOD+WQM(I,K,1))-       &
                !    K_NITRR(I,K)*TEMP_NITR**TT*WQM(I,K,1)*WQM(I,K,4)/                   &
                !    (KNITR+WQM(I,K,1)) * 64.0_SP / 14.0_SP-                          &
                !    DPP(I,K)*WQM(I,K,3)*32.0_SP / 12.0_SP-                                &
                !    SODD(I,K) / 0.7_SP * TEMP_SOD**TT +                                       &
                !    GPP(I,K)*(32.0_SP / 12.0_SP+48.0_SP*RATIO_NC*(1-PNH3G(I,K)) / 14.0_SP) * &
                !    WQM(I,K,3) - K_RESP1 * 32 * 24 * 0.001_SP * 2                
                !ELSE
                    MY_TEMP=  K_REAE(I,K)*TEMP_REAE**TT*(CS(I,K)-WQM(I,K,1))
                    MY_TEMP= K_REAE(I,K)*(CS(I,K)-WQM(I,K,1))
                    
                    RF(I,K,1) =KR1 * (CS(I,K)-WQM(I,K,1))-                     &
                    K_DEOX*TEMP_DEOX**TT*WQM(I,K,1)*WQM(I,K,2) / (KBOD+WQM(I,K,1))-   &
                    K_resp*Temp_resp**TT*WQM(I,K,3)*32.0_SP/12.0_SP-                                &
                    K_NITRR(I,K)*TEMP_NITR**TT*WQM(I,K,1)*WQM(I,K,4) /(KNITR+WQM(I,K,1))*64.0_SP/14.0_SP-&
                    SODD(I,K)/D(I)*TEMP_SOD**TT+                                      &
                    GPP(I,K)*(32.0_SP/12.0_SP+48.0_SP /14.0_SP * 1.17 * (1 - PNH3G(I,K))) *  WQM(I,K,3)
                !END IF

                RF(I,K,1) = RF(I,K,1) / DAY_SEC

                !------------- For carbonaceous biochemical oxygen demand -------------           
                RF(I,K,2) =  32.0_SP/12.0_SP*DPP(I,K)*WQM(I,K,3)-     &
                K_DEOX*TEMP_DEOX**TT*WQM(I,K,1)*WQM(I,K,2)/(KBOD+WQM(I,K,1))-         &
                 K_DENI*TEMP_DENI**TT*WQM(I,K,5)*KNO3/(KNO3+WQM(I,K,1))*5.0_SP/4.0_SP * 32.0_SP / 14.0_SP- &
                WSS2*(1-FD2)*WQM(I,K,2)/MAX(D(I),1.5_SP)  
                 
                RF(I,K,2) = RF(I,K,2) / DAY_SEC

                !-------------------------- For phytoplankton -------------------------
                RF(I,K,3) =                                                           &
                GPP(I,K)*WQM(I,K,3)-DPP(I,K)*WQM(I,K,3)-WSS3*WQM(I,K,3)/MAX(D(I),1.5_SP)

                RF(I,K,3) = RF(I,K,3) / DAY_SEC

                !----------------------------- For ammonia ----------------------------
                RF(I,K,4) =                                                           &
                DPP(I,K)*RATIO_NC*(1-F_ONN(I))*WQM(I,K,3)+                            &
                K_MINE1*TEMP_MINE1**TT*WQM(I,K,3)*WQM(I,K,6)/ (KMPC+WQM(I,K,3))  -                      & 
               GPP(I,K)*RATIO_NC*PNH3G(I,K)*WQM(I,K,3)-                              &
                K_NITRR(I,K)*TEMP_NITR**TT*WQM(I,K,1)*WQM(I,K,4)/  (KNITR+WQM(I,K,1))           

                RF(I,K,4) = RF(I,K,4) / DAY_SEC

                !----------------------- For nitrate and nitrite ----------------------
                RF(I,K,5) =                                                           &
                K_NITRR(I,K)*TEMP_NITR**TT*WQM(I,K,1)*WQM(I,K,4)/ (KNITR+WQM(I,K,1)) -      &
                GPP(I,K)*RATIO_NC*(1-PNH3G(I,K))* WQM(I,K,3) -                &
                K_DENI*TEMP_DENI**TT*WQM(I,K,5)*KNO3/ (KNO3+WQM(I,K,1))
            
                RF(I,K,5) = RF(I,K,5) / DAY_SEC

                !------------------------ For organic nitrogen ------------------------
                RF(I,K,6) =                                                           &
                DPP(I,K)*RATIO_NC*F_ONN(I)*WQM(I,K,3)-                                &
                K_MINE1*TEMP_MINE1**TT*WQM(I,K,3)*WQM(I,K,6)/ (KMPC+WQM(I,K,3)) - &
                WSS3*WQM(I,K,6)*(1-FD6)/MAX(D(I),1.5_SP)                            

                RF(I,K,6) = RF(I,K,6) / DAY_SEC

                !--------------------- For inorganic phosphorus -----------------------
                RF(I,K,7) =                                                           &
                DPP(I,K)*RATIO_PC*(1-FOP)*WQM(I,K,3)+                                 &
                K_MINE2*TEMP_MINE2**TT*WQM(I,K,3)*WQM(I,K,8)/                         &
                (KMPC+WQM(I,K,3))-                                                    &
                GPP(I,K)*RATIO_PC*WQM(I,K,3)                                            

                RF(I,K,7) = RF(I,K,7) / DAY_SEC

                !---------------------- For organic phosphorus ------------------------
                RF(I,K,8) =                                                           &
                DPP(I,K)*RATIO_PC*F_OPP(I)*WQM(I,K,3)-                                &
                K_MINE2*TEMP_MINE2**TT*WQM(I,K,3)*WQM(I,K,8)/                         &
                (KMPC+WQM(I,K,3))-                                                    &
                WSS3*WQM(I,K,8)*(1-FD8)/MAX(D(I),1.5_SP)

                RF(I,K,8) = RF(I,K,8) / DAY_SEC

            END DO
        END IF
    END DO

    !
    !--Calculate the Advection and Horizontal Diffusion Terms----------------------!
    !
    DO N1=1,NB
        DO K=1,KBM1
            PSPX  = 0.0_SP 
            PSPY  = 0.0_SP 
            PSPXD = 0.0_SP 
            PSPYD = 0.0_SP

            DO I=1,M
                DO J=1,NTSN(I)-1
                    I1=NBSN(I,J)
                    I2=NBSN(I,J+1)

#    if defined (WET_DRY)
                    IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 1)THEN
                        FFD=0.5_SP*(WQM(I,K,N1)+WQM(I2,K,N1) -WMEAN(I,K,N1)-WMEAN(I2,K,N1))
                        FF1=0.5_SP*(WQM(I,K,N1)+WQM(I2,K,N1))
                    ELSE IF(ISWETN(I1) == 1 .AND. ISWETN(I2) == 0)THEN
                        FFD=0.5_SP*(WQM(I1,K,N1)+WQM(I,K,N1) -WMEAN(I1,K,N1)-WMEAN(I,K,N1))
                        FF1=0.5_SP*(WQM(I1,K,N1)+WQM(I,K,N1))
                    ELSE IF(ISWETN(I1) == 0 .AND. ISWETN(I2) == 0)THEN
                        FFD=0.5_SP*(WQM(I,K,N1)+WQM(I,K,N1) -WMEAN(I,K,N1)-WMEAN(I,K,N1))
                        FF1=0.5_SP*(WQM(I,K,N1)+WQM(I,K,N1))
                    ELSE
                        FFD=0.5_SP*(WQM(I1,K,N1)+WQM(I2,K,N1)-WMEAN(I1,K,N1)-WMEAN(I2,K,N1))
                        FF1=0.5_SP*(WQM(I1,K,N1)+WQM(I2,K,N1))
                    END IF 
#    else	 
                    FFD=0.5_SP*(WQM(I1,K,N1)+WQM(I2,K,N1) -WMEAN(I1,K,N1)-WMEAN(I2,K,N1))
                    FF1=0.5_SP*(WQM(I1,K,N1)+WQM(I2,K,N1))
#    endif	 

                    PSPX(I)=PSPX(I)+FF1*(VY(I1)-VY(I2))
                    PSPY(I)=PSPY(I)+FF1*(VX(I2)-VX(I1))
                    PSPXD(I)=PSPXD(I)+FFD*(VY(I1)-VY(I2))
                    PSPYD(I)=PSPYD(I)+FFD*(VX(I2)-VX(I1))

                END DO
                PSPX(I)=PSPX(I)/ART2(I)
                PSPY(I)=PSPY(I)/ART2(I)
                PSPXD(I)=PSPXD(I)/ART2(I)
                PSPYD(I)=PSPYD(I)/ART2(I)
            END DO

            IF(K == KBM1)THEN
                DO I=1,M
                    PFPXB(I) = PSPX(I)
                    PFPYB(I) = PSPY(I)
                END DO
            END IF

            DO I=1,M
                VISCOFF(I)=VISCOFH(I,K)
            END DO

            IF(K == KBM1) THEN
                AH_BOTTOM(1:M) = HORCON*(FACT*VISCOFF(1:M) + FM1)
            END IF

            DO I=1,NCV_I
                IA=NIEC(I,1)
                IB=NIEC(I,2)
                XI=0.5_SP*(XIJE(I,1)+XIJE(I,2))
                YI=0.5_SP*(YIJE(I,1)+YIJE(I,2))

                DXA=XI-VX(IA)
                DYA=YI-VY(IA)
                DXB=XI-VX(IB)
                DYB=YI-VY(IB)

                FIJ1=WQM(IA,K,N1)+DXA*PSPX(IA)+DYA*PSPY(IA)
                FIJ2=WQM(IB,K,N1)+DXB*PSPX(IB)+DYB*PSPY(IB)

                WQM1MIN=MINVAL(WQM(NBSN(IA,1:NTSN(IA)-1),K,N1))
                WQM1MIN=MIN(WQM1MIN, WQM(IA,K,N1))
                WQM1MAX=MAXVAL(WQM(NBSN(IA,1:NTSN(IA)-1),K,N1))
                WQM1MAX=MAX(WQM1MAX, WQM(IA,K,N1))
                WQM2MIN=MINVAL(WQM(NBSN(IB,1:NTSN(IB)-1),K,N1))
                WQM2MIN=MIN(WQM2MIN, WQM(IB,K,N1))
                WQM2MAX=MAXVAL(WQM(NBSN(IB,1:NTSN(IB)-1),K,N1))
                WQM2MAX=MAX(WQM2MAX, WQM(IB,K,N1))
                IF(FIJ1 < WQM1MIN) FIJ1=WQM1MIN
                IF(FIJ1 > WQM1MAX) FIJ1=WQM1MAX
                IF(FIJ2 < WQM2MIN) FIJ2=WQM2MIN
                IF(FIJ2 > WQM2MAX) FIJ2=WQM2MAX

                UN=UVN(I,K)
#      if defined (SEMI_IMPLICIT)
                UN1=UVN1(I,K)
#      endif

                VISCOF=HORCON*(FACT*(VISCOFF(IA)+VISCOFF(IB))*0.5_SP + FM1)

                TXX=0.5_SP*(PSPXD(IA)+PSPXD(IB))*VISCOF
                TYY=0.5_SP*(PSPYD(IA)+PSPYD(IB))*VISCOF

                FXX=-DTIJ(I,K)*TXX*DLTYE(I)
                FYY= DTIJ(I,K)*TYY*DLTXE(I)

#      if !defined (SEMI_IMPLICIT)
                EXFLUX=-UN*DTIJ(I,K)* &
                ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP+FXX+FYY
#      else
                EXFLUX=-UN*DTIJ(I,K)* &
                ((1.0_SP+SIGN(1.0_SP,UN))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN))*FIJ1)*0.5_SP
                EXFLUX=(1.0_SP-IFCETA)*EXFLUX+IFCETA*(-UN1*DTIJ1(I,K)*((1.0_SP+SIGN(1.0_SP,UN1))*FIJ2+(1.0_SP-SIGN(1.0_SP,UN1))*FIJ1)*0.5_SP)+FXX+FYY
#      endif

                XFLUX(IA,K,N1)=XFLUX(IA,K,N1)+EXFLUX
                XFLUX(IB,K,N1)=XFLUX(IB,K,N1)-EXFLUX

                XFLUX_ADV(IA,K,N1)=XFLUX_ADV(IA,K,N1)+(EXFLUX-FXX-FYY)
                XFLUX_ADV(IB,K,N1)=XFLUX_ADV(IB,K,N1)-(EXFLUX-FXX-FYY)
            END DO
        ENDDO !!SIGMA LOOP
    ENDDO

    DO N1=1,NB
# if defined (MPDATA)
        WQM_FRESH = WQM(:,:,N1)

        IF(POINT_ST_TYPE == 'calculated') THEN
            IF(INFLOW_TYPE == 'node') THEN
                IF(NUMQBC > 0) THEN
                    DO J=1,NUMQBC
                        JJ=INODEQ(J)
                        STPOINT=WDIS(J,N1)
                        DO K=1,KBM1
                            IF(QDIS(J) > 0 )   then
                                XFLUX(JJ,K,N1)=XFLUX(JJ,K,N1) - QDIS(J)*VQDIST(J,K)*STPOINT   !/DZ(K)
                            else
                                XFLUX(JJ,K,N1)=XFLUX(JJ,K,N1) + QDIS(J)*VQDIST(J,K)*STPOINT
                                XFLUX(JJ,K,N1)= -XFLUX(JJ,K,N1)
                            endif
                            !IF(QDIS(J) < 0 )    XFLUX(JJ,K,N1)= -XFLUX(JJ,K,N1)
                        END DO
                    END DO
                END IF
            ELSE IF(INFLOW_TYPE == 'edge') THEN
                IF(NUMQBC > 0) THEN
                    DO J=1,NUMQBC
                        J1=N_ICELLQ(J,1)
                        J2=N_ICELLQ(J,2)
                        STPOINT=WDIS(J,N1) !!ASK LIU SHOULD THIS BE STPOINT1(J1)/STPOINT2(J2)
                        DO K=1,KBM1
                            XFLUX(J1,K,N1)=XFLUX(J1,K,N1)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT   !/DZ(K)
                            XFLUX(J2,K,N1)=XFLUX(J2,K,N1)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT   !/DZ(K)
                        END DO
                    END DO
                END IF
            END IF
        END IF


        ! The horizontal term of advection is neglected here
        DO K=1,KBM1
            DO I=1,M
                IF(ISONB(I) == 2) THEN
                    XFLUX(I,K,N1)=0.
                ENDIF
            END DO
        END DO

        ! Initialize variables of MPDATA
        WQM_S=0._SP
        WQM_SF=0._SP
        WWWS=0._SP
        WWWSF=0._SP
        DTWWWS=0._SP
        ZZZFLUX=0._SP
        BETA=0._SP
        BETAIN=0._SP
        BETAOUT=0._SP

        !!   first loop for vertical upwind
        !!   flux including horizontal and vertical upwind
        DO K=1,KBM1
            DO I=1,M
#    if defined (WET_DRY)
                IF(ISWETN(I)*ISWETNT(I) == 1) THEN
#    endif
                    IF(K == 1) THEN
                        TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*WQM(I,K,N1)   &
                        -(WTS(I,K+1)+ABS(WTS(I,K+1)))*WQM(I,K+1,N1) &
                        +(WTS(I,K)+ABS(WTS(I,K)))*WQM(I,K,N1)    
                    ELSE IF(K == KBM1) THEN
                        TEMP = +(WTS(I,K)-ABS(WTS(I,K)))*WQM(I,K-1,N1)     &
                        +(WTS(I,K)+ABS(WTS(I,K)))*WQM(I,K,N1)
                    ELSE
                        TEMP = -(WTS(I,K+1)-ABS(WTS(I,K+1)))*WQM(I,K,N1)   &
                        -(WTS(I,K+1)+ABS(WTS(I,K+1)))*WQM(I,K+1,N1) &
                        +(WTS(I,K)-ABS(WTS(I,K)))*WQM(I,K-1,N1)     &
                        +(WTS(I,K)+ABS(WTS(I,K)))*WQM(I,K,N1)
                    END IF
                    TEMP = 0.5_SP*TEMP 

                    IF(K == 1)THEN
                        WQMMAX = MAXVAL(WQM(NBSN(I,1:NTSN(I)),K,N1))
                        WQMMIN = MINVAL(WQM(NBSN(I,1:NTSN(I)),K,N1))
                        WQMMAX = MAX(WQMMAX,WQM(I,K+1,N1),WQM(I,K,N1),WQM_FRESH(I,K))
                        WQMMIN = MIN(WQMMIN,WQM(I,K+1,N1),WQM(I,K,N1),WQM_FRESH(I,K))
                    ELSEIF(K == KBM1)THEN
                        WQMMAX = MAXVAL(WQM(NBSN(I,1:NTSN(I)),K,N1))
                        WQMMIN = MINVAL(WQM(NBSN(I,1:NTSN(I)),K,N1))
                        WQMMAX = MAX(WQMMAX,WQM(I,K-1,N1),WQM(I,K,N1),WQM_FRESH(I,K))
                        WQMMIN = MIN(WQMMIN,WQM(I,K-1,N1),WQM(I,K,N1),WQM_FRESH(I,K))
                    ELSE
                        WQMMAX = MAXVAL(WQM(NBSN(I,1:NTSN(I)),K,N1))
                        WQMMIN = MINVAL(WQM(NBSN(I,1:NTSN(I)),K,N1))
                        WQMMAX = MAX(WQMMAX,WQM(I,K+1,N1),WQM(I,K-1,N1),WQM(I,K,N1),WQM_FRESH(I,K))
                        WQMMIN = MIN(WQMMIN,WQM(I,K+1,N1),WQM(I,K-1,N1),WQM(I,K,N1),WQM_FRESH(I,K))
                    END IF

                    ZZZFLUX(I,K) = TEMP*(DTI/DT(I))/DZ(I,K) + XFLUX(I,K,N1)/ART1(I)*(DTI/DT(I))/DZ(I,K) 
                    XXXX = ZZZFLUX(I,K)*DT(I)/DTFA(I)+WQM(I,K,N1)-WQM(I,K,N1)*DT(I)/DTFA(I) 

                    BETA(I,K)=0.5*(1.-SIGN(1.,XXXX)) * (WQMMAX-WQM(I,K,N1))/(ABS(XXXX)+1.E-10) &
                    +0.5*(1.-SIGN(1.,-XXXX)) * (WQM(I,K,N1)-WQMMIN)/(ABS(XXXX)+1.E-10)

                    WQM_SF(I,K)=WQM(I,K,N1)-MIN(1.,BETA(I,K))*XXXX

#    if defined (WET_DRY)
                END IF
#    endif
            END DO
        END DO  !! SIGMA LOOP

        !----------------------------------------------------------------------------------------
        NTERA = 4
        DO ITERA=1,NTERA   !! Smolaricizw Loop 
            IF(ITERA == 1)THEN
                WWWSF  = WTS
                WQM_S   = WQM_SF
                DTWWWS = DT
            ELSE
                WWWSF  = WWWS
                WQM_S   = WQM_SF
                DTWWWS = DTFA
            END IF
            DO K=2,KBM1
                DO I=1,M
                    TEMP=ABS(WWWSF(I,K))-DTI*(WWWSF(I,K))*(WWWSF(I,K))/DZ(I,K)/DTWWWS(I)
                    WWWS(I,K)=TEMP*(WQM_S(I,K-1)-WQM_S(I,K))/(ABS(WQM_S(I,K-1))+ABS(WQM_S(I,K))+1.E-14)

                    IF(TEMP < 0.0_SP .OR. WQM_S(I,K) == 0.0_SP)THEN 
                        WWWS(I,K)=0. 
                    END IF
                END DO 
            END DO
            DO I=1,M
                WWWS(I,1)=0.
            END DO

            DO I=1,M
                WQMMAX = MAXVAL(WQM(NBSN(I,1:NTSN(I)),1,N1))
                WQMMIN = MINVAL(WQM(NBSN(I,1:NTSN(I)),1,N1))
                WQMMAX = MAX(WQMMAX,WQM(I,2,N1),WQM(I,1,N1),WQM_FRESH(I,1))
                WQMMIN = MIN(WQMMIN,WQM(I,2,N1),WQM(I,1,N1),WQM_FRESH(I,1))

                TEMP=0.5*((WWWS(I,2)+ABS(WWWS(I,2)))*WQM_S(I,2))*(DTI/DTFA(I))/DZ(I,1)
                BETAIN(I,1)=(WQMMAX-WQM_S(I,1))/(TEMP+1.E-10)

                TEMP=0.5*((WWWS(I,1)+ABS(WWWS(I,1)))*WQM_S(I,1)-        &
                (WWWS(I,2)-ABS(WWWS(I,2)))*WQM_S(I,1))*(DTI/DTFA(I))/DZ(I,1)
                BETAOUT(I,1)=(WQM_S(I,1)-WQMMIN)/(TEMP+1.E-10)

                WWWSF(I,1)=0.5*MIN(1.,BETAOUT(I,1))*(WWWS(I,1)+ABS(WWWS(I,1))) +0.5*MIN(1.,BETAIN(I,1))*(WWWS(I,1)-ABS(WWWS(I,1)))
            END DO

            DO K=2,KBM1-1
                DO I=1,M
                    WQMMAX = MAXVAL(WQM(NBSN(I,1:NTSN(I)),K,N1))
                    WQMMIN = MINVAL(WQM(NBSN(I,1:NTSN(I)),K,N1))
                    WQMMAX = MAX(WQMMAX,WQM(I,K+1,N1),WQM(I,K-1,N1),WQM(I,K,N1),WQM_FRESH(I,K))
                    WQMMIN = MIN(WQMMIN,WQM(I,K+1,N1),WQM(I,K-1,N1),WQM(I,K,N1),WQM_FRESH(I,K))

                    TEMP=0.5*((WWWS(I,K+1)+ABS(WWWS(I,K+1)))*WQM_S(I,K+1)-  &
                    (WWWS(I,K)-ABS(WWWS(I,K)))*WQM_S(I,K-1))*(DTI/DTFA(I))/DZ(I,K)
                    BETAIN(I,K)=(WQMMAX-WQM_S(I,K))/(TEMP+1.E-10)

                    TEMP=0.5*((WWWS(I,K)+ABS(WWWS(I,K)))*WQM_S(I,K)-        &
                    (WWWS(I,K+1)-ABS(WWWS(I,K+1)))*WQM_S(I,K))*(DTI/DTFA(I))/DZ(I,K)
                    BETAOUT(I,K)=(WQM_S(I,K)-WQMMIN)/(TEMP+1.E-10)

                    WWWSF(I,K)=0.5*MIN(1.,BETAIN(I,K-1),BETAOUT(I,K))*(WWWS(I,K)+ABS(WWWS(I,K))) + 0.5*MIN(1.,BETAIN(I,K),BETAOUT(I,K-1))*(WWWS(I,K)-ABS(WWWS(I,K)))
                END DO
            END DO

            K=KBM1
            DO I=1,M
                WQMMAX = MAXVAL(WQM(NBSN(I,1:NTSN(I)),K,N1))
                WQMMIN = MINVAL(WQM(NBSN(I,1:NTSN(I)),K,N1))
                WQMMAX = MAX(WQMMAX,WQM(I,K-1,N1),WQM(I,K,N1),WQM_FRESH(I,K))
                WQMMIN = MIN(WQMMIN,WQM(I,K-1,N1),WQM(I,K,N1),WQM_FRESH(I,K))

                TEMP=0.5*((WWWS(I,K+1)+ABS(WWWS(I,K+1)))*WQM_S(I,K+1)-  &
                (WWWS(I,K)-ABS(WWWS(I,K)))*WQM_S(I,K-1))*(DTI/DTFA(I))/DZ(I,K)
                BETAIN(I,K)=(WQMMAX-WQM_S(I,K))/(TEMP+1.E-10)

                TEMP=0.5*((WWWS(I,K)+ABS(WWWS(I,K)))*WQM_S(I,K)-        &
                (WWWS(I,K+1)-ABS(WWWS(I,K+1)))*WQM_S(I,K))*(DTI/DTFA(I))/DZ(I,K)
                BETAOUT(I,K)=(WQM_S(I,K)-WQMMIN)/(TEMP+1.E-10)

                WWWSF(I,K)=0.5*MIN(1.,BETAIN(I,K-1),BETAOUT(I,K))*(WWWS(I,K)+ABS(WWWS(I,K))) + 0.5*MIN(1.,BETAIN(I,K),BETAOUT(I,K-1))*(WWWS(I,K)-ABS(WWWS(I,K)))
            END DO

            WWWS=WWWSF 

            DO K=1,KBM1
                DO I=1,M
#      if defined (WET_DRY)
                    IF(ISWETN(I)*ISWETNT(I) == 1) THEN
#      endif
                        IF(K == 1) THEN
                            TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*WQM_S(I,K)   &
                            -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*WQM_S(I,K+1) &
                            +(WWWS(I,K)+ABS(WWWS(I,K)))*WQM_S(I,K)
                        ELSE IF(K == KBM1) THEN
                            TEMP = +(WWWS(I,K)-ABS(WWWS(I,K)))*WQM_S(I,K-1)     &
                            +(WWWS(I,K)+ABS(WWWS(I,K)))*WQM_S(I,K)
                        ELSE
                            TEMP = -(WWWS(I,K+1)-ABS(WWWS(I,K+1)))*WQM_S(I,K)   &
                            -(WWWS(I,K+1)+ABS(WWWS(I,K+1)))*WQM_S(I,K+1) &
                            +(WWWS(I,K)-ABS(WWWS(I,K)))*WQM_S(I,K-1)     &
                            +(WWWS(I,K)+ABS(WWWS(I,K)))*WQM_S(I,K)
                        END IF
                        TEMP = 0.5_SP*TEMP
                        WQM_SF(I,K)=(WQM_S(I,K)-TEMP*(DTI/DTFA(I))/DZ(I,K)) 
#      if defined (WET_DRY)
                    END IF
#      endif
                END DO
            END DO  !! SIGMA LOOP
        END DO  !! Smolarvizw Loop
#  endif

# if ! defined(MPDATA)
        !--------------------------------------------------------------------
        !   The central difference scheme in vertical advection
        !--------------------------------------------------------------------
        DO K=1,KBM1
            DO I=1,M
#    if defined (WET_DRY)
                IF(ISWETN(I)*ISWETNT(I) == 1) THEN
#    endif
                    IF(K == 1) THEN                    
                        TEMP=-WTS(I,K+1)*(WQM(I,K,N1)*DZ(I,K+1)+WQM(I,K+1,N1)*DZ(I,K))/(DZ(I,K)+DZ(I,K+1))          
                    ELSE IF(K == KBM1) THEN
                        TEMP=WTS(I,K)*(WQM(I,K,N1)*DZ(I,K-1)+WQM(I,K-1,N1)*DZ(I,K))/  (DZ(I,K)+DZ(I,K-1))
                    ELSE
                        TEMP=WTS(I,K)*(WQM(I,K,N1)*DZ(I,K-1)+WQM(I,K-1,N1)*DZ(I,K))/ (DZ(I,K)+DZ(I,K-1))- &
                        WTS(I,K+1)*(WQM(I,K,N1)*DZ(I,K+1)+WQM(I,K+1,N1)*DZ(I,K))/ (DZ(I,K)+DZ(I,K+1))
                    END IF

                    IF(ISONB(I) == 2) THEN
                        !         XFLUX(I,K)=TEMP*ART1(I)/DZ(I,K)
                        XFLUX(I,K,N1)=TEMP*ART1(I)
                    ELSE
                        !         XFLUX(I,K)=XFLUX(I,K)+TEMP*ART1(I)/DZ(I,K)
                        XFLUX(I,K,N1)=XFLUX(I,K,N1)+TEMP*ART1(I)
                    END IF
#    if defined (WET_DRY)
                END IF
#    endif
            END DO
        END DO  !! SIGMA LOOP

        !
        !--Set Boundary Conditions-For Fresh Water Flux--------------------------------!
        !
        IF(POINT_ST_TYPE == 'calculated') THEN
            IF(INFLOW_TYPE == 'node') THEN
                IF(NUMQBC > 0) THEN
                    DO J=1,NUMQBC
                        JJ=INODEQ(J)
                        STPOINT=SDIS(J)
                        DO K=1,KBM1
                            !             XFLUX(JJ,K)=XFLUX(JJ,K) - QDIS(J)*VQDIST(J,K)*STPOINT/DZ(JJ,K)
                            XFLUX(JJ,K,N1)=XFLUX(JJ,K,N1) - QDIS(J)*VQDIST(J,K)*STPOINT

                            IF(QDIS(J) < 0 ) XFLUX(JJ,K,N1)= -XFLUX(JJ,K,N1)

                        END DO
                    END DO
                END IF
            ELSE IF(INFLOW_TYPE == 'edge') THEN
                IF(NUMQBC > 0) THEN
                    DO J=1,NUMQBC
                        J1=N_ICELLQ(J,1)
                        J2=N_ICELLQ(J,2)
                        STPOINT=SDIS(J) !!ASK LIU SHOULD THIS BE STPOINT1(J1)/STPOINT2(J2)
                        DO K=1,KBM1
                            !             XFLUX(J1,K)=XFLUX(J1,K)-   &
                            !                         QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT/DZ(J1,K)
                            !             XFLUX(J2,K)=XFLUX(J2,K)-   &
                            !                         QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT/DZ(J2,K)
                            XFLUX(J1,K,N1)=XFLUX(J1,K,N1)-QDIS(J)*RDISQ(J,1)*VQDIST(J,K)*STPOINT
                            XFLUX(J2,K,N1)=XFLUX(J2,K,N1)-QDIS(J)*RDISQ(J,2)*VQDIST(J,K)*STPOINT
                        END DO
                    END DO
                END IF
            END IF
        END IF

# endif

        !--Update -------------------------------------------------------------!
        DO I=1,M
#  if defined (WET_DRY)
            IF(ISWETN(I)*ISWETNT(I) == 1 )THEN
#  endif
                DO K=1,KBM1
#    if defined (MPDATA)     
                    !       SF1(I,K)=(S1(I,K)-XFLUX(I,K)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))*(DT(I)/DTFA(I))

                    !WQM_F(I,K,N1)=WQM(I,K,N1)-XFLUX(I,K,N1)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))* (DT(I)/DTFA(I)+RF(I,K,N1)*DTI
                    ! WQM_F(I,K,N1)=(WQM(I,K,N1)-XFLUX(I,K,N1)/ART1(I)*(DTI/DT(I)*DZ(I,K)))* (DT(I)/D(I)) + RF(I,K,N1)*DTI !/DTFA(I)
                    WQM_F(I,K,N1)=(WQM(I,K,N1)-XFLUX(I,K,N1)/ART1(I)*(DTI/(DT(I)*DZ(I,K))))* (DT(I)/DTFA(I)) + RF(I,K,N1)*DTI ! /DTFA(I)
#    else
                    WQM_F(I,K,N1)=WQM_SF(I,K) + RF(I,K,N1)*DTI/DTFA(I)             
#    endif              
                  IF(WQM_F(I,K,N1)<=0)THEN
                      WQM_F(I,K,N1)=WQ_DEFAULT(N1)
                  ENDIF
                  IF(WQM_F(I,K,N1) > WQ_MAX(N1))THEN
                        WQM_F(I,K,N1) = WQ_MAX(N1)
                  ENDIF
                END DO
#  if defined (WET_DRY)
            ELSE
                DO K=1,KBM1
                    WQM_F(I,K,N1)=WQM(I,K,N1)
                END DO
            END IF
#  endif
        END DO
    ENDDO
    RETURN

    END SUBROUTINE ADV_WASP