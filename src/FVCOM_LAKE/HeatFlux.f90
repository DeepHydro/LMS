    
!    subroutine NetLongWaveFlux(cloudCover)
!        USE ALL_VARS
!        IMPLICIT NONE
!        REAL(SP) :: emissivity_air,c_cloud,cloudCover,reflection_air
!
!        reflection_air = 0.03_SP
!        cloudCover = 0.1_SP
!        emissivity_water = 0.975_SP
!
!
!        return
!    end subroutine NetLongWaveFlux

!  MODEL = 1 is recommended
    subroutine SensibleHeat(PAIR,TAIR,TWATER, RH,UWIND,SensiableH,ZU,MODEL)
        USE ALL_VARS
        IMPLICIT NONE
        REAL(SP) :: TAIR,TWATER,RH,SensiableH,RHOA
        REAL(SP) :: VPAIR,VPWATER,PAIR
        REAL(SP) :: CLAT, LE,CSEN
        REAL(SP) :: UWIND,U10,USR,USR1,ZU
        REAL(SP) :: QWATER,QAIR,K
        REAL(SP) :: TWATER_K,TAIR_K,DELTA_T
        INTEGER:: MODEL

        CSEN =  1.41E-3_SP
        TWATER_K = TWATER + Kelvin
        TAIR_K = TAIR + Kelvin
    
        CALL RHOAIR(PAIR,TAIR_K - Kelvin,RH,RHOA) 
        DELTA_T = TAIR_K - TWATER_K

        CALL WindVelocity(UWIND,U10,ZU)
        USR=U10

        IF(MODEL == 0)THEN
            SensiableH =  RHOA * CSEN * 1003._SP * USR * DELTA_T
        ELSE IF(MODEL == 1)THEN
             if(TAIR > TWATER) then
                SensiableH =  RHOA * CSEN * 1003._SP * U10 * DELTA_T
            else
                 SensiableH =  RHOA * CSEN * 4186._SP * U10 * DELTA_T
            end if
        ELSE IF(MODEL == 2)THEN
             USR1 = 9.2 + 0.46 * U10**2
             SensiableH = 0.47 *  DELTA_T * USR1
        ENDIF

        RETURN
   end subroutine SensibleHeat

     subroutine LatentHeat(PAIR,TAIR,TWATER, RH,UWIND,LatentH,ZU,MODEL)
        USE ALL_VARS
        IMPLICIT NONE
        REAL(SP) :: PAIR,TAIR,TWATER,RH,LatentH,RHOA
        REAL(SP) :: UWIND,U10,USR,USR1,ZU
        REAL(SP) :: CLAT, LE,CSEN,K
        REAL(SP) :: QWATER,QAIR,QSAT26S
        REAL(SP) ::TAIR_K, TWATER_K
        REAL(SP) :: VPAIR,VPWATER, PAIR_C
       INTEGER:: MODEL

        CSEN =  1.41E-3_SP
        CLAT = 1.3E-3_SP
        LE = 2.2572E6_SP 
        K = 5418
        TWATER_K = TWATER + Kelvin
        TAIR_K = TAIR + Kelvin

        CALL RHOAIR(PAIR,TAIR,RH,RHOA) 
        

        CALL WindVelocity(UWIND,U10,ZU)
        USR=U10

        IF(MODEL == 0)THEN
            call VaporPressure(TAIR,VPAIR)
            call VaporPressure(TWATER,VPWATER)
            PAIR_C = VPAIR * RH * 0.01_SP
             LatentH = RHOA * CLAT * LE * USR * (PAIR_C - VPWATER)
        ELSE IF(MODEL == 1)THEN
             call VaporPressure(TAIR,VPAIR)
             call VaporPressure(TWATER,VPWATER)
             PAIR_C = VPAIR * RH * 0.01_SP
             LatentH = RHOA * CLAT * LE * USR * (PAIR_C * 0.2167 / TAIR - VPWATER * 0.2167 / TWATER)
        ELSE IF(MODEL == 2)THEN
             CALL QSAT(TWATER,PAIR,QWATER)
             CALL QSAT(TAIR,PAIR,QAIR)
             LatentH =  LE * CLAT * U10 * ( QAIR - QWATER)
       ELSE IF(MODEL == 3)THEN
               CALL QSAT(TAIR,PAIR,QSAT26S)
               QAIR=0.98_SP*QSAT26S/1000._SP
               !  SPECIFIC HUMIDITY OF AIR (G/KG)  
               CALL QSAT(TWATER,PAIR,QSAT26S)
               QWATER=(0.01_SP*RH)*QSAT26S/1000._SP
               LatentH = RHOA * CLAT * LE *  U10 *  ( QAIR - QWATER)
        ENDIF
        RETURN
    end subroutine  LatentHeat


   SUBROUTINE SensibleLatentHeat(PAIR,TAIR,TWATER, RH,UWIND,LatentH,SensiableH,ZU)
    USE ALL_VARS
    IMPLICIT NONE
    REAL(SP) :: TAIR,TWATER,RH,LatentH,SensiableH,RHOA
    REAL(SP) :: VPAIR,VPWATER,PAIR,PAIR_C
    REAL(SP) :: UWIND,DU,DTER,UG,U10,USR,UT
    REAL(SP) :: QWATER,QAIR,TWATER_K,TAIR_K
    REAL(SP) :: SensiableH1,LatentH1,USR1,QSAT26S
    REAL(SP) :: CLAT, LE,CSEN,K
    REAL(SP) ::ZU

   CSEN =  1.41E-3_SP
    CLAT = 1.3E-3_SP
    LE = 2.2572E6_SP 
    K = 5418

    TWATER_K = TWATER + Kelvin
    TAIR_K = TAIR + Kelvin

   DU=UWIND
   UG=0.5_SP
   DTER=0.3_SP
   UT=SQRT(DU*DU+UG*UG) 
   U10=UT*LOG(10._SP/1e-4)/LOG(ZU/1e-4)
!   USR=0.035_SP*U10
  USR =U10 !! 
  USR1 = 9.2 + 0.46 * U10**2

    call VaporPressure(TAIR,VPAIR)
    call VaporPressure(TWATER,VPWATER)

    CALL RHOAIR(PAIR,TAIR,RH,RHOA) 
    PAIR_C = VPAIR * RH * 0.01_SP

    SensiableH =  RHOA * CSEN * 1003._SP * USR * ( TWATER- TAIR)
    LatentH = RHOA * CLAT * LE * USR * (PAIR_C - VPWATER)

    SensiableH = 0.47 *  ( TWATER- TAIR) * USR1
    LatentH = USR1 * (PAIR_C - VPWATER)

    if(TAIR > TWATER) then
         SensiableH =  RHOA * CSEN * 1003._SP * U10 * ( TAIR - TWATER)
    else
         SensiableH =  RHOA * CSEN * 4186._SP * U10 * ( TAIR - TWATER)
    end if

    QWATER = EXP(K * ( 1.0_SP / Kelvin -  1.0_SP / TWATER_K)) * 0.2167 * 6.11 / TWATER_K
    QAIR = EXP(K * ( 1.0_SP / Kelvin -  1.0_SP / TAIR_K)) * RH * 0.01_SP * 0.2167 * 6.11 / TAIR_K
    
    CALL QSAT(TWATER,PAIR,QWATER)
    CALL QSAT(TAIR,PAIR,QAIR)

    LatentH =  LE * CLAT * U10 * ( QAIR - QWATER)
    LatentH = RHOA * CLAT * LE * USR * (PAIR * 0.2167 / TAIR - VPWATER * 0.2167 / TWATER)

    CALL QSAT(TAIR,PAIR,QSAT26S)
    
   QAIR=0.98_SP*QSAT26S/1000._SP

    !  SPECIFIC HUMIDITY OF AIR (G/KG)  
   CALL QSAT(TWATER,PAIR,QSAT26S)
   QWATER=(0.01_SP*RH)*QSAT26S/1000._SP

    LatentH = RHOA * CLAT * LE * 0.035_SP* U10 *  ( QAIR - QWATER)
   END SUBROUTINE SensibleLatentHeat

  
   !-----------------------------------------------------------
   SUBROUTINE VaporPressure(TAIR,VP)
      USE ALL_VARS 
    ! computes saturation Vapor Pressure
       IMPLICIT NONE
       REAL(SP) :: TAIR,VP,TEMP
       TEMP= (18.678_SP - TAIR/ 234.5_SP) * TAIR / (TAIR+273.16_SP-38.66_SP)
       VP=6.1121_SP * EXP(TEMP)
       RETURN
     END SUBROUTINE VaporPressure

          !-----------------------------------------------------------
       SUBROUTINE LongWaveFlux(TAIR,TWATER,CLOUDCOVER,lwFlux)
       USE ALL_VARS
    ! computes saturation specific humidity
       IMPLICIT NONE
   
       REAL(SP) :: TAIR,TWATER,CLOUDCOVER,lwFlux
       REAL(SP) :: e_water,e_sp,e_temp
        REAL(SP) :: vp, cc, tk
        call VaporPressure(TAIR,vp)
       e_water=0.975_SP
       e_sp=5.67032e-8_SP
       e_temp=e_water*e_sp
       tk=273.16
       cc = 1.0_sp - 0.72_sp*CLOUDCOVER
       lwFlux=e_temp * (TAIR+tk)**4 * (0.05_sp * SQRT(vp) - 0.39_sp) * cc - 4.0_sp * e_temp *  (TWATER+tk)**3 * (TWATER - TAIR)
       RETURN
   END SUBROUTINE LongWaveFlux

   !-----------------------------------------------------------
   ! P(hPa) T(C) RH(%)
   SUBROUTINE   RHOAIR(P,T_C,RH,RHOA)
       USE ALL_VARS 
       REAL(SP)::P,T_C,RGAS,Q,RH
       REAL(SP)::QSAT26S,RHOA
   
       RGAS=287.1_SP

       CALL QSAT(T_C,P,QSAT26S)
       Q=(0.01_SP*RH)*QSAT26S/1000._SP

       RHOA = P*100._SP/(RGAS*(T_C+Kelvin)*(1+0.61_SP*Q))
       return
   END SUBROUTINE RHOAIR

   !-----------------------------------------------------------
   SUBROUTINE QSAT(T_C,P,QS)
        USE ALL_VARS 
    ! computes saturation specific humidity
       IMPLICIT NONE
       REAL(SP) :: T_C,P,QS,ES
   
        ES=6.112_SP*EXP(17.502_SP*T_C/(T_C+241.0_SP)) *(1.0007_SP+3.46_SP*0.000001_SP*P)
        QS=ES*622/(P-.378_SP*ES)
   
       RETURN
   END SUBROUTINE QSAT
