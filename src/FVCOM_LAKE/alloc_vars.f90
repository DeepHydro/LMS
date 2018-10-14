!==============================================================================|
!    Allocate and Initialize Most Arrays                                       !
!==============================================================================|

   SUBROUTINE ALLOC_VARS

!==============================================================================!
   USE ALL_VARS
   IMPLICIT NONE
    
   INTEGER NCT,NDB,IERR
   REAL(SP) :: MEMTOT,MEMCNT
!==============================================================================!

   NDB = 1       !!GWC BASE THIS ON KIND
   NCT = NT*3
   MEMCNT = 0.

   IF(MSR)WRITE(IPT,*  )'!'
   IF(MSR)WRITE(IPT,*)'!                ALLOCATING MEMORY    '
   IF(MSR)WRITE(IPT,*  )'!'

!==============================================================================!
!  ALLOCATE:                                                                   !
!==============================================================================!

!--------------------------Grid Metrics---------------------------------------------!

   ALLOCATE(XC(0:NT))            ;XC   = ZERO   !!X-COORD AT FACE CENTER 
   ALLOCATE(YC(0:NT))            ;YC   = ZERO   !!Y-COORD AT FACE CENTER
   ALLOCATE(VX(0:MT))            ;VX   = ZERO   !!X-COORD AT GRID POINT
   ALLOCATE(VY(0:MT))            ;VY   = ZERO   !!Y-COORD AT GRID POINT
   ALLOCATE(ART(0:NT))           ;ART  = ZERO   !!AREA OF ELEMENT
   ALLOCATE(ART1(0:MT))          ;ART1 = ZERO   !!AREA OF NODE-BASE CONTROl VOLUME
   ALLOCATE(ART2(MT))            ;ART2 = ZERO   !!AREA OF ELEMENTS AROUND NODE

                  MEMCNT = MT*7*NDB + MEMCNT

!----------------Node, Boundary Condition, and Control Volume-----------------------!

   ! Commented by YT on 2017-5-12                
   !ALLOCATE(NV(0:NT,4))          
   NV       = 0  !!NODE NUMBERING FOR ELEMENTS
   ALLOCATE(NBE(0:NT,3))         ;NBE      = 0  !!INDICES OF ELMNT NEIGHBORS
   ALLOCATE(NTVE(0:MT))          ;NTVE     = 0 
   ALLOCATE(NTSN(MT))            ;NTSN     = 0 
   ALLOCATE(ISONB(0:MT))         ;ISONB    = 0  !!NODE MARKER = 0,1,2
   ALLOCATE(ISBCE(0:NT))         ;ISBCE    = 0 
   ALLOCATE(NIEC(NCT,2))         ;NIEC     = 0
   ALLOCATE(NTRG(NCT))           ;NTRG     = 0
   ALLOCATE(DLTXE(NCT))          ;DLTXE    = ZERO
   ALLOCATE(DLTYE(NCT))          ;DLTYE    = ZERO
   ALLOCATE(DLTXYE(NCT))         ;DLTXYE   = ZERO
   ALLOCATE(SITAE(NCT))          ;SITAE    = ZERO
   ALLOCATE(XIJE(NCT,2))         ;XIJE     = ZERO
   ALLOCATE(YIJE(NCT,2))         ;YIJE     = ZERO 
   ALLOCATE(N_ICELLQ(NUMQBC,2))  ;N_ICELLQ = 0

                 MEMCNT = NT*8 + MT*3  + NCT*3 + NCT*7*NDB  + MEMCNT

!----------------2-d arrays for the general vertical coordinate -------------------------------!

   ALLOCATE(Z(0:MT,KB))               ; Z    = ZERO    !!SIGMA COORDINATE VALUE 
   ALLOCATE(ZZ(0:MT,KB))              ; ZZ   = ZERO    !!INTRA LEVEL SIGMA VALUE
   ALLOCATE(DZ(0:MT,KB))              ; DZ   = ZERO    !!DELTA-SIGMA VALUE
   ALLOCATE(DZZ(0:MT,KB))             ; DZZ  = ZERO    !!DELTA OF INTRA LEVEL SIGMA 
   ALLOCATE(Z1(0:NT,KB))              ; Z1   = ZERO    !!SIGMA COORDINATE VALUE 
   ALLOCATE(ZZ1(0:NT,KB))             ; ZZ1  = ZERO    !!INTRA LEVEL SIGMA VALUE
   ALLOCATE(DZ1(0:NT,KB))             ; DZ1  = ZERO    !!DELTA-SIGMA VALUE
   ALLOCATE(DZZ1(0:NT,KB))            ; DZZ1 = ZERO    !!DELTA OF INTRA LEVEL SIGMA 
                     MEMCNT = MT*KB*4*NDB + NT*KB*4*NDB +MEMCNT

!---------------2-d flow variable arrays at elements-------------------------------!

   ALLOCATE(UA(0:NT))            ;UA        = ZERO  !!VERTICALLY AVERAGED X-VELOC
   ALLOCATE(VA(0:NT))            ;VA        = ZERO  !!VERTICALLY AVERAGED Y-VELOC
   ALLOCATE(UAF(0:NT))           ;UAF       = ZERO  !!UA FROM PREVIOUS RK STAGE 
   ALLOCATE(VAF(0:NT))           ;VAF       = ZERO  !!VA FROM PREVIOUS RK STAGE 
#  if !defined (SEMI_IMPLICIT)
   ALLOCATE(UARK(0:NT))          ;UARK      = ZERO  !!UA FROM PREVIOUS TIMESTEP 
   ALLOCATE(VARK(0:NT))          ;VARK      = ZERO  !!VA FROM PREVIOUS TIMESTEP 
   ALLOCATE(UARD(0:NT))          ;UARD      = ZERO  !!UA AVERAGED OVER EXTERNAL INT
   ALLOCATE(VARD(0:NT))          ;VARD      = ZERO  !!VA AVERAGED OVER EXTERNAL INT
#  endif
   ALLOCATE(COR(0:NT))           ;COR       = ZERO  !!CORIOLIS PARAMETER
   ALLOCATE(H1(0:NT))            ;H1        = ZERO  !!BATHYMETRIC DEPTH   
   ALLOCATE(D1(0:NT))            ;D1        = ZERO  !!DEPTH
   ALLOCATE(DT1(0:NT))           ;DT1       = ZERO  !!DEPTH  
   ALLOCATE(EL1(0:NT))           ;EL1       = ZERO  !!SURFACE ELEVATION
   ALLOCATE(ELF1(0:NT))          ;ELF1      = ZERO  !!SURFACE ELEVATION
   ALLOCATE(DTFA(0:MT))          ;DTFA      = ZERO  !!ADJUSTED DEPTH FOR MASS CONSERVATION
   ALLOCATE(ET1(0:NT))           ;ET1       = ZERO  !!SURFACE ELEVATION
#  if !defined (SEMI_IMPLICIT)
   ALLOCATE(ELRK1(0:NT))         ;ELRK1     = ZERO  !!SURFACE ELEVATION
#  endif
   ALLOCATE(CC_SPONGE(0:NT))     ;CC_SPONGE = ZERO  !!SPONGE DAMPING COEFFICIENT FOR MOMENTUM
                 MEMCNT = NT*17*NDB + MEMCNT

!---------------2-d flow variable arrays at nodes----------------------------------!

   ALLOCATE(H(0:MT))             ;H    = ZERO       !!BATHYMETRIC DEPTH   
   ALLOCATE(D(0:MT))             ;D    = ZERO       !!DEPTH   
   ALLOCATE(DT(0:MT))            ;DT   = ZERO       !!DEPTH   
   ALLOCATE(EL(0:MT))            ;EL   = ZERO       !!SURFACE ELEVATION
   ALLOCATE(ELF(0:MT))           ;ELF  = ZERO       !!SURFACE ELEVATION
   ALLOCATE(ET(0:MT))            ;ET   = ZERO       !!SURFACE ELEVATION
   ALLOCATE(EGF(0:MT))           ;EGF  = ZERO       !!SURFACE ELEVATION
#  if !defined (SEMI_IMPLICIT)
   ALLOCATE(ELRK(0:MT))          ;ELRK = ZERO       !!SURFACE ELEVATION
#  endif
               MEMCNT = MT*8*NDB + MEMCNT

!---------------surface/bottom/edge boundary conditions-----------------------------!

   ALLOCATE(CBC(0:NT))           ;CBC     = ZERO     !!BOTTOM FRICTION     
# if !defined (SEMI_IMPLICIT)
   ALLOCATE(WUSURF2(0:NT))       ;WUSURF2 = ZERO     !!SURFACE FRICTION FOR EXT
   ALLOCATE(WVSURF2(0:NT))       ;WVSURF2 = ZERO     !!SURFACE FRICTION FOR EXT
# endif
   ALLOCATE(WUBOT(0:NT))         ;WUBOT   = ZERO     !!BOTTOM FRICTION
   ALLOCATE(WVBOT(0:NT))         ;WVBOT   = ZERO     !!BOTTOM FRICTION
   ALLOCATE(WUSURF(0:NT))        ;WUSURF  = ZERO     !!SURFACE FRICTION FOR INT
   ALLOCATE(WVSURF(0:NT))        ;WVSURF  = ZERO     !!SURFACE FRICTION FOR INT
   ALLOCATE(UUWIND(0:NT))        ;UUWIND  = ZERO     !!SURFACE X-WIND 
   ALLOCATE(VVWIND(0:NT))        ;VVWIND  = ZERO     !!SURFACE Y-WIND 
   ALLOCATE(SWRAD(0:MT))         ;SWRAD   = ZERO     !!SURFACE INCIDENT RADIATION
   ALLOCATE(WTSURF(0:MT))        ;WTSURF  = ZERO 
#  if !defined (SEMI_IMPLICIT)
   ALLOCATE(QPREC2(0:MT))        ;QPREC2  = ZERO     !!SURFACE PRECIPITATION FOR EXT
   ALLOCATE(QEVAP2(0:MT))        ;QEVAP2  = ZERO     !!SURFACE EVAPORATION FOR EXT
#  endif
   ALLOCATE(QPREC3(0:MT))        ;QPREC3  = ZERO     !!SURFACE PRECIPITATION FOR INT
   ALLOCATE(QEVAP3(0:MT))        ;QEVAP3  = ZERO     !!SURFACE EVAPORATION FOR INT

#  if defined (ICE)
   ALLOCATE(T_AIR(0:MT))         ;T_AIR   = ZERO     !!BULK AIR TEMPERATURE
   ALLOCATE(RH_AIR(0:MT))        ;RH_AIR  = ZERO     !!RELATIVE HUMIDITY
   ALLOCATE(QA_AIR(0:MT))        ;QA_AIR  = ZERO     !!SPECIFIC HUMIDITY
   ALLOCATE(PA_AIR(0:MT))        ;PA_AIR  = ZERO     !!SURFACE PRESSURE
   ALLOCATE(DLW_AIR(0:MT))       ;DLW_AIR = ZERO     !!DOWNWARD LONGWAVE RADIATION
   ALLOCATE(DSW_AIR(0:MT))       ;DSW_AIR = ZERO     !!DOWNWARD SHORTWAVE RADIATION
   ALLOCATE(CLOUD(0:MT))         ;CLOUD   = ZERO     !!Total Cloud Cover
   ALLOCATE(WSSURF(0:MT))        ;WSSURF  = ZERO     !!sea surface salinity    flux
   ALLOCATE(WHSURF(0:MT))        ;WHSURF  = ZERO     !!sea surface heat        flux
   ALLOCATE(FWtSURF(0:MT))       ;FWtSURF = ZERO     !!sea surface fresh water flux
#  endif

!------------------------ GROUNDWATER T,S FLUX AT CURRENT TIM-------------------
   ALLOCATE(BFWDIS(IBFW))        ;BFWDIS  = ZERO     !!GROUNDWATER FLUX AT CURRENT TIME
   ALLOCATE(BFWTDIS(IBFW))       ;BFWTDIS = ZERO     !!GROUNDWATER TEMPERATURE FLUX AT CURRENT TIME
   ALLOCATE(BFWSDIS(IBFW))       ;BFWSDIS = ZERO     !!GROUNDWATER SALINITY FLUX AT CURRENT TIME
#  if !defined (SEMI_IMPLICIT)
   ALLOCATE(BFWDIS2(IBFW))       ;BFWDIS2  = ZERO     !!GROUNDWATER FLUX FOR EXT
!   ALLOCATE(BFWTDIS2(IBFW))      ;BFWTDIS2 = ZERO     !!GROUNDWATER TEMPERATURE FLUX AT CURRENT TIME FOR EXT
!   ALLOCATE(BFWSDIS2(IBFW))      ;BFWSDIS2 = ZERO     !!GROUNDWATER SALINITY FLUX AT CURRENT TIME FOR EXT
#  endif
   ALLOCATE(BFWDIS3(IBFW))       ;BFWDIS3  = ZERO     !!GROUNDWATER FLUX FOR INT
   ALLOCATE(BFWTDIS3(IBFW))      ;BFWTDIS3 = ZERO     !!GROUNDWATER TEMPERATURE FLUX AT CURRENT TIME FOR INT
   ALLOCATE(BFWSDIS3(IBFW))      ;BFWSDIS3 = ZERO     !!GROUNDWATER SALINITY FLUX AT CURRENT TIME FOR INT

!------------------------------------------------------------------------------
   ALLOCATE(QDIS(NUMQBC))        ;QDIS    = ZERO     !!FRESH WATER FLUX AT CURRENT TIME
#  if !defined (SEMI_IMPLICIT)
   ALLOCATE(QDIS2(NUMQBC))       ;QDIS2   = ZERO     !!FRESH WATER (EXT MODE NOT USED )  
#  endif
   ALLOCATE(TDIS(NUMQBC))        ;TDIS    = ZERO     !!FRESH WATER TEMP AT CURRENT TIME
   ALLOCATE(SDIS(NUMQBC))        ;SDIS    = ZERO     !!FRESH WATER SLNT AT CURRENT TIME
   ALLOCATE(QAREA(NUMQBC))       ;QAREA   = ZERO     !!AREA OF FLUX                    
   ALLOCATE(RDISQ(NUMQBC,2))     ;RDISQ   = ZERO     !!AREA OF FLUX                    
   ALLOCATE(ANGLEQ(NUMQBC))      ;ANGLEQ  = ZERO     !!FLUX VELOCITY                   
   ALLOCATE(VLCTYQ(NUMQBC))      ;VLCTYQ  = ZERO     !!FLUX ANGLE 
          MEMCNT = NT*9*NDB + MT*4*NDB + NUMQBC*9*NDB + MEMCNT

!-----------------------2-d flow fluxes---------------------------------------------!
   ALLOCATE(PSTX(0:NT))          ;PSTX  = ZERO       !!EXT MODE BAROTROPIC TERMS
   ALLOCATE(PSTY(0:NT))          ;PSTY  = ZERO       !!EXT MODE BAROTROPIC TERMS
   ALLOCATE(ADVUA(0:NT))         ;ADVUA = ZERO 
   ALLOCATE(ADVVA(0:NT))         ;ADVVA = ZERO
# if !defined (SEMI_IMPLICIT)
   ALLOCATE(ADX2D(0:NT))         ;ADX2D = ZERO
   ALLOCATE(ADY2D(0:NT))         ;ADY2D = ZERO
   ALLOCATE(DRX2D(0:NT))         ;DRX2D = ZERO
   ALLOCATE(DRY2D(0:NT))         ;DRY2D = ZERO
# endif
   ALLOCATE(TPS(0:NT))           ;TPS   = ZERO      !!WORKING ARRAY             
# if !defined (SEMI_IMPLICIT)
   ALLOCATE(ADVX(0:NT,KB))       ;ADVX  = ZERO 
   ALLOCATE(ADVY(0:NT,KB))       ;ADVY  = ZERO 
# endif
               MEMCNT = NT*9*NDB + NT*KB*2*NDB + MEMCNT


!---------------- internal mode   arrays-(element based)----------------------------!

   ALLOCATE(U(0:NT,KB))          ;U       = ZERO   !!X-VELOCITY
   ALLOCATE(V(0:NT,KB))          ;V       = ZERO   !!Y-VELOCITY
   ALLOCATE(UBETA(0:NT,KB))      ;UBETA   = ZERO   !X-VELOCITY temp time step without Coriolis
   ALLOCATE(VBETA(0:NT,KB))      ;VBETA   = ZERO   !Y-VELOCITY temp time step without Coriolis
   ALLOCATE(UBETA2D(0:NT))       ;UBETA2D = ZERO
   ALLOCATE(VBETA2D(0:NT))       ;VBETA2D = ZERO

   ALLOCATE(W(0:NT,KB))          ;W     = ZERO   !!VERTICAL VELOCITY IN SIGMA SYSTEM
   ALLOCATE(WW(0:NT,KB))         ;WW    = ZERO   !!Z-VELOCITY
   ALLOCATE(UF(0:NT,KB))         ;UF    = ZERO   !!X-VELOCITY FROM PREVIOUS TIMESTEP
   ALLOCATE(VF(0:NT,KB))         ;VF    = ZERO   !!Y-VELOCITY FROM PREVIOUS TIMESTEP
   ALLOCATE(WT(0:NT,KB))         ;WT    = ZERO   !!Z-VELOCITY FROM PREVIOUS TIMESTEP
   ALLOCATE(RHO(0:NT,KB))        ;RHO   = ZERO   !!DENSITY AT ELEMENTS
   ALLOCATE(RMEAN(0:NT,KB))      ;RMEAN = ZERO   !!MEAN INITIAL DENSITY AT ELEMENTS
   ALLOCATE(T(0:NT,KB))          ;T     = ZERO   !!TEMPERATURE AT ELEMENTS
   ALLOCATE(TMEAN(0:NT,KB))      ;TMEAN = ZERO   !!MEAN INITIAL TEMPERATURE AT ELEMENTS
   ALLOCATE(S(0:NT,KB))          ;S     = ZERO   !!SALINITY AT ELEMENTS
   ALLOCATE(SMEAN(0:NT,KB))      ;SMEAN = ZERO   !!MEAN INITIAL SALINITY AT ELEMENTS
               MEMCNT = NT*KB*13*NDB + MEMCNT

!-----------------------3d variable arrays-(node based)-----------------------------!

   ALLOCATE(T1(0:MT,KB))         ;T1     = ZERO  !!TEMPERATURE AT NODES               
   ALLOCATE(S1(0:MT,KB))         ;S1     = ZERO  !!SALINITY AT NODES               
   ALLOCATE(RHO1(0:MT,KB))       ;RHO1   = ZERO  !!DENSITY AT NODES               
   ALLOCATE(TF1(0:MT,KB))        ;TF1    = ZERO  !!TEMPERATURE FROM PREVIOUS TIME
   ALLOCATE(SF1(0:MT,KB))        ;SF1    = ZERO  !!SALINITY FROM PREVIOUS TIME 
   ALLOCATE(TMEAN1(0:MT,KB))     ;TMEAN1 = ZERO  !!MEAN INITIAL TEMP
   ALLOCATE(SMEAN1(0:MT,KB))     ;SMEAN1 = ZERO  !!MEAN INITIAL SALINITY 
   ALLOCATE(RMEAN1(0:MT,KB))     ;RMEAN1 = ZERO  !!MEAN INITIAL DENSITY 
   ALLOCATE(WTS(0:MT,KB))        ;WTS    = ZERO  !!VERTICAL VELOCITY IN SIGMA SYSTEM
   ALLOCATE(WTTS(0:MT,KB))       ;WTTS   = ZERO  !!WTS FROM PREVIOUS TIMESTEP        
   ALLOCATE(Q2(0:MT,KB))         ;Q2     = ZERO   !!TURBULENT KINETIC ENERGY AT NODES
   ALLOCATE(Q2L(0:MT,KB))        ;Q2L    = ZERO   !!TURBULENT KE*LENGTH AT NODES
#  if defined (GOTM)
   ALLOCATE(TKE(0:MT,KB))        ;TKE   = ZERO   !!TURBULENT KINETIC ENERGY AT ELEMENTS
   ALLOCATE(TKEF(0:MT,KB))       ;TKEF  = ZERO   !!TURBULENT KINETIC ENERGY AT ELEMENTS
   ALLOCATE(TEPS(0:MT,KB))       ;TEPS  = ZERO   !!TURBULENT KINETIC ENERGY AT ELEMENTS
   ALLOCATE(TEPSF(0:MT,KB))      ;TEPSF = ZERO   !!TURBULENT KE*LENGTH AT ELEMENTS
#  endif
   ALLOCATE(L(0:MT,KB))          ;L     = ZERO   !!TURBULENT LENGTH SCALE AT ELEMENTS
   ALLOCATE(KM(0:MT,KB))         ;KM    = ZERO   !!TURBULENT QUANTITY
   ALLOCATE(KH(0:MT,KB))         ;KH    = ZERO   !!TURBULENT QUANTITY
   ALLOCATE(KQ(0:MT,KB))         ;KQ    = ZERO   !!TURBULENT QUANTITY
   ALLOCATE(AAM(0:MT,KB))        ;AAM   = ZERO   !!??
   ALLOCATE(KM1(0:NT,KB))        ;KM1   = ZERO   !!TURBULENT QUANTITY AT ELEMENTS
               MEMCNT = MT*KB*18*NDB + NT*KB*NDB + MEMCNT
  
!---------------------------internal mode fluxes------------------------------------!

   ALLOCATE(DRHOX(0:NT,KB))      ;DRHOX = ZERO 
   ALLOCATE(DRHOY(0:NT,KB))      ;DRHOY = ZERO 
   ALLOCATE(Q2F(0:MT,KB))        ;Q2F   = ZERO 
   ALLOCATE(Q2LF(0:MT,KB))       ;Q2LF  = ZERO
               MEMCNT = NT*KB*2*NDB + MT*KB*2*NDB + MEMCNT

!------------shape coefficient arrays and control volume metrics--------------------!

   ALLOCATE(A1U(0:NT,4))         ;A1U   = ZERO
   ALLOCATE(A2U(0:NT,4))         ;A2U   = ZERO 
   ALLOCATE(AWX(0:NT,3))         ;AWX   = ZERO 
   ALLOCATE(AWY(0:NT,3))         ;AWY   = ZERO 
   ALLOCATE(AW0(0:NT,3))         ;AW0   = ZERO 
   ALLOCATE(ALPHA(0:NT))         ;ALPHA = ZERO
               MEMCNT = NT*4*2 + NT*3*3 + NT + MEMCNT

!-----salinity and temperature bottom diffusion condition/bottom depth gradients----!

   ALLOCATE(PHPN(MT))            ;PHPN      = ZERO 
   ALLOCATE(PFPXB(MT))           ;PFPXB     = ZERO
   ALLOCATE(PFPYB(MT))           ;PFPYB     = ZERO
   ALLOCATE(SITA_GD(MT))         ;SITA_GD   = ZERO 
   ALLOCATE(AH_BOTTOM(MT))       ;AH_BOTTOM = ZERO 
               MEMCNT = MT*5*NDB + MEMCNT

   ALLOCATE(VISCOFH(0:MT,KB))    ;VISCOFH = ZERO
               MEMCNT = MT*KB*NDB + MEMCNT

!---------------Coordinates of Center Pionts around the Nodes-----------------------!
   ALLOCATE(XCA(M))             ;XCA        = ZERO
   ALLOCATE(YCA(M))             ;YCA        = ZERO
   ALLOCATE(XCB(M))             ;XCB        = ZERO
   ALLOCATE(YCB(M))             ;YCB        = ZERO
   ALLOCATE(XCC(M,10))          ;XCC        = ZERO !ASSUMING THE MAXIUM OF NEIGHBORING ELEMENTS IS NOT MORE THAN 10
   ALLOCATE(YCC(M,10))          ;YCC        = ZERO !ASSUMING THE MAXIUM OF NEIGHBORING ELEMENTS IS NOT MORE THAN 10
   ALLOCATE(XCD(M,10))          ;XCD        = ZERO !ASSUMING THE MAXIUM OF NEIGHBORING ELEMENTS IS NOT MORE THAN 10
   ALLOCATE(YCD(M,10))          ;YCD        = ZERO !ASSUMING THE MAXIUM OF NEIGHBORING ELEMENTS IS NOT MORE THAN 10
   ALLOCATE(XCE(M))             ;XCE        = ZERO
   ALLOCATE(YCE(M))             ;YCE        = ZERO
   ALLOCATE(XCF(M))             ;XCF        = ZERO
   ALLOCATE(YCF(M))             ;YCF        = ZERO
   ALLOCATE(VAL_COS_VY(M))      ;VAL_COS_VY = ZERO

!-----special initialization which probably do nothing------------------------------!

   DT1(0) = 100.0_SP

!---------------report approximate memory usage-------------------------------------!

   MEMTOT = MEMCNT*4
#  if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_REDUCE(MEMCNT,MEMTOT,1,MPI_F,MPI_SUM,0,MPI_COMM_WORLD,IERR)
#  endif
   IF(MSR)WRITE(IPT,*)'!  # MBYTES TOTAL MEM    :',MEMTOT/1E+6
   IF(MSR .AND. .NOT.SERIAL )WRITE(IPT,*)'!  # AVERAGE MBYTES/PROC :',MEMTOT/(1E+6*NPROCS)
   
!---------------my variable definition-------------------------------------!
StefanBoltz = 5.67032e-8_sp
Kelvin = 273.16
R2D = 180 / PI
Language= 'zh-CN'
   RETURN
   END SUBROUTINE ALLOC_VARS
!==============================================================================|
