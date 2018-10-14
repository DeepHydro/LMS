!    module Mod_DotFVMProject
!    USE ALL_VARS
!    USE MOD_UTILS
!    USE MOD_INP
!    use READ_XML_PRIMITIVES
!    use XMLPARSE 
!    USE MOD_WQM, ONLY: WQM_ON,NB,BENWQM_KEY
!    USE MOD_WD, ONLY: WET_DRY_ON
!    USE MOD_HEATFLUX, ONLY: ZU,ZT,ZQ,C_HFX
!    USE MOD_OBCS, only:TYPE_TSOBC
!    !USE MOD_ASSIM, ONLY:SST_ASSIM , SST_METHOD, RAD_SST , GAMA_SST ,&
!    !    GALPHA_SST , ASTIME_WINDOW_SST,N_INFLU_SST,IAV_DAY 
!    implicit none
!
!    contains
!    
!    SUBROUTINE OpenProject    
!    CHARACTER(LEN=120) :: PrjName
!     logical :: errout=.true.
!    CALL get_command_argument(1, PrjName)
!   
!    call ReadPrjXmlFile(PrjName,errout)
!    END SUBROUTINE OpenProject
!
!    SUBROUTINE ReadPrjXmlFile(fname, errout)
!    
!    character(len=*), intent(in)   :: fname
!    logical, intent(out), optional   :: errout
!
!    type(XML_PARSE) :: info
!    logical :: error
!    logical :: endtag
!    logical:: has_var
!    
!    character(len=80)   :: tag
!    character(len=80)    :: starttag
!    character(len=80), dimension(1:2,1:20) :: attribs
!    character(len=200), dimension(1:10)  :: innerdata
!    integer :: nodata
!    integer :: noattribs
!    
!     real(SP), dimension(:), pointer  :: p_real_array
!     
!    errout=.false.
!    
!    call xml_open( info, fname, .true. )
!    call xml_options( info, report_errors=.true., ignore_whitespace=.true.)
!   
!   do
!      call xml_get( info, tag, endtag, attribs, noattribs, innerdata, nodata )
!      if ( xml_error(info) ) then
!         call STOPRUN('Error reading input file!')
!         error = .true.
!         return
!      endif
!      if ( endtag .and. tag .eq. starttag ) then
!         exit
!      endif
!      if ( endtag .and. noattribs .eq. 0 ) then
!         if ( xml_ok(info) ) then
!            cycle
!         else
!            exit
!         endif
!      endif
!      
!      !**********************Time Integration*************
!      select case( tag )
!          case('TimeStep')
!               call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,DTE, has_var)
!          case('Split')
!               call read_xml_integer(info, tag, endtag, attribs, noattribs, innerdata, nodata,ISPLIT, has_var)
!          case('RampUp')
!               call read_xml_integer(info, tag, endtag, attribs, noattribs, innerdata, nodata,IRAMP, has_var)
!          case('Cycles')
!               call read_xml_integer(info, tag, endtag, attribs, noattribs, innerdata, nodata,NSTEPS, has_var)
!      end select
!     !**********************InputOutput*************    
!      select case( tag )
!           case('CaseName') 
!               call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,CASENAME, has_var)
!          case('InputDir') 
!               call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,INPDIR, has_var)
!          case('OutputDir') 
!               call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,OUTDIR, has_var)
!          case('ScreenFile')
!               call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,INFOFILE, has_var)
!           case('SceenInteval')
!               call read_xml_integer(info, tag, endtag, attribs, noattribs, innerdata, nodata,IREPORT, has_var)
!               
!         case('OutputOn')
!               call read_xml_logical(info, tag, endtag, attribs, noattribs, innerdata, nodata,Output_ON, has_var)
!          case('OutputFormat')
!               call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,OUT_FORMAT, has_var)
!            case('OutputInteval')
!               call read_xml_integer(info, tag, endtag, attribs, noattribs, innerdata, nodata,IRECORD, has_var)           
!               
!         case('SMSOn')
!                call read_xml_logical(info, tag, endtag, attribs, noattribs, innerdata, nodata,SMS_ON, has_var)
!          case('SMSOutMode')
!                call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,SMS_OUT_Mode, has_var)  
!          case('SMSInteval')
!               call read_xml_integer(info, tag, endtag, attribs, noattribs, innerdata, nodata,IDMPSMS, has_var)              
!            
!          case('RestartInteval')
!               call read_xml_integer(info, tag, endtag, attribs, noattribs, innerdata, nodata,IRESTART, has_var)      
!          case('Restart')
!                call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,RESTART, has_var)  
!                
!          case('WQIniValueMode')
!                call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,InitialWQ_Mode, has_var)
!          end select
!          
!          !**********************Models*************    
!         select case( tag )
!               case('BarotropicModelOn')   
!                    call read_xml_logical(info, tag, endtag, attribs, noattribs, innerdata, nodata,BAROTROPIC, has_var)
!               case('SalinityModelOn')   
!                    call read_xml_logical(info, tag, endtag, attribs, noattribs, innerdata, nodata,SALINITY_ON, has_var)
!               case('TemperatureModelOn')   
!                    call read_xml_logical(info, tag, endtag, attribs, noattribs, innerdata, nodata,TEMP_ON, has_var)
!                case('WaterQualityModelOn')   
!                    call read_xml_logical(info, tag, endtag, attribs, noattribs, innerdata, nodata,WQM_ON, has_var)
!                 case('WetDryModelOn')   
!                    call read_xml_logical(info, tag, endtag, attribs, noattribs, innerdata, nodata,WET_DRY_ON, has_var)      
!         end select
!         
!        !**********************HydrodynamicsParas*************    
!         select case( tag )
!             case('VerticalLayers') 
!                call read_xml_integer(info, tag, endtag, attribs, noattribs, innerdata, nodata,KB, has_var)
!             
!            case('BtmDragCoef') 
!                call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,BFRIC, has_var)
!            case('BtmFrctDepLenScale') 
!                call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,Z0B, has_var)
!            case('BtmRoughModel') 
!                call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,BROUGH_TYPE, has_var)  
!           case('HorzDiffModel') 
!                call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,HORZMIX, has_var)     
!            case('HorzDiffCoef') 
!                call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,HORCON, has_var)    
!            case('HorzKeneticViso') 
!                call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,HPRNU, has_var)   
!                
!           case('VetcDiffModel') 
!                call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,VERTMIX, has_var)    
!            case('VetcDiffCoef') 
!                call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,UMOL, has_var)   
!            case('VetcKeneticViso') 
!                call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,VPRNU, has_var)   
!
!              case('BroClinPresGredn') 
!                call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,C_BAROPG, has_var)    
!              case('DensityModel') 
!                call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,CTRL_DEN, has_var)         
!              case('AdjustFlowStability')  
!                call read_xml_logical(info, tag, endtag, attribs, noattribs, innerdata, nodata,VERT_STAB, has_var)           
!               case('AdjustTempSalt')  
!                call read_xml_logical(info, tag, endtag, attribs, noattribs, innerdata, nodata,TS_FCT, has_var)     
!              end select
!              
!               !**********************Depthes*************    
!         select case( tag )
!                case('MinDepth')  
!                    call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,MIN_DEPTH, has_var)    
!               case('DepthAdjustFactor')  
!                    call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,DJUST, has_var)   
!              case('BaseSurfaceElevation')
!                    call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,BaseSurfaceElevation, has_var)
!               case('StandardDepthCounts') 
!                    call read_xml_integer(info, tag, endtag, attribs, noattribs, innerdata, nodata,KSL, has_var)  
!               case('StandardDepthes') 
!                    ALLOCATE(DPTHSL(KSL))
!                    ALLOCATE(p_real_array(KSL))
!                    call read_xml_real_array(info, tag, endtag, attribs, noattribs, innerdata, nodata,p_real_array, has_var) 
!                    DPTHSL=p_real_array
!                    deallocate(p_real_array)
!             end select
!
!             !=========Atmospheric Forcing=================!
!            select case( tag )
!                case('SurfaceHeatingType') 
!                    call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,H_TYPE, has_var)
!                case('MetoForcingType')  
!                    call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,M_TYPE, has_var)
!                case('WindType') 
!                    call read_xml_word(info, tag, endtag, attribs, noattribs, innerdata, nodata,WINDTYPE, has_var)
!                case('VetcHeating') 
!                    call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,RHEAT, has_var)
!                case('VetcHeatingLenScale1') 
!                    call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,ZETA1, has_var)        
!               case('VetcHeatingLenScale2') 
!                    call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,ZETA2, has_var)       
!               case('TimesHeatingInitiated') 
!                    call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,THOUR_HS, has_var)     
!                case('HeightOfWindSpeed') 
!                    call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,ZU, has_var)        
!               case('HeightOfAirTemp') 
!                    call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,ZT, has_var)       
!               case('HeightOfRelativeHumdt') 
!                    call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,ZQ, has_var)                           
!               end select    
!
!        !=========Lagrangian Tracking==================!        
!        !select case( tag )
!        !     case('LagrgTrackOn') 
!        !            call read_xml_logical(info, tag, endtag, attribs, noattribs, innerdata, nodata,LAG_ON, has_var)
!        !     case('LagrgInteval') 
!        !            call read_xml_real(info, tag, endtag, attribs, noattribs, innerdata, nodata,LAG_INTERVAL, has_var)
!        !     case('LagTrackOn') 
!        !            call read_xml_logical(info, tag, endtag, attribs, noattribs, innerdata, nodata,LAG_ON, has_var)     
!        !  end select
!            
!      nodata = 0
!      if ( .not. xml_ok(info) ) then
!        exit
!      endif
!   end do
!   
!       MEDM_ON=.FALSE.
!       AVGE_ON  = .FALSE.
!       Probe_Mode='Custom'
!
!        !=========Parameters Controlling Output of Average Data Fields=================!
!                                                                                
!        AVGE_ON  = .False.
!        INT_AVGE = 30 !18000
!        BEG_AVGE = 30 !72000
!        NUM_AVGE = 30 !2 
!       ! CDF_OUT_AVE = .False.
!        !CDF_VDP_AVE = el ua va
!
!        !=========Parameters Controlling Tidal Forcing==================================
!
!        S_TYPE = 'non-julian'
!        DELTT   = 720.
!
!        !=========Standard Depth Levels=================================================
!
!        KSL    =  11
!        !DPTHSL = 0.05 0.5 0.88 1.27 1.66 2.05 2.44 2.83 3.22 3.31 5
!
!
!        !============ Parameter Controlling Vertical Coordinate Distribution============
!        INDEX_VERCOR = 1
!
!        !============ Parameters Controlling Sigma Distribution=========================
!        !FOR CASE 1 (INDEX_VERCOR = 1). For detail, see set_sigma.F
!        P_SIGMA = 1.0
!
!        !============ The number of vertical layers ====================================
!        KB      = 3
!
!        !============ Parameters Controlling General Vertical Coordinate Distribution===
!        !FOR CASE 2 (INDEX_VERCOR = 2). For detail, see set_sigma.F
!        DU2 = 0.001
!        DL2 = 0.001
!
!        !FOR CASE 3 (INDEX_VERCOR = 3). For detail, see set_sigma.F
!        DUU   = 4.0
!        DLL   = 0.0
!        HMIN1 = 6.0
!        KU    = 2
!        KL    = 0
!
!        ALLOCATE(ZKU(KU)); ZKU=2.0_SP
!        ALLOCATE(ZKL(KL)); ZKL=2.0_SP
!        
!        !============ Parameters Controlling Time Series Output=========================
!
!      !  PROBE_ON = .False.
!
!        !============ Parameters Controlling Water Quality Module=======================
!
!        NB     = 8
!        BENWQM_KEY = .False.
!        
!        !============ Parmaeters Controlling SST Data Assimilation======================
!
!        !SST_ASSIM =.False. 
!        !SST_METHOD = 'OI'
!        !RAD_SST = 10000.
!        !GAMA_SST = 1.0
!        !GALPHA_SST = 3.e-3
!        !ASTIME_WINDOW_SST = .50
!        !N_INFLU_SST = 1
!        !IAV_DAY = 5
!
!        !====CURRENT DATA ASSIMILATION VARIABLES========================================
!
!        !CURRENT_ASSIM =.False. 
!        !CURRENT_METHOD = 'NG'
!        !RAD_CUR = 20000.
!        !GAMA_CUR = 1.0 
!        !GALPHA_CUR = 8.3e-3 
!        !ASTIME_WINDOW_CUR = .50
!        !N_INFLU_CUR = 1 
!
!        !====TEMPERATURE/SALINITY DATA ASSIMILATION VARIABLES===========================
!
!        !TS_ASSIM = .False.
!        !TS_METHOD = 'OI'
!        !RAD_TS = 30000.
!        !GAMA_TS = 1.0
!        !GALPHA_TS = 8.3e-3
!        !ASTIME_WINDOW_TS = 72.
!        !N_INFLU_TS = 1
!        
!        !==== Parameter Controlling Richardson # dep. dissipation correction============
!        SURFACEWAVE_MIX =.False.
!
!        !==== Parameter Controlling Open Boundary Temp/Salt Nudging=====================
!        TS_NUDGING_OBC = .False.
!        ALPHA_OBC = 0.
!
!        !==== Parameter controlling Open Boundary Temp/Salt Series Nudging===========
!        !TSOBC_ON = .False.
!        !ALPHA_SERIES_OBC = .0014
!
!        !=====VARIABLES CONTROLLING 2D MOMENTUM BALANCE CHECKING OUTOUT======
!        !OUT_BALANCE = .False.                  
!        !NUM_BALANCE = 4           !!sum of cell for 2d momentum balance output
!        !NO_CELL  =  11446  11212  11213  11447
!
!        !=====PARAMETER CONTROLLING THE TYPE OF Temp/Salt OBC=======
!        TYPE_TSOBC = 3                  
!
!        !==== Parameter controlling Tidal Open Boundary Output===========
!        !TIDE_INITIAL  = 14400
!        !TIDE_INTERVAL = 6
!
!        !==== Option for semi-implicit corriolis force
!!        ADCOR_ON =.False.
!        
!        !=====VARIABLES for SPECIFY DYE RELEASE                 
!        !DYE_ON = .False.
!        !IINT_SPE_DYE_B = 390962    !391051    
!        !IINT_SPE_DYE_E = 390970    !391140   
!        !KSPE_DYE = 15
!        !!K_SPECIFY = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 
!        !MSPE_DYE = 10
!        !!M_SPECIFY = 17536 17149 17150 16752 16751 16345 16344 15954 15953 15556 
!        !DYE_SOURCE_TERM = 1.0
!
!        !==Parameters of tidal amplitude and phase for one-D mode with biological model
!        !UMAX = 0.0 0.0 0.0 0.0 0.0 0.0
!        !PMAX = 0.0 0.0 0.0 0.0 0.0 0.0
!        
!        !==== Parameter Controlling Net Heat Flux Calculation in Model===========
!        C_HFX = .True.
!        
!        !==== Implicit Factor ===================================
!        !IFCETA = 0.55
!        !BEDF = 1.0
!        !KSTAGE_UV = 0
!        !KSTAGE_TE = 0
!        !KSTAGE_TS = 0
!        !MSTG = slow
!        
!        !==============================================================================|
!        !            SET PHYSICAL PARAMETERS                                           !
!        !==============================================================================|
!
!       KBM1 = KB-1 ; KBM2 = KB-2 
!       DTI=DTE*FLOAT(ISPLIT)
!       IEND=NSTEPS
!       DAYS=NSTEPS*DTI/24.0_SP/3600.0_SP
!       IINT = 0
!       IF(IREPORT == 0) IREPORT = IEND+2
!
!!==============================================================================|
!!            ERROR CHECKING                                                    !
!!==============================================================================|
!
!   IF(RESTART /= 'cold_start' .AND. RESTART /= 'hot_cold_s' .AND. &
!         RESTART /= 'hot_start')THEN
!     IF(MSR)WRITE(IPT,*) 'RESTART NOT CORRECT --->',RESTART   
!     IF(MSR)WRITE(IPT,*) 'SHOULD BE "cold_start","hot_cold_s", or "hot_start"'
!     STOP
!   END IF
!   IF(HORZMIX /= 'constant' .AND. HORZMIX /= 'closure') THEN
!     IF(MSR)WRITE(IPT,*) 'HORZMIX NOT CORRECT --->',HORZMIX  
!     IF(MSR)WRITE(IPT,*) 'SHOULD BE "constant" or "closure"'
!     STOP
!   END IF
!   IF(C_BAROPG /= 'sigma' .AND. C_BAROPG /= 's_levels') THEN
!     IF(MSR)WRITE(IPT,*) 'C_BAROPG NOT CORRECT --->',C_BAROPG
!     IF(MSR)WRITE(IPT,*) 'SHOULD BE "sigma" or "s_levels"'
!     STOP
!   END IF
!   IF(CTRL_DEN /= 'sigma-t' .AND. CTRL_DEN /= 'pdensity' .AND. &
!      CTRL_DEN /= 'sigma-t_stp') THEN
!     IF(MSR)WRITE(IPT,*) 'CTRL_DEN NOT CORRECT,--->',CTRL_DEN
!     IF(MSR)WRITE(IPT,*) 'SHOULD BE "sigma-t" , "pdensity", or "sigma-t_stp"'
!     STOP
!   END IF
!   IF(WINDTYPE /= 'stress' .AND. WINDTYPE /= 'speed') THEN
!     IF(MSR)WRITE(IPT,*) 'WINDTYPE NOT CORRECT,--->',WINDTYPE
!     IF(MSR)WRITE(IPT,*) 'SHOULD BE "stress" or "speed"'
!     STOP 
!   END IF
!   IF(M_TYPE /= 'uniform' .AND. M_TYPE /= 'non-uniform') THEN
!     IF(MSR)WRITE(IPT,*) 'M_TYPE NOT CORRECT,--->',M_TYPE
!     IF(MSR)WRITE(IPT,*) 'SHOULD BE "uniform" or "non-uniform"'
!     STOP 
!   END IF
!   IF(S_TYPE /= 'julian' .AND. S_TYPE /= 'non-julian') THEN
!     IF(MSR)WRITE(IPT,*) 'S_TYPE NOT CORRECT,--->',S_TYPE
!     IF(MSR)WRITE(IPT,*) 'SHOULD BE "julian" or "non-julian"'
!     STOP 
!   END IF
!   IF(BROUGH_TYPE /= 'orig' .AND. BROUGH_TYPE /= 'gotm' .AND. &
!      BROUGH_TYPE /= 'user_defined') THEN
!     IF(MSR)WRITE(IPT,*) 'BROUGH_TYPE NOT CORRECT,--->',BROUGH_TYPE
!     IF(MSR)WRITE(IPT,*) 'SHOULD BE "orig" or "gotm" or "user_defined"'
!     STOP 
!   END IF
!   IF(H_TYPE /= 'body_h' .AND. H_TYPE /= 'flux_h') THEN
!     IF(MSR)WRITE(IPT,*) 'H_TYPE NOT CORRECT,--->',H_TYPE
!     IF(MSR)WRITE(IPT,*) 'SHOULD BE "body_h" or "flux_h"'
!     STOP 
!   END IF
!
!   IF(KB > 200)THEN
!     WRITE(IPT,*)'KB EXCEEDS 200'
!     WRITE(IPT,*)'THIS WILL CAUSE ERROR IN SEVERAL READ STATEMENTS SINCE'
!     WRITE(IPT,*)'ASSOCIATED FORMAT STATEMENT ASSUMES MAX KB OF 200'
!     WRITE(IPT,*)'GREP CODE FOR READ AND 200 TO SEE'
!     call PSTOP
!   END IF
!    
!   !==============================================================================|
!!            REPORTING                                                         !
!!==============================================================================|
!   IF(MSR)WRITE(IPT,*)'!  # STD SALINITY LEVELS :',KSL
!   IF(MSR)WRITE(IPT,*)'!  # OF SIGMA LEVELS     :',KB
!
!!==============================================================================|
!!            SCREEN REPORT OF SET VARIABlES                                    !
!!==============================================================================|
!   IF(MSR)THEN
!   WRITE(IPT,*)'!  # DTE                 :',DTE
!   WRITE(IPT,*)'!  # ISPLIT              :',ISPLIT
!   WRITE(IPT,*)'!  # IRAMP               :',IRAMP
!   WRITE(IPT,*)'!  # NSTEPS              :',NSTEPS
!   WRITE(IPT,*)'!  # IRHO_MEAN           :',IRHO_MEAN
!   WRITE(IPT,*)'!  # RESTART             :',TRIM(RESTART)
!   WRITE(IPT,*)'!  # BFRIC               :',BFRIC
!   WRITE(IPT,*)'!  # MIN_DEPTH           :',MIN_DEPTH
!   WRITE(IPT,*)'!  # Z0B                 :',Z0B   
!   WRITE(IPT,*)'!  # HORZMIX             :',TRIM(HORZMIX)
!   WRITE(IPT,*)'!  # HORCON              :',HORCON
!   WRITE(IPT,*)'!  # HPRNU               :',HPRNU
!   WRITE(IPT,*)'!  # VERTMIX             :',TRIM(VERTMIX)
!   WRITE(IPT,*)'!  # UMOL                :',UMOL
!   WRITE(IPT,*)'!  # VPRNU               :',VPRNU
!   WRITE(IPT,*)'!  # C_BAROPG            :',TRIM(C_BAROPG)
!   WRITE(IPT,*)'!  # CTRL_DEN            :',TRIM(CTRL_DEN)
!   WRITE(IPT,*)'!  # H_TYPE              :',TRIM(H_TYPE   )
!   WRITE(IPT,*)'!  # WINDTYPE            :',TRIM(WINDTYPE   )
!   WRITE(IPT,*)'!  # DJUST               :',DJUST   
!   WRITE(IPT,*)'!  # ZETA1               :',ZETA1   
!   WRITE(IPT,*)'!  # ZETA2               :',ZETA2   
!   WRITE(IPT,*)'!  # RHEAT               :',RHEAT   
!   WRITE(IPT,*)'!  # THOUR_HS            :',THOUR_HS           
!   WRITE(IPT,*)'!  # IRESTART            :',IRESTART   
!   WRITE(IPT,*)'!  # IREPORT             :',IREPORT    
!   WRITE(IPT,*)'!  # IRECORD             :',IRECORD    
!   WRITE(IPT,*)'!  # IDMPSMS             :',IDMPSMS    
!   WRITE(IPT,*)'!  # KSL                 :',KSL      
!   WRITE(IPT,*)'!  # DPTHSL              :',DPTHSL   
!   WRITE(IPT,*)'!  # P_SIGMA             :',P_SIGMA   
!   WRITE(IPT,*)'!  # KB                  :',KB   
!   WRITE(IPT,*)'!  # DELTT               :',DELTT   
!   WRITE(IPT,*)'!  # CASETITLE           :',TRIM(CASETITLE)
!   WRITE(IPT,*)'!  # M_TYPE              :',TRIM(M_TYPE   )
!   WRITE(IPT,*)'!  # S_TYPE              :',TRIM(S_TYPE   )
!   WRITE(IPT,*)'!  # BROUGH_TYPE         :',TRIM(BROUGH_TYPE)
!   WRITE(IPT,*)'!  # OUTDIR              :',TRIM(OUTDIR   )
!   WRITE(IPT,*)'!  # INPDIR              :',TRIM(INPDIR   )
!   WRITE(IPT,*)'!  # INFOFILE            :',TRIM(INFOFILE )
!   IF(AVGE_ON)THEN
!     WRITE(IPT,*)'!  # FLOW AVGES OUTPUT   :  ACTIVE'
!     WRITE(IPT,*)'!  # START ITERATION     :',BEG_AVGE
!     WRITE(IPT,*)'!  # AVGING INTERVAL     :',INT_AVGE
!     WRITE(IPT,*)'!  # NUMBER OF INTERVALS :',NUM_AVGE
!   ELSE
!     WRITE(IPT,*)'!  # FLOW AVGES OUTPUT   :  INACTIVE'
!   END IF
!   IF(VERT_STAB)THEN
!     WRITE(IPT,*)'!  # CONVECTIVE OVERTURN :  ACTIVE'
!   ELSE
!     WRITE(IPT,*)'!  # CONVECTIVE OVERTURN :  INACTIVE'
!   END IF
!   IF(TS_FCT)THEN
!     WRITE(IPT,*)'!  # TEMP/SAL AVERAGING  :  ACTIVE'
!   ELSE
!     WRITE(IPT,*)'!  # TEMP/SAL AVERAGING  :  INACTIVE'
!   END IF
!   IF(BAROTROPIC)THEN
!     WRITE(IPT,*)'!  # BAROTROPIC RUN      :  ACTIVE'
!   END IF
!   IF(TEMP_ON .AND. .NOT. BAROTROPIC)THEN
!     WRITE(IPT,*)'!  # TEMPERATURE EQUATION:  ACTIVE'
!   END IF
!   IF(SALINITY_ON .AND. .NOT. BAROTROPIC)THEN
!     WRITE(IPT,*)'!  # SALINITY EQUATION   :  ACTIVE'
!   END IF
!   IF(SURFACEWAVE_MIX)THEN
!     WRITE(IPT,*)'!  # SURFACE WAVE MIXING :  ACTIVE'
!   END IF
!   IF(TS_NUDGING_OBC)THEN
!     WRITE(IPT,*)'!  # OBC TS NUDGING      :  ACTIVE'
!     WRITE(IPT,*)'!  # NUDGING  COEFF      :',ALPHA_OBC
!   END IF
!   
!   IF(SEDIMENT_ON)THEN
!     WRITE(IPT,*)'!  # SEDIMENT MODEL      :  ACTIVE'
!   END IF
!   IF(RESTART_SED)THEN
!     WRITE(IPT,*)'!  # SEDIMENT MODEL      :  HOT_STARTED'
!  ELSE
!     WRITE(IPT,*)'!  # SEDIMENT MODEL      :  COLD_STARTED'
!  END IF
!
!END IF
!
!
!!==============================================================================|
!!            FORMATS                                                           |
!!==============================================================================|
! 101  FORMAT(A10," = ",F10.6)
! 102  FORMAT(A10," = ",I10)
! 103  FORMAT(A10," = ",A25)
!1000  FORMAT (80a1)
!4000  FORMAT (3i10,1x,a10)
!5000  FORMAT (a10,2e10.3)
!6000  FORMAT (3(2x,a8),4x,a6)
!
!   RETURN
!    END SUBROUTINE ReadPrjXmlFile   
!
!    END MODULE
