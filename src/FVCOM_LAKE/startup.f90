!==============================================================================|
!   Begin Restart Run From Specified Time                                      |
!==============================================================================|

   SUBROUTINE STARTUP             

   USE ALL_VARS
   USE MOD_WQM

#  if defined (WET_DRY)
   USE MOD_WD
#  endif   
#  if defined (DYE_RELEASE)
   USE MOD_DYE
#  endif   

   USE BCS

   IMPLICIT NONE

   CHARACTER(LEN=120) :: FNAME
   CHARACTER(LEN=8)   :: RRKINP1
   CHARACTER(LEN=4)   :: RRKINP2
   CHARACTER(LEN=4)   :: ENKINP
!==============================================================================|
!--Set Water Depth-Using Bathymetry and Free Surface Elevation-----------------!

   CALL WATER_DEPTH

#  if defined (WET_DRY)
   IF(WET_DRY_ON) CALL SET_WD_DATA
#  endif
!
!--Set up Temperature, Salinity, and Turbulence Quantity Fields----------------!
! 
      
   IF((RESTART == 'cold_start').AND.(S_TYPE == 'non-julian'))THEN

     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    COLD_START'
     IF(MSR)WRITE(IPT,*)  '!  S_TYPE                :    NON-JULIAN'
     CALL INITIAL_TS
#    if defined (WET_DRY)
     IF(WET_DRY_ON) CALL SET_WD_DATA
#    endif
     CALL INITIAL_QQL
#    if defined (WATER_QUALITY)
     CALL INITIAL_WQM
#    endif

#    if defined (DYE_RELEASE)
     CALL INITIAL_DYE
#    endif

   ELSE IF((RESTART=='cold_start').AND.(S_TYPE=='julian'))THEN

     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    COLD_START'
     IF(MSR)WRITE(IPT,*)  '!  S_TYPE                :    JULIAN'
     CALL INITIAL_TS
     CALL INITIAL_UVEL
#    if defined (WET_DRY)
     IF(WET_DRY_ON) CALL SET_WD_DATA
#    endif
     CALL INITIAL_QQL
#    if defined (WATER_QUALITY)
     CALL INITIAL_WQM
#    endif

#    if defined (DYE_RELEASE)
     CALL INITIAL_DYE
#    endif

   ELSE IF((RESTART=='hot_cold_s').AND.(S_TYPE=='julian'))THEN
          
     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    HOT_COLD_S'
     IF(MSR)WRITE(IPT,*)  '!  S_TYPE                :    JULIAN'
     CALL HOT_START_DATA
     CALL INITIAL_TS
#    if defined (WATER_QUALITY)
     CALL INITIAL_WQM
#    endif

#    if defined (WET_DRY)
     FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_restart_wd.dat"
     IF(WET_DRY_ON) CALL WD_READ(FNAME)    
#    endif

#    if defined (DYE_RELEASE)
     CALL INITIAL_DYE
#    endif

   ELSE IF(RESTART == 'hot_start') THEN
     
     IF(MSR)WRITE(IPT,*)  '!  STARTUP TYPE          :    HOT_START'
     CALL HOT_START_DATA
     IF(MSR)WRITE(IPT,*)  '!  RESTART DATA          :    READ     '
#    if defined (WET_DRY)
     FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_restart_wd.dat"
     IF(WET_DRY_ON) CALL WD_READ(FNAME)    
#    endif

#    if defined (DYE_RELEASE)
     CALL INITIAL_DYE
#    endif

   ELSE
         
     PRINT*,'RESTAR AND S_TYPE DEFINITION NOT CORRECT'
     PRINT*,'RESTAR==',RESTART
     PRINT*,'S_TYPE==',S_TYPE
     CALL PSTOP
         
   END IF

   IF(SERIAL)RETURN

   RETURN
   END SUBROUTINE STARTUP
!==============================================================================|

   SUBROUTINE EXCHANGE_ALL 
!------------------------------------------------------------------------------|

   USE ALL_VARS
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif

#  if defined (DYE_RELEASE)
   USE MOD_DYE
#  endif   

   IMPLICIT NONE

   CALL RHO_MEAN

   RETURN
   END SUBROUTINE EXCHANGE_ALL
!==============================================================================|
