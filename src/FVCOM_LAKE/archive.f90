!==============================================================================|
!     ARCHIVE THE MODEL RESULTS                                                |
!==============================================================================|

   SUBROUTINE ARCHIVE
   USE ALL_VARS
   USE MOD_WQM
!==============================================================================|

#  if defined (WET_DRY)
   USE MOD_WD
#  endif   

#  if defined (SEDIMENT)
   USE MOD_SED
#  endif
   USE CONTROL  

   IMPLICIT NONE
   CHARACTER(LEN=8)   :: CH8
   CHARACTER(LEN=100) :: COPYFILE
   CHARACTER(LEN=120) :: FNAME
   CHARACTER(LEN=4)   :: CH4
!==============================================================================|

!
!--DUMP MEDM FILE EVERY IRECORD  ITERATIONS (NO DUMP IF IRECORD  = 0)----------!
!
IF(MEDM_ON)THEN
     IF(IRECORD /= 0)THEN
         IF(MOD(IINT,IRECORD) == 0)THEN
            IF(OUT_FORMAT == 'BINARY') THEN
                 CALL OUT_BINARY(IINT/IRECORD)
              ELSE IF (OUT_FORMAT == 'TXT') THEN
                CALL OUT_TXT(IINT/IRECORD)
              ENDIF
         END IF
     END IF
ENDIF

IF(Output_ON)THEN
     IF(IRECORD /= 0)THEN
         IF(MOD(IINT,IRECORD) == 0)THEN
            IF(OUT_FORMAT == 'BINARY') THEN
                    WRITE(IPT,*)  '!  DUMPING               :    BINARY FILE'    
                        CALL out_binary_cfd(0,IINT/IDMPSMS)
              ELSE IF (OUT_FORMAT == 'TXT') THEN
                     CALL OUT_TXT(IINT/IRECORD)
              ENDIF
         END IF
     END IF
ENDIF
!
!----------------------DUMP SMS FILE-------------------------------------------!
!
IF(SMS_ON) then
   IF(IDMPSMS /= 0)THEN
     IF(MOD(IINT,IDMPSMS) == 0)THEN
       IF(MSR) THEN
            WRITE(IPT,*)  '!  DUMPING               :    SMS FILE'    
                CALL OUT_SMS_ONE(IINT/IDMPSMS)
            IF(WQM_ON)THEN        
                CALL OUT_SMS_WQ(IINT/IDMPSMS)
            ENDIF
       ENDIF
     ENDIF
   END IF
ENDIF
   
!--DUMP RESTART FILE EVERY IRESTART ITERATIONS (NO DUMP IF IRESTART = 0)-------!
   IF(IRESTART /= 0)THEN
     IF(MOD(IINT,IRESTART) == 0)THEN
       IF(MSR)WRITE(IPT,*)  '!  DUMPING               :    RESTART FILE'
       CALL ARCRST

#      if defined(WET_DRY)
       IF(WET_DRY_ON) THEN
          WRITE(CH8,'(I8.8)') IINT
	       FNAME = 're_'//trim(CH8)//'_wd'
          CALL WD_DUMP(FNAME)
       ENDIF
#      endif       
#      if defined(SEDIMENT)
       IF(SEDIMENT_ON) CALL ARCHIVE_SED
#      endif       
     END IF
   END IF

!
!--FLOW FIELD AVERAGING AND OUTPUT OF FLOWFIELD AVERAGES-----------------------!
!
   IF(AVGE_ON) THEN
     IF(SMS_AVGE_ON) THEN 
        CALL out_sms_avge
     ELSE 
       CALL OUT_AVGE
     ENDIF
   ENDIF

   RETURN
   END SUBROUTINE ARCHIVE
!==============================================================================|
