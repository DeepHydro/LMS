!==============================================================================!
!   Open Input Files for Model Parameters and Output Files for Results         !
!==============================================================================!

   SUBROUTINE IOFILES
   USE ALL_VARS
   USE MOD_UTILS

   IMPLICIT NONE

   LOGICAL:: CHECK 
   CHARACTER(LEN=80)  :: TEMP,ISTR,OSTR
   INTEGER IERR,ISTAT
   
   IF(MSR)WRITE(IPT,*)'!                                                                !'
   IF(MSR)WRITE(IPT,*)'!                  OPENING FILES                                 !'
   IF(MSR)WRITE(IPT,*)'!                                                                !'

   INGRD =12
   INCWH =14
   INRIV =15
   INITS =16
   INJUL =18
   INHFX =19
   INWND =20
   INEVP =21
   INELF =22
   INUVF =23
   INCOR =26
   INBFW =27

   IOPRT =41
   IOPLT =42
   IOTSR =43

   IOSMSD=51
   IOSMSV=52
   IOSMST=53

   IGOTM =59

   INRES=54
   INJMP=55
   IREST=60

!---------CHECK EXISTENCE OF STANDARD FILES-AND OPEN---------------------------!
!
   ISTR = "./"//TRIM(INPDIR)//"/"//trim(casename)
   OSTR = "./"//TRIM(OUTDIR)//"/"//"out/"//trim(casename)
   CALL FOPEN(IOPRT, TRIM(OSTR)//'_prt.dat',"ofr")
   CALL FOPEN(INGRD, TRIM(ISTR)//'_grd.dat',"cfr")

   CALL FOPEN(INRIV, TRIM(ISTR)//'_riv.dat',"cfr")
   CALL FOPEN(INBFW, TRIM(ISTR)//'_bfw.dat',"cfr")

!
!-----------------INITIAL TEMPERATURE AND SALINITY-----------------------------!
!
   IF(RESTART /= 'hot_start') CALL FOPEN(INITS, TRIM(ISTR)//'_its.dat',"cfr")
!
!-----------------OPEN METEOROLOGICAL FORCING FILES----------------------------!
!
   IF(trim(M_TYPE) == 'uniform')THEN
     CALL FOPEN(INCWH, TRIM(ISTR)//'_mc.dat',"cfr")
   ELSE 
     CALL FOPEN(INHFX, TRIM(ISTR)//'_hfx.dat' ,"cur")
     CALL FOPEN(INWND, TRIM(ISTR)//'_wnd.dat' ,"cur")

     EVP_FLAG = .FALSE.
     INQUIRE(FILE=TRIM(ISTR)//'_evp.dat',EXIST=CHECK) 
     IF(CHECK)THEN
       EVP_FLAG = .TRUE.
       CALL FOPEN(INEVP, TRIM(ISTR)//'_evp.dat' ,"cur")
     END IF     
   END IF

!-----------------------FILES FOR ARCHIVING------------------------------------!
!
   CALL FOPEN(IOPLT,TRIM(OSTR)//'_plt.dat' ,"our")
   CALL FOPEN(IOTSR,TRIM(OSTR)//'_tsr.dat' ,"ofr")

!
!-----------------------DEPTH OUTPUT FOR SMS PLOT------------------------------!
!
   CALL FOPEN(IOSMSD,TRIM(OUTDIR)//"/sms/"//trim(casename)//"_dep.xy" ,"ofr")

   RETURN
   END SUBROUTINE IOFILES
!==============================================================================!




