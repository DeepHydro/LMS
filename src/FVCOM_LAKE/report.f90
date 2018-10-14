
!==============================================================================|
!     REPORT INITIAL INFORMATION                                               |
!==============================================================================|

   SUBROUTINE REPORT(INFO_STRING)

!==============================================================================|
   USE ALL_VARS
   USE MOD_CLOCK
#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE
   CHARACTER(LEN=*) :: INFO_STRING     !!INFORMATION STRING
   INTEGER :: E3TOT,ESTOT,IERR
   REAL(DP), DIMENSION(17) :: SBUF,RBUF1,RBUF2,RBUF3
   
   REAL(SP), ALLOCATABLE :: Q21(:,:),Q2L1(:,:),L1(:,:)
   REAL(SP), ALLOCATABLE :: KH1(:,:),KQ1(:,:)
   CHARACTER(LEN=13) :: TSTRING
#  if defined (GOTM)
   REAL(SP), ALLOCATABLE :: TKE1(:,:),TEPS1(:,:)
#  endif   
   INTEGER :: I,J,K
   
!==============================================================================|
   OPEN(IOLOG,FILE=TRIM(OUTDIR)//'_report.txt',POSITION='APPEND',STATUS='UNKNOWN')

   ALLOCATE(Q21(1:N,KBM1));   Q21   = 0.0_SP
   ALLOCATE(Q2L1(1:N,KBM1));  Q2L1  = 0.0_SP
   ALLOCATE(L1(1:N,KBM1));    L1    = 0.0_SP
   ALLOCATE(KH1(1:N,KBM1));   KH1   = 0.0_SP
   ALLOCATE(KQ1(1:N,KBM1));   KQ1   = 0.0_SP
#  if defined (GOTM)
   ALLOCATE(TKE1(1:N,KBM1));  TKE1  = 0.0_SP
   ALLOCATE(TEPS1(1:N,KBM1)); TEPS1 = 0.0_SP
#  endif   
   
   DO K=1,KBM1
     DO I=1,N
       DO J=1,3
         Q21(I,K)  = Q21(I,K)+Q2(NV(I,J),K)
         Q2L1(I,K) = Q2L1(I,K)+Q2L(NV(I,J),K)
         L1(I,K)   = L1(I,K)+L(NV(I,J),K)
         KH1(I,K)  = KH1(I,K)+KH(NV(I,J),K)
         KQ1(I,K)  = KQ1(I,K)+KQ(NV(I,J),K)
#  if defined (GOTM)
         TKE1(I,K) = TKE1(I,K)+TKE(NV(I,J),K)
	 TEPS1(I,K)= TEPS1(I,K)+TEPS(NV(I,J),K)
#  endif   
       END DO
       Q21(I,K)  = Q21(I,K)/3.0_SP
       Q2L1(I,K) = Q2L1(I,K)/3.0_SP
       L1(I,K)   = L1(I,K)/3.0_SP
       KH1(I,K)  = KH1(I,K)/3.0_SP
       KQ1(I,K)  = KQ1(I,K)/3.0_SP
#  if defined (GOTM)
         TKE1(I,K) = TKE1(I,K)/3.0_SP
	 TEPS1(I,K)= TEPS1(I,K)/3.0_SP
#  endif   
     END DO
   END DO    	    
      
   SBUF(1)  = SUM(DBLE(UA(1:N)))
   SBUF(2)  = SUM(DBLE(VA(1:N)))
   SBUF(3)  = SUM(DBLE(EL1(1:N)))
   SBUF(4)  = SUM(DBLE(H1(1:N)))
   SBUF(5)  = SUM(DBLE(U(1:N,1:KBM1)))
   SBUF(6)  = SUM(DBLE(V(1:N,1:KBM1)))
   SBUF(7)  = SUM(DBLE(S(1:N,1:KBM1)))
   SBUF(8)  = SUM(DBLE(T(1:N,1:KBM1)))
#  if defined (GOTM)
   SBUF(9)  = SUM(DBLE(TKE1(1:N,2:KBM1)))
   SBUF(10) = SUM(DBLE(TEPS1(1:N,2:KBM1)))
#  else
   SBUF(9)  = SUM(DBLE(Q21(1:N,2:KBM1)))
   SBUF(10) = SUM(DBLE(Q2L1(1:N,2:KBM1)))
#  endif
   SBUF(11) = SUM(DBLE(L1(1:N,1:KBM1)))
   SBUF(12) = SUM(DBLE(KM1(1:N,1:KBM1)))
   SBUF(13) = SUM(DBLE(KQ1(1:N,1:KBM1)))
   SBUF(14) = SUM(DBLE(KH1(1:N,1:KBM1)))
   SBUF(15) = SUM(DBLE(RHO(1:N,1:KBM1)))
   SBUF(16) = SUM(DBLE(D1(1:N)))
#  if !defined (WET_DRY)
   SBUF(17) = FLOAT(N)
#  else
   SBUF(17) = SUM(ISWETC(1:N))
#  endif


   RBUF1 = SBUF


   SBUF(1)  = MAXVAL(UA(1:N))
   SBUF(2)  = MAXVAL(VA(1:N))
   SBUF(3)  = MAXVAL(EL(1:M))
   SBUF(4)  = MAXVAL(H(1:M))
   SBUF(5)  = MAXVAL(U(1:N,1:KBM1))
   SBUF(6)  = MAXVAL(V(1:N,1:KBM1))
   SBUF(7)  = MAXVAL(S1(1:M,1:KBM1))
   SBUF(8)  = MAXVAL(T1(1:M,1:KBM1))
#  if defined (GOTM)
   SBUF(9)  = MAXVAL(TKE(1:M,2:KBM1))
   SBUF(10) = MAXVAL(TEPS(1:M,2:KBM1))
#  else
   SBUF(9)  = MAXVAL(Q2(1:M,1:KBM1))
   SBUF(10) = MAXVAL(Q2L(1:M,1:KBM1))
#  endif
   SBUF(11) = MAXVAL(L(1:M,1:KBM1))
   SBUF(12) = MAXVAL(KM(1:M,1:KBM1))
   SBUF(13) = MAXVAL(KQ(1:M,1:KBM1))
   SBUF(14) = MAXVAL(KH(1:M,1:KBM1))
   SBUF(15) = MAXVAL(RHO1(1:M,1:KBM1))
   SBUF(16) = MAXVAL(D(1:M))

   RBUF2 = SBUF
# if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_REDUCE(SBUF,RBUF2,16,MPI_DP,MPI_MAX,0,MPI_COMM_WORLD,IERR)
# endif

   SBUF(1)  = MINVAL(UA(1:N))
   SBUF(2)  = MINVAL(VA(1:N))
   SBUF(3)  = MINVAL(EL(1:M))
   SBUF(4)  = MINVAL(H(1:M))
   SBUF(5)  = MINVAL(U(1:N,1:KBM1))
   SBUF(6)  = MINVAL(V(1:N,1:KBM1))
   SBUF(7)  = MINVAL(S1(1:M,1:KBM1))
   SBUF(8)  = MINVAL(T1(1:M,1:KBM1))
#  if defined (GOTM)
   SBUF(9)  = MINVAL(TKE(1:M,2:KBM1))
   SBUF(10)  = MINVAL(TEPS(1:M,2:KBM1))
#  else
   SBUF(9)  = MINVAL(Q2(1:M,1:KBM1))
   SBUF(10)  = MINVAL(Q2L(1:M,1:KBM1))
#  endif
   SBUF(11) = MINVAL(L(1:M,1:KBM1))
   SBUF(12) = MINVAL(KM(1:M,1:KBM1))
   SBUF(13) = MINVAL(KQ(1:M,1:KBM1))
   SBUF(14) = MINVAL(KH(1:M,1:KBM1))
   SBUF(15) = MINVAL(RHO1(1:M,1:KBM1))
   SBUF(16) = MINVAL(D(1:M))

   RBUF3 = SBUF
# if defined (MULTIPROCESSOR)
   IF(PAR)CALL MPI_REDUCE(SBUF,RBUF3,16,MPI_DP,MPI_MIN,0,MPI_COMM_WORLD,IERR)
# endif

   IF(MSR)THEN
   IF(LEN_TRIM(INFO_STRING) /= 0)THEN
     CALL GETTIME(TSTRING,INT(DTI*IINT))
     WRITE(IPT,*  )'===================',TRIM(INFO_STRING),'======================'
     WRITE(IPT,102)TSTRING
     WRITE(IOLOG,* )'===================',TRIM(INFO_STRING),'======================'
     WRITE(IOLOG,102 )TSTRING
   END IF
   RBUF1(15) = (RBUF1(15)+NGL*KBM1)*1000.
   RBUF2(15) = (RBUF2(15)+1.)*1000.
   RBUF3(15) = (RBUF3(15)+1.)*1000.
   E3TOT = DBLE(NGL*KBM1)
   ESTOT = DBLE(NGL)
   
if(  Language == 'zh-CN')THEN
   WRITE(IPT,*  )'!变量名称               平均值       最大值     最小值'
   WRITE(IPT,100)'二维外模X方向速度:',RBUF1(1)/ESTOT,RBUF2(1),RBUF3(1)
   WRITE(IPT,100)'二维外模Y方向速度:',RBUF1(2)/ESTOT,RBUF2(2),RBUF3(2)
   WRITE(IPT,100)'自由水面水位:',RBUF1(3)/ESTOT,RBUF2(3),RBUF3(3)
   WRITE(IPT,100)'总水深:',RBUF1(4)/ESTOT,RBUF2(4),RBUF3(4)
   WRITE(IPT,100)'三维内模X方向速度:',RBUF1(5)/E3TOT,RBUF2(5),RBUF3(5)
   WRITE(IPT,100)'三维内模Y方向速度:',RBUF1(6)/E3TOT,RBUF2(6),RBUF3(6)
   WRITE(IPT,100)'污染物浓度:',RBUF1(7)/E3TOT,RBUF2(7),RBUF3(7)
   WRITE(IPT,100)'水温:',RBUF1(8)/E3TOT,RBUF2(8),RBUF3(8)

   WRITE(IOLOG,*  )'!变量名称                平均值      最大值     最小值'
   WRITE(IOLOG,100)'二维外模X方向速度:',RBUF1(1)/ESTOT,RBUF2(1),RBUF3(1)
   WRITE(IOLOG,100)'二维外模Y方向速度:',RBUF1(2)/ESTOT,RBUF2(2),RBUF3(2)
   WRITE(IOLOG,100)'自由水面水位:',RBUF1(3)/ESTOT,RBUF2(3),RBUF3(3)
   WRITE(IOLOG,100)'总水深:',RBUF1(4)/ESTOT,RBUF2(4),RBUF3(4)
   WRITE(IOLOG,100)'三维内模X方向速度:',RBUF1(5)/E3TOT,RBUF2(5),RBUF3(5)
   WRITE(IOLOG,100)'三维内模Y方向速度:',RBUF1(6)/E3TOT,RBUF2(6),RBUF3(6)
   WRITE(IOLOG,100)'污染物浓度:',RBUF1(7)/E3TOT,RBUF2(7),RBUF3(7)
   WRITE(IOLOG,100)'水温:',RBUF1(8)/E3TOT,RBUF2(8),RBUF3(8)
else    
    WRITE(IPT,*  )'!  QUANTITY              :     AVG           MAX         MIN'
   WRITE(IPT,100)'!  EXTERNAL UVEL         :',RBUF1(1)/ESTOT,RBUF2(1),RBUF3(1)
   WRITE(IPT,100)'!  EXTERNAL VVEL         :',RBUF1(2)/ESTOT,RBUF2(2),RBUF3(2)
   WRITE(IPT,100)'!  FREE SURFACE          :',RBUF1(3)/ESTOT,RBUF2(3),RBUF3(3)
   WRITE(IPT,100)'!  BATH DEPTH            :',RBUF1(4)/ESTOT,RBUF2(4),RBUF3(4)
   WRITE(IPT,100)'!  INTERNAL UVEL         :',RBUF1(5)/E3TOT,RBUF2(5),RBUF3(5)
   WRITE(IPT,100)'!  INTERNAL VVEL         :',RBUF1(6)/E3TOT,RBUF2(6),RBUF3(6)
   WRITE(IPT,100)'!  SALINITY              :',RBUF1(7)/E3TOT,RBUF2(7),RBUF3(7)
   WRITE(IPT,100)'!  TEMPERATURE           :',RBUF1(8)/E3TOT,RBUF2(8),RBUF3(8)

   WRITE(IOLOG,*  )'!  QUANTITY              :     AVG           MAX         MIN'
   WRITE(IOLOG,100)'!  EXTERNAL UVEL         :',RBUF1(1)/ESTOT,RBUF2(1),RBUF3(1)
   WRITE(IOLOG,100)'!  EXTERNAL VVEL         :',RBUF1(2)/ESTOT,RBUF2(2),RBUF3(2)
   WRITE(IOLOG,100)'!  FREE SURFACE          :',RBUF1(3)/ESTOT,RBUF2(3),RBUF3(3)
   WRITE(IOLOG,100)'!  BATH DEPTH            :',RBUF1(4)/ESTOT,RBUF2(4),RBUF3(4)
   WRITE(IOLOG,100)'!  INTERNAL UVEL         :',RBUF1(5)/E3TOT,RBUF2(5),RBUF3(5)
   WRITE(IOLOG,100)'!  INTERNAL VVEL         :',RBUF1(6)/E3TOT,RBUF2(6),RBUF3(6)
   WRITE(IOLOG,100)'!  SALINITY              :',RBUF1(7)/E3TOT,RBUF2(7),RBUF3(7)
   WRITE(IOLOG,100)'!  TEMPERATURE           :',RBUF1(8)/E3TOT,RBUF2(8),RBUF3(8)
endif  

#  if defined (GOTM)
   WRITE(IPT,100)'!  TURBULENT KE          :',RBUF1(9)/E3TOT,RBUF2(9),RBUF3(9)
   WRITE(IPT,100)'!  TURBULENT DISSIPATION :',RBUF1(10)/E3TOT,RBUF2(10),RBUF3(10)
#  ENDIF

    IF(Language == 'zh-CN')THEN
           WRITE(IPT,100)'紊动能 KE:',RBUF1(9)/E3TOT,RBUF2(9),RBUF3(9)
           WRITE(IPT,100)'紊动能 KE*L:',RBUF1(10)/E3TOT,RBUF2(10),RBUF3(10)
           WRITE(IPT,100)'紊动掺混长度:',RBUF1(11)/E3TOT,RBUF2(11),RBUF3(11)
           WRITE(IPT,100)'KM:',RBUF1(12)/E3TOT,RBUF2(12),RBUF3(12)
           WRITE(IPT,100)'KQ:',RBUF1(13)/E3TOT,RBUF2(13),RBUF3(13)
           WRITE(IPT,100)'KH:',RBUF1(14)/E3TOT,RBUF2(14),RBUF3(14)
           WRITE(IPT,100)'水体密度:',RBUF1(15)/E3TOT,RBUF2(15),RBUF3(15)
           WRITE(IPT,100)'水深:',RBUF1(16)/ESTOT,RBUF2(16),RBUF3(16)

           WRITE(IOLOG,100)'紊动能 KE:',RBUF1(9)/E3TOT,RBUF2(9),RBUF3(9)
           WRITE(IOLOG,100)'紊动能 KE*L:',RBUF1(10)/E3TOT,RBUF2(10),RBUF3(10)
           WRITE(IOLOG,100)'紊动掺混长度:',RBUF1(11)/E3TOT,RBUF2(11),RBUF3(11)
           WRITE(IOLOG,100)'KM:',RBUF1(12)/E3TOT,RBUF2(12),RBUF3(12)
           WRITE(IOLOG,100)'KQ:',RBUF1(13)/E3TOT,RBUF2(13),RBUF3(13)
           WRITE(IOLOG,100)'KH:',RBUF1(14)/E3TOT,RBUF2(14),RBUF3(14)
           WRITE(IOLOG,100)'水体密度:',RBUF1(15)/E3TOT,RBUF2(15),RBUF3(15)
           WRITE(IOLOG,100)'水深:',RBUF1(16)/ESTOT,RBUF2(16),RBUF3(16)


    ELSE 
          WRITE(IPT,100)'!  TURBULENT KE          :',RBUF1(9)/E3TOT,RBUF2(9),RBUF3(9)
           WRITE(IPT,100)'!  TURB KE*L             :',RBUF1(10)/E3TOT,RBUF2(10),RBUF3(10)

           WRITE(IOLOG,100)'!  TURBULENT KE          :',RBUF1(9)/E3TOT,RBUF2(9),RBUF3(9)
           WRITE(IOLOG,100)'!  TURB KE*L             :',RBUF1(10)/E3TOT,RBUF2(10),RBUF3(10)

           WRITE(IOLOG,100)'!  TURB LENGTH SCALE     :',RBUF1(11)/E3TOT,RBUF2(11),RBUF3(11)
           WRITE(IOLOG,100)'!  KM                    :',RBUF1(12)/E3TOT,RBUF2(12),RBUF3(12)
           WRITE(IOLOG,100)'!  KQ                    :',RBUF1(13)/E3TOT,RBUF2(13),RBUF3(13)
           WRITE(IOLOG,100)'!  KH                    :',RBUF1(14)/E3TOT,RBUF2(14),RBUF3(14)
           WRITE(IOLOG,100)'!  DENSITY               :',RBUF1(15)/E3TOT,RBUF2(15),RBUF3(15)
           WRITE(IOLOG,100)'!  DEPTH                 :',RBUF1(16)/ESTOT,RBUF2(16),RBUF3(16)

            WRITE(IPT,100)'!  TURB LENGTH SCALE     :',RBUF1(11)/E3TOT,RBUF2(11),RBUF3(11)
           WRITE(IPT,100)'!  KM                    :',RBUF1(12)/E3TOT,RBUF2(12),RBUF3(12)
           WRITE(IPT,100)'!  KQ                    :',RBUF1(13)/E3TOT,RBUF2(13),RBUF3(13)
           WRITE(IPT,100)'!  KH                    :',RBUF1(14)/E3TOT,RBUF2(14),RBUF3(14)
           WRITE(IPT,100)'!  DENSITY               :',RBUF1(15)/E3TOT,RBUF2(15),RBUF3(15)
           WRITE(IPT,100)'!  DEPTH                 :',RBUF1(16)/ESTOT,RBUF2(16),RBUF3(16)
ENDIF        

#  if defined (WET_DRY)
   WRITE(IPT,*  )'!  WET/DRY INFO          :   #WET       #DRY             %WET'
   WRITE(IOLOG,*  )'!  WET/DRY INFO          :   #WET       #DRY             %WET'
   IF(RBUF1(17) == FLOAT(NGL))THEN
        WRITE(IPT,*)'!  NO DRY POINTS          '
        WRITE(IOLOG,*)'!  NO DRY POINTS          '
   ELSE
        WRITE(IPT,101)'!  WET/DRY DATA          :',INT(RBUF1(17)),NGL-INT(RBUF1(17)),100.*RBUF1(17)/FLOAT(NGL)
         WRITE(IOLOG,101)'!  WET/DRY DATA          :',INT(RBUF1(17)),NGL-INT(RBUF1(17)),100.*RBUF1(17)/FLOAT(NGL)
   END IF
#  endif

   END IF

   CLOSE(IOLOG,STATUS='KEEP')
   DEALLOCATE(Q21,Q2L1,L1)
   DEALLOCATE(KH1,KQ1)
#  if defined (GOTM)
   DEALLOCATE(TKE1,TEPS1)
#  endif
   
   RETURN


101    FORMAT(1X,A26,2I12,F12.6)
#if defined (EN)       
100 FORMAT(1X,A26,3F12.6)
102 FORMAT('! Current time: ', A14)
#else    
100 FORMAT(A20,3F12.6)
102 FORMAT(' !当前时间: ', A14)
#endif    
   END SUBROUTINE REPORT 
!==============================================================================|

