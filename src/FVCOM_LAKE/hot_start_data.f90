
!==============================================================================|
!   READ IN RESTART DATA FILE AND RESTART                                      |
!==============================================================================|

   SUBROUTINE HOT_START_DATA      

!------------------------------------------------------------------------------|

   USE ALL_VARS

#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
#  if defined (DYE_RELEASE)
   USE MOD_DYE
#  endif   

#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
#  if defined (EQUI_TIDE)
   USE MOD_EQUITIDE
#  endif
#  if defined (ATMO_TIDE)
   USE MOD_ATMOTIDE
#  endif

#  if defined (ENKF_ASSIM)
   USE MOD_ENKF
#  endif

   IMPLICIT NONE
   INTEGER :: I,K,N1
   REAL(SP), DIMENSION(N) ::RTP
   CHARACTER CH4*4 
!==============================================================================|
!  NOTE: TO MAINTAIN COMPATIBILITY WITH PREVIOUS DOTFVM, ARRAY DTF1 (NO LONGER  !
!  USED) IS READ IN AS RTP                                                     !
!==============================================================================|

# if defined (ENKF_ASSIM)
   IF(IENS>=1) THEN    ! need change
      WRITE(CH4,'(i4.4)') IENS
      IF(ICYC>=ENKF_START/DELTA_ASS .AND. ICYC<=ENKF_END/DELTA_ASS+1) THEN
        OPEN(INRES,file=TRIM(OUTDIR)//'/anl/restart'//ch4//'.dat',form='unformatted',status='old')
      ELSE
        OPEN(INRES,file=TRIM(OUTDIR)//'/fct/restart'//ch4//'.dat',form='unformatted',status='old')
      ENDIF
   ENDIF   
# endif

   IF(SERIAL)THEN
     REWIND(INRES)
     READ(INRES) IINT
     READ(INRES) ((U(I,K),K=1,KB),I=0,N)
     READ(INRES) ((V(I,K),K=1,KB),I=0,N)
     READ(INRES) ((W(I,K),K=1,KB),I=0,N)
#    if defined (GOTM)
     READ(INRES) ((TKE(I,K),K=1,KB),I=0,M)
     READ(INRES) ((TEPS(I,K),K=1,KB),I=0,M)
#    else
     READ(INRES) ((Q2(I,K),K=1,KB),I=0,M)
     READ(INRES) ((Q2L(I,K),K=1,KB),I=0,M)
     READ(INRES) ((L(I,K),K=1,KB),I=0,M)
#    endif
     READ(INRES) ((S(I,K),K=1,KB),I=0,N)
     READ(INRES) ((T(I,K),K=1,KB),I=0,N)
     READ(INRES) ((RHO(I,K),K=1,KB),I=0,N)
     READ(INRES) ((TMEAN(I,K),K=1,KB),I=0,N)
     READ(INRES) ((SMEAN(I,K),K=1,KB),I=0,N)
     READ(INRES) ((RMEAN(I,K),K=1,KB),I=0,N)

     READ(INRES) ((S1(I,K),K=1,KB),I=1,M)
     READ(INRES) ((T1(I,K),K=1,KB),I=1,M)
     READ(INRES) ((RHO1(I,K),K=1,KB),I=1,M)
     READ(INRES) ((TMEAN1(I,K),K=1,KB),I=1,M)
     READ(INRES) ((SMEAN1(I,K),K=1,KB),I=1,M)
     READ(INRES) ((RMEAN1(I,K),K=1,KB),I=1,M)
     READ(INRES) ((KM(I,K),K=1,KB),I=1,M)
     READ(INRES) ((KH(I,K),K=1,KB),I=1,M)
     READ(INRES) ((KQ(I,K),K=1,KB),I=1,M)

     READ(INRES) (UA(I), I=0,N)
     READ(INRES) (VA(I), I=0,N)
     READ(INRES) (EL1(I), I=1,N)
     READ(INRES) (ET1(I), I=1,N)
     READ(INRES) (H1(I), I=1,N)
     READ(INRES) (D1(I), I=1,N)
     READ(INRES) (DT1(I), I=1,N)
!    READ(INRES) (DTF1(I), I=1,N)
     READ(INRES) (RTP(I), I=1,N)

     READ(INRES) (EL(I), I=1,M)
     READ(INRES) (ET(I), I=1,M)
     READ(INRES) (H(I), I=1,M)
     READ(INRES) (D(I), I=1,M)
     READ(INRES) (DT(I), I=1,M)

#    if defined (EQUI_TIDE)
     READ(INRES) (EL_EQI(I), I=1,M)
#    endif
#    if defined (ATMO_TIDE)
     READ(INRES) (EL_ATMO(I), I=1,M)
#    endif

#    if defined (WATER_QUALITY)
     DO N1=1,NB
       READ(INRES) ((WQM(I,K,N1),K=1,KB),I=1,M)
     END DO
#    endif

#    if defined (DYE_RELEASE)
       IF(IINT.GT.IINT_SPE_DYE_B) THEN
         READ(INRES) ((DYE(I,K),K=1,KB),I=1,M)
         READ(INRES) ((DYEMEAN(I,K),K=1,KB),I=1,M)
       ENDIF
#    endif

     CLOSE(INRES)
   ELSE
#  if defined (MULTIPROCESSOR)
     REWIND(INRES)
     READ(INRES) IINT
     CALL PREAD(INRES,U     ,LBOUND(U,1),    UBOUND(U,1),    N,NGL,KB,EGID(1),0,"U"     )
     CALL PREAD(INRES,V     ,LBOUND(V,1),    UBOUND(V,1),    N,NGL,KB,EGID(1),0,"V"     )
     CALL PREAD(INRES,W     ,LBOUND(W,1),    UBOUND(W,1),    N,NGL,KB,EGID(1),0,"W"     )
#    if defined (GOTM)
     CALL PREAD(INRES,TKE    ,LBOUND(TKE,1),   UBOUND(TKE,1),   M,MGL,KB,NGID(1),0,"TKE"    )
     CALL PREAD(INRES,TEPS   ,LBOUND(TEPS,1),  UBOUND(TEPS,1),  M,MGL,KB,NGID(1),0,"TEPS"   )
#    else
     CALL PREAD(INRES,Q2    ,LBOUND(Q2,1),   UBOUND(Q2,1),   M,MGL,KB,NGID(1),0,"Q2"    )
     CALL PREAD(INRES,Q2L   ,LBOUND(Q2L,1),  UBOUND(Q2L,1),  M,MGL,KB,NGID(1),0,"Q2L"   )
     CALL PREAD(INRES,L     ,LBOUND(L,1  ),  UBOUND(L,1),    M,MGL,KB,NGID(1),0,"L"   )
#    endif
     CALL PREAD(INRES,S     ,LBOUND(S,1),    UBOUND(S,1),    N,NGL,KB,EGID(1),0,"S"     )
     CALL PREAD(INRES,T     ,LBOUND(S,1),    UBOUND(T,1),    N,NGL,KB,EGID(1),0,"T"     )
     CALL PREAD(INRES,RHO   ,LBOUND(RHO,1),  UBOUND(RHO,1),  N,NGL,KB,EGID(1),0,"RHO"   )
     CALL PREAD(INRES,TMEAN ,LBOUND(TMEAN,1),UBOUND(TMEAN,1),N,NGL,KB,EGID(1),0,"TMEAN" )
     CALL PREAD(INRES,SMEAN ,LBOUND(SMEAN,1),UBOUND(SMEAN,1),N,NGL,KB,EGID(1),0,"SMEAN" )
     CALL PREAD(INRES,RMEAN ,LBOUND(RMEAN,1),UBOUND(RMEAN,1),N,NGL,KB,EGID(1),0,"RMEAN" )

     CALL PREAD(INRES,S1    ,LBOUND(S1,1),    UBOUND(S1,1),    M,MGL,KB,NGID,1,"S1"     )
     CALL PREAD(INRES,T1    ,LBOUND(T1,1),    UBOUND(T1,1),    M,MGL,KB,NGID,1,"T1"     )
     CALL PREAD(INRES,RHO1  ,LBOUND(RHO1,1),  UBOUND(RHO1,1),  M,MGL,KB,NGID,1,"RHO1"   )
     CALL PREAD(INRES,TMEAN1,LBOUND(TMEAN1,1),UBOUND(TMEAN1,1),M,MGL,KB,NGID,1,"TMEAN1" )
     CALL PREAD(INRES,SMEAN1,LBOUND(SMEAN1,1),UBOUND(SMEAN1,1),M,MGL,KB,NGID,1,"SMEAN1" )
     CALL PREAD(INRES,RMEAN1,LBOUND(RMEAN1,1),UBOUND(RMEAN1,1),M,MGL,KB,NGID,1,"RMEAN1" )
  
     CALL PREAD(INRES,KM  ,LBOUND(KM,1),UBOUND(KM,1),M,MGL,KB,NGID(1),1,"KM" )
     CALL PREAD(INRES,KH  ,LBOUND(KH,1),UBOUND(KH,1),M,MGL,KB,NGID(1),1,"KH" )
     CALL PREAD(INRES,KQ  ,LBOUND(KQ,1),UBOUND(KQ,1),M,MGL,KB,NGID(1),1,"KQ" )

     CALL PREAD(INRES,UA  ,LBOUND(UA,1), UBOUND(UA,1), N,NGL,1 ,EGID(1),0,"UA"  )
     CALL PREAD(INRES,VA  ,LBOUND(VA,1), UBOUND(VA,1), N,NGL,1 ,EGID(1),0,"VA"  )
     CALL PREAD(INRES,EL1 ,LBOUND(EL1,1),UBOUND(EL1,1),N,NGL,1 ,EGID(1),1,"EL1" )
     CALL PREAD(INRES,ET1 ,LBOUND(ET1,1),UBOUND(ET1,1),N,NGL,1 ,EGID(1),1,"ET1" )
     CALL PREAD(INRES,H1  ,LBOUND(H1,1), UBOUND(H1,1), N,NGL,1 ,EGID(1),1,"H1"  )
     CALL PREAD(INRES,D1  ,LBOUND(D1,1), UBOUND(D1,1), N,NGL,1 ,EGID(1),1,"D1"  )
     CALL PREAD(INRES,DT1 ,LBOUND(DT1,1),UBOUND(DT1,1),N,NGL,1 ,EGID(1),1,"DT1" )
     CALL PREAD(INRES,RTP ,LBOUND(RTP,1),UBOUND(RTP,1),N,NGL,1 ,EGID(1),1,"RTP" )

     CALL PREAD(INRES,EL  ,LBOUND(EL,1),UBOUND(EL,1),M,MGL,1 ,NGID,1,"EL"   )
     CALL PREAD(INRES,ET  ,LBOUND(ET,1),UBOUND(ET,1),M,MGL,1 ,NGID,1,"ET"   )
     CALL PREAD(INRES,H   ,LBOUND(H,1), UBOUND(H,1), M,MGL,1 ,NGID,1,"H"    )
     CALL PREAD(INRES,D   ,LBOUND(D,1), UBOUND(D,1), M,MGL,1 ,NGID,1,"D"    )
     CALL PREAD(INRES,DT  ,LBOUND(DT,1),UBOUND(DT,1),M,MGL,1 ,NGID,1,"DT"   )

#     if defined (EQUI_TIDE)
      CALL PREAD(INRES,EL_EQI,LBOUND(EL_EQI,1),UBOUND(EL_EQI,1),M,MGL,1 ,NGID,1,"EL_EQI"   )
#     endif
#     if defined (ATMO_TIDE)
      CALL PREAD(INRES,EL_ATMO,LBOUND(EL_ATMO,1),UBOUND(EL_ATMO,1),M,MGL,1 ,NGID,1,"EL_ATMO"   )
#     endif

#    if defined (WATER_QUALITY)
     DO N1 = 1, NB
       CALL PREAD(INRES,WQM(:,:,N1),LBOUND(WQM(:,:,N1),1),UBOUND(WQM(:,:,N1),1), &
                  M,MGL,KB,NGID,1,"WQM")
     END DO
#    endif

#    if defined (DYE_RELEASE)
     IF(IINT.GT.IINT_SPE_DYE_B) THEN
     CALL PREAD(INRES,DYE    ,LBOUND(DYE,1),    UBOUND(DYE,1),     M,MGL,KB,NGID,1,"DYE"     )
     CALL PREAD(INRES,DYEMEAN,LBOUND(DYEMEAN,1),UBOUND(DYEMEAN,1), M,MGL,KB,NGID,1,"DYEMEAN" )
     ENDIF
#    endif

     CLOSE(INRES)
#    endif
   END IF

!--Set Turbulent Macro-Scale
#  if defined(GOTM)
   L = .001 
   L(1:M,2:KBM1) = (.5544**3)*TKE(1:M,2:KBM1)**1.5/TEPS(1:M,2:KBM1)
#  endif

   CALL N2E3D(KM,KM1)

   RETURN
   END SUBROUTINE HOT_START_DATA
!==============================================================================|
