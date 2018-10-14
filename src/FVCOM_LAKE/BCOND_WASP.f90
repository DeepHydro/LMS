   SUBROUTINE BCOND_WASP     
!==============================================================================|
!   Set Boundary Conditions on Water Quality                                   |
!==============================================================================|

!------------------------------------------------------------------------------|
   USE MOD_WQM
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
   IMPLICIT NONE
   REAL(SP) :: T2D,T2D_NEXT,T2D_OBC,XFLUX2D,TMP
   INTEGER  :: I,J,K,J1,J11,J22,NCON2,N1
   INTEGER  ::I1,I2,V1,V2,V3,C
   REAL(SP) ::WQMMAX,WQMMIN,WQ_AVR
!------------------------------------------------------------------------------|


!
!--SET CONDITIONS FOR FRESH WATER INFLOW---------------------------------------|
!
   IF(POINT_ST_TYPE == 'specified') THEN
     IF(NUMQBC > 0) THEN
       IF(INFLOW_TYPE == 'node') THEN
         DO I=1,NUMQBC
           J11=INODEQ(I)
           DO K=1,KBM1
             DO N1=1,NB
               WQM_F(J11,K,N1) = WDIS(I,N1)
             END DO
           END DO
         END DO
       ELSE IF(INFLOW_TYPE == 'edge') THEN
         DO I=1,NUMQBC
           J11=N_ICELLQ(I,1)
           J22=N_ICELLQ(I,2)
           DO K=1,KBM1
             DO N1=1,NB
               WQM_F(J11,K,N1)=WDIS(I,N1)
               WQM_F(J22,K,N1)=WDIS(I,N1)
             END DO
           END DO
         END DO
       END IF
     END IF
  ELSE  IF(POINT_ST_TYPE == 'calculated') THEN
     IF(INFLOW_TYPE == 'node') THEN
       IF(NUMQBC > 0) THEN
          DO I=1,NUMQBC
            J11=INODEQ(I)
             IF(QDIS(I)>0)THEN
                   DO K=1,KBM1
                     DO N1=1,NB
                       WQM_F(J11,K,N1) = WDIS(I,N1)
                     END DO
                   END DO
             ELSE
              DO N1=1,NB
                WQ_AVR=0.0_SP
                 DO C=1,NTVE(J11)
                       I1=NBVE(J11,C)
                       V1 = NV(I1,1)
                       V2 = NV(I1,2)
                       V3 = NV(I1,3)
                       WQ_AVR = WQ_AVR + (WQM_F(V1,1,N1)+WQM_F(V2,1,N1)+WQM_F(V3,1,N1))/3.0_SP
                 ENDDO
                  DO K=1,KBM1
                       WQM_F(J11,K,N1) =  WQ_AVR / NTVE(J11)
                   END DO
            ENDDO          
            END IF
          ENDDO
       END IF
    END IF
  END IF

       
   IF(IOBCN > 0) THEN

!
!  SET WATER QUALITY CONDITIONS ON OUTER BOUNDARY
!
     DO N1=1,NB
       DO I=1,IOBCN
         J=I_OBC_N(I)
         J1=NEXT_OBC(I)
         T2D=0.0_SP
         T2D_NEXT=0.0_SP
         XFLUX2D=0.0_SP
         DO K=1,KBM1
           T2D=T2D+WQM(J,K,N1)*DZ(J,K)
           T2D_NEXT=T2D_NEXT+WQM_F(J1,K,N1)*DZ(J1,K)
           XFLUX2D=XFLUX2D+XFLUX_OBC(I,K)                             !*DZ(K)
         END DO
    
         IF(UARD_OBCN(I) > 0.0_SP) THEN
           TMP=XFLUX2D+T2D*UARD_OBCN(I)
           T2D_OBC=(T2D*DT(J)-TMP*DTI/ART1(J))/D(J)
           DO K=1,KBM1
             WQM_F(J,K,N1)=T2D_OBC+(WQM_F(J1,K,N1)-T2D_NEXT)
!             WQM_F(J,K,N1) = WQM_F(J1,K,N1)
           END DO

         DO K=1,KBM1
           WQMMAX = MAXVAL(WQM(NBSN(J,1:NTSN(J)),K,N1))
           WQMMIN = MINVAL(WQM(NBSN(J,1:NTSN(J)),K,N1))
         
           IF(K == 1)THEN
            WQMMAX = MAX(WQMMAX,(WQM(J,K,N1)*DZ(J,K+1)+WQM(J,K+1,N1)*DZ(J,K))/  &
	             (DZ(J,K)+DZ(J,K+1)))
            WQMMIN = MIN(WQMMIN,(WQM(J,K,N1)*DZ(J,K+1)+WQM(J,K+1,N1)*DZ(J,K))/  &
	             (DZ(J,K)+DZ(J,K+1)))
           ELSE IF(K == KBM1)THEN
            WQMMAX = MAX(WQMMAX,(WQM(J,K,N1)*DZ(J,K-1)+WQM(J,K-1,N1)*DZ(J,K))/  &
	             (DZ(J,K)+DZ(J,K-1)))
            WQMMIN = MIN(WQMMIN,(WQM(J,K,N1)*DZ(J,K-1)+WQM(J,K-1,N1)*DZ(J,K))/  &
	             (DZ(J,K)+DZ(J,K-1)))
           ELSE
            WQMMAX = MAX(WQMMAX,(WQM(J,K,N1)*DZ(J,K-1)+WQM(J,K-1,N1)*DZ(J,K))/  &
	             (DZ(J,K)+DZ(J,K-1)), &
                     (WQM(J,K,N1)*DZ(J,K+1)+WQM(J,K+1,N1)*DZ(J,K))/  &
		     (DZ(J,K)+DZ(J,K+1)))
            WQMMIN = MIN(WQMMIN,(WQM(J,K,N1)*DZ(J,K-1)+WQM(J,K-1,N1)*DZ(J,K))/  &
	             (DZ(J,K)+DZ(J,K-1)), &
                     (WQM(J,K,N1)*DZ(J,K+1)+WQM(J,K+1,N1)*DZ(J,K))/  &
		     (DZ(J,K)+DZ(J,K+1)))
           END IF
 
           IF(WQMMIN-WQM_F(J,K,N1) > 0.0_SP)WQM_F(J,K,N1) = WQMMIN
           IF(WQM_F(J,K,N1)-WQMMAX > 0.0_SP)WQM_F(J,K,N1) = WQMMAX

         END DO

         ELSE
           DO K=1,KBM1
               WQM_F(J,K,N1)=WQM(J,K,N1)
           END DO
         END IF
       END DO
     END DO !!OUTER LOOP OVER WQ VARIABLES


   END IF

!
!--SET BOUNDARY CONDITIONS-----------------------------------------------------|
!
   WQM(0,:,:) = 0.0_SP
!   DO K = 1,KBM1
!     DO N1 = 1,NB
!       WQM(0,K,N1) = 0.0_SP
!     END DO
!   END DO

   RETURN
   END SUBROUTINE BCOND_WASP