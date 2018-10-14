    !==============================================================================|
    !  FLUX CONTROL FOR SALINITY                                                        |
    !==============================================================================|

    SUBROUTINE FCT_WASP

    USE MOD_WQM

    !==============================================================================|
    USE ALL_VARS
    USE BCS
    USE MOD_OBCS
    IMPLICIT NONE
    REAL(SP):: SMAX,SMIN
    INTEGER :: I,J,K,N1
    !==============================================================================|

    DO N1=1,NB
        DO I=1,M
            IF(IOBCN > 0)THEN
                DO J=1,IOBCN
                    IF(I == I_OBC_N(J))GO TO 200
                END DO
            END IF  	 
            IF(NUMQBC > 0)THEN
                DO J=1,NUMQBC
                    IF(INFLOW_TYPE == 'node')THEN
                        IF(I == INODEQ(J))GO TO 200
                    END IF  
                    IF(INFLOW_TYPE == 'edge')THEN
                        IF(I == N_ICELLQ(J,1) .OR. I == N_ICELLQ(J,2))GO TO 200
                    END IF  
                END DO
            END IF
            DO K=1,KBM1
                SMAX = MAXVAL(WQM(NBSN(I,1:NTSN(I)),K,N1))
                SMIN = MINVAL(WQM(NBSN(I,1:NTSN(I)),K,N1))

                IF(K == 1)THEN
                    SMAX = MAX(SMAX,(WQM(I,K,N1)*DZ(I,K+1)+WQM(I,K+1,N1)*DZ(I,K))/  &
                    (DZ(I,K)+DZ(I,K+1)))
                    SMIN = MIN(SMIN,(WQM(I,K,N1)*DZ(I,K+1)+WQM(I,K+1,N1)*DZ(I,K))/  &
                    (DZ(I,K)+DZ(I,K+1)))
                ELSE IF(K == KBM1)THEN
                    SMAX = MAX(SMAX,(WQM(I,K,N1)*DZ(I,K-1)+WQM(I,K-1,N1)*DZ(I,K))/  &
                    (DZ(I,K)+DZ(I,K-1)))
                    SMIN = MIN(SMIN,(WQM(I,K,N1)*DZ(I,K-1)+WQM(I,K-1,N1)*DZ(I,K))/(DZ(I,K)+DZ(I,K-1)))
                ELSE
                    SMAX = MAX(SMAX,(WQM(I,K,N1)*DZ(I,K-1)+WQM(I,K-1,N1)*DZ(I,K))/  &
                    (DZ(I,K)+DZ(I,K-1)),                             &
                    (WQM(I,K,N1)*DZ(I,K+1)+WQM(I,K+1,N1)*DZ(I,K))/           &
                    (DZ(I,K)+DZ(I,K+1)))
                    SMIN = MIN(SMIN,(WQM(I,K,N1)*DZ(I,K-1)+WQM(I,K-1,N1)*DZ(I,K))/  &
                    (DZ(I,K)+DZ(I,K-1)),                             &
                    (WQM(I,K,N1)*DZ(I,K+1)+WQM(I,K+1,N1)*DZ(I,K))/           &
                    (DZ(I,K)+DZ(I,K+1)))
                END IF

                IF(SMIN-WQM(I,K,N1)  > 0.0_SP)WQM(I,K,N1) = SMIN
                IF(WQM(I,K,N1) -SMAX > 0.0_SP)WQM(I,K,N1)  = SMAX

            END DO
200         CONTINUE
        END DO
    ENDDO

    RETURN
    !#  endif
    END SUBROUTINE    FCT_WASP
