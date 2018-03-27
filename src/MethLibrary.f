       SUBROUTINE BIOVITERBII(ETAV,P,EMISSION,T,KTILDE,PATH,PSI)
       
       IMPLICIT NONE
       INTEGER T,KTILDE,IND,PATH(T),PSI(KTILDE,T),I,J,K,COUNTP
       DOUBLE PRECISION NORM,NORM1,NUMMAX
       DOUBLE PRECISION ETAV(KTILDE),EMISSION(KTILDE,T)
       DOUBLE PRECISION P(KTILDE,KTILDE),DELTA(KTILDE,T)
       DOUBLE PRECISION NDELTA(KTILDE,T),PDELTA(KTILDE,T)
       
       
       DO 202 I=1,KTILDE
          DELTA(I,1)=ETAV(I)+EMISSION(I,1)
          PSI(I,1)=0
 202   CONTINUE

       COUNTP=0
       DO 203 I=2,T
          DO 213 J=1,KTILDE
             NUMMAX=0.0
             NUMMAX=DELTA(1,I-1)+P(1,(J+COUNTP))
             IND=1
             DO 223 K=2,KTILDE
                IF((DELTA(K,I-1)+P(K,(J+COUNTP))).GT.NUMMAX) THEN
                   NUMMAX=(DELTA(K,I-1)+P(K,(J+COUNTP)))
                   IND=K
                ENDIF
 223         CONTINUE
             PSI(J,I)=IND
             DELTA(J,I)=NUMMAX+EMISSION(J,I)
 213      CONTINUE
          COUNTP=COUNTP+KTILDE
 203   CONTINUE
       
       NUMMAX=0.0
       NUMMAX=DELTA(1,T)
       IND=1
       DO 253 K=2,KTILDE
          IF(DELTA(K,T).GT.NUMMAX) THEN
             NUMMAX=DELTA(K,T)
             IND=K
           ENDIF
 253   CONTINUE
       
       PATH(T)=IND
       DO 263 K=T-1,1,-1
          PATH(K)=PSI(PATH(K+1),K+1)
 263   CONTINUE

       RETURN
       END 
	
       SUBROUTINE TRANSEMISI(MUK,NCOV,TOTALSEQ,KS,
     c COV,SEPSILON,T,PT,P,EMISSION)

       IMPLICIT NONE

       INTEGER T,KS,I,J,K,COUNTP,NCOV
       DOUBLE PRECISION PI,ELNSUB,COV(NCOV)
       DOUBLE PRECISION P(KS,KS*NCOV)
       DOUBLE PRECISION EMISSION(KS,T)
       DOUBLE PRECISION MUK(KS),PT(KS)
       DOUBLE PRECISION TOTALSEQ(T)
       DOUBLE PRECISION SEPSILON(KS)
       PARAMETER(PI=3.14159265358979)       

       COUNTP=0
       DO 700 I=1,NCOV
           DO 710 J=1,KS
               DO 720 K=1,KS
                   IF (J.EQ.K) THEN
                   P(J,(K+COUNTP))=ELNSUB(0.,(PT(K)+COV(I)))
                   ELSE
                   P(J,(K+COUNTP))=PT(K)+COV(I)
                   ENDIF
 720           CONTINUE
 710       CONTINUE
           COUNTP=COUNTP+KS
 700   CONTINUE



       DO 730 J=1,KS
          DO 740 K=1,T
              EMISSION(J,K)=EMISSION(J,K)+LOG(1/
     c        (SQRT(2*PI)*SEPSILON(J)))+(-0.5*((TOTALSEQ(K)-
     c        MUK(J))/(SEPSILON(J)))**2)
 740      CONTINUE
 730   CONTINUE
       RETURN
       END 





      DOUBLE PRECISION FUNCTION ELNSUM(X,Y)
	  
      IMPLICIT NONE
      DOUBLE PRECISION X,Y
      IF (X.GT.Y) THEN
         ELNSUM=X+LOG(1+EXP(Y-X))
      ELSE
         ELNSUM=Y+LOG(1+EXP(X-Y))
      ENDIF

      RETURN
      END


      DOUBLE PRECISION FUNCTION ELNSUB(X,Y)
	  
      IMPLICIT NONE
      DOUBLE PRECISION X,Y
      IF (X.GT.Y) THEN
         ELNSUB=X+LOG(1-EXP(Y-X))
      ELSE
         ELNSUB=Y+LOG(1-EXP(X-Y))
      ENDIF

      RETURN
      END

