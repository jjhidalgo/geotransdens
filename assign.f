       SUBROUTINE ASSIGN
     ;(  NPAR       ,NZPAR     ,IVPAR     ,INDPAR
     ;  ,PAR        ,PARC)
       
*****************************************************************
***   ASSIGN NEW VALUE OF PARAMETERS
*****************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION PAR(NPAR),IVPAR(NZPAR),PARC(NZPAR),INDPAR(NPAR)

       K=1
       DO I=1,NZPAR
          IF(IVPAR(I).NE.0) THEN
  
C-------------------------------Anti-log transformation
               IF(INDPAR(K).NE.0) THEN
                  PARC(I)=10**(DLOG10(PARC(I))+PAR(K))
               ELSE

C-------------------------------Assign new value of parameters
                  PARC(I)=PARC(I)+PAR(K)
               ENDIF
               K=K+1
          ENDIF
       ENDDO
       RETURN
       END 
