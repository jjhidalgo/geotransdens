      SUBROUTINE MODELSELEC
     ;(DETHESS,LIKS1,LIKS2,LIKS3,MAINF,NPAR,NUMTOBS, STAT)

********************************************************************************
*
* PURPOSE
*
* This subroutine calculate the value of 4 different model selection criteria: 
*   AIC, BIC, ARMA and Kashyap´s criterion. 
*
* THEORY
*   model selection criteria offer a way to compare the model fit of 
*   different models of the same system with different parametrizations.
*   When the qualtity of model fit is equal these criterions tend to select the 
*   model with the smallest number of parameters. 
*
* EQUATIONS
*
*   AIC = S + 2M 
*   BIC = S + M*ln(N)
*   ARMA = S+c*M*ln(ln(N))
*   d = S + M*ln(N/(2*pi)) + ln(det(F))
*
* where S (LIKS1,LIK2,LIKS3 according to different criteria, varying 
* alpha-param.) is the log likelihood function, F is the Fisher information 
* matrix, M is the number of optimized parameters, N is the number of 
* measurements and c is a constant factor of 2 or more (2 in this code)
* 
* HISTORY    LJS (First coding)  Dec 2002
*            AAR (Revision and formatting) Jan 2003
*
********************************************************************************

      IMPLICIT NONE
                                                   ! External variables: scalars
      INTEGER*4 NPAR, NUMTOBS,MAINF
      REAL*8 LIKS1, LIKS2, LIKS3, DETHESS
                                                   ! External variables: arrays
      REAL*8 STAT(40,11)
                                                   ! Internal variables: scalars
      REAL*8 AIC1, AIC2, AIC3, BIC1,BIC2,BIC3,ARMA1, ARMA2, ARMA3
     ;      ,D1, D2, D3, PI, TAUW1
      INTEGER*4 I
C____________________________________________________ AIC criteria (Akaike 1974)

      AIC1 = LIKS1 + 2*NPAR
      AIC2 = LIKS2 + 2*NPAR
      AIC3 = LIKS3 + 2*NPAR

C____________________________________________________ BIC criteria (Akaike 1977)

      IF (NUMTOBS .GT. 0) THEN
         BIC1 = LIKS1 + npar*LOG(REAL(NUMTOBS))
         BIC2 = LIKS2 + npar*LOG(REAL(NUMTOBS))
         BIC3 = LIKS3 + npar*LOG(REAL(NUMTOBS))
      ENDIF

C__________________________________________________________ Autoregressive model

      IF (NUMTOBS .GT. 1) THEN
          ARMA1= LIKS1 + 2*NPAR*LOG(LOG(REAL(NUMTOBS)))
          ARMA2= LIKS2 + 2*NPAR*LOG(LOG(REAL(NUMTOBS)))
          ARMA3= LIKS3 + 2*NPAR*LOG(LOG(REAL(NUMTOBS)))
      ENDIF

      PI=4D0* DATAN(1.D0)
      I=1
      DO WHILE (STAT(I,1) .EQ. 0)
          I=I+1
      ENDDO
      TAUW1=STAT(I,1)/STAT(I,11)

      IF (NUMTOBS .GT. 0 ) THEN
         D1=LIKS1 + NPAR*LOG(REAL(NUMTOBS)/(2D0*PI))
     ; + DETHESS-LOG(TAUW1)
         D2=LIKS2 + NPAR*LOG(REAL(NUMTOBS)/(2D0*PI)) 
     ; + DETHESS-LOG(TAUW1)
         D3=LIKS3 + NPAR*LOG(REAL(NUMTOBS)/(2D0*PI))
     ; + DETHESS-LOG(TAUW1)
      ENDIF

C_______________________________________________________________ Echoes results

      WRITE(MAINF, 1000)
 1000 FORMAT(//,28X,'MODEL SELECTION CRITERIA',/
     ;      ,28X,'===== ========= ========',//)

                                                                 ! Formulation 1
      WRITE(MAINF, 2000)
 2000 FORMAT(10X,' TAU PARAMETERS SUCH THAT LOG-LIKELIHOOD IS '
     ;           'MINIMIZED. FORMULATION 1',/
     ;      ,10X,' --- ---------- ---- ---- -------------- -- '
     ;           '---------- -------------',//)
      
      WRITE(MAINF,1200) AIC1
      IF (NUMTOBS .GT. 0) WRITE(MAINF,1300) BIC1
      IF (NUMTOBS .GT. 1) WRITE(MAINF,1400) ARMA1
      IF (NUMTOBS .GT. 0 ) WRITE(MAINF,1500) D1
 
                                                                 ! Formulation 2
      WRITE(MAINF, 3000)
 3000 FORMAT(/,10X,' TAU PARAMETERS CONSISTENT WITH WEIGHTING '
     ;           'PARAMETERS. FORMULATION 2',/ 
     ;      ,10X,  ' --- ---------- ---------- ---- --------- '
     ;           '----------- -------------',//)
      
      WRITE(MAINF,1200) AIC2
      IF (NUMTOBS .GT. 0) WRITE(MAINF,1300) BIC2
      IF (NUMTOBS .GT. 1) WRITE(MAINF,1400) ARMA2
      IF (NUMTOBS .GT. 0 ) WRITE(MAINF,1500) D2

                                                                 ! Formulation 3
      WRITE(MAINF, 4000)
 4000 FORMAT(/,10X,' TAU PARAMETERS SET TO 1 ',/
     ;      ,10X,  ' --- ---------- --- -- - ',//)
      
       WRITE(MAINF,1200) AIC3
      IF (NUMTOBS .GT. 0) WRITE(MAINF,1300) BIC3
      IF (NUMTOBS .GT. 1) WRITE(MAINF,1400) ARMA3
      IF (NUMTOBS .GT. 0) WRITE(MAINF,1500) D3

 1200 FORMAT(5X,' AKAIKE (1974)        = ',E10.3)
 1300 FORMAT(5X,' AKAIKE (1974) ET AL. = ',E10.3)
 1400 FORMAT(5X,' HANNAN (1980)        = ',E10.3)
 1500 FORMAT(5X,' KASHAP (1982)        = ',E10.3)

      RETURN 
      END
