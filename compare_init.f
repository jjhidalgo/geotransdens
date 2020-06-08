       SUBROUTINE COMPARE_INIT
     ;(INDFLTR  ,INDINIT  ,INTI     ,IOINV    ,IOPINIT  ,NCONVI
     ;,NINT     ,NUMITER  ,NUMNP    ,TICAL    ,TICALAN  ,TINC
     ;,TINCINI  ,TINCLAST ,BVAR     ,IBVAR    ,TIME     ,VAUX1
     ;,VCALIT     ,VCALAN   ,SOLUTION   ,VPREV1   ,VPREV2)  

********************************************************************************
*
* PURPOSE To choose a criterium for defining the way to initialize the state 
*         variable (pressure, head or conc.) at the current time step of the 
*         nonlinear (flow or transport) solution process.
*
* DESCRIPTION This subroutine computes 
*                 a) the sum (for all nodes of the grid),of the squares of the 
*             differences between state variable at the previous solution time 
*             and those values used for their initialization when: 
*                       I- extrapolation is used , or 
*                      II- information obtained from the previous
*                          inverse problem iteration is used . 
*             
*                b) using the L-2 norm, as the maximum error found.
*
*             Then, the way adopted for initializing at the current time step 
*             is that which obtains the lower error for the previous one, 
*             according to I and II.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  BVAR                                                                         
*  IBVAR                                                                        
*  TIME                   Observation times.                                    
*  VAUX1                                                                        
*  VCAL                                                                         
*  VCALAN                                                                       
*  VCALIT                                                                       
*  VPREV1                                                                       
*  VPREV2                                                                       
*
* EXTERNAL VARIABLES: SCALARS
*
*  INDINIT                                                                      
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  IOPINIT                Scalar. It can be IOPINITH or IOPINITC.
*                         Each of these values can be: 1) Extrap. prescribed=0
*                                                      2) Use inv. prob.info. =1
*                                                      3) Decide by yourself
*                                                         using MSE =2
*                                                      4) Decide by yourself
*                                                         using L-2 norm
*  NCONVI                                                                       
*  NINT                   Number of observation times                           
*  NUMITER                Current iteration in inverse problem process          
*  NUMNP                  Number of nodes                                       
*  TICAL                                                                        
*  TICALAN                                                                      
*  TINC                   Current time increment                                
*  TINCINI                Time increment used two time steps ago                
*  TINCLAST                                                                     
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  AVERAGES_VAR
*
* HISTORY: First coding: German Galarza (Nov-1997)
*          Revision: Andres Alcolea (Oct-1997)
*
********************************************************************************


       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       DIMENSION
     ;      IBVAR(NUMNP)      ,BVAR(NUMNP)        ,VAUX1(NUMNP)
     ;     ,VCALIT(NUMNP)     ,SOLUTION(NUMNP)    ,VCALAN(NUMNP)
     ;     ,VPREV1(NUMNP)     ,VPREV2(NUMNP)      ,TIME(NINT)

C______________________________Checks the way prescribed by the user for 
C______________________________initializing state variable.


C______________________________Extrapolation is prescribed or no inv. problem 
C______________________________information is available (only simulating or 
C______________________________code is performing first Marquardt iteration

       IF(IOPINIT.EQ.0.OR.IOINV.LE.0.OR.NUMITER.EQ.1) THEN
         INDINIT=0
         RETURN

C______________________________The use of inv. prob.information is prescribed

       ELSE IF (IOPINIT.EQ.1) THEN    
        INDINIT=1                      
        RETURN
       ENDIF

C______________________________If user has not prescribed any condition, code
C______________________________decides by itself.

C______________________________Initializes comparison criteria

       SUMEINTER=0D0
       SUMEFILE=0D0
       XMAXINTER=0D0
       XMAXFILE=0D0

C ------------------------------------Computes useful time values

       TIZERO=TIME(1)
       TKMS1=TIME(INTI)+TICAL
       TK=TIME(INTI)+TICALAN

C______________________________ Compute head/conc. values at the exact Kth time,
C______________________________ obtained in the previous inverse iteration

       CALL AVERAGES_VAR 
     ;(INDFLTR   ,NUMNP    ,TIZERO   ,TK       ,TKMS1    ,BVAR     
     ;,VAUX1     ,SOLUTION)                                            



C______________________________ Computes the quadratic mean error between the
C______________________________ extrap. solution and the current one or between
C______________________________ the initialized using inv. prob. info. and the
C______________________________ current one.

C______________________________ Process has reached convergence once.

       IF (NCONVI.EQ.1)THEN 

         DO I=1,NUMNP
           IF(IBVAR(I).EQ.1) GOTO 100
           VAR=VCALAN(I)

           IF(IOPINIT.EQ.2) THEN

             SUMEINTER=SUMEINTER+(VCALIT(I)-VAR)**2
             SUMEFILE=SUMEFILE+
     ;            ((SOLUTION(I)-VAUX1(I))-(VCALIT(I)-VAR))**2

           ELSE

             XABSINTER=DABS(VCALIT(I)-VAR)
             XABSFILE=DABS(SOLUTION(I)-VAUX1(I)-VCALIT(I)+VAR)

             IF(XABSINTER.GT.XMAXINTER) XMAXINTER=XABSINTER
             IF(XABSFILE.GT.XMAXFILE) XMAXFILE=XABSFILE

           END IF

100        CONTINUE
         END DO

C______________________________ Process has reached convergence twice.

       ELSE IF(NCONVI.EQ.2)THEN             !third time step

         DO I=1,NUMNP
           IF(IBVAR(I).EQ.1)GOTO 110
           VAR=VCALAN(I)+(VCALAN(I)-VPREV1(I))/TINCLAST*TINC  

           IF(IOPINIT.EQ.2) THEN

             SUMEINTER=SUMEINTER+(VCALIT(I)-VAR)**2
             SUMEFILE=SUMEFILE+
     ;            ((SOLUTION(I)-VAUX1(I))-(VCALIT(I)-VCALAN(I)))**2

           ELSE

             XABSINTER=DABS(VCALIT(I)-VAR)
             XABSFILE=DABS(SOLUTION(I)-VAUX1(I)-VCALIT(I)+VCALAN(I))

             IF(XABSINTER.GT.XMAXINTER) XMAXINTER=XABSINTER
             IF(XABSFILE.GT.XMAXFILE) XMAXFILE=XABSFILE

           END IF

110        CONTINUE
         END DO


C______________________________ Process has reached convergence three times or 
C______________________________ more.

       ELSE                     

         DO I=1,NUMNP
           IF(IBVAR(I).EQ.1)GOTO 120
           X2=TINCINI
           X3=X2+TINCLAST
           X4=X3+TINC
           Y1=VPREV2(I)
           Y2=VPREV1(I)
           Y3=VCALAN(I)
           DIVI=(X2*(Y3-Y1)-X3*(Y2-Y1))

           IF (DIVI.NE.0D0)THEN
             X0=(X2*X2*(Y3-Y1)-X3*X3*(Y2-Y1))/2D0/DIVI
             DIVI1=(2D0*X2*X0-X2*X2)
             IF (DIVI1.NE.0D0)THEN
               A=(Y1-Y2)/DIVI1
               Y0=-A*X0*X0+Y1
               VAR=Y0+A*(X4-X0)**2             
             ELSE
               VAR=(Y3-Y2)/TINCLAST*TINC+Y3    
             ENDIF
           ELSE
             VAR=Y3                            
           ENDIF

           IF(IOPINIT.EQ.2) THEN

             SUMEINTER=SUMEINTER+(VCALIT(I)-VAR)**2
             SUMEFILE=SUMEFILE+
     ;            (SOLUTION(I)-VAUX1(I)-VCALIT(I)+VCALAN(I))**2

           ELSE

             XABSINTER=DABS(VCALIT(I)-VAR)
             XABSFILE=DABS(SOLUTION(I)-VAUX1(I)-VCALIT(I)+VCALAN(I))

             IF(XABSINTER.GT.XMAXINTER) XMAXINTER=XABSINTER
             IF(XABSFILE.GT.XMAXFILE) XMAXFILE=XABSFILE

           END IF

 120       CONTINUE
         END DO
       ENDIF


C______________________________ Code decides to extrapolate.

       INDINIT=0

C______________________________ Code decides to use inv. prob. info.

       IF(IOPINIT.EQ.2) THEN
   
         IF(SUMEFILE.LT.SUMEINTER) INDINIT=1

       ELSE

         IF(XMAXFILE.LT.XMAXINTER) INDINIT=1

       END IF


       RETURN
       END

