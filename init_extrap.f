       SUBROUTINE INIT_EXTRAP
     ;(IDIMFNT  ,INDFLTR  ,INTI       ,NCONVI     ,NINT     
     ;,NPARNP   ,NTYPAR   ,NUMNP      ,NZPAR      ,THETA    ,TICALAN  
     ;,TINC     ,TINCINI  ,TINCLAST   ,TINTERVOBS ,CFPARNP  ,FNT      
     ;,IBVAR    ,INORPAR  ,IXPARNP    ,NFTPAR     ,PARC     ,VCALAN   
     ;,VCALIT   ,VPREV1   ,VPREV2)  

********************************************************************************
*
* PURPOSE To assign initial values to the state variable, at the begining of the
*         Newton-Raphson iterative process, related to the nonlinear flow or 
*         transport equation solution, by extrapolation of previous solutions
*
*
* DESCRIPTION The algorithm extrapolates the variable (head, pressure or solute
*             concentracion at a particular time, from solutions obtained for 
*             previous times
*
* EXTERNAL VARIABLES: ARRAYS
*
*  CFPARNP                Array containing node coefinient of node j            
*                         corresponding to INpar index parameter zone.          
*  FNT                    Array containing time functions values                
*  IBVAR                  Array containing IBCOD or IBTCO
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IXPARNP                Array containing zone number of each node j,          
*                         corresponding to INpar index parameter zone.          
*  NFTPAR                 Vector containing time function number at every       
*                         parameter zone                                        
*  PARC                   Vector containing calculated values for all           
*                         parameters                                            
*  VCALAN                                                                       
*  VCALIT                                                                       
*  VPREV1                                                                       
*  VPREV2                                                                       
*
* EXTERNAL VARIABLES: SCALARS
*
*  DTIMBVAR               Scalar. It is equal to DTIMBFLU or DTIMBTRA
*  IDIMFNT                First dimension of array FNT, it coincides with       
*                         NFNT if the latter is not zero                        
*  INDFLTR                Indicator variable (0=flow, 1=transport)
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  NCONVI                 Number of convergences reached.
*  NINT                   Number of observation times                           
*  NPARNP                 Number of nodal parameters in current problem         
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMNP                  Number of nodes                                       
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  THETA                  Temporal weighting parameter for flow or transport 
*                         state variable
*  TICALAN                                                                      
*  TINC                   Current time increment                                
*  TINCINI                Time increment used two time steps ago                
*  TINCLAST                                                                     
*  TINTERVOBS             Length of the current observation time interval       
*
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  EQUAL_ARRAY                                                                  
*  TEMPCOEFF                                                                    
*
* HISTORY: First coding: German Galarza (Dec-1997)
*          Revision:     Andres Alcolea (Oct-1998)
*
********************************************************************************



       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       DIMENSION
     ;    IBVAR(NUMNP)   ,IXPARNP(NUMNP,NPARNP)  ,CFPARNP(NUMNP,NPARNP)
     ;   ,PARC(NZPAR)    ,VCALIT(NUMNP)          ,VCALAN(NUMNP)
     ;   ,VPREV1(NUMNP)  ,VPREV2(NUMNP)          ,INORPAR(NTYPAR)
     ;   ,NFTPAR(NZPAR)  


C______________________________ Has not converged yed. k=0

       IF (NCONVI.EQ.0)THEN         

         CALL EQUAL_ARRAY (VCALIT,VCALAN,NUMNP)

C______________________________ Has converged once. k=1 . Linear interpolation

       ELSE IF(NCONVI.EQ.1)THEN     

         DO I=1,NUMNP
           VCALIT(I)=VCALAN(I)+(VCALAN(I)-VPREV1(I))/TINCLAST*TINC
         END DO

C______________________________ Has converged twice or more. Linear or quadratic

       ELSE                               

         DO I=1,NUMNP

           IF(IBVAR(I).EQ.1)GOTO 100
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
               VCALIT(I)=Y0+A*(X4-X0)**2            
             ELSE
               VCALIT(I)=(Y3-Y2)/TINCLAST*TINC+Y3   
             ENDIF
           ELSE
             VCALIT(I)=Y3                           
           ENDIF

100        CONTINUE

         END DO

       ENDIF



C______________________________Only for nodes with prescribed head or conc.
C______________________________boundary condition 

       DO I=1,NUMNP

         IF (IBVAR(I).EQ.1)THEN

C______________________________Flow equation

           IF(INDFLTR.EQ.0) THEN    

             NZVAR=IXPARNP(I,2)
             VAR=PARC(INORPAR(9)+NZVAR)*CFPARNP(I,2)
             NFVAR=NFTPAR(INORPAR(9)+NZVAR)

C______________________________Transport equation

           ELSE              


             NZVAR=IXPARNP(I,8)
             VAR=PARC(INORPAR(18)+NZVAR)*CFPARNP(I,8)
             NFVAR=NFTPAR(INORPAR(18)+NZVAR)

           END IF


           DTIMBVAR=(TICALAN+THETA*TINC)/TINTERVOBS
           IF(NFVAR.NE.0) VAR=VAR
     ;           *TEMPCOEFF(DTIMBVAR,IDIMFNT,INTI,NFVAR,NINT,FNT)



           VCALIT(I)=VAR

          

         ENDIF

       END DO

       RETURN
       END
