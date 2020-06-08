       SUBROUTINE INIT_INV_PROB_INFO
     ;(THETA     ,IDIMFNT  ,INDFLTR  ,INTI     ,NINT     ,NPARALG
     ;,NPARNP    ,NUMNP    ,NZPAR    ,TICAL    ,TICALAN  ,TINC
     ;,TINTERVOBS,BVAR     ,CFPARNP  ,FNT      ,IBVAR    ,INORPAR
     ;,IXPARNP   ,NFTPAR   ,PARC     ,TIME     ,VAUX1    ,VCALAN
     ;,VCALIT)  

********************************************************************************
*
* PURPOSE  To assign initial values to the state variable, at the begining of 
*          the Newton-Raphson iterative process, associated to the nonlinear 
*          flow or transport equation solution, based on the previous inverse 
*          problem iteration results.
*
*
* DESCRIPTION The algoritm extracts the solution obtained at a particular 
*             computation time, during the previous inverse problem iteration. 
*             This one is stored in the vector which contains the heads or 
*             concentrations in the Newton-Raphson process.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  BVAR                   Array containing BFLU or BTRA
*  CFPARNP                Array containing node coefinient of node j            
*                         corresponding to INpar index parameter zone.          
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
*  PAR_DIR                Array containing all real direct problem              
*                         parameters                                            
*  TIME                   Observation times.                                    
*  VAUX1                                                                        
*  VCALAN                                                                       
*  VCALIT                                                                       
*
* EXTERNAL VARIABLES: SCALARS
*
*  FNT                    Array containing time functions values                
*  IDIMFNT                First dimension of array FNT, it coincides with       
*                         NFNT if the latter is not zero                        
*  INDFLTR                Indicator variable. =0-->flow ; =0-->transport
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  NINT                   Number of observation times                           
*  NPARNP                 Number of nodal parameters in current problem         
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMNP                  Number of nodes                                       
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  TICAL                                                                        
*  TICALAN                                                                      
*  TINC                   Current time increment                                
*  TINTERVOBS             Length of the current observation time interval       
*
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  AVERAGES_HEAD                                                                
*  TEMPCOEFF                                                                    
*
* HISTORY: First coding: German Galarza (Dec-1997)
*          Revision: Andres Alcolea (Oct-1998)
*
********************************************************************************

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       DIMENSION
     ;       IBVAR(NUMNP)    ,IXPARNP(NUMNP,NPARNP)    ,NFTPAR(NZPAR)
     ;      ,BVAR(NUMNP)     ,CFPARNP(NUMNP,NPARNP)    ,VAUX1(NUMNP)
     ;      ,VCALIT(NUMNP)   ,PARC(NZPAR)              ,TIME(NINT)
     ;      ,VCALAN(NUMNP)   ,INORPAR(NPARALG)

C ---------------------------------------- compute the magnitude of the time to be solved

       TIZERO=TIME(1)
       TKMS1=TIME(INTI)+TICAL
       TK=TIME(INTI)+TICALAN

C ---------------------------------------- recover the solution obtained during the previous
C ---------------------------------------- inverse problem iteration

       CALL AVERAGES_VAR
     ;(INDFLTR   ,NUMNP    ,TIZERO   ,TK       ,TKMS1    ,BVAR
     ;,VAUX1     ,VCALIT)

C ---------------------------------------- Carry out the initialization

       DO I=1,NUMNP
         VCALIT(I)=VCALAN(I)+(VCALIT(I)-VAUX1(I))
       END DO

C ---------------------------------------- Corriges nodes with prescribed head or concentration
C ---------------------------------------- Boundary condition 

       DO I=1,NUMNP

         IF (IBVAR(I).EQ.1)THEN

           IF(INDFLTR.EQ.0) THEN        !   Flow equation

             NZVAR=IXPARNP(I,2)
             VAR=PARC(INORPAR(9)+NZVAR)*CFPARNP(I,2)
             NFVAR=NFTPAR(INORPAR(9)+NZVAR)

           ELSE                    ! Transport equation

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
