       SUBROUTINE AVERAGES_VAR
     ;(INDFLTR   ,NUMNP    ,TIZERO   ,TK       ,TKMS1    ,BVAR     
     ;,VAUX1     ,SOLUTION)                                            

********************************************************************************
*
* PURPOSE  Recover heads or concentration values at times K and K+1, computed
*          during the previous inverse problem iteration, for all nodes of the 
*          finite element grid.
*
* DESCRIPTION Examines the file containing the flow or transport equation
*             solution at the previous inverse iteration. Based on this 
*             information, the algorithm computes state variable at two 
*             particular times, K and K+1.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  BVAR                                                                         
*  VAUX1                                                                        
*  SOLUTION                                                                       
*
* EXTERNAL VARIABLES: SCALARS
*
*  INDFLTR
*  NUMNP                  Number of nodes                                       
*  TIZERO                                                                       
*  TK                                                                           
*  TKMS1                                                                        
*
* INTERNAL VARIABLES: SCALARS
*
*  I                                                                            
*  INDICT                                                                       
*  NORDEN                                                                       
*  NORDENK                                                                      
*  NORDENKMS1                                                                   
*  TEMPS                                                                        
*  TICOLD                                                                       
*  TICOMP                                                                       
*  XPONDK                                                                       
*  XPONDKMS1                                                                    
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ZERO_ARRAY                                                                   
*
* HISTORY: First coding : German Galarza (Dec-1997)
*          Revision : Andres Alcolea (Oct-1998)
*
********************************************************************************

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       DIMENSION BVAR(NUMNP)  ,VAUX1(NUMNP)  ,SOLUTION(NUMNP)


C______________________________ Initializes the auxiliar vector BVAR
C______________________________ (BFLU O BTRA)

       CALL ZERO_ARRAY(BVAR,NUMNP)

C______________________________Rewind the file containing the direct problem
C______________________________solution indices at the previous inverse 
C______________________________iteration

       IUNIT=52+4*INDFLTR
       REWIND(IUNIT)

C______________________________Founds current time.

C______________________________Initializes the time position indicator INDICT
C______________________________and the absolute value of the time whose solution
C______________________________is required.

       INDICT=0
       TICOMP=TK

C______________________________K corresponds to the first simulation time order

       IF (DABS(TK-TIZERO).LT.1E-10)THEN     
         NORDENK=2
         XPONDK=0D0
         INDICT=1
         TICOMP=TKMS1
         TICOLD=TK
       ENDIF

C______________________________Reads time and its corresponding order number

 100   READ(IUNIT)TEMPS,NORDEN
 
C______________________________Time read is great or equal than that looked for

       IF (TEMPS.GE.TICOMP)THEN     

C______________________________Kth time has not been found yet (so it is found now)

         IF (INDICT.EQ.0)THEN 
           NORDENK=NORDEN

C______________________________Computes the ponderation factor for Kth time

           XPONDK=1D0-((TEMPS-TICOMP)/(TEMPS-TICOLD))        
           INDICT=1                                        

C______________________________Updates the time searched
  
           TICOMP=TKMS1                                    
           IF (TKMS1.LE.TEMPS)THEN   

C______________________________(K+1)-th time is also less than that read

             NORDENKMS1=NORDENK

C______________________________Computes the ponderation factor for (K+1)th time

             XPONDKMS1=1D0-((TEMPS-TICOMP)/(TEMPS-TICOLD))  
             GOTO 200
           ENDIF
           TICOLD=TEMPS

         ELSE IF (INDICT.EQ.1) THEN 

C______________________________Kth time has found already (so it is found the 
C______________________________(k+1)th now)

           NORDENKMS1=NORDEN

C______________________________Computes the ponderation factor for (K+1)th time

           XPONDKMS1=1D0-((TEMPS-TICOMP)/(TEMPS-TICOLD)) 
           GOTO 200
         ENDIF

       ELSE

C______________________________Storage the last absolute time found

         TICOLD=TEMPS                                 

       ENDIF

       GOTO 100


C______________________________Read solutions at the previous inverse iteration
C______________________________from a direct access file

 200   READ(IUNIT-2,REC=NORDENK-1) BVAR 
       READ(IUNIT-2,REC=NORDENK) VAUX1  

C______________________________Computes the solution at kth time

       DO I=1,NUMNP
         VAUX1(I)=XPONDK*VAUX1(I)+(1D0-XPONDK)*BVAR(I)   
       END DO

C______________________________Read solutions at the previous inverse iteration
C______________________________from a direct access file :K+1 solution

       READ(IUNIT-2,REC=NORDENKMS1-1) BVAR  
       READ(IUNIT-2,REC=NORDENKMS1) SOLUTION  

C______________________________Computes the solution at (k+1)th time

       DO I=1,NUMNP
         SOLUTION(I)=XPONDKMS1*SOLUTION(I)+(1D0-XPONDKMS1)*BVAR(I)  
       END DO

       RETURN       
       END


                                                                        
