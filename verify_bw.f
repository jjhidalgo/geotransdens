       SUBROUTINE VERIFY_BW
     ;(   LNNDEL     ,KXX        ,NUMEL      ,NBAND       ,LMXNDL
     ;   ,IERROR     ,IOWAR      ,MAINF      ,IUGRID      ,FILENAME)
********************************************************************************
* PURPOSE 
*
*    This subroutine verifyes if input bandwidth is the calculated band with
* 
* DESCRIPTION
*
*    This subroutine must be called when all elements conectivities have been
*    read and interpolated. True bandwidth is calculated as the maximum 
*    difference between two node order of any element. So that, there´s a first 
*    loop that starts from 1 to NUMEL elements of the curent problem. For each
*    element, there are two loops that goes throug each element nodes, 
*    calculating IL (=difference between two node order of current element). 
*    LMAX is actualized with the maximum IL calculated in any element. 
*    Finally, LMAX value is assigned to NVBAN (calculated bandwith).
* 
* EXTERNAL VARIABLES: ARRAYS 
*
*    LNNDEL(NE)           Number of nodes for element NE                            
*    KXX(I,NE)            Node number for element node I
*                         (I is in countedclockwise order) 
*    FILENAME             Array containing all input and output filenames
*
* EXTERNAL VARIABLES: SCALARS
*
*    NUMEL                Number of elements
*    NBAND                Input value for bandwidth
*    LMXNDL               Maximum number of nodes per element
*    IERROR               Number of errors detected reading data
*    IOWAR                Allows printing data in MAIN file
*                     
* INTERNAL VARIABLES: SCALARS 
*
*    NVBAN                True bandwidth (Calculated)
*    N1                   Node number
*    N2                   Node number
*    IL                   Absolute difference between two node numbers
*                         of an element
*    LMAX                 Maximum value of IL for every pair of nodes of
*                         any element    
*
* SUBROUTINES REFERENCED
*
*    ERROR                Writes the current error message and error number 
*                         on MAIN FILE.
* HISTORY
*
*     SCR 15-abr-1997     Firt coding
*
********************************************************************************
     
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION KXX(LMXNDL,NUMEL),LNNDEL(NUMEL) 
       CHARACTER FILENAME(20)*20
 
       NVBAN=1
      
       DO NE=1,NUMEL
       
         LMAX=0
         DO I=1,LNNDEL(NE)-1
            N1=KXX(I,NE)
            DO J=I+1,LNNDEL(NE)
               N2=KXX(J,NE)
               IL=IABS(N1-N2)
               IF (IL.GT.LMAX) LMAX=IL
            END DO
         END DO
        
         IF (LMAX.GT.NBAND) THEN 
           
            CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                 'INCORRECT BANDWITH'
     ;                 ,0,0,IUGRID,1,3.01)       
                  
            WRITE(MAINF,3000) NE
 3000       FORMAT(/,2X,'BANDWIDTH EXCEEDED IN ELEMENT NUMBER ',I5)
 
         END IF 
       
         IF (LMAX.GT.NVBAN) NVBAN=LMAX

       END DO


        
       IF (NVBAN.GT.NBAND) THEN 
       
         WRITE(MAINF,3100) NVBAN,NBAND
 3100    FORMAT(/,'INPUT BANDWITH IS LOWER',
     ;            ' THAN TRUE BANDWIDTH',
     ;          //,5X,'TRUE BANDWIDTH  =',I5,/
     ;            5X,'INPUT BANDWIDTH =',I5)
     
       ELSE 
         
         IF (NVBAN.LT.NBAND) WRITE(MAINF,3200) NVBAN,NBAND
 
 
 3200    FORMAT(/' INPUT BANDWIDTH IS GREATHER',
     ;           ' THAN TRUE BANDWIDTH',
     ;          //,5X,'TRUE BANDWIDTH  =',I5,/
     ;               5X,'INPUT BANDWIDTH =',I5)
         

       END IF    
       
       RETURN
       
       END
