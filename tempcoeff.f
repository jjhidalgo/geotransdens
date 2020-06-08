       REAL*8 FUNCTION TEMPCOEFF (DTIM,IDIMFNT,INTI,NFTIME,NINT,FNT)

*****************************************************************************
*
* PURPOSE
*      Computes the temporal coefficient of a given parameter
*
* DESCRIPTION
*
*      Makes a linear interpolation between two consecutive time function 
*      values as 
*                 (1-DTIM)*FNT(NFTIME,INTI)+DTIM*FNT(NFTIME,INTI+1)
*
*      Where 0 <= DTIM <= 1 accounts for the relative location of the current 
*      time inside the interval given by TIME(INTI) and TIME(INTI+1), ie, 
*      between two consecutive observation times. 
*                  
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FNT                    Array containing time functions values                
*
* EXTERNAL VARIABLES: SCALARS
*
*  DTIM                   Current time for the computation of current time
*                         function value (counted since the beginning of the
*                         current observation interval, not since the beginning
*                         of the problem) divided by the length of the
*                         observation interval (interval between two
*                         consecutive observation times). Used only to make a
*                         linear interpolation of time function values at 
*                         computation times.
*  IDIMFNT                Used for dimensioning array FNT, it coincides with    
*                         NFNT if the latter is not zero                       
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  NFTIME                 Time function number of the current parameter zone
*  NINT                   Number of observation times                           
*
* HISTORY
*
*     AMS        1988     First coding
*     AMS      4-1998     Revision and common elimination
*
*****************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)

       DIMENSION FNT(IDIMFNT,NINT)

C------------------------- FIRST (and unique) EXECUTABLE STATEMENT.

       TEMPCOEFF=(1.D0-DTIM)*FNT(NFTIME,INTI) + DTIM*FNT(NFTIME,INTI+1)

       RETURN
       END
