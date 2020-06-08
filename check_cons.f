      SUBROUTINE CHECK_CONSIST
     ;(MAINF,NTYPAR,NZPAR,INORPAR2,IVPAR,STAT,TYPENAME)

********************************************************************************
*
* PURPOSE
*
*   Check wether all estimated parameters have estimation weights. 
*
* EXTERNAL VARIABLES: ARRAYS
*
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  STAT                   Contains statistical parameters and properties
*                         for all parameter types and for all observation types. 
*  TYPENAME               Array containing  the names of the state var. and 
*                         parameter types in the same order as OBSCLASS and STAT
*
* EXTERNAL VARIABLES: SCALARS
*
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.        
*
* INTERNAL VARIABLES: SCALARS  
*
*   ESTIMPAR              logical, =true when a parameter is estimated
*   I                     dummy counter
*   J                     dummy counter
*   NOF                   index of starting location of data of device in vobs 
*   NOL                   index of last entry of data of device in vobs


* HISTORY: LJS (Dec. 2002): First coding                                  
*          AAR (Jan. 2003): Revision and formatting
*
********************************************************************************

 
      IMPLICIT NONE

      INTEGER*4 NTYPAR, NZPAR, MAINF
      INTEGER*4 INORPAR2(NTYPAR+1), IVPAR(NZPAR) 
      REAL*8 STAT(40,11)
      CHARACTER*10 TYPENAME(40)
      
      INTEGER*4 I,J,NOF,NOL
      LOGICAL ESTIMPAR

C____________________________________________ Check consistency lambda and ivpar

      DO I =1, NTYPAR

         NOF = INORPAR2(I)+1
         NOL = INORPAR2(I+1)

         ESTIMPAR = .FALSE.

         DO J = NOF,NOL
            IF (IVPAR(J).GT. 0) ESTIMPAR = .TRUE.
         ENDDO

         IF (ESTIMPAR .AND. STAT(I+10,2) .EQ. 0.D0)
     ;      WRITE(MAINF,100) TYPENAME(I+10)

      ENDDO

 100  FORMAT(//,5X,' WARNING: ',A10, ' IS ESTIMATED AND WEIGHTING '
     ;             'PARAMETER IS ZERO',/)

      RETURN
      END
