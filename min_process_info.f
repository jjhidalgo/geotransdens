       SUBROUTINE MIN_PROCESS_INFO
     ;(FNEW     ,FOLD      ,IN1      ,IOWPI    ,ITER1    ,MAINF     
     ;,MESS     ,NITERF1   ,NITERF2  ,NUMITER  ,XMAXIM)

********************************************************************************
*
* PURPOSE Displays the results of the optimization process
*
* DESCRIPTION 
*
* - Step 0: Declaration of variables
* - Step 1: Writes main header
* - Step 2: Writes reason of stopping
* - Step 3: Writes maximum increment of parameters if process converged due to a
*           little increment in the parameters
* - Step 4: Writes last variation of the obj. fun. if process converged due to a 
*           little increment in the objective function
* - Step 5: Writes Marquardt's process history
* - Step 6: Writes header of last better parameters           
*
* EXTERNAL VARIABLES: ARRAYS
*
* INTERNAL VARIABLES: ARRAYS
*
* EXTERNAL VARIABLES: SCALARS
*
*  FNEW                   Objective function value in the current iteration     
*  FOLD                   Objective function computed value in last iteration   
*  IN1                    Flag for stopping
*  IOWPI                  Allows writing the results of the minim. process
*  ITER1                  Actual iteration
*  MAINF                  Main output file unit
*  MESS                   Message of stopping
*  NITERF1                Number of failed iterations (total)
*  NITERF2                Number of failed iterations (since last good one)
*  NUMITER                Actual iteration
*  XMAXIM                 Maximum increment of unknown parameters
*
* INTERNAL VARIABLES: SCALARS
*
*  GNORM                  Gradient norm (as read)
*  ISUMFO                 Flag (as read)
*  IT                     Dummy counter
*  NUMITER_AUX            Iteration number (as read)
*  OBJCON                 Objective function of concentrations (as read)
*  OBJGEO                 Objective function of parameters (as read)
*  OBJHED                 Objective function of head levels (as read)
*  OBJPAR                 Objective function of geophysical information (as read)
*  PHI                    Coeff. of goodness of quad. approach (as read)
*  XMARQ                  Marquardt's parameter (as read)
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
* HISTORY: AMS,GGL,JCR: First coding (???)
*          AAR: Header inclusion (Aug-2003)
*          JHG: Added geophysical observations (May-2020)
*     
********************************************************************************

C_______________________ Step 0: Declaration of variables

       IMPLICIT NONE
                                                              ! Integer external
       INTEGER*4 IN1,MAINF,ITER1,NITERF1,NITERF2,NUMITER,IOWPI
                                                                 ! Real external
       REAL*8 XMAXIM,FNEW,FOLD
                                                              ! Integer internal
       INTEGER*4 IT,ISUMFO,NUMITER_AUX
                                                                 ! Real internal
       REAL*8 XMARQ,OBJHED,OBJCON,OBJPAR,OBJGEO,GNORM,PHI
                                                                     ! Characters
       CHARACTER*50 MESS
        
C_______________________ Step 1: Writes main header

       IF (IN1.NE.99) THEN

         WRITE (MAINF,100)
 100     FORMAT
     ;   (//,9X,'  M I N I M I Z A T I O N   I N F O R M A T I O N')   
           
C_______________________ Step 2: Writes reason of stopping

        WRITE (MAINF,200) MESS,ITER1,NITERF1,NITERF2,NUMITER
 200     FORMAT(/' THE PROGRAM WAS STOP BECAUSE:',/,A50,/,
     ;   ' INDICATORS CORRESPONDING TO NEXT PARAMETERS',
     ;   ' (LAST GOOD ITERATION)',/,
     ;   ' WILL BE FIND AT INFORMATION OF ',I3,'th ITERATION',//,
     ;   ' FALL ITERATIONS SINCE LAST GOOD ITERATION..............',I5,/,
     ;   ' FALL TOTAL ITERATIONS..................................',I5,/,
     ;   ' ORDER OF LAST ITERATION................................',I5)
  
C_______________________ Step 3: Writes maximum increment of parameters
C_______________________         if process converged due to a little increment 
C_______________________         in the parameters

         IF (IN1.EQ.3) WRITE(MAINF,300) XMAXIM*100
 300     FORMAT(/,' MAXIM INCREMENT OF PARAMETERS........',E10.3,'%')
         
C_______________________ Step 4: Writes last variation of the obj. fun.
C_______________________         if process converged due to a little increment 
C_______________________         in the objective function

         IF (IN1.EQ.5) WRITE(MAINF,400) (FNEW-FOLD)/FOLD*100
 400     FORMAT(/,
     ;     ' LAST VARIATION OF OBJECTIVE FUNCTION........',E10.3,'%')   
    
C_______________________ Step 5: Writes Marquardt's process history

         IF (IOWPI.NE.0) THEN

           REWIND(71)
           DO IT=1,11111
             READ (71,END=999) NUMITER_AUX,ISUMFO,FNEW,OBJHED,OBJCON
     ;              ,OBJPAR,OBJGEO,XMARQ,PHI,XMAXIM,GNORM

             IF (NUMITER_AUX.EQ.1) WRITE(MAINF,500)
 500     FORMAT(///,' ITER BAD      FNEW    OBJHED    OBJCON    OBJPAR',
     ;      '    OBJGEO     XMARQ  PHI(k-1)    XMAXIM     GNORM',/,
     ;              ' ---- ---      ----    ------    ------    ------',
     ;      '     -----     -----  --------    ------     -----',//)

             WRITE (MAINF,600) NUMITER_AUX,ISUMFO,FNEW,OBJHED,OBJCON
     ;              ,OBJPAR,OBJGEO,XMARQ,PHI,XMAXIM,GNORM

 600         FORMAT(1P,I5,I4,9E10.3) 
           END DO

         END IF ! IOWPI.NE.0

C_______________________ Step 6: Writes header of last better parameters           

 999     WRITE (MAINF,700)
 700     FORMAT(///,31X,'PARAMETERS HISTORY',/,
     ;              31X,'========== =======')

       END IF ! IN1.NE.99

       RETURN
       END
