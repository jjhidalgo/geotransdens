      SUBROUTINE WRI_PARAM_HISTORY
     ;(IOWIT   ,ISOT   ,MAINF     ,NTYPAR   ,NZPAR
     ;,INORPAR ,IVPAR  ,NZONE_PAR ,PARGOOD)

********************************************************************************
*     
* PURPOSE Writes values of the estimated/interpolated parameters
*
* DESCRIPTION 
*
* - Step 0: Declaration of variables
* - Step 1: Writes history of ESTIMATED parameters.
*   - Step 1.1: Writes main header and reads param. set (IOWIT<>0) or last 
*             PARGOOD (IOWIT=0) is written
*   - Step 1.2: Writes values. Transmissivity, due to anisotropy, requires an 
*               special treatment
*
* EXTERNAL VARIABLES: ARRAYS
*
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IVPAR                  Vector containing estimation index for all            
*                         parameters  
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters                                            
*  PARAUX                 Here, will be read. External due to workspace
*
* EXTERNAL VARIABLES: SCALARS

*  IOWIT                  If 1, estimated parameters history through the
*                         optimization problem are written. If 0, only
*                         last estimated parameters are written
*  ISOT                   Anisotropy index
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*
* INTERNAL VARIABLES: ARRAYS
*
* INTERNAL VARIABLES: SCALARS
* 
*  ITER                   Dummy counter of iterations
*  ITYPE                  Dummy counter of types of parameters
*  NUMITER                Iteration number, as read
*  VAR_TYPE               Parameter type name                                
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  WRI_PAR_HIST_AUX       Writes the values of the estimated parameters of 
*                         a given type
*
* HISTORY
*
*     AAR        9-2003   First coding
*
********************************************************************************

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 IOWIT,MAINF,NTYPAR,NZPAR,ISOT
     ;         ,NZONE_PAR(NTYPAR),IVPAR(NZPAR),INORPAR(NTYPAR)
                                                                 ! Real external
      REAL*8 PARGOOD(NZPAR)
                                                              ! Integer internal
      INTEGER*4 ITER,ITYPE,NUMITER,ISUMFO
                                                                    ! Characters
      CHARACTER VAR_TYPE(14)*16
      DATA VAR_TYPE /'TRANSMISSIVITY  ','STORAGE COEFF.  '
     ;              ,'AREAL RECHARGE  ','PRESCRIBED HEAD '
     ;              ,'PRESCRIBED FLOW ','LEAKAGE COEFF.  '
     ;              ,'MOLECULAR DIFF. ','LONGIT. DISPERS.'
     ;              ,'TRANSV. DISPERS.','POROSITY        '
     ;              ,'LINEAR DECAY    ','RETARDATION     '
     ;              ,'EXTERNAL CONC.  ','GENERIC PARAM.  '/


C_______________________ Step 1: Writes history of parameters.

      REWIND (70)
      DO ITER=1,11111

C_______________________ Step 1.1: Writes main header and reads param. 
C_______________________           set (IOWIT<>0) or writes last PARGOOD 
C_______________________           (IOWIT=0)

        IF (IOWIT.NE.0) THEN
          READ (70,END=999) ISUMFO,NUMITER,PARGOOD
          IF (ISUMFO.EQ.0) THEN
              WRITE (MAINF,2000) NUMITER
          ELSE
              WRITE (MAINF,2001) NUMITER
          END IF
        ELSE
            IF (ITER.GT.1) RETURN     ! Only once, with last estimated param.
            WRITE (MAINF,1000)
        END IF

 1000   FORMAT (/,5X,' LAST ESTIMATED PARAMETERS',/
     ;           ,5X,' ==== ========= ==========',/)
 2000   FORMAT(//,' GOOD ITERATION:',I5,/,' ==== ==========',//,12X,
     ;             'TYPE    ZONE  ISOZ         VALUE',/)
 2001   FORMAT(//,' BAD ITERATION:',I5,/,' === ==========',//,12X,
     ;             'TYPE    ZONE  ISOZ         VALUE',/)
        
C_______________________ Step 1.2: Writes values. Transmissivity, due to 
C_______________________           anisotropy requires an special treatment

        IF (NZONE_PAR(1).NE.0) CALL WRI_PAR_HIST_AUX
     ;(1           ,ISOT      ,MAINF   ,NZONE_PAR(1)
     ;,VAR_TYPE(1) ,IVPAR(1)  ,PARGOOD(1))
        
        DO ITYPE=2,14

           IF (NZONE_PAR(ITYPE).NE.0) CALL WRI_PAR_HIST_AUX
     ;(0                ,1                  ,MAINF  
     ;,NZONE_PAR(ITYPE) ,VAR_TYPE(ITYPE)    ,IVPAR(INORPAR(ITYPE+5)+1)
     ;,PARGOOD(INORPAR(ITYPE+5)+1))

        END DO ! ITYPE=2,14

      END DO ! ITER=1,11111

 999  RETURN
      END

********************************************************************************
********************************************************************************
********************************************************************************

      SUBROUTINE WRI_PAR_HIST_AUX
     ;(IOP_TRA  ,ISOT  ,MAINF  ,NZVAR  ,VARTYPE  ,IVPAR ,PARGOOD)

********************************************************************************
*
* PURPOSE Writes values of the estimated/interpolated parameter of a given type
*
* DESCRIPTION 
*
* - Step 0: Declaration of variables
* - Step 1: Loop over zones of parameters
*   - Step 1.1.1: Defines anisotropy character
*   - Step 1.1.2: Writes values
*
* EXTERNAL VARIABLES: ARRAYS
*
*  IVPAR                  Vector containing estimation index for all            
*                         parameters  
*  PARGOOD                Array containing estimated ZONAL parameters
*
* EXTERNAL VARIABLES: SCALARS
*
*  CHARISOZ               Value of anisotropy index (writing purposes)
*  IOP_TRA                If 0, type is not transmissivity
*  ISOT                   Anisotropy index
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NZVAR                  Total number of zones of a given param. type
*  VARTYPE                Parameter type name                                
*
* INTERNAL VARIABLES: ARRAYS
*
* INTERNAL VARIABLES: SCALARS
*
*  ICOMPO                 Dummy counter of anisotropy components
*  IZPAR                  Dummy counter of zonal parameters
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
* HISTORY
*
*     AAR        9-2003   First coding
*
********************************************************************************

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 NZVAR,MAINF,IOP_TRA,ISOT
     ;         ,IVPAR(NZVAR,ISOT)
                                                                 ! Real external
      REAL*8 PARGOOD(NZVAR*ISOT)
                                                              ! Integer internal
      INTEGER*4 IZVAR,ICOMPO
                                                                    ! Characters
      CHARACTER VARTYPE*16,CHARISOZ*1

C_______________________ Step 1: Loop over anisotropy components (only sense
C_______________________         for transmissivity. 

      DO ICOMPO=1,ISOT


C_______________________ Step 1.1: Loop over zones of parameters

        DO IZVAR=1,NZVAR


          IF (IVPAR(IZVAR,ICOMPO).NE.0) THEN

C_______________________ Step 1.1.1: Defines anisotropy character

             IF (IOP_TRA.NE.0) THEN
               WRITE(CHARISOZ,'(I1)') ICOMPO
             ELSE
               CHARISOZ='-'
             END IF
 
C_______________________ Step 1.1.2: Writes values

             WRITE(MAINF,1000) 
     ;         VARTYPE,IZVAR,CHARISOZ,PARGOOD((ICOMPO-1)*NZVAR+IZVAR)
 1000        FORMAT(A16,3X,I5,5X,A1,2X,G16.8)

          END IF ! IVPAR(IZVAR,ICOMPO).NE.0
        END DO ! IZVAR=1,NZVAR
      END DO ! ICOMPO=1,ISOT

      RETURN
      END
