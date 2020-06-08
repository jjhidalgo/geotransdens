       SUBROUTINE LEC_CFE
     ; (IERROR   ,INARR    ,INARRT   ,INCOE    ,INCRD    ,INDFM
     ; ,INDSP    ,INFOD    ,INPOR    ,INPWR    ,INSTG    ,INTRA
     ; ,IOEQT    ,IOFLSAT  ,IOTRS    ,IOWAR    ,IUPAR    ,MAINF
     ; ,NPAREL   ,NROW     ,NUMEL    ,NZARR    ,NZCOE    ,NZCRD
     ; ,NZDFM    ,NZDSP    ,NZFOD    ,NZPOR    ,NZSTG    ,NZTRA
     ; ,CFPAREL  ,FILENAME)

*****************************************************************************
* PURPOSE
*     Reads element coefficients of flow and/or transport equations
*
* DESCRIPTION
*     Reads element coefficients of flow and/or transport equations
*
* EXTERNAL VARIABLES: ARRAYS
*
*  CFPAREL                Array containing node coefficient of element j       
*                         corresponding to INpar index parameter zone.          
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  INARR                  Index for areal recharge                              
*                         in array variables (LXPAREL and CFPAREL)              
*  INARRT                 Index for transient areal recharge                    
*                         in array variables (LXPAREL and CFPAREL)              
*  INCOE                  Index for external concentration (elements)           
*                         in array variables (LXPAREL and CFPAREL)              
*  INCRD                  Index for retardation coefficient                     
*                         in array variables (LXPAREL and CFPAREL)              
*  INDFM                  Index for molecular diffusion                         
*                         in array variables (LXPAREL and CFPAREL)              
*  INDSP                  Index for disperssiviy                                
*                         in array variables (LXPAREL and CFPAREL)              
*  INFOD                  Index for first order decay coefficient               
*  INPOR                  Index for porosity                                    
*                         in array variables (LXPAREL and CFPAREL)              
*  INPWR                  Allows writing on MAIN FILE                           
*  INSTG                  Index for storage coefficient                         
*                         in array variables (LXPAREL and CFPAREL)              
*  INTRA                  Index for transmissivity                              
*                         in array variables (LXPAREL and CFPAREL)              
*  IOEQT                  Type of problem to be solved                          
*  IOFLSAT                Indicates the possibility that one part of the domain 
*                         reaches unsaturated state.                            
*  IOTRS                  Flow regime                                           
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUPAR                  Unit number of file PAR                               
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NPAREL                 Number of element parameters in current problem       
*  NROW                   Current record number                                 
*  NUMEL                  Number of elements                                    
*  NZARR                  Number of areal recharge zones                        
*  NZCOE                  Number of external concentration zones                
*  NZCRD                  Number of retardation coefficient zones               
*  NZDFM                  Number of molecular difusion zones                    
*  NZDSP                  Number of dispersivity zones                          
*  NZFOD                  Number of zones of first order decay                  
*  NZPOR                  Number of porosity zones                              
*  NZSTG                  Number of storage coefficient zones                   
*  NZTRA                  Number of transmissivity zones                        
*
* INTERNAL VARIABLES: SCALARS
*
*  DARR                   Steady state recharge default value.                  
*  DARRT                  Transient recharge default value.                     
*  DCOE                   Steady state external concentration default value.    
*  DCRD                   Retardation coefficient default value.                
*  DDFM                   Molecular diffusion default value.                    
*  DDSP                   Dispersivity default value.                           
*  DFOD                   First order decay default value.
*  DPOR                   Porosity default value.                               
*  DSTG                   Storage default value. 
*  DTRA                   Transmissivity default value.
*  LEAUX                  Auxiliar string where the last read line of the       
*                         current input file is stored                          
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ASS_REAL_VAL           Assigns the read or the default value
*  ERROR                  Writes the current error message and error number     
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
*****************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*100 LEAUX,LEEL,FILENAME(20)*20

       DIMENSION CFPAREL(NUMEL,NPAREL)

C------------------------- FIRST EXECUTABLE STATEMENT.

       IF (IOEQT.NE.2) THEN

*_______________________Reads flow element's coefficients

*_______________________Reads defaults (CARD C2.1)

          LEAUX=LEEL(FILENAME,IUPAR,MAINF,NROW,INPWR)          
          READ(LEAUX,1000,ERR=9000) DTRA,DSTG,DARR,DARRT
 1000     FORMAT (5F10.0)
          IF (INPWR.NE.0) THEN 
             WRITE(MAINF,3000) 
 3000        FORMAT(////,
     ;              10X,'ELEM.COEFF. FOR FLOW EQUATION',/,
     ;              10X,'=============================',/)
             WRITE(MAINF,3100) DTRA,DSTG,DARR,DARRT
 3100        FORMAT(/,5X,'DEFAULT VALUES ',/
     ;                5X,'--------------',/
     ;       5X,'TRANSMISIVITY .......... =',G10.3,/
     ;       5X,'STORAGE COEFF. ......... =',G10.3/
     ;       5X,'RECHARGE COEFF(S.S). ... =',G10.3/
     ;       5X,'RECHARGE COEFF(TRANS) .. =',G10.3/)   
          END IF
            
*_______________________Assign default values

              IF (DTRA.NE.0. AND. NZTRA.NE.0) THEN 
                DO L=1,NUMEL
                  CFPAREL(L,INTRA)=DTRA
                END DO
              END IF 

              IF (DSTG.NE.0) THEN
                IF (IOTRS.EQ.0) THEN
                  CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                      'ST. STATE CONDITIONS AND NON-ZERO '//
     ;                      'STORAGE ELEMENT COEFFICIENT',
     ;                       NROW,0,IUPAR,0,0.00)
                ELSE
                  IF (NZSTG.NE.0) THEN 
                    DO L=1,NUMEL
                      CFPAREL(L,INSTG)=DSTG 
                    END DO   
                  END IF
                END IF
              END IF

              IF (DARR.NE.0 .AND. NZARR.NE.0) THEN
               
                DO L=1,NUMEL 
                  CFPAREL(L,INARR)=DARR
                END DO
                 
              END IF
                 
              IF (DARRT.NE.0 .AND. NZARR.NE.0) THEN

                DO L=1,NUMEL
                  CFPAREL(L,INARRT)=DARRT
                END DO
                
              END IF   

             NOLD=0

*_______________________Starts loop reading element coefficient (card C2.2)

C------------------------- Initializes some auxiliar variables

       AUXST=0.D0
       AUXAR=0.D0 
       AUXART=0.D0

         NE=0
               DO WHILE (NE.LT.NUMEL)  

                LEAUX=LEEL(FILENAME,IUPAR,MAINF,NROW,INPWR)          
                READ(LEAUX,1100,ERR=9100) NE,TRA,STG,ARR,ARRT
 1100           FORMAT (I5,4F10.0)   
                
                IF (NOLD.EQ.0 .AND. ((NE.NE.0 .AND. INPWR.NE.0) .OR. 
     .          INPWR.GT.1))  WRITE(MAINF,3200) 
 3200     FORMAT(//'   ELEM   TRANSM     STOR      ARR       ARRT  '/)

                IF (NE.GT.NUMEL .OR. NE.LE.0) THEN
                   CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                       'ELEMENT NUMBER IS OUT OF '//
     .                       'BOUNDS WHEN READING TRANSP. PARAM.'//
     .                       'ELEM. COEFF. ',NROW,1,IUPAR,1,5.01)

*_______________________Checks increasing order of elements

                ELSE

                   IF (NE.LE.NOLD)
     ;                CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;                'ELEMENT OUT OF SEQUENCE (FLOW ELEM. COEFF)',
     ;                NROW,1,IUPAR,1,5.02)
                   IF (NZTRA.NE.0)
     ;                CALL ASS_REAL_VAL(TRA,DTRA,CFPAREL(NE,INTRA))
                
                   IF (NZARR.NE.0)
     ;                CALL ASS_REAL_VAL(ARR,DARR,CFPAREL(NE,INARR))
                      
                   IF (NZARR.NE.0)
     ;                CALL ASS_REAL_VAL(ARRT,DARRT,CFPAREL(NE,INARRT))

                   IF (IOTRS.NE.0) THEN
                     IF (NZSTG.NE.0)
     ;                 CALL ASS_REAL_VAL(STG,DSTG,CFPAREL(NE,INSTG))

                   ELSE 

                      IF (STG.NE.0)
     ;                  CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                  'ST. STATE CONDITIONS AND NON-ZERO '//
     ;                  'STORAGE ELEMENT COEFFICIENT ',
     ;                  NROW,3,IUPAR,0,0.00)
                   END IF              
                END IF

*_______________________Writes missing elements

                IF (INPWR.GT.1) THEN
                   DO N=NOLD+1,NE-1
                      IF (IOTRS.NE.0) THEN
                         AUXST=CFPAREL(N,INSTG)
                      END IF
                      IF (NZARR.NE.0) AUXAR=CFPAREL(N,INARR)
                      IF (NZARR.NE.0) AUXART=CFPAREL(N,INARRT)
                       WRITE(MAINF,3300) N,CFPAREL(N,INTRA),AUXST,
     ;                 AUXAR,AUXART
 3300                  FORMAT (I5,4G10.3)   
                   END DO
                END IF

                   NOLD=NE
                
*_______________________Writes last element read                

                IF (INPWR.NE.0) THEN  
                   IF (IOTRS.NE.0) THEN
                      AUXST=CFPAREL(NE,INSTG)
                   END IF
                   IF (NZARR.NE.0) AUXAR=CFPAREL(NE,INARR)
                   IF (NZARR.NE.0) AUXART=CFPAREL(NE,INARRT)
                   WRITE(MAINF,3300) NE,CFPAREL(NE,INTRA),AUXST,AUXAR,
     ;                               AUXART
                END IF
             END DO    !SEARCH NEXT ELEMENT
        END IF

        IF (IOEQT.NE.1) THEN

*_______________________Reads elements coefficients for transport 

*_______________________Reads default values (CARD C3.1)
                             
          LEAUX=LEEL(FILENAME,IUPAR,MAINF,NROW,INPWR)
          READ(LEAUX,1200,ERR=9200) DDSP,DDFM,DPOR,DCRD,DCOE,DFOD
 1200     FORMAT (6F10.0)    
          IF (INPWR.NE.0) THEN 
             WRITE(MAINF,3400) 
 3400        FORMAT(////,10X,'ELEM.COEFF. FOR TRANSP. EQUATION',/
     ;                   10X,'================================',/)
             WRITE(MAINF,3500) DDSP,DDFM,DPOR,DCRD,DCOE,DFOD
 3500        FORMAT(/,5X,'DEFAULT VALUES',/
     ;                5X,'--------------',/
     ;       5X,'DISPERSIVITY ......... =',G10.3,/
     ;       5X,'MOLECULAR DIFFUSION .. =',G10.3,/
     ;       5X,'POROSITY ............. =',G10.3,/
     ;       5X,'RETARD. COEFFICIENT .. =',G10.3,/
     ;       5X,'EXTERNAL CONCENT. .... =',G10.3,/
     ;       5X,'FIRST ORDER DEC. COEF. =',G10.3)
          END IF

C------------------------- Assign default values

              IF (DDSP.NE.0 .AND. NZDSP.NE.0) THEN
                
                DO L=1,NUMEL
                  CFPAREL(L,INDSP)=DDSP        
                END DO

              END IF

              IF (DDFM.NE.0 .AND. NZDFM.NE.0) THEN
                
                DO L=1,NUMEL
                  CFPAREL(L,INDFM)=DDFM  
                END DO

              END IF
               
              IF (DPOR.NE.0. AND .NZPOR.NE.0) THEN

                  DO L=1,NUMEL
                  CFPAREL(L,INPOR)=DPOR
                  END DO
                    
              END IF

              IF (DCRD.NE.0 .AND. NZCRD.NE.0) THEN

                DO L=1,NUMEL
                  CFPAREL(L,INCRD)=DCRD
                END DO

              END IF       
      
              IF (DCOE.NE.0. AND .NZCOE.NE.0) THEN

                DO L=1,NUMEL
                  CFPAREL(L,INCOE)=DCOE
                END DO

              END IF   
              
               IF (DFOD.NE.0. AND .NZFOD.NE.0) THEN

                DO L=1,NUMEL
                  CFPAREL(L,INFOD)=DFOD
                END DO

              END IF   

*______________________Starts loop reading elem. coeff. for transp. (card C3.2)

          NOLD=0
          NE=0

C------------------------- Initializes some auxiliar variables

       AUXDS=0.D0
       AUXDF=0.D0
       AUXCR=0.D0
       AUXFO=0.D0

              DO WHILE (NE.LT.NUMEL) 
                LEAUX=LEEL(FILENAME,IUPAR,MAINF,NROW,INPWR)          
                READ(LEAUX,1300,ERR=9300) NE,DSP,DFM,POR,CRD,COE,FOD
 1300           FORMAT (I5,6F10.0) 

C------------------------- Writes header on MAIN file
 
               IF (NOLD.EQ.0 .AND. ((NE.NE.0 .AND. INPWR.NE.0) .OR. 
     .         INPWR.GT.1))  WRITE(MAINF,3600) 
 3600    FORMAT(//
     ;'   ELEM    DSP       DFM       POR       CRD      COE      FOD'/)
                   
*_______________________Checks if element numbers are in bounds

                     IF (NE.GT.NUMEL .OR. NE.LE.0) THEN
                        CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;                    'ELEMENT NUMBER IS OUT OF '//
     ;                    'BOUNDS WHEN READING TRANSP. PARAM.'//
     ;                    'ELEM. COEFF. ',NROW,1,IUPAR,1,5.03)

*_______________________Checks increasing order of elements

                     ELSE
                        IF (NE.LE.NOLD) 
     ;                    CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;                                'ELEMENT OUT '//
     ;                                'OF SEQUENCE (FLOW ELEM. COEFF) ',
     ;                                 NROW,1,IUPAR,1,5.04)    

*_______________________Assign default values through subroutine ASS_REAL_VAL

                        IF (NZDSP.NE.0)
     ;                    CALL ASS_REAL_VAL (DSP,DDSP,CFPAREL(NE,INDSP))

                        IF (NZDFM.NE.0)                      
     ;                    CALL ASS_REAL_VAL (DFM,DDFM,CFPAREL(NE,INDFM))

                        IF (NZPOR.NE.0)                           
     ;                    CALL ASS_REAL_VAL (POR,DPOR,CFPAREL(NE,INPOR))

                        IF (NZCRD.NE.0)
     ;                    CALL ASS_REAL_VAL (CRD,DCRD,CFPAREL(NE,INCRD))
                           
                        IF (NZFOD.NE.0)
     ;                    CALL ASS_REAL_VAL (FOD,DFOD,CFPAREL(NE,INFOD))
        
                        IF (NZCOE.NE.0)
     ;                    CALL ASS_REAL_VAL (COE,DCOE,CFPAREL(NE,INCOE))
                                   
                     END IF

*_______________________Writes missing nodes

                   IF (INPWR.GT.1) THEN
                      DO N=NOLD+1,NE-1
                         IF (NZDSP.NE.0) AUXDS=CFPAREL(N,INDSP)
                         IF (NZDFM.NE.0) AUXDF=CFPAREL(N,INDFM)
                         IF (NZCRD.NE.0) AUXCR=CFPAREL(N,INCRD)
                         IF (NZFOD.NE.0) AUXFO=CFPAREL(N,INFOD)
                         WRITE(MAINF,3700) N,AUXDS,AUXDF,
     ;                                     CFPAREL(N,INPOR),
     ;                                     AUXCR,CFPAREL(N,INCOE),AUXFO
 3700                    FORMAT (I5,6G10.3)  

                      END DO
                   END IF

*_______________________Writes last node read

                   IF (INPWR.NE.0) THEN
                      IF (NZDSP.NE.0) AUXDS=CFPAREL(NE,INDSP)
                      IF (NZDFM.NE.0) AUXDF=CFPAREL(NE,INDFM)
                      IF (NZCRD.NE.0) AUXCR=CFPAREL(NE,INCRD)
                      IF (NZFOD.NE.0) AUXFO=CFPAREL(NE,INFOD)
                      WRITE(MAINF,3700) NE,AUXDS,AUXDF,CFPAREL(NE,INPOR)
     ;                                 ,AUXCR,CFPAREL(NE,INCOE),AUXFO
                   END IF

                   NOLD=NE
              END DO   !Next element
       
       ELSE IF(IOFLSAT.NE.0)THEN  

*_______________________For unsaturated flow have to read default values of 
*_______________________coefficient of porosity (card C3.1) 
                                   
          LEAUX=LEEL(FILENAME,IUPAR,MAINF,NROW,INPWR)
          READ(LEAUX,1200,ERR=9200) DDSP,DDFM,DPOR,DCRD,DCOE,DFOD

          IF(INPWR.NE.0)WRITE(MAINF,3800)DPOR
 3800     FORMAT(/,' DEFAULT FOR POROSITY COEFFICIENT=',E10.4,/)

*_______________________Default values
          
          IF (NZPOR. NE.0) THEN  
            DO L1=1,NUMEL
              CFPAREL(L1,INPOR)=DPOR
            END DO
          END IF

          NE=0

          DO WHILE(NE.LT.NUMEL)

            LEAUX=LEEL(FILENAME,IUPAR,MAINF,NROW,INPWR)
            READ(LEAUX,1300,ERR=9300) NE,DSP,DFM,POR,CRD,COE,FOD
            
              IF(NE.LE.0.OR.NE.GT.NUMEL) THEN
                 CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                     'ELEMENT NUMBER IS OUT OF '//
     .                     'BOUNDS WHEN READING TRANSP. PARAM.'//
     .                     'ELEM. COEFF. ',NROW,1,IUPAR,1,5.01)

              ELSE
  
                 CFPAREL(NE,INPOR)=POR

              END IF

          END DO !NEXT ELEMENT

          IF (INPWR.NE.0)THEN
            WRITE(MAINF,3900)
 3900       FORMAT('    ELEMENT   POROSITY COEFFICIENT')
         
            DO L1=1,NUMEL

              WRITE(MAINF,4000)L1,CFPAREL(L1,INPOR)
 4000         FORMAT(4X,I5,7X,E12.5)

            END DO
          ENDIF

       ENDIF

       RETURN

 9000  CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;            'GENERIC FORTRAN ERROR WHEN READING FLOW DEFAULTS ',
     .            0,0,IUPAR,1,5.05)

       RETURN
       
 9100  CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;            'GENERIC FORTRAN ERROR WHEN READING FLOW '//
     ;            'ELEMENT COEFF. ',NROW,0,IUPAR,1,5.06)

       RETURN
       
 9200  CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;            'GENERIC FORTRAN ERROR WHEN READING TRANSPORT '//
     ;            'DEFAULTS ',NROW,0,IUPAR,1,5.07)

       RETURN

 9300  CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;             'GENERIC FORTRAN ERROR WHEN READING TRANSPORT '//
     ;             'ELEMENT COEFF. ',NROW,0,IUPAR,1,5.08)
       
       END                                                 

