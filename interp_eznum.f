      SUBROUTINE INTERP_EZNUM
     ; (INPWR    ,IOEQT    ,IOFLSAT  ,IOTRS    ,LDARR    ,LDARRT
     ; ,LDCOE    ,LDCRD    ,LDDFM    ,LDDSP    ,LDFOD    ,LDPOR
     ; ,LDSTG    ,LDTRA    ,MAINF    ,N        ,NE       ,NOLD
     ; ,NUMEL    ,NZARR    ,NZCOE    ,NZCRD    ,NZDFM    ,NZDSP
     ; ,NZFOD    ,NZPOR    ,NZSTG    ,NZTRA    ,LDIM     ,LTYPE
     ; ,LXARR    ,LXARRT   ,LXCOE    ,LXCRD    ,LXDFM    ,LXDSP
     ; ,LXFOD    ,LXPOR    ,LXSTG    ,LXTRA)

*****************************************************************************
* PURPOSE 
*
*     This subroutine interpolates zone numbers for missing elements in the 
*     input file
*
* DESCRIPTION
*
*     This subroutine starts with a loop that goes through all missing
*     elements between NOLD and NE. Each missing element zone number is
*     assigned through subrutine ASS_EVAL: if  a default value for a given
*     parameter zone number equals 0, current element parameter zone number
*     is assigned that corresponding to NOLD element, otherwise is assigned
*     the default value LDPAR.
* 
* EXTERNAL VARIABLES: ARRAYS 
*     
*  LDIM                   Vector containing fisical dimension of j-th element   
*  LTYPE                  Vector containing type for element j                  
*  LXARR                  Areal recharge (steady) zone number at a given element
*  LXARRT                 Areal recharge (trans.) zone number at a given element
*  LXCOE                  External concentration zone number at a given element 
*  LXCRD                  Retardation zone number at a given element            
*  LXDFM                  Molecular diffusion zone number at a given element    
*  LXDSP                  Dispersivity zone number at a given element           
*  LXFOD                  First order decay zone number at a given element      
*  LXPOR                  Porosity zone number at a given element               
*  LXSTG                  Storage coefficient zone number at a given element    
*  LXTRA                  Transmissivity zone number at a given element         
*
* EXTERNAL VARIABLES: SCALARS  
*
*  INPWR                  Allows writing on MAIN FILE                           
*  IOEQT                  Type of problem to be solved                          
*  IOFLSAT                Indicates the possibility that one part of the domain 
*                         reaches unsaturated state.                            
*  IOTRS                  Flow regime                                           
*  LDTRA                  Default value for transmisivity zone number
*  LDARR                  Default value for areal recharge
*  LDARRT                 Default value for transient areal recharge
*  LDSTG                  Default value for storage coeficient    
*  LDDSP                  Default value for dispersivity
*  LDDFM                  Default value for matrix diffusion
*  LDCOE                  Default value for external concentration
*  LDPOR                  Default value for porosity
*  LDCRD                  Default value for retardation coefficient
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NUMEL                  Number of elements                                    
*  NZARR                  Number of areal recharge zones                        
*  NZCOE                  Number of external concentration zones                
*  NZCRD                  Number of retardation Coefficient zones               
*  NZDFM                  Number of molecular difusion zones                    
*  NZDSP                  Number of dispersivity zones                          
*  NZFOD                  Number of zones of first order decay                  
*  NZPOR                  Number of porosity zones                              
*  NZSTG                  Number of storage Coefficient zones                   
*  NZTRA                  Number of transmissivity zones                        
*  NOLD                   Last element read before the jump 
*  N                      NOLD+1
*  NE                     Last element read.
*
* INTERNAL VARIABLES: SCALARS 
*
*  LOLTRA                 Transmissivity zone number. Auxiliary var.
*  LOLARR                 Areal recharge zone number. Auxiliary var.
*  LOLARRT                Transient areal recharge zone number. Auxiliary var.
*  LOLSTG                 Storage coefficient zone number. Auxiliary var.
*  LOLDSP                 Dispersivity zone number. Auxiliary var.
*  LOLDFM                 Matrix diffusion zone number. Auxiliary var.
*  LOLCOE                 External concentration zone number. Auxiliary var.
*  LOLPOR                 Porosity zone number. Auxiliary var.
*  LOLCRD                 Retardation coeficient zone number. Auxiliary var.
* 
* SUBROUTINES REFERENCED
*
*  ASS_EVAL               Assigns the zone number of one element
*  LDIMEN                 Computes the dimension of a given element         
*
* HISTORY
*
*     SCR 15-apr-1997     First coding
*     AMS      1-1998     Revision
*     AMS      4-1999     Correction of the dimension of LXFOD array
*     AMS     10-2000     Inicialization of variables LOL*, and correction
*                         of an erroneous assignment
*
*****************************************************************************
        
       IMPLICIT REAL*8 (A-H,O-Z)       
        
       DIMENSION LDIM(NUMEL),LTYPE(NUMEL),
     ;      LXTRA(NUMEL),LXARR(NUMEL),LXARRT(NUMEL),
     ;      LXSTG(NUMEL),LXDSP(NUMEL),LXDFM(NUMEL),LXCOE(NUMEL),
     ;      LXPOR(NUMEL),LXCRD(NUMEL),LXFOD(NUMEL)
        
C------------------------- FIRST EXECUTABLE STATEMENT.

C------------------------- Initializes auxiliary variables

       LOLTRA=  0
       LOLSTG=  0
       LOLARR=  0
       LOLARRT= 0
       LOLDSP=  0
       LOLDFM=  0
       LOLPOR=  0
       LOLCRD=  0
       LOLCOE=  0
       LOLFOD=  0

       DO NN=N,NE-1 !N=NOLD+1

*_______________________Put in the zone numbers
*_______________________if defaults are not zero

          IF (IOEQT.NE.2) THEN
             CALL ASS_EVAL
     ;            (NZTRA,LDTRA,LXTRA(NOLD),LXTRA(NN))

C------------------------- Stores element dimension of every transmissivity 
C------------------------- zone. Array LDIM is used temporary for this function.
C------------------------- Array LDIM will be redefined later for its 
C------------------------- definitely use (to store the dimension of every 
C------------------------- element.

             NZTR=LXTRA(NN)
             IF (NZTR.NE.0) LDIM(NZTR)=LDIMEN(LTYPE(NN))

             CALL ASS_EVAL
     ;            (NZARR,LDARR,LXARR(NOLD),LXARR(NN))
                
             IF (IOTRS.NE.0) THEN

                CALL ASS_EVAL (NZSTG,LDSTG,LXSTG(NOLD),LXSTG(NN))

             ENDIF
             CALL ASS_EVAL (NZARR,LDARRT,LXARRT(NOLD),LXARRT(NN))
          ENDIF

          IF (IOEQT.NE.1) THEN
         
             CALL ASS_EVAL (NZDSP,LDDSP,LXDSP(NOLD),LXDSP(NN))

             CALL ASS_EVAL (NZDFM,LDDFM,LXDFM(NOLD),LXDFM(NN))

             CALL ASS_EVAL (NZPOR,LDPOR,LXPOR(NOLD),LXPOR(NN))

             CALL ASS_EVAL (NZCRD,LDCRD,LXCRD(NOLD),LXCRD(NN))

             CALL ASS_EVAL (NZCOE,LDCOE,LXCOE(NOLD),LXCOE(NN))

             CALL ASS_EVAL (NZFOD,LDFOD,LXFOD(NOLD),LXFOD(NN))
     
          ELSE IF(IOFLSAT.NE.0)THEN
         
             CALL ASS_EVAL (NZPOR,LDPOR,LXPOR(NOLD),LXPOR(NN))      

          ENDIF

*_______________________Defines auxiliary variables to write
*_______________________last interpolated element zone numbers

          IF (INPWR.NE.0) THEN  
             IF (IOEQT.NE.2) THEN
                LOLTRA=LXTRA(NN)
                IF (NZARR.NE.0) LOLARR=LXARR(NN)
                IF (IOTRS.NE.0) THEN
                   LOLSTG=LXSTG(NN)
                END IF
                IF (NZARR.NE.0) LOLARRT=LXARRT(NN)
             END IF
             IF (IOEQT.NE.1) THEN
                IF (NZDSP.NE.0) LOLDSP=LXDSP(NN)
                IF (NZDFM.NE.0) LOLDFM=LXDFM(NN)
                IF (NZPOR.NE.0) LOLPOR=LXPOR(NN)
                IF (NZCRD.NE.0) LOLCRD=LXCRD(NN)
                IF (NZCOE.NE.0) LOLCOE=LXCOE(NN)
                IF (NZFOD.NE.0) LOLFOD=LXFOD(NN)
             ELSE IF(IOFLSAT.NE.0)THEN
                IF (NZPOR.NE.0) LOLPOR=LXPOR(NN)
             ENDIF

*_______________________Writes last interpolated element zone numbers
*_______________________on MAIN file

             WRITE (MAINF,3000) NN,LOLTRA,LOLSTG,LOLARR,LOLARRT,
     ;                 LOLDSP,LOLDFM,LOLPOR,LOLCRD,LOLCOE,LOLFOD
 3000        FORMAT(11I5,2F10.0)     
          ENDIF   
             
       ENDDO !NEXT ELEMENT
        
       RETURN
       END
