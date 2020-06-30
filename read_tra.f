       SUBROUTINE READ_TRA
     ;(IDIMWGT   ,IERROR    ,INPWR      ,IOINV      ,IOLGTRA    ,IOPBLI
     ;,IOWAR     ,IPARDET   ,ISOT       ,IOTIM      ,IUPAR      ,IZTRAC
     ;,MAINF     ,MXGRPZN   ,NFNL       ,NGROUP_ZN  ,NPAR       ,NPFNL      
     ;,NUMEL     ,NZTRA     ,WEIGHT     ,FILENAME   ,INDPAR     ,INTRAC     
     ;,IOPTLOG   ,IOPT_GS   ,IPNT_END   ,IPNT_PAR   ,IPNT_START ,ISOZ      
     ;,IVTRAGRP  ,LDIM      ,NFNLTRA    ,NFTTRA     ,STTRA      ,TRAC  
     ;,TRAM      ,TRAZ      ,WGT_TRA    ,WGT_UNK    ,NROW       ,PARNAME
     &,INAME)
                         
*****************************************************************************
* PURPOSE
*      Reads and checks zonal values of transmissivity
*
* DESCRIPTION
*      Reads and checks zonal values of transmissivity
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  INDPAR                 Array of 0's and 1's. 0 means that the parameter is   
*                         estimated aritmethically, otherwise logarithmically.  
*  INTRAC                 Array containing the location of first transmissivity 
*                         zone (for each tensor component) in array             
*                         variables PARC, PARM, STPAR ... minus 1               
*  ISOZ                   Anisotropy of every transmissivity zone               
*  IVTRA                  Transmissivity estimation index (0=no estimation)     
*  LDIM                   Vector containing fisical dimension of j-th element   
*                         Temporarily this array stores the dimension 
*                         of one element of each transmissivity zone (it is 
*                         assumed that all elements belonging to the same 
*                         transmissivity zone have the same dimension)
*  NFNLTRA                Transmissivity nonlinear functions                    
*  NFTTRA                 Transmissivity time functions                         
*  STTRA                  Transmissivity standard deviation                     
*  TRAC                   Computed Transmissivity zonal values                  
*  TRAM                   Prior information Transmissivity zonal values         
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  IOINV                  Inverse problem option                                
*  IOLGTRA                Transmissivity log-scaling index                      
*  IOTRS                  Flow regime                                           
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IPAR                   Number of parameters to estimate                                                      
*  ISOT                   Maximum hydraulic conductivity tensor anisotropy      
*                         degree in the problem.                                
*  IUPAR                  Unit number of file PAR                               
*  IZTRAC                 NZTRA*MAX(IODIM,ISOT). Used to dimension TRAC array
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NPAR                   Total number of parameters to be estimated            
*  NPFNL                  Counts the number of nonlinear functions read
*  NROW                   Current record number                                 
*  NUMEL                  Number of elements                                    
*  NZTRA                  Number of transmissivity zones                        
*
* INTERNAL VARIABLES: SCALARS
*
*  I                      Dummy counter
*  ICOMPO                 Pointer
*  IDIM                   Dimension of the current transmissivity zone                                                      
*  IGRP1                  Read group of zones
*  ISOZ1                  Read anisotropy
*  IVTRA1                 Read transmissivity estimation index
*  IVVAR                  Not used                                                      
*  j                      Dummy counter
*  LEAUX                  Auxiliar string where the last read line of the       
*                         current input file is stored                          
*  N                      Dummy counter of zones of actual param. type
*  NFNLTRA1               Read nonlinear transmissivity function number
*  NFTTRA1                Read time transmissivity function number
*  NROW                   Current record number                                 
*  NZ                     Zone number as read
*  STTRA1                 Read standard deviation of transmissivity 
*  TRAM1                  Read transmissivity prior information value of the 
*                         current zone
*  TRAZ1                  Read transmissivity initial value of the current zone
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  CHECK_PAR              Checks transmissivity zonal data
*  ERROR                  Writes the current error message and error number     
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
* HISTORY
*
*     SCR      5-1997     First version
*     AMS      1-1998     Revision. Addition of header
*     AMS     12-1998     Corrects update of IVTRA array 
*     AAR      7-2003     Revision
*
*****************************************************************************

       IMPLICIT NONE
                                                             ! Integer external
       INTEGER*4 IZTRAC,NZTRA,NUMEL,IDIMWGT,NPAR,INPWR,MAINF,IUPAR
     ;          ,IERROR,IOWAR,ISOT,IOINV,IOPBLI,IOTIM,MXGRPZN,NFNL,NPFNL
     ;          ,IOLGTRA,IPARDET,NGROUP_ZN
     ;          ,NFTTRA(IZTRAC),ISOZ(NZTRA),LDIM(NUMEL)
     ;          ,NFNLTRA(IZTRAC),INDPAR(NPAR),INTRAC(6)
     ;          ,IVTRAGRP(IZTRAC),IPNT_START(IZTRAC)
     ;          ,IPNT_END(IZTRAC),IOPT_GS(MXGRPZN,20)
     ;          ,IPNT_PAR(IZTRAC*IDIMWGT),IOPTLOG(IZTRAC)
                                                                ! Real external
       REAL*8 TRAZ(IZTRAC),STTRA(IZTRAC),TRAM(NPAR),WGT_UNK(NPAR)
     ;       ,WGT_TRA(IZTRAC*IDIMWGT),TRAC(NPAR),WEIGHT
                                                             ! Integer internal
       INTEGER*4 N,NROW,NZ,ISOZ1,IVTRA1,NFNLTRA1,NFTTRA1,IGRP1,IDIM,J,I
     ;          ,ICOMP, INAME
                                                                ! Real internal
       REAL*8 TRAZ1,TRAM1,STTRA1
                                                                   ! Characters
       CHARACTER FILENAME(18)*20,LEEL*100,LEAUX*100,PARNAME(NPAR)*4        
       
C------------------------- FIRST EXECUTABLE STATEMENT.
C------------------------- Writes title in MAIN file

       IF (INPWR.NE.0) WRITE(MAINF,3000)
 3000  FORMAT(////20X,'  ZONES INFORMATION',/20X,17('*'))


       IF (INPWR.NE.0) WRITE(MAINF,3100)
 3100  FORMAT(///,10X,'TRANSMISIVITY   ZONES',/,
     ;        10X,'-------------------',/,
     ;' ZONE ISOZ   COMPUTED    IV   ST.DESV.    MEASURED  NFNL',
     ;'  NFT IGRP')

C------------------------- Starts transmissivity zones loop

       DO N=1,NZTRA

C------------------------- Reads first tensor component
          IGRP1 = 0
          LEAUX=LEEL(FILENAME,IUPAR,MAINF,NROW,INPWR)          
          READ(LEAUX,1000,ERR=9000) NZ,ISOZ1,TRAZ1,IVTRA1,
     ;                             STTRA1,TRAM1,NFNLTRA1,NFTTRA1,IGRP1

 1000     FORMAT (2I5,F10.0,I5,2F10.0,3I5)                      
 
C------------------------- Checks transmissivity zone number sequence

          IF (NZ.NE.N) CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;       'TRANSMISSIVITY ZONE NUMBER IS OUT OF SEQUENCE ',
     ;        NROW,1,IUPAR,1,6.01)
 
C------------------------- Checks transmissivity anisotrpy degree and
C------------------------- assigns value to missing tensor components

          IF (ISOZ1.GT.ISOT) CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;       'ZONAL ISOTROPY INDEX IS TOO LARGE ',NROW,2,IUPAR,1,6.02)

          ISOZ(NZ)=ISOZ1 

C------------------------- Temporarily LDIM array stores the dimension
C------------------------- of one element of each transmissivity zone (it is
C------------------------- assumed that all elements belonging to the same
C------------------------- transmissivity zone have the same dimension)

          IDIM=LDIM(NZ)

C------------------------- Checks transmissivity read parameters 
C------------------------- through CHEK_PAR subroutine. 

          CALL CHECK_PAR 
     ;(6.02     ,IERROR  ,IGRP1     ,IOINV   ,IOPBLI   ,IOTIM    
     ;,IOWAR    ,IUPAR   ,IVTRA1    ,MAINF   ,MXGRPZN  ,NFNL    
     ;,NFNLTRA1 ,NFTTRA1 ,NGROUP_ZN ,NPFNL   ,STTRA1  
     ;,'   TRANSMISSIVITY   '       ,FILENAME ,NROW)

          IF (ISOZ1.EQ.1) THEN 

C------------------------- Assigns the same values to all tensor component 
C------------------------- (isotropic case). TRAZs are always assigned. If later
C------------------------- due to geostat. considerations, they are calculated
C------------------------- these "read" values will be substituted

             DO J=1,IDIM
                ICOMP=INTRAC(J)+NZ
*                IF (IOPT_GS(IGRP1,2).EQ.0) TRAZ(ICOMP)=TRAZ1
                TRAZ(ICOMP)=TRAZ1   ! Always assigned. 
                NFTTRA(ICOMP)=NFTTRA1
                NFNLTRA(ICOMP)=NFNLTRA1
                IVTRAGRP(ICOMP)=IGRP1
                IOPTLOG(ICOMP)=IOLGTRA

  ! If inverse problem is solved and param. is estimated deterministically

             ENDDO

             IF (IOINV.GT.0 .AND.
     ;             IVTRA1.NE.0 .AND. IOPT_GS(IGRP1,2).EQ.0) THEN
                ICOMP=INTRAC(1)+NZ
                IPARDET=IPARDET+1
                IF (IPARDET.GT.NPAR) CALL ERROR 
     ;(IERROR,IOWAR,MAINF,FILENAME,
     ;'NUMBER OF PARAMETERS TO BE ESTIMATED (PAR FILE) IS GREATER'//
     ;' THAN THE DEFINED AT DIM FILE',NROW,0,IUPAR,2,0.07)

                IPNT_START(ICOMP)=(ICOMP-1)*IDIMWGT+1
                IPNT_END(ICOMP)=IPNT_START(ICOMP)
                IPNT_PAR((ICOMP-1)*IDIMWGT+1)=IPARDET
                WGT_TRA((ICOMP-1)*IDIMWGT+1)=1.0D0
                IF(TRAM1.EQ.0D0) TRAM1=TRAZ1

                IF (IOLGTRA.EQ.1) THEN  ! Log-estimation
                   TRAM(IPARDET)=DLOG10(TRAM1)
                   TRAC(IPARDET)=DLOG10(TRAZ1)
                ELSE
                   TRAM(IPARDET)=TRAM1
                   TRAC(IPARDET)=TRAZ1
                END IF

                WGT_UNK(IPARDET)=WEIGHT

                IF (STTRA1.LE.0D0) THEN
                  CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;      'ESTIMATED ZONAL PARAM. AND STANDARD DEVIATION IS ZERO '//
     ;      'OR NEGATIVE. IT IS SET TO 1.0',NROW,0,IUPAR,1,6.16)
                   STTRA1=1.0D0
                END IF

                STTRA(ICOMP)=STTRA1
                INDPAR(IPARDET)=0
                IF (IOLGTRA.EQ.1) INDPAR(IPARDET)=IOLGTRA

c------------Stores parameter name. Used to write PSH and PSC.
                INAME = INAME +1
                WRITE(PARNAME(INAME),11) 'TX',N
   11           FORMAT(A2,I2)

             END IF

          ELSE    ! Anysotropic case

             ICOMP=INTRAC(1)+NZ
             IF (IOPT_GS(IGRP1,2).EQ.0) TRAZ(ICOMP)=TRAZ1
             NFTTRA(ICOMP)=NFTTRA1
             NFNLTRA(ICOMP)=NFNLTRA1
             IVTRAGRP(ICOMP)=IGRP1
             IOPTLOG(ICOMP)=IOLGTRA
             IF (IOINV.GT.0 .AND. IVTRA1.NE.0 .AND. 
     ;           IOPT_GS(IGRP1,2).EQ.0) THEN  
 
*____________________ Parameter is estimated deterministically

                IPARDET=IPARDET+1
                IF (IPARDET.GT.NPAR) CALL ERROR 
     ;(IERROR,IOWAR,MAINF,FILENAME,
     ;'NUMBER OF PARAMETERS TO BE ESTIMATED (PAR FILE) IS GREATER'//
     ;' THAN THE DEFINED AT DIM FILE',NROW,0,IUPAR,2,0.07)

                IPNT_START(ICOMP)=(ICOMP-1)*IDIMWGT+1
                IPNT_END(ICOMP)=IPNT_START(ICOMP)
                IPNT_PAR((ICOMP-1)*IDIMWGT+1)=IPARDET
                WGT_TRA((ICOMP-1)*IDIMWGT+1)=1.0D0
                IF(TRAM1.EQ.0D0) TRAM1=TRAZ1

                IF (IOLGTRA.EQ.1) THEN  ! Log-estimation
                   TRAM(IPARDET)=DLOG10(TRAM1)
                   TRAC(IPARDET)=DLOG10(TRAZ1)
                ELSE
                   TRAM(IPARDET)=TRAM1
                   TRAC(IPARDET)=TRAZ1
                END IF

                WGT_UNK(IPARDET)=WEIGHT

                IF (STTRA1.LE.0D0) THEN
                  CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;      'ESTIMATED ZONAL PARAM. AND STANDARD DEVIATION IS ZERO '//
     ;      'OR NEGATIVE. IT IS SET TO 1.0',NROW,0,IUPAR,1,6.16)
                   STTRA1=1.0D0
                END IF

                STTRA(ICOMP)=STTRA1
                INDPAR(IPARDET)=0
                IF (IOLGTRA.EQ.1) INDPAR(IPARDET)=IOLGTRA
             END IF

          ENDIF 

          IF (INPWR.NE.0) WRITE(MAINF,3200)
     ;       N,ISOZ(NZ),TRAZ(NZ),IVTRA1,STTRA1,TRAM1,NFNLTRA(NZ)
     ;      ,NFTTRA(NZ),IGRP1
 3200     FORMAT (1P,2I5,G12.5,I5,2G12.5,3I5)

C------------------------- Anisotropy

          IF (ISOZ1.GT.1) THEN                               
             DO I=2,ISOZ1

C------------------------- Reads the I-th tensor component

                LEAUX=LEEL(FILENAME,IUPAR,MAINF,NROW,INPWR)          
                READ(LEAUX,1100,ERR=9000) TRAZ1,IVTRA1,STTRA1,
     ;                                 TRAM1,NFNLTRA1,NFTTRA1,IGRP1
 1100           FORMAT(10X,F10.0,I5,F10.0,F10.0,3I5)        

C------------------------- Checks transmissivity read parameters 
C------------------------- through CHEK_PAR subroutine

                CALL CHECK_PAR 
     ;(6.02     ,IERROR  ,IGRP1     ,IOINV   ,IOPBLI   ,IOTIM    
     ;,IOWAR    ,IUPAR   ,IVTRA1    ,MAINF   ,MXGRPZN  ,NFNL    
     ;,NFNLTRA1 ,NFTTRA1 ,NGROUP_ZN ,NPFNL   ,STTRA1  
     ;,'   TRANSMISSIVITY   '       ,FILENAME ,NROW)

C------------------------- Sets values of the current tensor component

                ICOMP=INTRAC(I)+NZ
                IF (IOPT_GS(IGRP1,2).EQ.0) TRAZ(ICOMP)=TRAZ1
                NFTTRA(ICOMP)=NFTTRA1
                NFNLTRA(ICOMP)=NFNLTRA1
                IVTRAGRP(ICOMP)=IGRP1
                IOPTLOG(ICOMP)=IOLGTRA
                IF (IOINV.GT.0. AND.
     ;             IVTRA1.NE.0 .AND. IOPT_GS(IGRP1,2).EQ.0) THEN
                  IPARDET=IPARDET+1
                  IF (IPARDET.GT.NPAR) CALL ERROR 
     ;(IERROR,IOWAR,MAINF,FILENAME,
     ;'NUMBER OF PARAMETERS TO BE ESTIMATED (PAR FILE) IS GREATER'//
     ;' THAN THE DEFINED AT DIM FILE',NROW,0,IUPAR,2,0.07)

                  IPNT_START(ICOMP)=(ICOMP-1)*IDIMWGT+1
                  IPNT_END(ICOMP)=IPNT_START(ICOMP)
                  IPNT_PAR((ICOMP-1)*IDIMWGT+1)=IPARDET
                  WGT_TRA((ICOMP-1)*IDIMWGT+1)=1.0D0
                  IF(TRAM1.EQ.0D0) TRAM1=TRAZ1

                  IF (IOLGTRA.EQ.1) THEN  ! Log-estimation
                     TRAM(IPARDET)=DLOG10(TRAM1)
                     TRAC(IPARDET)=DLOG10(TRAZ1)
                  ELSE
                     TRAM(IPARDET)=TRAM1
                     TRAC(IPARDET)=TRAZ1
                  END IF

                  WGT_UNK(IPARDET)=WEIGHT

                  IF (STTRA1.LE.0D0) THEN
                    CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;      'ESTIMATED ZONAL PARAM. AND STANDARD DEVIATION IS ZERO '//
     ;      'OR NEGATIVE. IT IS SET TO 1.0',NROW,0,IUPAR,1,6.16)
                     STTRA1=1.0D0
                  END IF

                  STTRA(ICOMP)=STTRA1
                  INDPAR(IPARDET)=0
                  IF (IOLGTRA.EQ.1) INDPAR(IPARDET)=IOLGTRA

C-----------------Stores parameter name. Used to write PSH and PSC
                  INAME = INAME +1
                  IF (I.EQ.2) THEN
                      WRITE(PARNAME(INAME),11) 'TY',N
                  ELSE IF(I.EQ.3) THEN
                      WRITE(PARNAME(INAME),11) 'TY',N
                  END IF

                END IF

                IF (INPWR.NE.0) WRITE(MAINF,3300)
     ;       TRAZ1,IVTRA1,STTRA1,TRAM1,NFNLTRA(NZ)
     ;      ,NFTTRA(NZ),IGRP1
 3300           FORMAT(1P,10X,G12.5,I5,2G12.5,3I5)

             ENDDO  !Next transmissivity component

C------------------------- If anisotropy is smaller than problem dimension, 
C------------------------- last component is assigned

             IF (ISOZ1.EQ.2 .AND. IDIM.EQ.3) THEN

C------------------------- Inverse problem variables are not set, because 
C------------------------- Tzz is not a new parameter (it is Tyy)

                ICOMP=INTRAC(3)+NZ
                IF (IOPT_GS(IGRP1,2).EQ.0) TRAZ(ICOMP)=TRAZ1
                NFTTRA(ICOMP)=NFTTRA1
                NFNLTRA(ICOMP)=NFNLTRA1

             ENDIF
          ENDIF   ! ISOZ1 > 1
       ENDDO      !Next transmissivity zone

       RETURN
       
 9000  CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;      'GENERIC FORTRAN ERROR WHEN READING '//
     ;      'TRANSMISSIVITY ZONE PARAM.',NROW,0,IUPAR,1,6.16)
       RETURN        
       END 
