       SUBROUTINE READ_PAR
     ;(ERNUM    ,IDIMWGT  ,IERROR     ,INPWR      ,IOINV      ,IOLGVAR
     ;,IOPBLI   ,IOTIM    ,IPARDET    ,IOWAR      ,ISTART     ,IUPAR    
     ;,MAINF    ,MXGRPZN  ,NFNL       ,NGROUP_ZN  ,NPAR       ,NPFNL      
     ;,NZVAR    ,VAR      ,WEIGHT     ,FILENAME   ,INDPAR     ,IOPTLOG
     ;,IOPT_GS  ,IPNT_END ,IPNT_PAR   ,IPNT_START ,IVVARGRP   ,NFNLVAR
     ;,NFTVAR   ,STVAR    ,VARC       ,VARM       ,VARZ       ,WGT_UNK    
     ;,WGT_VAR  ,NROW)
                          
*****************************************************************************
*
* PURPOSE
*      Reads and checks zonal values of a given parameter type except 
*      transmissivities,due to the anisotropy case
*
* DESCRIPTION
*
* 1) Writes main header
* 2) Loop over zones of actual parameter type
*   2.1) Reads all information
*   2.2) Checks sequential number of zones and all supplied information
*   2.3) Zonal parameter, non linear and time function are always assigned
*   2.4) If parameter is estimated deterministically, assigns pointers to 
*        IPNT_PAR and some useful data. If parameter is estim. geost. or 
*        interpolated assignations will be done elsewhere
*   2.5) Writes parameters once assigned
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
*  IOPTLOG                Part of IVPAR (column 4), containing log-estimation 
*                         option for all zones related to this parameter type
*  IOPT_GS                Vector contaning geost. options of all groups of zones
*  IPNT_END               Part of column 2 of array IVPAR, related to actual 
*                         parametyer type. Contains the last useful position 
*                         at array IPNT_PAR to be used on the parameterization
*                         of actual zonal parameter
*  IPNT_PAR               Array contaning pointers to arrays DLT_PAR and WGT_PAR
*                         for defining linear combinations of unknowns
*  IPNT_START             Part of column 1 of array IVPAR, related to actual 
*                         parametyer type. Contains the first useful position 
*                         at array IPNT_PAR to be used on the parameterization
*                         of actual zonal parameter
*  ITYPEVAR               Part of column 2 of array IVPAR, related to actual 
*                         parametyer type. Contains index of parameter type
*  IVVARGRP               Part of column 4 of array IVPAR, related to actual 
*                         parametyer type. Contains the group of zones to which
*                         actual zone belongs to
*  NFNLVAR                Array contaning non linear functions indexes of actual 
*                         parameter type
*  NFTVAR                 Array containing time functions indexes of actual 
*                         parameter type
*  STVAR                  Array containing standard deviations of actual 
*                         parameter type
*  VARC                   Calcululated values of unknown parameters
*  VARM                   Prior information zonal values of actual parameter type
*  VARZ                   Zonal values of actual parameter type
*  WGT_UNK                Array containing WEIGHTS (for objective function
*                         calculation) of problems unknowns
*  WGT_VAR                Weights defining the linear combination of unknowns
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  ERNUM                  Error message code
*  IDIMWGT                Used to dimension WGT_PAR
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IOINV                  Inverse problem option                                
*  IOLGVAR                Log-estimation option of actual param. type
*  IOPBLI                 If zero, linear  problem, otherwise nonlinear.
*  IOTIM                  Problem regime (0, SS, 1, transient with given 
*                         initial cond., 2 transient with SS initial cond.)
*  IPARDET                Counter of deterministically estimated parameters
*  ISTART                 Starting point at some arrays
*  IUPAR                  Unit number of file PAR                               
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  MXGRPZN                Used to dimension IOPT_GS
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NGROUP_ZN              Number of groups of zones
*  NPAR                   Total number of parameters to be estimated            
*  NPFNL                  Counts the number of nonlinear functions read
*  NZVAR                  Number of zones of actual type of param.
*  VAR                    String containing name of actual parameter type
*  WEIGHT                 Parameters objective function weight related to actual
*                         type of parameter
*
* INTERNAL VARIABLES: SCALARS
*
*  IGRP1                  Read group of zones
*  IVVAR1                 Read estimation index
*  LEAUX                  Auxiliar string where the last read line of the       
*                         current input file is stored                          
*  N                      Dummy counter of zones of actual param. type
*  NFNLVAR1               Read nonlinear function number
*  NFTVAR1                Read time function number
*  NROW                   Current record number                                 
*  NZ                     Zone number as read
*  STVAR1                 Read standard deviation of a parameter               `
*  VARM1                  Read prior information value of the current zone
*  VARZ1                  Read initial value of the current zone
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  CHECK_PAR              Checks zonal data
*  ERROR                  Writes the current error message and error number     
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
*
* HISTORY
*
*     SCR      5-1997     First version
*     AMS      1-1998     Revision. Addition of header
*     AAR      7-2003     Revision.
*
*****************************************************************************

*______________________ 0) Declaration of variables

       IMPLICIT NONE
                                                            ! Integer external
       INTEGER*4 NZVAR,NPAR,IDIMWGT,INPWR,MAINF,IUPAR,IERROR,IOWAR,IOINV
     ;          ,IOPBLI,IOTIM,MXGRPZN,NFNL,NPFNL,IOLGVAR,IPARDET
     ;          ,ISTART,NGROUP_ZN 
     ;          ,NFNLVAR(NZVAR),INDPAR(NPAR),NFTVAR(NZVAR)
     ;          ,IVVARGRP(NZVAR),IPNT_START(NZVAR),IPNT_END(NZVAR)
     ;          ,IOPT_GS(MXGRPZN,20)
     ;          ,IPNT_PAR(NZVAR*IDIMWGT),IOPTLOG(NZVAR)
                                                               ! Real external
       REAL*8 VARZ(NZVAR),STVAR(NZVAR),VARM(NPAR),WGT_UNK(NPAR)
     ;       ,VARC(NPAR),WGT_VAR(NZVAR*IDIMWGT),WEIGHT
       REAL*4 ERNUM
                                                            ! Integer internal
       INTEGER*4 N,NROW,NZ,IVVAR1,NFNLVAR1,NFTVAR1,IGRP1
                                                               ! Real internal
       REAL*8 VARZ1,STVAR1,VARM1
                                                                  ! Characters
       CHARACTER*20 FILENAME(18),LEEL*100,LEAUX*100,VAR        

*_______________________ 1) Writes main header

       IF (INPWR.NE.0) WRITE(MAINF,3000) VAR
 3000  FORMAT(//10X,A20,' ZONES',/,
     ; 10X,'--------------------------',/,
     ; ' ZONE   COMPUTED    IV   ST.DESV.    MEASURED  NFNL  NFT IGRP')

*_______________________ 2) Loop over zones of actual parameter type

       DO N=1,NZVAR
                   
*_______________________ 2.1) Reads all information
          IGRP1 = 0
          LEAUX=LEEL(FILENAME,IUPAR,MAINF,NROW,INPWR)          
          READ(LEAUX,1000,ERR=9000) NZ,VARZ1,IVVAR1,STVAR1,
     ;                             VARM1,NFNLVAR1,NFTVAR1,IGRP1
 1000     FORMAT(I5,F10.0,I5,F10.0,F10.0,3I5) !READ    
 3100     FORMAT(1P,I5,G12.5,I5,2G12.5,3I5)   !WRITE
         
*_______________________ 2.2) Checks sequential number of zones and all 
*_______________________      supplied information

          IF (NZ.NE.N)          
     ;       CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;       VAR//' ZONE NUMBER IS OUT OF SEQUENCE',
     ;       NROW,1,IUPAR,1,ERNUM)  
    
          CALL CHECK_PAR 
     ;(ERNUM    ,IERROR  ,IGRP1     ,IOINV   ,IOPBLI   ,IOTIM    
     ;,IOWAR    ,IUPAR   ,IVVAR1    ,MAINF   ,MXGRPZN  ,NFNL    
     ;,NFNLVAR1 ,NFTVAR1 ,NGROUP_ZN ,NPFNL   ,STVAR1   ,VAR      
     ;,FILENAME  ,NROW)

*_______________________ 2.3) Non linear and time function are always assigned
*_______________________      Zonal parameter is assigned if estimation is
*_______________________      deterministical or if it is not estimated

         IF (IOPT_GS(IGRP1,2).EQ.0) VARZ(NZ)=VARZ1
         NFTVAR(NZ)=NFTVAR1
         NFNLVAR(NZ)=NFNLVAR1

*_______________________ 2.4) If parameter is estimated deterministically, 
*_______________________      assigns pointers to IPNT_PAR and some useful data
*_______________________      If parameter is estim. geost. or interpolated,
*_______________________      assignations will be done elsewhere

         IVVARGRP(NZ)=IGRP1
         IOPTLOG(NZ)=IOLGVAR

                        ! Inverse problem solved and deterministical estimation
         IF (IOINV.GT.0. 
     ;       AND. IVVAR1.NE.0 .AND. IOPT_GS(IGRP1,2).EQ.0) THEN   

           IPARDET=IPARDET+1
           IF (IPARDET.GT.NPAR)  CALL ERROR 
     ;(IERROR,IOWAR,MAINF,FILENAME,
     ;'NUMBER OF PARAMETERS TO BE ESTIMATED (PAR FILE) IS GREATER'//
     ;' THAN THE DEFINED AT DIM FILE',NROW,0,IUPAR,2,ERNUM+0.07)


           IPNT_START(NZ)=ISTART+(NZ-1)*IDIMWGT
           IPNT_END(NZ)=IPNT_START(NZ)

           IPNT_PAR((NZ-1)*IDIMWGT+1)=IPARDET
           WGT_VAR((NZ-1)*IDIMWGT+1)=1.0D0

           IF (VARM1.EQ.0D0) VARM1=VARZ1
           IF (IOLGVAR.EQ.1) THEN  ! Log-estimation
              VARM(IPARDET)=DLOG10(VARM1)
              VARC(IPARDET)=DLOG10(VARZ1)
           ELSE
              VARM(IPARDET)=VARM1
              VARC(IPARDET)=VARZ1
           END IF
           WGT_UNK(IPARDET)=WEIGHT
           
                               ! Checks STVAR value. Echoes a warning if is zero

           IF (STVAR1.LE.0.0D0) THEN

             CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;       'PARAMETER ESTIMATED DETERMINISTICALLY AND STANDARD'//
     ;       ' DEVIATION IS ZERO OR NEGATIVE. IT IS SET TO ONE',
     ;       NROW,6,IUPAR,1,ERNUM+0.01)
             
             STVAR1=1.0D0

           END IF ! STVAR1.EQ.0.0D0

           STVAR(NZ)=STVAR1
           INDPAR(IPARDET)=0
           IF (IOLGVAR.EQ.1) INDPAR(IPARDET)=IOLGVAR

         END IF

*_______________________ 2.5) Writes parameters once assigned and checked

          IF (INPWR.NE.0) WRITE(MAINF,3100)
     ;       N,VARZ(NZ),IVVAR1,STVAR1,VARM1,NFNLVAR(NZ),NFTVAR(NZ)
     ;      ,IVVARGRP(NZ)

       ENDDO ! N=1,NZVAR

       RETURN
       
 9000  CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;      'GENERIC FORTRAN ERROR WHEN READING '//
     ;      VAR//' ZONE PARAM.',NROW,0,IUPAR,1,ERNUM+0.07)
        
       RETURN
       END            
