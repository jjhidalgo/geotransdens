      SUBROUTINE CONDITIONAL_SIMULATION
     ;(NPOINTS       ,NDATA       ,MAINF
     ;,MXNPRIM_GS  ,IDIMCROSS_GS ,ISEED_SC    ,IACTSIMUL
     ;,IDEBUG        ,KTYPE_GR    ,MXKRIG_GR    ,MXROT_GR    ,MXSBX_GR
     ;,MXSBY_GR      ,MXSBZ_GR    ,MXDISC_GS    ,MXSB_GR     ,NDISC_GR
     ;,IDIMIVARIO_GS ,NEXDR_GR    ,NMXP_GR      ,NDMIN_GR    ,NOCT_GR
     ;,IDIMVAR_GS    ,MXCLOSE_GS  ,MXMEASPP_GS  ,MXZONPP_GS  ,MAXSAM_GR
     ;,NUMSB_GS      ,ICROSSCOV_GS,IFLAG_SIMUL  ,IPOLDRIFT_GR,IVARIO_GR   
     ;,ISUPBL_GS     ,SIMVALUE    ,X_SIM        ,Y_SIM       ,Z_SIM,MEAS  
     ;,TRIM_GR       ,X_DATA      ,Y_DATA       ,Z_DATA
     ;,RND_PATH      ,SEQPOSDAT   ,SEQDATA      ,EXDRIFT_POINT1
     ;,EXDRIFT_POINT2,EXDRIFT_POINT3,EXDRIFT_POINT4          ,VARIO_GR    
     ;,SEARCH_GR     ,VSTATS_GR   ,SUPBL_GR     ,KRISOL_GS  ,CLOSESAM_GS 
     ;,KRISYS_GR     ,CROSSCOV_GS ,POSDIS_GR    ,ROTMAT_GR   ,KRIGAUX_GS  
     ;,POSDISAUX_GS  ,IDIMDATASC_GR,EXDRIFT_MEAS1,EXDRIFT_MEAS2
     ;,EXDRIFT_MEAS3 ,EXDRIFT_MEAS4,SEQDATA_SEC ,SEQPOSDAT_SEC)

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 NPOINTS,NDATA,MAINF,MXNPRIM_GS
     ;         ,IDIMCROSS_GS,ISEED_SC,IACTSIMUL,IDEBUG,KTYPE_GR
     ;         ,MXKRIG_GR,MXROT_GR,MXSBX_GR,MXSBY_GR,MXSBZ_GR,MXDISC_GS
     ;         ,MXSB_GR,NDISC_GR,IDIMIVARIO_GS,NEXDR_GR,NMXP_GR,NDMIN_GR
     ;         ,NOCT_GR,IDIMVAR_GS,MXCLOSE_GS,MXZONPP_GS,MAXSAM_GR
     ;         ,IDIMDATASC_GR,MXMEASPP_GS
     ;         ,NUMSB_GS(MXSB_GR),IFLAG_SIMUL(NPOINTS)
     ;         ,ICROSSCOV_GS(MXNPRIM_GS,NPOINTS),IPOLDRIFT_GR(9)
     ;         ,IVARIO_GR(IDIMIVARIO_GS,2),ISUPBL_GS(MXSB_GR*8,3)

                                                                 ! Real external
      REAL*8 SIMVALUE(NPOINTS),X_SIM(NPOINTS),Y_SIM(NPOINTS)
     ;      ,Z_SIM(NPOINTS),MEAS(NDATA),TRIM_GR(8),X_DATA(NDATA)
     ;      ,Y_DATA(NDATA),Z_DATA(NDATA)
     ;      ,RND_PATH(NPOINTS,2),SEQPOSDAT(IDIMDATASC_GR,3)
     ;      ,SEQDATA(IDIMDATASC_GR,5),EXDRIFT_MEAS1(NDATA)
     ;      ,EXDRIFT_MEAS2(NDATA),EXDRIFT_MEAS3(NDATA)
     ;      ,EXDRIFT_MEAS4(NDATA)
     ;      ,EXDRIFT_POINT1(NPOINTS),EXDRIFT_POINT2(NPOINTS)
     ;      ,EXDRIFT_POINT3(NPOINTS),EXDRIFT_POINT4(NPOINTS)
     ;      ,VARIO_GR(IDIMIVARIO_GS,8),SEARCH_GR(11)
     ;      ,VSTATS_GR(IDIMVAR_GS,4),SUPBL_GR(6),KRISOL_GS(MXKRIG_GR,3)
     ;      ,CLOSESAM_GS(MXCLOSE_GS,2),KRISYS_GR(MXKRIG_GR*MXKRIG_GR)
     ;      ,CROSSCOV_GS(IDIMCROSS_GS,MXZONPP_GS,2)
     ;      ,POSDIS_GR(NDISC_GR,3,NPOINTS),ROTMAT_GR(3,3,MXROT_GR)
     ;      ,KRIGAUX_GS(MAXSAM_GR,8),POSDISAUX_GS(NDISC_GR,3)
     ;      ,SEQPOSDAT_SEC(IDIMDATASC_GR,3),SEQDATA_SEC(IDIMDATASC_GR,5)
                                                              ! Integer internal
      INTEGER*4 IPOINT,IDUMMY,IDPOINT,NUMDATA,ISEED(13)
                                                                 ! Real internal
      REAL*8 R_DUMMY,ACORNI,CMEAN,GAUPROB,CVARIANCE

C_______________________ Step 1: Initializes values of simulated points FLAG

      DO IPOINT=1,NPOINTS
        IFLAG_SIMUL(IPOINT)=-1
      END DO

C_______________________ Step 2: Work out a random path (see footnote 1)

C_______________________ Step 2.1: Initialisation of auxiliar array for random
C_______________________           number generation. To avoid repeated 
C_______________________           sequences, input seed is multiplied per 
C_______________________           number of actual simulation

      CALL ZERO_ARRAY_I(ISEED,13)

*** OJO2 PARA VERIFICACION DE DOS SIMULACIONES CONDICIONADAS IDENTICAS
*      ISEED(1)=ISEED_SC
      ISEED(1)=ISEED_SC*IACTSIMUL !orig estaba así en el de Andrés.
*** FIN DEL OJO2

      DO IDUMMY=1,5
        R_DUMMY=ACORNI(ISEED)  ! ACORNI will make ISEED to change 5 times
      END DO

C_______________________ Step 2.2: Calculates random numbers

      DO IPOINT=1,NPOINTS
        RND_PATH(IPOINT,2)=ACORNI(ISEED)                  ! Random numbers (0,1)
        RND_PATH(IPOINT,1)=IPOINT                 ! Will contain the random path
      END DO      

C_______________________ Step 2.3: Sorts RND_PATH(i,2) in increasing order, 
C_______________________           accordingly to RND_PATH(i,1) values. At the
C_______________________           end of this step, RND_PATH(i,1) contains the
C_______________________           random path

      CALL SORTEM 
     ;(1,NPOINTS,RND_PATH(1,2),1,RND_PATH(1,1),R_DUMMY,R_DUMMY,R_DUMMY
     ;,R_DUMMY,R_DUMMY,R_DUMMY)

C_______________________ Step 3: Loads measurements in sequential updatable data 
C_______________________         arrays

      CALL ZERO_ARRAY (SEQPOSDAT_SEC,IDIMDATASC_GR*3)   ! Initialization
      CALL ZERO_ARRAY (SEQDATA_SEC,IDIMDATASC_GR*5)
      CALL ZERO_ARRAY (SEQPOSDAT,IDIMDATASC_GR*3)
      CALL ZERO_ARRAY (SEQDATA,IDIMDATASC_GR*5)

      CALL EQUAL_ARRAY(SEQPOSDAT_SEC(1,1),X_DATA,NDATA) ! Coord. of measurements
      CALL EQUAL_ARRAY(SEQPOSDAT_SEC(1,2),Y_DATA,NDATA)   
      CALL EQUAL_ARRAY(SEQPOSDAT_SEC(1,3),Z_DATA,NDATA)  

      CALL EQUAL_ARRAY(SEQDATA_SEC(1,1),MEAS,NDATA)  ! Data
     
      CALL EQUAL_ARRAY(SEQDATA_SEC(1,2),EXDRIFT_MEAS1,NDATA)  ! Drifts
      CALL EQUAL_ARRAY(SEQDATA_SEC(1,3),EXDRIFT_MEAS2,NDATA)
      CALL EQUAL_ARRAY(SEQDATA_SEC(1,4),EXDRIFT_MEAS3,NDATA)
      CALL EQUAL_ARRAY(SEQDATA_SEC(1,5),EXDRIFT_MEAS4,NDATA)

      NUMDATA=NDATA    ! Number of data for first point to be simulated. 
                       ! In step 4, this variable is incremented (+1) sequentially.

C_______________________ Step 4: Loop over points to be simulated

      DO IPOINT=1,NPOINTS

C_______________________ Step 4.1: Selects point to be simulated from RND_PATH
C_______________________           and verifies that this point has not been
C_______________________           already simulated/assigned 

        IDPOINT=INT(RND_PATH(IPOINT,1)+0.5)

        IF (IDEBUG.NE.0) THEN
           IF (IPOINT.EQ.1) WRITE (696,2000) IACTSIMUL
 2000      FORMAT(' VALORES EN LA SIMULACION CONDICIONADA:',I5,/
     ;           ,' ======= == == ========== =============',//
     ;           ,'  PTO       GAUPROB         MEDIA       SIGMA      '
     ;            ' VALOR',/,'  ===       =======         =====       '
     ;            '=====       =====')
        END IF
     

        IF (IFLAG_SIMUL(IDPOINT).EQ.-1) THEN

C_______________________ Step 4.2: Calculates the conditional mean and st. dev. 
C_______________________           Recovers data to be updated, as they change
C_______________________           at KRIGING routine

           CALL EQUAL_ARRAY (SEQPOSDAT(1,1),SEQPOSDAT_SEC(1,1),NUMDATA)
           CALL EQUAL_ARRAY (SEQPOSDAT(1,2),SEQPOSDAT_SEC(1,2),NUMDATA)
           CALL EQUAL_ARRAY (SEQPOSDAT(1,3),SEQPOSDAT_SEC(1,3),NUMDATA)
           CALL EQUAL_ARRAY (SEQDATA(1,1),SEQDATA_SEC(1,1),NUMDATA)
           CALL EQUAL_ARRAY (SEQDATA(1,2),SEQDATA_SEC(1,2),NUMDATA)
           CALL EQUAL_ARRAY (SEQDATA(1,3),SEQDATA_SEC(1,3),NUMDATA)
           CALL EQUAL_ARRAY (SEQDATA(1,4),SEQDATA_SEC(1,4),NUMDATA)
           CALL EQUAL_ARRAY (SEQDATA(1,5),SEQDATA_SEC(1,5),NUMDATA)

           IDUMMY=1  ! Not used, only to avoid passing a lot of unuseful argts
           R_DUMMY=0.D0

           CALL KRIGING
     ;(IDEBUG                ,IDIMCROSS_GS        ,0       
     ;,0                     ,0                   ,KTYPE_GR      
     ;,MAINF                 ,MXKRIG_GR           ,MXROT_GR       
     ;,MXSB_GR               ,MXDISC_GS           ,NUMDATA          
     ;,NDISC_GR              ,IVARIO_GR(1,1)      ,1         
     ;,NEXDR_GR              ,NMXP_GR             ,NDMIN_GR      
     ;,NOCT_GR               ,IDUMMY              ,MXSBX_GR      
     ;,MXSBY_GR              ,MXSBZ_GR            ,VARIO_GR(1,1)         
     ;,SEARCH_GR(1)          ,SEARCH_GR(9)        ,SEARCH_GR(10)    
     ;,SEARCH_GR(11)         ,SEARCH_GR(3)        ,SEARCH_GR(4)     
     ;,VSTATS_GR(1,2)        ,SUPBL_GR(1)         ,SUPBL_GR(4)      
     ;,SUPBL_GR(2)           ,SUPBL_GR(5)         ,SUPBL_GR(3)      
     ;,SUPBL_GR(6)           ,VARIO_GR(1,4)       ,VARIO_GR(1,5)    
     ;,VARIO_GR(1,6)         ,VARIO_GR(1,7)       ,VARIO_GR(1,8)    
     ;,CLOSESAM_GS(1,1)      ,KRISOL_GS(1,1)      ,KRISOL_GS(1,2)       
     ;,KRISYS_GR             ,CROSSCOV_GS(1,IDPOINT,1)  ,SEQDATA(1,1)         
     ;,SEQDATA(1,2)          ,SEQDATA(1,3)        ,SEQDATA(1,4)
     ;,SEQDATA(1,5)          ,IPOLDRIFT_GR(1)     ,IVARIO_GR(1,2)      
     ;,ISUPBL_GS(1,1)        ,ISUPBL_GS(1,2)      ,ISUPBL_GS(1,3)      
     ;,NUMSB_GS              ,ICROSSCOV_GS(1,IDPOINT),POSDIS_GR(1,1,1)
     ;,VARIO_GR(1,2)         ,ROTMAT_GR           ,EXDRIFT_POINT1  
     ;,EXDRIFT_POINT2        ,EXDRIFT_POINT3      ,EXDRIFT_POINT4
     ;,VARIO_GR(1,3)         ,CLOSESAM_GS(1,2)    ,TRIM_GR(1)          
     ;,KRIGAUX_GS(1,4)       ,KRIGAUX_GS(1,5)     ,KRISOL_GS(1,3)      
     ;,CROSSCOV_GS(1,IDPOINT,2),SEQPOSDAT(1,1)    ,X_SIM(IDPOINT)   
     ;,KRIGAUX_GS(1,1)       ,POSDISAUX_GS(1,1)   ,SEQPOSDAT(1,2)
     ;,Y_SIM(IDPOINT)        ,KRIGAUX_GS(1,2)     ,POSDISAUX_GS(1,2)  
     ;,SEQPOSDAT(1,3)        ,Z_SIM(IDPOINT)      ,KRIGAUX_GS(1,3)     
     ;,POSDISAUX_GS(1,3)     ,CMEAN               ,SEQPOSDAT_SEC(1,1)
     ;,SEQPOSDAT_SEC(1,2)    ,SEQPOSDAT_SEC(1,3)  ,MXNPRIM_GS
     ;,MXMEASPP_GS           ,CVARIANCE)

C_______________________ Step 4.3: Draw a N(0,1) random nuber (factor of ST DEV)
 
          CALL NORMAL(ISEED(1),1.0D0,GAUPROB,0.0D0)
     

C_______________________ Step 4.4: Changes ISEED for next point

          DO IDUMMY=1,5
            R_DUMMY=ACORNI(ISEED)
          END DO

C_______________________ Step 4.5: Assigns final value to current point. If 
C_______________________ we are conditioning zones to measurements, we assign
C_______________________ directly to the value of the measurement that falls 
C_______________________ within current zone

          SIMVALUE(IDPOINT)=GAUPROB*DSQRT(CVARIANCE)+CMEAN

*** tal vez esto vaya fuera; si se queda, deberia asignar los pesos de puntos piloto
*** de este pixel a 0.0. Esa es la unica manera de que se respeten a capon las medidas
*** Tal cual esta sin asignar los pesos de pp a 0, solo se respetan las medidas para la SC inicial

*          IF (ICONDZONES.EQ.1) THEN   ! Conditioning of zones to meas.
*             IF (IZONMEAS_GS(IDPOINT).NE.0) THEN   ! Zone contains a measurement
*                SIMVALUE(IDPOINT)=MEAS(IZONMEAS_GS(IDPOINT))
*             END IF
*          END IF

*** hasta aqui

C_______________________ Step 4.6: Adds simulated point to the set of data for
C_______________________           further simulated points

          NUMDATA=NUMDATA+1
          SEQPOSDAT_SEC(NUMDATA,1)=X_SIM(IDPOINT)
          SEQPOSDAT_SEC(NUMDATA,2)=Y_SIM(IDPOINT)
          SEQPOSDAT_SEC(NUMDATA,3)=Z_SIM(IDPOINT)
          SEQDATA_SEC(NUMDATA,1)=SIMVALUE(IDPOINT)          
          SEQDATA_SEC(NUMDATA,2)=EXDRIFT_POINT1(IDPOINT)
          SEQDATA_SEC(NUMDATA,3)=EXDRIFT_POINT2(IDPOINT)
          SEQDATA_SEC(NUMDATA,4)=EXDRIFT_POINT3(IDPOINT)
          SEQDATA_SEC(NUMDATA,5)=EXDRIFT_POINT4(IDPOINT)
          IFLAG_SIMUL(IDPOINT)=1
          IF (IDEBUG.NE.0) WRITE(696,2100) IDPOINT,GAUPROB,CMEAN
     ;                   ,DSQRT(CVARIANCE),SIMVALUE(IDPOINT)

        END IF ! IFLAG_SIMUL(IDPOINT).EQ.-1

      END DO ! IPOINT=1,NPOINTS

 2100 FORMAT(I5,4(3X,E10.4))
      RETURN
      END

********************************************************************************
****************** GAUSSIAN RANDOM NUMBER GENERATOR ****************************
********************************************************************************

      SUBROUTINE NORMAL(NSEED,STD,XGEN,XMEAN)

********************************************************************************
*
* PURPOSE Generates a random number XGEN from a Gaussian distribution with mean
*         XMEAN and standard deviation STD with integer seed NSEED. Should not
*         repeat numbers if NSEED varies from realization to realization.
*
* EXTERNAL VARIABLES: SCALARS
*
*  NSEED              Integer seed to generate random number
*  STD                Standard deviation of Gaussian distribution
*  XGEN               Generated random number
*  XMEAN              Mean of Gaussian distribution
*
* INTERNAL VARIABLES: SCALARS
*
*  I                  Dummy counter
*
********************************************************************************

      IMPLICIT NONE
      REAL*8 XGEN,XMEAN,STD 
      REAL*4 RAN            
      INTEGER*4 NSEED,I
      XGEN=0.D0
      DO I=1,12
         XGEN=XGEN+DBLE(RAN(NSEED))
      END DO
      XGEN=XMEAN+STD*(XGEN-6.D0)
      RETURN
      END
