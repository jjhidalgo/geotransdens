      SUBROUTINE COMP_OBS_GEO
     ;     (CCAL     ,ACTH     ,CFPAREL  ,DERC     ,DVOBS    ,IDIMDERC
     ;     ,IFLAGS   ,INEW     ,INEWT    ,INORPAR  ,IODEVICE ,IOINV
     ;     ,IOPLC    ,IOSMTP   ,IPBTP    ,LXPAREL  ,NDEVGEO  ,NDEVS
     ;     ,NFLAGS   ,NOOBSIT  ,NPAR     ,NPAREL   ,NPBMX    ,NPBTP 
     ;     ,NTYPAR   ,NUMEL    ,NUMNP    ,NUMTIT   ,NUMTOBS  ,NZPAR
     ;     ,PARZ     ,TABSOLUT,TIT      ,VJAC      ,VOBSC)
***********************************************************************
*     PURPOSE
*     
*     Computes values which correspond to the gephysical observations
*     
*     DESCRIPTION
*     
*     The subroutine is entered once for each simulation time. Initially,
*     it is looked into, for each device, whether or not it is necessary to
*     compute a value which corresponds to an observation.
*     If this is the case, gimli is called to obtain the value of
*     the observation (RESCALIT) and its dereivative with respect to 
*     concentration (DRESDC).
*     The observation is stored in VOBSC and the derivative in DVOBS.      
*     
*     WARNING: It won't work if more than one transport problem.
*     
*     EXTERNAL VARIABLES: ARRAYS
*     
*     CCAL                   Computed concentration at every node                  
*     DERC                   Nodal concentration derivatives with respect to       
*     estimated parameters.                                 
*     DVOBS                                                                        
*     INDEXNOD               Index relating nodes                                  
*     IODEVICE               Column 1: Data type                                   
*     Column 2: Status for calc. of obs.                    
*     Column 3: Method of spat. integr.                     
*     Column 4: Method of temp. integr.                     
*     Column 5: Number of integr. time                      
*     NOOBSIT                Observation number to which an integration time       
*     belongs to                                            
*     RESCAL                 Computed resistivities at every node                          
*     TIT                    Integration time                                      
*     TOBS                   Time of observation                                   
*     VJAC                   Jacobian matrix                                       
*     VOBSC                  Value of simulated value corresponding to observation 
*     
*     INTERNAL VARIABLES: ARRAYS
*     
*     
*     EXTERNAL VARIABLES: SCALARS
*     
*     INEW                                                                         
*     IOLD                                                                         
*     NDEVS                                                                        
*     NPAR                   Total number of parameters to be estimated            
*     NUMNP                  Number of nodes                                       
*     NUMTIT                 Total number of integration times                     
*     NUMTNOD                Total number of nodes used for calculating obs.       
*     NUMTOBS                Total number of observations                          
*     TABSOLUT               Current absolut computation time                      
*     TINC                   Current time increment                                
*     
*     INTERNAL VARIABLES: SCALARS
*     
*     ND                                                                           
*     NOF                                                                          
*     
*     FUNCTIONS AND SUBROUTINES REFERENCED
*     
*     COMP_OBS_AUX                                                                 
*     
*     HISTORY
*     
*     JHG       5-2020     First coding.
*     
***********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER*4::NDEVGEO

      DIMENSION TIT(NUMTIT),IODEVICE(NDEVS+1,10)
     ;     ,CCAL(NUMNP,NPBTP)
     ;     ,DERC(NUMNP,NPAR,IDIMDERC,NPBTP)
     ;     ,IFLAGS(NFLAGS) 
     ;     ,VOBSC(NUMTOBS+NDEVS),VJAC(NUMTOBS,NPAR)
     ;     ,NOOBSIT(NUMTIT),DVOBS(NPAR),CFPAREL(NUMEL,NPAREL)
     ;     ,PARZ(NZPAR),LXPAREL(NUMEL,NPAREL,NPBMX), ACTH(NUMEL)
     ;     ,INORPAR(NTYPAR)

C-------------------------Internal variables
      INTEGER*4::IOCALGEO
      REAL*8::RESCAL(NDEVGEO),DRESDC(NUMNP,NDEVGEO),POROSITY(NUMEL)

C-------------------------First executable statment.

      IF (IFLAGS(3).EQ.1) CALL IO_SUB('COMP_OBS_GEO',0)

      IOCALGEO = 1
C-------------------------Checks if any device has geophysical observations
C-------------------------at the current solution time.
      DO ND=1,NDEVS
         DVOBS = 0D0
         TNX=TIT(IODEVICE(ND,2)) ! Next integration time

C------------------------If not geophysical data or not used, go to next device.
         IF (IODEVICE(ND,1).NE.6   .OR. 
     &        IODEVICE(ND,2).LT.0 .OR.
     &        IODEVICE(ND,10).EQ.0) THEN
            
            CYCLE

         END IF

C------------------------Checks device problem. 
         IF (IOSMTP.EQ.0) THEN
            IPROB=1             !always problem 1 if simultaneous
         ELSE
            IPROB=IODEVICE(ND,9) ! Flow/tpt. problem 
         ENDIF
         
         IF (IOSMTP.NE.0 .OR. IPBTP.EQ.IODEVICE(ND,9)) THEN
            TNX=TIT(IODEVICE(ND,2)) ! Next int. time
         ELSE
            CYCLE 
         END IF

         IF ((DABS(TABSOLUT-TNX).LE.1.E-15*(TABSOLUT+TNX)/2D0 .OR.
     ;        TABSOLUT-TNX.GT.1.0E-15).AND. IODEVICE(ND,2).GT.0) THEN

C------------------------Geophysical observations are computed
C------------------------the first time a device with geophysical
C------------------------observations is found.
            IF (IOCALGEO.EQ.1) THEN

                POROSITY = CFPAREL(1:NUMEL,7)
     &                     *PARZ(INORPAR(15)+LXPAREL(1:NUMEL,7,IPROB))
     &                     *ACTH(1:NUMEL)

                CALL   CALC_GIMLI(CCAL(1,IPROB)
     &                           ,DERC(1,1,INEWT,IPROB)        ,DRESDP
     &                           ,'HYDROMESH.vtk'    ,IFLAGS   ,NDEVGEO
     &                           ,NFLAGS             ,NPAR     ,NUMELEM
     &                           ,NUMNP              ,POR      ,RESCAL
     &                           ,TABSOLUT)

                
               IOCALGEO = 0     !Do not calc. geo. obs. again this time step.             
               !RESCAL = 1D0     !for verification
               !DRESDC = 1D0
            END IF
            
            NO=NOOBSIT(IODEVICE(ND,2)) ! Observation number


            VOBSC(NO)=RESCAL(ND)

            IF (IOINV.GT.0) THEN
               DO NP=1,NPAR
                  VJAC(NO,NP)=DOT_PRODUCT(DRESDC(:,ND)
     ;                                   ,DERC(:,NP,INEW,IPBT))
               ENDDO
            END IF
            
            IODEVICE(ND,2)=IODEVICE(ND,2)+1 ! Next integr. time

*_____________Number of integr. times not exceeded

            IF(NOOBSIT(IODEVICE(ND,2)).LT.IODEVICE(ND+1,8) .AND.        &
     &           IODEVICE(ND,2).LE.NUMTIT) THEN
               TNX=TIT(IODEVICE(ND,2)) ! Next int. time
            ELSE
               IODEVICE(ND,2)=-IODEVICE(ND,2) ! "Finished" index
               write(777,*) 'DEVICE FINISHED'
               write(777,*)
            ENDIF

         ENDIF                  !(DABS(TABSOLUT-TNX).LE.1.E-15*(TABSOLUT+TNX)/2D0 ...
      END DO                    !NDEVS

      IF(IFLAGS(3).EQ.1) CALL IO_SUB('COMP_OBS_GEO',1)

      RETURN
      END SUBROUTINE COMP_OBS_GEO
