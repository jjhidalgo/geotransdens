      SUBROUTINE COMP_OBS
     ;(IAOBS    ,INEW     ,IOINV    ,IOLD     ,IOPLC    ,IOPLH
     ;,NDEVS    ,NPAR     ,NPARF    ,NPBFL    ,NPBTP    ,NUMNP
     ;,NUMTIT   ,NUMTNOD  ,NUMTOBS  ,TABSOLUT ,TINC     ,CCAL
     ;,CCALAN   ,DERC     ,DERH     ,DVOBS    ,HCAL     ,HCALAN
     ;,INDEXNOD ,IODEVICE ,NOOBSIT  ,TIT      ,TOBS     ,VJAC
     ;,VOBSC    ,WTOBSN   ,WTOBST   ,INEWT    ,IOLDT    ,IOSMFL
     ;,IOSMTP   ,IPBFL    ,IPBTP    ,TIME     ,NINT     ,IDIMDERH
     ;,IDIMDERC ,IFLAGS   ,NFLAGS)
***********************************************************************
* PURPOSE
*
* Computes values which correspond to the observations
*
* DESCRIPTION
*
* The subroutine is entered once for each simulation time. Initially,
* it is looked into, for each device, whether or not it is necessary to
* compute a value which corresponds to an observation - or a part of it
* in the case of observations over a time interval. If this is the
* case, values of nodal heads and concentrations are weighted spatially
* and temporally to obtain the value corresponding to an observation or
* a part of it.
* Similarly, the "observation" Jacobian (the derivative of the
* observations with respect to the parameters to be estimated) is
* determined by applying the same weights to the "nodal" Jacobians (the
* derivative of nodal heads and concentration with respect to the
* parameters to be estimated).
*
* EXTERNAL VARIABLES: ARRAYS
*
*  CCAL                   Computed concentration at every node                  
*  CCALAN                 Computed concentrations in the previous time step.    
*  DERC                   Nodal concentration derivatives with respect to       
*                         estimated parameters.                                 
*  DERH                   Nodal head derivatives with respect to estimated      
*                         flow parameters.                                      
*  DVOBS                                                                        
*  HCAL                   Computed heads at every node                          
*  HCALAN                 Head level at previous time                           
*  INDEXNOD               Index relating nodes                                  
*  IODEVICE               Column 1: Data type                                   
*                         Column 2: Status for calc. of obs.                    
*                         Column 3: Method of spat. integr.                     
*                         Column 4: Method of temp. integr.                     
*                         Column 5: Number of integr. time                      
*  NOOBSIT                Observation number to which an integration time       
*                         belongs to                                            
*  RESCAL                 Computed resistivities at every node                          
*  TIT                    Integration time                                      
*  TOBS                   Time of observation                                   
*  VJAC                   Jacobian matrix                                       
*  VOBSC                  Value of simulated value corresponding to observation 
*  WTOBSN                 Weight for node used to calculate observation         
*  WTOBST                 Weight for integration time                           
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  INEW                                                                         
*  IOLD                                                                         
*  NDEVS                                                                        
*  NPAR                   Total number of parameters to be estimated            
*  NUMNP                  Number of nodes                                       
*  NUMTIT                 Total number of integration times                     
*  NUMTNOD                Total number of nodes used for calculating obs.       
*  NUMTOBS                Total number of observations                          
*  TABSOLUT               Current absolut computation time                      
*  TINC                   Current time increment                                
*
* INTERNAL VARIABLES: SCALARS
*
*  ND                                                                           
*  NOF                                                                          
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  COMP_OBS_AUX                                                                 
*
* HISTORY
*
*     CK      11-1999     First coding
*
***********************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION TIT(NUMTIT),IODEVICE(NDEVS+1,9),TOBS(NUMTOBS,2)
     ;         ,HCAL(NUMNP,NPBFL),HCALAN(NUMNP,NPBFL)
     ;         ,CCAL(NUMNP,NPBTP),CCALAN(NUMNP,NPBTP)
     ;         ,DERH(NUMNP,NPARF,IDIMDERH,NPBFL)
     ;         ,DERC(NUMNP,NPAR,IDIMDERC,NPBTP)
     ;         ,IFLAGS(NFLAGS) 
     ;         ,VOBSC(NUMTOBS+NDEVS),VJAC(NUMTOBS,NPAR),TIME(NINT)


C-------------------------First executable statment.

      IF (IFLAGS(3).EQ.1) CALL IO_SUB('COMP_OBS',0)

      
      DO ND=1,NDEVS             ! Checks for each device if necessary to calc.
        NOF=IODEVICE(ND,8)                        ! First obs. of current device NOT USED!!!

C______________________________ Datatype: Head level. Computed head or VJAC or
C______________________________ both must be calculated

          IF (IODEVICE(ND,1).EQ.1) THEN
             IF (IOSMFL.EQ.0) THEN
                IPROB=1
              ELSE
                IPROB=IODEVICE(ND,9)  ! Flow/tpt. problem 
             ENDIF

C------------------------- Objective function must be computed only 
C------------------------- when the actual flow problem number coincides with 
C------------------------- the problem number of the device OR when the flow
C------------------------- problems are simultaneous, because in this case the
C------------------------- variable ISOLEQ(.,3) losses its significance

             IF (IOSMFL.NE.0 .OR. IPBFL.EQ.IODEVICE(ND,9)) 
     ;          CALL COMP_OBS_AUX
     ;(INEW     ,IOINV    ,IOLD     ,IOPLH    ,ND       ,NDEVS
     ;,NPARF    ,NUMNP    ,NUMTIT   ,NUMTNOD  ,NUMTOBS  ,TABSOLUT
     ;,TINC,DERH(1,1,1,IPROB),DVOBS ,INDEXNOD ,IODEVICE ,NOOBSIT
     ;,TIT      ,TOBS     ,HCAL(1,IPROB),HCALAN(1,IPROB),VJAC,VOBSC
     ;,WTOBSN   ,WTOBST    ,TIME     ,NINT)

*             write(82,*) ' !!',IOSMFL,IPBFL,IODEVICE(ND,9),nd
*             write(82,*)  ' obsc',vobsc
          ENDIF

C______________________________ Datatype: Concentration. Computed head or VJAC 
C______________________________ or both must be calculated

          IF (IODEVICE(ND,1).EQ.2) THEN
             IF (IOSMTP.EQ.0) THEN
                IPROB=1
              ELSE
                IPROB=IODEVICE(ND,9)  ! Flow/tpt. problem 
             ENDIF

C------------------------- Objective function must be computed only 
C------------------------- when the actual flow problem number coincides with 
C------------------------- the problem number of the device OR when the flow
C------------------------- problems are simultaneous, because in this case the
C------------------------- variable ISOLEQ(.,4) losses its significance

             IF (IOSMTP.NE.0 .OR. IPBTP.EQ.IODEVICE(ND,9)) THEN

                 CALL COMP_OBS_AUX
     ;(INEWT    ,IOINV   ,IOLDT     ,IOPLC    ,ND       ,NDEVS
     ;,NPAR     ,NUMNP   ,NUMTIT    ,NUMTNOD  ,NUMTOBS  ,TABSOLUT
     ;,TINC,DERC(1,1,1,IPROB),DVOBS ,INDEXNOD ,IODEVICE ,NOOBSIT
     ;,TIT      ,TOBS    ,CCAL(1,IPROB)       ,CCALAN(1,IPROB)
     ;,VJAC     ,VOBSC   ,WTOBSN    ,WTOBST   ,TIME     ,NINT)

             END IF

          ELSE IF(IODEVICE(ND,1).EQ.3) THEN         ! Datatype: Water Contents
            CONTINUE
          ELSE IF(IODEVICE(ND,1).EQ.4) THEN             ! Datatype: Water Flux
            CONTINUE
          ELSE IF(IODEVICE(ND,1).EQ.5) THEN            ! Datatype: Solute Flux
            CONTINUE
          ELSE IF(IODEVICE(ND,1).EQ.6) THEN            ! Datatype: Geophysical
            CONTINUE !Computed in comp_obs_geo
          ENDIF                 ! Data type

      ENDDO                                                        ! Next device

      IF (IOPLH.EQ.2.OR.IOPLC.EQ.2) 
     ;   WRITE(IAOBS) TABSOLUT,(VOBSC(J),J=NUMTOBS+1,NUMTOBS+NDEVS)

      IF(IFLAGS(3).EQ.1) CALL IO_SUB('SIM_JAC',1)

      RETURN
      END
