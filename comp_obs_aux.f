      SUBROUTINE COMP_OBS_AUX
     ;(INEW     ,IOINV    ,IOLD     ,IOPL     ,ND       ,NDEVS    
     ;,NPAR     ,NUMNP    ,NUMTIT   ,NUMTNOD  ,NUMTOBS  ,TABSOLUT 
     ;,TINC     ,DERVAR   ,DVOBS    ,INDEXNOD ,IODEVICE ,NOOBSIT  
     ;,TIT      ,TOBS     ,VARNEW   ,VAROLD   ,VJAC     ,VOBSC    
     ;,WTOBSN   ,WTOBST   ,TIME     ,NINT)

*********************************************************
*
* PURPOSE
*
*
* DESCRIPTION
*
* EXTERNAL VARIABLES: ARRAYS
*
*  DERVAR                 Derivatives of state var. with respect to par.        
*  DVOBS                  Derivatives of observations with respect to param.
*  INDEXNOD               Index relating nodes                                  
*  IODEVICE               Column 1: Data type                                   
*                         Column 2: Status for calc. of obs.                    
*                         Column 3: Method of spat. integr.                     
*                         Column 4: Method of temp. integr.                     
*                         Column 5: Number of integr. time                      
*  NOOBSIT                Observation number to which an integration time       
*                         belongs to                                            
*  TIT                    Integration time                                      
*  TOBS                   Time of observation                                   
*  VARNEW                 Value at present time (T)                             
*  VAROLD                 Value at prev. time (T-DT)                            
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
*  INEW                   Index related to DERH, DERC management
*  IOLD                   Index related to DERH, DERC management
*  ND                     Numer of current device
*  NDEVS                  Number of devices
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
*  I                                                                            
*  NFN                                                                          
*  NLN                                                                          
*  NNOD                                                                         
*  NO                                                                           
*  NOF                                                                          
*  NOL                                                                          
*  NP                                                                           
*  NUMFIN                                                                       
*  TE                                                                           
*  TNX                                                                          
*  TOB                                                                          
*  VOBS0                                                                        
*  VOBS1                                                                        
*  VOBS2                                                                        
*  WT                                                                           
*  WT1                                                                          
*  WT2                                                                          
*  WT3                                                                          
*  WT4                                                                          
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ZERO_ARRAY                                                                   
*
* HISTORY: CK    11-1999 First coding
*          AAR   07-2001 Revision and header inclusion
*
********************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION WTOBSN(NUMTNOD),TIT(NUMTIT),INDEXNOD(NUMTNOD),
     ;NOOBSIT(NUMTIT),VARNEW(NUMNP),VAROLD(NUMNP),VJAC(NUMTOBS,NPAR),
     ;IODEVICE(NDEVS+1,10),TOBS(NUMTOBS,2),VOBSC(NUMTOBS+NDEVS),
     ;WTOBST(NUMTOBS),DERVAR(NUMNP,NPAR,2),DVOBS(NPAR,3),TIME(NINT)

C______________________________ Identifies and initialises some variables

      NFN=IODEVICE(ND,7)                        ! First node constituting device
      NLN=IODEVICE(ND+1,7)-1                     ! Last node constituting device

      VOBS0=0D0
      VOBS2=0D0
      IF (IOINV.GT.0) CALL ZERO_ARRAY (DVOBS,3*NPAR)

!-----Summation over nodes that belong to the device.
!-----Computes the observation values at current (VOBS2) and previous (VBOS0) solution times.
!-----This is done for all devices and used for PLT when IOPL=2.
!-----Then VOBS0 and VBOS1 are used to compute observations at the observation times.
!-----Becasue it is written in VOBSC in a position after NUMTOBS.
      
      DO I=NFN,NLN                                ! Spacial summation over nodes

        NNOD=INDEXNOD(I)                                  ! Index relating nodes
        VOBS2=VOBS2+WTOBSN(I)*VARNEW(NNOD)                           ! Value (T)

        IF (IODEVICE(ND,2).GT.0) THEN
          VOBS0=VOBS0+WTOBSN(I)*VAROLD(NNOD)                        ! Value (T-DT)
          IF (IOINV.GT.0) THEN
            DO NP=1,NPAR
              DVOBS(NP,1)=DVOBS(NP,1)+WTOBSN(I)*DERVAR(NNOD,NP,INEW) ! Der. (T-DT)
              DVOBS(NP,2)=DVOBS(NP,2)+WTOBSN(I)*DERVAR(NNOD,NP,IOLD) ! Der. (T)
            ENDDO
          END IF
        END IF
      ENDDO

C______________________________ Assigns value at TABSOLUT
!-----VOBSC contains the observations in ALL TIMES. Used to write PLT and maybe more things.
!-----IOPL=1 write observation times; IOPL=2 write all solution times.
      IF (IOPL.EQ.2) VOBSC(NUMTOBS+ND)=VOBS2

C______________________________ Subroutine ends for actual device if all its 
C______________________________ measurement have been used
!-----IODEVICE(ND,10) is the iniital (read in OBS file) value of IODEVIDE(ND,2).
!----- Initially IODEVICE(ND,2) is the starting integration time (position in TIT, NOOBSIT).
!-----If 0, it means that the device is not used for calibration.
!-----During the algorithm  IODEVICE(ND,2) is incremented by 1 to go to the next integration time.
!-----When the number of integration times is exceeded IODEVICE(ND,2) changes to a negative value.
!-----IODEVICE(ND,2) negative makes the algorithm stop.
!-----The original value of  IODEVICE(ND,2) is recovered in SIM_JAC from  IODEVICE(ND,10).
!-----I do not know why an auxiiar variable was not used.
        IF (IODEVICE(ND,2).LT.0 .OR. IODEVICE(ND,10).EQ.0) RETURN

*_______________________Intervals not def. by simulation
!----- Time integration point (IODEVICE(ND,4) = 0) or average in interval (IODEVICE(ND,4) = 2).
      IF(IODEVICE(ND,4).EQ.0 .OR. IODEVICE(ND,4).EQ.2) THEN

        TNX=TIT(IODEVICE(ND,2))                          ! Next integration time

C______________________________ The condition in the while looop
C______________________________ avoids problems with rounding errors.IODEVICE(ND,4).EQ.0
C______________________________ OLD version: (TNX.LE.TABSOLUT .AND. IODEVICE(ND,2).GT.0)
!------Checks if the solution time correspond to observation time.

        DO WHILE ((DABS(TABSOLUT-TNX).LE.1.E-15*(TABSOLUT+TNX)/2D0 .OR.
     ;      TABSOLUT-TNX.GT.1.0E-15).AND. IODEVICE(ND,2).GT.0)

          IF(TABSOLUT.EQ.TIME(1)) THEN         ! T=0 => All weight on obs. at T=0
            WT=1D0
          ELSE
            WT=(TNX-TABSOLUT+TINC)/TINC
          ENDIF

*_______________________Integration time values of obs. and derivative 
*_______________________Linear interpolation


          VOBS1=VOBS2*WT+VOBS0*(1-WT)
          IF (IOINV.GT.0) THEN
            DO NP=1,NPAR
              DVOBS(NP,3)=DVOBS(NP,2)*WT+DVOBS(NP,1)*(1-WT)
            ENDDO
          END IF

*_______________________Checking for passing of interval ends (weights reduced)

          NO=NOOBSIT(IODEVICE(ND,2))                        ! Observation number

          IF(IODEVICE(ND,4).EQ.2 .AND. (TNX.EQ.TOBS(NO,1) .OR.
     ;      TNX.EQ.TOBS(NO,2))) THEN
            VOBS1=VOBS1/2.0                           ! Weight reduction at ends
            IF (IOINV.GT.0) THEN
              DO NP=1,NPAR
                DVOBS(NP,3)=DVOBS(NP,3)/2.0           ! Weight reduction at ends
              ENDDO
            END IF
          ENDIF

          VOBSC(NO)=VOBSC(NO)+VOBS1*WTOBST(NO)

          IF (IOINV.GT.0) THEN
            DO NP=1,NPAR
              VJAC(NO,NP)=VJAC(NO,NP)+DVOBS(NP,3)*WTOBST(NO)
            ENDDO
          END IF

          IODEVICE(ND,2)=IODEVICE(ND,2)+1                    ! Next integr. time

*_____________Number of integr. times not exceeded

          IF(NOOBSIT(IODEVICE(ND,2)).LT.IODEVICE(ND+1,8) .AND.
     ;       IODEVICE(ND,2).LE.NUMTIT) THEN
            TNX=TIT(IODEVICE(ND,2))                             ! Next int. time
          ELSE
            IODEVICE(ND,2)=-IODEVICE(ND,2)                    ! "Finished" index
          ENDIF
        ENDDO !while TNX.LE.TABSOLUT

*_______________________Intervals def. by simulation

      ELSE

        NOF=IODEVICE(ND,8)
        NOL=IODEVICE(ND+1,8)-1
        NUMFIN=0                             ! Init. of counter of finished obs.

        DO NO=NOF,NOL                           ! Check for contr. to more obs.?
          TOB=TOBS(NO,1)
          TE=TOBS(NO,2)

*_______________________Necessary to calculate contribution?

          IF(TABSOLUT.GE.TOB .AND. (TABSOLUT-TINC).LT.TE) THEN

*_______________________Observation not ended

            IF(TABSOLUT.LT.TE) THEN

*_______________________Observation initiated or continued

              IF((TABSOLUT-TINC).GE.TOB) THEN

                WT=WTOBST(NO)*TINC/2.0
                VOBSC(NO)=VOBSC(NO)+WT*(VOBS0+VOBS2)
                IF (IOINV.GT.0) THEN
                  DO NP=1,NPAR
                    VJAC(NO,NP)=VJAC(NO,NP)+WT*(DVOBS(NP,1)+DVOBS(NP,2))
                  ENDDO
                END IF

*_______________________Observation initiated

              ELSE IF(TABSOLUT.GE.TOB .AND. (TABSOLUT-TINC).LT.TOB) THEN

                WT1=WTOBST(NO)*(TABSOLUT-TOB)*0.5
                WT2=(TOB-TABSOLUT+TINC)/TINC
                WT3=WT1*(1+WT2)
                WT4=WT1*(1-WT2)
                VOBSC(NO)=WT3*VOBS2+WT4*VOBS0

                IF (IOINV.GT.0) THEN
                  DO NP=1,NPAR
                    VJAC(NO,NP)=WT3*DVOBS(NP,2)+WT4*DVOBS(NP,1)
                  ENDDO
                END IF

              ENDIF

*_______________________Observation ended (and maybe also initiated)

*_______________________Observation initiated and ended

            ELSE IF((TABSOLUT-TINC).LT.TOB) THEN

              WT1=WTOBST(NO)*(TE-TOB)
              WT2=(TE+TOB-TABSOLUT-TABSOLUT)/(TINC+TINC)
              WT3=-WT1*WT2
              WT4=WT1+WT1*WT2
              VOBSC(NO)=WT3*VOBS0+WT4*VOBS2
              IF (IOINV.GT.0) THEN
                DO NP=1,NPAR
                  VJAC(NO,NP)=WT3*DVOBS(NP,1)+WT4*DVOBS(NP,2)
                ENDDO
              END IF
              NUMFIN=NUMFIN+1               ! Update of counter of finished obs.

*_______________________Observation ended

            ELSE

              WT1=WTOBST(NO)*(TE-TABSOLUT+TINC)
              WT2=(TE-TABSOLUT+TINC)/(TINC+TINC)
              WT4=WT1*WT2
              WT3=WT1-WT4
              VOBSC(NO)=VOBSC(NO)+WT3*VOBS0+WT4*VOBS2
              IF (IOINV.GT.0) THEN
                DO NP=1,NPAR
                  VJAC(NO,NP)=VJAC(NO,NP)
     ;                       +WT3*DVOBS(NP,1)+WT4*DVOBS(NP,2)
                ENDDO
              END IF
              NUMFIN=NUMFIN+1               ! Update of counter of finished obs.

            ENDIF

*_______________________Observation already ended

          ELSE IF((TABSOLUT-TINC).GE.TE) THEN
            NUMFIN=NUMFIN+1                 ! Update of counter of finished obs.
          ENDIF

        ENDDO

*_____________________Checks whether all obs. (and thereby device) are finished

        IF(NUMFIN.EQ.(NOL-NOF+1)) IODEVICE(ND,2)=-1

      ENDIF

      RETURN
      END
