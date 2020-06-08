      SUBROUTINE KRIGING
     ;(IDEBUG         ,IDIMCROSS     ,IO_CROSS       ,IO_SAMPLE      
     ;,IO_WEIGHTS     ,KRIGTYPE      ,MAINF          ,MAXEQ          
     ;,MAXROT         ,MAXSB         ,MXDISC_GS      ,NDATA          
     ;,NDISC          ,NESTED        ,NESTIM         ,NEXDRIFTS      
     ;,NMAXKRIG       ,NMINKRIG      ,NOCTANT        ,NRESTRI
     ;,NSUP_BL_X      ,NSUP_BL_Y     ,NSUP_BL_Z      ,NUGGET         
     ;,RADIUS         ,SEARCH_ANGLE1 ,SEARCH_ANGLE2  ,SEARCH_ANGLE3  
     ;,SEARCH_ANIS1   ,SEARCH_ANIS2  ,SKMEAN         ,XMNSUP         
     ;,XSIZSUP        ,YMNSUP        ,YSIZSUP        ,ZMNSUP         
     ;,ZSIZSUP        ,ANGLE1        ,ANGLE2         ,ANGLE3         
     ;,ANIS1          ,ANIS2         ,CLOSE          ,COVESTMEAS     
     ;,COVESTMEAS_SEC ,COVMEASMEAS   ,CROSS_COV      ,DATA           
     ;,EXTATR1        ,EXTATR2       ,EXTATR3        ,EXTATR4        
     ;,IDRIF          ,ITYPEVARIO    ,IXSBTOSR       ,IYSBTOSR       
     ;,IZSBTOSR       ,NDATA_SB      ,NSAMPLE        ,OFF_DISC       
     ;,RANGE          ,ROTMAT        ,SECEST1        ,SECEST2        
     ;,SECEST3        ,SECEST4       ,SILL_NUGG      ,TMP            
     ;,TRIM_LIMIT     ,VARKRIG       ,VEXTKRIG       ,WEIGHT         
     ;,WEIGHTS        ,XDATA         ,XESTIM         ,XKRIG          
     ;,XOFF           ,YDATA         ,YESTIM         ,YKRIG          
     ;,YOFF           ,ZDATA         ,ZESTIM         ,ZKRIG          
     ;,ZOFF           ,ESTIM_VAL     ,XDATASECU      ,YDATASECU
     ;,ZDATASECU      ,MXNPRIM_GS    ,MXMEASPP_GS    ,XVARKRIG)

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 NESTED,MAXROT,KRIGTYPE,NDATA,NEXDRIFTS,MXNPRIM_GS
     ;         ,NSUP_BL_X,NSUP_BL_Y,NSUP_BL_Z,MAXSB,NESTIM,NDISC,IDEBUG
     ;         ,MXDISC_GS,MAINF,NOCTANT,NMAXKRIG,NMINKRIG,MAXEQ
     ;         ,IO_SAMPLE,IO_CROSS,IDIMCROSS,IO_WEIGHTS,MXMEASPP_GS
     ;         ,ITYPEVARIO(NESTED),IDRIF(9),NDATA_SB(MAXSB)
     ;         ,IXSBTOSR(8*MAXSB),IYSBTOSR(8*MAXSB),IZSBTOSR(8*MAXSB)
     ;         ,NSAMPLE(MXNPRIM_GS,NESTIM)
                                                                 ! Real external
      REAL*8 NUGGET,SEARCH_ANGLE1,SEARCH_ANGLE2,SEARCH_ANGLE3
     ;      ,SEARCH_ANIS1,SEARCH_ANIS2,RADIUS,SKMEAN
     ;      ,XMNSUP,XSIZSUP,YMNSUP,YSIZSUP,ZMNSUP,ZSIZSUP
     ;      ,ROTMAT(MAXROT,3,3),ANGLE1(NESTED),ANGLE2(NESTED)
     ;      ,ANGLE3(NESTED),ANIS1(NESTED),ANIS2(NESTED),TMP(NDATA)
     ;      ,SILL_NUGG(NESTED),XDATA(NDATA),YDATA(NDATA),ZDATA(NDATA)
     ;      ,DATA(NDATA),XESTIM(NESTIM),YESTIM(NESTIM),ZESTIM(NESTIM)
     ;      ,OFF_DISC(MXDISC_GS,3,NESTIM),XOFF(NDISC),YOFF(NDISC)
     ;      ,ZOFF(NDISC),RANGE(NESTED),SECEST1(NESTIM),SECEST2(NESTIM)
     ;      ,SECEST3(NESTIM),SECEST4(NESTIM),TRIM_LIMIT(8),CLOSE(NDATA)
     ;      ,XKRIG(NMAXKRIG),YKRIG(NMAXKRIG),ZKRIG(NMAXKRIG)
     ;      ,VARKRIG(NMAXKRIG),VEXTKRIG(NMAXKRIG,4),EXTATR1(NDATA)
     ;      ,EXTATR2(NDATA),EXTATR3(NDATA),EXTATR4(NDATA)
     ;      ,COVMEASMEAS(MAXEQ*MAXEQ),COVESTMEAS(MAXEQ)
     ;      ,CROSS_COV(IDIMCROSS,NESTIM),COVESTMEAS_SEC(MAXEQ)
     ;      ,WEIGHT(MAXEQ),WEIGHTS(IDIMCROSS,NESTIM),ESTIM_VAL(NESTIM)
     ;      ,XDATASECU(NDATA),YDATASECU(NDATA),ZDATASECU(NDATA)
                                                              ! Integer internal
      INTEGER*4 NSBTOSR,NRESTRI,IESTIM,IDISC,JDISC,NEQUA,NRESTRIUNIV
     ;         ,INF_OCTANTS,NCLOSE,NACCEPT_SAMPLES,ISAMPLE,ID_SAMPLE
     ;         ,JSAMPLE,IROWCOVMEASMEAS,IROWCROSSCOV,IEXT,ISINGULAR
     ;         ,IDATA,IRESTRI
                                                                 ! Real internal
      REAL*8 XMEANKRIG,XVARKRIG,XEST,YEST,ZEST,CMAX,COVA,BLOCK_COVA
     ;      ,SCALE_OK,SCALE_UNIV,SCALE_EXDRIFT,DX,DY,DZ,DIST,ESTIMATION
     ;      ,VARESTIM,COVMAX
     ;      ,MEAN_DRIFT(9),EXTDRESTIM(4)

C_______________________ Step 1: Set up of rotation matrices considering 
C_______________________         anisotropy of search ellipsoid and all nested 
C_______________________         structures defining variogram. Checks if
C_______________________         universal kriging is done and, in this case,
C_______________________         calculates the scale factor. Calculates maximum
C_______________________         covariance.

      CALL ROTMAT_UNIVKRIG_COVMAX
     ;(COVMAX        ,KRIGTYPE      ,MAXROT        ,NESTED    
     ;,NRESTRIUNIV   ,NUGGET        ,RADIUS        ,SCALE_UNIV 
     ;,SEARCH_ANGLE1 ,SEARCH_ANGLE2 ,SEARCH_ANGLE3 ,SEARCH_ANIS1
     ;,SEARCH_ANIS2  ,ANGLE1        ,ANGLE2        ,ANGLE3
     ;,ANIS1         ,ANIS2         ,IDRIF         ,ITYPEVARIO
     ;,ROTMAT        ,SILL_NUGG)

C_______________________ Step 2: Set up of super block search grid. Defines
C_______________________         regular grid of superblocks and the offsets of
C_______________________         superblocks that should be used for a given 
C_______________________         estimation point centered at (0,0,0)

      CALL SET_SUPERBLOCK_GRID
     ;(NDATA     ,XDATA   ,YDATA     ,ZDATA     ,DATA     ,TMP     
     ;,NEXDRIFTS ,EXTATR1 ,EXTATR2   ,EXTATR3   ,NDATA_SB ,NSUP_BL_X
     ;,XMNSUP    ,XSIZSUP ,NSUP_BL_Y ,YMNSUP    ,YSIZSUP  ,NSUP_BL_Z
     ;,ZMNSUP    ,ZSIZSUP)

      CALL WHICH_SUPER_BLOCK
     ;(NSUP_BL_X ,XSIZSUP ,NSUP_BL_Y ,YSIZSUP       ,NSUP_BL_Z ,ZSIZSUP
     ;,NESTED+1  ,MAXROT  ,ROTMAT    ,RADIUS*RADIUS ,NSBTOSR   ,IXSBTOSR
     ;,IYSBTOSR  ,IZSBTOSR)

C_______________________ Step 3: Calculates number of restrictions of the 
C_______________________         extended matrix of kriging system

      NRESTRI = 1                                                ! OK
      NRESTRI = NRESTRI + NRESTRIUNIV                            ! UK
      IF (KRIGTYPE.EQ.3) NRESTRI = NRESTRI + NEXDRIFTS           ! KED
      IF (KRIGTYPE.EQ.0) NRESTRI = 0                             ! SK
      IF (KRIGTYPE.EQ.2) NRESTRI = NRESTRI + 1                   ! KLVM

C_______________________ Step 4: Big loop over estimation points. 
C_______________________         Initialization of accumulators

      XMEANKRIG=0.0D0
      XVARKRIG=0.0D0

      DO IESTIM=1,NESTIM

C_______________________ Step 4.1: Identifies coordinates of estimation point
C_______________________           and offsets of block discretization points
C_______________________           and initializes kriging system arrays

         XEST=XESTIM(IESTIM)
         YEST=YESTIM(IESTIM)
         ZEST=ZESTIM(IESTIM)

         IF (NDISC.GT.1) THEN
           DO IDISC=1,NDISC
             XOFF(IDISC)=OFF_DISC(IDISC,1,IESTIM)
             YOFF(IDISC)=OFF_DISC(IDISC,2,IESTIM)
             ZOFF(IDISC)=OFF_DISC(IDISC,3,IESTIM)
           END DO ! IDISC=1,NDISC
         END IF ! NDISC.GT.1

         CALL ZERO_ARRAY (COVMEASMEAS,MAXEQ*MAXEQ)
         CALL ZERO_ARRAY (COVESTMEAS,MAXEQ)
         CALL ZERO_ARRAY (WEIGHTS(1,IESTIM),IDIMCROSS)
         CALL ZERO_ARRAY (CLOSE,NDATA)
         CALL ZERO_ARRAY_I (NSAMPLE(1,IESTIM),NMAXKRIG)

C_______________________ Step 4.2: Calculates block covariance and the scale 
C_______________________           factor for ordinary kriging)

         CALL COVARIANCE
     ;(XOFF(1) ,YOFF(1) ,ZOFF(1) ,XOFF(1) ,YOFF(1)    ,ZOFF(1)
     ;,1       ,NESTED  ,NESTED  ,NUGGET  ,ITYPEVARIO ,SILL_NUGG
     ;,RANGE   ,1       ,MAXROT  ,ROTMAT  ,CMAX       ,COVA)

         BLOCK_COVA = COVA
         SCALE_OK = BLOCK_COVA   

         IF (NDISC.GT.1) THEN ! Block kriging 

            BLOCK_COVA = 0.0D0
            DO IDISC=1,NDISC
               DO JDISC=1,NDISC

                 CALL COVARIANCE
     ;(XOFF(IDISC) ,YOFF(IDISC) ,ZOFF(IDISC) ,XOFF(JDISC) ,YOFF(JDISC)
     ;,ZOFF(JDISC) ,1           ,NESTED      ,NESTED      ,NUGGET  
     ;,ITYPEVARIO  ,SILL_NUGG   ,RANGE       ,1           ,MAXROT  
     ;,ROTMAT      ,CMAX        ,COVA)

                 IF (IDISC.EQ.JDISC) COVA = COVA - NUGGET
                 BLOCK_COVA = BLOCK_COVA +COVA


               END DO ! JDISC=1,NDISC
            END DO ! IDISC=1,NDISC

            BLOCK_COVA = BLOCK_COVA / FLOAT(NDISC*NDISC)

         END IF ! NDISC.GT.1

C_______________________ Step 4.3: Calculates scaled mean values of polynomial 
C_______________________           drifts. These values are part of the 
C_______________________           kriging system independent term. Only if 
C_______________________           universal kriging is done
 
         IF (NRESTRIUNIV.NE.0) CALL CALC_MEAN_POLDRIFT
     ;(NDISC    ,SCALE_UNIV    ,IDRIF    ,MEAN_DRIFT    ,XOFF
     ;,YOFF     ,ZOFF)

C_______________________ Step 4.4: Assigns external drifts (if KED/KLVM is done) 
C_______________________           and calculates external drifts scaling factor
 
         IF (NEXDRIFTS.GT.0) CALL EXT_ATR_ASSIGN
     ;(COVMAX          ,IESTIM          ,MAINF          ,NEXDRIFTS    
     ;,SCALE_EXDRIFT   ,SECEST1(IESTIM) ,SECEST2(IESTIM),SECEST3(IESTIM)  
     ;,SECEST4(IESTIM) ,EXTDRESTIM      ,TRIM_LIMIT(5))

C_______________________ Step 4.5: Looks for closest measurements. Stores their 
C_______________________           ID. Saves position and meas. values at
C_______________________           auxiliar arrays. Checks if there are enough 
C_______________________           samples to estimate block.

         CALL SEARCH_SAMPLES
     ;(XEST     ,YEST      ,ZEST     ,RADIUS*RADIUS   ,NESTED+1 
     ;,MAXROT   ,ROTMAT    ,NSBTOSR  ,IXSBTOSR        ,IYSBTOSR 
     ;,IZSBTOSR ,NOCTANT   ,NDATA    ,XDATA           ,YDATA
     ;,ZDATA    ,TMP       ,NDATA_SB ,NSUP_BL_X       ,XMNSUP
     ;,XSIZSUP  ,NSUP_BL_Y ,YMNSUP   ,YSIZSUP         ,NSUP_BL_Z
     ;,ZMNSUP   ,ZSIZSUP   ,NCLOSE   ,CLOSE           ,INF_OCTANTS)

         NACCEPT_SAMPLES=0
         DO ISAMPLE=1,NCLOSE
            IF (NACCEPT_SAMPLES.LT.NMAXKRIG) THEN  ! Still not enough samples

               NACCEPT_SAMPLES = NACCEPT_SAMPLES + 1
               ID_SAMPLE = INT(CLOSE(ISAMPLE)+0.5)
                                                    ! Works ALWAYS with offsets
               XKRIG(NACCEPT_SAMPLES)= XDATA(ID_SAMPLE) - XEST
               YKRIG(NACCEPT_SAMPLES)= YDATA(ID_SAMPLE) - YEST
               ZKRIG(NACCEPT_SAMPLES)= ZDATA(ID_SAMPLE) - ZEST
               VARKRIG(NACCEPT_SAMPLES) = DATA(ID_SAMPLE)
               VEXTKRIG(NACCEPT_SAMPLES,1) = EXTATR1(ID_SAMPLE)
               VEXTKRIG(NACCEPT_SAMPLES,2) = EXTATR2(ID_SAMPLE)
               VEXTKRIG(NACCEPT_SAMPLES,3) = EXTATR3(ID_SAMPLE)
               VEXTKRIG(NACCEPT_SAMPLES,4) = EXTATR4(ID_SAMPLE)

               IF (IO_SAMPLE.NE.0) THEN
                  DO IDATA=1,NDATA

                     DX=DABS(XDATA(ID_SAMPLE)-XDATASECU(IDATA))
                     DY=DABS(YDATA(ID_SAMPLE)-YDATASECU(IDATA))
                     DZ=DABS(ZDATA(ID_SAMPLE)-ZDATASECU(IDATA))

                     IF (DX.LT.1D-6.AND.DY.LT.1D-6.AND.DZ.LT.1D-6) 
     ;                     NSAMPLE(ISAMPLE,IESTIM)=IDATA

                  END DO ! IDATA=1,NDATA
               END IF ! IO_SAMPLE.NE.0

            END IF ! NACCEPT_SAMPLES.LT.NDMAX
         END DO ! ISAMPLE=1,NCLOSE

         IF (NACCEPT_SAMPLES.LT.NMINKRIG) THEN 
           WRITE(6,1000) NACCEPT_SAMPLES,NMINKRIG,IESTIM
           WRITE(MAINF,1000) NACCEPT_SAMPLES,NMINKRIG,IESTIM
 1000      FORMAT(/,' ERROR: ENCOUNTERED A LOCATION TO BE ESTIMATED'
     ;              ' WHERE THERE ARE ONLY ',I5,' DATA. YOU ARE'
     ;              ' REQUESTING A MINIMUM OF ',I5,' .LOCATION: ',I5,
     ;              ' CRITICAL STOP')  
           STOP
         END IF ! NACCEPT_SAMPLES.LT.NMINKRIG

         IF (NACCEPT_SAMPLES.LT.NRESTRI) THEN  
           WRITE(6,1100) NACCEPT_SAMPLES,NRESTRI,IESTIM
           WRITE(MAINF,1100) NACCEPT_SAMPLES,NRESTRI,IESTIM
 1100      FORMAT(/,' ERROR: ENCOUNTERED A LOCATION TO BE ESTIMATED'
     ;              ' WHERE THERE ARE ONLY ',I5,' DATA. YOU ARE'
     ;              ' REQUESTING ',I5,' RESTRICTIONS DUE TO DRIFT'
     ;              ' TERMS.LOCATION: ',I5,' CRITICAL STOP')  
           STOP
         END IF ! NACCEPT_SAMPLES.LT.NMINKRIG

C_______________________ Step 4.6: Builds and solves kriging system (if there
C_______________________           are enough samples)
C_______________________            Calculates kriging matrix and independent
C_______________________            term (scalars in this case), solves the 
C_______________________            system, updates accumulators and goes to
C_______________________            step 4.7. Also, cross-covariance matrix is
C_______________________            saved, for a posteriori covariance matrix
C_______________________            calculation purposes

         NEQUA = NACCEPT_SAMPLES + NRESTRI

                         ! Calculates kriging matrix (measurements covariances)

         DO ISAMPLE = 1,NACCEPT_SAMPLES 

            DO JSAMPLE = 1,NACCEPT_SAMPLES 
 
               CALL COVARIANCE
     ;(XKRIG(ISAMPLE)  ,YKRIG(ISAMPLE)  ,ZKRIG(ISAMPLE) ,XKRIG(JSAMPLE)
     ;,YKRIG(JSAMPLE)  ,ZKRIG(JSAMPLE)  ,1              ,NESTED       
     ;,NESTED          ,NUGGET          ,ITYPEVARIO     ,SILL_NUGG   
     ;,RANGE           ,1               ,MAXROT         ,ROTMAT      
     ;,CMAX            ,COVA)
 
               COVMEASMEAS(NEQUA*(ISAMPLE-1)+JSAMPLE) = COVA
               COVMEASMEAS(NEQUA*(JSAMPLE-1)+ISAMPLE) = COVA

            END DO ! JSAMPLE = 1,NACCEPT_SAMPLES 
         END DO ! ISAMPLE = 1,NACCEPT_SAMPLES 
            
                                ! Calculates kriging matrix (OK unbias portion)

         IF (NEQUA.GT.NACCEPT_SAMPLES) THEN
            DO ISAMPLE=1,NACCEPT_SAMPLES
               COVMEASMEAS(NEQUA*(ISAMPLE-1)+NACCEPT_SAMPLES+1) = 
     ;                          SCALE_OK
               COVMEASMEAS(NEQUA*NACCEPT_SAMPLES+ISAMPLE) = 
     ;                          SCALE_OK
            END DO ! ISAMPLE=1,NACCEPT_SAMPLES
         END IF ! NEQUA.GT.NACCEPT_SAMPLES

                              ! Calculates independent term (measurements part)

         DO ISAMPLE=1,NACCEPT_SAMPLES
          
            IF (NDISC.EQ.1) THEN

               CALL COVARIANCE
     ;(XKRIG(ISAMPLE) ,YKRIG(ISAMPLE)  ,ZKRIG(ISAMPLE)  ,XOFF(1)  
     ;,YOFF(1)        ,ZOFF(1)         ,1               ,NESTED       
     ;,NESTED         ,NUGGET          ,ITYPEVARIO      ,SILL_NUGG   
     ;,RANGE          ,1               ,MAXROT          ,ROTMAT      
     ;,CMAX           ,COVESTMEAS(ISAMPLE))
                                               ! Stores cross-covariance matrix
               IF (IO_CROSS.NE.0) 
     ;            CROSS_COV(NSAMPLE(ISAMPLE,IESTIM),IESTIM) = 
     ;                COVESTMEAS(ISAMPLE)

             ELSE
                  
               COVESTMEAS(ISAMPLE)=0.0D0
               DO IDISC=1,NDISC
                  CALL COVARIANCE
     ;(XKRIG(ISAMPLE)  ,YKRIG(ISAMPLE)  ,ZKRIG(ISAMPLE)  ,XOFF(IDISC) 
     ;,YOFF(IDISC)     ,ZOFF(IDISC)     ,1               ,NESTED       
     ;,NESTED          ,NUGGET          ,ITYPEVARIO      ,SILL_NUGG   
     ;,RANGE           ,1               ,MAXROT          ,ROTMAT      
     ;,CMAX            ,COVA)

                  COVESTMEAS(ISAMPLE) = COVESTMEAS(ISAMPLE) + COVA
                  DX=XKRIG(1)-XOFF(IDISC)
                  DY=YKRIG(1)-YOFF(IDISC)
                  DZ=ZKRIG(1)-ZOFF(IDISC)
                  DIST= DX*DX + DY*DY + DZ*DZ  
                  IF (DIST.LE.1.0D-6) COVESTMEAS(ISAMPLE) = 
     ;                                COVESTMEAS(ISAMPLE) - NUGGET

               END DO ! IDISC=1,NDISC
               
               COVESTMEAS(ISAMPLE)=COVESTMEAS(ISAMPLE) / FLOAT(NDISC)

                                               ! Stores cross-covariance matrix
               IF (IO_CROSS.NE.0) 
     ;            CROSS_COV(NSAMPLE(ISAMPLE,IESTIM),IESTIM) = 
     ;               COVESTMEAS(ISAMPLE)

            END IF ! NDISC.EQ.1

         END DO ! ISAMPLE=1,NACCEPT_SAMPLES
            
                              ! Calculates independent term (OK unbias portion)

         IF (NEQUA.GT.NACCEPT_SAMPLES) THEN
            COVESTMEAS(NACCEPT_SAMPLES+1) = SCALE_OK
            IF (IO_CROSS.NE.0) CROSS_COV(MXMEASPP_GS+1,IESTIM)=SCALE_OK
         END IF ! NEQUA.GT.NACCEPT_SAMPLES

                                        ! Restrictions due to universal kriging

         IROWCOVMEASMEAS = NACCEPT_SAMPLES + 1
         IROWCROSSCOV = MXMEASPP_GS + 1

         IF (NRESTRIUNIV.GT.0) CALL CALDRIF
     ;(IDIMCROSS      ,IESTIM         ,IO_CROSS         ,IROWCOVMEASMEAS
     ;,IROWCROSSCOV   ,MAXEQ          ,NACCEPT_SAMPLES  ,NEQUA            
     ;,NESTIM         ,NMAXKRIG       ,SCALE_UNIV       ,COVESTMEAS 
     ;,COVMEASMEAS    ,CROSS_COV      ,IDRIF            ,MEAN_DRIFT 
     ;,XKRIG          ,YKRIG          ,ZKRIG)
                           
                                          ! Restrictions due to external drifts

         IF (KRIGTYPE.EQ.3 .OR. KRIGTYPE.EQ.2) THEN
            DO IEXT=1,NEXDRIFTS
               IROWCOVMEASMEAS = IROWCOVMEASMEAS + 1
               IROWCROSSCOV = IROWCROSSCOV + 1
               DO ISAMPLE=1,NACCEPT_SAMPLES
                  COVMEASMEAS(NEQUA*(IROWCOVMEASMEAS-1)+ISAMPLE) =
     ;               VEXTKRIG(ISAMPLE,IEXT)*SCALE_EXDRIFT
                  COVMEASMEAS(NEQUA*(ISAMPLE-1)+IROWCOVMEASMEAS) =
     ;               VEXTKRIG(ISAMPLE,IEXT)*SCALE_EXDRIFT
               END DO ! ISAMPLE=1,NACCEPT_SAMPLES
               COVESTMEAS(IROWCOVMEASMEAS) =
     ;            EXTDRESTIM(IEXT)*SCALE_EXDRIFT
                                               ! Stores cross-covariance matrix
                  IF (IO_CROSS.NE.0) 
     ;               CROSS_COV(IROWCROSSCOV,IESTIM) = 
     ;                  EXTDRESTIM(IEXT)*SCALE_EXDRIFT

            END DO ! IEXT=1,NEXDRIFTS
         END IF ! KRIGTYPE.EQ.3

                       ! Saves independent term (variance calculation purposes)

         CALL EQUAL_ARRAY(COVESTMEAS_SEC,COVESTMEAS,NEQUA)

                                                        ! Solves kriging system

         CALL KTSOL
     ;(NEQUA  ,1  ,1  ,COVMEASMEAS  ,COVESTMEAS  ,WEIGHT  ,ISINGULAR  
     ;,MAXEQ)

         IF (ISINGULAR.NE.0) THEN
            WRITE(MAINF,1200) IESTIM
            WRITE(6,1200) IESTIM
 1200       FORMAT(/,' ERROR: SINGULAR KRIGING MATRIX ESTIMATING'
     ;               ' LOCATION: ',I5,' . CRITICAL STOP')
            STOP
         END IF ! ISINGULAR.NE.0

                                ! Calculates estimation and estimation variance
 
         ESTIMATION=0.D0
         VARESTIM=BLOCK_COVA

         IF (KRIGTYPE.EQ.2) SKMEAN = EXTDRESTIM(1) ! KLVM

         DO ISAMPLE=1,NEQUA
            VARESTIM=VARESTIM-WEIGHT(ISAMPLE)*COVESTMEAS_SEC(ISAMPLE)
            IF (ISAMPLE.LE.NACCEPT_SAMPLES) THEN
               IF (KRIGTYPE.EQ.0 .OR. KRIGTYPE.EQ.2) THEN
                  ESTIMATION = ESTIMATION + 
     ;               WEIGHT(ISAMPLE)*(VARKRIG(ISAMPLE)-SKMEAN)
               ELSE
                  ESTIMATION = ESTIMATION + 
     ;               WEIGHT(ISAMPLE)*VARKRIG(ISAMPLE)
               END IF ! KRIGTYPE.EQ.0 .OR. KRIGTYPE.EQ.2
            END IF ! ISAMPLE.LE.NACCEPT_SAMPLES
         END DO ! ISAMPLE=1,NEQUA
         IF (KRIGTYPE.EQ.0 .OR. KRIGTYPE.EQ.2) 
     ;      ESTIMATION = ESTIMATION + SKMEAN
            ESTIM_VAL(IESTIM) = ESTIMATION

                                                     ! Writes DEBUG information
         IF (IDEBUG.NE.0) CALL WRITE_DEBUG_INFO
     ;(IESTIM            ,WEIGHT(NACCEPT_SAMPLES+1)   ,NACCEPT_SAMPLES
     ;,XEST              ,YEST                        ,ZEST    
     ;,NSAMPLE(1,IESTIM) ,VARKRIG                     ,WEIGHT
     ;,XKRIG             ,YKRIG                       ,ZKRIG)

C_______________________ Step 4.7: Saves kriging weights

         IF (IO_WEIGHTS.NE.0) THEN
            DO ISAMPLE=1,NACCEPT_SAMPLES
               WEIGHTS(NSAMPLE(ISAMPLE,IESTIM),IESTIM) = WEIGHT(ISAMPLE)
            END DO ! ISAMPLE=1,NACCEPT_SAMPLES
            IF (NRESTRI.GT.0) THEN
               IRESTRI=1
               DO ISAMPLE=NACCEPT_SAMPLES+1,NEQUA
                  WEIGHTS(MXMEASPP_GS+IRESTRI,IESTIM) = WEIGHT(ISAMPLE)
                  IRESTRI = IRESTRI + 1
               END DO
            END IF
         END IF ! IO_WEIGHT.NE.0

C_______________________ Step 4.8: End of big loop over estimation points. 
C_______________________           Updates acumulators and goes to next one

         XMEANKRIG = XMEANKRIG + ESTIMATION
         XVARKRIG = XVARKRIG + VARESTIM
        
      END DO !IESTIM=1,NESTIM

C_______________________ Step 5: Calculates global mean and variance

      XMEANKRIG = XMEANKRIG / FLOAT(NESTIM)
      XVARKRIG = XVARKRIG / FLOAT(NESTIM)
      IF (IDEBUG.NE.0) WRITE(666,1300) XMEANKRIG,XVARKRIG
 1300 FORMAT(/,' KRIGING MEAN:     ',E10.4,/
     ;        ,' KRIGING VARIANCE: ',E10.4)

      RETURN
      END

