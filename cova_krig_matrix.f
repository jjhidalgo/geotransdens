      SUBROUTINE COVA_KRIG_MATRIX
     ;(IDIMCOV_GS    ,IDIMCROSS    ,IOCOVPAR    ,IODEBUG   ,MAXROT        
     ;,MXDISC_GS     ,NDISC        ,NESTED      ,NESTIM    ,NRESTRI       
     ;,NUGGET        ,NZEROS       ,COVA_KRIG   ,COVPAR    ,CROSSCOV      
     ;,ITYPEVARIO    ,MAT_DIFF     ,OFFSET      ,RANGE     ,ROTMAT    
     ;,SILL_NUGG     ,WEIGHTS      ,XESTIM      ,YESTIM    ,ZESTIM)

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                             ! Integer external
      INTEGER*4 NESTIM,IDIMCROSS,NDISC,MXDISC_GS,NRESTRI
     ;         ,NESTED,MAXROT,NZEROS,IDIMCOV_GS,IODEBUG,IOCOVPAR
     ;         ,ITYPEVARIO(NESTED)
                                                                ! Real external
      REAL*8 NUGGET,CROSSCOV(IDIMCROSS,NESTIM),WEIGHTS(IDIMCROSS,NESTIM)
     ;      ,XESTIM(NESTIM),YESTIM(NESTIM),ZESTIM(NESTIM)
     ;      ,OFFSET(MXDISC_GS,3,NESTIM),SILL_NUGG(NESTED),RANGE(NESTED)
     ;      ,ROTMAT(MAXROT,3,3),COVPAR(IDIMCOV_GS)
     ;      ,COVA_KRIG(NESTIM*(NESTIM+1)/2)
     ;      ,MAT_DIFF(NESTIM*(NESTIM+1)/2)
                                                             ! Integer internal
      INTEGER*4 I,J,K,IDISC,JDISC,ILAST,IPOSCOVA_KRIG,IROWOLD,IPOSCOVPAR
                                                                ! Real internal
      REAL*8 AUXSUM,CROSS,LAMBDA,X1,Y1,Z1,X2,Y2,Z2,COVAPRIORIJ
     ;      ,XGAUSSI,YGAUSSI,ZGAUSSI,XGAUSSJ,YGAUSSJ,ZGAUSSJ,CMAX,COVA

C_______________________ Step 1: Multiplies extended cross-covariance matrix
C_______________________         times kriging weights matrix

      CALL ZERO_ARRAY(MAT_DIFF,NESTIM*(NESTIM+1)/2)

      DO I=1,NESTIM
         DO J=1,I
            AUXSUM=0.0D0
            IPOSCOVA_KRIG = I * (I+1) / 2 - (I - J)
            DO K=1,IDIMCROSS
                  CROSS = CROSSCOV(K,I)
                  LAMBDA = WEIGHTS(K,J)
                  AUXSUM = AUXSUM + CROSS * LAMBDA
            END DO ! K=1,IDIMCROSS
            MAT_DIFF(IPOSCOVA_KRIG) = MAT_DIFF(IPOSCOVA_KRIG) + AUXSUM
 
         END DO ! J=1,NESTIM
      END DO ! I=1,NESTIM

C_______________________ Step 2: Loop over estimated points/blocks (ROWS), 
C_______________________         identifying coordinates

      IROWOLD=0
      DO I=1,NESTIM
  
         X1=XESTIM(I)
         Y1=YESTIM(I)
         Z1=ZESTIM(I)

C_______________________ Step 2.1: For each row, runs over all columns
C_______________________           identifying coordinates

         DO J=1,I

            X2=XESTIM(J)
            Y2=YESTIM(J)
            Z2=ZESTIM(J)

C_______________________ Step 2.1.2: Calculates a priori covariance (I,J)

            COVAPRIORIJ = 0.0D0

            DO JDISC=1,NDISC

               IF (NDISC.GT.1) THEN
                  XGAUSSJ = X2 + OFFSET(JDISC,1,J)
                  YGAUSSJ = Y2 + OFFSET(JDISC,2,J)
                  ZGAUSSJ = Z2 + OFFSET(JDISC,3,J)
               ELSE
                  XGAUSSJ = X2 
                  YGAUSSJ = Y2 
                  ZGAUSSJ = Z2 
               END IF ! NDISC.GT.1

               DO IDISC=1,NDISC

                  IF (NDISC.GT.1) THEN
                     XGAUSSI = X1 + OFFSET(IDISC,1,I)
                     YGAUSSI = Y1 + OFFSET(IDISC,2,I)
                     ZGAUSSI = Z1 + OFFSET(IDISC,3,I)
                  ELSE
                     XGAUSSI = X1
                     YGAUSSI = Y1
                     ZGAUSSI = Z1
                  END IF ! NDISC.GT.1

                  CALL COVARIANCE
     ;(XGAUSSI     ,YGAUSSI     ,ZGAUSSI      ,XGAUSSJ  ,YGAUSSJ
     ;,ZGAUSSJ     ,1           ,NESTED       ,NESTED   ,NUGGET  
     ;,ITYPEVARIO  ,SILL_NUGG   ,RANGE        ,1        ,MAXROT  
     ;,ROTMAT      ,CMAX        ,COVA)

                  COVAPRIORIJ = COVAPRIORIJ + COVA

               END DO ! IDISC=1,NDISC

            END DO ! JDISC=1,NDISC

            COVAPRIORIJ = COVAPRIORIJ / FLOAT(NDISC*NDISC)

C_______________________ Step 1.1.3: Calculates positions at COVPAR
C_______________________             (global covariance matrix of parameters)
C_______________________             and COVA_KRIG (kriging covariance matrix)

            IPOSCOVA_KRIG=I*(I+1)/2-(I-J)
            IF (I.GT.IROWOLD) THEN
               IPOSCOVPAR=I*(I+1)/2-(I-J)+I*NZEROS
               ILAST=IPOSCOVPAR
               IROWOLD=I
            ELSE
               IPOSCOVPAR = ILAST + 1
               ILAST = ILAST + 1
            END IF ! I.GT.IROWOLD

C_______________________ Step 1.1.4: Saves covariance

            IF (IOCOVPAR.NE.0) 
     ;         COVPAR(IPOSCOVPAR) = COVAPRIORIJ-MAT_DIFF(IPOSCOVA_KRIG)
            COVA_KRIG(IPOSCOVA_KRIG) = 
     ;         COVAPRIORIJ - MAT_DIFF(IPOSCOVA_KRIG)

         END DO ! J=1,I

      END DO ! I=1,NESTIM+NRESTRI

C_______________________ Step 2: Echoes kriging covariance matrix

      IF (IODEBUG.NE.0) THEN
         WRITE(666,1000)
 1000    FORMAT(//,' KRIGING COVARIANCE MATRIX',/
     ;            ,' ======= ========== ======',/)

         WRITE(666,1100) COVA_KRIG
 1100    FORMAT(7E10.4)
      END IF ! IODEBUG.NE.0

      RETURN
      END 


