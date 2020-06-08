      SUBROUTINE COMVEL
     &          (AREA     ,BUOYANCY ,COORD    ,DBUOYANCY
     &          ,DPARELDC ,DPARELDH ,DVDC     ,DVDH
     &          ,GP_COORD ,GRADLOC  ,GRAVEL   ,GRDFF
     &          ,HAUX1    ,IDIMQ    ,IODENS   ,IODIM
     &          ,IOFLLI   ,ISOZ     ,KXX      ,LDIM
     &          ,IFLAGS   ,LINMET   ,LMXNDL   ,LNNDEL
     &          ,LTYPE    ,LXPAREL  ,MAXPG    ,NPAREL
     &          ,NPPEL    ,NTYPEL   ,NUMEL    ,NUMNP
     &          ,NZTRA    ,PAREL    ,POINTWEIGHT
     &          ,NFLAGS   ,QXYZ     ,VD       ,XNORVD)
********************************************************************************
*
* PURPOSE
*
*    Computes Darcy's velocity, its norm and other related values
*
* DESCRIPTION
*
*    Computes Darcy's velocity, its norm and other related values
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COORD                  Nodal coordinates                                     
*  GRAVEL                 Projection of gravity at each element
*  GRDFF                  Array containing the product between interpolation    
*                         functions integrals and interp. functions gradient    
*  HCALIT                   head.
*  IFLAGS                 Array with different writing options. Used mainly for
*                         debugging.                                            
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LDIM                   Vector containing the dimension of each element       
*  LNNDEL                 Number of nodes at every element                      
*  QXYZ                   Products between the different components of          
*                         Darcy's velocity divided by its norm                  
*  VD                     Darcy's velocity                                      
*  XNORVD                 Euclidean norm of Darcy's velocity                    
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMQ                  Used to dimension array QXYZ                          
*  IODIM                  Maximum dimension of any element included             
*                         in the problem                                        
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFLAGS                 Maximum number of allowed flags                       
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*
* INTERNAL VARIABLES: SCALARS
*
*  L                      Counter index of elements                                                      
*  NNUD                   Number of nodes of the current element                
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  CALHEAD                Transforms pressure to head level
*
* HISTORY
*
*     AMS        1988     First coding
*     JCA      5-1998     Includes pressures and ETRA array
*     AMS     12-1998     Revision and small modifications
*     AMS      8-2000     Correction of array IND
*
********************************************************************************

       IMPLICIT NONE
C---------------------------  EXTERNAL VARIABLES: SCALARS

      INTEGER*4::IDIMQ,IOCALCDEVF,IOCALCDEVT,IODENS,IODIM,K
     &          ,LMXNDL,MAXPG,NPAREL,NPPEL,NTYPEL,NUMEL,NFLAGS
     &          ,NUMNP,NZTRA,IOCALCDEV,IOFLLI,gravel(numel,iodim)


C---------------------------  EXTERNAL VARIABLES: ARRAYS

      REAL*8::AREA(NUMEL),COORD(NUMNP,3)
     &       ,DPARELDC(NPPEL,NUMEL),DPARELDH(NPPEL,NUMEL)
     &       ,DVDC(LMXNDL,IODIM,NUMEL),DVDH(LMXNDL,IODIM,NUMEL)
     &       ,GRDFF(IODIM, LMXNDL,NUMEL),HAUX1(NUMNP) 
     &       ,GRADLOC(IODIM,LMXNDL,MAXPG),PAREL(NUMEL,NPPEL)
     &       ,POINTWEIGHT(MAXPG,NTYPEL),QXYZ(IDIMQ, NUMEL)
     &       ,VD(IODIM, NUMEL),XNORVD(NUMEL),GP_COORD(6,8,3)
     &       ,BUOYANCY(IODIM,LMXNDL,NUMEL)
     &       ,DBUOYANCY(IODIM,LMXNDL*LMXNDL,NUMEL)
    
      INTEGER*4::ISOZ(NZTRA),IFLAGS(NFLAGS)
     &          ,KXX(LMXNDL, NUMEL),LDIM(NUMEL)
     &          ,LNNDEL(NUMEL),LTYPE(NUMEL),LXPAREL(NUMEL,NPAREL)
     &          ,LINMET(3,2)

C---------------------------  INTERNAL VARIABLES: SCALARS

      INTEGER*4::I,L,LD, LANI, NZONE,ISZ,ISMAX,IS,IELEMTYP,NNUD

      REAL*8::AREAL, XNOR

C--------------------------- INTERNAL VARIABLES: ARRAYS

      REAL*8::TRACT(9),DVDH_L(LMXNDL,IODIM),DVDC_L(LMXNDL,IODIM)
     &       ,VEL(IODIM)

C--------------------------- Determine wether derivatives should be calculated
C--------------------------- Calculate derivatives to concentration if we have a coupled problem solved with 
C--------------------------- Newtons method or a coupled problem with parameter estimation.
C--------------------------- calculate derivatives to head if we solve flow with newtons method or 
C--------------------------- if we do parameter estimation.

      IOCALCDEVT = 0


      IF (IODENS.EQ.1) THEN

          IF (LINMET(3,2).EQ.2 .OR. LINMET(2,2).EQ.2) THEN

              IOCALCDEVT = 1

          END IF !LINMET(3,2).EQ.2) .OR. LINMET(2,2).EQ.2

      END IF !IODENS.EQ.1

      IOCALCDEVF = 0

      IF (IODENS.EQ.1 .AND. LINMET(3,2).EQ.2) THEN

          IOCALCDEVF = 1

      END IF !IODENS.EQ.1 AND. LINMET(3,2).EQ.2
      

      IOCALCDEV = MAX(IOCALCDEVF,IOCALCDEVT)
C--------------------------- step 1: loop over elements

      DO L = 1, NUMEL

          NNUD = LNNDEL(L)
          NZONE=LXPAREL(L,1)         ! T zone of current element
          ISZ=ISOZ(NZONE)            ! Anisotropy degree
          LD = LDIM(L)
          LANI = LD*(LD + 1)/2   ! Maximum anisotropy of the current element                                         
          ISMAX=MAX(ISZ,LD)  
          IELEMTYP=LTYPE(L)       
          AREAL=AREA(L)
          

C--------------------------- step 2: identify transmissivity components

          DO IS=1,ISMAX
              TRACT(IS)=PAREL(L,IS)        
              IF (IOCALCDEVF.EQ.1 .AND. IS.GT.1. 
     &            .AND. (ISZ.LT.2.OR.IS.NE.3) ) THEN      

                   DPARELDH(IS,L)=DPARELDH(IS-1,L)

              END IF

              IF (IOCALCDEVT.EQ.1 .AND. IS.GT.1. 
     &            .AND. (ISZ.LT.2.OR.IS.NE.3) ) THEN

                  DPARELDC(IS,L) = DPARELDC(IS-1,L)
              END IF

          END DO !IS=1,ISMAX


C--------------------------- step 3: calculate velocity and its derivatives 

          VEL = 0D0
          DVDH_L = 0D0
          DVDC_L = 0D0

          CALL COMVEL_ELEMENT
     &        (LD        ,IOCALCDEVF,IODIM     ,ISMAX   
     &        ,LMXNDL    ,NPPEL     ,NNUD      ,NUMEL    ,NUMNP  ,L       
     &        ,KXX       ,DPARELDH  ,GRDFF    ,HAUX1  ,DVDH_L      
     &        ,VEL       ,TRACT)


C--------------------------- Computes buoyancy term contribution to velocity
C--------------------------- in density dependent flow.

          IF (IODENS.EQ.1) THEN

              CALL COMP_VEL_BUOYANCY
     &            (AREA     ,BUOYANCY ,COORD    ,DBUOYANCY,DVDC_L
     &            ,GP_COORD ,GRADLOC  ,IODIM    ,IOCALCDEV,ISOZ
     &            ,KXX      ,L        ,LDIM     ,LMXNDL   ,LNNDEL
     &            ,LTYPE    ,LXPAREL  ,MAXPG    ,NPAREL   ,NPPEL
     &            ,NUMEL    ,NUMNP    ,NZTRA    ,PAREL    ,POINTWEIGHT
     &            ,VEL)

          ELSE IF (IODENS.EQ.0 .AND. IOFLLI.EQ.1) THEN
          
              CALL COMP_VEL_GRAV
     &            (AREA       ,DVDH       ,GRAVEL     ,IOCALCDEVF
     &            ,IODIM      ,ISOZ       ,L          ,LDIM
     &            ,LMXNDL     ,LNNDEL     ,LXPAREL    ,NPAREL
     &            ,NPPEL      ,NUMEL      ,NZTRA      ,VEL)
     
          END IF !IODENS.EQ.1

C--------------------------- Stores velocity and its derivatives.

          VD(1:LD,L) = VEL(1:LD)

          IF (IOCALCDEVF.EQ.1) THEN
              DVDH(1:NNUD,1:LD,L) = DVDH_L(1:NNUD,1:LD)
          END IF

          IF (IOCALCDEVT.EQ.1) THEN
              DVDC(1:NNUD,1:LD,L) = DVDC_L(1:NNUD,1:LD)
          END IF


    
C--------------------------- step 4: calculate euclidean norm

          XNOR = 0D0

          DO I=1,LD

             XNOR = XNOR + VD(I, L) * VD(I, L)

          END DO !I=1,LD

          XNOR = DSQRT(XNOR)


C--------------------------- step 5:If Darcy's flow is too small, then the norm and 
C--------------------------- component's products are set to zero

          IF (XNOR.LT.1D-25) THEN

              XNORVD(L) = 0D0

              QXYZ(1:IDIMQ, L) = 0D0

          ELSE

C--------------------------- step 6: the norm and component'sproducts are stored in
C--------------------------- arrays XNORVD and QXYZ

              XNORVD(L) = XNOR

              DO I = 1,LD        !Terms q_i^2/|q|
                  QXYZ(I, L) = VD(I, L) * VD(I, L) / XNOR        !DIAGONAL
              END DO

              DO I = LD + 1, LANI      !Terms ( q_i * q_j ) /|q|
                  QXYZ(I, L) = VD(MIN(I - LD + 1, 3), L) *     ! NON-DIAGONAL
     &                      VD(MAX(I - 4, 1), L) / XNOR
              END DO

          ENDIF    !XNOR .LT. 1.D-25 Euclidian norm of Darcy's flow too small
      END DO !L=!,NUMEL

C------------------------- Writes velocities and their cross products when 
C------------------------- IFLAGS(5) is set to 1

       IF (IFLAGS(5).EQ.1) THEN
          WRITE(8, *) ' VELOCITY'
          DO L = 1, NUMEL
             WRITE(8, 2000) L,(VD(K, L), K = 1, IODIM)
          ENDDO
 2000     FORMAT(I7,3E21.10)
       ENDIF

      END SUBROUTINE COMVEL
