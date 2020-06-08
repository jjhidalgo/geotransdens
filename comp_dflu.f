      SUBROUTINE COMP_DFLU
     &          (AREA     ,BETAC    ,CAUX1    ,CREF     ,DENSITY
     &          ,DENSREF  ,DFLU     ,DFLUDFLU ,DFLUDTRA ,DPARELDC
     &          ,DPARELDH ,EPSFLU   ,EPSTRA   ,HAUX2    ,IDIMDFLU
     &          ,IFLAGS   ,IOCALCDEVF         ,IOCALCDEVT
     &          ,IODENS   ,IOVRWC   ,KXX      ,LMXNDL   ,LNNDEL
     &          ,LTYPE    ,MAINF    ,NFLAGS   ,NPPEL    ,NUMEL
     &          ,NUMNP    ,PAREL    ,THETAT)



********************************************************************************
*
* PURPOSE Calculation and assembly of storavity matrix DFLU and its contribution
*         to the non linear problem derivatives matrix DFLUDFLU
*
* DESCRIPTION Performs two calculations: 
*                          1) DFLU=STGC*AREA/NNUD (ALWAYS)
*                          2) DFLUDFLU(I)=EPSFLU*DERSTGC*AREA*DeltaH/(TINC*NNUD)
*             Subroutine can be summarized as follows
*            
*             - Step 0: Declaration of variables
*             - Step 1: Initialises array DFLU
*             - Step 2: Begins main loop over mesh elements
*               - Step 2.1: Loop over element nodal points
*                 - Step 2.1.1: Assigns DFLU component
*                 - Step 2.1.2: Non linear case. Adds the contribution to 
*                   DFLUDFLU
*             - Step 3: Debugging purposes. Echoes DFLU to MAINF
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  DFLUDFLU               Contains the derivatives of AFLU and DFLU matrices    
*                         with respect to head. This term appears on the left   
*                         hand side of both flow direct and inverse problems,   
*                         as well as on the right hand side of flow inverse     
*                         problem.                                              
*  DERSTGH                Derivatives of storavity wrt. head level/pressure
*  DFLU                   Matrix of finite elements equations for flow          
*                         problem related to storage term.                      
*  HCALAN                 Head level at previous time                           
*  HCALIT                 Computed heads in last iteration                      
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*  PAREL                  Parameter values at every element and current time    
*                         for all nodal parameters (each value is computed as   
*                         the product of up to four terms:                      
*                           elem. coeff*zonal value*time funct.*nonl. funct. )  
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMDADFLU             Used to dimension array DFLUDFLU     
*  IDIMDFLU               Used to dimension array DFLU                          
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  LMXNDL                 Maximum number of nodes per element                   
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  TINC                   Current time increment                                
*
* INTERNAL VARIABLES: SCALARS
*
*  IDDF                   Used to dimension array DERAD
*  INUD                   Dummy counter of element nodal points
*  I                      Unimportant dummy counter
*  J                      Unimportant dummy counter
*  L                      Dummy counter of elements                 
*  NNUD                   Number of nodes of the current element                
*  NODE                   Identifies current node
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ZERO_ARRAY                                                                   
*
* HISTORY : First coding                                 German Galarza (Nov-97)
*           New code                                     AAR (Oct-98)
*           Inclusion of PAREL                           AAR (Feb-02)
*
********************************************************************************

C-------------------- Step 0: Declaration of variables

      IMPLICIT NONE


C-------------------- EXTERNAL VARIABLES: SCALARS

      INTEGER*4::IDIMDFLU ,IOCALCDEVF         ,IOCALCDEVT,IODENS
     &          ,IOVRWC   ,LMXNDL   ,MAINF    ,NFLAGS   ,NPPEL
     &          ,NUMEL    ,NUMNP

      REAL*8::CREF,DENSREF,EPSFLU,EPSTRA,THETAT

      REAL*8::DENS

C-------------------- EXTERNAL VARIABLES: ARRAYS

      INTEGER*4::IFLAGS(NFLAGS) ,KXX(LMXNDL,NUMEL)
     &          ,LNNDEL(NUMEL)  ,LTYPE(NUMEL)

      REAL*8::AREA(NUMEL)                 ,CAUX1(NUMNP)
     &       ,DENSITY(NUMEL)              ,DFLU(NUMEL,IDIMDFLU)
     &       ,HAUX2(NUMNP)                ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL)
     &       ,DFLUDFLU(NUMEL,LMXNDL*LMXNDL)
     &       ,DPARELDC(NPPEL,NUMEL)       ,DPARELDH(NPPEL,NUMEL)
     &       ,PAREL(NUMEL,NPPEL)

C-------------------- INTERNAL VARIABLES: SCALARS

      REAL*8::AREAL,AREALN,DFLUI,BETAC,DDFLUDH,DDFLUDC,DENSNODE,DENSL
     &       ,HI,STG

      INTEGER*4::L,NNUD,I,IJ_POS,J,LTYP,IOCALCDEV,INODE


C-------------------- Step 1: Initialises array DFLU

      DFLU = 0D0
      
   
      IOCALCDEV = MAX(IOCALCDEVF,IOCALCDEVT)



C-------------------- Step 2: Begins main loop over mesh elements

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          LTYP = LTYPE(L)
          AREAL = AREA(L)
          AREALN = AREAL/NNUD

          DENSNODE = 1D0
          DENSL = 1D0

          STG = PAREL(L,7)*AREALN

          DO I=1,NNUD

              INODE= KXX(I,L)

C-------------------- Step 2.1: Assigns DFLU component       

              SELECT CASE(IOVRWC)


              CASE (0,1) ! Constant or elementwise

                  IF (IODENS.EQ.1) THEN

                      DENSL = DENSITY(L)

                  END IF !IODENS.EQ.1

                  DFLUI = STG*DENSL

              CASE (2) ! Nodewise

                  IF (IODENS.EQ.1) THEN

                      DENSNODE = DENS(DENSREF,BETAC,CAUX1(INODE),CREF)

                  END IF !IODENS.EQ.1

                  DFLUI = STG*DENSNODE

              END SELECT !IOVRWC

              DFLU(L,I) = DFLUI

C-------------------- Computes derivatives.

              IF (IOCALCDEV.NE.0) THEN

                  HI = HAUX2(INODE)

                  SELECT CASE (IOVRWC)

                  CASE (0,1) ! Constant or elementwise

                      DDFLUDH = DENSL*EPSFLU*DPARELDH(7,L)*AREALN*HI

                      DDFLUDC = (EPSTRA*DPARELDC(7,L)*AREAL
     &                           + THETAT*BETAC*DFLUI/NNUD)*HI

                  CASE (2) ! Nodewise

                      DDFLUDH = DENSNODE*EPSFLU*DPARELDH(7,L)*AREAL*HI

                      DDFLUDC = (EPSTRA*DPARELDC(7,L)*AREAL
     &                           + THETAT*BETAC*DFLUI)*HI

                  END SELECT !IOVRWC


                  DO J=1,NNUD

                      IJ_POS = (I-1)*NNUD + J

C-------------------- If nodewise, derivatives are null except when I.EQ.J

                      IF (I.EQ.J .OR. IOVRWC.NE.2) THEN

                          IF (IOCALCDEVF.EQ.1) THEN

                              DFLUDFLU(L,IJ_POS) = DFLUDFLU(L,IJ_POS)
     &                                            + DDFLUDH

                          END IF !IOCALCDEVF

                          IF (IOCALCDEVT.EQ.1) THEN 

                              DFLUDTRA(L,IJ_POS) = DFLUDTRA(L,IJ_POS)
     &                                            + DDFLUDC

                          END IF !IOCALCDEVT

                      END IF !I.EQ.J .OR. IOVRWC.NE.2

                  END DO !J=1,NNUD

              END IF !IOCALCDEV.NE.0

          END DO !I=1,NNUD

      END DO  !L=!,NUMEL

C-------------------- Step 3: Debugging purposes. Echoes DFLU to MAINF

      IF (IFLAGS(14).EQ.1) THEN

         WRITE(MAINF,'(A20)') ' DFLU EN COMP_DFLU'

         DO I=1,NUMEL

             DO J=1,IDIMDFLU

                WRITE(MAINF,'(2I5,E20.13)') I,J,DFLU(I,J)

             END DO !J=1,IDIMDFLU

         END DO !I=1,NUMEL

      END IF !IFLAGS(14).EQ.1

      END SUBROUTINE COMP_DFLU
