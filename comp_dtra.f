      SUBROUTINE COMP_DTRA
     &          (ACTH     ,AREA     ,BETAC    ,CAUX1    ,CAUX2
     &          ,CCALAN   ,CCALIT   ,CREF     ,DENSITY  ,DENSREF
     &          ,DTRA     ,DTRADFLU ,DTRADTRA ,DWDH     ,IDIMDTRA
     &          ,IODENS   ,IOINV    ,IOVRWC   ,ITPTVAR  ,KXX
     &          ,LINMET   ,LMXNDL   ,LNNDEL   ,NUMEL    ,NUMNP
     &          ,RETARD   ,THETAT   ,WATVOL   ,WSPECHEAT)

********************************************************************************
*
* PURPOSE
*
*  Manages the computation of DTRA matrix and its derivatives.
*
*
* DESCRIPTION
*
*  Manages the computation of DTRA matrix by calling
*  the routines that use WATVOL by nodes or elements as needed.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  COMP_DTRA_I            Computes DTRA with WATVOL by nodes
*  COMP_DTRA_L            Computes DTRA with WATVOL by elements
*
* HISTORY
*
*     AMS      7-2002     First coding
*     JHG      5-2003     Changes in calls to COMP_DTRA_I
*                         and COMP_DTRA_L.
*******************************************************************************

       IMPLICIT NONE
 
C------------------------- External

      INTEGER*4::IDIMDTRA,IODENS,IOINV,IOVRWC,ITPTVAR,LMXNDL,NUMEL,NUMNP

      REAL*8::BETAC,CREF,DENSREF,THETAT,WSPECHEAT

      REAL*8::DENS

      INTEGER*4::KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)

      REAL*8::ACTH(NUMEL)    ,AREA(NUMEL)   ,CCALAN(NUMNP)
     &       ,CCALIT(NUMNP)  ,CAUX1(NUMNP)  ,CAUX2(NUMNP)
     &       ,DENSITY(NUMEL) ,DTRA(NUMEL,IDIMDTRA)
     &       ,DTRADFLU(NUMEL,LMXNDL*LMXNDL)
     &       ,DTRADTRA(NUMEL,LMXNDL*LMXNDL)
     &       ,DWDH(MAX(1,(IOVRWC-1)*2*LMXNDL),NUMEL)
     &       ,RETARD(NUMEL)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)
      
       INTEGER*4::LINMET(3,2)

       IF (IOVRWC.LE.1) THEN

          CALL COMP_DTRA_L
     &        (ACTH     ,AREA     ,BETAC    ,CCALAN   ,CCALIT
     &        ,CREF     ,DENSREF  ,DTRA     ,IDIMDTRA ,IODENS
     &        ,IOVRWC   ,ITPTVAR  ,KXX      ,LMXNDL   ,LNNDEL
     &        ,NUMEL    ,NUMNP    ,RETARD   ,THETAT   ,WATVOL
     &        ,WSPECHEAT)

       ELSE

          CALL COMP_DTRA_I
     &        (ACTH     ,AREA     ,BETAC    ,CCALAN   ,CCALIT
     &        ,CREF     ,DENSREF  ,DTRA     ,IDIMDTRA ,IODENS
     &        ,IOVRWC   ,ITPTVAR  ,KXX      ,LMXNDL   ,LNNDEL
     &        ,NUMEL    ,NUMNP    ,RETARD   ,THETAT   ,WATVOL
     &        ,WSPECHEAT)

       END IF !IOVRWC.LE.1
          
       IF (LINMET(3,2).EQ.2 .OR. IOINV.EQ.1) THEN

          CALL COMP_DER_DTRA
     &        (AREA     ,BETAC    ,CAUX1    ,CAUX2    ,CREF
     &        ,DENSITY  ,DENSREF  ,DTRA     ,DTRADFLU ,DTRADTRA
     &        ,DWDH     ,IDIMDTRA ,IOVRWC   ,ITPTVAR  ,KXX
     &        ,LMXNDL   ,LNNDEL   ,NUMEL    ,NUMNP    ,THETAT
     &        ,WATVOL)

       END IF !LINMET(3,2).EQ.2 .OR. IOINV.EQ.1

       END SUBROUTINE COMP_DTRA

*****************************************************************************
*****************************************************************************

      SUBROUTINE COMP_DTRA_L
     &          (ACTH     ,AREA     ,BETAC    ,CCALAN   ,CCALIT
     &          ,CREF     ,DENSREF  ,DTRA     ,IDIMDTRA ,IODENS
     &          ,IOVRWC   ,ITPTVAR  ,KXX      ,LMXNDL   ,LNNDEL
     &          ,NUMEL    ,NUMNP    ,RETARD   ,THETAT   ,WATVOL
     &          ,WSPECHEAT)

     
    
********************************************************************************
*
* PURPOSE
*
*      Computes DTRA matrix its derivatives w. r. t. state variable 
*      with WATVOL by elements
*
* DESCRIPTION
*
*      Computes DTRA matrix and its derivatives w. r. t. state variables
*      with WATVOL by elements.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  CNST                   Interpolation functions gradient for a given element  
*                         nodes
*  DTRA                   Matrix of finite elements equations for transport
*                         problem related to storage.
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing the type of each element            
*  RETARD                 Array containing the retardation of every element
*  WATVOL                 Array containing the water content of every element
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMDTRA               Used to dimension array DTRA (second dimension).
*                             IDIMDTRA=LMXNDL if consistent scheme.
*                             IDIMDTRA=LMXNDL * (LMXNDL-1)/2 if non-consistent scheme.
*  IOCNST                 Scheme for mass storage term in transport problem
*                             IOCNST = 1 Lumped scheme (diagonal matrix).
*                             IOCNST = 0 Consistent scheme.
*  LMXNDL                 Maximum number of nodes per element                   
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*
* INTERNAL VARIABLES: SCALARS
*
*  NNUD                   Number of nodes of the current element                
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ZERO_ARRAY                                                                   
*
* HISTORY
*
*     JCA      5-1998     First coding
*     AMS      3-1999     Revision (inclusion of WATVOL)
*     JHG      5-2003     Elimination of bandwise assembly.
*                         New retard calculation.
********************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMDTRA,IODENS,IOVRWC,ITPTVAR,LMXNDL,NUMEL,NUMNP

      REAL*8::BETAC,CREF,DENSREF,THETAT,WSPECHEAT

      REAL*8::DENS

      INTEGER*4::KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)

      REAL*8::ACTH(NUMEL)             ,AREA(NUMEL),CCALAN(NUMNP)
     &       ,CCALIT(NUMNP)
     &       ,DTRA(NUMEL,IDIMDTRA)    ,RETARD(NUMEL)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)     

C------------------------- Internal

      REAL*8::AREAL,CAVG,RCRR,RET,RHO,TH1,WTV

      INTEGER*4::I,INODE,L,NNUD

C------------------------- First executable statement

      TH1 = 1D0 - THETAT

C------------------------- Starts loop by elements

      DO L=1,NUMEL

          NNUD = LNNDEL(L)

C------------------------- Computes density at k+1-thetat

          IF (IODENS.GT.0) THEN

              CAVG = 0D0

              DO I=1,NNUD
                INODE = KXX(I,L)
                CAVG = CAVG + TH1*CCALIT(INODE) + THETAT*CCALAN(INODE)
              END DO !I=1,NNUD

              CAVG = CAVG/NNUD
              RHO = DENS(DENSREF,BETAC,CAVG,CREF)

          ELSE

              RHO = 1D0

          END IF !(IODENS.GT.0)

          RET = RETARD(L)*ACTH(L)
          AREAL = AREA(L)
          WTV = WATVOL(1,L,3)

          IF (ITPTVAR.EQ.0) THEN !solute transport

              RCRR = RHO*(WTV + RET)*AREAL

          ELSE !energy transport

              IF (IODENS.EQ.1) THEN

                  RCRR = (RHO*WTV +RET/WSPECHEAT)*AREAL

              ELSE
C------------------------- If density is constant, water density appears
C------------------------- in this term (see TRANSDENS Gu’a r‡pida, 3.7).
                  RCRR = (RHO*WTV +RET/(DENSREF*WSPECHEAT))*AREAL

              END IF !IODENS.EQ.1

          END IF !ITPTVAR.EQ.0

          DTRA(L,1:NNUD) = DTRA(L,1:NNUD) + RCRR/NNUD

      END DO !L=1,NUMEL

      END SUBROUTINE COMP_DTRA_L

*****************************************************************************
*****************************************************************************

      SUBROUTINE COMP_DTRA_I
     &          (ACTH     ,AREA     ,BETAC    ,CCALAN   ,CCALIT
     &          ,CREF     ,DENSREF  ,DTRA     ,IDIMDTRA ,IODENS
     &          ,IOVRWC   ,ITPTVAR  ,KXX      ,LMXNDL   ,LNNDEL
     &          ,NUMEL    ,NUMNP    ,RETARD   ,THETAT   ,WATVOL
     &          ,WSPECHEAT)

********************************************************************************
*
* PURPOSE
*
*      Computes DTRA matrix with WATVOL by nodes
*
* DESCRIPTION
*
*      Computes DTRA matrix  with WATVOL by nodes. 
*                       THERE IS NO CONSISTENT SCHEME. ONLY LUMPED.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  CNST                   Interpolation functions gradient for a given element  
*                         nodes
*  DTRA                   Matrix of finite elements equations for transport
*                         problem related to storage.
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing the type of each element            
*  RETARD                 Array containing the retardation of every element
*  WATVOL                 Array containing the water content of every element
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMDTRA               Used to dimension array DTRA (second dimension)       
*  IOCNST                 Scheme for mass storage term in transport problem     
*  LMXNDL                 Maximum number of nodes per element                   
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*
* INTERNAL VARIABLES: SCALARS
*
*  NNUD                   Number of nodes of the current element                
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ZERO_ARRAY                                                                   
*
* HISTORY
*
*     JCA      5-1998     First coding
*     AMS      3-1999     Revision (inclusion of WATVOL)
*     JHG      5-2003     Elimination of bandwise assembly.
*                         New retard calculation.
********************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMDTRA,IODENS,IOVRWC,ITPTVAR,LMXNDL,NUMEL,NUMNP

      REAL*8::BETAC,CREF,DENSREF,THETAT,WSPECHEAT

	REAL*8::DENS

      INTEGER*4::KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)

      REAL*8::ACTH(NUMEL)     ,AREA(NUMEL)           ,CCALAN(NUMNP)
     &       ,CCALIT(NUMNP)   ,DTRA(NUMEL,IDIMDTRA)  ,RETARD(NUMEL)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3)
        
C------------------------- Internal

      INTEGER*4::I,INODE,L,NNUD

      REAL*8::AREAL,CONCTH1,DENSNODE,RCRR,RET,TH1,WTVI

C------------------------- First executable statement


      TH1 = 1D0 - THETAT

C------------------------- Starts loop by elements

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          RET = RETARD(L)*ACTH(L)
          AREAL = AREA(L)
          
          DO I=1,NNUD

              WTVI = WATVOL(I,L,3) ! k+1-Theta
              INODE=KXX(I,L)

	        IF (IODENS.EQ.1) THEN !density at k+1-thetat
                  CONCTH1 = TH1*CCALIT(INODE) + THETAT*CCALAN(INODE)
                  DENSNODE = DENS(DENSREF,BETAC,CONCTH1,CREF)
	        ELSE

	            DENSNODE = 1D0

	        END IF !IODENS.EQ.1


              IF (ITPTVAR.EQ.0) THEN !solute transport
                           
                  RCRR = DENSNODE*(WTVI+RET)*AREAL

              ELSE

C------------------------- If density is constant, water density appears
C------------------------- in this term (see TRANSDENS Gu’a r‡pida, 3.7).
                  RCRR = (DENSNODE*WTVI + RET/(DENSREF*WSPECHEAT))*AREAL

              END IF !ITPTVAR.EQ.0

              DTRA(L,I) = DTRA(L,I) + RCRR/NNUD
      
          END DO !I=1,NNUD

      END DO !L=1,NUMEL

      END SUBROUTINE COMP_DTRA_I
