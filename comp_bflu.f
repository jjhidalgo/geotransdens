      SUBROUTINE COMP_BFLU
     &          (AREA     ,BETAC    ,BFLU     ,BUOYANCY ,CAUX1
     &          ,CONCFLOW ,COORD    ,CREF     ,DBFLUDFLU,DBFLUDTRA
     &          ,DBUOYANCY,DENSITY  ,DENSREF  ,DFLUDFLU ,DFLUDTRA
     &          ,DNODALRH ,DPARELDH  ,EPSFLU
     &          ,GP_COORD ,GRADLOC  ,HAUX1    ,HCALAN   ,IBCOD
     &          ,IBTCO    ,INDSSTR  ,IOCALCDEVF         ,IOCALCDEVT
     &          ,IOCONSRC ,IODENS   ,IODIM    ,IOFLLI   ,IONEWT   ,ISOZ
     &          ,KXX      ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE
     &          ,LXPAREL  ,MAXPG    ,NPAREL   ,NPPEL    ,NPPNP
     &          ,NTYPAR   ,NTYPEL   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &          ,NZTRA    ,PAREL    ,PARNP    ,POINTWEIGHT
     &          ,THETAF   ,THETAT   ,gravel   ,grdff)


********************************************************************************
*
* PURPOSE Computes the RHS term of flow equation (sinks/sources term). Will be
*         corrected because of presc. head and mixed boundary conditions
*
* DESCRIPTION Subroutine is quite simple. Can be summarized as follows:
*            
*             - Step 0: Declaration of variables
*             - Step 1: Initialization of RHS of flow equation: Array BFLU
*             - Step 2: Adds areal recharge contribution. 
*             - Step 3: Adds prescribed flow contribution
*             - Step 4: Adds gravity flow if liquid pressure is the state 
*                       variable
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  BFLU                   Right hand side of flow discretized equation.         
*  FLOWGRAV               Array containing gravity flow for all nodes
*  IBCOD                  Flow boundary condition index                         
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*  LXPAREL                Array containing zone numbers for a given             
*                         element parameter                                     
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters                                            
*  PAREL                  Parameter values at every element and current time    
*                         for all nodal parameters (each value is computed as   
*                         the product of up to four terms:                      
*                           elem. coeff*zonal value*time funct.*nonl. funct. )  
*  PARNP                  Parameter values at every node and current time for   
*                         all nodal parameters (each value is computed as the   
*                         product of up to four terms:                          
*                           nodal coeff*zonal value*time funct.*nonl. funct. )  
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IOPRHED                Indicates whether the flow state variable state is    
*                         preasure (set to 1) or head (set to 0)                
*  LMXNDL                 Maximum number of nodes per element                   
*  NPAREL                 Number of element parameters in current problem       
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
*  NPPNP                  Total number of parameters by nodes (not confuse      
*                         with NPARNP, because in this casethere is no          
*                         difference between a given parameter in steady or tr.)
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*
* INTERNAL VARIABLES: SCALARS
*
*  I                      Dummy counter of nodal points           
*  IB                     Boundary condition of a given node
*  L                      Dummy counter of elements                             
*  NNUD                   Number of nodes of the current element                
*  NODE                   Nodal point ident.                    
*  NZONEARR               Recharge zone to which element L belongs to
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ZERO_ARRAY                                                                   
*
* HISTORY: First coding                         German Galarza (Nov-97)
*          Revision and header inclusion        AAR (Nov-00)
*          PARNP,PAREL inclusion                AAR (Feb-02)
*
********************************************************************************

C---------------------------  Step 0: Declaration of variables

      IMPLICIT NONE

C---------------------------   EXTERNAL VARIABLES: SCALARS

      REAL*8::BETAC,DENSREF,CREF,THETAF,EPSFLU,THETAT

      INTEGER*4::INDSSTR,IOCONSRC,IODIM,IOFLLI,IONEWT,NUMNP,NUMEL,NPPEL
     &          ,NTYPAR,NPAREL,IODENS,NPPNP,LMXNDL,MAXPG,NZTRA,NTYPEL

C---------------------------   EXTERNAL VARIABLES: ARRAYS

      REAL*8::AREA(NUMEL),BFLU(NUMNP),BUOYANCY(IODIM,LMXNDL,NUMEL)
     &       ,CAUX1(NUMNP),CONCFLOW(NUMNP),COORD(NUMNP,3)
     &       ,DBFLUDFLU(NUMNP),DBFLUDTRA(NUMNP)
     &       ,DBUOYANCY(IODIM,LMXNDL*LMXNDL,NUMEL)
     &       ,DENSITY(NUMEL),DNODALRH(NUMNP,4)
     &       ,DFLUDFLU(NUMEL,LMXNDL*LMXNDL)
     &       ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL)
     &       ,DPARELDH(NPPEL,NUMEL)
     &       ,GP_COORD(6,8,IODIM)
     &       ,GRADLOC(IODIM,LMXNDL,MAXPG)
     &       ,HAUX1(NUMNP),HCALAN(NUMNP)
     &       ,PAREL(NUMEL,NPPEL),PARNP(NUMNP,NPPNP)
     &       ,POINTWEIGHT(MAXPG,NTYPEL)
     &       ,gravel(numel,iodim),grdff(iodim,lmxndl,numel)

      INTEGER*4 NZONE_PAR(NTYPAR),LXPAREL(NUMEL,NPAREL),LNNDEL(NUMEL)
     &         ,IBCOD(NUMNP),IBTCO(NUMNP),ISOZ(NZTRA),KXX(LMXNDL,NUMEL)
     &         ,LDIM(NUMEL),LTYPE(NUMEL)


C---------------------------   INTERNAL VARIABLES: SCALARS

      REAL*8::ALFAX,AREAL,CAUD,CEXT,CNODE,DENS,DENSREC,DENSNODE,DALFADH
     &       ,DHEADDH,HEAD,RECHARGE,THTF1

      INTEGER*4::I,IB,IBT,IOCALCDEV,IOCALCDEVF,IOCALCDEVT,L,NZONEARR
     &           ,NNUD,NZONECEXT


C--------------------------- Step 1: Initialization of RHS of flow
C---------------------------         equation: Array BFLU

      BFLU = 0D0
      IOCALCDEV = 0
C--------------------------- Initialization of variables related
C--------------------------- to the computation of derivatives.

      IOCALCDEV = MAX(IOCALCDEVF, IOCALCDEVT)

      IF (IOCALCDEVF.GT.0) THEN
          DBFLUDFLU = 0D0
      END IF

      IF (IOCALCDEVT.GT.0) THEN
          DBFLUDTRA = 0D0
      END IF

C--------------------------- Step 2: Adds areal recharge contribution

      IF (NZONE_PAR(3).NE.0) THEN ! Recharge zones are considered

          DO L=1,NUMEL

              NZONEARR = LXPAREL(L,3+INDSSTR)
              IF (NZONEARR.NE.0) THEN

                  AREAL = AREA(L)
                  NNUD = LNNDEL(L)
                  RECHARGE = PAREL(L,8)
                  NZONECEXT = LXPAREL(L,10)

                  IF (NZONECEXT.GT.0 .AND. RECHARGE.GT.0) THEN

                      CEXT = PAREL(L,15)

                  ELSE

                      CEXT = 0D0

                  END IF !NZONECEXT.GT.0 .AND. RECHARGE.GT.0

                  IF (IODENS.EQ.1) THEN

                      DENSREC =  DENS(DENSREF,BETAC,CEXT,CREF)

                  ELSE

                      DENSREC = 1.D0

                  END IF !IODENS.EQ.1

C--------------------------- All nodes in the element receive the same contribution.
C--------------------------- Recharge does not contributes to the derivative.

                  BFLU(KXX(1:NNUD,L)) = BFLU(KXX(1:NNUD,L))
     &                                + DENSREC*RECHARGE*AREAL/NNUD

              END IF ! L belongs to a recharge zone
          END DO ! L=1,NUMEL
      END IF !NZONE_PAR(3).NE.0 Recharge zones are considered


C--------------------------- Step 3: Boundary conditions.

      THTF1 = THETAF - 1D0

      DO I=1,NUMNP

          IB = IBCOD(I)
          IBT = IBTCO(I)

C--------------------------- Prescribed flow.

          IF (IB.EQ.2) THEN

              CAUD = PARNP(I,2)

              IF (IODENS.EQ.1) THEN

                  IF (CAUD.GT.0) THEN

                      CNODE = PARNP(I,4) !external concentration

                  ELSE

                      CNODE = CAUX1(I) !internal concentration

                  END IF !PARNP(NODE,2).GT.0

                  DENSNODE = DENS(DENSREF,BETAC,CNODE,CREF)

              ELSE

                  DENSNODE = 1.D0

              END IF !IODENS.EQ.1

              BFLU(I) = BFLU(I) + DENSNODE*CAUD

              IF (IOCALCDEVF.GT.0) THEN

                  DBFLUDFLU(I) = DBFLUDFLU(I)
     &                          + EPSFLU*DENSNODE*DNODALRH(I,2)

              END IF

              IF (IOCALCDEVT.GT.0 .AND. CAUD.LT.0) THEN

                  DBFLUDTRA(I) = DBFLUDTRA(I)
     &                          + THETAT*BETAC*DENSNODE*CAUD

              END IF

          END IF !IB.EQ.2 Flow


C--------------------------- leakage,mixed 

          IF (IB.GT.2) THEN 

              ALFAX = PARNP(I,3)
              HEAD = PARNP(I,1)
              CAUD = ALFAX*(HEAD - HAUX1(I))

              IF (IB .EQ. 4) THEN

                  CAUD = CAUD + PARNP(I,2)

              END IF !IB .EQ. 4

              DHEADDH = EPSFLU*DNODALRH(I,1)
              DALFADH = EPSFLU*DNODALRH(I,3)


              IF (IODENS.EQ.1) THEN

                  IF (CAUD.GT.0) THEN

                      CNODE = PARNP(I,4) !external concentration

                  ELSE

                      CNODE = CAUX1(I)

                  END IF !CAUD.GT.0

                  DENSNODE = DENS(DENSREF,BETAC,CNODE,CREF)

              ELSE

                  DENSNODE = 1.D0

              END IF !IODENS.EQ.1

              IF (IONEWT.EQ.0) Then

                  BFLU(I) = BFLU(I)
     &                     + DENSNODE*ALFAX*(HEAD + THTF1*HCALAN(I))

              ELSE

                  BFLU(I) = BFLU(I) + DENSNODE*ALFAX*(HEAD-HAUX1(I))

              END IF !IONEWT.EQ.0

C--------------------------- Presc. flow contribution for mixed b. c.
C--------------------------- (the same for Newton and Picard).

              IF (IB .EQ. 4) THEN

                  BFLU(I) = BFLU(I) + DENSNODE*PARNP(I,2)

              END IF !IB .EQ. 4

              IF (IOCALCDEVF.GT.0) THEN

                  DBFLUDFLU(I) = DBFLUDFLU(I)
     &                         + DENSNODE*(ALFAX*(DHEADDH-THETAF)
     &                                     +DALFADH*(HEAD-HAUX1(I)))

C--------------------------- Presc. flow contribution for mixed b. c.
                  IF (IB .EQ. 4) THEN

                      DBFLUDFLU(I) = DBFLUDFLU(I)
     &                              + EPSFLU*DENSNODE*DNODALRH(I,2)

                  END IF !IB .EQ. 4

              END IF !IOCALCDEVF.GT.0


              IF (IOCALCDEVT.GT.0 .AND. CAUD.LT.0) THEN

                  DBFLUDTRA(I) = DBFLUDTRA(I) 
     &                     + THETAT*BETAC*DENSNODE*ALFAX*(HEAD-HAUX1(I))

C--------------------------- Presc. flow contribution for mixed b. c.
                  IF (IB .EQ. 4) THEN

                      DBFLUDTRA(I) = DBFLUDTRA(I)
     &                              + THETAT*BETAC*DENSNODE*PARNP(I,2)

                  END IF !IB .EQ. 4

              END IF !IOCALCDEVT.GT.0 .AND. CAUD.LT.0

          END IF !IB.GT.2


          IF (IOCONSRC.EQ.1) THEN

              IF (IB.NE.1) THEN

                  IF (IBT.EQ.1 .OR. IBT.EQ.4 .OR. IBT.EQ.5) THEN

                      BFLU(I) = BFLU(I) + CONCFLOW(I)

                  END IF !IBT.EQ.1 .OR. IBT.EQ.4 .OR. IBT.EQ.5

              END IF !IB.NE.1

          END IF !IOCONSRC.EQ.1

      END DO !I=1,NUMNP


      IF (IODENS.EQ.1) THEN

          CALL COMP_BFLU_BUOYANCY
     &        (AREA     ,BETAC    ,BFLU     ,BUOYANCY ,COORD
     &        ,DBUOYANCY,DENSITY  ,DFLUDTRA ,GP_COORD ,GRADLOC
     &        ,IOCALCDEV,IODIM    ,ISOZ     ,KXX      ,LDIM
     &        ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXPAREL  ,MAXPG
     &        ,NPAREL   ,NPPEL    ,NTYPEL   ,NUMEL    ,NUMNP
     &        ,NZTRA    ,PAREL    ,POINTWEIGHT        ,THETAT)

      ELSE

          IF (IOFLLI.EQ.1) THEN

              CALL COMP_FLOWGRAV
     &            (AREA       ,BFLU       ,DFLUDFLU   ,DPARELDH
     &            ,GRAVEL     ,GRDFF      ,IOCALCDEVF ,IODIM
     &            ,KXX        ,LDIM       ,LMXNDL     ,LNNDEL
     &            ,NUMNP      ,NPPEL      ,NUMEL      ,PAREL)

          END IF !IOFLLI.EQ.1

      END IF ! IODENS.EQ.1

      END SUBROUTINE COMP_BFLU
