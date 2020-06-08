      SUBROUTINE COMP_DER_FLOW
     &          (AFLU     ,BETAC    ,CAUDAL   ,CFLU     ,CREF
     &          ,DENSREF  ,DFLU     ,DFLUDFLU ,DFLUDTRA ,DNODALRH
     &          ,DQDFLU   ,DQDTRA   ,EPSFLU   ,EPSTRA   ,HAUX1
     &          ,IBCOD    ,IDIMAFLU ,IDIMDFLU ,IOATK    ,IODENS
     &          ,ISOLFL   ,KXX      ,LINMET   ,LMXNDL   ,LNNDEL
     &          ,NPPNP    ,NUMEL    ,NUMNP    ,PARNP    ,TINC)
     
********************************************************************************
*
* PURPOSE
*
*      Computes derivatives of nodal flow at every node w. r. t.state varables
*      at k or k+1
*
* DESCRIPTION
*
*      Computes derivatives of nodal flow at every node w. r. t.state varables
*      at k or k+1
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AFLU                   Matrix of finite elements equations for flow problem
*                         No boundary conditions are included on it.
*  CAUDAL                 Input/output flow at every node.
*  DFLU                   Matrix of finite elements equations for flow
*                         problem related to storage term.
*  HAUX1                  Array containing heads, ponderated by THETAF time
*                         weight. Is equal to THETAF*HCAL+(1-THETAF)*HCALAN
*  IBCOD                  Flow boundary condition index
*  PARNP                  Parameter values at every node and current time for
*                         all nodal parameters (each value is computed as the
*                         product of up to four terms:
*                            coeff*zonal value*time funct.*nonl. funct. )
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMDFLU               Used to dimension array DFLU
*  ISOLFL                 If 1, steady flow has been solved. If 2, transient.
*  IOCNSF                 Scheme for storage term in flow problem
*                                 1.Lumped
*                                 2.Consistent.
*                         numbers of two nodes belonging to the same element)
*  NUMNP                  Number of nodes
*  THETAF                 Time weighting parameter for flow problems
*  TINC                   Current time increment
*
* INTERNAL VARIABLES: SCALARS
*
*  I                      Nodal counter
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  MUL_SS
*  MUL_TT
*
* HISTORY
*
*     JHG     10-2003     First coding.
********************************************************************************


      IMPLICIT NONE

      INTEGER*4::IDIMAFLU ,IDIMDFLU ,IOATK    ,IODENS   ,ISOLFL   
     &          ,LMXNDL   ,NPPNP    ,NUMEL    ,NUMNP  

      INTEGER*4:: IBCOD(NUMNP)  ,KXX(LMXNDL,NUMEL)
     &           ,LNNDEL(NUMEL) ,LINMET(3,2)

      REAL*8::BETAC    ,CREF     ,DENSREF  ,EPSFLU   ,EPSTRA   ,TINC    

      REAL*8::AFLU(NUMEL,IDIMAFLU)  ,CAUDAL(NUMNP)
     &       ,CFLU(NUMEL,IDIMDFLU)
     &       ,DFLU(NUMEL,IDIMDFLU)  ,DFLUDFLU(NUMEL,LMXNDL*LMXNDL)
     &       ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL),DNODALRH(NUMNP,4)
     &       ,DQDFLU(NUMEL,LMXNDL*LMXNDL)  ,DQDTRA(NUMEL,LMXNDL*LMXNDL)
     &       ,HAUX1(NUMNP)
     &       ,PARNP(NUMNP,NPPNP)


      CALL COMP_DQDTRA
     &    (BETAC    ,CAUDAL   ,CFLU     ,CREF     ,DENSREF
     &    ,DFLUDTRA ,DQDTRA   ,EPSTRA   ,IBCOD    ,IDIMDFLU
     &    ,IOATK    ,IODENS   ,ISOLFL   ,KXX      ,LMXNDL
     &    ,LNNDEL   ,NUMEL    ,NPPNP    ,NUMNP    ,PARNP
     &    ,TINC)   
          
      IF (LINMET(3,2).EQ.2) THEN

          CALL COMP_DQDFLU
     &        (AFLU     ,BETAC    ,CAUDAL   ,CREF     ,DENSREF
     &        ,DFLU     ,DFLUDFLU ,DNODALRH ,DQDFLU   ,EPSFLU
     &        ,HAUX1    ,IBCOD    ,IDIMAFLU ,IDIMDFLU ,IOATK
     &        ,IODENS   ,ISOLFL   ,KXX      ,LMXNDL   ,LNNDEL
     &        ,NPPNP    ,NUMEL    ,NUMNP    ,PARNP    ,TINC)

      END IF


      END SUBROUTINE COMP_DER_FLOW

C***********************************************************************
C***********************************************************************
C***********************************************************************
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE COMP_DQDTRA
     &          (BETAC    ,CAUDAL   ,CFLU     ,CREF     ,DENSREF
     &          ,DFLUDTRA ,DQDTRA   ,EPSTRA   ,IBCOD    ,IDIMDFLU
     &          ,IOATK    ,IODENS   ,ISOLFL   ,KXX      ,LMXNDL
     &          ,LNNDEL   ,NUMEL    ,NPPNP    ,NUMNP    ,PARNP
     &          ,TINC)   

********************************************************************************
*
* PURPOSE
*
*      Computes derivatives of nodal flow at every node w. r. t.state varables
*      at k or k+1
*
* DESCRIPTION
*
*      Computes derivatives of nodal flow at every node w. r. t.state varables
*      at k or k+1
*
********************************************************************************
      IMPLICIT NONE

      INTEGER*4::I        ,IDIMDFLU ,INODE    ,IOATK    ,IODENS
     &          ,ISOLFL   ,J        ,JNODE    ,L        ,LMXNDL
     &          ,NNUD     ,NPPNP    ,NUMEL    ,NUMNP    ,IIPOS
     &          ,IJPOS    ,JIPOS
     

      REAL*8::BETAC,CEXT,CREF,DENS,DENSEXT,DENSREF,EPSTRA,TINC
     &       ,TINCINV,EPST,ATKFACTOR

      INTEGER*4:: IBCOD(NUMNP),KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)

      REAL*8::CAUDAL(NUMNP),CFLU(NUMEL,IDIMDFLU)
     &       ,DFLUDTRA(NUMEL,LMXNDL*LMXNDL),DQDTRA(NUMEL,LMXNDL*LMXNDL)
     &       ,PARNP(NUMNP,NPPNP)
  
    

C------------------------ Initialization
C------------------------ DQDTRA is computed by adding the nodal values
C------------------------ so, the initialization to zero is mandatory.
      
      DQDTRA = 0.D0

      IF (TINC.NE.0) TINCINV = 1D0/TINC

C------------------------ Needed so that DQDTRA can be converted into
C------------------------ derivative respecto to conc at k or k+1

      IF (IOATK.EQ.1) THEN

          EPST = (1D0 - EPSTRA)/EPSTRA
          ATKFACTOR = -1D0

	ELSE

          EPST = 1D0
          ATKFACTOR = 1D0

      END IF !IOATK.EQ.1

C----------------------------------------------
C------------------------ Prescribed head -----
C----------------------------------------------

      DO L=1,NUMEL

          NNUD = LNNDEL(L) 

          DO I=1,NNUD-1

              INODE=KXX(I,L)

              IF (IBCOD(INODE).EQ.1 .AND. CAUDAL(INODE).GT.0) THEN

                  IIPOS = (I-1)*NNUD + I
                  DQDTRA(L,IIPOS) = DQDTRA(L,IIPOS)
     &                            + EPST*DFLUDTRA(L,IIPOS)

              
                  IF (ISOLFL.GT.1 .AND. IODENS.EQ.1) THEN
                      
                      DQDTRA(L,IIPOS) = DQDTRA(L,IIPOS)
     &                                + CFLU(L,I)*TINCINV*ATKFACTOR

                  END IF !(ISOLFL.GT.1)

             END IF ! IBCOD(INODE).EQ.1 .AND. CAUDAL(INODE).GT.0

              DO J=I+1,NNUD

                  JNODE = KXX(J,L)
        
                  IF (IBCOD(INODE).EQ.1 .AND. CAUDAL(INODE).GT.0) THEN

                     IJPOS = (I-1)*NNUD + J
                     DQDTRA(L,IJPOS) = DQDTRA(L,IJPOS)
     &                               + EPST*DFLUDTRA(L,IJPOS)

                  END IF

                  IF (IBCOD(JNODE).EQ.1 .AND. CAUDAL(JNODE).GT.0) THEN

                     JIPOS = (J-1)*NNUD + I
                     DQDTRA(L,JIPOS) = DQDTRA(L,JIPOS)
     &                               + EPST*DFLUDTRA(L,JIPOS)

                  ENDIF 

              END DO !J=I+1,NNUD

          END DO !I=1,NNUD-1

C------------------------ Last Diagonal term

          INODE = KXX(NNUD,L)
          IIPOS = (NNUD-1)*NNUD + NNUD
          
          IF (IBCOD(INODE).EQ.1.AND. CAUDAL(INODE).GT.0) THEN
          
              DQDTRA(L,IIPOS) = DQDTRA(L,IIPOS) + EPST*DFLUDTRA(L,IIPOS)
          
         
              IF (ISOLFL.GT.1 .AND. IODENS.EQ.1) THEN
                 
                  DQDTRA(L,IIPOS) = DQDTRA(L,IIPOS)
     &                            + CFLU(L,NNUD)*TINCINV*ATKFACTOR

              END IF !IODENS.EQ.1

          END IF !IBCOD(INODE).EQ.1 .AND. IODENS.EQ.1

C------------------------ So far, the derivatives are those 
C------------------------ of (rhoi*·Qi), NOT those of Qi.
C------------------------ Notice that: 
C------------------------ @(rhoi*Qi)/@cj = rhoi* @(Qi)/@cj since @(rhoi*)/@cj = 0
C------------------------ Then, @(Qi)/@cj = 1/rhoi* @(rhoi*Qi)/@cj 
C------------------------ Now, the rhoi* factor is removed.

          IF (IODENS.EQ.1) THEN

              DO I=1,NNUD

                  INODE = KXX(I,L)

                  IF (IBCOD(INODE).EQ.1.AND. CAUDAL(INODE).GT.0) THEN

                      DO J=1,NNUD

                          IJPOS = (I-1)*NNUD + J

                          CEXT = PARNP(INODE,4)
                          DENSEXT = DENS(DENSREF,BETAC,CEXT,CREF)

                          DQDTRA(L,IJPOS) = DQDTRA(L,IJPOS) / DENSEXT
            
                      END DO  !J=1,NNUD

                  END IF ! IBCOD(INODE).EQ.1.AND. CAUDAL(INODE).GT.0

              END DO !I=1,NNUD

          END IF !IODENS.EQ.1

      END DO !L=1,NUMEL


      END SUBROUTINE COMP_DQDTRA
C***********************************************************************
C***********************************************************************
C***********************************************************************
C***********************************************************************
C***********************************************************************
C***********************************************************************
      SUBROUTINE COMP_DQDFLU
     &          (AFLU     ,BETAC    ,CAUDAL   ,CREF     ,DENSREF
     &          ,DFLU     ,DFLUDFLU ,DNODALRH ,DQDFLU   ,EPSFLU
     &          ,HAUX1    ,IBCOD    ,IDIMAFLU ,IDIMDFLU ,IOATK
     &          ,IODENS   ,ISOLFL   ,KXX      ,LMXNDL   ,LNNDEL
     &          ,NPPNP    ,NUMEL    ,NUMNP    ,PARNP    ,TINC)

      IMPLICIT NONE

      INTEGER*4::IDIMAFLU ,IDIMDFLU ,IOATK    ,IODENS   ,ISOLFL   
     &          ,LMXNDL   ,NPPNP    ,NUMEL    ,NUMNP
     &          ,IIPOS    ,IJPOS    ,JIPOS    ,L        ,NNUD
     &          ,I        ,INODE    ,J        ,JNODE
     &          ,IA_IJ    ,JJPOS

      INTEGER*4:: IBCOD(NUMNP)  ,KXX(LMXNDL,NUMEL)
     &           ,LNNDEL(NUMEL)

      REAL*8::BETAC    ,CEXT,CREF,DENS,DENSEXT     ,DENSREF
     &       ,TINC,TINCINV,DERCAUD,DHEADDH,DQBDH,DALFDH,QB,ALFA
     &       ,EPSAUX,EPSAUX2,EPSFLU,ATKFACTOR


      REAL*8::AFLU(NUMEL,IDIMAFLU)  ,CAUDAL(NUMNP)
     &       ,DFLU(NUMEL,IDIMDFLU)  ,DFLUDFLU(NUMEL,LMXNDL*LMXNDL)
     &       ,DNODALRH(NUMNP,4)
     &       ,DQDFLU(NUMEL,LMXNDL*LMXNDL),HAUX1(NUMNP)
     &       ,PARNP(NUMNP,NPPNP)
     
      INTEGER*4,ALLOCATABLE::IOUSED(:)

C------------------------Inicialitation.
C------------------------ DQDFLU is computed by adding the nodal values
C------------------------ so, the initialization to zero is mandatory.      
      DQDFLU = 0.D0

      IF (TINC.NE.0D0) TINCINV = 1D0/TINC

C------------------------ Needed so that DQDFLU can be converted into
C------------------------ derivative respecto to h at k or k+1

      IF (IOATK.EQ.1) THEN

          ATKFACTOR = -1D0
	    EPSAUX = (1D0-EPSFLU)/EPSFLU
	    EPSAUX2 = (1d0-EPSFLU)

      ELSE

          ATKFACTOR = 1D0
	    EPSAUX = 1D0
	    EPSAUX2 = EPSFLU

      END IF !IOATK.EQ.1
      
C----------------------------------------------
C------------------------ Prescribed head -----
C----------------------------------------------

      DO L=1,NUMEL

          NNUD=LNNDEL(L) 

          DO I=1,NNUD-1

              INODE=KXX(I,L)

              IF (IBCOD(INODE).EQ.1 .AND. CAUDAL(INODE).GT.0) THEN

                  IIPOS = (I-1)*NNUD + I

                  DQDFLU(L,IIPOS) = DQDFLU(L,IIPOS)
     &                            + EPSAUX*DFLUDFLU(L,IIPOS)

C-------------Transient model

                  IF (ISOLFL.GT.1) THEN

                      DQDFLU(L,IIPOS) = DQDFLU(L,IIPOS)
     &                                + DFLU(L,I)*TINCINV*ATKFACTOR

                  END IF !(ISOLFL.GT.1)

              END IF !IBCOD(INODE).EQ.1 .AND. CAUDAL(INODE).GT.0

c-------------Nondiagonal terms

              DO J=I+1,NNUD

                  JNODE = KXX(J,L)
                  IA_IJ = (I-1)*NNUD + J - I*(I+1)/2
                  IJPOS = (I-1)*NNUD + J
                  JIPOS = (J-1)*NNUD + I
                  JJPOS = (J-1)*NNUD + J

                  
                  IF (IBCOD(INODE).EQ.1 .AND. CAUDAL(INODE).GT.0) THEN

                       DQDFLU(L,IJPOS) = DQDFLU(L,IJPOS) 
     &                                 + EPSAUX*DFLUDFLU(L,IJPOS)
     &                                 + EPSAUX2*AFLU(L,IA_IJ)

                       DQDFLU(L,IIPOS) = DQDFLU(L,IIPOS) 
     &                                 - EPSAUX2*AFLU(L,IA_IJ)
                  ENDIF

                  IF (IBCOD(JNODE).EQ.1 .AND. CAUDAL(JNODE).GT.0) THEN

                      DQDFLU(L,JIPOS) = DQDFLU(L,JIPOS)
     &                                + EPSAUX*DFLUDFLU(L,JIPOS)
     &                                + EPSAUX2*AFLU(L,IA_IJ)

                      DQDFLU(L,JJPOS) = DQDFLU(L,JJPOS) 
     &                                - EPSAUX2*AFLU(L,IA_IJ)
                  END IF


              END DO !J=I+1,NNUD


          END DO !I=1,NNUD-1

C------------------------ Last Diagonal term

          INODE=KXX(NNUD,L)
          IIPOS =(NNUD-1)*NNUD + NNUD
          
          IF (IBCOD(INODE).EQ.1 .AND. CAUDAL(INODE).GT.0) THEN

              DQDFLU(L,IIPOS) =DQDFLU(L,IIPOS)+EPSAUX2*DFLUDFLU(L,IIPOS)

              IF (ISOLFL.GT.1) THEN

                  DQDFLU(L,IIPOS) = DQDFLU(L,IIPOS)
     &                            + DFLU(L,NNUD)*TINCINV*ATKFACTOR
              ENDIF ! ISOLFL.GT.1

          END IF ! IBCOD(INODE).EQ.1 .AND. CAUDAL(INODE).GT.0

C------------------------ So far, the derivatives are those 
C------------------------ of (rhoi*·Qi), NOT those of Qi.
C------------------------ Notice that:
C------------------------ @(rho*Q)/@h = rhoi* @(Qi)/@hj since @(rhoi*)/@hj = 0
C------------------------ Then, @(Qi)/@hj = 1/rhoi* @(rhoi*Qi)/@hj 
C------------------------ Now, the rhoi* factor is removed.

          IF (IODENS.EQ.1) THEN

              DO I=1,NNUD

                  INODE = KXX(I,L)

                  IF (IBCOD(INODE).EQ.1.AND. CAUDAL(INODE).GT.0) THEN

                      DO J=1,NNUD

                          IJPOS = (I-1)*NNUD + J

                          CEXT = PARNP(INODE,4)
                          DENSEXT = DENS(DENSREF,BETAC,CEXT,CREF)

                          DQDFLU(L,IJPOS) = DQDFLU(L,IJPOS) / DENSEXT
            
                      END DO  !J=1,NNUD

                  END IF ! IBCOD(INODE).EQ.1.AND. CAUDAL(INODE).GT.0


              END DO !I=1,NNUD

          END IF !IODENS.EQ.1
    
      END DO !L=1,NUMEL

C----------------------------------------------
C------------------------ Leakage and flow ----
C----------------------------------------------

C------------------------ Now the derivatives related to
C------------------------ the other boundary conditions.

C------------------------ Since these boundary conditions are nodewise
C------------------------ and DQDFLU is an elementwise matrix, the
C------------------------ contribution to the derivative is stored in the
C------------------------ first element containing the node.
C------------------------ The IOUSED vector is used to remember if the
C------------------------ node was found in other element before.

      ALLOCATE(IOUSED(NUMNP))
      IOUSED = 0

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
      
          DO I=1,NNUD

              INODE = KXX(I,L)
                              
              IF (CAUDAL(INODE).GT.0 .AND. IBCOD(INODE).GT.1
     &           .AND. IOUSED(INODE).EQ.0) THEN

                  IOUSED(INODE) = 1

                  SELECT CASE (IBCOD(INODE))

                      CASE (2) !Flow

                          DERCAUD = DNODALRH(INODE,2) !Der. of prsc. flow wrt. h


                      CASE(3,4) !Leakage and flow

                          DHEADDH = DNODALRH(INODE,1) !Der. of prsc. head wrt. h
                          DQBDH = DNODALRH(INODE,2)   !Der. of prsc. flow wrt. h
                          DALFDH = DNODALRH(INODE,3)  !Der. of leak coeff.wrt. h
                          QB = DNODALRH(INODE,4)      !Presc. flow
                          ALFA = PARNP(INODE,3)
                      
                          DERCAUD = DALFDH*(PARNP(INODE,1)-HAUX1(INODE))
     &                            + ALFA*(DHEADDH - 1D0)


                  END SELECT !IBCOD(INODE)

                 IIPOS = (I-1)*NNUD + I
                 DQDFLU(L,IIPOS) = DQDFLU(L,IIPOS) + EPSAUX2*DERCAUD


              END IF !CAUDAL(I).GT.0 .AND. IBCOD(I).GT.1 .AND. IOUSED(I).EQ.0

            END DO !INODE=1,NNUD

      END DO !L=1,NUMEL
      
      DEALLOCATE(IOUSED)

      END SUBROUTINE COMP_DQDFLU
