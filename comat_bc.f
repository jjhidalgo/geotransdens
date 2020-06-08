      SUBROUTINE COMAT_BC 
     &          (ADSC     ,BETAC    ,CAUDAL   ,CAUX1    ,CREF
     &          ,DENSREF  ,DQDFLU   ,DQDTRA   ,EPSTRA   ,IAD
     &          ,IADD     ,IADN     ,IBTCO    ,IDSC_COLS,IDSC_ROWS
     &          ,IOCOUPLED,IOASSEMB ,IOATK    ,IODENS   ,IONEWT
     &          ,ITYPADSC ,KXX      ,LMXNDL   ,LNNDEL   ,MAXNB
     &          ,MAXNN    ,NBAND    ,NPPNP    ,NUMEL    ,NUMNP
     &          ,PARNP    ,THETAT,AREA,NZONE_PAR,PAREL,NPPEL,NTYPAR)

********************************************************************************
*
* PURPOSE
*
*     Set in ADSC mass flow boundary conditions
*
* DESCRIPTION
*
*     Set in ADSC mass flow boundary conditions. When CAUDAL(I)>0, it is 
*     added to the diagonal of ADSC.
*
*
*      CAUDAUX = RHO*Q, se ensambla multiplicado por THETAT.
*
*          Ensambled to diagonal of ADSC when Picard's method is used.
*
*
*      CAUDAUX2 = d/dh Q*RHOext*(Wext - Wk+th) =
*               = dQdh*THETA*RHOext*(Wext - Wk+th)
*
*      CAUDAUX3  = d/dw Q*RHOext*(Wext - Wk+th)
*                = dQdw*RHOext*(Wext - Wk+th) - THETA*Q*RHOext
*
*          Emsambled to the blocks 3 and 4 of ADSC when Newton's method is used
*          multiplied by -EPSILON.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  ATRA                   Matrix of finite elements equations for transport     
*                         problem. No boundary conditions are included.         
*  CAUDAL                 Input/output flow at every node.                      
*  IBTCO                  Transport boundary condition index                    
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IOCOUPLED              1 If solving a the coupled density-dependent system
*                           This implies the use on Newton's method. Then vectors
*                           IAD,IADD and IADN are IAD_D, IADD_D and IADN_D
*  NUMNP                  Number of nodes                                       
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ADD_TO_DIAG
*
* HISTORY
*
*     AMS        1988     First coding
*     AMS      3-1999     Revision
*     JHG      8-2003     Density dependent flow
*                         Elementwise storage of ATRA.
*
********************************************************************************

      IMPLICIT NONE

      INTEGER*4 IB,ITYPADSC,IDSC_COLS,IDSC_ROWS,LMXNDL,MAXNB
     &          ,MAXNN,NUMEL,NUMNP,NBAND,IONEWT,I,NPPNP,L,NNUD
     &          ,J,IODENS,INODE,IOCOUPLED,IOASSEMB,IOATK,NPPEL,NTYPAR

      INTEGER*4 IBTCO(NUMNP),IAD(MAXNB,MAXNN),IADD(MAXNN),IADN(MAXNN),
     &           KXX(LMXNDL,NUMEL),LNNDEL(NUMEL),NZONE_PAR(NTYPAR)

      REAL*8 THETAT,DENSEXT,EPSTRA
     &       ,DENS,BETAC,DENSREF,CREF,CEXT


      REAL*8::AREA(NUMEL),ADSC(IDSC_ROWS,IDSC_COLS),CAUDAL(NUMNP)
     &       ,CAUX1(NUMNP),
     &        DQDTRA(NUMEL,LMXNDL*LMXNDL),DQDFLU(NUMEL,LMXNDL*LMXNDL),
     &        PAREL(NUMEL,NPPEL),PARNP(NUMNP,NPPNP)

      REAL*8::AREALN,CAUD,CNOD,CREC,DENSREC,EPST,FACT,FACT2,RECH
      INTEGER*4::II_POS,IJ_POS

      INTEGER*4,ALLOCATABLE::IOUSED(:)
      REAL*8,ALLOCATABLE::CAUDAUX(:),CAUDAUX3(:,:),CAUDAUX2(:,:)



C------------------------- If Picard method is used, THETA*RHOext*Qi
C------------------------- must be added to the diagonal of ADSC

      IF (IONEWT.EQ.0) THEN

          ALLOCATE(CAUDAUX(NUMNP))
          CAUDAUX=0D0

          DO I=1,NUMNP
      
              IB = IBTCO(I)
              CAUD = CAUDAL(I)

              IF (IB.EQ.2 .OR. IB.EQ.3) THEN

                  IF (CAUD.GT.0) THEN

                      CEXT = PARNP(I,4)

                      IF (IODENS.EQ.1) THEN
                          DENSEXT = DENS(DENSREF,BETAC,CEXT,CREF)
                      ELSE
                          DENSEXT = 1D0
                      END IF !IODENS.EQ.1

                      CAUDAUX(I) = DENSEXT*CAUD

                  ELSE

                      CAUDAUX(I) = 0D0

                  END IF !CAUDAL(I).GT.0

              END IF !IB.EQ.2 .OR. IB.EQ.3

          END DO !I=1,NUMNP

          IF (IOCOUPLED.EQ.1) THEN

              CALL ASSEMBLE_SHUFFLED
     &(THETAT       ,1               ,ITYPADSC  ,1         ,NUMNP   
     &,IDSC_COLS    ,IDSC_ROWS       ,4             ,LMXNDL    ,NUMEL     
     &,MAXNB        ,MAXNN           ,CAUDAUX       ,ADSC      ,IAD           
     &,IADD         ,IADN            ,KXX           ,LNNDEL    ,NBAND)
     
          ELSE

              CALL ASSEMBLE
     &(THETAT    ,1            ,NUMNP     ,IDSC_COLS  ,IDSC_ROWS
     &,1         ,ITYPADSC     ,LMXNDL    ,MAXNB      ,NBAND
     &,MAXNN     ,NUMEL        ,CAUDAUX   ,ADSC       ,IAD
     &,IADD      ,IADN         ,KXX       ,LNNDEL)

          END IF !IOCOUPLED.EQ.1

C------------------------- Memory used by CAUDAUX is released.

          DEALLOCATE(CAUDAUX)

      END IF !IONEWT.EQ.0

C------------------------- The following part is only for Newton

      IF (IONEWT.EQ.1) THEN

          IF (IOATK.EQ.1) THEN

              EPST = (1D0 - EPSTRA)

          ELSE

              EPST = EPSTRA

          END IF !IOATK.EQ.1

          ALLOCATE(CAUDAUX3(NUMEL,LMXNDL*LMXNDL))
          CAUDAUX3 = 0D0

          IF (IOCOUPLED.EQ.1) THEN
              ALLOCATE(CAUDAUX2(NUMEL,LMXNDL*LMXNDL))
              CAUDAUX2 = 0D0
          END IF !IOCOUPLED.EQ.1


          ALLOCATE (IOUSED(NUMNP))
          IOUSED = 0

          DO L=1,NUMEL
              
              NNUD=LNNDEL(L)
              AREALN = AREA(L)/NNUD
	        RECH = PAREL(L,8)*AREALN

            IF (RECH.GT.0) THEN

                CREC = PAREL(L,15)

            ELSE

                CREC = 0D0

            END IF !RECH.GT.0

            IF (IODENS.EQ.1) THEN

                  DENSREC = DENS(DENSREF,BETAC,CREC,CREF)

            ELSE

                DENSREC = 1D0

            END IF !IODENS.EQ.1

              DO I=1,NNUD

                  INODE = KXX(I,L)
                  IB = IBTCO(INODE)
                  CAUD = CAUDAL(INODE)
                  CNOD = CAUX1(INODE)
                  II_POS = (I-1)*NNUD + I

C------------------------- Derivative of recharge contribution w. r. t.
C------------------------- concentration is always added

                  IF (NZONE_PAR(3).NE.0 .AND. IOCOUPLED.EQ.1) THEN

                      FACT = EPST*RECH*DENSREC
                      CAUDAUX3(L,II_POS) = CAUDAUX3(L,II_POS) - FACT

                  END IF !NZONE_PAR(3).NE.0


C------------------------- Only nodes with mass flow boundary condition
C------------------------- and incoming flow account for this contribution.

                  IF (CAUD.GT.0 .AND. (IB.EQ.2 .OR. IB.EQ.3)) THEN

                      CEXT = PARNP(INODE,4)

                      IF (IODENS.EQ.1) THEN

                          DENSEXT = DENS(DENSREF,BETAC,CEXT,CREF)

                      ELSE

                          DENSEXT = 1D0

                      END IF !IODENS.EQ.1


                      IF (IOCOUPLED.EQ.1 .AND. IOUSED(INODE).LT.1) THEN

                          FACT = EPST*CAUD*DENSEXT

                          CAUDAUX3(L,II_POS) = CAUDAUX3(L,II_POS)
     &                                       - FACT

C------------------------- Avoids repetition of the contribution to this node
C------------------------- during elements loop.

                          IOUSED(INODE) = 1 
                                            
                      END IF !IOCOUPLED.EQ.1

                      FACT2 = DENSEXT*(CEXT - CNOD)

                      DO J=1,NNUD

                          IJ_POS = (I-1)*NNUD + J

                          IF (IOCOUPLED.EQ.1) THEN

                              CAUDAUX2(L,IJ_POS) = CAUDAUX2(L,IJ_POS) 
     &                                          + DQDFLU(L,IJ_POS)*FACT2
                          END IF !IOCOUPLED.EQ.1

                          CAUDAUX3(L,IJ_POS) = CAUDAUX3(L,IJ_POS)
     &                                       + DQDTRA(L,IJ_POS)*FACT2

                      END DO !J=1,NNUD

                  END IF !CAUDAL(I).GT.0 .AND. (IB.EQ.2 .OR. IB.EQ.3)
              END DO !I=1,NNUD
          
          END DO ! L=1,NUMEL

C------------------------- Memory used by IOUSED is released.

          DEALLOCATE(IOUSED)

          IF (IOASSEMB.EQ.1) THEN

              IF (IOCOUPLED.EQ.1) THEN

C------------------------- Contibution to block 4 (dTdT)

                  CALL ASSEMBLE_SHUFFLED
     &(-1D0         ,4               ,ITYPADSC      ,LMXNDL*LMXNDL,NUMEL
     &,IDSC_COLS    ,IDSC_ROWS       ,4             ,LMXNDL       ,NUMEL     
     &,MAXNB        ,MAXNN           ,CAUDAUX3      ,ADSC        ,IAD           
     &,IADD         ,IADN            ,KXX           ,LNNDEL      ,NBAND)

C------------------------- Contibution to block 3 (dTdF)

                  CALL ASSEMBLE_SHUFFLED
     &(-1D0         ,4               ,ITYPADSC      ,LMXNDL*LMXNDL,NUMEL
     &,IDSC_COLS    ,IDSC_ROWS       ,3             ,LMXNDL       ,NUMEL     
     &,MAXNB        ,MAXNN           ,CAUDAUX2      ,ADSC        ,IAD           
     &,IADD         ,IADN            ,KXX           ,LNNDEL      ,NBAND)

     
C------------------------- Memory used by CAUDAUX2 is released.

                  DEALLOCATE(CAUDAUX2)

              ELSE

                  CALL ASSEMBLE
     ;(-1D0      ,LMXNDL*LMXNDL,NUMEL     ,IDSC_COLS  ,IDSC_ROWS  
     ;,4         ,ITYPADSC     ,LMXNDL    ,MAXNB      ,NBAND  
     ;,NUMNP     ,NUMEL        ,CAUDAUX3  ,ADSC   ,IAD      
     ;,IADD      ,IADN         ,KXX       ,LNNDEL)


              END IF !IOCOUPLED.EQ.1

          ELSE

              
              DQDTRA(:,:) = CAUDAUX3(:,:)

              IF (IOCOUPLED.EQ.1) THEN

                  DQDFLU(:,:) = CAUDAUX2(:,:)

                  DEALLOCATE(CAUDAUX2)

              END IF !IOCOUPLED.EQ.1

          END IF !IOASSEMB.EQ.1

C------------------------- Memory used by CAUDAUX3 is released.

          DEALLOCATE(CAUDAUX3)

      END IF !IONEWT.EQ.1

      
      END SUBROUTINE COMAT_BC 

