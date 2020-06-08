      SUBROUTINE CALC_WATVOL
     &          (HCALAN   ,HCALIT   ,INTI     ,IOVRWC   ,KXX
     &          ,LMXNDL   ,LNNDEL   ,NPBTP    ,NPPEL    ,NUMEL
     &          ,NUMNP    ,PAREL    ,THETAT   ,WATVOL)

********************************************************************************
*
* PURPOSE
*
*  Computes volume of water per unit of aquifer.
*
*
* DESCRIPTION
*
*  Manages de the computation of WATVOL vector in times k+1 and k+theta
*
*
*
* EXTERNAL VARIABLES: ARRAYS
*
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).
*  LNNDEL                 Number of nodes at every element
*  LXPAREL                Array containing zone numbers for a given             
*                         element parameter
*  CFPAREL                Array containing node coefinient of element j         
*                         corresponding to INpar index parameter zone.          
*  HCALIT                 Computed heads at every node                          
*  HCALAN                 Head level at previous time
*  PARC                   Vector containing calculated values for all           
*                         parameters
*  PAR_DIR                Array containing all real direct problem              
*                         parameters                                            
*  PAREL                  Parameter values at every element and current time    
*                         for all nodal parameters (each value is computed as   
*                         the product of up to four terms:                      
*                           elem. coeff*zonal value*time funct.*nonl. funct. )  
*  WATVOL                 Array containing the water content of every element
*                         The array stores the water content in time k and
*                         k+theta.
*                             WATVOL(#,#,1) --> Time k.
*                             WATVOL(#,#,2) --> Time k+theta
*                             WATVOL(#,#,3) --> Time k+1-theta
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  LMXNDL                 Maximum number of nodes per element                   
*  NUMEL                  Number of elements                                    
*  IOVRWC                 Equal to IOPTS(31).
*                             1.WATVOL calculated elementwise
*                             2.WATVOL calculated nodewise.
*  NUMNP                  Number of nodes                                       
*  NPAREL                 Number of element parameters in current problem       
*  NPARALG                Maximum number of algorithm parameters, including     
*                         minimization and simulation (used for dimensioning)   
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  THETAT                 Time weighting parameter for transport problems       
*



* INTERNAL VARIABLES: SCALARS
*
*  H_AVG
*  I
*  J
*  K
*  L                      Current element
*  NNUD                   Number of nodes of the current element                
*
*

* HISTORY
*
*     JHG      6-2003     First coding.
*                 
*******************************************************************************
      IMPLICIT NONE
      
C------------------------- External       

      INTEGER::INTI, LMXNDL, NUMEL, IOVRWC,NPBTP,NPPEL, NUMNP


      INTEGER::KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)

      REAL*8::THETAT


      REAL*8::HCALAN(NUMNP) ,HCALIT(NUMNP) ,PAREL(NUMEL,NPPEL)
     &       ,WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3,NPBTP)


C------------------------- Internal

      INTEGER*4::I,J,L,NNUD
      REAL*8::H_AVG, THT1

C------------------------- First executable estatment

      IF (INTI.NE.0) THEN


          THT1 = 1D0 -THETAT

**********
C________WATVOL EN K+THETA (antes en COMP_ELEM_VAR)
**********

C------------------------- Details of computation
C------------------------- 
C------------------------- W(k+theta) = (1-theta)*W(k) + theta*W(k+1)
C-------------------------
C------------------------- W(k+theta) = (1-theta)*W(k) + theta*(W(k) + DeltaW)
C-------------------------
C------------------------- Simplifying
C-------------------------
C-------------------------  W(k+theta) = W(k) + theta*(DeltaW)
C-------------------------


C------------------------- Water volume 

C------------------------- IF (IOVRWC.EQ.0) THEN Constant in time

          IF (IOVRWC.EQ.1) THEN

C------------------------- Variable in time by elements

C------------------------- Head average increment at current element

              DO L=1,NUMEL

                  H_AVG = 0D0
                  NNUD = LNNDEL(L)

                  DO J=1,NNUD
                      I = KXX(J,L)
                      H_AVG = H_AVG + (HCALIT(I)-HCALAN(I))
                  END DO !J=1,NNUD

C------------------------- Water volume of the element is updated

                  WATVOL(1,L,2,1:NPBTP) = WATVOL(1,L,1,1:NPBTP)
     &                                   + PAREL(L,7)*H_AVG*THETAT/NNUD

                  WATVOL(1,L,2,1:NPBTP) = WATVOL(1,L,1,1:NPBTP)
     &                                   + PAREL(L,7)*H_AVG*THT1/NNUD
              END DO !L=1,NUMEL

          ELSE IF (IOVRWC.EQ.2) THEN

C------------------------- Variable in time by nodes

              DO L=1,NUMEL

C------------------------- Water volume of the element is updated

                  NNUD = LNNDEL(L)

                      DO J=1,NNUD

                          I = KXX(J,L)
                          WATVOL(J,L,2,1:NPBTP) = WATVOL(J,L,1,1:NPBTP)
     &                         + THETAT*PAREL(L,7)*(HCALIT(I)-HCALAN(I))

                          WATVOL(J,L,3,1:NPBTP) = WATVOL(J,L,1,1:NPBTP)
     &                           + THT1*PAREL(L,7)*(HCALIT(I)-HCALAN(I))

                      END DO !J=1,NNUD

              END DO !L=1,NUMEL

          END IF !IOVRC.EQ.1,2

      END IF !INTI.NE.0

      END SUBROUTINE CALC_WATVOL
