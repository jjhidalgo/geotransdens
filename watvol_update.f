      SUBROUTINE WATVOL_UPDATE
     &          (HCAL     ,HCALAN   ,IOVRWC   ,KXX
     &          ,LMXNDL   ,LNNDEL   ,NPBTP    ,NPPEL    ,NUMEL
     &          ,NUMNP    ,PAREL    ,WATVOL)
      
********************************************************************************
*
* PURPOSE
*
*  Updates WATVOL for the next time step.
*
*
* DESCRIPTION
*
*  Updates WATVOL for the next time step.
*  at the end of time k+1, WATVOL(#,#,1) must store the value of WATVOL in k+1
*  so that it can be used as W in k for the next iteration.
*
*     End of time step k+1:
*
*             WATVOL(#,#,1) = WATVOL(in k)
*
*     Update:
*
*             WATVOL(#,#,1) = WATVOL(in k+1)
*
*     Begining of time k+2
*
*             WATVOL(#,#,1) = WATVOL(in k+1) (i.e. previous time step)
*
*
* EXTERNAL VARIABLES: ARRAYS
*
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).
*  LNNDEL                 Number of nodes at every element
*  ACTH                   Aquifer thickness of every element. Cross sectional   
*                         area for 1-D elements, thickness for 2-D elements.
*  HCAL                   Computed heads at every node                          
*  HCALAN                 Head level at previous time
*  PAREL                  Parameter values at every element and current time    
*                         for all nodal parameters (each value is computed as   
*                         the product of up to four terms:                      
*                           elem. coeff*zonal value*time funct.*nonl. funct. )  
*  WATVOL                 Array containing the water content of every element
*                         The array stores the water content in time k and
*                         k+theta.
*                             WATVOL(#,#,1) --> Time k.
*                             WATVOL(#,#,2) --> Time k+theta
*
*
* EXTERNAL VARIABLES: SCALARS
*
*                         number INTI and observation time number INTI+1        
*  LMXNDL                 Maximum number of nodes per element                   
*  NUMEL                  Number of elements                                    
*  IOVRWC                 Equal to IOPTS(31).
*                             1.WATVOL calculated elementwise
*                             2.WATVOL calculated nodewise.
*  NUMNP                  Number of nodes                                       
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
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

      INTEGER::LMXNDL, NUMEL, IOVRWC,NPBTP, NPPEL,  NUMNP


      INTEGER::KXX(LMXNDL,NUMEL),
     &         LNNDEL(NUMEL)


      REAL*8::HCAL(NUMNP),
     &        HCALAN(NUMNP),
     &        PAREL(NUMEL,NPPEL),
     &        WATVOL(MAX(1,(IOVRWC-1)*LMXNDL),NUMEL,3,NPBTP)


C------------------------- Internal

      INTEGER*4::I,J,L,NNUD

      REAL*8::H_AVG

C------------------------- First executable estatment


**********
C________WATVOL EN K+THETA (antes en COMP_ELEM_VAR)
**********

C------------------------- Water volume 

C-------------------------  IF (IOVRWC.EQ.0) THEN Constant in time

      IF (IOVRWC.EQ.1) THEN

C------------------------- Variable in time by elements

C------------------------- Head average increment at current element

          DO L=1,NUMEL

              H_AVG = 0D0
              NNUD = LNNDEL(L)
              DO J=1,NNUD
                  I = KXX(J,L)
                  H_AVG = H_AVG + (HCAL(I)-HCALAN(I))
              END DO !J=1,NNUD

C------------------------- Water volume of the element is updated

              WATVOL(1,L,1,1:NPBTP) = WATVOL(1,L,1,1:NPBTP)
     &                               + PAREL(L,7)*H_AVG/NNUD
          END DO !L=1,NUMEL

      ELSE IF (IOVRWC.EQ.2) THEN

C------------------------- Variable in time by nodes

          DO L=1,NUMEL

C------------------------- Water volume of the element is updated

              NNUD = LNNDEL(L)

                  DO J=1,NNUD
                      I = KXX(J,L)
                      WATVOL(J,L,1,1:NPBTP) = WATVOL(J,L,1,1:NPBTP)
     &                                + PAREL(L,7)*(HCAL(I)-HCALAN(I))
                      END DO !J=1,NNUD

          END DO !L=1,NUMEL
      END IF !IOVRC.EQ.1,2

      END SUBROUTINE WATVOL_UPDATE
