       SUBROUTINE COMP_VOL_NOD 
     ;(LMXNDL   ,NUMEL    ,NUMNP    ,ACTH   ,AREA     ,KXX      
     ;,LDIM     ,LNNDEL   ,VOLNOD)  

********************************************************************************
*
* PURPOSE Compute the volume related to a given node
*
* DESCRIPTION  To do that, a loop over elements is done, identifying its nodes
*              and computing the proportional part of the element size 
*              (i.e. lenght for 1d elements, area for 2D and volume for 3D)
*              One must take care with 2D models when solving tpt. eq. and
*              and multiply the last entity time aquifer thickness, which is 
*              stored on vector BTRA (by now).
*
* EXTERNAL VARIABLES: ARRAYS
*
*  ACTH                   Aquifer thickness of every element. Cross section for 
*                         1-D elements, thickness for 2-D elements.             
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  VOLNOD                 Nodal volume
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LDIM                   Vector containing the dimension of each element
*  LNNDEL                 Number of nodes at every element                      
*
* EXTERNAL VARIABLES: SCALARS
*
*  LMXNDL                 Maximum number of nodes per element                   
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*
* INTERNAL VARIABLES: SCALARS
*
*  AR                                                                           
*  I                                                                            
*  L                                                                            
*  NI                                                                           
*  NNUD                   Number of nodes of the current element                
*  WIDE                                                                         
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
* HISTORY: First coding ?
*          Revision: AAR (Feb-00)
*
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION KXX(LMXNDL,NUMEL),VOLNOD(NUMNP),LNNDEL(NUMEL),
     ;      AREA(NUMEL),ACTH(NUMEL),LDIM(NUMEL)

C______________________________ Initializes vector VOLNOD

       DO I=1,NUMNP
         VOLNOD(I)=0D0
       END DO

C______________________________ Loop over mesh elements

       DO L=1,NUMEL
         NNUD=LNNDEL(L)
         AR=AREA(L)/NNUD
         IF (LDIM(L).LT.3) THEN
           VOLEL=AR*ACTH(L)
         ELSE
           VOLEL=AR
         END IF

C______________________________ Loop over element nodes

         DO I=1,NNUD
           NI=KXX(I,L)
           VOLNOD(NI)=VOLNOD(NI)+VOLEL
         END DO       ! Next node
       END DO         ! Next element
       RETURN
       END
