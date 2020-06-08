      SUBROUTINE ASSEMBLE_SYM_INTO_SPARSE
     &(A            ,A_DSC        ,FACTOR       ,IAD_S       ,IADN_S         
     &,KXX          ,IA_COLS      ,LMXNDL       ,LNNDEL      ,NB        
     &,NN           ,NUMEL)


***********************************************************************
*   PURPOSE
*   To add the content of a symetric matrix to an assembled matrix in
*   sparse storage. 
*
*   DESCRIPTION
*   Step 1: loop over elements
*   Step 2: loop over element entries
*   Step 3: placing element entries
*
*
***********************************************************************

      IMPLICIT NONE
* EXTERNAL VARIABLES: SCALARS
      INTEGER*4 NUMEL,IA_COLS,LMXNDL,NB,NN
      REAL*8 FACTOR

* EXTERNAL VARIABLES: ARRAYS
      INTEGER  IAD_S(NB,NN),IADN_S(NN),KXX(LMXNDL,NUMEL)
     ;        ,LNNDEL(NUMEL)
      REAL*8 A(NUMEL,IA_COLS),A_DSC(NB,NN)
     ;               

* INTERNAL VARIABLES: SCALARS
      REAL*8  VALUE                     
      INTEGER*4 J,L,IGNODE1,IGNODE2,ILNODE1,ILNODE2,ICOL
     ;          ,NNUD,IMAXPOINTER

* INTERNAL VARIABLES: ARRAYS
      INTEGER*4 IPOINTER(21,2)

      DATA IPOINTER(1,1) /1/, IPOINTER(1,2) /1/
     ;    ,IPOINTER(2,1) /2/, IPOINTER(2,2) /1/
     ;    ,IPOINTER(3,1) /2/, IPOINTER(3,2) /2/
     ;    ,IPOINTER(4,1) /3/, IPOINTER(4,2) /1/
     ;    ,IPOINTER(5,1) /3/, IPOINTER(5,2) /2/
     ;    ,IPOINTER(6,1) /3/, IPOINTER(6,2) /3/
     ;    ,IPOINTER(7,1) /4/, IPOINTER(7,2) /1/
     ;    ,IPOINTER(8,1) /4/, IPOINTER(8,2) /2/
     ;    ,IPOINTER(9,1) /4/, IPOINTER(9,2) /3/
     ;    ,IPOINTER(10,1) /4/, IPOINTER(10,2) /4/
     ;    ,IPOINTER(11,1) /5/, IPOINTER(11,2) /1/
     ;    ,IPOINTER(12,1) /5/, IPOINTER(12,2) /2/
     ;    ,IPOINTER(13,1) /5/, IPOINTER(13,2) /3/
     ;    ,IPOINTER(14,1) /5/, IPOINTER(14,2) /4/
     ;    ,IPOINTER(15,1) /5/, IPOINTER(15,2) /5/
     ;    ,IPOINTER(16,1) /6/, IPOINTER(16,2) /1/
     ;    ,IPOINTER(17,1) /6/, IPOINTER(17,2) /2/
     ;    ,IPOINTER(18,1) /6/, IPOINTER(18,2) /3/
     ;    ,IPOINTER(19,1) /6/, IPOINTER(19,2) /4/
     ;    ,IPOINTER(20,1) /6/, IPOINTER(20,2) /5/
     ;    ,IPOINTER(21,1) /6/, IPOINTER(21,2) /6/

      
C____________________________________Step 1: loop over elements
      DO L=1,NUMEL
         NNUD = LNNDEL(L)
         IMAXPOINTER = (NNUD*(NNUD+1))/2

C____________________________________Step 2: loop over element entries
         
         
         DO J=1,IMAXPOINTER
            
            !identifying nodes
            ILNODE1 = IPOINTER(J,1)    
            ILNODE2 = IPOINTER(J,2) 
            IGNODE1 = KXX(ILNODE1,L)
            IGNODE2 = KXX(ILNODE2,L)  
            VALUE=A(L,J)
            VALUE=VALUE*FACTOR

C____________________________________Step 3: placing element entries
            !placing element (gnode1,gnode2)
            CALL FIND(IGNODE1,IGNODE2,ICOL,IAD_S,IADN_S,NB,NN)
            A_DSC(ICOL,IGNODE1)=A_DSC(ICOL,IGNODE1)+VALUE
            
            
            !if it is not a diagonal element, place gnode2,gnode1
            IF (IGNODE1.NE.IGNODE2) THEN
               !placing element (gnode2,gnode1)
               CALL FIND(IGNODE2,IGNODE1,ICOL,IAD_S,IADN_S,NB,NN)
               A_DSC(ICOL,IGNODE1)=A_DSC(ICOL,IGNODE1)+VALUE
            ENDIF 

        ENDDO          ! loop over element entries
      ENDDO            ! loop over elements
          
      RETURN
      END