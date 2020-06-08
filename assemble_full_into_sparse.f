      SUBROUTINE ASSEMBLE_FULL_INTO_SPARSE
     ;(A           ,A_DSC       ,FACTOR     ,IA_COLS    ,IAD_S
     ;,IADN_S      ,KXX         ,LNNDEL     ,LMXNDL     ,NB       
     ;,NN          ,NUMEL)
 

************************************************************************
*   PURPOSE
*   To add the content of a nonsymetric matrix to an assembled matrix in
*   sparse storage.
*
*   DESCRIPTION
*   Step 1: initiate 3 do loops (over elements, over node-i and over node-j)
*   Step 2: Saving matrix entry (l,i,j)
*
*
***********************************************************************



      IMPLICIT NONE
* EXTERNAL VARIABLES: SCALARS
      INTEGER*4 NUMEL,IA_COLS,LMXNDL,NB,NN
      DOUBLE PRECISION FACTOR

* EXTERNAL VARIABLES: ARRAYS
      INTEGER*4 KXX(LMXNDL,NUMEL),LNNDEL(NUMEL),IAD_S(NB,NN),IADN_S(NN)
      DOUBLE PRECISION A(NUMEL,IA_COLS),A_DSC(NB,NN)
     ;            

* INTERNAL VARIABLES: SCALARS
      DOUBLE PRECISION  VALUE                     
      INTEGER*4 I,J,L,INODE,JNODE,ICOL,NNUD,POS

* INTERNAL VARIABLES: ARRAYS

    
C_________________________________Step1: initiate loops over node combinations of elements      
      DO L=1,NUMEL
         NNUD=LNNDEL(L)
         DO I=1,NNUD
            INODE=KXX(I,L)
            DO J=1,NNUD
               JNODE=KXX(J,L)
                   
               POS=(I-1)*NNUD+J

C__________________________________Step 2: Saving matrix entry (l,i,j)
               VALUE=A(L,POS)
               VALUE = VALUE * FACTOR

               CALL FIND(INODE,JNODE,ICOL,IAD_S,IADN_S,NB,NN)
               A_DSC(ICOL,INODE)=A_DSC(ICOL,INODE)+VALUE


            ENDDO
         ENDDO
      ENDDO

      RETURN
      END