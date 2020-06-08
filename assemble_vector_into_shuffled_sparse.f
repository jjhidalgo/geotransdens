      SUBROUTINE ASSEMBLE_VECTOR_NODE_INTO_SHUFFLED_SPARSE
     ;(FACTOR      ,IA_ROWS     ,LBLOCK_NR      ,NN        ,NB
     ;,A           ,A_DSC       ,IAD_D          ,IADN_D) 

*******************************************************************
*  PURPOSE
*  To assemble an element-wise stored matrix without diagonal into a sparse matrix,
*  the rows and columns of which have been shuffled to decrease bandwidth.
*
*   EXTERNAL VARIABLES, SCALARS
* NUMEL
* LNMXDL
* LBLOCK_NR: the block of the non-shufled matrix where the element-wise matrix
*            should go. 
*            if 1, the block of derivatives of flow eqn to head
*            if 2, the block of derivatives of flow eqn to mass fraction
*            if 3, the block of derivatives of tpt eqn to head
*            if 4, the block of derivatives of tpt eqn to mass fraction
* LNNUDP
* IDIM
******************************************************************

      IMPLICIT NONE

      ! EXTERNAL VARIABLES, SCALARS
      INTEGER*4 LBLOCK_NR,IA_ROWS,NN,NB
      REAL*8 FACTOR 

      ! EXTERNAL VARIABLES, ARRAYS
      INTEGER*4 IAD_D(NB,NN),IADN_D(NN)
      REAL*8 A(IA_ROWS,1), A_DSC(NB,NN) 

      ! INTERNAL VARIABLES, SCALARS
      INTEGER*4 IROW,ICOL, KIJ,INODE
      REAL*8 A_COMP


      DO INODE = 1,IA_ROWS        !Loop over components of vector

          A_COMP=A(INODE,1)
          A_COMP=A_COMP*FACTOR


          SELECT CASE(LBLOCK_NR)

              CASE (1)  !block dfludflu
                  IROW = (INODE-1)*2+1
                  ICOL = (INODE-1)*2+1

              CASE (2)  !block dfludtra
                  IROW = (INODE-1)*2+1
                  ICOL = INODE*2

              CASE (3)  !block dtradflu
                  IROW = INODE*2
                  ICOL = (INODE-1)*2+1

              CASE(4)   !block dtradtra
                  IROW = INODE*2
                  ICOL = INODE*2

          END SELECT

          CALL FIND( IROW,ICOL,KIJ,IAD_D,IADN_D,NB,NN)
          A_DSC(KIJ,IROW) =  A_DSC(KIJ,IROW) + A_COMP

      ENDDO 


      RETURN
      END
*********************************************************************
      SUBROUTINE ASSEMBLE_VECTOR_ELEM_INTO_SHUFFLED_SPARSE
     ;(FACTOR    ,IA_COLS  ,IA_ROWS       ,LBLOCK_NR        ,NN        
     ;,NB        ,A        ,A_DSC         ,IAD_D            ,IADN_D
     ;,LNNDEL    ,KXX) 

*******************************************************************
*  PURPOSE
*  To assemble an element-wise stored matrix without diagonal into a sparse matrix,
*  the rows and columns of which have been shuffled to decrease bandwidth.
*
*   EXTERNAL VARIABLES, SCALARS
* NUMEL
* LNMXDL
* LBLOCK_NR: the block of the non-shufled matrix where the element-wise matrix
*            should go. 
*            if 1, the block of derivatives of flow eqn to head
*            if 2, the block of derivatives of flow eqn to mass fraction
*            if 3, the block of derivatives of tpt eqn to head
*            if 4, the block of derivatives of tpt eqn to mass fraction
* LNNUDP
* IDIM
******************************************************************

      IMPLICIT NONE

      ! EXTERNAL VARIABLES, SCALARS
      INTEGER*4 LBLOCK_NR,IA_ROWS,IA_COLS,NB,NN


      REAL*8 FACTOR 

      ! EXTERNAL VARIABLES, ARRAYS
      INTEGER*4 KXX(IA_COLS,IA_ROWS), LNNDEL(IA_ROWS),IAD_D(NB,NN)
     ;          ,IADN_D(NN)
      REAL*8 A(IA_ROWS,IA_COLS), A_DSC(NB,NN) 

      ! INTERNAL VARIABLES, SCALARS
      INTEGER*4 IROW,ICOL, KIJ,I,INODE,L,NNUD,NUMEL
      REAL*8 A_COMP
 


      NUMEL = IA_ROWS
      DO L=1,NUMEL                                           !Loop over elements 

          NNUD = LNNDEL(L)

          DO I =1,NNUD                                !Loop over nodes of element
              INODE = KXX(I,L)

              A_COMP=A(L,I)
              A_COMP=A_COMP*FACTOR

           
              SELECT CASE(LBLOCK_NR)

                  CASE (1)  !block dfludflu
                      IROW = (INODE-1)*2+1
                      ICOL = (INODE-1)*2+1
                    
                  CASE (2)  !block dfludtra
                      IROW = (INODE-1)*2+1
                      ICOL = INODE*2

                  CASE (3)  !block dtradflu
                      IROW = INODE*2
                      ICOL = (INODE-1)*2+1

                  CASE(4)   !block dtradtra
                      IROW = INODE*2
                      ICOL = INODE*2
         
              END SELECT

              CALL FIND( IROW,ICOL,KIJ,IAD_D,IADN_D,NB,NN)
              A_DSC(KIJ,IROW) =  A_DSC(KIJ,IROW) + A_COMP

         ENDDO                                       !Loop over nodes of element
      ENDDO                                                  !loop over elements


      RETURN
      END
