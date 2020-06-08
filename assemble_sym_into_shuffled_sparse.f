      SUBROUTINE ASSEMBLE_SYM_INTO_SHUFFLED_SPARSE
     ;(FACTOR     ,IA_COLS   ,IA_ROWS    ,LBLOCK_NR  ,LMXNDL    
     ;,NB         ,NN        ,NUMEL      ,A          ,A_DSC     
     ;,IAD_D      ,IADN_D    ,KXX        ,LNNDEL)
          
*******************************************************************
*  PURPOSE
*  To assemble an element-wise stored symetric matrix into a sparse matrix,
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
* IDIM
******************************************************************

      
      IMPLICIT NONE

! EXTERNAL VARIABLES, SCALARS
      INTEGER*4 NUMEL,LMXNDL,LBLOCK_NR,IA_COLS,IA_ROWS,NB,NN    
      REAL*8 FACTOR

! EXTERNAL VARIABLES, ARRAYS
      INTEGER*4 KXX(LMXNDL,NUMEL),LNNDEL(NUMEL),IAD_D(NB,NN)
     ;          ,IADN_D(NN)
      REAL*8 A(IA_ROWS,IA_COLS),A_DSC(NB,NN) 

! INTERNAL VARIABLES, SCALARS
      INTEGER*4 IROW,ICOL, KIJ, KJI, L,NNUD,I,J,INODE
     ;          ,JNODE,LPOS_A
      REAL*8   A_ELEM





      DO L=1,NUMEL                         !loop over elements
          NNUD=LNNDEL(L)
          DO I=1,NNUD                      !loop over nodes of element
              INODE=KXX(I,L)
              DO J= I,NNUD                 !loop over nodes of element
                  JNODE= KXX(J,L)

                  !retrieving matrix-element ij
                  LPOS_A=(I-1)*NNUD+J
                  A_ELEM= A(L,LPOS_A)
                  A_ELEM=A_ELEM*FACTOR


                  SELECT CASE(LBLOCK_NR)

                      CASE (1)  !block dfludflu
                          IROW = (INODE-1)*2+1
                          ICOL = (JNODE-1)*2+1

                      CASE (2)  !block dfludtra
                          IROW = (INODE-1)*2+1
                          ICOL = JNODE*2

                      CASE (3)  !block dtradflu
                          IROW = INODE*2
                          ICOL = (JNODE-1)*2+1

                      CASE(4)   !block dtradtra
                          IROW = INODE*2
                          ICOL = JNODE*2
                  END SELECT

                  CALL FIND( IROW,ICOL,KIJ,IAD_D,IADN_D,NB,NN)
                  A_DSC(IROW,KIJ) =  A_DSC(IROW,KIJ) + A_ELEM

                  IF (J.NE.I) THEN
                      CALL FIND( ICOL,IROW,KJI,IAD_D,IADN_D,NB,NN)
                      A_DSC(ICOL,KJI) =  A_DSC(IROW,KJI) + A_ELEM
                  ENDIF

              ENDDO     !J
          ENDDO   !I
      ENDDO  !L

      RETURN
      END
