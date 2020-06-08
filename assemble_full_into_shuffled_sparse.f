      SUBROUTINE ASSEMBLE_FULL_INTO_SHUFFLED_SPARSE
     ;(FACTOR    ,IA_COLS    ,IA_ROWS   ,NUMEL     ,LMXNDL     
     ;,LBLOCK_NR ,NB         ,NN        ,IAD_D     ,IADN_D     
     ;,KXX       ,LNNDEL     ,A         ,A_DSC) 

*******************************************************************
*  PURPOSE
*  To assemble an element-wise stored matrix into a sparse matrix,
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
      INTEGER*4 IA_COLS,IA_ROWS,NUMEL,LMXNDL,LBLOCK_NR,NB,NN
      REAL*8 FACTOR 

! EXTERNAL VARIABLES, ARRAYS
      INTEGER*4 IAD_D(NB,NN),IADN_D(NN)
     ;         ,KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)
      REAL*8 A(IA_ROWS,IA_COLS),A_DSC(NB,NN) 

! INTERNAL VARIABLES, SCALARS
      INTEGER*4 IROW,ICOL, KIJ, KJI,L,I,J,NNUD,INODE,JNODE,LPOS_A,KII
     ;          ,IJROW,IJCOL,JIROW,JICOL
      REAL*8   A_ELEM_IJ,A_ELEM_JI,A_ELEM_II  





      DO L=1,NUMEL                         !loop over elements
          NNUD=LNNDEL(L)
          DO I=1,NNUD-1                    !loop over nodes of element
              INODE=KXX(I,L)
              DO J= I+1,NNUD               !loop over nodes of element
                  JNODE= KXX(J,L)

                  !retrieving matrix-element ij
                  LPOS_A=(I-1)*NNUD+J
                  A_ELEM_IJ= A(L,LPOS_A)

                  !retrieving matrix-element ji
                  LPOS_A=(J-1)*NNUD+I
                  A_ELEM_JI= A(L,LPOS_A)


                  A_ELEM_IJ=A_ELEM_IJ*FACTOR
                  A_ELEM_JI=A_ELEM_JI*FACTOR


                 SELECT CASE(LBLOCK_NR)

                     CASE (1)  !block dfludflu
                         IJROW = (INODE-1)*2+1
                         IJCOL = (JNODE-1)*2+1
                         JIROW = (JNODE-1)*2+1
                         JICOL = (INODE-1)*2+1

                     CASE (2)  !block dfludtra
                         IJROW = (INODE-1)*2+1
                         IJCOL = JNODE*2
                         JIROW = (JNODE-1)*2+1
                         JICOL = INODE*2


                     CASE (3)  !block dtradflu
                         IJROW = INODE*2
                         IJCOL = (JNODE-1)*2+1
                         JIROW = JNODE*2
                         JICOL = (INODE-1)*2+1


                     CASE(4)   !block dtradtra
                         IJROW = INODE*2
                         IJCOL = JNODE*2
                         JIROW = JNODE*2
                         JICOL = INODE*2


                    END SELECT

                   CALL FIND(IJROW,IJCOL,KIJ,IAD_D,IADN_D,NB,NN)
                   A_DSC(KIJ,IJROW) =  A_DSC(KIJ,IJROW) + A_ELEM_IJ

                   CALL FIND(JIROW,JICOL,KJI,IAD_D,IADN_D,NB,NN)
                   A_DSC(KJI,JIROW) =  A_DSC(KJI,JIROW) + A_ELEM_JI

               ENDDO     !J
          ENDDO   !I

c___________________________Diagonal elements
          DO I = 1,NNUD

              INODE =KXX(I,L)
              !retrieving matrix-element ij
              LPOS_A=(I-1)*NNUD+I
              A_ELEM_II= A(L,LPOS_A)

              A_ELEM_II=A_ELEM_II*FACTOR


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

            CALL FIND(IROW,ICOL,KII,IAD_D,IADN_D,NB,NN)
            A_DSC(KII,IROW) =  A_DSC(KII,IROW) + A_ELEM_II

        ENDDO !I

      ENDDO  !L


      RETURN
      END
