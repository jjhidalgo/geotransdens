      SUBROUTINE ASSEMBLE_CROSSTERM_INTO_SHUFFLED_SPARSE
     ;(FACTOR    ,IA_COLS    ,IA_ROWS   ,NUMEL     ,LMXNDL     
     ;,LBLOCK_NR ,NB         ,NN        ,IAD_D     ,IADD_D
     ;,IADN_D    ,KXX        ,LNNDEL    ,A         ,A_DSC) 
      IMPLICIT NONE

 ! EXTERNAL VARIABLES, SCALARS
      INTEGER*4 IA_COLS,IA_ROWS,NUMEL,LMXNDL,LBLOCK_NR,NB,NN
      REAL*8 FACTOR 

! EXTERNAL VARIABLES, ARRAYS
      INTEGER*4 IAD_D(NB,NN),IADN_D(NN),IADD_D(NN)
     ;         ,KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)
      REAL*8 A(IA_ROWS,IA_COLS),A_DSC(NB,NN) 


! INTERNAL VARIABLES, SCALARS
      INTEGER*4 INODE, JNODE,L,I,J,NNUD,KIJ,KJI,KII
     ;         ,IJROW,JIROW,IJCOL,JICOL,IIROW
      REAL*8 ELEM_IJ,ELEM_JI,ELEM_II
      
! INTERNAL VARIABLES, ARRAYS
      REAL*8, DIMENSION(:), ALLOCATABLE:: ROWVECTOR


      !allocate rowvector
      ALLOCATE(ROWVECTOR(LMXNDL))

      !loop over elements
      DO L=1,NUMEL
        NNUD=LNNDEL(L)
 
        !loop over nodes of element
        !fill local vector with row and calculate missing element
        ROWVECTOR=0D0
        DO I=2,NNUD
            ROWVECTOR(1) = ROWVECTOR(1)-A(L,I-1)
            ROWVECTOR (I) = A(L,I-1)
        ENDDO

        !assemble matrix to nondiagonal terms
        DO I=1,NNUD-1
            INODE = KXX(I,L)
            DO J=I,NNUD
                JNODE = KXX(J,L)

                !search the values to be assembled at ij and ji
                ELEM_IJ = ROWVECTOR(J)*FACTOR
                ELEM_JI = ROWVECTOR(I)*FACTOR

                !searching rows and columns corresponding to elements ij and ji
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
                   CALL FIND(JIROW,JICOL,KJI,IAD_D,IADN_D,NB,NN)

                  !and assemble 
                  A_DSC(KIJ,IJROW)= A_DSC(KIJ,IJROW) + ELEM_IJ 
                  A_DSC(KJI,JIROW)= A_DSC(KJI,JIROW) + ELEM_JI

            ENDDO
        ENDDO

        !assemble matrix to diagonal terms
        DO I=1,NNUD
             ELEM_II= ROWVECTOR(I)*FACTOR
             INODE = KXX(I,L)

             SELECT CASE(LBLOCK_NR)

              CASE (1)  !block dfludflu
                  IIROW = (INODE-1)*2+1

              CASE (2)  !block dfludtra
                  IIROW = (INODE-1)*2+1

              CASE (3)  !block dtradflu
                  IIROW = INODE*2

              CASE(4)   !block dtradtra
                  IIROW = INODE*2
           END SELECT

           KII = IADD_D(IIROW)
           A_DSC(KII,IIROW)= A_DSC(KII,INODE)+ ELEM_II

        ENDDO

      ENDDO !loop over elements
      
      !Finally destroy local array
      DEALLOCATE(ROWVECTOR)

      RETURN
      END
