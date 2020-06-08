      SUBROUTINE ASSEMBLE_CROSSTERM_INTO_SHUFFLED_BAND_NO_SYM
     ;          (FACTOR    ,IA_COLS   ,IA_ROWS    ,IDSC_COLS  ,IDSC_ROWS
     ;          ,LMXNDL    ,NBAND     ,LBLOCK_NR  ,A          ,A_DSC
     ;          ,LNNDEL    ,KXX)


********************************************************************************
*
* PURPOSE  To assemble a matrix whose rows are equal and the first element is
*          equal to the sum of the rest of th elements with the opposite sing
*          into a shuffled non-symmetric banded matrix.
*
*
* DESCRIPTION The subroutine assembles the elements (multiplied by a factor)
*            of a 'row-wise' matrix into shuffled symmetric banded one.
*
*
*        The 'row-wise' matrix is a matrix whiose rows are all equal.
*             Furthermore, the elemetn A1J = -SUMJ(AIJ), J>2.
*             It must be stored in a vector in the following
*             fashion:
*
*                     | 11  12  13 |
*                     | 11  12  13 |   =>  |12, 13|
*                     | 11  12  13 |
*
*            Then, The element (I,J) is located the position
*            (J - 1) in the vector (I > 1)
*
* EXTERNAL VARIABLES: ARRAYS
*
*  A                        Origin matrix (full)
*  A_DSC                    Destiny matrix (non-symmetric banded).
*  LBLOCK_NR                The block of the non-shufled matrix where the element-wise matrix
*                        should go.
*                            if 1, the block of derivatives of flow eqn to head
*                            if 2, the block of derivatives of flow eqn to mass fraction
*                            if 3, the block of derivatives of tpt eqn to head
*                            if 4, the block of derivatives of tpt eqn to mass fraction
*
* EXTERNAL VARIABLES: SCALARS
*
*  FACTOR                    Factor that multiplies the matrix to be assembled.
*                        For time-weighted schemes pourposes (theta-scheme).
*
*  IA_COLS                Number of columns (first dimension) of A matrix.
*  IA_ROWS                Number of rows (second dimension) of A matrix.
*  IDSC_COLS                Number of columns (first dimension) of A_DSC matrix.
*  IDSC_ROWS                Number of rowss (second dimension) of A_DSC matrix.
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
* HISTORY: First coding: JHG (11-2004)
*
********************************************************************************



        IMPLICIT NONE

C-------------------- External

      INTEGER::IA_COLS  ,IA_ROWS  ,IDSC_COLS,IDSC_ROWS,LMXNDL
     &        ,NBAND    ,LBLOCK_NR

      Real*8::FACTOR

      INTEGER*4::LNNDEL(IA_ROWS),KXX(LMXNDL,IA_ROWS)

      REAL*8::A(IA_ROWS,IA_COLS), A_DSC(IDSC_ROWS,IDSC_COLS)

C-------------------- Internal

      INTEGER::I        ,IDIAG    ,IJ_COL   ,IJ_ROW   ,J
     &        ,JI_COL   ,JI_ROW   ,KEXT     ,KEXT_IND ,KEXTBAND
     &        ,KINT     ,KINT_IND ,KINTBAND ,L        ,LMDIAG
     &        ,NNUD     ,NODE1    ,NODEI    ,NODEJ


      REAL*8,ALLOCATABLE:: ROWVECTOR(:)

C--------------------
C--------------------

      ALLOCATE(ROWVECTOR(LMXNDL))

      LMDIAG = NBAND * 2 + 2

      DO L=1,IA_ROWS

          ROWVECTOR=0D0

          DO I=2,NNUD

              ROWVECTOR(1) = ROWVECTOR(1)-A(L,I-1)
              ROWVECTOR (I) = A(L,I-1)

          ENDDO !I=2,NNUD

          NNUD=LNNDEL(L)

          DO I=1,NNUD-1 

              NODEI= KXX(I,L)  

              DO J=I+1,NNUD

                  NODEJ=KXX(J,L)

                  SELECT CASE (LBLOCK_NR)

                      CASE(1)
                          IJ_ROW = (NODEI - 1) * 2 + 1
                          IJ_COL = (NODEJ - 1) * 2 + 1
                          JI_ROW=  (NODEJ - 1) * 2 + 1
                          JI_COL = (NODEI - 1) * 2 + 1

                      CASE(2)
                          IJ_ROW = (NODEI - 1) * 2 + 1
                          IJ_COL = NODEJ * 2
                          JI_ROW = (NODEJ - 1) * 2 + 1
                          JI_COL = NODEI * 2 

                      CASE(3)
                          IJ_ROW =  NODEI * 2
                          IJ_COL=  (NODEJ - 1) * 2 + 1
                          JI_ROW = (NODEJ * 2)
                          JI_COL = (NODEI - 1) * 2 + 1

                      CASE(4)
                          IJ_ROW = NODEI * 2
                          IJ_COL = NODEJ * 2
                          JI_ROW = NODEJ * 2
                          JI_COL = NODEI * 2

                  END SELECT


                  IF (NODEI .LT. NODEJ) THEN
                      KINT = IJ_ROW
                      KEXT = JI_ROW
                      KINTBAND = LMDIAG + IABS(IJ_ROW - IJ_COL)
                      KEXTBAND = LMDIAG - IABS(JI_ROW - JI_COL)
                      KINT_IND =  J
                      KEXT_IND = I
                  ELSE
                      KINT = JI_ROW
                      KEXT = IJ_ROW
                      KINTBAND = LMDIAG + IABS(JI_ROW - JI_COL)
                      KEXTBAND = LMDIAG - IABS(IJ_ROW - IJ_COL)
                      KINT_IND = I
                      KEXT_IND = J
                  ENDIF


                  A_DSC(KINT,KINTBAND) = A_DSC(KINT,KINTBAND)
     ;                                 + ROWVECTOR(KINT_IND)*FACTOR
                  A_DSC(KEXT,KEXTBAND) = A_DSC(KEXT,KEXTBAND)
     ;                                 + ROWVECTOR(KEXT_IND)*FACTOR
              END DO ! J=I+1,NNUD

          END DO ! I=1,NNUD-1


C------------------------- Diagonal.
          DO I=1,NNUD

              SELECT CASE (LBLOCK_NR)
                  CASE(1)
                      NODE1 =    (KXX(I,L) - 1) * 2 + 1
                      IDIAG = LMDIAG
                  CASE(2)
                      NODE1 = (KXX(I,L) - 1) * 2 + 1
                      IDIAG = LMDIAG + 1
                  CASE(3)
                      NODE1 =    KXX(I,L) * 2
                      IDIAG = LMDIAG - 1
                  CASE(4)
                      NODE1 =    KXX(I,L) * 2
                      IDIAG = LMDIAG
              END SELECT


              A_DSC(NODE1,IDIAG) = A_DSC(NODE1,IDIAG)
     &                           + ROWVECTOR(I)*FACTOR

          END DO !I=1,NNUD

      END DO ! L=1,IA_ROWS

      END SUBROUTINE ASSEMBLE_CROSSTERM_INTO_SHUFFLED_BAND_NO_SYM
