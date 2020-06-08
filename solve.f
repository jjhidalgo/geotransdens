      SUBROUTINE SOLVE
     ;(IADSC_ROWS    ,IADSC_COLS  ,IB_ROWS       ,IDESC       ,IDETAIL
     ;,IDIMWORK      ,IDIRECT     ,IEQN          ,INTI        ,IPLAN_B
     ;,IPRECOND      ,ISYMETRIC   ,ITERM         ,IWALGO
     ;,MAINF         ,NBAND1        ,NBF         ,NITMAX
     ;,NORTHMAX      ,NORTHSTART  ,RTWOTOL       ,RMAXTOL     ,SMAXTOL
     ;,A_DSC         ,A_DSCF      ,B             ,IAD         ,IADD
     ;,IADN          ,IAFD        ,IAFDD         ,IAFDN       ,WORK
     ;,SOLUTION)


**************************************************************************                        
*
*     PURPOSE
*     to solve systems of the kind Ax=b where a is a matrix and b and x vectors,
*     x being the vector iof unknowns. The routine can solve using direct methods
*     or using iterative methods.
*
*     NUEVAS VARIABLES:
C
c       ISYMETRIC: indicates wether the system matrix is symetric
c       IDIRECT: if 1, indicates that a direct solving method is desired. 
c                 otherwise an iterative solver is used. 
C       IPLAN_B:  This variable specifies what WATSOLV needs to do when 
c                 it does not converge.
c                 if IPLAN_B = 0,    watsolv will give up.
c                 If IPLAN_B = 1, watsolv will change the algorithm to CGSTAB 
c                 if IPLAN_B = 2, watsolv will increase the number of orthogonal 
c                 vectors used by GMRES.  
C       IPRECOND  specifies wether matrix preconditioning is required
C       NORTHSTART  specifies the number of orthogonal vectors to start with (for GMRES)
C       NORTHMAX  specifies the largest allowed  number of orthogonal vectors to start with (for GMRES)
c       IEQN  replaces the old variable IND 
C


      IMPLICIT NONE
C     EXTERNAL VARIABLES: SCALARS
      INTEGER*4  IADSC_ROWS,IADSC_COLS, IB_ROWS,IDIRECT, ISYMETRIC
     ;           , IDESC,NBAND1,IWALGO
     ;           ,MAINF,IEQN, IPRECOND, ITERM, NORTH, NORTHMAX
     ;           ,NORTHSTART, NITMAX, IDETAIL, IPLAN_B, INTI
     ;           ,NBF,IDIMWORK
      REAL*8 RTWOTOL,  RMAXTOL, SMAXTOL

C     EXTERNAL VARIABLES: ARRAYS
      INTEGER*4 IAD(IADSC_ROWS,IADSC_COLS),IADD(IB_ROWS),IADN(IB_ROWS)
     ;        ,IAFD(NBF,IB_ROWS),IAFDD(IB_ROWS)
     ;        ,IAFDN(IB_ROWS)
      REAL*8 A_DSC(IADSC_ROWS,IADSC_COLS),A_DSCF(NBF,IB_ROWS)
     ;       ,B(IB_ROWS)
     ;       ,WORK(IDIMWORK),SOLUTION(IB_ROWS)
       

C     INTERNAL VARIABLES: SCALARS
      CHARACTER*30 NOM
      CHARACTER*10 SUBNAME
      INTEGER*4 IER, NSWAPS,IJOB
      REAL*8 D1, D2
      LOGICAL SOLVE_FLAG



      IF (IEQN .EQ.1)THEN
        NOM='FLOW EQ.'
      ELSEIF (IEQN .EQ. 2) THEN
        NOM='TRANSPORT EQ.'
      ELSEIF (IEQN .EQ. 3) THEN
        NOM = 'COUPLED FLOW AND TRANSPORT EQN'
      ENDIF


C_______________________________________________________________
C__________________________PART 1_______________________________
C_________________________direct methods________________________
C_______________________________________________________________

      IF (IDIRECT .EQ.1 ) THEN

C____________________________________________PART 1A:    symetric matrix
        IF (ISYMETRIC .EQ. 1) THEN

C_____________set error indicator to 0
          IER = 0
C______________________________ Coefficient matrix has been factorized 
C                               previously. Therefore, only system solution is
C                               required
          IF (IDESC.EQ.1) THEN

            CALL LUELPB
     ;          (A_DSC    ,B          ,IB_ROWS   ,NBAND1 -1    ,IB_ROWS
     ;          ,SOLUTION)

C______________________________ Cholesky factorization and solution
          ELSE IF (IDESC.EQ.2) THEN

            CALL EQUAL_ARRAY(SOLUTION,B,IB_ROWS)
            SUBNAME ='LEQ1PB'

            CALL LEQ1PB
     ;          (A_DSC     ,IB_ROWS ,NBAND1 -1  ,IB_ROWS ,SOLUTION
     ;          ,IB_ROWS  ,1        ,D1      ,D2
     ;          ,IER)


C______________________________ Only factorization is required 
          ELSE

            SUBNAME ='LUDAPB'

            CALL LUDAPB
     ;          (A_DSC    ,IB_ROWS   ,NBAND1 -1  ,IB_ROWS      ,A_DSC
     ;          ,IB_ROWS  ,D1         ,D2        ,IER)

          ENDIF !IDESC.EQ.1,2...

C________________________________ PART 1B: non symetric matrix
        ELSE

          CALL EQUAL_ARRAY(SOLUTION,B,IB_ROWS)
          SUBNAME ='LEQT1B'

          IF (IDESC.EQ.2) IJOB = 0
          IF (IDESC.EQ.1) IJOB = 2

          CALL LEQT1B
     ;        (A_DSC     ,IB_ROWS     ,NBAND1-1    ,NBAND1-1   ,IB_ROWS
     ;        ,SOLUTION  ,1           ,IB_ROWS     ,IJOB       ,WORK
     ;        ,IER)

          IF (IER.EQ.129) THEN
            WRITE(6,*) ' AT INTI=',INTI,
     ;            ' TRANSPORT MATRIX IS ALGORITHMICALLY SINGULAR'
            STOP ' SINGULAR TRANSPORT'
          ENDIF

        ENDIF

C________________________________ PART 1C: error handling
        IF (IER .NE. 0) THEN

          IF (IER.EQ.129)  WRITE(MAINF,100 ) NOM,SUBNAME, INTI
  100     FORMAT(' ERROR IN SOLUTION OF ',/,A30,' IN SUBROUTINE',/,
     ;             A10,' OF DIRECT SOLVER',/,'TIME INTERVAL: ',I5)
          RETURN

       ENDIF !ISYMETRIC .EQ. 1
C_______________________________________________________________
C__________________________PART 2_______________________________
C______________________iterative methods________________________
C_______________________________________________________________
      ELSE


C_____________________________________if preconditioning is requiered 
           

        IF (IPRECOND .EQ. 1) CALL W_FACTOR
     ;(A_DSC    ,IB_ROWS   ,IADSC_ROWS,IAD    ,IADN
     ;,A_DSCF   ,IAFD      ,IAFDD     ,IAFDN  ,NBF)
         


        NORTH= NORTHSTART
        NSWAPS = 0
        !while the problem has not converged and while we have permis-
       ! sion to increase NORTH
        ITERM=4
        IF (IWALGO.EQ.1) SOLVE_FLAG = .FALSE.
        IF (IWALGO.EQ.2) SOLVE_FLAG = .TRUE.

        DO WHILE (ITERM .GT. 3 .AND.
     ;       ((NORTH .LE. NORTHMAX .AND. NSWAPS .LT.2)
     ;        .AND. (IPLAN_B .NE. 0)))

             !solve the system

          CALL WATSOLV
     ;(IB_ROWS  ,A_DSC       ,B           ,SOLUTION       ,NORTH
     ;,NITMAX   ,IDETAIL     ,SOLVE_FLAG  ,RTWOTOL        ,RMAXTOL
     ;,SMAXTOL  ,IAD         ,IADD        ,IADN           ,IADSC_ROWS
     ;,IAFD     ,IAFDD       ,IAFDN       ,NBF            ,A_DSCF
     ;,ITERM)



          IF (ITERM.GT.3) THEN

            !if the alternative plan is to change to cgstab
            IF (IPLAN_B .EQ. 1 .AND. NSWAPS.EQ.0) THEN
              SOLVE_FLAG = .NOT. SOLVE_FLAG
              NSWAPS = NSWAPS +1
            ENDIF


            !if the alternative plan is to increase north
            IF (IPLAN_B .EQ. 2 ) THEN
              NORTH =NORTH + 1
            ENDIF

          ENDIF !ITERM.GT.3
        ENDDO !WHILE (ITERM .GT. 3 .AND. ...

C________________________________ PART 2B: error handling

        IF (ITERM .GT. 3) THEN
          WRITE(MAINF,100 )NOM,SUBNAME,INTI

          RETURN
        ENDIF

      ENDIF !IDIRECT .EQ.1...

      RETURN
      END SUBROUTINE SOLVE
