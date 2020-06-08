      SUBROUTINE DERQ_GEN
     &         (IAD      ,IADN     ,IDIMDFLU ,INDSSTR  ,IPAR     ,NFLAGS
     &         ,NPAR     ,NPPNP    ,NUMNP    ,AFLU     ,CAUX1    ,DERC
     &         ,DERH     ,DERH1    ,DERH2    ,DFLU     ,IBCOD    ,IBTCO
     &         ,IFLAGS   ,PARNP    ,CAUDAL   ,IDIMAFLU ,ITYPAFLU 
     &         ,ITYPDFLU ,NUMEL    ,LMXNDL   ,KXX      ,LNNDEL)

********************************************************************************
*
* PURPOSE
*
*      Computes part of the derivatives of nodal flow w.r.t. flow parameters
*
* DESCRIPTION
*
*      The derivatives of nodal flow w.r.t. flow parameters is made of up two
*      parts. First part is common for all parameters. Second part is only 
*      for transmissivity and storage. 
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AFLU                   Matrix of finite elements equations for flow problem
*                         No boundary conditions are included on it.
*  DFLU                   Matrix of finite elements equations for flow
*                         problem related to storage term.
*  HAUX1                  Array containing heads, ponderated by THETAF time
*                         weight. Is equal to THETAF*HCAL+(1-THETAF)*HCALAN
*  HAUX2                  Array containing difference of heads in two
*                         consecutives times. Is equal to HCAL-HCALAN/TIME STEP
*  IBCOD                  Flow boundary condition index
*  PARNP                  Parameter values at every node and current time for
*                         all nodal parameters (each value is computed as the
*                         product of up to four terms:
*                            coeff*zonal value*time funct.*nonl. funct. )
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMDFLU               Used to dimension array DFLU
*  INDSSTR                Current problem state. If 1, transient state,
*                         if 0, steady-state
*  NBAND                  Half Bandwith (maximum difference between the
*                         numbers of two nodes belonging to the same element)
*  NBAND1                 Used to dimension. It is equal to NBAND+1
*  NPPNP                  Number of nodal parameters
*  NUMNP                  Number of nodes
*  THETAF                 Time weighting parameter for flow problems
*  TINC                   Current time increment
*
* INTERNAL VARIABLES: SCALARS
*
*  D_CAUD                 Derivative of nodal flow w.r.t. parameter IPAR
*                         (only implicit dependence through heads)
*  I                      Nodal counter
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  MUL_SS
*  MUL_TT
*
* HISTORY
*
*     AMS      3-2002     First coding
*
********************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMAFLU,IDIMDFLU,INDSSTR
     &          ,IPAR,ITYPAFLU,ITYPDFLU,LMXNDL,NFLAGS
     &          ,NPAR,NPPNP,NUMEL,NUMNP
      

      INTEGER*4::IAD(*),IADN(*),IBCOD(NUMNP),IBTCO(NUMNP)
     &           ,IFLAGS(NFLAGS),KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)

      REAL*8::AFLU(NUMEL,IDIMAFLU),CAUDAL(NUMNP),CAUX1(NUMNP)
     &       ,DERC(NUMNP,NPAR),DERH(NUMNP,NPAR)
     &       ,DERH1(NUMNP),DERH2(NUMNP),DFLU(NUMEL,IDIMDFLU)
     &       ,PARNP(NUMNP,NPPNP)
     

C------------------------- Internal

      INTEGER*4::I,IB

      REAL*8::CEXT,CNODE,DERQ

      REAL*8,ALLOCATABLE::D_CAUD(:)


C------------------------- FIRST EXECUTABLE STATEMENT.

      IF (IFLAGS(30).EQ.1) THEN

          WRITE(77,*) ' DERIVADAS DE D_CAUD (NO PARTE DIRECTA)' 

      END IF !IFLAGS(30).EQ.1


      ALLOCATE (D_CAUD(NUMNP))
       
      D_CAUD(:) = 0.D0


C------------------------- Estas cuentas sólo sirven para nivel fijo pero como
C------------------------- prod_mat_vec trabaja con toda la matriz, se hace
C------------------------- de todas formas y luego se anula en los nudos donde
C------------------------- no hay condición de nivel fijo.

      IF (ANY(IBCOD(:).EQ.1)) THEN

          IF (INDSSTR.NE.0) THEN

              CALL PROD_MAT_VEC
     &            (1D0      ,IAD      ,IADN     
     &            ,IDIMAFLU ,NUMEL    ,NUMNP    ,1
     &            ,ITYPAFLU ,1        ,LMXNDL   ,NUMEL    ,NUMNP
     &            ,KXX      ,LNNDEL   ,AFLU     ,D_CAUD
     &            ,DERH1)

              CALL PROD_MAT_VEC
     &            (1D0      ,IAD       ,IADN
     &            ,IDIMDFLU ,NUMEL    ,NUMNP    ,1
     &            ,ITYPDFLU ,1        ,LMXNDL   ,NUMEL    ,NUMNP
     &            ,KXX      ,LNNDEL   ,DFLU     ,D_CAUD
     &            ,DERH2)

          ELSE

	        CALL PROD_MAT_VEC
     &            (1D0      ,IAD      ,IADN     
     &            ,IDIMAFLU ,NUMEL    ,NUMNP    ,1
     &            ,ITYPAFLU ,1        ,LMXNDL   ,NUMEL    ,NUMNP
     &            ,KXX      ,LNNDEL   ,AFLU     ,D_CAUD
     &            ,DERH(1,IPAR))

          END IF !INDSSTR.NE.0
	END IF !ANY(IBCOD.EQ.1)

              
C------------------------- Cross over nodes

      DO I=1,NUMNP

          DERQ = 0D0

C------------------------- If there is in-flow and mass flow boundary condition...

          IF (CAUDAL(I).GT.0 .AND.
     &        (IBTCO(I).EQ.2 .OR.IBTCO(I).EQ.3) ) THEN

	        IB = IBCOD(I)

              SELECT CASE(IB)

C------------------------- Prescribed head boundary condition 
              
                  CASE(1) 

                      DERQ = D_CAUD(I)

                      IF (IFLAGS(30).EQ.1) THEN
                          WRITE(77,'(I5,2E23.16)') I,CAUDAL(I),DERQ
                      END IF

C------------------------- Leakage boundary condition 

                  CASE(3,4)
                
                      IF (INDSSTR.NE.0) THEN
                          DERQ = -PARNP(I,3)*DERH1(I)
                      ELSE
                          DERQ = -PARNP(I,3)*DERH(I,IPAR)
                      END IF !INDSSTR.NE.0


                      IF (IFLAGS(30).EQ.1) THEN

                          WRITE(77,'(2E23.16)') CAUDAL(I),DERQ

                      END IF !IFLAGS(30).EQ.1


              END SELECT !IB

C------------------------- Adds to inverse problem RHS

              CEXT = PARNP(I,4)
	        CNODE = CAUX1(I)
           
              DERC(I,IPAR) = DERC(I,IPAR) + DERQ*(CEXT-CNODE)

          END IF ! CAUDAL(I)>0 .AND. IBTCO=2,3


      END DO !I=1,NUMNP

      DEALLOCATE (D_CAUD)

      END SUBROUTINE DERQ_GEN
