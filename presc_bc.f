      SUBROUTINE PRESC_BC
     &          (A          ,A_DSC      ,B          ,IA_COLS
     &          ,IA_DSC_COLS,IA_DSC_ROWS,IADD_S     ,IADN_S
     &          ,IBCD       ,INDCHANGES ,INDFLTR    ,INDSSTR
     &          ,IONEWT     ,ISPARSE    ,ITYPADSC   ,KXX
     &          ,LMXNDL     ,LNNDEL     ,NBAND1     ,NPPNP
     &          ,NUMEL      ,NUMNP      ,PARNP      ,THETA
     &          ,VCALIT)


********************************************************************************
*
* PURPOSE  Performs all the computations related to presc. head boundary
*          conditions
*
* DESCRIPTION Updates RHS of flow equation (BFLU) and sets to zero the
*             row corresponding to the presc. head node (also the column, for
*             linear problems). Sets to 1.0 the diagonal term. Substracts
*             Aij*Hj for i=1,numnp, from BFLU to mantain the symmetric struc. of
*             AFLUDSC. The subroutine can be summarized as follows:
*             
*             - Step 0: Declaration of variables
*             - Step 1: Begins main loop over nodal points
*             - Step 2: Sets to zero current node row and column. Diagonal term
*                       is set to 1d0
*               - Step 2.1: Row is ALWAYS initialised
*               - Step 2.2: Initialises column whenever problem is currently in
*                           linear steady-state or, whitin time loop (transient
*                           state) if there are matrix changes (so that, must
*                           be assembled again)
*               - Step 2.3: In any case, sets to 1d0 the diagonal term
*               - Step 2.4: Assigns the value of prescribed head at the RHS of
*                           flow equation
*                 - Step 2.4a: Substracts Aij*Hj due to the symmetric structure
*                              of AFLUDSC
*
* EXTERNAL VARIABLES: ARRAYS
*
*  A                      Matrix of finite elements equations   
*                         No boundary conditions are included on it.            
*  A_DSC                  Coefficient matrix of flow system (2.15), most        
*                         often in its decomposed form.                         
*  B                      Right hand side of discretized equation.         
*  IBCD                   Boundary condition index                         
*  PARNP                  Parameter values at every node and current time for   
*                         all nodal parameters (each value is computed as the   
*                         product of up to four terms:                          
*                           nodal coeff*zonal value*time funct.*nonl. funct. )  
*
* INTERNAL VARIABLES: ARRAYS
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMAFLU               Used to dimension array AFLU                          
*  INDCHANGES             Indicates (<>0) whether or not flow matrix must be 
*                         assembled, for instance, because of parameters vary 
*                         with time
*  INDSSTR                Current problem state. If 1, transient state,         
*                         if 0, steady-state                                    
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  NBAND1                 Used to dimension. It is equal to NBAND+1             
*  NPPNP                  Total number of parameters by nodes (not confuse      
*                         with NPARNP, because in this casethere is no          
*                         difference between a given parameter in steady or tr.)
*  NUMNP                  Number of nodes                                       
*  THETA                 Time weighting parameter for flow problems            
*
* INTERNAL VARIABLES: SCALARS
*
*  I                      Nodal points dummy counter
*  I1                     Unimportant dummy counter                            
*  J                      Unimportant dummy counter                             
*  JI                     Unimportant dummy counter                             
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
* HISTORY: First coding                   German Galarza (Nov-97)
*          Revision and header inclusion  AAR (Nov-00)
*          Use of PARNP.                  AAR (Feb-02)
*          Adaptation for density dependent flow LJS
*
********************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IA_COLS    ,IA_DSC_COLS,IA_DSC_ROWS,INDFLTR
     &          ,INDCHANGES ,INDSSTR    ,IONEWT     ,ISPARSE
     &          ,ITYPADSC   ,LMXNDL     ,NBAND1     ,NPPNP
     &          ,NUMEL      ,NUMNP

      REAL*8::THETA

      INTEGER*4::IBCD(NUMNP)   ,IADN_S(NUMNP)
     &          ,IADD_S(NUMNP) ,KXX(LMXNDL,NUMEL)
     &          ,LNNDEL(NUMEL)

      REAL*8::A(NUMEL, IA_COLS)  ,A_DSC(IA_DSC_ROWS,IA_DSC_COLS)
     &       ,B(NUMNP)
     &       ,PARNP(NUMNP,NPPNP) ,VCALIT(NUMNP)

C------------------------- Internal
      INTEGER*4::I        ,INDPAR   ,INPOS    ,J        ,K        ,KNODE
     &          ,L        ,N        ,NNODE   ,NNUD


C------------------------- First executable estatement.

C------------------------- Step 0: Establish which row of PARNP to use:
C------------------------- that of flow or that of transport

      IF (INDFLTR.EQ.0) INDPAR = 1
      IF (INDFLTR.EQ.1) INDPAR = 4

C------------------------- Step 1: Begins main loop over nodal points

      DO I=1,NUMNP

          IF (IBCD(I).EQ.1) THEN ! Prescribed boundary condition

C------------------------- Step 2: Sets to zero current node row and
C------------------------- column. Diagonal term is set to 1d0.

              IF (INDCHANGES. NE.0 .OR. INDSSTR.EQ.0) THEN
                  CALL ZERO_ROW
     &            (IA_DSC_ROWS,IA_DSC_COLS,ISPARSE,NUMNP,I,IADN_S,A_DSC)


                  IF (ISPARSE.EQ.0) THEN

                      A_DSC(I,NBAND1) = 1D0

                  ELSE

                      A_DSC(IADD_S(I),I) = 1D0

                  END IF !ISPARSE.EQ.0

              END IF !INDCHANGES. NE.0


C------------------------- Step 3: set RHS term to the right value
C------------------------- Variable increment if Newton method is uses
C------------------------- or varable value otherwise.

              IF (IONEWT.EQ.1) THEN

                  B(I) = PARNP(I,INDPAR) - VCALIT(I)

              ELSE

                  B(I) = PARNP(I,INDPAR)

              END IF !IONEWT.EQ.0                                 

C------------------------- Step 5: Initialises column. Initializing
C------------------------- columns only for symmetric banded storage.
              IF (ITYPADSC.EQ.7) THEN

                  IF (INDCHANGES.NE.0 .OR. INDSSTR.EQ.0) THEN

                      DO J=1,MIN0(NBAND1-1,NUMNP-I)
                          A_DSC(I+J,NBAND1-J) = 0D0
                      END DO !J=1,MIN0(NBAND1-1,NUMNP-I)

                  END IF !INDCHANGES.NE.0

C------------------------- Step 6: Substracts theta*Aij*Hj due to the
C------------------------- symmetric structure of AFLUDSC.

                  DO L=1,NUMEL

                      NNUD = LNNDEL(L)

                      DO K=1,NNUD

                          KNODE = KXX(K,L)

C------------------------- If knode is the node with the bound. cond.

                          IF (KNODE.EQ.I) THEN

C------------------------- Loop over the nodes of l 	 	

                              DO N=1,NNUD

                                  NNODE = KXX(N,L)

C------------------------- Do not carry out operation for diagonal elements
C------------------------- nor for them with prescribed bound. cond.

                                  IF (N.NE.K .AND. IBCD(NNODE).NE.1)THEN

C------------------------- get position of element Ank
                                      IF (N.GT.K) THEN
                                         INPOS = (K-1)*NNUD+N-K*(K+1)/2
                                      ELSE IF(K.GT.N)THEN
                                         INPOS = (N-1)*NNUD+K-N*(N+1)/2
                                      END IF

                                      B(NNODE) = B(NNODE)
     &                                     - THETA*A(L,INPOS)*PARNP(I,1)

                                  END IF !N.NE.K .AND. IBCD(NNODE).NE.1
                              END DO !N=1,NNUD
                          END IF !KNODE.EQ.I
                      END DO !K=1,NNUD   
                  END DO !L=1,NUMEL
              END IF !ITYPADSC.EQ.7 .AND. INDCHANGES. NE.0 (band. symm.)
          END IF ! IBCD(I).EQ.1   
      END DO ! I=1.NUMNP

      END SUBROUTINE PRESC_BC
