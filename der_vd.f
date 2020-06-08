       SUBROUTINE DER_VD
     ; (IODIM    ,LMXNDL   ,MAINF    ,NFLAGS   ,NPAREL   ,NPPEL
     ; ,NUMEL    ,NUMNP    ,NZTRA    ,DER_H    ,DVDP     ,GRDFF
     ; ,IFLAGS   ,ISOZ     ,KXX      ,LDIM     ,LNNDEL   ,LXPAREL
     ; ,PAREL    ,VD)

********************************************************************************
*
* PURPOSE
*
*    Computes derivatives of Darcy's velocity for the current flow parameter
*
* DESCRIPTION
*
*    Computes derivatives of Darcy's velocity for the current flow parameter
*
* EXTERNAL VARIABLES: ARRAYS
*
* EXTERNAL VARIABLES: ARRAYS
*
*  DER_H                  Derivatives of nodal heads w.r.t current parameter
*  DVDP                   Derivatives of Darcy's velocity with respect to       
*                         flow parameters. In some cases it is dimensioned      
*                         as (IODIM,NUMEL) to reduce storage.                   
*  GRDFF                  Array containing the product between interpolation    
*                         functions integrals and interp. functions gradient    
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  ISOZ                   Anisotropy of every transmissivity zone               
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LDIM                   Vector containing the dimension of each element       
*  LNNDEL                 Number of nodes at every element                      
*  LXPAREL                Array containing zone numbers for a given             
*                         element parameter                                     
*  PAREL                  Parameter values at every element and current time    
*                         for all nodal parameters (each value is computed as   
*                         the product of up to four terms:                      
*                           elem. coeff*zonal value*time funct.*nonl. funct. )  
*
* INTERNAL VARIABLES: ARRAYS
*
*  IND                    Array used to simplify the coding of 
*                         conductivity matrix times the gradient of FEM 
*                         interpolation functions (shape functions)
*
* EXTERNAL VARIABLES: SCALARS
*
*  IODIM                  Maximum dimension of the problem                      
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFLAGS                 Maximum number of allowed flags                       
*  NPAREL                 Number of element parameters in current problem       
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NZTRA                  Number of transmissivity zones                        
*
* INTERNAL VARIABLES: SCALARS
*
*  IST                    Maximum between ISOZ and LDIM
*  L                      Counter index of elements                                                      
*  NNUD                   Number of nodes of the current element                
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*   None
*
* HISTORY
*
*     AMS      1-2002     First coding
*
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION PAREL(NUMEL,NPPEL),GRDFF(IODIM, LMXNDL, NUMEL), 
     ;    IFLAGS(NFLAGS), ISOZ(NZTRA), KXX(LMXNDL, NUMEL),
     ;    LDIM(NUMEL), LNNDEL(NUMEL), LXPAREL(NUMEL,NPAREL), 
     ;    DVDP(IODIM, NUMEL), DER_H(NUMNP) ,VD(IODIM, NUMEL)

       DIMENSION IND(3, 3, 3)

C------------------------- Array IND is used to simplify the coding of 
C------------------------- conductivity matrix times the gradient of FEM 
C------------------------- interpolation functions (shape functions)

       DATA ((IND(I,J,3),I=1,3),J=1,3)/1,4,5,4,3,6,5,6,2/
       DATA ((IND(I,J,2),I=1,2),J=1,2)/1,3,3,2/,IND(1,1,1)/1/

C------------------------- FIRST EXECUTABLE STATEMENT.

C------------------------- Cross over elements

       DO L = 1, NUMEL
          NNUD = LNNDEL(L)
          IDIM = LDIM(L)

          IST=MAX(ISOZ(LXPAREL(L,1)),IDIM)

C------------------------- Computes velocities at every element

          DO I = 1, IDIM
             S = 0.D0
             DO J = 1, NNUD 
                S2 = 0.D0
                DO K = 1, IDIM
                   IF (IND(I, K, IDIM).LE.IST) S2 = S2 + 
     ;                PAREL(L, IND(I, K, IDIM)) * GRDFF(K, J, L)
                ENDDO
                S = S - DER_H(KXX(J, L)) * S2
             ENDDO
*            if (l.eq.1) 
*    ;          WRITE(MAINF,*) '-->',i,(DER_H(KXX(J, L)),j=1,nnud),s
             DVDP(I, L) = S
          ENDDO
      

       ENDDO       ! Elements

C------------------------- Writes derivatives of velocities 

       IF (IFLAGS(33).EQ.1) THEN
          
          DO K=1,IODIM
             WRITE(MAINF,2100) ' COMPONENTE ',K,' DE LA VELOCIDAD'
             DO L = 1, NUMEL
                WRITE(MAINF, 2000) VD(K,L),DVDP(K, L)
             ENDDO
          ENDDO
 2000     FORMAT(2E23.16)
 2100     FORMAT(A,I5,A)
       ENDIF

       RETURN
       END
