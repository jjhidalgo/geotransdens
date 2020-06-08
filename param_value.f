       SUBROUTINE PARAM_VALUE
     ;(CFPARAM   ,DERSCC   ,DERSCH    
     ;,DTIM     ,EPS       ,IDIMFNT   ,INDSSTR  ,INTI      
     ;,IOFML    ,L         ,LMXNDL   
     ;,NFLAGS   ,NFNL      ,NFNLPARAM ,NFTIP     ,NFTPARAM ,NINT     
     ;,NNUD     ,NPARALG   ,NUMEL     ,NUMNP     ,NZPRG    ,PARAMC    
     ;,PARAMVAL ,PRGC      ,FNT        
     ;,VCALAN   ,VCALIT    ,IFLAGS    ,IPAR_DIR  ,KXX      ,NFNLPRG   
     ;,PARACD   ,XPARAM)
*****************************************************************************
* PURPOSE 
*
*   Computes a parameter value at the current node/element
* 
*
* DESCRIPTION
*
*   Computes a parameter value and its derivative with respect to the state 
*   variable (heads/concentrations). These values are computed as:
*
*        par=parC*CFpar*FTpar*FNLpar
*        d(par)/dV=parC*CFpar*FTpar*d(FNLpar)/dV
*
*   where V    is the state variable,
*   par        is the full value of the parameter  (PARAMVAL on output)
*   parC       is the zonal contribution (PARAMC argument on input)
*   CFpar      is the coefficient contribution
*   FTpar      is the time function contribution 
*   FNLpar     is the nonlinear function contribution 
*   d(par)/dV  is the derivative of the full value of the parameter with 
*              respect to the state variable V (DERVEC if 1<=IOCPAR <=3 
*              or DERSC otherwise on output)
*
* The parameter can be a parameter computed by node (1) or by element (2):
*
*    (1) The first argument, L, is then the current node value, and the 
*        second argument must be 1 (NNUD)
*
*    (2) L is the current element and NNUD its number of nodes. 
*        In this case, alternative formulations are enabled 
*        to compute the non-linear contribution to the
*        parameter value. These options are controlled by the IOFML
*        internal variable.
*
*              IOFML:  Controls the type of temporal interpolation
*
* The formulation is:
*
*    Internal Nodal Ponderation: Nonlinear contribution 
*    is computed at the average of the state variable.
*
*            Some notations:
*                     Vm, k    is the average value of state variable V
*                              between the nodes of the current element and 
*                              evaluated at time step k
*                     FNL      Nonlinear contribution 
*                     Vm,k+eps is defined as 
*                              eps*Vm,k+1 + (1-eps)* Vm,k
*
*    IOFML.EQ.1: Internal interpolation in time for FNL and d(FNL)/dV
*
*        FNL         is evaluated as   FNL[Vm,k+eps]              
*        d(FNL)/dV   is evaluated as   d(FNL)/dV [Vm,k+eps] 
*
*    IOFML.EQ.2: External interpolation in time for FNL. d(FNL)/dV is 
*                evaluated at time step k+1
*
*        FNL        is evaluated as   eps*FNL[Vm,k+1]+(1-eps)*FNL[Vm,k]
*        d(FNL)/dV  is evaluated as   d(FNL)/dV [Vm,k+1]
* 
*    IOFML.EQ.3:  FNL and d(FNL)/dV are evaluated at time step k+1
*
*         FNL       is evaluated as   FNL[Vm,k+1]
*        d(FNL)/dV  is evaluated as   d(FNL)/dV [Vm,k+1]
*
*
* 
*            Some notations:
*                        Vi, k    is the nodal value of state variable V at 
*                                 node i and at time step k
*
*    IOFML.EQ.1: Internal interpolation in FNL in time 
*
*        FNL        is evaluated as   MEAN ( FNL[Vi,k+eps] )  
*                           The MEAN is computed with a call to PARMEDIO
*
*        d(FNL/dV)  at node i is evaluated as   d(FNL/dV) [Vi,k+eps]
* 
*
*    IOFML.EQ.2: External interpolation in time 
*
*        FNL       is evaluated as  
*                  eps * MEAN ( FNL[Vi,k+1] ) + (1-eps)* MEAN ( FNL[Vi,k+1] )
*                           Each MEAN is computed with a call to PARMEDIO
*
*        d(FNL/dV) at node i  is evaluated as   d(FNL/dV) [Vi,k+1]     
*
*
*    IOFML.EQ.3: Internal interpolation in MEAN and external in FNL in time
*
*        FNL       is evaluated as
*                  MEAN ( eps*FNL[Vi,k+1]+(1-eps)*FNL[Vi,k] )
*
*                  MEAN is computed with a call to PARMEDIO
*
*       d(FNL/dV)  at node i  is evaluated as 
*                   eps*d(FNL/dV) [Vi,k+1] + (1-eps)*d(FNL/dV) [Vi,k]
*
* If the current parameter is linear, it is computed as:
*
*         par=PARC*CFPAR*NFTPAR
*         d(par)/dV=0
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FNT                    Array containing time function values                     
*  VCALAN                 Array containing computed state variable values 
*                         at last time step
*  VCALIT                 Array containing computed state variable values at 
*                         the last Newton iteration
*  KXX                    Node numbers of every element (counterclockwise
*                         order)
*  NFNLPRG                Generic parameter zone numbers for every nonlinear   `
*                         function
*  PARACD                 Agreement parameter values
*  PRGC                   Array of generic parameter values
*  XPARAM                 Auxiliar array containing parameters associated to 
*                         non-linear functions.
*
* EXTERNAL VARIABLES: SCALARS
*
*  CFPARAM                Element or Nodal coefficient of the current parameter
*  DERSCC                 Derivative of a given parameter with respect to the
*                         cconcentration
*  DERSCH                 Derivative of a given parameter with respect to the
*                         head level
*  DTIMEF                 Current time for the computation of flow time
*                         functions (counted since the beginning of the
*                         current observation interval, not since the beginning
*                         of the problem) divided by the length of the
*                         observation interval (interval between two
*                         consecutive observation times). Used only to make a
*                         linear interpolation of time function values.
*  EPS                    Time weighting parameter for nonlinear problems
*  IDIMFNT                First dimension of array FNT, it coincides with
*                         NFNT if the latter is not zero
*  INTI                   Observation time number such that the current
*                         computation time lies in between observation time
*                         number INTI and observation time number INTI+1
*  IOFML                  Time interpolation option that controls the 
*                         contribution of the nonlinear function
*  INDSSTR                Current problem state. If 1, transient state, 
*                         if 0, steady-state
*  L                      Current node/element
*  LMXNDL                 Maximum number of nodes per element
*  MAINF                  Unit number for MAIN FILE
*  NDENSFUNC              Number of the function describing the density
*                         dependency of the parameter
*  NFNL                   Total number of non-linear functions required
*                         in the problem
*  NFNLPARAM              Non-linear function number of the current parameter
*  NFTIP                  Type of the non-linear function of the current 
*                         parameter
*  NFTPARAM               Time function number of the current parameter 
*  NINT                   Number of observation times
*  NNUD                   Number of nodes of the current element. 1 if L is a
*                         node number. 
*  NUMEL                  Total number of elements                                                      
*  NUMNP                  Total number of nodes                                                      
*  NZPRG                  Total number of generic parameter zones
*  PARAMC                 Zonal value of the current parameter
*  PARAMVAL               Returned value of the current parameter
*
* INTERNAL VARIABLES: SCALARS
*
*  CMEDNEW                Mean Average of the concentration computed 
*                         at k+1 time step      
*  CMEDOLD                Mean Average of the concentration computed
*                         at k time step
*  DPARAMNL1              Non-linear functions derivative respect to the state
*                         variable                                                      
*  DPARNL                 Derivative of the non-linear function contribution 
*                         with respect to the state variable
*  DYDZ                   
*  HMEDNEW                Mean Average of the head computed 
*                         at k+1 time step      
*  HMEDOLD                Mean Average of the head computed
*                         at k time step
*  PARAM1                 Averame for parameter value returned from 
*                         parmedio subroutine                                                      
*  PARAM2                 The same as PARAM1                                                      
*  PARAMNL                Non-linear contribution computed in FUNNOLI
*  PARAMNL1               Non-linear contribution computed in FUNNOLI                                                      
*  PARAMNL2               Non-linear contribution computed in FUNNOLI                              
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  FUNNOLI                Computes non-linear function value contribution
*                         and its derivative with respect to the variable state
*  PARMEDIO               Computes different types of mean of a set of values
*
* HISTORY
*
*     SCR 20-jul-1997     First codding
*     AMS        5-98     Revision and correction of some mistakes
*     LJS     11-2003     Adaptation for variable density
*
*****************************************************************************
     
      IMPLICIT NONE
C EXTERNAL VARIABLES: SCALARS
      INTEGER*4 IDIMFNT,INDSSTR,INTI,IOFML
     ;         ,L,LMXNDL,NFLAGS,NFNL,NFNLPARAM,NFTIP
     ;         ,NFTPARAM,NINT,NNUD,NPARALG,NUMEL,NUMNP,NZPRG  
      REAL*8 CFPARAM, DTIM,EPS,PARAMC,PARAMVAL
     ;      ,PRGC,DERSCH,DERSCC
C EXTERNAL VARIABLES: ARRAYS
      INTEGER*4 IFLAGS(NFLAGS),IPAR_DIR(NPARALG)
     ;         ,KXX(LMXNDL,NUMEL),NFNLPRG(8,NFNL)  
      REAL*8 PARACD(3,NFNL),VCALAN(NUMNP),VCALIT(NUMNP)
     ;      ,XPARAM(8),FNT(IDIMFNT,NINT)
C INTERNAL VARIABLES: SCALARS
      INTEGER*4 I
      REAL*8 VMED,VMEDNEW,VMEDOLD
     ;      ,PARNL,PARAMNL1,PARAMNL2
     ;      ,DPARAMNL1H,DPARAMNL2H,DPARAMNL1C,DPARAMNL2C,FT1,TEMPCOEF

*       IF(IFLAGS(3).EQ.1)CALL IO_SUB('PARAM_VALUE',0)

C------------------------- Computes non linear contribution to parameter value
C------------------------- if the problem is transient and the current 
C------------------------- parameter has nonlinear function

       IF (INDSSTR.NE.0 .AND. NFNLPARAM.NE.0) THEN
       
          VMEDOLD = 0D0
          VMEDNEW = 0D0

         IF (NNUD.EQ.1) THEN 

C------------------------- Nodal parameter 
      
             VMEDOLD = VCALAN(L)
             VMEDNEW = VCALIT(L)

         ELSE

C------------------------- Parameter computed by elements. Average value of
C------------------------- state variable at the current element

            DO I=1,NNUD            


                 VMEDOLD = VMEDOLD + VCALAN(KXX(I,L))
                 VMEDNEW = VMEDNEW + VCALIT(KXX(I,L))

            END DO !I=1,NNUD

            VMEDOLD = VMEDOLD/NNUD
            VMEDNEW = VMEDNEW/NNUD


         ENDIF

C------------------------- Still internal nodal ponderation

         IF (IOFML.EQ.1) THEN

C------------------------- Internal interpolation in time'

            VMED = EPS*VMEDNEW+(1.D0-EPS)*VMEDOLD

            CALL FUNNOLI
     &(DPARAMNL1H    ,DPARAMNL1C    ,IPAR_DIR(7)   ,VMED
     &,NFLAGS        ,NFNL          ,NFNLPARAM     ,NFTIP
     &,NZPRG         ,PARNL         ,IFLAGS        ,NFNLPRG
     ;,PARACD        ,PRGC          ,XPARAM)  


            PARAMVAL=CFPARAM*PARAMC*PARNL
            DERSCH=DPARAMNL1H*PARAMC*CFPARAM/NNUD
            DERSCC=DPARAMNL1C*PARAMC*CFPARAM/NNUD

         ELSE  IF (IOFML.EQ.2 .OR. IOFML.EQ.3) THEN

              CALL FUNNOLI 
     &(DPARAMNL1H    ,DPARAMNL1C    ,IPAR_DIR(7)   ,VMEDOLD
     &,NFLAGS        ,NFNL          ,NFNLPARAM     ,NFTIP
     &,NZPRG         ,PARAMNL1      ,IFLAGS        ,NFNLPRG  
     &,PARACD        ,PRGC          ,XPARAM)  

               CALL FUNNOLI
     &(DPARAMNL2H    ,DPARAMNL2C    ,IPAR_DIR(7)   ,VMEDNEW
     &,NFLAGS        ,NFNL          ,NFNLPARAM     ,NFTIP
     &,NZPRG         ,PARAMNL2      ,IFLAGS        ,NFNLPRG  
     &,PARACD        ,PRGC          ,XPARAM)  

            IF (IOFML.EQ.2) THEN 
 

C------------------------- External interpolation in time

               PARAMVAL=( EPS*PARAMNL2+(1-EPS)*PARAMNL1 )*
     ;                                      CFPARAM*PARAMC

C------------------------- Derivative computed at time k+1

               DERSCH=CFPARAM*PARAMC*DPARAMNL2H/NNUD
               DERSCC=CFPARAM*PARAMC*DPARAMNL2C/NNUD


            ELSEIF (IOFML.EQ.3) THEN 

C------------------------- Nonlinear contribution and its derivative computed 
C------------------------- at time k+1

               PARAMVAL=CFPARAM*PARAMC*PARAMNL2
               DERSCH=CFPARAM*PARAMC*DPARAMNL2H/NNUD
               DERSCC=CFPARAM*PARAMC*DPARAMNL2C/NNUD

            ENDIF           ! IOFML = 2 or 3
         ENDIF              ! IOFML = 1

C------------------------- Linear parameter

       ELSE               ! When INDSSTR = 0 .OR. NFNLPARAM = 0

          PARAMVAL=CFPARAM*PARAMC

C------------------------- Derivatives set to zero (no sense in linear case)

          DERSCH=0D0 
          DERSCC=0D0
       ENDIF              !INDSSTR.NE.0 .AND. NFNLPARAM.NE.0

C------------------------- Computes time function contribution to parameter 
C------------------------- value

       IF (INDSSTR.NE.0 .AND. NFTPARAM.NE.0) THEN

          FT1=FNT(NFTPARAM,INTI)
          TEMPCOEF=(FNT(NFTPARAM,INTI+1)-FT1)*DTIM+FT1 
          PARAMVAL=PARAMVAL*TEMPCOEF

       END IF

*       IF(IFLAGS(3).EQ.1)CALL IO_SUB('PARAM_VALUE',1)

       RETURN
       END
