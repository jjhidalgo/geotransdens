       SUBROUTINE DER_PARAM 
     ; (CFPARAM  ,DTIM     ,EPS      ,IDIMFNT            ,INDSSTR 
     ; ,INPRGC   ,INTI     ,IOCAP    ,IOCPAR   ,IOFML    ,IP 
     ; ,L        ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLPARAM
     ; ,NFNLTIP  ,NFTPAR   ,NINT     ,NNUD               ,NPTOT 
     ; ,NPZON    ,NUMEL    ,NUMNP    ,NZPAR    ,NZPRG    ,PARAMC 
     ; ,DERIV    ,FNT      ,IFLAGS   ,INDEX              ,IPOS 
     ; ,IVPAR    ,KXX      ,NFNLPRG  ,PARACD   ,PRGC     ,VCAL 
     ; ,VCALAN   ,XPARAM   ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR) 

***************************************************************************** 
* 
* EXTERNAL VARIABLES: ARRAYS 
* 
*  DERIV                  Its first component is the value of the current 
*                         nonlinear function. The rest of components are 
*                         the derivatives of the current nonlinear function 
*                         with respect to all the estimated parameters of this 
*                         function 
*  FNT                    Array containing time function values 
*  VCALAN                 Array containing computed state variable values 
*                         at last time step 
*  VCAL                   Array containing computed state variable values at 
*                         the current time step 
*  KXX                    Node numbers of every element (counterclockwise 
*                         order) 
*  IPOS                   Location order of the estimated parameter 
*  NFNLPRG                Generic parameter zone numbers for every nonlinear   ` 
*                         function 
*  PARACD                 Agreement parameter values 
*  PRGC                   Array of generic parameter values 
*  XPARAM                 Auxiliar array containing parameters associated to 
*                         non-linear functions. 
* 
* INTERNAL VARIABLES: ARRAYS 
* 
*  DERIV1                 Array containing the derivatives of the current 
*                         parameter, p, with respect to all the estimated 
*                         parameters which p is dependent to at previous time 
*  DERIV2                 Array containing the derivatives of the current 
*                         parameter, p, with respect to all the estimated 
*                         parameters which p is dependent to at current time 
*  DERIV                  Array containing the derivatives of the current 
*                         parameter, p, with respect to all the estimated 
*                         parameters which p is dependent to at required time 
* 
* EXTERNAL VARIABLES: SCALARS 
* 
*  CFPARAM                Element or Nodal coefficient of the current parameter 
*  DTIM                   Current time for the computation of time 
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
*  IOCPAR                 Nodal ponderation option to compute the contribution 
*                         of the nonlinear function to the value of the 
*                         current parameter 
*  IOFML                  Time interpolation option that controls the 
*                         contribution of the nonlinear function 
*  INDSSTR                Current problem state. If 1, transient state, 
*                         if 0, steady-state 
*  L                      Current node/element 
*  LMXNDL                 Maximum number of nodes per element 
*  NFNL                   Total number of non-linear functions required 
*                         in the problem 
*  NFNLPARAM              Non-linear function number of the current parameter 
*  NFTIP                  Type of the non-linear function of the current 
*                         parameter 
*  NFTPARAM               Time function number of the current parameter 
*  NINT                   Number of observation times 
*  NNUD                   Number of nodes of the current element. 1 if L is a 
*                         node number. 
*  NPTOT                  Total number of estimated parameters on the current 
*                         nonlinear function 
*  NUMEL                  Total number of elements 
*  NUMNP                  Total number of nodes 
*  NZPRG                  Total number of generic parameter zones 
*  PARAMC                 Zonal value of the current parameter 
* 
* INTERNAL VARIABLES: SCALARS 
* 
*  VMEDNEW                Mean Average of the state variable computed 
*                         at k+1 time step 
*  VMEDOLD                Mean Average of the state variable computed 
*                         at k time step 
* 
* FUNCTIONS AND SUBROUTINES REFERENCED 
* 
*  FUNNOLI_IN             Computes derivatives of the current non-linear 
*                         function with respect to estimated parameters which 
*                         are on its definition 
*  PARMEDIO               Computes different types of mean of a set of values 
* 
* HISTORY 
* 
*     AMS       11-98     First coding 
* 
***************************************************************************** 

       IMPLICIT DOUBLE PRECISION (A-H,O-Z) 

       DIMENSION 
     ;    NFNLPRG(8,NFNL)   ,PARACD(3,NFNL)   ,DERIV(NPAR) 
     ;   ,DERIV1(12)        ,DERIV2(12)       ,IPOS(NPAR) 
     ;   ,FNT(IDIMFNT,NINT) ,VCALAN(NUMNP)    ,VCAL(NUMNP) 
     ;   ,KXX(LMXNDL,NUMEL) ,IFLAGS(NFLAGS)   ,INDEX(12) 
     ;   ,TEMPCOEF(12)      ,NFTPAR(NZPAR)    ,CFPARAM(12) 
     ;   ,XPARAM(12)        ,NFNLTIP(NFNL)    ,IVPAR(NZPAR,4) 
     ;   ,WGT_PAR(NZPAR*IDIMWGT)              ,IPNT_PAR(NZPAR*IDIMWGT) 
  

C------------------------- Computes time function contribution to parameter 
C------------------------- value 
  

       IF (INDSSTR.NE.0) THEN 
  

          DO I=1,NPZON 

             IF (NFTPAR(INDEX(I)).NE.0) THEN 
                FT1=FNT(NFTPAR(INDEX(I)),INTI) 
                TEMPCOEF(I)=(FNT(NFTPAR(INDEX(I)),INTI+1)-FT1)*DTIM+FT1 
                XPARAM(I)=TEMPCOEF(I)*CFPARAM(I) 
             ELSE 
                XPARAM(I)=CFPARAM(I) 
             ENDIF 

          ENDDO 
       ELSE 
          XPARAM(1)=CFPARAM(1)                    ! "Steady state" 
       ENDIF 

       IF (IP.GT.0) THEN 

C------------------------- DERIVAUX is the derivative of P w.r.t zonal parameter 

          DERIVAUX=XPARAM(1) 

C------------------------- Derivative of Pz w.r.t the parameters it depends on 

          KONT=1 
          JJ=INDEX(1) 
          DO I=IVPAR(JJ,1),IVPAR(JJ,2) 
             DERIV(KONT)=WGT_PAR(I)*DERIVAUX
             IF (IVPAR(JJ,4).EQ.1) 
     ;           DERIV(KONT)=DERIV(KONT)*PARAMC*DLOG(10.0D0)
             IPOS(KONT)=IPNT_PAR(I) 
             KONT=KONT+1 
          ENDDO 
          NPTOT=IVPAR(JJ,2)-IVPAR(JJ,1)+1 

C------------------------- First position in array DERIV for the derivatives 
C------------------------- w.r.t. generic parameters 

          INI=NPTOT+1 

       ELSE 
          NPTOT=0 
          INI=1 
       ENDIF 
  

       IF (NFNLPARAM.NE.0) THEN 
  

C------------------------- Computes non linear contribution to parameter value 
C------------------------- if the problem is transient and the current 
C------------------------- parameter has nonlinear function 

          NFTIP=NFNLTIP(NFNLPARAM) 

C------------------------- If external nodal ponderation is required 

          IF (IOCPAR.GE.1 .AND. IOCPAR.LE.3) THEN 

             RETURN         !!!!!!!!!OJO NO FUNCIONA ESTA PARTE 

          ELSE              ! IOCPAR.NE. 1, 2, 3 
  
C------------------------- Internal Nodal ponderation 

             VMEDOLD=0D0 
             VMEDNEW=0D0 
             IF (NNUD.EQ.1) THEN 

C------------------------- Nodal parameter 

                VMEDOLD=VCALAN(L) 
                VMEDNEW=VCAL(L) 
             ELSE 

C------------------------- Parameter computed by elements. Average value of 
C------------------------- state variable at the current element 

                DO I=1,NNUD 
                   VMEDOLD=VMEDOLD+VCALAN(KXX(I,L)) 
                   VMEDNEW=VMEDNEW+VCAL(KXX(I,L)) 
                ENDDO 

                VMEDOLD=VMEDOLD/NNUD 
                VMEDNEW=VMEDNEW/NNUD 
             ENDIF 

C------------------------- Still internal nodal ponderation 

             IF (IOFML.EQ.1) THEN 

C------------------------- Internal interpolation in time 
  
                VMED=EPS*VMEDNEW+(1-EPS)*VMEDOLD 

                CALL FUNNOLI_IN 
     ;(INDEX(1)  ,INDEX(2)  ,INPRGC   ,IOCAP    ,NFLAGS   ,NFNL 
     ;,NFNLPARAM ,NFTIP     ,NPTOT    ,NZPAR    ,NZPRG    ,PARAMC 
     ;,VMED      ,DERIV(INI),IFLAGS   ,IPOS(INI),IVPAR    ,NFNLPRG 
     ;,PARACD    ,PRGC      ,XPARAM) 

             ELSE  IF (IOFML.EQ.2) THEN 

                CALL FUNNOLI_IN 
     ;(INDEX(1)  ,INDEX(2)  ,INPRGC   ,IOCAP    ,NFLAGS   ,NFNL 
     ;,NFNLPARAM ,NFTIP     ,NPTOT    ,NZPAR    ,NZPRG    ,PARAMC 
     ;,VMEDOLD   ,DERIV1    ,IFLAGS   ,IPOS(INI),IVPAR    ,NFNLPRG 
     ;,PARACD    ,PRGC      ,XPARAM) 

                CALL FUNNOLI_IN 
     ;(INDEX(1)  ,INDEX(2)  ,INPRGC   ,IOCAP    ,NFLAGS   ,NFNL 
     ;,NFNLPARAM ,NFTIP     ,NPTOT    ,NZPAR    ,NZPRG    ,PARAMC 
     ;,VMEDNEW   ,DERIV2    ,IFLAGS   ,IPOS(INI),IVPAR    ,NFNLPRG 
     ;,PARACD    ,PRGC      ,XPARAM) 

  
                DO I=1,NPTOT 
                   DERIV(INI+I-1)=EPS*DERIV2(I)+(1-EPS)*DERIV1(I) 
                ENDDO 
  
             ENDIF              ! IOFML = 1 
          ENDIF                 ! 1 <= IOCPAR <= 3 

          NPTOT=NPTOT+INI-1  ! Includes all unk. param. (generics and the rest)

       ENDIF              ! NFNLPARAM.NE.0 

       RETURN 
       END 
