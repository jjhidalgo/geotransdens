      SUBROUTINE FUNNOLI_IN
     ;(IND1      ,IND2     ,INPRGC   ,IOCAP    ,NFLAGS    ,NFNL      
     ;,NFNLPARAM ,NFTIP    ,NPTOT    ,NZPAR    ,NZPRG     ,PARAMC    
     ;,X         ,DERIV    ,IFLAGS   ,IPOS     ,IVPAR     ,NFNLPRG   
     ;,PARACD    ,PRGC     ,XPARAM)  

*******************************************************************************
*
* PURPOSE  Computing the derivatives of non-linear functions with respect to
*          all concerning physical and generic parameters 
*
*
* DESCRIPTION  This subroutine includes derivatives functions concerning to 
*              transport and flow phaenomena. All of them are related to the 
*              computation of derivatives od matrix coefficients (AFLU-DFLU-
*              BLU-ATRA-DTRA-BTRA) with respect to physical and generic 
*              parameters
*              For each kind of non-linear function, several computations are
*              done. If a coefficient is defined as:
*                        
*                           a(i,j)=W*CFz*FTz*Pz*FNL(Gi,Ps,Pz,state var.)
*              
*              then computes:
*
*                    1)  d(a(i,j))              d(FNL)
*                        --------- =CFz*(FNL+Pz--------)
*                           dPz                  dPz
*
*                    2)  d(a(i,j))          d(FNL)
*                        --------- =CFz*Pz*--------   ; for all i 
*                           dGi              dGi
*
*                    3)  d(a(i,j))          d(FNL)
*                        --------- =CFz*Pz*--------   ; for all s
*                           dPs               dPs
*
*              The product VALUE times FTz is perform in DER_PARAM, and the
*              last one (*W) in the main subroutine (DERTRA ,DERSTG...)
*              Expression(2) is perform for all of the generic parameters appearing
*              in a explicit form in the expression of FNL and for those that
*              appear in a implicit form. For an example, relative permeability of
*              Van Genuchten only depends (in a explicit form) on LAMBDA, minimum
*              and maximum saturation degrees, but also, it depends on bubble 
*              pressure trough the saturation degree.
*              On the other hand, expression (3) is perform for all of the non-
*              generic parameters appearing in the definition of the non linear 
*              function. For example, retention curves are computed as:
*
*                                    Porosity          dSw
*                            Sw+ ----------------*--------------
*                                Specific storage  d(state var.)
*
*              In this case, Ps would be the porosity
*              Be careful with the expression of the type:
*
*                      dFNL
*                     ------- because FNL is a function of (CFz*Pz*FTz)
*                       dPz   
*
*              so, we must compute:
*
*                      dFNL                    dFNL
*                     -------= CFz*FTz*-----------------
*                       dPz              dPz (complete)
*                     
* 
*
* EXTERNAL VARIABLES: ARRAYS
*
*  DERIV                  Its first component is the value of the current 
*                         nonlinear function. The rest of components are 
*                         the derivatives of the current nonlinear function 
*                         with respect to all the estimated parameters of this
*                         function
*  IFLAGS                 Array with different writing options. Used mainly for
*                         debugging.
*  IPOS                   Location order of the estimated parameter
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  NFNLPRG                Generic parameter zone number for every nonlinear
*                         function
*  PARACD                 Agreement parameters
*  PRGC                   Array containing generic parameters (subset of PARC)
*  XPARAM                 Array containing porosity divided by storavity on the
*                         first component and the Porosity on the second one 
*
* EXTERNAL VARIABLES: SCALARS
*
*  CFPARAM                Element or Nodal coefficient of the current parameter
*  IND1                   Zone number of the current parameter     
*  IND2                   Zone number of the non-generic parameter
*  IOCAP                  Parameterization scheme. If =0 --> capacitive
*                                                  else  --> conservative
*  NFLAGS                 Maximum number of allowed flags                       
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NFNLPARAM              Absolut order number of the current non-linear 
*                         function                            
*  NFTIP                  Non linear kind of function:
*                            (From 1 to 20)  Transmissivity or hyd. conductivity
*                            (From 21 to 40) Storage / retention curves
*                            (From 41 to 45) Leakage
*                            .
*                            .    To be completed...
*                            .
*  NPTOT                  Total number of estimated parameters on the current 
*                         non-linear function plus one.    
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  NZPRG                  Total number of generic parameter zones               
*  PARAMC                 Zonal value of the current parameter
*  IOCAP                  Parameterization scheme. If =0 --> capacitive
*                                                  else  --> conservative
*  X                      Independent variable (state variable at time=k+1)
*
* HYSTORY: First coding: Andres Alcolea (March-1999)
*          Structure;    Agustin Medina (Jan-February 1999)
*
********************************************************************************                                                                       

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION XPARAM(8)         ,NFNLPRG(8,NFNL)    ,PRGC(NZPRG) 
     ;         ,IFLAGS(NFLAGS)    ,PARACD(3,NFNL)     ,IVPAR(NZPAR)
     ;         ,IPOS(12)          ,DERIV(12)

      IF(IFLAGS(3).NE.0) CALL IO_SUB('FUNNOLI_IN',0)


      IF (NFTIP.EQ.1) THEN

C______________________________ Unconfined aquifers
C______________________________ Asymptotic free aquifer transmissivity

        XBASE=XPARAM(2)                 ! Bottom of the aquifer
        DELTA=PARACD(1,NFNLPARAM)       ! Agreement parameter
        CFFT=XPARAM(1)

C______________________________ Computation of non-linear function value
C______________________________ It is stored in DERIV(1). Only if zonal paramet`
C______________________________ is going to be estimated. In this case Tij(zone`

        IC1=1                   ! Estimated parameter index counter
        IPOS(IC1)=IVPAR(IND1)   ! IND1 is the zonal parameter zone number

        IF (IPOS(IC1).NE.0) THEN

          IF (X.GE.(XBASE+DELTA)) THEN    ! Linear zone (right)

            Y=(X-XBASE)
            DYDX=1D0
  
          ELSE                            ! Exponential zone (left)
  
            Y=DELTA*DEXP((X-XBASE-DELTA)/DELTA)
            DYDX=Y/DELTA

          ENDIF

C______________________________ In any case, non-linear function depends in a
C______________________________ explicit way on the zonal parameter Tij(zone k)

          DFNLDPZ=0D0

          DERIV(IC1)=CFFT*(Y+PARAMC*DFNLDPZ)
          IC1=IC1+1                                    ! Updates counter

        END IF

C______________________________ This type of non-linear function does not
C______________________________ depend on any zonal or generic parameter.
C______________________________ So, we have finished dirty work

        NPTOT=IC1-1

      ELSE IF (NFTIP.EQ.2) THEN

C______________________________ Unconfined aquifer
C______________________________ Parallel free aquifer transmissivity
C______________________________ X= Current head level
C______________________________ Y= Saturated thickness

        XBASE=XPARAM(2)                  ! Bottom of the aquifer (Card B1.2)
        DELTA=PARACD(1,NFNLPARAM)          ! 0.5*Length of the transicion zone
        X1=(X-XBASE)                     ! Relative head level
        CFFT=XPARAM(1)

C______________________________ Computation of non-linear function value
C______________________________ It is stored in DERIV(1). Only if zonal paramet`
C______________________________ is going to be estimated. In this case Tij(zone`

        IC1=1                   ! Estimated parameter index counter
        IPOS(IC1)=IVPAR(IND1)   ! IND1 is the zonal parameter zone number

        IF (IPOS(IC1).NE.0) THEN

          IF (X1.GT.DELTA)THEN             ! Linear zone (right)
   
            Y=(X-XBASE)
            DYDX=1D0
   
          ELSE IF (X1.GE.-DELTA)THEN        ! Transition zone (agreement) 
   
            TMIN=PRGC(NFNLPRG(1,NFNLPARAM))       ! Treshold value for B(h)
            E=DELTA
  
C______________________________ Computes the coeff. of the agreement 
  
            D=TMIN/(E*(DEXP(E)+DEXP(-E))-(DEXP(E)-DEXP(-E)))
            B=(1D0-D*(DEXP(E)+DEXP(-E)))/2D0
            A=(D*DEXP(-E)+B)/2D0/E           
            C=TMIN-D*DEXP(-E)+B*E-A*E*E
  
            Y=A*X1*X1+B*X1+C+D*DEXP(X1)
            DYDX=2*A*X1+B+D*DEXP(X1)
  
          ELSE
  
C______________________________ Constant zone (left)
  
            TMIN=PRGC(NFNLPRG(1,NFNLPARAM))       ! Treshold value for B(h)
            Y=TMIN
            DYDX=0D0
  
          ENDIF

C______________________________ In any case, non-linear function depends in a
C______________________________ explicit way on the zonal parameter Tij(zone k)

          DFNLDPZ=0D0

          DERIV(IC1)=CFFT*(Y+PARAMC*DFNLDPZ)
          IC1=IC1+1                                    ! Updates counter

        END IF
C______________________________ Computation of the derivatives of FNL with 
C______________________________ respect to all of its generic parameters 
C______________________________ (lowest). It means that current non-linear
C______________________________ function could depend on the generic parameters 
C______________________________ of another non-linear one.

C______________________________ First (and only) generic parameter: TMIN

        IPOS(IC1)=IVPAR(INPRGC+NFNLPRG(1,NFNLPARAM))
        IF (IPOS(IC1).NE.0) THEN

          IF (X1.GT.DELTA)THEN              ! Linear zone (right)

            DERYDERB=0D0

          ELSE IF (X1.GE.-DELTA)THEN        ! Transition zone (agreement) 
             
            AUX=DEXP(DELTA)*(1D0-DELTA)-DEXP(-1D0*DELTA)*(1D0+DELTA)
            DERDDERB=-1D0/AUX
            DERADERB=(DEXP(-1D0*DELTA)-DEXP(DELTA))*
     ;               DERDDERB/(4D0*DELTA)
            DERGDERB=-0.5D0*(DEXP(-1D0*DELTA)+DEXP(DELTA))*DERDDERB
            DERCDERB=-1D0*DEXP(DELTA)*DERDDERB-DELTA*DELTA*DERADERB
     ;               -DELTA*DERGDERB
    
            DERYDERB=DERADERB*X1*X1+DERGDERB*X1+DERCDERB
     ;               +DERDDERB*DEXP(X1)

          ELSE                              ! Constant zone (left)

            DERYDERB=1D0

          END IF
            
          DERIV(IC1)=CFFT*PARAMC*DERYDERB            
          IC1=IC1+1

        END IF

        NPTOT=IC1-1

      ELSE IF (NFTIP.EQ.3) THEN

C______________________________ Broadbridge and White relative hyd. conduct.

        CONTINUE                ! Coming soon...

      ELSE IF (NFTIP.EQ.4) THEN

C______________________________ Relative hydraulic conductivity of Van Genuchten

        XLAMBDA=PRGC(NFNLPRG(1,NFNLPARAM))  ! Van Genuchten exponent
        SMIN=PRGC(NFNLPRG(2,NFNLPARAM))     ! Minimum saturation degree
        SMAX=PRGC(NFNLPRG(3,NFNLPARAM))     ! Maximum saturation degree
        PCAP=PRGC(NFNLPRG(4,NFNLPARAM))     ! Bubble pressure 
        DELTA=PARACD(1,NFNLPARAM)           ! Length of agreement func. int.
        PERM=XPARAM(1)*PARAMC               ! COnductivity
        CFFT=XPARAM(1)                      ! (CFTRA*FT) of K_r

C______________________________ Computation of non-linear function value
C______________________________ It is stored in DERIV(1). Only if zonal parameter
C______________________________ is going to be estimated. In this case Tij(zone k)

        IC1=1                   ! Estimated parameter index counter
        IPOS(IC1)=IVPAR(IND1)   ! IND1 is the zonal parameter zone number

        IF (IPOS(IC1).NE.0) THEN

          IF (X.GE.0D0)THEN                       

C______________________________ Saturated state

            Y=1D0                               ! Hyd. Cond. = sat. hyd. cond.
            DYDX=0D0
         
          ELSE                                  

C______________________________ Unsaturated state.
C______________________________ Computes some auxiliar terms
C______________________________ and Sw and dSw/dh

            RELP=(-X/PCAP)**(1D0/(1D0-XLAMBDA))          
            XMULT=XLAMBDA*(SMIN-SMAX)/(1D0-XLAMBDA)/PCAP*RELP**XLAMBDA 

            SATW=SMIN+(SMAX-SMIN)*(1D0+RELP)**(-XLAMBDA) ! Sat. degree
            DERSATW=XMULT*(1D0+RELP)**(-XLAMBDA-1D0)     ! First order deriv.
            

            IF (SATW.GT.(SMAX-DELTA))THEN    

C______________________________ Zone of agreement function. The agreement func.
C______________________________ is the sum of a parabole and an exponential term

              SDELTA=1D0-DELTA/(SMAX-SMIN)                 ! Auxiliar term
              RADDELTA=SQRT(SDELTA)                        ! Auxiliar term
              SDPOT=1D0-SDELTA**(1D0/XLAMBDA)              ! Auxiliar term

C______________________________ Computes some auxiliar terms:
C______________________________ K1 is K (hyd. cond.) for SATW=SMAX-DELTA and
C______________________________ K2 is the deriv. of K with respect to SATW when
C______________________________ SATW=SMAX-DELTA

              K1=RADDELTA*(1D0-SDPOT**XLAMBDA)**2D0     
              K2=(1D0/(2D0*RADDELTA)*(1D0-SDPOT**XLAMBDA)**2D0+
     ;           2D0*RADDELTA*(1D0-SDPOT**XLAMBDA)*SDPOT**(XLAMBDA-1D0)
     ;           *SDELTA**((1D0-XLAMBDA)/XLAMBDA))/(SMAX-SMIN) 


              EX=DEXP(DELTA)
              D=(1D0-DELTA/2D0*K2-K1)/(EX*(1D0-DELTA/2D0)-DELTA/2D0-1D0)
              C=K1-D
              B=K2-D
              A=(D*(1D0-EX)-K2)/2D0/DELTA

              S=SATW+DELTA-SMAX             ! Substitution
              Y=A*S**2D0+B*S+C+D*DEXP(S)    ! Value of K
              DKDS=2D0*A*S+B+D*DEXP(S)      ! Derivative of K with respect to S
              DYDX=DKDS*DERSATW             ! Derivative of K with respect to X

            ELSE                            

C______________________________ Zone of van Genuchten function

              SATEF=(SATW-SMIN)/(SMAX-SMIN)        ! Auxiliar term
              RAD=DSQRT(SATEF)                     ! Auxiliar term
              SPOT=1D0-SATEF**(1D0/XLAMBDA)        ! Auxiliar term  

              Y=RAD*(1-SPOT**XLAMBDA)**2D0
              DKDS=(1D0/(2D0*RAD)*(1D0-SPOT**XLAMBDA)**2D0+
     ;        2D0*RAD*(1D0-SPOT**XLAMBDA)*SPOT**(XLAMBDA-1D0)*
     ;        SATEF**((1D0-XLAMBDA)/XLAMBDA))/(SMAX-SMIN)
              DYDX=DKDS*DERSATW

            ENDIF

          ENDIF

C______________________________ In any case, non-linear function depends in a 
C______________________________ explicit way on the zonal parameter Tij(zone k)

          DFNLDPZ=0D0

          DERIV(IC1)=CFFT*(Y+PARAMC*DFNLDPZ)   
          IC1=IC1+1                                    ! Updates counter

        END IF

C______________________________ Computation of the derivatives of FNL with 
C______________________________ respect to all of its generic parameters 
C______________________________ (lowest). It means that current non-linear
C______________________________ function could depend on the generic parameters 
C______________________________ of another non-linear one.

        IF (X.LT.0D0) RELP=(-1D0*X/PCAP)**(1D0/(1D0-XLAMBDA))          
        RELSUC=-1D0*X/PCAP

C______________________________ Due to the appearance of DER(K)/DER(SATEFF) in
C______________________________ the whole develope of this part of the subrou.
C______________________________ , it is always computed

        SATW=SMIN+(SMAX-SMIN)*(1D0+RELP)**(-1.D0*XLAMBDA) ! Sat. degree
        IF(X.GE.0D0) SATW=SMAX
        XMU=(SATW-SMIN)/(SMAX-SMIN)

        IF(X.GE.0D0) THEN
          DKDMU=0D0
        ELSE IF (XMU.LE.1D0-DELTA) THEN     ! Left zone: Van Genuchten function
          AUX=1D0-XMU**(1D0/XLAMBDA)
          DKDMU=2.D0*XMU**((2D0-XLAMBDA)/(2D0*XLAMBDA))
     ;         *(1D0-AUX**XLAMBDA)
     ;         *AUX**(XLAMBDA-1D0)+0.5D0*(1D0-AUX**XLAMBDA)
     ;         *(1D0-AUX**XLAMBDA)/DSQRT(XMU)
        ELSE                               ! Agreement function
          SDELTA=1D0-DELTA/(SMAX-SMIN)                 ! Auxiliar term
          RADDELTA=SQRT(SDELTA)                        ! Auxiliar term
          SDPOT=1D0-SDELTA**(1D0/XLAMBDA)              ! Auxiliar term
          AUX=XMU+DELTA-1D0                            ! Auxiliar term

C______________________________ Computes some auxiliar terms:
C______________________________ K1 is K (hyd. cond.) for SATW=SMAX-DELTA and
C______________________________ K2 is the deriv. of K with respect to SATW when
C______________________________ SATW=SMAX-DELTA

          K1=RADDELTA*(1D0-SDPOT**XLAMBDA)**2D0     
          K2=(1D0/(2D0*RADDELTA)*(1D0-SDPOT**XLAMBDA)**2D0+
     ;       2D0*RADDELTA*(1D0-SDPOT**XLAMBDA)*SDPOT**(XLAMBDA-1D0)
     ;       *SDELTA**((1D0-XLAMBDA)/XLAMBDA))/(SMAX-SMIN) 


          EX=DEXP(DELTA)
          D=(1D0-DELTA/2D0*K2-K1)/(EX*(1D0-DELTA/2D0)-DELTA/2D0-1D0)
          B=K2-D
          A=(D*(1D0-EX)-K2)/2D0/DELTA

          DKDMU=2D0*A**AUX+B+D*DEXP(AUX)
      
        END IF

C______________________________ First generic parameter: XLAMBDA

        IPOS(IC1)=IVPAR(INPRGC+NFNLPRG(1,NFNLPARAM))
        IF (IPOS(IC1).NE.0) THEN

          IF (X.LT.0D0) THEN       ! Non-linear function 

            U=1D0+RELSUC**(1D0/(1D0-XLAMBDA))       ! Sat. degree base (aux)
            V=1D0-XMU**(1D0/XLAMBDA)                ! Rel .perm base (aux)

            DERUDERLAM=RELSUC**(1D0/(1D0-XLAMBDA))*DLOG(RELSUC)
     ;                /((1D0-XLAMBDA)*(1D0-XLAMBDA))
            DERVDERLAM=XMU**(1D0/XLAMBDA)*DLOG(XMU)/(XLAMBDA*XLAMBDA)

            DERMUDERLAM=-1D0*XLAMBDA*U**(-1D0*XLAMBDA-1D0)
     ;                 *DERUDERLAM-U**(-1D0*XLAMBDA)*DLOG(U)
          
            DKDLAM=2D0*DSQRT(XMU)*(V**XLAMBDA-1D0)
     ;            *(XLAMBDA*V**(XLAMBDA-1D0)*DERVDERLAM
     ;            +V**XLAMBDA*DLOG(V))

            DERIV(IC1)=PERM*(DKDMU*DERMUDERLAM+DKDLAM)

          ELSE                           ! Constant, therefore, it does not
                                         ! depend on generic parameters
            DERIV(IC1)=0.D0

          END IF
          IC1=IC1+1

        END IF

C______________________________ Second generic parameter: minimum saturation degree

         
        IPOS(IC1)=IVPAR(INPRGC+NFNLPRG(2,NFNLPARAM))
        IF (IPOS(IC1).NE.0) THEN
          DERIV(IC1)=0.D0        
          IC1=IC1+1
        END IF

C______________________________ Third generic parameter: maximum saturation degree

        IPOS(IC1)=IVPAR(INPRGC+NFNLPRG(3,NFNLPARAM))
        IF (IPOS(IC1).NE.0) THEN
          DERIV(IC1)=0.D0
          IC1=IC1+1
        END IF

C______________________________ Fourth generic parameter: bubble pressure  

        IPOS(IC1)=IVPAR(INPRGC+NFNLPRG(4,NFNLPARAM))
        IF (IPOS(IC1).NE.0) THEN

          DMUDH=(RELSUC/PCAP)*(XLAMBDA/(1.D0-XLAMBDA))*
     ;           RELSUC**(XLAMBDA/(1D0-XLAMBDA))*
     ;           (1.D0+RELSUC**(1.D0/(1D0-XLAMBDA)))**(-XLAMBDA-1.D0)


          DERIV(IC1)=PERM*DKDMU*DMUDH
          IC1=IC1+1
        END IF

C______________________________ This type of non-linear function does not 
C______________________________ depend on any zonal parameter, such as porosity.
C______________________________ So, we have finished dirty work        

        NPTOT=IC1-1

      ELSE IF (NFTIP.EQ.5) THEN

C______________________________ Brooks & Corey relative hydraulic conductivity 

        CONTINUE                ! Coming soon...

      ELSE IF (NFTIP.EQ.21) THEN

C______________________________ Two layer storage non-linear function
C______________________________ Free aquifers
C______________________________ X= Current head level
C______________________________ Y= Factor such us S=STGC*FACTOR

        SM=PARACD(1,NFNLPARAM)          ! Minimum relative value of STGC
        DELTA=PARACD(2,NFNLPARAM)       ! Agreement parameter
        XBASE=XPARAM(2)                 ! Bottom of the aquifer
        CFFT=XPARAM(1)                  ! CFPAREL*FNT

C______________________________ Computation of non-linear function value
C______________________________ It is stored in DERIV(1). Only if zonal paramet`
C______________________________ is going to be estimated. In this case Tij(zone`

        IC1=1                   ! Estimated parameter index counter
        IPOS(IC1)=IVPAR(IND1)   ! IND1 is the zonal parameter zone number

        IF (IPOS(IC1).NE.0) THEN
          
          IF (X.GE.(XBASE+DELTA)) THEN    

C______________________________ Constant zone (right)

            Y=1D0
            DYDX=0D0

          ELSE IF(X.LT.(XBASE-DELTA)) THEN        

C______________________________ Constant zone (left)

            Y=SM
            DYDX=0D0

          ELSE                                     

C______________________________ Transition zone

            X1=(X-XBASE)/DELTA
            IF(X1.LT.0D0) FACTOR=-1D0
            Y=(1D0-SM)*(3D0*X1-FACTOR*DABS(X1)**3)/4+(1+SM)/2
            DYDX=(1D0-SM)*3D0*(1D0-X1*X1)/(4*DELTA)

          ENDIF

C______________________________ In any case, non-linear function depends in a
C______________________________ explicit way on the zonal parameter Tij(zone k)

          DFNLDPZ=0D0

          DERIV(IC1)=CFFT*(Y+PARAMC*DFNLDPZ)
          IC1=IC1+1                                    ! Updates counter

        END IF

C______________________________ This type of non-linear function does not
C______________________________ depend on any zonal or generic parameter.
C______________________________ So, we have finished dirty work

        NPTOT=IC1-1

      ELSE IF (NFTIP.EQ.22) THEN

C______________________________ Quasi-linear retention curve

        CONTINUE

      ELSE IF (NFTIP.EQ.23) THEN

C______________________________ Broadbridge & White retention curve 

        CONTINUE                ! Coming soon...
      
      ELSE IF (NFTIP.EQ.24) THEN

C______________________________ Van Genuchten retention curve
C______________________________ Defines formula parameters
         

        STOR=XPARAM(1)*PARAMC              ! (CF*FT*Pz) Ss
        CFFTSS=XPARAM(1)                   ! (CF*FT) Ss
        CFFTPOR=XPARAM(2)                  ! (CF*FT) Por 
        POR=XPARAM(2)*XPARAM(3)            ! (CF*FT*Pz) Por 
        XLAMBDA=PRGC(NFNLPRG(1,NFNLPARAM)) ! Van Genuchten exponent
        SMIN=PRGC(NFNLPRG(2,NFNLPARAM))    ! Residual saturation
        SMAX=PRGC(NFNLPRG(3,NFNLPARAM))    ! Maximum saturation degree
        PCAP=PRGC(NFNLPRG(4,NFNLPARAM))    ! Bubble pressure

C______________________________ Computation of non-linear function value
C______________________________ It is stored in DERIV(1). Only if zonal parameter
C______________________________ is going to be estimated. In this case Tij(zone k)

        IC1=1                   ! Estimated parameter index counter
        IPOS(IC1)=IVPAR(IND1)   ! IND1 is the zonal parameter zone number

        IF (IPOS(IC1).NE.0) THEN

C______________________________Computations for saturated state

          IF (X.GE.0D0)THEN
   
            SATW=SMAX

C______________________________ Computations for unsaturated state
C______________________________ Computes the saturation degree and its first
C______________________________ and second order derivative with respect to
C______________________________ the state variable.

          ELSE
   
            RELP=(-X/PCAP)**(1D0/(1D0-XLAMBDA))   ! Auxiliar scalar
            SATW=SMIN+(SMAX-SMIN)*(1D0+RELP)**(-XLAMBDA)    

          ENDIF

C______________________________ In this case, non-linear function depends on 
C______________________________ the zonal parameter

          DERIV(IC1)=CFFTSS*SATW    
          IC1=IC1+1                  ! Updates counter

        END IF

        RELSUC=-1D0*X/PCAP           ! All cases
        DERSWDMU=SMAX-SMIN           ! All cases

C______________________________Effective saturation is used to estimate Smin
C______________________________and Smax

        IF (X.GE.0D0)THEN
          SATW=SMAX
        ELSE
          RELP=(-X/PCAP)**(1D0/(1D0-XLAMBDA))   ! Auxiliar scalar
          SATW=SMIN+(SMAX-SMIN)*(1D0+RELP)**(-XLAMBDA)
        ENDIF

        XMU=(SATW-SMIN)/(SMAX-SMIN)

C______________________________ First generic parameter: XLAMBDA

        IPOS(IC1)=IVPAR(INPRGC+NFNLPRG(1,NFNLPARAM))
        IF (IPOS(IC1).NE.0) THEN

          IF (X.LT.0.D0) THEN
            U=1D0+RELSUC**(1D0/(1D0-XLAMBDA))       ! Sat. degree base (aux)

            DERUDERLAM=RELSUC**(1D0/(1D0-XLAMBDA))*DLOG(RELSUC)
     ;                /((1D0-XLAMBDA)*(1D0-XLAMBDA))

            DERMUDERLAM=-1D0*XLAMBDA*U**(-1D0*XLAMBDA-1D0)
     ;                 *DERUDERLAM-U**(-1D0*XLAMBDA)*DLOG(U)
          
            U=RELSUC**(1D0/(1D0-XLAMBDA))
            DUDLAM=U*DLOG(RELSUC)/((1D0-XLAMBDA)**2)
            TERM1=DLOG(U)-DLOG(1D0+U)+DUDLAM*((XLAMBDA-U)/(U*(1D0+U)))
            TERM1=TERM1*XLAMBDA-1D0/(XLAMBDA-1D0)
            COM=(SMAX-SMIN)*U**XLAMBDA*(1D0+U)**(-1D0*XLAMBDA-1D0)
     ;         /PCAP/(XLAMBDA-1D0)
            DERLAM=TERM1*COM
            DERIV(IC1)=STOR*DERSWDMU*DERMUDERLAM-POR*DERLAM
          ELSE
            DERIV(IC1)=0.D0
          END IF
          IC1=IC1+1
        END IF


C__________________________ Second generic parameter: minimum saturation degree

         
        IPOS(IC1)=IVPAR(INPRGC+NFNLPRG(2,NFNLPARAM))
        IF (IPOS(IC1).NE.0) THEN

          DERSWDERSMIN=1.D0-XMU
          U=RELSUC**(1D0/(1D0-XLAMBDA))
          DERSMIN=XLAMBDA/(1D0-XLAMBDA) *(1D0+U)**(-1D0*XLAMBDA-1D0)
     ;           *U**XLAMBDA/PCAP
          DERIV(IC1)=STOR*DERSWDERSMIN-POR*DERSMIN
          IC1=IC1+1
        END IF

C___________________________ Third generic parameter: maximum saturation degree

        IPOS(IC1)=IVPAR(INPRGC+NFNLPRG(3,NFNLPARAM))
        IF (IPOS(IC1).NE.0) THEN

          DERSWDERSMAX=XMU
          U=RELSUC**(1D0/(1D0-XLAMBDA))

          DERSMAX=XLAMBDA/(XLAMBDA-1D0) *(1D0+U)**(-1D0*XLAMBDA-1D0)
     ;           *U**XLAMBDA/PCAP
          DERIV(IC1)=STOR*DERSWDERSMAX-POR*DERSMAX
          IC1=IC1+1
        END IF

C______________________________ Fourth generic parameter: bubble pressure  

        IPOS(IC1)=IVPAR(INPRGC+NFNLPRG(4,NFNLPARAM))
        IF (IPOS(IC1).NE.0) THEN

          DMUDH=(RELSUC/PCAP)*(XLAMBDA/(1.D0-XLAMBDA))*
     ;           RELSUC**(XLAMBDA/(1D0-XLAMBDA))*
     ;           (1.D0+RELSUC**(1.D0/(1D0-XLAMBDA)))**(-XLAMBDA-1.D0)
         

          A=(SMAX-SMIN)*XLAMBDA/(XLAMBDA-1.0D0)
          U=RELSUC**(1D0/(1D0-XLAMBDA))  
          DERUDERH=(1D0/(1D0-XLAMBDA))*U**XLAMBDA*(-1D0*RELSUC/PCAP)
          AUX=(1D0+U)**((-1D0)*(XLAMBDA+1D0))*U**XLAMBDA/PCAP
          AUX2=(XLAMBDA-U)*DERUDERH/(U*(1D0+U))-1D0/PCAP
          DERH=A*AUX*AUX2

          DERIV(IC1)=STOR*DERSWDMU*DMUDH-POR*DERH
          IC1=IC1+1
        END IF
C______________________________ This kind of non-linear function also depends on other
C______________________________ non-generic parameters, which are not the zonal parameter
C______________________________ so an additional term must be added (in this case, porosity)

        IPOS(IC1)=IVPAR(IND2)         ! IND2 indicates the number of the porosity zone
        IF(IPOS(IC1).NE.0) THEN       ! Computation of the first derivative of sat. degree
                                      ! with respect to cap. pressure

          IF (X.GE.0D0)THEN

            DYDX=0D0            ! Right vertical branch (not on figure A1.9.4)

C______________________________ Computations for unsaturated state
C______________________________ Computes the first and second order derivative 
C______________________________ with respect to the state variable.

          ELSE

            RELP=(-X/PCAP)**(1D0/(1D0-XLAMBDA))   !RELP and MULT are auxiliar var.
            XMULT=XLAMBDA*(SMIN-SMAX)/(1D0-XLAMBDA)/PCAP*RELP**XLAMBDA

            IF(IOCAP.EQ.0) THEN            ! Capacitive scheme
              U=RELSUC**(1D0/(1D0-XLAMBDA))
              DYDX=(SMAX-SMIN)*(XLAMBDA/(XLAMBDA-1D0))
     ;           *(1D0+U)**(-XLAMBDA-1D0)*U**XLAMBDA/PCAP
            ELSE                           ! Conservative scheme
              CONTINUE
            END IF
          END IF  

          DERIV(IC1)=CFFTPOR*DYDX
          IC1=IC1+1

        END IF

        NPTOT=IC1-1

      ELSE IF (NFTIP.EQ.25) THEN

C______________________________ Brooks&Corey retention curve 

        CONTINUE                      ! Coming soon...

      ELSE IF (NFTIP.EQ.41) THEN

C______________________________ Non-linear leakage (aquifer-river interaction)

        CONTINUE                      ! Coming soon...

      END IF

      IF(IFLAGS(3).NE.0) CALL IO_SUB('FUNNOLI_IN',1)

      RETURN
      END 
