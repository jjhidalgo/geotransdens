       SUBROUTINE FUNNOLI
     ;(DYDH     ,DYDC     
     ;,IOCAP    ,H
     ;,NFLAGS   ,NFNL     ,NFNLPARAM,NFTIP    ,NZPRG          
     ;,Y        ,IFLAGS   ,NFNLPRG  
     ;,PARACD   ,PRGC     ,XPARAM)  
*******************************************************************************
*
* PURPOSE  Computing a particular function in a specific point of its domain and
*          the first derivative respect to the state variable.
*
*
* DESCRIPTION  This subroutine includes functions concerning to transport and 
*              flow phaenomena. All of them are related to the calculation of a
*              physical parameter of the flow or transport equations as:
*                  
*                       CFPARAM*PARZON*FTIME*FNONLIN
*
*              where FNONFLIN is computed here.
*              Moreover, all of them are computed using the state variable of 
*              the corresponding equation (or a function of it) as the 
*              independent variable.
*              It deals with transmissivity, storage, retention curves and 
*              leakage non linear functions (Oct-98), but will be increased...
*              There's a particular case; user can distinguish between two 
*              parameterization schemes when using a retention curve:
*                - Capacitive: comp. Ss*(Sw+DERSATWDERPSI), where the last term 
*                              is computed analitically.
*                - Conservative: comp. Ss*(Sw+DERSATWDERPSI), where the last 
*                                term is computed numerically.
*
*
* EXTERNAL VARIABLES: ARRAYS
*
*  NFNLPRG                Generic parameter zone number for every nonlinear
*                         function
*  PARACD                 Agreement parameters
*  PRGC                   Array containing generic parameters (subset of PARC)
*  XPARAM                 Array containing porosity divided by storavity on the
*                         first component and the Porosity on the second one 
*
* EXTERNAL VARIABLES: SCALARS
*
*  DYDX                   Derivative of the function respect to the independent 
*                         variable at time=k+1
*  DYDXOLD                Derivative of the function respect to the independent 
*                         variable at time=k
*  IOCAP                  Parameterization scheme. If =0 --> capacitive
*                                                  else  --> conservative
*  X                      Independent variable (state variable at time=k+1)
*  XOLD                   Independent variable (state variable at time=k)
*  Y                      Non linear function to be evaluated at time=k+1
*  YOLD                   Non linear function to be evaluated at time=k
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
*
*  NZPRG                  Total number of generic parameter zones               
*
*
* HYSTORY: First coding: German Galarza (May-98)
*
*          Re-coding: Andres Alcolea (Oct-98)
*                     Christen Knudby (Oct-98) ! The fifth element 
*                     Agustin Medina (Oct-98)  
*
********************************************************************************                                                                       

      IMPLICIT NONE 

C  EXTERNAL VARIABLES: SCALARS
      INTEGER*4 IOCAP,NFLAGS,NFNL,NFNLPARAM,NFTIP,NZPRG

      REAL*8 DYDH,DYDC,H,Y

C  EXTERNAL VARIABLES: ARRAYS
      INTEGER*4 IFLAGS(NFLAGS) ,NFNLPRG(8,NFNL)
      REAL*8 PARACD(3,NFNL),PRGC(NZPRG),XPARAM(8)

C INTERNAL VARIABLES: SCALARS
      INTEGER*4 ITMAX,I
      REAL*8 XBASE,DELTA,X1,TMIN,E,D,B,A,C1,THETAN,THETAS,XLAMS,XMUMIN
     ;      ,TOLER,XMU,PSI,XMUINI,PSIINI,XMUEND, PSIEND,SLOPE,ALFA
     ;      ,XMUREL,XLAMBDA,SMIN,SMAX,PCAP,RELP,PSIMIN,DEXP,XMULT
     ;      ,SATW,DERSATW,SDELTA,RADDELTA,SDPOT,K1,K2,EX,S,DKDS,SATEF
     ;      ,RAD,DSQRT,SPOT,PPAR,SM,FACTOR,XMOICAP,CAPPRE,DER2SATW
     ;      ,G,RELPOR,POR,THETA,CAPSUC,HEAD,DERHEAD,QMAX,DEPTH,X2
     ;      ,PREF,AUX,AGREE,BETA

         
      IF(IFLAGS(3).EQ.1) CALL IO_SUB('FUNNOLI',0)
      
      DYDC=0D0

      IF (NFTIP.EQ.1) THEN

C______________________________ Unconfined aquifers
C______________________________ Asymptotic free aquifer transmissivity
C______________________________ X= Current head level
C______________________________ Y= Saturated thickness

        XBASE=XPARAM(1)                 ! Bottom of the aquifer
        DELTA=PARACD(1,NFNLPARAM)       ! Agreement parameter
         
        IF (H.GE.(XBASE+DELTA)) THEN    ! Linear zone (right)
         
          Y=(H-XBASE)
          DYDH=1D0
          
        ELSE                            ! Exponential zone (left)

          Y=DELTA*DEXP((H-XBASE-DELTA)/DELTA)
          DYDH=Y/DELTA

        ENDIF


      ELSE IF (NFTIP.EQ.2) THEN

C______________________________ Unconfined aquifer
C______________________________ Parallel free aquifer transmissivity
C______________________________ X= Current head level
C______________________________ Y= Saturated thickness

        XBASE=XPARAM(1)                  ! Bottom of the aquifer (Card B1.2)
        DELTA=PARACD(1,NFNLPARAM)          ! 0.5*Length of the transicion zone
        X1=(H-XBASE)                     ! Relative head level

        IF (X1.GT.DELTA)THEN             ! Linear zone (right)

          Y=(H-XBASE)
          DYDH=1D0

        ELSE IF(X1.GE.-DELTA)THEN        ! Transition zone (agreement) 

          TMIN=PRGC(NFNLPRG(1,NFNLPARAM))       ! Treshold value for B(h)
          E=DELTA

C______________________________ Computes the coeff. of the agreement 

          D=TMIN/(E*(DEXP(E)+DEXP(-E))-(DEXP(E)-DEXP(-E)))
          B=(1D0-D*(DEXP(E)+DEXP(-E)))/2D0
          A=(D*DEXP(-E)+B)/2D0/E           
          C1=TMIN-D*DEXP(-E)+B*E-A*E*E

          Y=A*X1*X1+B*X1+C1+D*DEXP(X1)
          DYDH=2*A*X1+B+D*DEXP(X1)

        ELSE

C______________________________ Constant zone (left)

          TMIN=PRGC(NFNLPRG(1,NFNLPARAM))       ! Treshold value for B(h)
          Y=TMIN
          DYDH=0D0

        ENDIF


      ELSE IF (NFTIP.EQ.3) THEN

C______________________________Broadbridge and White relative hyd. conductivity 

C______________________________Assigns some useful parameters

        THETAN=PRGC(NFNLPRG(1,NFNLPARAM))   ! Moisture (non-sat. state)
        THETAS=PRGC(NFNLPRG(2,NFNLPARAM))   ! Moisture (sat. state)
        C1=PRGC(NFNLPRG(3,NFNLPARAM))        ! Parameter C
        XLAMS=PRGC(NFNLPRG(4,NFNLPARAM))    ! Capillary lenght
        XMUMIN=PARACD(1,NFNLPARAM)          ! Minimum value of specific sat.
        TOLER=PARACD(2,NFNLPARAM)           ! Convergence criteria
        DELTA=PARACD(3,NFNLPARAM)           ! Agreement parameter
        ITMAX=10000                         ! Maximum number od iterations

        IF (H.GE.0D0) THEN

C______________________________Saturated state

          Y=1D0
          DYDH=0D0

        ELSE

C______________________________Non-saturated. It performs the calculation blocks
C______________________________The first one concerns to compute the water 
C______________________________content related to current pressure head. 
C______________________________The second one computes the unsaturated relative
C______________________________hydraulic conductivity and, at the third one 
C______________________________computes its derivative with respect to p. head.


C______________________________First block. Iterative method. Point allocation


C______________________________Checks the coherence of the formulation

          PSIMIN=-XLAMS*((1D0/C1)*DLOG((C1-XMUMIN)/(XMUMIN*(C1-1D0)))
     ;                      +(1D0-XMUMIN)/XMUMIN)

          IF (PSIMIN.GE.H) THEN
            WRITE(*,10) NFNLPARAM 
            STOP
          END IF


          I=-14

          DO WHILE (I.LE.-1)

            XMU=10D0**I
            PSI=-XLAMS*((1D0/C1)*DLOG((C1-XMU)/(XMU*(C1-1D0)))
     ;                  +(1D0-XMU)/XMU)      

C______________________________Solution has been allocated

            IF(PSI.GE.H . AND. PSIMIN.LT.H) THEN

              XMUINI=XMUMIN
              PSIINI=PSIMIN
              XMUEND=XMU
              PSIEND=PSI
              GOTO 100
            
            ELSE

C______________________________The process must continue

              XMUMIN=XMU
              PSIMIN=PSI

            END IF


          END DO

C______________________________Once the solution is allocated, begins iteration
C______________________________cycle between XMUINI and XMUEND.

  100      I=0

          DO WHILE (I.LE.ITMAX)

            I=I+1
            SLOPE=(PSIEND-PSIINI)/(XMUEND-XMUINI)
            XMU=(H-PSIINI)/SLOPE+XMUINI
            PSI=-XLAMS*((1D0/C1)*DLOG((C1-XMU)/(XMU*(C1-1D0)))
     ;                  +(1D0-XMU)/XMU)

C______________________________Checks tolerance

            IF (DABS(PSI-H).LE.TOLER) THEN
              I=ITMAX
            END IF

          END DO            

        END IF

        IF (DABS(PSI-H).GT.TOLER) THEN
          WRITE(*,20) 
          STOP
        END IF


C______________________________Once MU has been found, let's compute Kr and its
C______________________________derivative with respect to pressure head.
C______________________________There are two main zones:the formula and that of
C______________________________the agreement function, as a sum of a parabolic
C______________________________function and a exponential one.

        IF (XMU.LE.1D0-DELTA) THEN

C______________________________Broadbridge & White non-linear function

          Y=XMU*XMU*(C1-1D0)/(C1-XMU)
          DYDH=(C1-1D0)*XMU*XMU*XMU*(2D0*C1-XMU)/(XLAMS*C1*(C1-XMU))

        ELSE IF (XMU.GT.1D0) THEN

C______________________________What a piece of shit!!!!

          WRITE(*,*)'PROBLEMS WITH BROADBRIDGE AND WHITE:MU>1.I STOP'
          STOP

        ELSE

C______________________________Agreement function. Y=AX**2+BX+D+E*EXP(X)

          ALFA=1D0-(C1-1D0)*(1D0-DELTA)*(1D0-DELTA)/(C1+DELTA-1D0)
          BETA=(1D0-DELTA)*(2D0*C1+DELTA-1D0)*(C1-1D0)/
     ;         ((C1+DELTA-1D0)*(C1+DELTA-1D0))

          E=(2D0*ALFA-DELTA*BETA)/
     ;      (2D0*(DEXP(DELTA)-1D0-DELTA)-DELTA*(DEXP(DELTA)-1D0))
          D=-1D0*E          
          B=BETA-E
          A=(ALFA-E*(DEXP(DELTA)-1D0-DELTA)-BETA*DELTA)/(DELTA*DELTA)

          XMUREL=XMU+DELTA-1D0
          Y=1D0-ALFA+A*XMUREL*XMUREL+B*XMUREL+D+E*DEXP(XMUREL)
          DYDH=(2D0*A*XMUREL+B+E*DEXP(XMUREL))*XMU*XMU*(C1-XMU)/
     ;         (XLAMS*C1)

        END IF

      ELSE IF (NFTIP.EQ.4) THEN

C______________________________Relative hydraulic conductivity of Van Genuchten

C______________________________Identifies formula parameters

        XLAMBDA=PRGC(NFNLPRG(1,NFNLPARAM))  ! Van Genuchten exponent
        SMIN=PRGC(NFNLPRG(2,NFNLPARAM))     ! Minimum saturation degree
        SMAX=PRGC(NFNLPRG(3,NFNLPARAM))     ! Maximum saturation degree

        PCAP=PRGC(NFNLPRG(4,NFNLPARAM)) ! Capillary head (???)

        DELTA=PARACD(1,NFNLPARAM)           ! Length of agreement func. int.

        IF (H.GE.0D0) THEN

C______________________________ Saturated state

          Y=1D0                               ! Hyd. Cond. = sat. hyd. cond.
          DYDH=0D0

        ELSE                                  

C______________________________ Unsaturated state.
C______________________________ Computes some auxiliar terms
C______________________________ and Sw and dSw/dh

          RELP=(-H/PCAP)**(1D0/(1D0-XLAMBDA))          
          XMULT=XLAMBDA*(SMAX-SMIN)/(1D0-XLAMBDA)/PCAP*RELP**XLAMBDA

          SATW=SMIN+(SMAX-SMIN)*(1D0+RELP)**(-XLAMBDA) ! Sat. degree
          DERSATW=XMULT*(1D0+RELP)**(-XLAMBDA-1D0)     ! First order deriv.


          IF (SATW.GT.(SMAX-DELTA))THEN    

C______________________________Zone of agreement function. The agreement func.
C______________________________is the sum of a parabole and an exponential term

            SDELTA=1D0-DELTA/(SMAX-SMIN)                 ! Auxiliar term
            RADDELTA=SQRT(SDELTA)                        ! Auxiliar term
            SDPOT=1D0-SDELTA**(1D0/XLAMBDA)              ! Auxiliar term

C______________________________Computes some auxiliar terms:
C______________________________K1 is K (hyd. cond.) for SATW=SMAX-DELTA and
C______________________________K2 is the deriv. of K with respect to SATW when
C______________________________SATW=SMAX-DELTA

            K1=RADDELTA*(1D0-SDPOT**XLAMBDA)**2D0     
            K2=(1D0/(2D0*RADDELTA)*(1D0-SDPOT**XLAMBDA)**2D0+
     ;         2D0*RADDELTA*(1D0-SDPOT**XLAMBDA)*SDPOT**(XLAMBDA-1D0)
     ;         *SDELTA**((1D0-XLAMBDA)/XLAMBDA))/(SMAX-SMIN) 


            EX=DEXP(DELTA)
            D=(1D0-DELTA/2D0*K2-K1)/(EX*(1D0-DELTA/2D0)-DELTA/2D0-1D0)
            C1=K1-D
            B=K2-D
            A=(D*(1D0-EX)-K2)/2D0/DELTA

            S=SATW+DELTA-SMAX             ! Substitution
            Y=A*S**2D0+B*S+C1+D*DEXP(S)    ! Value of K
            DKDS=2D0*A*S+B+D*DEXP(S)      ! Derivative of K with respect to S
            DYDH=DKDS*DERSATW             ! Derivative of K with respect to X

          ELSE                            

C______________________________ Zone of van Genuchten function

            SATEF=(SATW-SMIN)/(SMAX-SMIN)        ! Auxiliar term
            RAD=DSQRT(SATEF)                     ! Auxiliar term
            SPOT=1D0-SATEF**(1D0/XLAMBDA)        ! Auxiliar term  

            Y=RAD*(1-SPOT**XLAMBDA)**2D0
            DKDS=(1D0/(2D0*RAD)*(1D0-SPOT**XLAMBDA)**2D0+
     ;      2D0*RAD*(1D0-SPOT**XLAMBDA)*SPOT**(XLAMBDA-1D0)*
     ;      SATEF**((1D0-XLAMBDA)/XLAMBDA))/(SMAX-SMIN)
            DYDH=DKDS*DERSATW

          ENDIF

        ENDIF

      ELSE IF (NFTIP.EQ.5) THEN

C______________________________Brooks & Corey relative hydraulic conductivity 

C______________________________Identifies some useful parameters

        PCAP=PRGC(NFNLPRG(1,NFNLPARAM))   ! Capillar pressure
        PPAR=PRGC(NFNLPRG(2,NFNLPARAM))   ! Parameter P > 3
        DELTA=PARACD(1,NFNLPARAM)         ! Length of agreement func. int.

        IF(H.GE.0D0) THEN

C______________________________Saturated state

          Y=1D0
          DYDH=0D0

        ELSE

C______________________________Non-saturated state

          XMU=(-1D0*PCAP/H)**(2D0/(PPAR-3D0))

          IF (XMU.LE.1D0-DELTA) THEN

C______________________________Brooks & Corey rel. perm. function

            Y=XMU**PPAR


            DYDH=2D0*PPAR*XMU**(PPAR-1D0)*(-PCAP)*
     ;           (-1D0*PCAP/H)**((5D0-PPAR)/(PPAR-3D0))/(PPAR-3D0)

          ELSE

C______________________________Agreement. Parabolic-exponential

            ALFA=(1D0-DELTA)**PPAR
            BETA=PPAR*(1D0-DELTA)**(PPAR-1)
            D=(2D0*ALFA-BETA*DELTA)/
     ;        (2D0*(DEXP(DELTA)-1D0-DELTA)-DELTA*(DEXP(DELTA)-1D0))
            C1=-1D0*D
            B=BETA-D
            A=(ALFA-E*(DEXP(DELTA)-1D0-DELTA)-BETA*DELTA)/(DELTA*DELTA)

            XMUREL=XMU+DELTA-1D0

            Y=1D0-ALFA+A*XMUREL*XMUREL+B*XMUREL+C1+D*DEXP(XMUREL)
            DYDH=(2D0*A*XMUREL+B+D*DEXP(XMUREL))*2D0*(-1D0*PCAP)*
     ;           (-1D0*PCAP/H)**((5D0-PPAR)/(PPAR-3D0))/(PPAR-3D0)

          END IF

        END IF

      ELSE IF (NFTIP.EQ.21) THEN

C______________________________ Two layer storage non-linear function
C______________________________ Free aquifers
C______________________________ X= Current head level
C______________________________ Y= Factor such us S=STGC*FACTOR

        SM=PARACD(1,NFNLPARAM)          ! Minimum relative value of STGC
        DELTA=PARACD(2,NFNLPARAM)       ! Agreement parameter

        IF (H.GE.(XBASE+DELTA)) THEN    

C______________________________ Constant zone (right)

          Y=1D0
          DYDH=0D0

        ELSE IF(H.LT.(XBASE-DELTA)) THEN        

C______________________________ Constant zone (left)

          Y=SM
          DYDH=0D0

        ELSE                                     

C______________________________ Transition zone

          X1=(H-XBASE)/DELTA
          IF(X1.LT.0D0) FACTOR=-1D0
          Y=(1D0-SM)*(3D0*X1-FACTOR*DABS(X1)**3)/4+(1+SM)/2
          DYDH=(1D0-SM)*3D0*(1D0-X1*X1)/(4*DELTA)

        ENDIF

      ELSE IF (NFTIP.EQ.22) THEN

C______________________________Quasi-linear retention curve
C______________________________Defines some useful parameters

        XMOICAP=PRGC(NFNLPRG(1,NFNLPARAM))   ! Moisture capacity
        SMIN=PRGC(NFNLPRG(2,NFNLPARAM))      ! Residual saturation
        SMAX=PRGC(NFNLPRG(3,NFNLPARAM))      ! Maximum saturation
        CAPPRE=PRGC(NFNLPRG(4,NFNLPARAM))    ! Capillar pressure
        DELTA=PARACD(1,NFNLPARAM)            ! Agreement parameter         
        PREF=CAPPRE-(SMAX-SMIN)-DELTA/2D0    ! Reference point 

C______________________________For head levels that don't verifies
C______________________________sat./non-sat. state, a previous operation
C______________________________is done: porosity and specific storage coeff.
C______________________________have been computed and stored on XPARAM.
C______________________________Terms: SATW= Saturation degree
C                                     DERSATW= Deriv. of SATW with respect to
C                                              state variable
C                                     DER2SATW= Second deriv.
C                                     Y= Term multypling to INC(PSY)
C                                     DYDX= Its deriv with respect to state v.


C______________________________Covering all zones from left to right
C______________________________Left zone: Non-saturated state

        IF (H.LE.PREF) THEN

          SATW=SMIN

          IF(IOCAP.EQ.0) THEN            ! Capacitive scheme

            DERSATW=0D0
          
          ELSE                           ! Conservative scheme

            CONTINUE

          END IF


          Y=SATW
          DYDH=0D0

C______________________________Parabolic agreement (left)

        ELSE IF (H.GT.PREF.AND.H.LE.(PREF+DELTA)) THEN

          X1=H-PREF                             ! Relative head
          SATW=XMOICAP*X1*X1/(2D0*DELTA)+SMIN

          IF(IOCAP.EQ.0) THEN            ! Capacitive scheme

            DERSATW=X1*XMOICAP/DELTA
            DER2SATW=XMOICAP/DELTA

          ELSE                           ! Conservative scheme

            CONTINUE

          END IF

          Y=SATW+XPARAM(1)*DERSATW
          DYDH=DERSATW+XPARAM(1)*DER2SATW

C______________________________Linear zone (middle)

        ELSE IF (H.GT.(PREF+DELTA).AND.H.LT.(CAPPRE-DELTA)) THEN
        
          SATW=XMOICAP*(H-PREF)+SMIN-XMOICAP*DELTA/2D0

          IF(IOCAP.EQ.0) THEN            ! Capacitive scheme

            DERSATW=XMOICAP

          ELSE                           ! Conservative scheme

            CONTINUE

          END IF

          Y=SATW+XPARAM(1)*DERSATW
          DYDH=DERSATW

C______________________________Parabolic agreement (right)

        ELSE IF (H.GT.(CAPPRE-DELTA).AND.H.LT.CAPPRE) THEN

C______________________________Computes agreement coeff.

          X1=H-(CAPPRE-DELTA)
          G=SMAX-XMOICAP*DELTA
          D=-2D0*G/
     ;      (DELTA*(DEXP(DELTA)-1D0)-2D0*(DEXP(DELTA)-1D0-DELTA))
          A=(G-D*(DEXP(DELTA)-1D0-DELTA))/(DELTA*DELTA)
          C1=-1D0*D
          B=XMOICAP-D
          SATW=A*X1*X1+B*X1+C1+D*DEXP(X1)+SMAX-XMOICAP*DELTA

          IF(IOCAP.EQ.0) THEN            ! Capacitive scheme

            DERSATW=2D0*A*X1+B+D*DEXP(X1)
            DER2SATW=2D0*A+D*DEXP(X1)

          ELSE                           ! Conservative scheme

            CONTINUE

          END IF

          Y=SATW+XPARAM(1)*DERSATW
          DYDH=DERSATW+XPARAM(1)*DER2SATW


C______________________________Right zone: Saturated state

        ELSE           

          SATW=SMAX

          IF(IOCAP.EQ.0) THEN            ! Capacitive scheme

            DERSATW=0D0

          ELSE                           ! Conservative scheme

            CONTINUE

          END IF

          Y=SATW
          DYDH=0D0

        END IF


      ELSE IF (NFTIP.EQ.23) THEN

C______________________________Broadbridge & White retention curve 

C______________________________Assigns some useful parameters

        THETAN=PRGC(NFNLPRG(1,NFNLPARAM))   ! Moisture (non-sat. state)
        THETAS=PRGC(NFNLPRG(2,NFNLPARAM))   ! Moisture (sat. state)
        C1=PRGC(NFNLPRG(3,NFNLPARAM))        ! Parameter C
        XLAMS=PRGC(NFNLPRG(4,NFNLPARAM))    ! Capillary lenght
        XMUMIN=PARACD(1,NFNLPARAM)          ! Minimum value of specific sat.
        TOLER=PARACD(2,NFNLPARAM)           ! Convergence criteria
        RELPOR=XPARAM(1)                    ! Porosity/Specific Storage
        POR=XPARAM(2)                       ! Porosity
        ITMAX=10000                         ! Maximum number od iterations

        IF (H.GE.0D0) THEN

C______________________________Saturated state

          Y=1D0
          DYDH=0D0

        ELSE

C______________________________Non-saturated. It performs the calculation blocks
C______________________________The first one concerns to compute the water 
C______________________________content related to current pressure head. 
C______________________________The second one computes the unsaturated relative
C______________________________hydraulic conductivity and, at the third one 
C______________________________computes its derivative with respect to p. head.
          

C______________________________First block. Iterative method. Point allocation


C______________________________Checks the coherence of the formulation

          PSIMIN=-XLAMS*((1D0/C1)*DLOG((C1-XMUMIN)/(XMUMIN*(C1-1D0)))
     ;                      +(1D0-XMUMIN)/XMUMIN)

          IF (PSIMIN.GE.H) THEN
            WRITE(*,10) NFNLPARAM 
            STOP
          END IF

   10     FORMAT('VALUE OF PARACD(2,',I5,') IS GREATER THAN ALLOWED')

          I=-14

          DO WHILE (I.LE.-1)

            XMU=10D0**I
            PSI=-XLAMS*((1D0/C1)*DLOG((C1-XMU)/(XMU*(C1-1D0)))
     ;                  +(1D0-XMU)/XMU)      

C______________________________Solution has been allocated

            IF (PSI.GE.H . AND. PSIMIN.LT.H) THEN

              XMUINI=XMUMIN
              PSIINI=PSIMIN
              XMUEND=XMU
              PSIEND=PSI
              GOTO 1000
            
            ELSE

C______________________________The process must continue

              XMUMIN=XMU
              PSIMIN=PSI

            END IF


          END DO

C______________________________Once the solution is allocated, begins iteration
C______________________________cycle between XMUINI and XMUEND.

 1000     I=0

          DO WHILE (I.LE.ITMAX)

            I=I+1
            SLOPE=(PSIEND-PSIINI)/(XMUEND-XMUINI)
            XMU=(H-PSIINI)/SLOPE+XMUINI
            PSI=-XLAMS*((1D0/C1)*DLOG((C1-XMU)/(XMU*(C1-1D0)))
     ;                  +(1D0-XMU)/XMU)

C______________________________Checks tolerance

            IF (DABS(PSI-H).LE.TOLER) THEN
              I=ITMAX
            END IF

          END DO            

        END IF

        IF (DABS(PSI-H).GT.TOLER) THEN
          WRITE(*,20) 
          STOP
        END IF

   20   FORMAT('TOLERANCE QUITE SMALL. I"M SORRY, I MUST STOP')

C______________________________Once MU has been forund, let's compute the 
C______________________________derivatives of Sw with respect to PSI, knowing 
C______________________________that Sw=MOISTURE/POROSITY

        DELTA=THETAS-THETAN
        THETA=XMU*DELTA+THETAN
        SATW=THETA/POR

        IF (IOCAP.EQ.0) THEN

C______________________________Capacitive scheme

          DERSATW=DELTA*XMU*XMU*(C1-XMU)/(XLAMS*C1*POR)
          DER2SATW=(2D0*C1-3D0*XMU)*XMU*XMU*XMU*DELTA*(C1-XMU)
     ;            /(POR*XLAMS*XLAMS*C1*C1)

        ELSE

C______________________________Conservative scheme

          CONTINUE

        END IF

        Y=SATW+RELPOR*DERSATW
        DYDH=DERSATW+RELPOR*DER2SATW


      ELSE IF (NFTIP.EQ.24) THEN

C______________________________Van Genuchten retention curve
C______________________________Defines formula parameters
         
        RELPOR=XPARAM(1)                   ! Porosity divided by storativity
        XLAMBDA=PRGC(NFNLPRG(1,NFNLPARAM)) ! Van Genuchten exponent
        SMIN=PRGC(NFNLPRG(2,NFNLPARAM))    ! Residual saturation
        SMAX=PRGC(NFNLPRG(3,NFNLPARAM))    ! Maximum saturation degree
        CAPSUC=PRGC(NFNLPRG(4,NFNLPARAM))  ! Capillary suction

C______________________________Computations for saturated state

        IF (H.GE.0D0)THEN

          Y=SMAX              ! Only the elastic term is active
          DYDH=0D0            ! Right vertical branch (not on figure A1.9.4)
          XPARAM(1)=SMAX      ! Saturation degree is returned in this array
          RETURN

C______________________________Computations for unsaturated state
C______________________________Computes the saturation degree and its first
C______________________________and second order derivative with respect to
C______________________________the state variable.

        ELSE

          RELP=(-H/PCAP)**(1D0/(1D0-XLAMBDA))
          XMULT=XLAMBDA*(SMAX-SMIN)/(1D0-XLAMBDA)/PCAP*RELP**XLAMBDA
           
          SATW=SMIN+(SMAX-SMIN)*(1D0+RELP)**(-XLAMBDA)    

          IF(IOCAP.EQ.0) THEN            ! Capacitive scheme

            DERSATW=XMULT*(1D0+RELP)**(-XLAMBDA-1D0)
            DER2SATW=-XMULT/(1D0-XLAMBDA)/PCAP*RELP**XLAMBDA*
     ;      (XLAMBDA/RELP*(1D0+RELP)**(-XLAMBDA-1D0)-
     ;      (XLAMBDA+1D0)*(1D0+RELP)**(-XLAMBDA-2D0))

          ELSE                           ! Conservative scheme

            CONTINUE

          END IF


          Y=RELPOR*DERSATW+SATW           ! Comments.....?
          DYDH=(RELPOR*DER2SATW+DERSATW)

          XPARAM(1)=SATW    ! XPARAM is assigned the value of the sat. deg.

        ENDIF

      ELSE IF (NFTIP.EQ.25) THEN

C______________________________Brooks&Corey retention curve 

C______________________________Identifies some useful parameters

        PCAP=PRGC(NFNLPRG(1,NFNLPARAM))     ! Capillar pressure
        PPAR=PRGC(NFNLPRG(2,NFNLPARAM))     ! Parameter P > 3
        THETAN=PRGC(NFNLPRG(3,NFNLPARAM))   ! Moisture (non-sat.)
        THETAS=PRGC(NFNLPRG(4,NFNLPARAM))   ! Moisture (sat.)
        RELPOR=XPARAM(1)                    ! Porosity/Specific Storage
        POR=XPARAM(2)                       ! Porosity

        IF (H.GE.0D0) THEN

C______________________________Saturated state

          Y=1D0
          DYDH=0D0

        ELSE

C______________________________Non-saturated state

          XMU=(-1D0*PCAP/H)**(2D0/(PPAR-3D0))
          THETA=THETAN+(THETAS-THETAN)*XMU          
          SATW=THETA/POR

C______________________________Computes the derivatives od saturation degree 

          IF (IOCAP.EQ.0) THEN

C______________________________Capacitive scheme

            DERSATW=(THETAS-THETAN)*2D0*(-1D0*PCAP)*
     ;             (-1D0*PCAP/H)**((5D0-PPAR)/(PPAR-3D0))/
     ;             (POR*(PPAR-3D0))

            DER2SATW=(THETAS-THETAN)*2D0*PCAP*PCAP*(5D0-PPAR)*
     ;             (-1D0*PCAP/H)**((8D0-2D0*PPAR)/(PPAR-3D0))/
     ;             ((PPAR-3D0)*(PPAR-3D0)*POR)

          ELSE

C______________________________Conservative scheme

            CONTINUE

          END IF

C______________________________Finally computes the output values

          Y=SATW+RELPOR*DERSATW
          DYDH=DERSATW+RELPOR*DER2SATW

        END IF

      ELSE IF (NFTIP.EQ.41.OR.NFTIP.EQ.42) THEN

C______________________________Sink-source limited condition expressed as 
C______________________________function of the maximum sink/source
C______________________________Input description:

        HEAD=XPARAM(1)         ! Input: Prescribed head value
        DERHEAD=XPARAM(2)      ! Input: Derivative of p. head with r. t. head
        ALFA=XPARAM(3)         ! Input: CF*ALFX*F(t)   (linear leakage)

        IF (NFTIP.EQ.41) THEN

          QMAX=PRGC(NFNLPRG(1,NFNLPARAM)) ! Maximum value of leached flow
          DEPTH=QMAX/ALFA
            
        ELSE

          DEPTH=PRGC(NFNLPRG(1,NFNLPARAM)) 

        END IF

        DELTA=PARACD(1,NFNLPARAM) ! Half lenght of the agreement parabole
        X1=HEAD-DEPTH-DELTA
        X2=HEAD-DEPTH+DELTA

        IF (H.LE.X1) THEN        ! Constant flow (left zone)

          XPARAM(4)=ALFA*DEPTH
          XPARAM(5)=0D0        
          Y=0D0
          DYDH=0D0

        ELSE IF (H.GE.X2) THEN   ! Linear (right zone)

          XPARAM(4)=0D0        
          XPARAM(5)=0D0        
          Y=1D0
          DYDH=0D0

        ELSE                     ! Agreement zone

          AGREE=DEPTH+DELTA
          AUX=1.D0/(4D0*DELTA)

          Y=AUX*( H-HEAD+2*AGREE )
          DYDH=AUX*( 1-DERHEAD)

          XPARAM(4)=ALFA*( DEPTH-AUX*AGREE*AGREE )
          XPARAM(5)=0D0

        END IF

      ENDIF

      IF(IFLAGS(3).EQ.1) CALL IO_SUB('FUNNOLI',1)

      RETURN

      END SUBROUTINE FUNNOLI
