      SUBROUTINE DERTRA 
     &          (AREA     ,BIBI     ,BUOYANCY ,CFPAREL  ,COORD
     &          ,DENSITY  ,DERH     ,DTIM     ,EPSFLU   ,FNT
     &          ,GP_COORD ,GRADLOC  ,HAUX1    ,HBASE    ,HCALAN
     &          ,HCALIT   ,IDIMBB   ,IDIMDENS ,IDIMDERH ,IDIMFNT
     &          ,IFLAGS   ,INDSSTR  ,INEW     ,INORPAR  ,INTI
     &          ,IOCAP    ,IOCTRA   ,IODENS   ,IODIM    ,IOFLLI
     &          ,IOFMLF   ,ISOZ     ,IVPAR    ,KXX      ,LDIM
     &          ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXPAREL  ,MAXPG
     &          ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &          ,NFTPAR   ,NINT     ,NPAREL   ,NPARF    ,NTYPAR
     &          ,NTYPEL   ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR
     &          ,NZTRA    ,PARACD   ,PARC     ,POINTWEIGHT
     &          ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,NPAR     ,IPOS
     &          ,DERIV)


*****************************************************************************
*
* PURPOSE
*
*     Compute the contribution of the derivatives of AFLU matrix with respect 
*     to transmissivities to the RHS of the derivatives of h  wrt parameters
*
* DESCRIPTION
*
*     Compute the contribution of the derivatives of AFLU matrix with respect
*     to transmissivities to the RHS of the derivatives of h  wrt parameters
*     More precisely, adds to the RHS the expression 
*
*                              d AFLU                 k+1             k
*                         ----------------   ( theta*h   + (1-theta)*h  )
*                         d Transmissivity
*
*
* EXTERNAL VARIABLES: ARRAYS
*
*  BIBI                   Array containing the product of interpolation         
*                         functions gradient, for a given element               
*  CFPAREL                Array containing node coefinient of element j         
*                         corresponding to INpar index parameter zone.          
*  DERH                   Nodal head derivatives with respect to estimated      
*                         flow parameters.                                      
*  FNT                    Array containing time functions values                
*  HAUX1                  Array containing HEADS, ponderated by THETAF          
*                         time factor                                           
*  HBASE                  Bottom level of the aquifer
*  HCALIT                 Computed heads at every node                          
*  HCALAN                 Head level at previous time                           
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  ISOZ                   Anisotropy of every transmissivity zone               
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LDIM                   Vector containing the dimension of each element       
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing the type of each element            
*  LXPAREL                Array containing zone numbers for a given             
*                         element parameter                                     
*  NFNLPAR                Vector containing non-linear function order           
*                         afecting every parameter at each zone.                
*  NFNLPRG                Generic parameter zone number for every nonlinear     
*                         function                                              
*  NFNLTIP                Type of non-linear function                           
*  NFTPAR                 Vector containing time function number at every       
*                         parameter zone                                        
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters                                            
*  PARACD                 Agreement parameters                                  
*  PARC                   Vector containing calculated values for all           
*                         parameters                                            
*
* INTERNAL VARIABLES: ARRAYS
*
*  DERIV                  Derivatives of the current parameter with respect
*                         to its zonal value plus its generic parameters
*  IPOS                   Location order for estimated parameters
*  XPARAM                 Array used to store auxiliar values for the 
*                         computation of the non-linear function
*
* EXTERNAL VARIABLES: SCALARS
*
*  DTIM                   Value to compute the time function value at the 
*                         current time (k+theta)
*  EPSFLU                 Time weighting parameter for nonlinear flow problems  
*  IDIMBB                 Used to dimension array BIBI                          
*  IDIMFNT                First dimension of array FNT, it coincides with       
*                         NFNT if the latter is not zero                        
*  INDSSTR                Current problem state. If 1, transient state,         
*                         if 0, steady-state                                    
*  INEW                   Index to locate the RHS of the derivatives of head 
*                         at current time with respect to parameters inside 
*                         array DERH
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  IOCTRA                 Option that defines the way of computing the          
*                         nonlinear function term of the transmissivity value   
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  IOFMLF                 Flow Formulation number                               
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFLAGS                 Maximum number of allowed flags                       
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NINT                   Number of observation times                           
*  NPAREL                 Number of element parameters in current problem       
*  NPARF                  Number of transient parameters to be estimated        
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  NZTRA                  Number of transmissivity zones                        
*
* INTERNAL VARIABLES: SCALARS
*
*  ISOT1                  Transmissivity anisotropy of the current element
*  ISOTTR                 Maximum possible anisotropy for the current element
*  IZON                   Transmissivity zone number of the current element
*  JJ                     Location of transmissivity values of the current 
*                         transmissivity zone in arrays PARC, NFNLPAR, etc
*  L                      Current element
*  LD                     Dimension of the current element
*  LTYPE1                 Type of the current element
*  M                      Counter location in array BIBI
*  NCNF                   Nonlinear function number of the current 
*                         transmissivity zone 
*  NNUD                   Number of nodes of the current element                
*  NPTOT                  Total number of parameters to be estimated in the 
*                         current transmissivity zone (includes zonal plus 
*                         generic parameters)
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  DER_PARAM              Computes the derivatives of the parameter with 
*                         respect to its zonal value plus its generic parameters
*
* HISTORY
*
*     SCR      ?-1997     First coding
*     AMS     12-1998     Full modification
*
*****************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMBB   ,IDIMDENS ,IDIMDERH ,IDIMFNT  ,INDSSTR
     &          ,INEW     ,INTI     ,IOCAP    ,IOCTRA   ,IODENS
     &          ,IODIM    ,IOFLLI   ,IOFMLF   ,LMXNDL   ,MAXPG
     &          ,NFLAGS   ,NFNL     ,NINT     ,NPAREL   ,NPARF
     &          ,NTYPAR   ,NTYPEL   ,NUMEL    ,NUMNP    ,NZPAR
     &          ,NZTRA    ,IDIMWGT  ,IPNT_PAR ,NPAR

      REAL*8::DTIM   ,EPSFLU


      INTEGER*4::IFLAGS(NFLAGS),INORPAR(NTYPAR)   ,ISOZ(NZTRA)
     &          ,IVPAR(NZPAR)  ,KXX(LMXNDL,NUMEL) ,LDIM(NUMEL)
     &          ,LNNDEL(NUMEL) ,LTYPE(NUMEL)      ,LXPAREL(NUMEL,NPAREL)
     &          ,NFNLPAR(NZPAR),NFNLPRG(8,NFNL)   ,NFTPAR(NZPAR)
     &          ,NFNLTIP(NFNL) ,NZONE_PAR(NTYPAR)

      
       REAL*8::AREA(NUMEL)                  ,BIBI(IDIMBB,NUMEL)
     &        ,BUOYANCY(IODIM,LMXNDL,NUMEL) ,CFPAREL(NUMEL,NPAREL)
     &        ,COORD(NUMNP,3)               ,DENSITY(IDIMDENS)
     &        ,DERH(NUMNP,NPARF,IDIMDERH)   ,FNT(IDIMFNT,NINT)
     &        ,GRADLOC(IODIM,LMXNDL,MAXPG)  ,GP_COORD(6,8,IODIM)
     &        ,HAUX1(NUMNP)                 ,HBASE(NUMEL)
     &        ,HCALAN(NUMNP)                ,HCALIT(NUMNP)
     &        ,PARACD(3,NFNL)               ,PARC(NZPAR)
     &        ,POINTWEIGHT(MAXPG,NTYPEL)    ,WGT_PAR(IDIMWGT)  

C------------------------- Internal

      INTEGER*4::I1     ,I2     ,IC     ,IP     ,ISOT1  ,ISOTTR ,ITRA
     &          ,IZON   ,JJ     ,K1     ,K2     ,L      ,LD     ,LTYPE1
     &          ,M      ,NCNF   ,NNUD   ,NPTOT


      REAL*8::AF     ,AREAL  ,AREAN  ,H1     ,H2


      INTEGER*4  INDEX(12)   ,IPOS(NPAR)

      REAL*8::CFPARAM(12)    ,DERIV(NPAR)    ,DERIVTRA(6)  ,XPARAM(8)

C------------------------- FIRST EXECUTABLE STATEMENT
C------------------------- Cross all elements


       DO L=1,NUMEL

           NNUD = LNNDEL(L)
           LTYPE1 = LTYPE(L)
           LD = LDIM(L)
           IZON = LXPAREL(L,1)
           ISOT1 = ISOZ(IZON)
           ISOTTR = LD*(LD+1)/2
           AREAL = AREA(L)
           AREAN = AREAL/NNUD
           
      
C------------------------- Go through all transmissivity tensor components

           DO ITRA=1,ISOT1

               IP = IVPAR(INORPAR(ITRA)+IZON)


C------------------------- Derivatives have to be computed if IP>0 or 
C------------------------- even when the zonal parameter is not estimated, it
C------------------------- may depend on some generic parameters that may have 
C------------------------- to be estimated

               IF (IP.GT.0 .OR. IOFLLI.NE.0) THEN

                   M = 0

C------------------------- Computes derivatives 

                   JJ = INORPAR(ITRA)+IZON            ! Auxiliar variables

                   IF (IOFLLI.NE.0) THEN

                      NCNF = NFNLPAR(JJ) 
                      XPARAM(2) = HBASE(L)

                   ELSE

                       NCNF = 0

                   END IF !IOFLLI.NE.0

                   INDEX(1)=JJ

                   CFPARAM(1) = CFPAREL(L,1)

                   CALL DER_PARAM
     ; (CFPARAM  ,DTIM     ,EPSFLU   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     ; ,INTI     ,IOCAP    ,IOCTRA   ,IOFMLF   ,IP       ,L        
     ; ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     ; ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,1        ,NUMEL    
     ; ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT      
     ; ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG  
     ; ,PARACD   ,PARC(INORPAR(19)+1),HCALIT   ,HCALAN   ,XPARAM
     ; ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)


C------------------------- When no estimated parameters

                   IF (NPTOT.NE.0) THEN

C------------------------- Computes the value of Aflu_ij at each node of the 
C------------------------- current element

                       DO K1=1,NNUD-1

                           I1 = KXX(K1,L)
                           H1 = HAUX1(I1)

                           DO K2=K1+1,NNUD

                               I2 = KXX(K2,L)
                               H2 = HAUX1(I2)
                               AF = BIBI(M+ITRA,L)

                               IF (ISOT1.EQ.1 .AND. LTYPE1.NE.1)THEN

                                   AF = AF + BIBI(M+2,L)
                                   IF (LD.EQ.3) AF = AF + BIBI(M+3,L)

                               ELSEIF (ISOT1.EQ.2 .AND. LTYPE1.GT.1.AND.
     &                                 ITRA.EQ.1 .AND. LD.GT.2) THEN

                                   AF = AF + BIBI(M+3,L)

                               END IF

C------------------------- Computes the constant part of each derivative

                               AF = AF*DENSITY(L)*(H2-H1)
C------------------------- Updates the inverse problem RHS. There are two 
C------------------------- updates because we compute only half part of matrix
C------------------------- AFLU due to its symmetry

                               DO IC=1,NPTOT

                                   DERH(I1,IPOS(IC),INEW) =
     &                               DERH(I1,IPOS(IC),INEW)-AF*DERIV(IC)

                                   DERH(I2,IPOS(IC),INEW) =
     &                               DERH(I2,IPOS(IC),INEW)+AF*DERIV(IC)

                               END DO

                               M = M + ISOTTR

                           END DO   !K2
                       END DO     !K1


C------------------------- Buoyancy term contribution 
C------------------------- when Density dependent flow

                      IF (IODENS.EQ.1) THEN

C------------------------- Builds the derivative of transmissivity tensor
C------------------------- w. r. t. to the current direction of anisotropy

                          CALL DERTRA_TENSOR(DERIVTRA,ISOT1,ITRA,LD)

C------------------------- Contribution of buoyancy term.

                          CALL DERTRA_BUOYANCY
     &                        (AREA     ,BUOYANCY ,COORD    ,DENSITY
     &                        ,DERH     ,DERIV    ,DERIVTRA ,GP_COORD
     &                        ,GRADLOC  ,IDIMDERH ,INEW     ,IPOS
     &                        ,IODIM    ,ISOZ     ,KXX      ,L
     &                        ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE
     &                        ,LXPAREL  ,MAXPG    ,NPAREL   ,NPARF
     &                        ,NPTOT    ,NTYPEL  ,NUMEL     ,NUMNP
     &                        ,NZTRA    ,POINTWEIGHT)

                      END IF !IODENS.EQ.1


                   END IF       !NPTOT .NE. 0
               END IF        !IP .GT. 0  .OR.  IOFLLI .NE. 0
           END DO           !ITRA (Anisotropy)
       END DO              !Elements


       END SUBROUTINE DERTRA
