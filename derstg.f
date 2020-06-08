      SUBROUTINE DERSTG 
     &          (AREA     ,BETAC    ,CAUX1    ,CAUX2    ,CFPAREL
     &          ,CREF     ,DENSITY  ,DENSREF  ,DERH     ,DTIM
     &          ,EPSFLU   ,FNT      ,HAUX2    ,HBASE    ,HCALAN
     &          ,HCALIT   ,IDIMDERH ,IDIMFNT  ,IFLAGS   ,INDSSTR
     &          ,INEW     ,INORPAR  ,INTI     ,IOCSTG   ,IODENS
     &          ,IOFLLI   ,IOFMLF   ,IOVRWC   ,IVPAR    ,KXX
     &          ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXPAREL  ,NFLAGS
     &          ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR
     &          ,NINT     ,NPAREL   ,NPARF    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD   ,PARC
     &          ,THETAT   ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,NPAR
     &          ,IPOS     ,DERIV)

*****************************************************************************
*
* PURPOSE
*
*     Compute the contribution of the derivatives of DFLU matrix with respect
*     to storage coeff. to the RHS of the derivatives of h  wrt parameters
*
* DESCRIPTION
*
*     Compute the contribution of the derivatives of AFLU matrix with respect
*     to storage coeff. to the RHS of the derivatives of h  wrt parameters
*     More precisely, adds to the RHS the expression
*
*                     1           d DFLU          k+1     k
*                  ------  * ----------------  ( h    -  h  )
*                  Delta t   d Storage coeff.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  CFPAREL                Array containing node coefinient of element j         
*                         corresponding to INpar index parameter zone.          
*  CNST                   Interpolation functions gradient for a given element  
*                         nodes                                                 
*  DERH                   Nodal head derivatives with respect to estimated      
*                         flow parameters.                                      
*  HAUX2                  Array containing diference of HEADS in two            
*                         consecutives times.                                   
*  HBASE                  Bottom level of the aquifer
*  HCAL                   Computed heads at every node                          
*  HCALAN                 Head level at previous time                           
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
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
*  FNT                    Array containing time functions values                
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
*  IOCNSF                 Scheme for storage term in flow problem               
*  IOCSTG                 Nodal ponderation option to compute the contribution  
*                         of the nonlinear function to the storage value        
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
*
* INTERNAL VARIABLES: SCALARS
*
*  IZON                   Storage coeff. zone number of the current element
*  JJ                     Location of storage coeff. values of the current
*                         storage coeff. zone in arrays PARC, NFNLPAR, etc
*  L                      Current element
*  LTYPE1                 Type of the current element
*  NCNF                   Nonlinear function number of the current
*                         storage coeff. zone
*  NNUD                   Number of nodes of the current element
*  NPTOT                  Total number of parameters to be estimated in the
*                         current storage coeff. zone (includes zonal plus
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

      INTEGER*4::IDIMDERH ,IDIMFNT  ,INDSSTR  ,INEW     ,INTI
     &          ,IOCAP    ,IOCSTG   ,IODENS   ,IOFLLI   ,IOFMLF
     &          ,IOVRWC   ,LMXNDL   ,NFLAGS   ,NFNL     ,NINT
     &          ,NPAREL   ,NPARF    ,NTYPAR   ,NUMEL    ,NUMNP
     &          ,NZPAR    ,IDIMWGT  ,IPNT_PAR ,NPAR

      REAL*8::BETAC    ,CREF     ,DENSREF  ,DTIM     ,EPSFLU   ,THETAT

      REAL*8::DENS
          
      REAL*8::AREA(NUMEL)               ,CAUX1(NUMNP)  ,CAUX2(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL)     ,DENSITY(NUMEL)
     &       ,DERH(NUMNP,NPARF,IDIMDERH),FNT(IDIMFNT,NINT),HAUX2(NUMNP)
     &       ,HBASE(NUMEL)              ,HCALAN(NUMNP)    ,HCALIT(NUMNP)
     &       ,PARACD(3,NFNL)            ,PARC(NZPAR)  ,WGT_PAR(IDIMWGT)

      INTEGER*4::IFLAGS(NFLAGS)   ,INORPAR(NTYPAR),IVPAR(NZPAR)
     &          ,KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)  ,LXPAREL(NUMEL,NPAREL)
     &          ,LTYPE(NUMEL)     ,NFNLPAR(NZPAR) ,NFNLTIP(NFNL)
     &          ,NFNLPRG(8,NFNL)  ,NFTPAR(NZPAR)  ,NZONE_PAR(NTYPAR)


C------------------------- Internal

      INTEGER*4::I        ,IC       ,INODE    ,IP       ,IZON
     &          ,JJ       ,L        ,LTYPE1   ,NCNF     ,NNUD
     &          ,NPTOT    ,NPZON

      REAL*8::AREALN   ,BETAREALN,CNODE    ,DH_AVG

      INTEGER*4::IPOS(12)  ,INDEX(NPAR)

      REAL*8::CFPARAM(12)  ,DERIV(NPAR)  ,DH(LMXNDL)
     &       ,RHO(LMXNDL)  ,XPARAM(8)

      


C------------------------- FIRST EXECUTABLE STATEMENT.


C------------------------- Cross over all elements

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          LTYPE1 = LTYPE(L)

          IZON = LXPAREL(L,2)
          IP = IVPAR(INORPAR(7)+IZON)
	    AREALN = AREA(L)/NNUD

C------------------------- Computes integration factors to save operations

          IF (IODENS.EQ.1) THEN

              BETAREALN = BETAC*AREALN

C------------------------- Computes change of head  multiplied by THETAT
C------------------------- WTV(k+th) = WTV(k) + STG*THETAT*Delta_H
              IF (IOVRWC.EQ.0) THEN

                  RHO(:) = DENSITY(L)
	            DH(:) = 0D0

              ELSE IF (IOVRWC.EQ.1) THEN

	            DO I=1,NNUD

	                INODE = KXX(I,L)
	                DH_AVG = DH_AVG + HCALIT(INODE) - HCALAN(INODE)
	                RHO(I) = DENSITY(L)

	            END DO !I=1,NNUD

                  DH_AVG = DH_AVG /NNUD

	            DH(:) = DH_AVG*THETAT

	        ELSE IF (IOVRWC.EQ.2) THEN

                  DO I=1,NNUD

	                INODE = KXX(I,L)
	                CNODE = CAUX1(INODE)
	                DH(I) = (HCALIT(INODE) - HCALAN(INODE))*THETAT
	                RHO(I) = DENS(DENSREF,BETAC,CNODE,CREF)

	            END DO !I=1,NNUD

	        END IF !IOVRWC.EQ.1,2

          ELSE

              BETAREALN = AREALN
	        RHO(:) = 1D0

          END IF !IODENS.EQ.1

C------------------------- Derivatives have to be computed if IP>0 or
C------------------------- even when the zonal parameter is not estimated, it
C------------------------- may depend on some generic parameters that may have
C------------------------- to be estimated

          IF (IP.GT.0 .OR. IOFLLI.NE.0) THEN

              JJ = INORPAR(7)+IZON            ! Auxiliar variables

              NPZON = 1

              IF (IOFLLI.NE.0) THEN

                  NCNF = NFNLPAR(JJ)

                  IF (NCNF.NE.0) THEN

                      NPZON = 2
                      INDEX(2) = INORPAR(15) + LXPAREL(L,7)
                      CFPARAM(2) = CFPAREL(L,7)
                      XPARAM(1) = PARC(JJ)*CFPAREL(L,2)
                      XPARAM(2) = HBASE(L)
                      XPARAM(3) = PARC(INDEX(2))

                  END IF !NCNF.NE.0

              ELSE

                  NCNF = 0

              END IF !IOFLLI.NE.0
               
              INDEX(1) = JJ

              CFPARAM(1) = CFPAREL(L,2)

              CALL DER_PARAM
     & (CFPARAM  ,DTIM     ,EPSFLU   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     & ,INTI     ,IOCAP    ,IOCSTG   ,IOFMLF   ,IP       ,L        
     & ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     & ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,NPZON    ,NUMEL    
     & ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT      
     & ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG  
     & ,PARACD   ,PARC(INORPAR(19)+1),HCALIT   ,HCALAN   ,XPARAM
     ; ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

C------------------------- Updates the inverse problem RHS.                 
C------------------------- Only Lumped scheme

              IF (NPTOT.NE.0) THEN

                  DO I=1,NNUD

                      INODE = KXX(I,L)

                      DO IC=1,NPTOT

                          DERH(INODE,IPOS(IC),INEW) = 
     &                        DERH(INODE,IPOS(IC),INEW)
     &                       -RHO(I)*AREALN*HAUX2(INODE)*DERIV(IC)

                          IF (IODENS.EQ.1 .AND. IOVRWC.GT.0) THEN

C------------------------- Contribution of CFLU term when water content is not
C------------------------- constant

                              DERH(INODE,IPOS(IC),INEW) = 
     &                        DERH(INODE,IPOS(IC),INEW)
     &                        -RHO(I)*BETAREALN*DERIV(IC)
     &                         *DH(I)*CAUX2(INODE)

                          END IF !IODENS.EQ.1

                      END DO !IC=1,NPTOT

                  END DO !I=1,NNUD

              END IF ! NPTOT .NE. 0

          END IF ! IP.GT.0 .OR. IOFLLI.NE.0

      END DO ! L=1,NUMEL 


      END SUBROUTINE DERSTG 
