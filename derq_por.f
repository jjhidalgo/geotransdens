      SUBROUTINE DERQ_POR
     &          (ACTH     ,AREA     ,BETAC    ,CAUDAL   ,CAUX1
     &          ,CAUX2    ,CCALAN   ,CCALIT   ,CFPAREL  ,CREF
     &          ,DENSITY  ,DENSREF  ,DERC     ,DTIM     ,EPSTRA
     &          ,FNT      ,IBCOD    ,IBTCO    ,IDIMDERC ,IDIMFNT
     &          ,IFLAGS   ,INDSSTR  ,INEW     ,INORPAR  ,INTI
     &          ,IOCAP    ,IOFMLT   ,IOTRLI   ,IOVRWC   ,IVPAR
     &          ,KXX      ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXPOR
     &          ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &          ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   ,NPPNP
     7          ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR
     &          ,PARACD   ,PARC     ,PARNP    ,IDIMWGT  ,WGT_PAR
     &          ,IPNT_PAR ,IPOS     ,DERIV)

********************************************************************************
*
* PURPOSE
*
*   Computes the derivative of nodal flow w.r.t. porosity in prescribed head
*   nodes with mass flow transport boundary condition.
*
* DESCRIPTION
*
*   Computes the derivative of nodal flow w.r.t. porosity in prescribed head
*   nodes with mass flow transport boundary condition.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  CAUDAL                 Input/output flow at every node.                      
*  CAUX1                  Array containing concentrations, ponderated by THETAT
*                         time factor
*  CNST                   CNST(i,j,k) is the integral of the product of         
*                         interpolation functions i and j in an element of      
*                         type k divided by the AREA of the element. It is used 
*                         only in consistent scheme.                            
*  DERC                   Nodal concentration derivatives with respect to       
*                         estimated parameters.                                 
*  HAUX2                  Array containing difference of heads in two           
*                         consecutives times. Is equal to HCAL-HCALAN/TIME STEP 
*  IBCOD                  Flow boundary condition index                         
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing the type of each element            
*  LXSTG                  Storage coefficient zone number at a given element    
*  PARC                   Vector containing calculated values for all
*                         parameters                                            
*  PAREL                  Parameter values at every element and current time    
*                         for all nodal parameters (each value is computed as   
*                         the product of up to four terms:                      
*                           elem. coeff*zonal value*time funct.*nonl. funct. )  
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
*  IOCNSF                 Scheme for storage term in flow problem               
*  LMXNDL                 Maximum number of nodes per element                   
*  NPAR                   Total number of parameters to be estimated            
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
*  NPPNP                  Total number of parameters by nodes (not confuse      
*                         with NPARNP, because in this casethere is no          
*                         difference between a given parameter in steady or tr.)
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
*  DER_STOR               Derivative of storage w.r.t. storage (zone)
*  LTY                    Type of element
*  NNUD                   Number of nodes of the current element                
*  NZ                     Storage zone number
*  SX                     Derivative of flow at node I w.r.t. storage (explicit 
*                         dependence)
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
* HISTORY
*
*     AMS      3-2002     First coding (starting from TRANSIN-II)
*
********************************************************************************

      IMPLICIT NONE

C-------------------- External

      INTEGER*4::IDIMDERC ,IDIMFNT  ,INDSSTR  ,INEW     ,INTI
     &          ,IOCAP    ,IOFMLT   ,IOTRLI   ,IOVRWC   ,LMXNDL
     &          ,NFLAGS   ,NFNL     ,NINT     ,NPAR     ,NPAREL
     &          ,NTYPAR   ,NUMEL    ,NUMNP    ,NPPNP    ,NZPAR
     &          ,IDIMWGT  ,IPNT_PAR

    
      INTEGER*4::IBCOD(NUMNP)    ,IBTCO(NUMNP)       ,IFLAGS(NFLAGS)
     &          ,INORPAR(NTYPAR) ,IVPAR(NZPAR)       ,KXX(LMXNDL,NUMEL)
     &          ,LNNDEL(NUMEL)   ,LTYPE(NUMEL)       ,LXPOR(NUMEL)
     &          ,NFNLPAR(NZPAR)  ,NFNLPRG(8,NFNL)    ,NFNLTIP(NFNL)
     &          ,NFTPAR(NZPAR)   ,NZONE_PAR(NTYPAR)

      REAL*8::ACTH(NUMEL)        ,AREA(NUMEL)
     &       ,CAUDAL(NUMNP)      ,CAUX1(NUMNP)
     &       ,CAUX2(NUMNP)       ,CCALAN(NUMNP)
     &       ,CCALIT(NUMNP)      ,CFPAREL(NUMEL,NPAREL)
     &       ,DENSITY(NUMEL)     ,DERC(NUMNP,NPAR,IDIMDERC)
     &       ,FNT(IDIMFNT,NINT)  ,PARACD(3,NFNL)
     &       ,PARC(NZPAR)   ,PARNP(NUMNP,NPPNP)    ,WGT_PAR(IDIMWGT)
     &        

      REAL*8::BETAC,CREF,DENSREF,DTIM,EPSTRA

      REAL*8::DENS

C-------------------- Internal

      INTEGER*4::I      ,IB     ,IBT    ,INODE  ,IP     ,IPAR   ,JJ
     &          ,L      ,LTY    ,NCNF   ,NNUD   ,NPTOT  ,NZ

      REAL*8::ACTHL    ,BETA_AREALN        ,CEXT_CI  ,CEXTNODE ,CNODE
     &       ,DENSNODE ,DENSPOR  ,POR

      INTEGER*4::INDEX(12)   ,IPOS(NPAR)
      REAL*8::CFPARAM(12)    ,DERIV(NPAR)  ,XPARAM(8)


C-------------------- First executable statement

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          NZ = LXPOR(L)
	    JJ = INORPAR(15)+NZ
          IP = IVPAR(JJ)
          BETA_AREALN =BETAC*AREA(L)/NNUD
	    ACTHL = ACTH(L)
          LTY = LTYPE(L)

	    IF (IP.NE.0 .OR. IOTRLI.NE.0) THEN

C------------------------- Derivatives of storage

              IF (IOTRLI.NE.0) THEN

          	    NCNF = NFNLPAR(JJ)

	        ELSE

          	    NCNF = 0

	        END IF !IOTRLI.NE.0

              INDEX(1) = JJ
              CFPARAM(1) = CFPAREL(L,7)

C------------------------- Derivatives of porosity if density is not constant

              CALL DER_PARAM
     &    (CFPARAM  ,DTIM     ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,0        ,IOFMLT   ,IP       ,L
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,1        ,NUMEL
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     &    ,PARACD   ,PARC(INORPAR(19)+1),CCALIT   ,CCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

              IF (NPTOT.GT.0) THEN
          
                  DO I=1,NNUD

                      INODE = KXX(I,L)
	                IB = IBCOD(INODE)
	                IBT = IBTCO(INODE)

	                IF (IB.EQ.1 .AND. CAUDAL(INODE).GT.0
     &                   .AND. (IBT.EQ.2 .OR. IBT.EQ.3)) THEN

                          CEXTNODE = PARNP(INODE,4)
                          CNODE = CAUX1(INODE)

                          DENSNODE = DENS(DENSREF,BETAC,CNODE,CREF)

                          CEXT_CI = CEXTNODE - CNODE

	                    IF (IOVRWC.EQ.2) THEN

	                        DENSPOR = DENSNODE

	                    ELSE

	                        DENSPOR = DENSITY(L)

	                    END IF !IOVRWC.LT.2

C------------------------- POR = @(DENSEXT�Q)/@por, then, it is not
C------------------------- necessary to multiply CEXT_CI by DENSEXT.

                          POR = DENSPOR*BETA_AREALN*ACTHL*CEXT_CI
     &                          *CAUX2(INODE)

C------------------------- CFLU contribution and correction of flow
C------------------------- if density is variable

C------------------------- Contribution to RHS.

                          DO IPAR=1,NPTOT

C-------------------- Derivatives w. r. t. porosity parameters
C-------------------- Only when density dependent flow

                              DERC(INODE,IPOS(IPAR),INEW) =
     &                                       DERC(INODE,IPOS(IPAR),INEW)
     &                                       + DERIV(IPAR)*POR

                          END DO !IPAR=1,NPTOT

                      END IF !IBCOD(I).EQ.1 .AND. CAUDAL(I).GT.0

	            END DO ! I=1,NNUD

              END IF ! NPTOT.EQ.0

          END IF !IP.NE.0 .OR. IOTRLI.NE.0

      END DO !L=1,NUMEL

      END SUBROUTINE DERQ_POR
