      SUBROUTINE DERCOE_ARR_FLU
     &          (AREA     ,BETAC    ,CCALAN   ,CCALIT   ,CFPAREL
     &          ,CREF     ,DENSREF  ,DERH     ,DTIM     ,EPSTRA
     &          ,FNT      ,IBCOD    ,IDIMDERH ,IDIMFNT  ,IFLAGS
     &          ,INARR    ,INDSSTR  ,INEW     ,INORPAR  ,INTI
     &          ,IOCAP    ,IOFMLT   ,IOTRLI   ,IVPAR    ,KXX
     &          ,LMXNDL   ,LNNDEL   ,LXPAREL  ,NFLAGS   ,NFNL
     &          ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT
     &          ,NPAR     ,NPAREL   ,NPPEL    ,NPZON    ,NTYPAR
     &          ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD
     &          ,PARC     ,PAREL    ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR
     &          ,IPOS     ,DERIV)

*******************************************************************************
*
* PURPOSE
*
*   Computes flow RHS of derivatives with respect to the external concentration
*   assciated to areal recharge.
*
* DESCRIPTION
*
*   Computes flow RHS of derivatives with respect to the external concentration
*   assciated to areal recharge.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  CAUDAL                 Input/output flow at every node.                      
*  DERC                   Nodal concentration derivatives with respect to       
*                         estimated parameters.                                 
*  IBCOD                  Flow boundary condition index                    
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  IXPARNP                Array containing zone number of each node j,          
*                         corresponding to INpar index parameter zone.          
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*  LXPAREL                Array containing zone number of each element j,       
*                         corresponding to INpar index parameter zone.          
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
*  INARR                  Index for areal recharge                              
*                         in array variables (LXPAREL and CFPAREL)              
*  INCON                  Index for external concentration                      
*                         in array variables (IXPARNP and CFPARNP)              
*  INEW                   Location in array DERC
*  LMXNDL                 Maximum number of nodes per element                   
*  NPAR                   Total number of parameters to be estimated            
*  NPAREL                 Number of element parameters in current problem       
*  NPARNP                 Number of nodal parameters in current problem         
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
*  IBT                    Nodal transport boundary condition
*  L                      Current element
*  NNUD                   Number of nodes of the current element                
*  NZA                    Areal recharge zone at current element
*  NZC                    External concentration zone at current element
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
* HISTORY
*
*     JHG     12-2005     First coding
*
*******************************************************************************

      IMPLICIT NONE
     
C------------------------- External

      INTEGER*4::IDIMDERH ,IDIMFNT  ,INARR    ,INDSSTR  ,INEW
     &          ,INTI     ,IOCAP    ,IOFMLT   ,IOTRLI   ,LMXNDL
     &          ,NFLAGS   ,NFNL     ,NINT     ,NPAR     ,NPAREL
     &          ,NPPEL    ,NPZON    ,NTYPAR   ,NUMEL    ,NUMNP
     &          ,NZPAR    ,IDIMWGT  ,IPNT_PAR


      REAL*8::BETAC    ,CREF      ,DENSREF  ,DTIM     ,EPSTRA

      REAL*8::DENS

      INTEGER*4::IBCOD(NUMNP)         ,IFLAGS(NFLAGS)   ,INORPAR(NTYPAR)
     &          ,IVPAR(NZPAR)         ,KXX(LMXNDL,NUMEL),LNNDEL(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR)   ,NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL)        ,NFTPAR(NZPAR)  ,NZONE_PAR(NTYPAR)

      REAL*8::AREA(NUMEL)          ,CCALAN(NUMNP)    ,CCALIT(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL),DERH(NUMNP,NPAR,IDIMDERH)
     &       ,FNT(IDIMFNT,NINT)    ,PARC(NZPAR)      ,PAREL(NUMEL,NPPEL)
     &       ,PARACD(3,NFNL)       ,WGT_PAR(IDIMWGT)
       
C------------------------- Internal

      INTEGER*4::IB   ,IC   ,IP   ,JJ   ,K    ,KNODE    ,L    ,NCNF
     &          ,NNUD ,NPTOT,NZA  ,NZC

      REAL*8::AREALN   ,CEXT    ,DENSEXT  ,FACTOR   ,RECHRG   ,VAUX

      INTEGER*4::INDEX(12)  ,IPOS(NPAR)

      REAL*8::CFPARAM(12)   ,DERIV(NPAR)   ,XPARAM(8)

C------------------------- First executable statement


      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          NZA = LXPAREL(L,INARR)
          NZC = LXPAREL(L,10)
          RECHRG = PAREL(L,8)

          IF (NZC.NE.0 .AND. NZA.NE.0 .AND. RECHRG.GT.0) THEN

              JJ = INORPAR(18) + NZC
              IP = IVPAR(JJ)

              IF (IP.NE.0 .OR. IOTRLI.NE.0) THEN

                  IF (IOTRLI.NE.0) THEN
                      NCNF = NFNLPAR(JJ)
                  ELSE
                      NCNF = 0
                  END IF !IOTRLI.NE.0

                  CFPARAM(1) = CFPAREL(L,10)
                  INDEX(1) = JJ
              
                  CALL DER_PARAM
     &    (CFPARAM  ,DTIM     ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,0        ,IOFMLT   ,IP       ,L
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,NPZON    ,NUMEL
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     &    ,PARACD   ,PARC(INORPAR(18)+1),CCALIT   ,CCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

                  IF (NPTOT.NE.0) THEN

                      NNUD = LNNDEL(L)
                      AREALN = AREA(L)/NNUD

                      CEXT = PAREL(L,15)

                      DENSEXT = DENS(DENSREF,BETAC,CEXT,CREF)

                      FACTOR = BETAC*DENSEXT*RECHRG*AREALN


                      DO IC=1,NPTOT

                      	VAUX = FACTOR*DERIV(IC)

                          DO K=1,NNUD

                              KNODE = KXX(K,L)
	                        IB = IBCOD(KNODE)

C------------------------- This contribution must not be added to nodes with
C------------------------- prescribed head.

	                        IF (IB.NE.1) THEN

                                  DERH(KNODE,IPOS(IC),INEW) =
     &                                  DERH(KNODE,IPOS(IC),INEW) + VAUX

                              END IF !IB.NE.1

                          END DO !K=1,NNUD

                      END DO !IC=1,NPTOT

                  END IF !NPTOT.NE.0

              END IF !IP.NE.0 .OR. IOTRLI.NE.0

          END IF !NZC.NE.0 .AND. NZA.NE.0

      END DO !L=1,NUMEL

      END SUBROUTINE DERCOE_ARR_FLU
