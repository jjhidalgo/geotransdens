      SUBROUTINE DERCOE_NOD_FLU
     &          (BETAC    ,CCALAN   ,CCALIT   ,CFPARNP  ,CREF
     &          ,DENSREF  ,DERH     ,DTIM     ,EPSTRA   ,FNT
     &          ,HAUX1    ,IBCOD    ,IDIMDERH ,IDIMFNT  ,IFLAGS
     &          ,INDSSTR  ,INEW     ,INORPAR  ,INTI     ,IOCAP
     &          ,IOFMLT   ,IOTRLI   ,IVPAR    ,IXPARNP  ,KXX
     &          ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG
     &          ,NFNLTIP  ,NFTPAR   ,NINT     ,NODE     ,NPAR
     &          ,NPARNP   ,NPPNP    ,NPZON    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD   ,PARC
     &          ,PARNP    ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,IPOS
     &          ,DERIV)

*******************************************************************************
*
* PURPOSE
*
*   Computes RHS of derivatives with respect to external concentration
*   associated to nodal fluxes at a given node.
*   Only when density dependent flow is solved.
*
* DESCRIPTION
*
*   Computes the RHS of derivatives with respect to external concentration
*   associated to nodal fluxes at a given node.
*   Only when density dependent flow is solved.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  CAUDAL                 Input/output flow at every node.                      
*  DERC                   Nodal concentration derivatives with respect to       
*                         estimated parameters.                                 
*  IBTCO                  Transport boundary condition index                    
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
*     AMS      3-2002     First coding (starting from TRANSIN-II)
*
*******************************************************************************

      IMPLICIT NONE
     
C------------------------- External

      INTEGER*4::IDIMDERH ,IDIMFNT  ,INDSSTR  ,INEW     ,INTI
     &          ,IOCAP    ,IOFMLT   ,IOTRLI   ,LMXNDL   ,NFLAGS
     &          ,NFNL     ,NINT     ,NODE     ,NPAR     ,NPARNP
     &          ,NPPNP    ,NPZON    ,NTYPAR   ,NUMEL    ,NUMNP
     &          ,NZPAR    ,IDIMWGT  ,IPNT_PAR


      REAL*8::BETAC    ,CREF      ,DENSREF  ,DTIM     ,EPSTRA

      REAL*8::DENS

      INTEGER*4::IBCOD(NUMNP)   ,IFLAGS(NFLAGS)       ,INORPAR(NTYPAR)
     &          ,IVPAR(NZPAR)   ,IXPARNP(NUMNP,NPARNP),KXX(LMXNDL,NUMEL)
     &          ,NFNLPAR(NZPAR) ,NFNLPRG(8,NFNL)      ,NFNLTIP(NFNL)
     &          ,NFTPAR(NZPAR)  ,NZONE_PAR(NTYPAR)

      REAL*8::CCALAN(NUMNP)         ,CCALIT(NUMNP)
     &       ,CFPARNP(NUMNP,NPARNP) ,DERH(NUMNP,NPAR,IDIMDERH)
     &       ,FNT(IDIMFNT,NINT)     ,HAUX1(NUMNP)
     &       ,PARC(NZPAR)           ,PARNP(NUMNP,NPPNP)
     &       ,PARACD(3,NFNL)        ,WGT_PAR(IDIMWGT)

       
C------------------------- Internal

      INTEGER*4::IB   ,IC   ,IP   ,JJ   ,NCNF ,NNUD ,NPTOT,NZC

      REAL*8::ALFAX    ,CAUD     ,CEXT    ,DENSEXT  ,FACTOR   ,HEAD

      INTEGER*4::INDEX(12)  ,IPOS(NPAR)

      REAL*8::CFPARAM(12)   ,DERIV(NPAR)   ,XPARAM(8)

C------------------------- First executable statement

      IB = IBCOD(NODE)
	CAUD = 0D0

	IF (IB.EQ.2) THEN

	    CAUD = PARNP(NODE,2)

	END IF !IB.EQ.2

	IF (IB.EQ.3 .OR. IB.EQ.4) THEN

          ALFAX = PARNP(NODE,3)
          HEAD = PARNP(NODE,1)
          CAUD = ALFAX*(HEAD - HAUX1(NODE))

	END IF !IB.EQ.3 OR. IB.EQ.4

	IF (IB.EQ.4) THEN
          
          IF (PARNP(NODE,2).GT.0) THEN

	        IF (CAUD.GT.0) THEN

	            CAUD = CAUD + PARNP(NODE,2)

	        ELSE

	            CAUD = PARNP(NODE,2)

	        END IF !CAUD.GT.0

	    END IF !PARNP(NODE,2).GT.0

	END IF !IB.EQ.4


C------------------------- If the node has boundary condition...
C------------------------- (CAUD.GT.0, otherwise there is no 'external'
C------------------------- concentration).
      IF (IB.GT.1 .AND. CAUD.GT.0) THEN 

          NZC = IXPARNP(NODE,7+INDSSTR)
          JJ = INORPAR(18) + NZC
          IP = IVPAR(JJ)
      
          IF (IP.NE.0 .OR. IOTRLI.NE.0) THEN

              IF (IOTRLI.NE.0) THEN
                  NCNF = NFNLPAR(JJ)
              ELSE
                    NCNF = 0
              END IF !IOTRLI.NE.0

              CFPARAM(1) = CFPARNP(NODE,7+INDSSTR)
              INDEX(1) = JJ
          
              CALL DER_PARAM
     &    (CFPARAM  ,DTIM     ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,0        ,IOFMLT   ,IP       ,NODE
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,NPZON    ,NUMEL
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     &    ,PARACD   ,PARC(INORPAR(18)+1),CCALIT   ,CCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

              IF (NPTOT.GT.0) THEN

                  CEXT = PARNP(NODE,4)
                  DENSEXT = DENS(DENSREF,BETAC,CEXT,CREF)

                  FACTOR =  BETAC*DENSEXT*CAUD

                  DO IC=1,NPTOT

                      DERH(NODE,IPOS(IC),INEW) = 
     &                       DERH(NODE,IPOS(IC),INEW) + FACTOR*DERIV(IC)

                  END DO !IC=1,NPTOT

              END IF !NPTOT.GT.0

          END IF !IP.NE.0

      END IF !IBT.GT.1 .AND. CAUDAL.GT.0

      END SUBROUTINE DERCOE_NOD_FLU
