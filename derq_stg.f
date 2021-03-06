      SUBROUTINE DERQ_STG
     &          (AREA     ,BETAC    ,CAUDAL   ,CAUX1    ,CAUX2
     &          ,CFPAREL  ,CREF     ,DENSITY  ,DENSREF  ,DERC
     &          ,DTIM     ,EPSFLU   ,FNT      ,HAUX1    ,HAUX2
     &          ,HCALAN   ,HCALIT   ,HINI     ,IBCOD    ,IBTCO
     &          ,IDIMDERC ,IDIMFNT  ,IFLAGS   ,INDSSTR  ,INEW
     &          ,INORPAR  ,INTI     ,IOCAP    ,IOCSTG   ,IODENS
     &          ,IOFLLI   ,IOFMLF   ,IOVRWC   ,IVPAR    ,KXX
     &          ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXSTG    ,NFLAGS
     &          ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR
     &          ,NINT     ,NPAR     ,NPAREL   ,NPPNP    ,NPZON
     &          ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR
     &          ,PARACD   ,PARC     ,PARNP    ,IDIMWGT  ,WGT_PAR
     &          ,IPNT_PAR ,IPOS     ,DERIV)

********************************************************************************
*
* PURPOSE
*
*   Computes the derivative of nodal flow w.r.t. storage (explicit dependence)
*
* DESCRIPTION
*
*   Computes the derivative of nodal flow w.r.t. storage (explicit dependence)
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
*  IODERPCALC             Used to avoid calling DER_PARAM more than once per element.
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
     &          ,IOCAP    ,IOCSTG   ,IODENS   ,IOFLLI   ,IOFMLF
     &          ,IOVRWC   ,LMXNDL   ,NFLAGS   ,NFNL     ,NINT
     &          ,NPAR     ,NPAREL   ,NPPNP    ,NPZON    ,NTYPAR
     &          ,NUMEL    ,NUMNP    ,NZPAR    ,IDIMWGT

      REAL*8::BETAC    ,CREF     ,DENSREF  ,DTIM     ,EPSFLU

      REAL*8::DENS   

      INTEGER*4::IBCOD(NUMNP)   ,IBTCO(NUMNP)     ,IFLAGS(NFLAGS)
     &          ,INORPAR(NTYPAR),IVPAR(NZPAR)     ,KXX(LMXNDL,NUMEL)
     &          ,LNNDEL(NUMEL)  ,LTYPE(NUMEL)     ,LXSTG(NUMEL)
     &          ,NFNLPAR(NZPAR) ,NFNLTIP(NFNL)    ,NFNLPRG(8,NFNL)
     &          ,NFTPAR(NZPAR)  ,NZONE_PAR(NTYPAR)  ,IPNT_PAR

      REAL*8::AREA(NUMEL)               ,CAUDAL(NUMNP)
     &       ,CAUX1(NUMNP)              ,CAUX2(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL)     ,DENSITY(NUMEL)
     &       ,DERC(NUMNP,NPAR,IDIMDERC) ,FNT(IDIMFNT,NINT)
     &       ,HAUX1(NUMNP)              ,HAUX2(NUMNP)
     &       ,HCALAN(NUMNP)             ,HCALIT(NUMNP)
     &       ,HINI(NUMNP)               ,PARACD(3,NFNL)
     &       ,PARC(NZPAR)   ,PARNP(NUMNP,NPPNP)    ,WGT_PAR(IDIMWGT)



C-------------------- Internal

      INTEGER*4::I      ,IB     ,IBT    ,INODE  ,IP     ,IPAR   ,JJ
     &          ,K      ,KNODE  ,L      ,LTY    ,NCNF   ,NNUD   ,NPTOT
     &          ,NZ

      REAL*8::AREALN   ,CAUD     ,CEXT_CI  ,CEXTNODE ,CNODE    ,DELTHAVG
     &       ,DENSEXT  ,DENSNODE ,DENSSTG  ,DERWTV   ,STG      ,STGCFLU

      INTEGER*4::INDEX(12)   ,IPOS(NPAR)
      REAL*8::CFPARAM(12)    ,DERIV(NPAR)  ,XPARAM(8)



C-------------------- First executable statement

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          AREALN = AREA(L)/NNUD
          LTY = LTYPE(L)
          NZ = LXSTG(L)
          JJ = INORPAR(7)+NZ
          IP = IVPAR(JJ)

          IF (IP.NE.0 .OR. IOFLLI.NE.0) THEN

C------------------------- Derivatives of storage

              IF (IOFLLI.NE.0) THEN

                  NCNF = NFNLPAR(JJ)

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


              IF (NPTOT.GT.0) THEN

                  DO I=1,NNUD

                      INODE = KXX(I,L)
	                IB = IBCOD(INODE)
	                IBT = IBTCO(INODE)
	                CAUD = CAUDAL(INODE)

C------------------------- Only nodes with prescribed head have this
C------------------------- contribution (Q = Q(S,K,etc.), throug DFLU).

                      IF (IB.EQ.1 .AND. CAUD.GT.0
     &                   .AND.(IBT.EQ.2 .OR. IBT.EQ.3))  THEN

C------------------------- Initilization of nodewise auxiliar variables

                          CEXTNODE = PARNP(INODE,4)
                          CNODE = CAUX1(INODE)


                          IF (IODENS.EQ.1) THEN

                              DENSEXT=DENS(DENSREF,BETAC,CEXTNODE,CREF)
                              DENSNODE=DENS(DENSREF,BETAC,CNODE,CREF)

                              CEXT_CI = DENSEXT*(CEXTNODE - CNODE)

	                        IF (IOVRWC.LT.2) THEN

	                            DENSSTG = DENSITY(L)

	                        ELSE

	                            DENSSTG = DENSNODE

	                        END IF !IOVRWC.LT.2

                          ELSE

                              CEXT_CI = CEXTNODE - CNODE

	                        DENSSTG = 1D0
	                        DENSEXT = 1D0
	                        DENSNODE = 1D0

                          END IF !IODENS.EQ.1

C------------------------- DFLU contribution.

                          STG = DENSSTG*AREALN*CEXT_CI*HAUX2(INODE)

C------------------------- The derivative w. r. t. parameters are
C------------------------- computed only once per element.

C------------------------- Contribution to RHS.


                          DO IPAR=1,NPTOT

C-------------------- Derivatives w. r. t. storage parameters.

                              DERC(INODE,IPOS(IPAR),INEW) =
     &                                       DERC(INODE,IPOS(IPAR),INEW)
     &                                       + DERIV(IPAR)*STG


C------------------------- Derivative of WATVOL in CFLU w. r. t. storage
C------------------------- Only if density dependent flow.

                              IF (IODENS.EQ.1) THEN

	                            SELECT CASE(IOVRWC)

	                            CASE(0)

	                                DERWTV = 0D0

	                            CASE(1)

                                      DO K=1,NNUD

	                                    KNODE = KXX(K,L)
	                                    DELTHAVG = DELTHAVG
     &                                        + HAUX1(KNODE)-HINI(KNODE)

                                      END DO !K=1,NNUD

                                      DELTHAVG = DELTHAVG/NNUD
	                                DERWTV = DERIV(IPAR)*DELTHAVG

	                            CASE(2)

	                                KNODE = KXX(K,L)
	                                DERWTV = DERIV(IPAR)*
     &                                      (HAUX1(INODE) - HINI(INODE))

	                            END SELECT !IOVRWC

C------------------------- CFLU contribution.
C------------------------- Q_CFLU = CFLU*deltaC /DENSEXT 

                              STGCFLU = DENSSTG*AREALN*CEXT_CI*BETAC
     &                                 *DERWTV*CAUX2(INODE)/DENSEXT

C-------------------- Derivatives w. r. t. storage paremters (CFLU)

                              DERC(INODE,IPOS(IPAR),INEW) =
     &                                     DERC(INODE,IPOS(IPAR),INEW)
     &                                     + STGCFLU
                              END IF !IODENS.EQ.1


                              IF (IFLAGS(30).EQ.1) THEN
                                  WRITE(77,10) I,L,IP,DERIV(IPAR),STG
   10                             FORMAT (3I6,2E20.10)
                              END IF !IFLAGS(30).EQ.1


                          END DO !IPAR=1,NPTOT

                      END IF !IB.EQ.1 .AND. CAUD.GT.0 .AND. ...

                  END DO !I=1,NNUD

              END IF !NPTOT.GT.0

          END IF !IP.NE.0

      END DO !L=1,NUMEL

      END SUBROUTINE DERQ_STG
