      SUBROUTINE DERCOE
     &          (AREA     ,BETAC    ,CAUDAL   ,CAUX1    ,CCALAN
     &          ,CCALIT   ,CFPAREL  ,CFPARNP  ,CREF     ,DENSREF
     &          ,DERC     ,DTIM     ,EPSTRA   ,FNT      ,IBCOD
     &          ,IBTCO    ,IDIMDERC ,IDIMFNT  ,IFLAGS   ,INARR
     &          ,INCON    ,INDSSTR  ,INEW     ,INORPAR  ,INTI
     &          ,IOCAP    ,IODENS   ,IOFMLT   ,IOTRLI   ,IVPAR
     &          ,IXPARNP  ,KXX      ,LMXNDL   ,LNNDEL   ,LXPAREL
     &          ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &          ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   ,NPARNP
     &          ,NPPEL    ,NPPNP    ,NPZON    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD   ,PARC
     &          ,PAREL    ,PARNP    ,THETAT   ,TINC     ,TINTERVOBS
     ;          ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV
     &          ,ITPTVAR  ,WSPECHEAT)

*******************************************************************************
*
* PURPOSE
*
*   Computes RHS of derivatives with respect to external concentration
*
* DESCRIPTION
*
*   Computes the RHS of derivatives with respect to external concentration
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

      INTEGER*4::IDIMDERC ,IDIMFNT  ,INARR    ,INCON    ,INDSSTR
     &          ,INEW     ,INTI     ,IOCAP    ,IODENS   ,IOFMLT
     &          ,IOTRLI   ,LMXNDL   ,NFLAGS   ,NFNL     ,NINT
     &          ,NPAR     ,NPAREL   ,NPARNP   ,NPPEL    ,NPPNP
     &          ,NPZON    ,NTYPAR   ,NUMEL    ,NUMNP    ,NZPAR
     &          ,IDIMWGT  ,IPNT_PAR ,ITPTVAR


      REAL*8::BETAC    ,CREF      ,DENSREF  ,DTIM     ,EPSTRA
     &       ,THETAT  ,TINC       ,TINTERVOBS         ,WSPECHEAT

      REAL*8::DENS

      INTEGER*4::IBCOD(NUMNP)   ,IBTCO(NUMNP)         ,IFLAGS(NFLAGS)
     &          ,INORPAR(NTYPAR)
     &          ,IVPAR(NZPAR)   ,IXPARNP(NUMNP,NPARNP),KXX(LMXNDL,NUMEL)
     &          ,LNNDEL(NUMEL)  ,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR)
     &          ,NFNLPRG(8,NFNL),NFNLTIP(NFNL)        ,NFTPAR(NZPAR)
     &          ,NZONE_PAR(NTYPAR)

      REAL*8::AREA(NUMEL)            ,CAUDAL(NUMNP)
     &       ,CAUX1(NUMNP)           ,CCALAN(NUMNP)
     &       ,CCALIT(NUMNP)          ,CFPAREL(NUMEL,NPAREL)
     &       ,CFPARNP(NUMNP,NPARNP)  ,DERC(NUMNP,NPAR,IDIMDERC)
     &       ,FNT(IDIMFNT,NINT)      ,PARC(NZPAR)
     &       ,PAREL(NUMEL,NPPEL)     ,PARNP(NUMNP,NPPNP)
     &       ,PARACD(3,NFNL)         ,WGT_PAR(IDIMWGT)


C------------------------- Internal

      INTEGER*4::I    ,IBF  ,IBT  ,IC   ,IP   ,JJ   ,K    ,L    ,NCNF
     &          ,NNUD ,NPTOT,NZA  ,NZC

      REAL*8::AREALN   ,CAUD     ,CEXT     ,CEXT_CINT,CNOD     ,CREC
     &       ,DENSEXT  ,DENSREC  ,DTIMAUX  ,FACT     ,FACTOR   ,RECHRG
     &       ,VAUX     ,VAUX2

      INTEGER*4::INDEX(12)  ,IPOS(NPAR)

      REAL*8::CFPARAM(12)   ,DERIV(NPAR)   ,XPARAM(8)

C------------------------- First executable statement

C------------------------- Areal Recharge term

      DO L=1,NUMEL

          NZA = LXPAREL(L,INARR)
          NZC = LXPAREL(L,10)
          RECHRG = PAREL(L,8)
          NNUD = LNNDEL(L)

          IF (NZC.NE.0 .AND. NZA.NE.0 .AND. RECHRG.GT.0D0) THEN

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

                      IF (IODENS.EQ.1) THEN !It's already been checked
                                            !that recharge is positive.
                          CREC = PAREL(L,15)

                          DENSREC = DENS(DENSREF,BETAC,CREC,CREF)

                          FACT = DENSREC*RECHRG*AREALN

                      ELSE

                          FACT = RECHRG*AREALN

                      END IF !IODENS.EQ.1

                      DO IC=1,NPTOT


                          DO K=1,NNUD

                              I = KXX(K,L)
                              IBT = IBTCO(I)

C------------------------- This contribution must not be added to nodes
C------------------------- with prescribed concentration, since DERC has
C------------------------- been set to zero for those nodes before
C------------------------- entering this subroutine.

                              IF (IBT.NE.1) THEN

                                  CNOD = CAUX1(I)
                                  IBF = IBCOD(I)
                                  CAUD = CAUDAL(I)

C------------------------- Contribution of DENSREC*REC(CEXT-CNOD) term.

                                  VAUX = FACT*(1D0 + BETAC*(CREC-CNOD))
     &                                   *DERIV(IC)

                                  DERC(I,IPOS(IC),INEW) =
     &                                      DERC(I,IPOS(IC),INEW) + VAUX

C------------------------- Implicit contribution to nodes with 
C------------------------- prescribed head through recharge (from nodal
C------------------------- flux computation)...

                                  IF (IBF.EQ.1 .AND. CAUD.GT.0
     &                               .AND.IODENS.EQ.1) THEN

                                      CNOD = CAUX1(I)
                                      CEXT = PARNP(I,4)

                                      VAUX2 =BETAC*DENSREC*RECHRG*AREALN
     &                                       *(CEXT - CNOD)

                                      DERC(I,IPOS(IC),INEW) =
     &                                             DERC(I,IPOS(IC),INEW)
     &                                             - DERIV(IC)*VAUX2

                                  END IF !IBF.EQ.1.AND. ...

                              END IF !IBT.NE.1

                          END DO !K=1,NNUD

                      END DO !IC=1,NPTOT

                  END IF !NPTOT.NE.0

              END IF !IP.NE.0 .OR. IOTRLI.NE.0

          END IF !NZC.NE.0 .AND. NZA.NE.0

      END DO !L=1,NUMEL

C------------------------- Boundary terms

      DO I=1,NUMNP

          IBT = IBTCO(I)
          IBF = IBCOD(I)
          CAUD = CAUDAL(I)
          FACTOR = 0D0

C------------------------- If the node has boundary condition...

          IF (IBT.GT.0 .AND. IBT.LT.5) THEN

              NZC = IXPARNP(I,INCON)
              JJ = INORPAR(18)+NZC
              IP = IVPAR(JJ)

C------------------------- the derivative of the parameter is computed
C------------------------- if needed

              IF (IP.NE.0 .OR. IOTRLI.NE.0) THEN

                  IF (IOTRLI.NE.0) THEN
                      NCNF = NFNLPAR(JJ)
                  ELSE
                      NCNF = 0
                  END IF !IOTRLI.NE.0

                  CFPARAM(1) = CFPARNP(I,7+INDSSTR)
                  INDEX(1) = JJ

C------------------------- Derivatives in nodes with prescribed
C------------------------- concentration boundary condition must be
C------------------------- computed at k+1.

                  DTIMAUX = DTIM + (1D0-THETAT)*TINC/TINTERVOBS 

                  CALL DER_PARAM
     &    (CFPARAM  ,DTIM     ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,0        ,IOFMLT   ,IP       ,L
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,1        ,NUMEL
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     &    ,PARACD   ,PARC(INORPAR(18)+1),CCALIT   ,CCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

                  IF (NPTOT.GT.0) THEN

                      SELECT CASE(IBT)

                          CASE(1) ! Prescribed concentration

C------------------------- Sets previous contributions to zero.

                              DERC(I,1:NPAR,INEW) = 0D0

                              FACTOR = 1D0

                          CASE(2,3) ! Mass flux

                              IF (CAUD.GT.0) THEN

                                  CEXT = PARNP(I,4)
                                  CNOD = CAUX1(I)
                                  CEXT_CINT = CEXT - CNOD

                                  IF (IODENS.EQ.1) THEN

                                      DENSEXT =
     &                                     DENS(DENSREF,BETAC,CEXT,CREF)

                                      FACTOR = DENSEXT*CAUD
     &                                         *(1D0 + BETAC*CEXT_CINT)

C------------------------- Implicit contribution to nodes with
C------------------------- prescribed head through nodal flux.

                                      IF (IBF.EQ.1) THEN

                                          FACTOR = FACTOR
     &                                     -BETAC*CAUD*DENSEXT*CEXT_CINT

                                      END IF !IBF.EQ.1

                                  ELSE

                                      FACTOR = CAUD

                                  END IF !IODENS.EQ.1


                              ELSE

                                      FACTOR = 0D0

                              END IF !CAUDAL(I).GT.0

                          CASE(4) !Input mass

                              FACTOR = 1D0

                          CASE(5) ! Concentration leakage

                              IF (ITPTVAR.EQ.0) THEN !Solute transport

                                  FACTOR = PARNP(I,6)

                              ELSE IF (ITPTVAR.EQ.1) THEN !heat transport

                                  FACTOR = PARNP(I,6)/WSPECHEAT

                              END IF !ITPTVAR.EQ.0,1


                      END SELECT !IBT

                      DO IC=1,NPTOT

                          DERC(I,IPOS(IC),INEW) = DERC(I,IPOS(IC),INEW)
     &                                          + FACTOR*DERIV(IC)
                      END DO !IC=1,NPTOT

                  END IF !NPTOT.GT.0

              END IF !IP.NE.0

          END IF !IBT.GT.0 .AND. IBT.LT.5

      END DO !I=1,NUMNP

      END SUBROUTINE DERCOE
