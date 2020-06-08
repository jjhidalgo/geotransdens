      SUBROUTINE DERQ_TRA 
     &          (AREA     ,BIBI     ,BUOYANCY ,CAUDAL   ,CAUX1
     &          ,CFPAREL  ,COORD    ,DENSITY  ,DERC     ,DTIM
     &          ,EPSFLU   ,FNT      ,GP_COORD ,GRADLOC  ,HAUX1
     &          ,HBASE    ,HCALAN   ,HCALIT   ,IBCOD    ,IDIMBB
     &          ,IDIMFNT  ,IFLAGS   ,INDSSTR  ,INORPAR  ,INTI
     &          ,IOCAP    ,IOCTRA   ,IODENS   ,IODIM    ,IOFLLI
     &          ,IOFMLF   ,ISOZ     ,IVPAR    ,KXX      ,LDIM
     &          ,LMXNDL   ,LNNDEL   ,LTYPE    ,LXPAREL  ,MAXPG
     &          ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &          ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   ,NPPNP
     &          ,NTYPAR   ,NTYPEL   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &          ,NZPAR    ,NZTRA    ,PARACD   ,PARC     ,PARNP
     &          ,POINTWEIGHT        ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR
     &          ,IPOS     ,DERIV)

*******************************************************************************
*
* PURPOSE
*
*    Computes derivative of flow w.r.t transmissivity (only direct dependence)
*
* DESCRIPTION
*
*    Computes derivative of flow w.r.t transmissivity (only direct dependence)
*    The indirect dependence (through heads) is computed in subroutine DER_VD
*    Another difference with DER_VD is in the fact that DER_VD computes the
*    derivative of VD w.r.t a UNIQUE zonal parameter, while this routine
*    includes the direct dependence w.r.t ALL transmissivities. OBSERVE that 
*    this dependence is only necessary in type 1 boundary condition nodes
*
* EXTERNAL VARIABLES: ARRAYS
*
*  BIBI                   Array containing the product of the gradient of          
*                         interpolation functions, for a given element               
*  CAUDAL                 Input/output flow at every node.                      
*  CAUX1                  Array containing concentrations, ponderated by THETAT 
*                         time factor                                           
*  DERC                   Nodal concentration derivatives with respect to       
*                         estimated parameters.                                 
*  HAUX1                  Array containing heads, ponderated by THETAF time     
*                         weight. Is equal to THETAF*HCAL+(1-THETAF)*HCALAN     
*  IBCOD                  Flow boundary condition index                         
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  ISOZ                   Anisotropy of every transmissivity zone               
*  IVPAR                  Vector containing estimation index for all            
*                         parameters                                            
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*  LTYPE                  Vector containing the type of each element            
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
*  IND                    Auxiliar array used to simplify the product           
*                         T x gradient of h (used to compute Darcy's velocity)  
*                         Basically allows to transform the vector storage      
*                         of T in matrix storage (IODIM x IODIM)                
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMBB                 Used to dimension array BIBI. Is equal to IDIMQ times 
*                         the maximum possible anisotropy of the problem        
*  IODIM                  Maximum dimension of the problem                      
*  LMXNDL                 Maximum number of nodes per element                   
*  NPAR                   Total number of parameters to be estimated            
*  NPAREL                 Number of element parameters in current problem       
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
*  NZTRA                  Number of transmissivity zones                        
*
* INTERNAL VARIABLES: SCALARS
*
*  IC                     Counter to position in array BIBI
*  IP                     Number of parameter estimation
*  ISZL                   Anisotropy of current element
*  LTY                    Type of current element
*  NNUD                   Number of nodes of the current element                
*  NZT                    Transmissivity zone number of current element
*  PAR_ZON                Zonal transmissivity value at zone NZT
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*
*   None
*
* HISTORY
*
*     AMS      4-2002     First coding
*
*******************************************************************************

      IMPLICIT NONE

C--------------------  External

      INTEGER*4::IDIMBB   ,IDIMFNT  ,INDSSTR  ,INTI     ,IOCAP
     &          ,IOCTRA   ,IODENS   ,IODIM    ,IOFLLI   ,IOFMLF
     &          ,LMXNDL   ,MAXPG    ,NFLAGS   ,NFNL     ,NINT
     &          ,NPAR     ,NPAREL   ,NPPNP    ,NTYPAR   ,NTYPEL
     &          ,NUMEL    ,NUMNP    ,NZPAR    ,NZTRA

      REAL*8::DTIM     ,EPSFLU

      INTEGER*4::IBCOD(NUMNP)         ,IFLAGS(NFLAGS) ,INORPAR(NTYPAR)
     &          ,ISOZ(NZTRA)          ,IVPAR(NZPAR)   ,KXX(LMXNDL,NUMEL)
     &          ,LDIM(NUMEL)          ,LNNDEL(NUMEL)  ,LTYPE(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL),NFNLPAR(NZPAR) ,NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL)        ,NFTPAR(NZPAR)  ,NZONE_PAR(NTYPAR)
     ;          ,IDIMWGT  ,IPNT_PAR

      REAL*8::AREA(NUMEL),BIBI(IDIMBB,NUMEL)
     &       ,BUOYANCY(IODIM,LMXNDL,NUMEL),CAUDAL(NUMNP),CAUX1(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL),COORD(NUMNP,IODIM)
     &       ,DENSITY(NUMEL),DERC(NUMNP,NPAR)
     &       ,GP_COORD(6,8,IODIM),GRADLOC(IODIM,LMXNDL,MAXPG)
     &       ,FNT(IDIMFNT,NINT),HAUX1(NUMNP),HBASE(NUMEL)
     &       ,HCALAN(NUMNP),HCALIT(NUMNP)
     &       ,PARACD(3,NFNL),PARC(NZPAR),PARNP(NUMNP,NPPNP)
     &       ,POINTWEIGHT(6,8)    ,WGT_PAR(IDIMWGT)
     &       
     
C--------------------  Internal

      INTEGER*4::I      ,IANIS  ,IC     ,INODE  ,IP     ,IPAR   ,ISZL
     &          ,ITRA   ,J      ,JJ     ,JNODE  ,L      ,LD     ,LTY
     &          ,NCNF   ,NNUD   ,NPTOT  ,NZT

      REAL*8::CEXT_CI  ,CEXTNODE ,CNODE    ,D_TRA
     &       ,HJHI     ,PAR_ZON  ,SX       ,SX_BYNCY

      INTEGER*4::IND(6,5,6),INDEX(12),IPOS(NPAR)


      REAL*8::CFPARAM(12),DERIV(NPAR),DERIVTRA(6),XPARAM(8)

      DATA ((IND(I,J,6),J=1,5),I=1,6)/0,6,12,18,24,  0,30,36,42,48,
     &       6,30,54,60,66,  12,36,54,72,78,  18,42,60,72,84,  
     &       24,48,66,78,84/
      DATA ((IND(I,J,4),J=1,3),I=1,4)/0,6,12,0,18,24,6,18,30,12,24,30/
      DATA ((IND(I,J,3),J=1,3),I=1,4)/0,3,6,0,9,12,3,9,15,6,12,15/
      DATA ((IND(I,J,5),J=1,3),I=1,3)/0,0,3,0,0,6,3,6,0/
      DATA ((IND(I,J,2),J=1,3),I=1,3)/0,0,3,0,0,6,3,6,0/
      DATA ((IND(I,J,1),J=1,1),I=1,1)/0/

C--------------------  First Executable Statement

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          NZT = LXPAREL(L,1)
          ISZL = ISOZ(NZT)
          LTY = LTYPE(L)
          LD = LDIM(L)

          DO I=1,NNUD

              INODE = KXX(I,L)
              

C------------------------- Contribution only in nodes with prescribed
C------------------------- head (IBCOD=1) and out-flow

              IF (IBCOD(INODE).EQ.1 .AND. CAUDAL(INODE).GT.0) THEN

                  CEXTNODE = PARNP(INODE,4)
                  CNODE = CAUX1(INODE)
                  CEXT_CI = CEXTNODE - CNODE

                  DO ITRA=1,ISZL

                      IP = IVPAR(NZT+INORPAR(ITRA))

                      IF (IP.GT.0 .OR. IOFLLI.NE.0) THEN

                          JJ = INORPAR(ITRA)+NZT

                          IF (IOFLLI.NE.0) THEN

                              NCNF = NFNLPAR(JJ)
                              XPARAM(2) = HBASE(L)

                          ELSE

                              NCNF = 0

                          END IF !IOFLLI.NE.0

                          INDEX(1) = JJ

                          CFPARAM(1) = CFPAREL(L,1)

                          PAR_ZON = PARC(INORPAR(ITRA)+NZT)

                          CALL DER_PARAM
     &    (CFPARAM  ,DTIM     ,EPSFLU   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,IOCTRA   ,IOFMLF   ,IP       ,L        
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,1        ,NUMEL    
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT      
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG  
     &    ,PARACD   ,PARC(INORPAR(19)+1),HCALIT   ,HCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)
                      
                          DO IPAR=1,NPTOT

                              D_TRA=DERIV(IPAR)

                              SX = 0D0

                              DO J=1,NNUD

                                  IF (I.NE.J) THEN

                                      IC = IND(I,J,LTY) 
                                      JNODE = KXX(J,L)
                                      HJHI = HAUX1(JNODE) - HAUX1(INODE)
                                  
C------------------------- Takes isotropy into account (if ISZL < IODIM, 
C------------------------- Txx is repeated in transmissivity tensor)

                                      IF(ITRA.EQ.1 .AND.
     &                                   ISZL.LT.IODIM) THEN

                                          DO IANIS=0,IODIM-ISZL

                                              SX = SX + D_TRA
     &                                       *BIBI(IC+ITRA+IANIS,L)*HJHI
     
                                          END DO !IANIS=0,IODIM-ISZL

                                      ELSE

C------------------------- If ITRA <> 1 or anisotropy is larger or equal to 
C------------------------- element's dimension, transmissivity appears only once

                                          SX = SX
     &                                       +D_TRA*BIBI(IC+ITRA,L)*HJHI

                                      END IF !ITRA.EQ.1 .AND. ISZL.LT.IODIM

                                  END IF !I.NE.J

                              END DO ! J=1,NNUD

                              IF (IODENS.EQ.1) THEN

C------------------------- ATRA matrix is multiplied by
C------------------------- the average density of the element

                                  SX = DENSITY(L)*SX

C------------------------- Buoyancy contribution to density-dependent flow

                                  CALL DERTRA_TENSOR
     &                                (DERIVTRA,ISZL,ITRA,LD)

                                  DERIVTRA = DERIV(IPAR)*DERIVTRA

                                  CALL DERQ_TRA_BUOYANCY
     &              (AREA     ,BUOYANCY ,COORD    ,DENSITY  ,DERIVTRA
     &              ,GP_COORD ,GRADLOC  ,IODIM    ,ISOZ     ,KXX
     &              ,L        ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE
     &              ,LXPAREL  ,MAXPG    ,I        ,NPAREL   ,NTYPEL
     &              ,NUMEL    ,NUMNP    ,NZTRA    ,POINTWEIGHT
     &              ,SX_BYNCY)

                                  SX = SX + SX_BYNCY

                              END IF !IODENS.EQ.1

C--------------------------- Adds to inverse problem RHS 
C--------------------------- The fluid flow computed already contains
C--------------------------- the external density contribution, so it
C--------------------------- is only necessary to multiply by the
C--------------------------- internal and external concentration
C--------------------------- difference

                              DERC(INODE,IPOS(IPAR))=
     &                               DERC(INODE,IPOS(IPAR)) + SX*CEXT_CI

C-------------------------    WRITE(78,*) ' QT',INODE,L,SX

                          END DO !IPAR=1,NPTOT

                      END IF !IP.GT.0 .OR. IOFLLI.NE.0

                  END DO ! ITRA=1,ISZL (element anisotropy)

              END IF ! IBCOD(INODE).EQ.1 .AND. CAUDAL(INODE).GT.0

          END DO ! I=1,NNUD   (nodes in the current element)

      END DO ! L=1,NUMEL  (elements)


      END SUBROUTINE DERQ_TRA
