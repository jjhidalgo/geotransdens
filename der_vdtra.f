      SUBROUTINE DER_VD_TRA
     &          (AREA     ,BUOYANCY ,CCAL     ,CFPAREL  ,COORD
     &          ,DAT_VD   ,DENSITY  ,DER_VD   ,DTIM     ,EPSFLU
     &          ,FNT      ,GP_COORD ,GRADLOC  ,GRDFF    ,HAUX1
     &          ,HBASE    ,HCALAN   ,HCALIT   ,IDIMDQ   ,IDIMFNT
     &          ,IFLAGS   ,INDSSTR  ,INORPAR  ,INTI     ,IOCAP
     &          ,IOCTRA   ,IODENS   ,IODIM    ,IOFLLI   ,IOFMLF
     &          ,ISOZ     ,IVPAR    ,KXX      ,LDIM     ,LMXNDL
     &          ,LNNDEL   ,LTYPE    ,LXPAREL  ,MAINF    ,MAXPG
     &          ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &          ,NFTPAR   ,NINT     ,NPAR     ,NPAREL   ,NTYPAR
     &          ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR    ,NZTRA
     &          ,PARACD   ,PARC     ,POINTWEIGHT        ,RHS_IP
     ;          ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV)


*****************************************************************************
*
* PURPOSE
*
*     Computes the derivative of ATRA matrix w.r.t transmissivity (direct 
*     dependence)
*
* DESCRIPTION
*
*     In the case of the derivatives of ATRA w.r.t transmissivity, it is
*     divided in two parts. First part is the indirect dependence of VD 
*     w.r.t. T (dependence through head), that is an analogous dependence 
*     as for the rest of parameters, and this is computed with two calls
*     DER_VD and DER_ATRA_FP_VD. The second part, direct dependence of VD
*     w.r.t. T is added in this subroutine. The first part of this routine
*     is quite similar to DER_VD, and the second part is similar to 
*     DER_ATRA_FP_VD
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  CCAL                   Computed concentration at every node                  
*  COORD                  Nodal coordinates                                     
*  DAT_VD                 Derivatives of ATRA w.r.t. the components of Darcy's  
*                         velocity (qx, qy, qz).                                
*  DER_VD                 Derivatives of Darcy's velocity with respect to       
*                         flow parameters. In some cases it is dimensioned      
*                         as (IODIM,NUMEL) to reduce storage.                   
*  GRAVEL                 Projection of gravity at each element                 
*  GRDFF                  Array containing the product between interpolation    
*                         functions integrals and interp. functions gradient    
*  HCAL                   Computed heads at every node                          
*  HEAD                   Auxiliar array used to store nodal heads              
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
*  LXPAREL                Array containing zone numbers for a given             
*                         element parameter                                     
*  PARC                   Vector containing calculated values for all           
*                         parameters                                            
*  PAREL                  Parameter values at every element and current time    
*                         for all nodal parameters (each value is computed as   
*                         the product of up to four terms:                      
*                           elem. coeff*zonal value*time funct.*nonl. funct. )  
*  RHS_IP                                                                       
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
*  IDIMDQ                 Used to dimension array DAT_VD (second dimension). It 
*                         is the number of different terms (ATRA)ij varying     
*                         i and  j in the local numeration of one element.      
*                         Is equal to LMXNDL*(LMXNDL-1)/2                       
*  IODIM                  Maximum dimension of the problem                      
*  IOPRHED                Indicates whether the flow state variable state is    
*                         preasure (set to 1) or head (set to 0)                
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFLAGS                 Maximum number of allowed flags                       
*  NPAR                   Total number of parameters to be estimated            
*  NPAREL                 Number of element parameters in current problem       
*  NPPEL                  Total number of parameters by elements (not confuse   
*                         with NPAREL, because in this case, different          
*                         anisotropy terms are treated separatedly)             
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
*  IP                     Index of current parameter to be estimated
*  ISZL                   Anisotropy of the current element
*  LD                     Dimension of current element
*  NNUD                   Number of nodes of the current element                
*  NZT                    Transmissivity zone of current element
*  PAR_ZON                Zonal value of transmissivity of current element
*  ZERO                   Constant set to zero                                  
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  CALHEAD                Transforms pressure to head level                     
*
*  HISTORY
*
*     AMS      1-2002     First coding
*
*****************************************************************************

       IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMDQ   ,IDIMFNT  ,INDSSTR  ,INTI     ,IOCAP
     &          ,IOCTRA   ,IODENS   ,IODIM    ,IOFLLI   ,IOFMLF
     &          ,LMXNDL   ,MAINF    ,MAXPG    ,NFNL     ,NFLAGS
     &          ,NINT     ,NPAR     ,NPAREL   ,NTYPAR
     &          ,NUMEL    ,NUMNP    ,NZPAR    ,NZTRA    ,IDIMWGT

      REAL*8:: DTIM,EPSFLU    ,WGT_PAR(IDIMWGT)

      INTEGER*4::IFLAGS(NFLAGS)                 ,ISOZ(NZTRA)
     &          ,INORPAR(NTYPAR),IVPAR(NZPAR)   ,KXX(LMXNDL, NUMEL)
     &          ,LDIM(NUMEL)    ,LNNDEL(NUMEL)  ,LTYPE(NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL)          ,NFNLPAR(NZPAR)
     &          ,NFNLPRG(8,NFNL),NFNLTIP(NFNL)  ,NFTPAR(NZPAR)
     &          ,NZONE_PAR(NTYPAR)

      REAL*8::AREA(NUMEL)                  ,BUOYANCY(IODIM,LMXNDL,NUMEL)
     &       ,CCAL(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL)        ,COORD(NUMNP,3)
     &       ,DAT_VD(IODIM,IDIMDQ,NUMEL)   ,DENSITY(NUMEL)
     &       ,DER_VD(IODIM)
     &       ,FNT(IDIMFNT,NINT)            ,GP_COORD(6,8,IODIM)
     &       ,GRADLOC(IODIM,LMXNDL,MAXPG)  ,GRDFF(IODIM,LMXNDL,NUMEL)
     &       ,HAUX1(NUMNP)                 ,HBASE(NUMEL)
     &       ,HCALAN(NUMNP)                ,HCALIT(NUMNP)
     &       ,PARACD(3,NFNL)               ,PARC(NZPAR)
     &       ,POINTWEIGHT(6,8)             ,RHS_IP(NUMNP,NPAR)
     

C------------------------- Internal

      INTEGER*4::I    ,IC   ,ILD  ,IND_IK_LD  ,INODE,IP   ,IPAR ,ISZL
     &          ,ISZL1,ISZL2,ITRA ,J    ,JJ   ,JNODE,K    ,L    ,LD
     &          ,NCNF ,NNUD ,NZT  ,NPTOT      ,IPNT_PAR

      REAL*8::AREALN,PAR_ZON,S,S1,S2,T1,T2,T3,ZERO

      INTEGER*4::IND(3, 3, 3),INDEX(12),IPOS(NPAR)

      REAL*8::CFPARAM(12),DERIV(NPAR),DERIVTRA(6),XPARAM(8)

C------------------------- Array IND is used to simplify the coding of 
C------------------------- conductivity matrix times the gradient of FEM 
C------------------------- interpolation functions (shape functions)

      DATA ((IND(I,J,3),I=1,3),J=1,3)/1,4,5,4,2,6,5,6,3/
      DATA ((IND(I,J,2),I=1,2),J=1,2)/1,3,3,2/,IND(1,1,1)/1/
      DATA ZERO /0.D0/

C------------------------- FIRST EXECUTABLE STATEMENT.

C------------------------- Cross over elements

      DO L=1,NUMEL

          NNUD = LNNDEL(L)
          LD = LDIM(L)
          NZT = LXPAREL(L,1)
          ISZL = ISOZ(NZT)
          ISZL1 = ISZL + 1
          ISZL2 = ISZL + 2
          AREALN = AREA(L)/NNUD

          DO ITRA=1,ISZL

              IP = IVPAR(INORPAR(ITRA)+NZT)

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

C------------------------- Computes partial derivatives of velocities w.r.t.
C------------------------- transmissivity (only direct dependence)

                      DO IPAR=1,NPTOT

                          DO I=1,LD

                          S = ZERO

                          DO J=1,NNUD 

                              JNODE = KXX(J, L)

                              S2 = ZERO

                              DO K=1,LD

                                  IND_IK_LD = IND(I,K,LD)

                                  IF ((ISZL2.EQ.LD .AND. IND_IK_LD.LE.3)
     &                            .OR.
     &                                (ISZL1.EQ.LD .AND. IND_IK_LD.LE.2)
     &                            .OR.
     &                             (ISZL.GE.LD  .AND. IND_IK_LD.EQ.ITRA)
     &                            ) THEN

                                      S2 = S2 + DERIV(IPAR)*GRDFF(K,J,L)

                                  END IF !ISZL...

                              END DO !K=1,LD

                              S = S - HAUX1(JNODE)*S2

                          END DO !J=1,NNUD

                          DER_VD(I) = S

                      END DO !I=1,LD
                  
                      IF (IODENS.EQ.1) THEN

C------------------------- Builds the derivative of transmissivity tensor
C------------------------- w. r. t. to the current direction of anisotropy

                          CALL DERTRA_TENSOR(DERIVTRA,ISZL,ITRA,LD)

                          DERIVTRA = DERIV(IPAR)*DERIVTRA

C------------------------- Contribution of buoyancy term.

                          CALL DER_VDTRA_BUOYANCY
     &                   (AREA     ,BUOYANCY ,COORD    ,DERIVTRA ,DER_VD
     &                   ,GP_COORD ,GRADLOC  ,IODIM    ,KXX      ,L
     &                   ,LDIM     ,LMXNDL   ,LNNDEL   ,LTYPE    ,MAXPG
     &                   ,NUMEL    ,NUMNP    ,POINTWEIGHT)

                      END IF !IODENS.EQ.1

                      IF (IFLAGS(33).EQ.1) THEN

                         WRITE(MAINF,10) ' VDT',L,(DER_VD(ILD),ILD=1,LD)

   10                    FORMAT(A4,I5,3E20.13)

                      END IF !IFLAGS(33).EQ.1


C------------------------- Assembles Vd derivatives with DAT_VD to compute
C------------------------- derivatives of ATRA w.r.t. transmissivity
C------------------------- (only direct dependence)

C------------------------- Since the whole ATRA matrix is multiplied by the
C------------------------- density of the element, it is better to multyply
C------------------------- the velocity vector by the density  now to save
C------------------------- operations instead of multiply inside one of the
C------------------------- following loops.

                      IF (IODENS.EQ.1) THEN

                          DER_VD = DENSITY(L)*DER_VD

                      END IF !IODENS.EQ.1

                      IC = 1

                      DO I=1,NNUD-1
      
                          INODE = KXX(I,L)

                          DO J=I+1,NNUD

                              JNODE = KXX(J,L)
                              T1 = ZERO
                              T2 = ZERO
                              T3 = ZERO

                              DO K=1,LD

                                  T1 = T1 +DAT_VD(K,IC,L)*DER_VD(K)
                                  T2 = T2 +AREALN*GRDFF(K,J,L)*DER_VD(K)
                                  T3 = T3 +AREALN*GRDFF(K,I,L)*DER_VD(K)

                              END DO !K=1,LD
                           
                              S1 = (T1 + T2)*(CCAL(JNODE) - CCAL(INODE))
                              S2 = (T1 + T3)*(CCAL(INODE) - CCAL(JNODE))

                              RHS_IP(INODE,IPOS(IPAR)) =
     &                                     RHS_IP(INODE,IPOS(IPAR)) - S1

                              RHS_IP(JNODE,IPOS(IPAR)) =
     &                                     RHS_IP(JNODE,IPOS(IPAR)) - S2

                              IC = IC + 1

                          END DO !J=I+1,NNUD

                      END DO ! I=I,NNUD-1

                  END DO !IPAR=1,NPTOT

              END IF ! IP.NE.0

          END DO ! ITRA=1,ISOZ

      END DO !L=1,NUMEL

C------------------------- Writes derivatives of velocities w.r.t. T
C------------------------- (only direct dependence)
C-----IF (IFLAGS(23).EQ.1) THEN
C-----     DO L = 1, NUMEL
C-----        WRITE(MAINF, 2000) (VD(K, L), K = 1, IDIM)
C-----     ENDDO
C-2000     FORMAT(6E21.12)
C-----  ENDIF

      END SUBROUTINE DER_VD_TRA
