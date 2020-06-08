      SUBROUTINE DERQ_ARR
     &          (AREA     ,BETAC    ,CAUDAL   ,CAUX1    ,CFPAREL
     &          ,CREF     ,DENSREF  ,DERC     ,DTIM     ,EPSFLU
     &          ,FNT      ,HCALAN   ,HCALIT   ,IBCOD    ,IBTCO
     &          ,IDIMDERC ,IDIMFNT  ,IFLAGS   ,INDSSTR  ,INEW
     &          ,INORPAR  ,INTI     ,IOCAP    ,IODENS   ,IOFLLI
     &          ,IOFMLF   ,IVPAR    ,KXX      ,LMXNDL   ,LNNDEL
     &          ,LXARR    ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG
     &          ,NFNLTIP  ,NFTPAR   ,NINT     ,NPAR     ,NPAREL
     &          ,NPPEL    ,NPPNP    ,NPZON    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZONE_PAR,NZPAR    ,PARACD   ,PARC
     &          ,PAREL    ,PARNP    ,IDIMWGT  ,WGT_PAR  ,IPNT_PAR
     &          ,IPOS     ,DERIV)

********************************************************************************
*
* PURPOSE
*
*  Derivative of nodal flow w.r.t. areal recharge (implicit dependence)
*
*
* DESCRIPTION
*
*  Derivative of nodal flow w.r.t. areal recharge (implicit dependence). This
*  dependence is added to the implicit dependence
*
*  Only for prescribed head nodes with mass flux boundary condition.
*
*  The contribution of recharge to CAUDAL is:
*
*         Qr = ( DESNRECH*RECH*AREA(L)/NNUD ) / DENSEXT
*
*  It is divided by DENSEXT in orer to obtrain a real flow (not a mass flow).
*
*  So, the contribution to DERC is @Qr/@p * DENSEXT*(CEXT- CINT).
*
* EXTERNAL VARIABLES: ARRAYS
*
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  CAUDAL                 Input/output flow at every node.                      
*  CAUX1                  Array containing concentrations, ponderated by THETAT 
*                         time factor                                           
*  DERC                   Nodal concentration derivatives with respect to       
*                         estimated parameters.                                 
*  IBCOD                  Flow boundary condition index                         
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
*  LXARR                  Areal recharge (steady) zone number at a given element
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
*  LMXNDL                 Maximum number of nodes per element                   
*  NFLAGS                 Maximum number of allowed flags                       
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
*  DER_ARR                Derivative of areal recharge w.r.t zonal recharge
*  NNUD                   Number of nodes of the current element                
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
* HISTORY
*
*     AMS      3-2002     First coding (starting from TRANSIN-II)
*
********************************************************************************

       IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMDERC ,IDIMFNT  ,INDSSTR  ,INEW     ,INTI
     &          ,IOCAP    ,IODENS   ,IOFLLI   ,IOFMLF   ,LMXNDL
     &          ,NFLAGS   ,NFNL     ,NINT     ,NPAR     ,NPAREL
     &          ,NPPEL    ,NPPNP    ,NPZON    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZPAR    ,IDIMWGT  ,IPNT_PAR

	REAL*8::BETAC  ,CREF   ,DENSREF,DTIM   ,EPSFLU

      REAL*8::DENS

      INTEGER*4::IBCOD(NUMNP)        ,IBTCO(NUMNP)    ,IFLAGS(NFLAGS)
     &          ,INORPAR(NTYPAR)     ,IVPAR(NZPAR)    ,KXX(LMXNDL,NUMEL)
     &          ,LNNDEL(NUMEL)       ,LXARR(NUMEL)    ,NFNLPAR(NZPAR)
     &          ,NFNLPRG(8,NFNL)     ,NFNLTIP(NFNL)   ,NFTPAR(NZPAR)
     &          ,NZONE_PAR(NTYPAR) 


      REAL*8::AREA(NUMEL)               ,CAUDAL(NUMNP)
     &       ,CAUX1(NUMNP)              ,CFPAREL(NUMEL,NPAREL)
     &       ,DERC(NUMNP,NPAR,IDIMDERC) ,FNT(IDIMFNT,NINT)
     &       ,HCALAN(NUMNP)             ,HCALIT(NUMNP)
     &       ,PARACD(3,NFNL)            ,PARC(NZPAR)
     & ,PAREL(NUMEL,NPPEL)  ,PARNP(NUMNP,NPPNP)    ,WGT_PAR(IDIMWGT)



      
C------------------------- Internal

      INTEGER*4::I      ,IB     ,IBT    ,INODE  ,INARR  ,IP     ,IPAR
     &          ,JJ     ,L      ,NCNF   ,NZARR  ,NNUD   ,NPTOT

	REAL*8::AREALN   ,CAUD     ,CEXT     ,CINT     ,CREC     ,DENSREC
     &       ,FACTOR   ,RECHARGE

	INTEGER*4::INDEX(12),IPOS(NPAR)

      REAL*8::CFPARAM(12),DERIV(NPAR),XPARAM(8)
      
C------------------------- First executable statement    

      INARR = 3 + INDSSTR

      DO L=1,NUMEL

          NZARR = LXARR(L)
          NNUD = LNNDEL(L)
	    JJ = INORPAR(8)+NZARR
          IP = IVPAR(JJ)
          RECHARGE = PAREL(L,8)
          AREALN = AREA(L)/NNUD


          IF (NZARR.NE.0 .AND. (IP.NE.0 .OR. IOFLLI.NE.0)) THEN

              IF (IOFLLI.NE.0) THEN
                  NCNF = NFNLPAR(JJ)
	        ELSE
	            NCNF = 0
	        END IF !IOFLLI.NE.0

              INDEX(1) = JJ
	        CFPARAM(1) = CFPAREL(L,INARR)

              CALL DER_PARAM
     &    (CFPARAM  ,DTIM     ,EPSFLU   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,0        ,IOFMLF   ,IP       ,L
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,NPZON    ,NUMEL
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     &    ,PARACD   ,PARC(INORPAR(19)+1),HCALIT   ,HCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)


              IF (NPTOT.GT.0) THEN
              
                  DO I=1,NNUD

                      INODE = KXX(I,L)
	                IB = IBCOD(INODE)
	                IBT = IBTCO(INODE)
	                CAUD = CAUDAL(INODE)

	                IF (IB.EQ.1 .AND. CAUD.GT.0 .AND.
     &                   (IBT.EQ.2 .OR. IBT.EQ.3) ) THEN

                          CEXT = PARNP(INODE,4)
                          CINT = CAUX1(INODE)

                          IF (RECHARGE.GT.0) THEN

	                        CREC = PAREL(L,15)
	                        
	                    ELSE

                              CREC = 0D0

	                    END IF !RECHARGE.GT.0

	                    IF (IODENS.EQ.1) THEN

                              DENSREC = DENS(DENSREF,BETAC,CREC,CREF)

                          ELSE

                              DENSREC = 1D0

                          END IF !IODENS.EQ.1


	                    FACTOR = AREALN*DENSREC*(CEXT - CINT)

                          DO IPAR=1,NPTOT

                              DERC(INODE,IPOS(IPAR),INEW) =
     &                                       DERC(INODE,IPOS(IPAR),INEW)
     &                                       - DERIV(IPAR)*FACTOR

                             IF (IFLAGS(30).EQ.1) THEN

                                  WRITE(77,10)
     &                            'Q_ARR',INODE,L,IPOS(IPAR),DERIV(IPAR)

   10                             FORMAT (A,3I6,E20.10)

	                        END IF !IFLAGS(30).EQ.1

                          END DO !IPAR=1,NPTOT

                      END IF !IBCOD(INODE).EQ.1. ...
             
                  END DO !I=1,NNUD

	        END IF !NPTOT.GT.0

          END IF !NZ.NE.0 .AND. IP.NE.0

      END DO !L=1,NUMEL

      END SUBROUTINE DERQ_ARR
