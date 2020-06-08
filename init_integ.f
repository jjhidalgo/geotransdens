      SUBROUTINE INIT_INTEG 
     ;(IDIMAFLU ,IDIMBB   ,IDIMCOV  ,IDIMCROSS_GS ,IDIMDFLU   ,IDIMDQ
     ;,IDIMDTRA ,IDIMFNT  ,IDIMHESS ,IDIMQ    ,IDIMVAR_GS ,IOCNSF
     ;,IOCNST   ,IODENS   ,IODIM   ,IOEQT    ,IOFLLI     ,IOFLSAT
     ;,IOINV    ,IOPART   ,IORTS    ,IOTRLI   ,IOTRS
     ;,IPARTRA  ,ISOT     ,LMXNDL   ,MAINF     
     ;,MAXNEOP
     ;
     ;,NBAND    ,NBAND1   ,NBAND2   ,NBANDCOV   ,NDEVS
     ;,NFNT     ,NINT
     ;,NPAR     ,NPAREL   ,NPARNP   ,NTYPAR     ,NUMEL
     ;,NUMNP    ,NUMTIT   ,NUMTNOD  ,NUMTOBS  ,NZALF
     ;,NZARR    ,NZCHP    ,NZCOE    ,NZCRD    ,NZDFM 
     ;,NZDMT    ,NZDSP    ,NZFOD    ,NZPAR    ,NZPOR
     ;,NZPRG    ,NZQQP    ,NZSTG    ,NZTRA    ,INORPAR
     ;,IDIMWORK ,IOSMTP   ,IOSMFL   ,NPBTP    ,NPBFL
     ;!nuevos
     ;,IDIMCFLU ,IDIMATRA
     ;,IOSPARSE ,IAFLUDSC_ROWS, IAFLUDSC_COLS,IATRADSC_ROWS      
     ;,IATRADSC_COLS  ,IA_COUPL_ROWS ,IA_COUPL_COLS
     ;,ICAN_CN  ,IALW_CN  ,IPAR_DIR ,NPARALG, ITYPAFLU, ITYPBFLU
     ;,ITYPCFLU ,ITYPDFLU ,ITYPATRA ,ITYPDTRA
     ;,ITYPBTRA ,ITYPFLUDSC,ITYPTRADSC,ITYPCOUPLDSC, ITYPDERIV
     ;,IDIMDENS ,IOVRWC ,IDIMGRAVEL,NZCLK
     &,IDIMDERH,IDIMDERC
     ; ,IDIMIVARIO_GS  ,IDIMWGT        ,IDIMZONPP_GS   ,IOINV_GS
     ; ,MXCLOSE_GS     ,MXDISC_GS      ,MXGRPZN        ,MXKRIG_GS
     ; ,MXLINCMB       ,MXMEASPP_GS    ,MXNPP_GS       ,MXNPRIM_GS
     ; ,MXNST_GS       ,MXNVAR_GS      ,MXNZON_GS      ,MXROT_GS
     ; ,MXSAM_GS       ,MXSB_GS        ,MXSC_GS        ,MXVGM_GS
     ; ,MXZONPP_GS     ,NGROUP_ZN      ,IOPT_GS        ,IO_KG_GS     
     ; ,IDIMDATASC_GS)

********************************************************************************
*
* PURPOSE Defines some integer variables
*
* DESCRIPTION Defines some variables used to dimension several arrays,
*             and indexes to the location of the different parameters in arrays
*             PARZ, PARC, PARM, IVPAR, NFTPAR, STPAR and FNTPAR. Also, calculates
*             main dimensions of arrays related to geostatistics
*
* EXTERNAL VARIABLES: ARRAYS
*
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR 
*
* INTERNAL VARIABLES: ARRAYS
*
*  INTRAC                 Array containing the location of first transmissivity 
*                         zone (for each tensor component) in array             
*                         variables PARC, PARM, STPAR ... minus 1               
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMAFLU               Used to dimension array AFLU                          
*  IDIMBB                 Used to dimension array BIBI                          
*  IDIMCOV                Used to dimension array COVINV
*  IDIMCROSS              Used to dimesnion CROSSCOV_GS
*  IDIMDADFLU             Used to dimension array DERADFLU
*  IDIMDFLU               Used to dimension array DFLU                          
*  IDIMDQ                 Used to dimension array DADQ (second dimension)       
*  IDIMDTRA               Used to dimension array DTRA (second dimension)       
*  IDIMFNT                First dimension of array FNT, it coincides with       
*                         NFNT if the latter is not zero                        
*  IDIMHESS               Used to dimension array HESS. It is equal to          
*                         NPAR*(NPAR+1)/2                                       
*  IDIMQ                  Used to dimension array QXYZ                          
*  IDIMIVARIO_GS          Used to dimension IVARIO_GS, VARIO_GS
*  IDIMVAR_GS             NVAR_GS*NVAR_GS+NEXDR_GS
*  IDIMWGT                Used to dimension WGT_PAR, IPNT_PAR
*  IDIMWORK               Used to dimension workspace WORK
*  IDIMZONPP_GS           Maximum among maximum number of pilot points and 
*  IOCNSF                 Scheme for storage term in flow problem               
*  IOCNST                 Scheme for mass storage term in transport problem     
*  IODIM                  Maximum dimension of any element included             
*                         in the problem                                        
*  IOEQT                  Type of problem to be solved                          
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  IOFLSAT                Indicates the possibility that one part of the domain 
*                         reaches unsaturated state.                            
*  IOINV_GS               Geostatistical inverse problem option
*  IOINV                  Inverse problem option                                
*  IOPART                 If non zero, writes partition and indexing variables
*  IOPT_GS                General options for inverse problem. 
*                         Each row contains information of a given group of zones
*                         On each row:
*                         - Column 1: Identifier of parameter type
*                         - Column 2: Estimation option. (0: zones belonging to
*                                     this group are estimated deterministically;
*                                     1: geostatistically;2:linear combination)
*                         - Column 3: Number of variables used for krig./cokrig.
*                         - Column 4: Number of measurements used for k/ck
*                         - Column 5: Number of pilot points used for k/ck
*                         - Column 6: Pilot points drawing option
*                                     0: Read and fixed
*                                     1: Completely Random and variable trough
*                                        optimization process
*                                     2: Layed on a random regular mesh and 
*                                        variable (that mesh) trough opt. proc.
*                         - Column 7: Number of zones in this group
*                         - Column 8: Density of pilot points (X-direction)
*                         - Column 9: Density of pilot points (Y-direction)
*                         - Column 10: Density of pilot points (Z-direction)
*                         - Column 11; Number of terms defining linear combination
*                         - Column 12. Flow/transport problem
*  IO_KG_GS               Kriging dimensions and options for geost. inv. prob. 
*                         Each row contains information of a given group of zones 
*                         Only sense if group is estimated geost. On each row:
*                         - Column 1: Kriging/Cokriging type
*                                     0: Simple kriging
*                                     1: Residual kriging
*                                     2: Kriging with locally varying mean
*                                     3: Kriging with external drift (up to four)
*                                     4: Simple cokriging
*                                     5: Standardized ordinary cokriging
*                                     6: Traditional ordinary cokriging
*                         - Column 2: Minimum number of samples to be considered 
*                                     in the kriging estimation. Point/Zone will 
*                                     remain unestimated if number of closest 
*                                     samples is less than this number
*                         - Column 3: Maximum number of primary variable samples 
*                                     used for kriging
*                         - Column 4: Maximum number of secondary samples used for 
*                                     kriging (all secondary variable samples are
*                                     jointly considered)
*                         - Column 5: Number of samples retain per octant. If 0, 
*                                     octant search is not performed
*                         - Column 6: Number of super-blocks in X-direction (super-
*                                     block search mesh)
*                         - Column 7: Number of super-blocks in Y-direction (super-
*                                     block search mesh)
*                         - Column 8: Number of super-blocks in Z-direction (super-
*                                     block search mesh)
*                         - Column 9: Number of discretization points of a zone in
*                                     X-direction.
*                         - Column 10: Number of discretization points of a zone in
*                                      Y-direction.
*                         - Column 11: Number of discretization points of a zone in
*                                      Z-direction.
*                                      Note: If product of last three numbers is 1, 
*                                            then point kriging is performed
*                         - Column 12: Number of nested structures defining the most
*                                      complex variogram of this group of zones 
*                                      (unique variogram if kriging is performed)
*                         - Column 13: Number of external drift terms
*                         - Column 14: Option of conditional/unconditional simulation
*                                      - 0: Conditional Estimation
*                                      - 1: Conditional simulation
*                                      - 2: Unconditional simulation
*                         - Column 15: Seed for conditional/unconditional simulations
*                         - Column 16: Number of conditional simulations
*  IORTS                  Transport regime                                      
*  IOTRLI                 Idem to IOFLLI, in the case of transport.             
*  IOTRS                  Flow regime                                           
*  IPARTRA                Used for dimensioning, it equals                      
*                         MAX(1,NPARTRA*(NPARTRA-1)/2)                          
*  ISOT                   Maximum hydraulic conductivity tensor anisotropy      
*                         degree in the problem.                                
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  MAXNEOP                Used to reserve some space. It is equal to            
*                         MAX (NUMEL,NUMNP,NDEVS,NPAR)
*  NBAND                  Half Bandwith (maximum difference between the         
*                         numbers of two nodes belonging to the same element)   
*  NBAND1                 Used to dimension. It is equal to NBAND+1             
*  NBAND2                 Used to dimension. It is equal to 2*NBAND+1           
*  NBANDCOV               Band width of the inverse covariance matrix           
*  NDEVS                  Number of observation Devices
*  NFNT                   Number of time functions used for describing time     
*                         dependence of all transient parameters                
*  NGROUP_ZN              Number of groups of zones
*  NINT                   Number of observation times                           
*  NMEAS                  Number of sampling locations
*  NPAR                   Total number of parameters to be estimated            
*  NPAREL                 Number of element parameters in current problem       
*  NPARNP                 Number of nodal parameters in current problem         
*  NPARTRA                Number of transmissivity zones to be estimated        
*                         with non-diagonal prior information covariance matrix 
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NUMTIT                 Total number of integration times                     
*  NUMTNOD                Total number of nodes used for calculating obs.       
*  NUMTOBS                Total number of observations                          
*  NVAR                   Number of varioables (prim.+all secondary)
*  NZALF                  Number of leakage zones                               
*  NZARR                  Number of areal recharge zones                        
*  NZCHP                  Number of prescribed head zones                       
*  NZCOE                  Number of external concentration zones                
*  NZCRD                  Number of retardation Coefficient zones               
*  NZDFM                  Number of molecular difusion zones                    
*  NZDMT                  Number of matrix diffusion zones                      
*  NZDSP                  Number of dispersivity zones                          
*  NZFOD                  Number of zones of first order decay                  
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  NZPOR                  Number of porosity zones                              
*  NZPRG                  Total number of generic parameter zones               
*  NZQQP                  Number of prescribed flow zones                       
*  NZSTG                  Number of storage Coefficient zones                   
*  NZTRA                  Number of transmissivity zones                        
*
* INTERNAL VARIABLES: SCALARS
*
*  I                      Counter dummy variable
*  INALFC                 Location of first leakage zone in array               
*                         variables PARC, PARM, STPAR ... minus 1               
*  INARRC                 Location of first areal recharge zone in array        
*                         variables PARC, PARM, STPAR ... minus 1               
*  INCHPC                 Location of first prescribed head zone in array       
*                         variables PARC, PARM, STPAR ... minus 1               
*  INCOEC                 Location of first external concent. zone in array     
*                         variables PARC, PARM, STPAR ... minus 1               
*  INCRDC                 Location of first retardation coeff. zone in array    
*                         variables PARC, PARM, STPAR ... minus 1               
*  INDFMC                 Location of first molecular diffusion zone in array   
*                         variables PARC, PARM, STPAR ... minus 1               
*  INDSLC                 Location of first longitudinal disp. zone in array    
*                         variables PARC, PARM, STPAR ... minus 1               
*  INDSTC                 Location of first transversal disp. zone in array     
*                         variables PARC, PARM, STPAR, etc.                     
*  INFODC                 Location of first first order decay zone in array     
*                         variables PARC, PARM, STPAR ... minus 1               
*  INPORC                 Location of first porosity zone in array              
*                         variables PARC, PARM, STPAR ... minus 1               
*  INPRGC                 Location of first generic parameter zone in array     
*                         variables PARC, PARM, STPAR ... minus 1               
*  INQQPC                 Location of first prescribed flow zone in array       
*                         variables PARC, PARM, STPAR ... minus 1               
*  INSTGC                 Location of first storage coefficient zone in array   
*                         variables PARC, PARM, STPAR ... minus 1               
*  IZSUM                  Auxiliar variable to compute the indexes to locate 
*                         the parameters in arrays PARC,PARM, etc.
*
* HISTORY: AMS    3-1997    First coding
*          AMS    1-1998    Revision
*          AAR    2-2001    Revision and inclusion of observations related var.
*          AAR    4-2002    Revision and inclusion of geoestat. related var.
*          JHG	6-2003	  Change of IDIMDTRA for elementwise storage of DTRA.
*
*  I_TYPE_A               Type of A matrix.
*                         The possible types of matrices are:
*                             1 -->   nodewise Vector.
*                             2 -->   elementwise vector
*                             3 -->   derivative type matrix
*                             4 -->   Full matrix (elementwise).
*                             5 -->   Symmetric matrix (elementwise).
*                             6 -->   Symmetric matrix without diagonal (elementwise).
*                                     (It is supposed tha the sum of the terms
*                                     out of the diagonal is equal to the diagonal with
*                                     the opposite sign).
*                             7 -->   Symmetric banded matrix.
*                             8 -->   Non symmetric banded matrix.
*                             9 -->   Sparse matrix (as the one used by WatSolve).
*
********************************************************************************

       IMPLICIT REAL*8(A-H,O-Z)

       DIMENSION  INTRAC(6),INORPAR(NTYPAR),IOPT_GS(MXGRPZN,20)
     ;          ,IO_KG_GS(MXGRPZN,16),IPAR_DIR(NPARALG)

C------------------------- Matrix Type constants

      INTEGER*4,PARAMETER::
     &                     ITYP_NODE_VECTOR    = 1
     &                    ,ITYP_ELEM_VECTOR    = 2
     &                    ,ITYP_DERIV_MAT      = 3
     &                    ,ITYP_FULL_MAT       = 4
     &                    ,ITYP_ELEM_SYM_MAT   = 5
     &                    ,ITYP_SYM_NODIAG_MAT = 6
     &                    ,ITYP_SYM_BAND_MAT    = 7
     &                    ,ITYP_NOSYM_BAND_MAT = 8
     &                    ,ITYP_SPARSE_MAT     = 9

C------------------------- FIRST EXECUTABLE STATEMENT.
C------------------------- Definition of some variables used for dimensioning

      MAXNEOP = MAX (NUMEL,NUMNP,NDEVS,NPAR)
      NBAND1 = NBAND+1
      NBAND2 = 2*NBAND+1
      IDIMQ = IODIM*(IODIM+1)/2
      IDIMDQ = (LMXNDL*LMXNDL-LMXNDL)/2
      IDIMHESS = NPAR*(NPAR+1)/2
      IDIMCOV = (2*NUMTOBS - NBANDCOV + 1)*NBANDCOV/2
      NPBFLDIM = MAX(IOSMFL,NPBFL)

C------------------------- Initialization of DERH and DERC 3rd dimension
C------------------------- to avoid "subscript out of range errors"

      IDIMDERH=2
      IDIMDERC=2

      IDIMWORK = 0

C------------------------- Not used variables or not initilized
      IDIMDCROSS = 0


      IF (IOEQT.NE.1 .OR. IOFLLI.NE.0) THEN

          IDIMWORK = (NBAND1+2)*NUMNP*MAX(1,IOSMTP*NPBTP,NPBFLDIM)

      END IF

C------------------------- +2 space for auxiliar arrays
C------------------------- DERH1 and DERH2 used in DERQ_GEN module
C------------------------- Also computes IDIMDERH e IDIMDERC
      IF (IOINV.GT.0) THEN

          IDIMWORK = MAX(IDIMWORK,IDIMHESS,3*NPAR)

          IF (IOINV.EQ.1 .OR. IOINV.EQ.3) THEN

              IF (IODENS.GT.0 .OR. IOFLLI.GT.0) THEN

                  IDIMDERH = 3

               END IF !IODENS.GT.0 .OR. IOFLLI.GT.0

          END IF !IOINV.EQ.1 .OR. IOINV.EQ.3

          IF (IOINV.GE.2) THEN

              IF (IODENS.GT.0 .OR. IOTRLI.GT.0) THEN

                  IDIMDERC = 3
              END IF !IODENS.GT.0 .OR. IOTRLI.GT.0

          END IF !IOINV.GT.2
          

      END IF !IOINV.GT.0

C--------------- Initializes nr of columns of varios matrix types

      IDIMSYM = LMXNDL * (LMXNDL+1) / 2     !symmetric matrix
      IDIMFULL = LMXNDL * LMXNDL            !full matrix
      IDIMSYMND = (LMXNDL -1) * LMXNDL / 2  !symm. mat. without diagonal


C------------------------- Elementwise flow matrices


C------------------------- AFLU
C------------------------- Symm. matrix without diagonal elementwise.

      IDIMAFLU = IDIMSYMND
      ITYPAFLU = ITYP_SYM_NODIAG_MAT
                      
C------------------------- DFLU
C------------------------- Elementwise vector

      IDIMDFLU = LMXNDL
      ITYPDFLU = ITYP_ELEM_VECTOR
C------------------------- CFLU
C------------------------- Elementwise vector

      IDIMCFLU = LMXNDL*IODENS
      ITYPCFLU = ITYP_ELEM_VECTOR

C------------------------- BFLU
C------------------------- Nodewise vector

      ITYPBFLU = ITYP_NODE_VECTOR

C------------------------- DENSITY

      IF (IOVRWC.LT.2) THEN  !dimension of array density
          IDIMDENS = NUMEL
      ELSE
          IDIMDENS = NUMNP
      END IF !IOVRWC.EQ.1
c-parche-provisional-hasta-que-la-densidad-vaya-por-nudos
          IDIMDENS = NUMEL
c-fin-parche

C------------------------- Elementwise transport matrices

C------------------------- ATRA
C------------------------- Full matrix elementwise.

      IDIMATRA = IDIMFULL
      ITYPATRA = ITYP_FULL_MAT

C------------------------- DTRA
C------------------------- Elementwise vector

      IDIMDTRA = LMXNDL
      ITYPDTRA = ITYP_ELEM_VECTOR

C------------------------- BTRA
C------------------------- Nodewise vector

      ITYPBTRA = ITYP_NODE_VECTOR

C------------------------- System matrices
 
c------------------------- Check wether we always use coupled newton solver 
c------------------------- and wether we can use coupled newton solver

      IALW_CN = 0  
      ICAN_CN = 0

      IF (IODENS.EQ.1) THEN
          IF (IPAR_DIR(10).NE.0 .AND. IPAR_DIR(13).NE.0) THEN
              ICAN_CN=1
          END IF

          IF (IPAR_DIR(14).EQ.0 .AND. IPAR_DIR(9).EQ.0
     &       .AND.IPAR_DIR(11).EQ.0 .AND. IPAR_DIR(12).EQ.0) THEN
              IALW_CN=1
         END IF
      END IF !IODENS.EQ.1

C------------------------- AFLUDSC - ATRADSC

      IF (IOSPARSE.EQ.0) THEN  !Banded storage 

          IAFLUDSC_ROWS = NUMNP

          IATRADSC_ROWS = NUMNP
          IATRADSC_COLS = NBAND2

          ITYPTRADSC = ITYP_NOSYM_BAND_MAT
                        
          IF (IODENS.EQ.1 .OR. IOFLLI.EQ.1) THEN

              IAFLUDSC_COLS = NBAND2
              ITYPFLUDSC = ITYP_NOSYM_BAND_MAT

          ELSE

              IAFLUDSC_COLS = NBAND1
              ITYPFLUDSC = ITYP_SYM_BAND_MAT

         END IF !IODENS.EQ.1 .OR. IOFLLI.EQ.0


      ELSE !Sparse storage

          IAFLUDSC_ROWS = IPAR_DIR(19)
          IAFLUDSC_COLS = NUMNP
          IATRADSC_ROWS = IPAR_DIR(19)
          IATRADSC_COLS = NUMNP

          ITYPFLUDSC = ITYP_SPARSE_MAT
          ITYPTRADSC = ITYP_SPARSE_MAT

      END IF !IOSPARSE.EQ.0


C------------------------- ACOUPLEDDSC

      IA_COUPL_ROWS = 0
      IA_COUPL_COLS = 0

      IF (IODENS.EQ.1 .AND. ICAN_CN.EQ.1) THEN

          IF(IOSPARSE.EQ.0) THEN      !Banded

              IA_COUPL_ROWS = 2*NUMNP
              IA_COUPL_COLS = (NBAND2*2+1)
              ITYPCOUPLDSC = 8

C------------------------- If we use coupled newton, WORK might need
C------------------------- more space. Tthe needed space is
C------------------------- NUMNP * (NBAND +1) (for a generic matrix).
C------------------------- Coupled bandwidth is 2*NBAND +1
C------------------------- Coupled number of nodes is 2*NUMNP
C------------------------- Then the space needed in the coupled case is
C------------------------- 2*NUMNP * (2*NBAND+1 +1)
         
              I_SPACE_NEEDED = 2*(NUMNP)*(2*NBAND+1 + 1)
              IDIMWORK = MAX(IDIMWORK,I_SPACE_NEEDED)

          ELSE   !Sparse

              IA_COUPL_ROWS = IPAR_DIR(19)*2
              IA_COUPL_COLS = 2*NUMNP
              ITYPCOUPLDSC = 9

          END IF !IOSPARSE.EQ.0
      ELSE

          ITYPCOUPLDSC = 0

      END IF !IODENS.EQ.1 .AND. ICAN_CN.EQ.1

C------------------------- Derivatives matrices.	

      ITYPDERIV = 7

C------------------------- BIBI

      IDIMBB=1
      IF (IODIM.EQ.2) THEN
          IDIMBB=9
          IF (LMXNDL.EQ.4) IDIMBB=18
      ELSE IF (IODIM.EQ.3) THEN
          IDIMBB=36
          IF (LMXNDL.EQ.6) IDIMBB=90
      ENDIF


C------------------------- PARC, STPAR, etc

      IF (IOEQT.NE.2) THEN        ! Only flow or flow plus transport
          NZPAR = NZTRA*MAX(ISOT,IODIM)+NZSTG+NZARR+NZCHP+NZQQP+NZALF+
     ;          NZPRG+NZCLK
          IF (IOEQT.EQ.3) THEN
              NZPAR= NZPAR+2*NZDSP+NZPOR+NZDFM+NZCRD+NZFOD+NZCOE
          ELSE

C------------------------- Porosity has to be explicitly added when no 
C------------------------- transport is solved and unsaturated flow is required

              IF (IOFLSAT.NE.0) NZPAR=NZPAR+NZPOR      
          END IF !IOEQT.EQ.3
      ELSE

          NZPAR=2*NZDSP+NZPOR+NZDFM+NZCRD+NZFOD+NZCOE+NZPRG

      END IF !IOEQT.NE.2

C------------------------- Assigns NPAREL and NPARNP

      NPAREL = 11
      NPARNP = 10

C------------------------- Assigns index for zonal arrays: PARC, IVPAR, etc

      INTRAC(:)=0

      DO I=2,MAX (ISOT,IODIM)
          INTRAC(I)= INTRAC(I-1)+NZTRA
      END DO

      INSTGC = INTRAC( MAX (ISOT,IODIM))
      IZSUM = NZTRA
      IF (IOTRS.NE.0) INSTGC=INSTGC+IZSUM
      IF (NZSTG.NE.0) IZSUM=NZSTG
      INARRC=INSTGC
      IF (NZARR.NE.0) THEN
          INARRC=INARRC+IZSUM
          IZSUM=NZARR
      END IF
      INCHPC=INARRC
      IF (NZCHP.NE.0) THEN
          INCHPC=INCHPC+IZSUM
          IZSUM=NZCHP
      END IF
      INQQPC=INCHPC
      IF (NZQQP.NE.0) THEN
          INQQPC=INQQPC+IZSUM
          IZSUM=NZQQP
      END IF
      INALFC=INQQPC
      IF (NZALF.NE.0)  THEN
          INALFC=INALFC+IZSUM
          IZSUM=NZALF
      END IF
      INDSLC=INALFC
      IF (NZDSP.NE.0)  THEN
          INDSLC=INDSLC+IZSUM
          IZSUM=NZDSP
      END IF
      INDSTC=INDSLC
      IF (NZDSP.NE.0)  THEN
          INDSTC=INDSTC+IZSUM
          IZSUM=NZDSP
      END IF
      INDFMC=INDSTC
      IF (NZDFM.NE.0)  THEN
          INDFMC=INDFMC+IZSUM
          IZSUM=NZDFM
      END IF
      INPORC=INDFMC
      IF (NZPOR.NE.0)  THEN
          INPORC=INPORC+IZSUM
          IZSUM=NZPOR
      END IF
      INFODC=INPORC
      IF (NZFOD.NE.0)  THEN
          INFODC=INFODC+IZSUM
          IZSUM=NZFOD
      END IF
      INCRDC=INFODC
      IF (NZCRD.NE.0)  THEN
          INCRDC=INCRDC+IZSUM
          IZSUM=NZCRD
      END IF
      INCOEC=INCRDC
      IF (NZCOE.NE.0)  THEN
          INCOEC=INCOEC+IZSUM
          IZSUM=NZCOE
      END IF

      INPRGC=INCOEC
      IF (NZPRG.NE.0)  THEN
          INPRGC=INPRGC+IZSUM
          IZSUM = NZPRG
      END IF

      INCLK=INPRGC
      IF (NZCLK.NE.0)  THEN
          INCLK=INPRGC+IZSUM
      END IF

C_______________________ Geostat. inv. prob.
C_______________________ Initialization

       MXLINCMB = 0
       IDIMZONPP_GS = 0
       MXVGM_GS=0
       MXNST_GS=0
       MXROT_GS=0
       MXMEASPP_GS=0
       MXNPP_GS=0
       IDIMVAR_GS=0
       MXNZON_GS=0
       MXZONPP_GS = 0
       IDIMPOSDIS_GS=0
       MXNVAR_GS=0
       MXNPRIM_GS=0
       IDIMCROSS = 0 !Not used?
       IDIMCROSS_GS=0
       MXSAM_GS=0
       MXKRIG_GS=0
       MXCLOSE_GS=0
       IDIMIVARIO_GS=0
       IDIMBLPP_GS=0
       MXDISC_GS=0
       MXSB_GS=0
       MAX_1=0
       MAX_2=0
       MXSC_GS=1

       IF (IOINV_GS.NE.0) THEN                             ! Geostat. inv. prob.

         MXSC_GS=1
C_______________________ Calculates maximum dimensions among geological forms


         DO I=1,NGROUP_ZN
                                              ! Identifies useful variables
           NVAR_GR=IOPT_GS(I,3)    
           NMEAS_GR=IOPT_GS(I,4)   
           NPP_GR=IOPT_GS(I,5)     
           NZON_GR=IOPT_GS(I,7)    

           NMXP_GR=IO_KG_GS(I,3)   
           NMXS_GR=IO_KG_GS(I,4)
           MXNST_GR=IO_KG_GS(I,12)
           NEXDR_GR=IO_KG_GS(I,13)
           MXSC_GR=IO_KG_GS(I,16)
           IF (MXSC_GR.EQ.0) MXSC_GR=1
           NDISC_GR=IO_KG_GS(I,9)*IO_KG_GS(I,10)*IO_KG_GS(I,11)
           MXSB_GR=IO_KG_GS(I,6)*IO_KG_GS(I,7)*IO_KG_GS(I,8)

                                             ! Calculates maximum dimension
           
           IF (NVAR_GR*NVAR_GR.GE.MXVGM_GS) MXVGM_GS=NVAR_GR*NVAR_GR

           IF (MXNST_GR.GE.MXNST_GS) MXNST_GS=MXNST_GR

           IF (NMEAS_GR+NPP_GR.GE.MXMEASPP_GS) 
     ;                                   MXMEASPP_GS=NMEAS_GR+NPP_GR

           IF (NEXDR_GR.GE.MXNEXDR_GS) MXNEXDR_GS=NEXDR_GR

           IF (NPP_GR.GE.MXNPP_GS) MXNPP_GS=NPP_GR

           IF (NVAR_GR.GE.MXNVAR_GS) MXNVAR_GS=NVAR_GR

           IF (NZON_GR.GE.MXNZON_GS) MXNZON_GS=NZON_GR

           IF (NDISC_GR.GE.MXDISC_GS) MXDISC_GS=NDISC_GR

           IF (NMXP_GR.GE.MXNPRIM_GS) MXNPRIM_GS=NMXP_GR

           
           IF (NMXP_GR+NMXS_GR.GE.MXSAM_GS) MXSAM_GS=NMXP_GR+NMXS_GR

           IF (NMXP_GR+NMXS_GR+11.GE.MAX_1) MAX_1=NMXP_GR+NMXS_GR+11
           IF ( (NMXP_GR+NMXS_GR)*NVAR_GR+NVAR_GR.GE.MAX_2) 
     ;                       MAX_2=(NMXP_GR+NMXS_GR)*NVAR_GR+NVAR_GR
           MAX_MAX=MAX0(MAX_1,MAX_2)
           MAX_1=0
           MAX_2=0
           MXKRIG_GS=MAX0(MXKRIG_GS,MAX_MAX)

           IF (NMEAS_GR+NPP_GR+NZON_GR.GE.MXCLOSE_GS) 
     ;                            MXCLOSE_GS=NMEAS_GR+NPP_GR+NZON_GR

           IF (NVAR_GR*NVAR_GR*MXNST_GR.GE.IDIMIVARIO_GS)
     ;                    IDIMIVARIO_GS=NVAR_GR*NVAR_GR*MXNST_GR


           IF (NMXP_GR*NZON_GR.GE.IDIMZONPP_GS) 
     ;                    IDIMZONPP_GS=NMXP_GR*NZON_GR


           IF (MXSB_GR.GE.MXSB_GS) MXSB_GS=MXSB_GR

           IF (MXSC_GR.GE.MXSC_GS) MXSC_GS=MXSC_GR


         END DO

         MXROT_GS=MXNST_GS*MXVGM_GS+1
         IDIMVAR_GS=MXNVAR_GS+4
         MXZONPP_GS=MAX0(MXNPP_GS,MXNZON_GS)
         IDIMWORK=MAX0(IDIMWORK,MXKRIG_GS*(MXKRIG_GS)/2)
         IDIMWORK=MAX0(IDIMWORK,2*IDIMHESS)
         IDIMDATASC_GS=MXZONPP_GS+MXMEASPP_GS
         IDIMCROSS_GS=IDIMDATASC_GS
       END IF

       IDIMWGT=MAX0(MXLINCMB,MXNPP_GS)
       IF (IDIMWGT.EQ.0) IDIMWGT=1   ! Necessary for dimensioning

C------------------------- Writes all integer variables used for dimensioning 
C------------------------- or indexing

      IF (IOPART.NE.0) THEN
          WRITE(MAINF,2999)
          WRITE(MAINF,3000) 
     ;     'IDIMAFLU',IDIMAFLU,         'IDIMBB',IDIMBB,
     ;     'IDIMCOV',IDIMCOV,           'IDIMCROSS',IDIMCROSS,
     ;     'IDIMDENS',IDIMDENS,
     ;     'IDIMDFLU',IDIMDFLU,         'IDIMDQ',IDIMDQ,
     ;     'IDIMDTRA',IDIMDTRA,         'IDIMFNT',IDIMFNT,
     ;     'IDIMHESS',IDIMHESS,         'IDIMQ',IDIMQ,
     ;     'IOCNSF',IOCNSF,             'IOCNST',IOCNST,
     ;     'IODIM',IODIM,
     ;     'IOEQT',IOEQT,               'IOFLLI',IOFLLI,
     ;     'IOFLSAT',IOFLSAT,           'IOINV',IOINV,
     ;     'IORTS',IORTS,               'IOTRLI',IOTRLI,
     ;     'IOTRS',IOTRS,               'IPARTRA',IPARTRA,
     ;     'ISOT',ISOT,                 'LMXNDL',LMXNDL,
     ;     'MAINF',MAINF,               'MAXNEOP',MAXNEOP,
     ;     'NBAND',NBAND,               'NBAND1',NBAND1,
     ;     'NBAND2',NBAND2,             'NBANDCOV',NBANDCOV,
! Geoestadisticos.............
     ;     'MXCLOSE_GS',MXCLOSE_GS,     'MXDISC_GS',MXDISC_GS,
     ;     'MXGRPZN',MXGRPZN,           'MXKRIG_GS',MXKRIG_GS,
     ;     'MXLINCMB',MXLINCMB,         'MXMEASPP_GS',MXMEASPP_GS,
     ;     'MXNPP_GS',MXNPP_GS,         'MXNPRIM_GS',MXNPRIM_GS,
     ;     'MXNST_GS',MXNST_GS,         'MXNVAR_GS',MXNVAR_GS,
     ;     'MXNZON_GS',MXNZON_GS,       'MXROT_GS',MXROT_GS,
     ;     'MXSAM_GS',MXSAM_GS,         'MXSB_GS',MXSB_GS,
     ;     'MXSC_GS',MXSC_GS,           'MXVGM_GS',MXVGM_GS,
     ;     'MXZONPP_GS',MXZONPP_GS,     'NGROUP_ZN',NGROUP_ZN,
     ;     'IDIMCROSS_GS',IDIMCROSS_GS, 'IDIMIVARIO_GS',IDIMIVARIO_GS,
     ;     'IDIMVAR_GS',IDIMVAR_GS,     'IDIMWGT',IDIMWGT, 
     ;     'IDIMZONPP_GS',IDIMZONPP_GS, 'IOINV_GS',IOINV_GS,
     ;     'NDEVS',NDEVS,               'NFNT',NFNT,
     ;     'NINT',NINT,                 'NPAR',NPAR,
     ;     'NPAREL',NPAREL,             'NPARNP',NPARNP,
     ;     'NTYPAR',NTYPAR,
     ;     'NUMEL',NUMEL,               'NUMNP',NUMNP,
     ;     'NUMTIT',NUMTIT,             'NUMTNOD',NUMTNOD,
     ;     'NUMTOBS',NUMTOBS,           'NZALF',NZALF,
     ;     'NZARR',NZARR,               'NZCHP',NZCHP,
     ;     'NZCOE',NZCOE,               'NZCRD',NZCRD,
     ;     'NZDFM',NZDFM,               'NZDMT',NZDMT,
     ;     'NZDSP',NZDSP,               'NZFOD',NZFOD,
     ;     'NZPAR',NZPAR,               'NZPOR',NZPOR,
     ;     'NZPRG',NZPRG,               'NZQQP',NZQQP,
     ;     'NZSTG',NZSTG,               'NZTRA',NZTRA,
     &     'NZCLK',NZCLK,                'INCLK',INCLK,
     ;     'INALFC',INALFC,             'INARRC',INARRC,
     ;     'INCHPC',INCHPC,             'INCOEC',INCOEC,
     ;     'INCRDC',INCRDC,             'INDFMC',INDFMC,
     ;     'INDSLC',INDSLC,             'INDSTC',INDSTC,
     ;     'INFODC',INFODC,             'INPORC',INPORC,
     ;     'INPRGC',INPRGC,             'INQQPC',INQQPC,
     ;     'INSTGC',INSTGC,             'INTRAC(1)',INTRAC(1),
     ;     'INTRAC(2)',INTRAC(2),       'INTRAC(3)',INTRAC(3),
     ;     'INTRAC(4)',INTRAC(4),       'INTRAC(5)',INTRAC(5),
     ;     'INTRAC(6)',INTRAC(6),
     ;     'IDIMCFLU',IDIMCFLU,         'IDIMATRA',IDIMATRA,      
     ;     'ICAN_CN',ICAN_CN,           'IALW_CN',IALW_CN,
     ;     'IOSPARSE',IOSPARSE,         'IAFLUDSC_ROWS',IAFLUDSC_ROWS,           
     ;     'IAFLUDSC_COLS',IAFLUDSC_COLS,'IACOUPL_ROWS',IA_COUPL_ROWS,
     ;     'IACOUPL_COLS',IA_COUPL_COLS, 'IATRADSC_ROWS',IATRADSC_ROWS,
     ;     'IATRADSC_COLS',IATRADSC_COLS


      END IF !IOPART.NE.0

 2999  FORMAT(///,21X,'INDEXING OR DIMENSIONING VARIABLES ',/,
     ;            21X,'======== == ============ =========',//)

 3000  FORMAT(2(A13,I10))

C--------------------------------------Write matrix types
3010   FORMAT(///,21X,'MATRIX TYPES',/,
     ;            21X,'====== =====',//)

3020  FORMAT
     ;('MATRIX A OF FLOW......................=',I5,/,
     ;'MATRIX C OF FLOW......................=',I5,/,
     ;'MATRIX D OF FLOW......................=',I5,/,
     ;'VECTOR B OF FLOW......................=',I5,/,
     ;'MATRIX A OF TRANSPORT.................=',I5,/,
     ;'MATRIX D OF TRANSPORT.................=',I5,/,
     ;'VECTOR B OF TRANSPORT.................=',I5,/,
     ;'SYSTEM MATRIX FLOW....................=',I5,/,
     ;'SYSTEM MATRIX TRANSPORT...............=',I5,/,
     ;'SYSTEM MATRIX COUPLED.................=',I5,/,
     ;'ALL EXISTING DERIVATIVES MATRICES.....=',I5,/)



      IF (IOPART.NE.0) THEN
          WRITE (MAINF,3010)
          WRITE(MAINF,3020) ITYPAFLU, ITYPCFLU, ITYPDFLU
     ;    , ITYPBFLU, ITYPATRA, ITYPDTRA, ITYPBTRA 
     ;    ,ITYPFLUDSC, ITYPTRADSC, ITYPCOUPLDSC, ITYPDERIV
      END IF

C------------------------- Groups indexes in array INORPAR

       INORPAR(1) =INTRAC(1)
       INORPAR(2) =INTRAC(2)
       INORPAR(3) =INTRAC(3)
       INORPAR(4) =INTRAC(4)
       INORPAR(5) =INTRAC(5)
       INORPAR(6) =INTRAC(6)
       INORPAR(7) =INSTGC
       INORPAR(8) =INARRC
       INORPAR(9) =INCHPC
       INORPAR(10)=INQQPC
       INORPAR(11)=INALFC
       INORPAR(12)=INDSLC
       INORPAR(13)=INDSTC
       INORPAR(14)=INDFMC
       INORPAR(15)=INPORC
       INORPAR(16)=INFODC
       INORPAR(17)=INCRDC
       INORPAR(18)=INCOEC
       INORPAR(19)=INPRGC
*	 INORPAR(20)=INAGE
       INORPAR(21)=INCLK

C------------------------- Gravity vector

      IF (IODENS.EQ.1 .OR. IOFLLI.EQ.1) THEN
          IDIMGRAVEL=NUMEL
      ELSE
          IDIMGRAVEL=0
      END IF

      END SUBROUTINE INIT_INTEG 
