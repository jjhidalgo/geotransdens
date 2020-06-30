      SUBROUTINE PRINCIPAL
     ;(!--------------------------------------------------dimensioning variables 
     ;IAFLUDSC_COLS         ,IAFLUDSC_ROWS          ,IATRADSC_COLS  
     ;,IATRADSC_ROWS       ,IA_COUPLED_DSC_COLS    ,IA_COUPLED_DSC_ROWS
     ;,IDIMAFLU            ,IDIMDFLU               ,IDIMCFLU
     ;,IDIMATRA
     &,IDIMDTRA            ,IDIMDQ     ,IDIMBB                             
     ;,IDIMCOV             ,IDIMFNT                ,IDIMGRAVEL ,IDIMHESS               
     ;,IDIMQ               ,IDIMWORK   ,IDIMDERC    ,IDIMDERH                
     ;,NBAND               ,NBAND1                 ,NBANDCOV                 
     ;,NDEVS               ,NFLAGS                 ,NUMTOBS,NBLCVP
     ;!------------------------------------------some new dimensioning variables
     ;,MAXNB    ,maxnbf   ,NPBMX        ,MAXPG      ,NTYPEL
     ;,NZPRG    ,NUMTIT   ,NUMTNOD 
     ;,NZTRA
     ;,IDIMDENS ,NMAXF    ,NMAXT
     ;!-------------------------------------------------------------io-variables
     ;,IOCNSF   ,IOCNST   ,IOCONSRC  ,IOCRITRAP    ,IODENS_INI ,IODIM    
     ;,IOEQT    ,IOFLLI   ,IOFLSAT   ,IOFMLF       ,IOFMLT
     ;,IOINV    ,IOPRHED  ,IOPINITC ,IOPINITH     ,IORTS     ,IOSMFL
     &,IOSMTP   ,IOSPARSE ,IODIRECT
     ;,IOTRLI   ,IOTRS    ,IOVAR        ,ITPTVAR
     ;
     ;!-------------------------------------------------------------------arrays 
     ;,FOBJ_WGT  ,IFLAGS    ,INORPAR   ,IOLG_PAR  ,IOPTS     ,IOWRITE  
     ;,IPAR_DIR  ,IPAR_INV  ,LINMET    
     &,NZONE_PAR ,PAR_DIR   ,PAR_INV   ,PAR_WGT
     ;!--------------------------------------------------------general variables	 
     ;,FILENAME ,IERROR   ,ISOT      ,LMXNDL   ,MAINF    ,NFNL     
     ;,NFNT     ,NINT     ,NOPTS     ,NPAR     ,NPARALG  ,NPAREL   
     ;,NPARF    ,NPARFPRG ,NPARNP    ,NPBFL    ,NPBTP    
     ;,NPPEL    ,NPPNP    ,NSTAT     ,NTDMT    ,NTYPAR   ,NUMEL    
     ;,NUMNP    ,NWRITE   ,NZPAR    
     ;!-------------------------------------------geostatistical problem scalars
     ;,MXROT_GS      ,MXMEASPP_GS
     ;,MXNPP_GS    ,IDIMVAR_GS  ,MXNZON_GS     ,MXDISC_GS   
     ;,MXNVAR_GS   ,MXNPRIM_GS  ,IDIMCROSS_GS  ,MXSAM_GS
     ;,MXKRIG_GS   ,MXCLOSE_GS  ,IDIMIVARIO_GS ,IDIMZONPP_GS 
     ;,MXSB_GS     ,IDIMWGT     ,IOINV_GS
     ;,MXGRPZN     ,NGROUP_ZN     ,IOPT_GS
     ;,IO_KG_GS    ,MXSC_GS     ,IDIMDATASC_GS
     ;!-------------------------------------------------------------Matrix types 
     ;,ITYPAFLU  ,ITYPCFLU     ,ITYPDFLU
     &,ITYPATRA  ,ITYPDTRA     ,ITYPFLUDSC
     &,ITYPTRADSC,ITYPACOUPLDSC
     ;!--------------------------------------------------------------real arrays
     ;,PARC     ,PARM      ,STPAR   ,CFPAREL
     ;,CFPARNP  ,PARNP     ,PAREL   ,COVPAR 
     ;,SOURCE   ,GRAV    ,GRAVEL
     ;,AFLU     ,BFLU      ,HCALAN
     ;,DFLU     ,ALFA      ,HAUX1   ,HAUX2
     ;,CFLU     ,AFLUDSC   ,AFLUDSCF
     ;,DERH     ,BM_ND_FL  ,BM_ZN_FL      
     ;,solution              
     ;!------------------------------------------ Inverse problem general arrays
     ;,GRAD     ,HESSAUX   ,HESS    ,DLT_PAR
     ;,PARAUX    ,COVINV  ,VJAC
     ;!---------------------------Coordinates and finite element integral arrays
     ;
     ;,COORD    ,BIBI      ,VOLNOD   ,AREA
     ;,GRDFF
     ;!-------------------------------------Auxiliar arrays and Real time arrays
     ;,DAT_VD   ,WORK      ,FNT      ,TIME     
     ;,CAUDAL   ,CONCFLOW
     ;!------------------------------------ Darcy velocity and related variables
     ;,VD       ,QXYZ      ,XNORVD 
     ;,dqdflu   ,dqdtra    ,dvdh
     ;,dvdc
     ;!-----------------------------------------------Consistent velocity arrays
     ;,POINTWEIGHT ,GRADLOC, BUOYANCY,DBUOYANCY
     ;,GP_COORD
     ;!---------------------------------------------------------transport arrays
     ;,WATVOL   ,ACTH               
     ;,DERC     ,CAUX1     ,CAUX2     ,DVDP        
     ;,ATRA     ,ATRADSC   ,ATRADSCF  ,BTRA       
     ;,CCALAN   ,DTRA      ,BM_ZN_TT  ,BM_ND_TT  ,DWDH
     ;
     ;!------------------------------------------------------ Observation arrays
     ;,VOBS     ,VOBSC     ,DVOBS   ,TOBS	
     ;,TIT      ,BUDAT     ,EXTNBU  ,WTOBSBU                       
     ;,WTOBSU   ,WTOBSN    ,WTOBST                          
     ;!---------------------------------------------------------- FLOW nonlinear
     ;,HBASE    ,HINI      ,HPREV1    ,DFLUDFLU ,DBFLUDFLU
     ;,HPREV2   ,DNODALRH  ,DPARELDH
     ;!------------------------------------------------------Transport nonlinear
     ;,CBASE    ,CPREV1     ,CPREV2     
     ;,CCALIT   ,DTRADTRA
     ;!-----------------------------------------------Coupled flow and transport
     ;,A_COUPL_DSC ,A_COUPL_DSCF,DFLUDTRA ,DTRADFLU ,DPARELDC
     ;,BCOUPLED    ,DBFLUDTRA,DER_VISC
     ;!------------------------- Common non linear variables (flow or transport
     ;,HCALIT   ,DTMXDS      ,PARACD
     ;
     ;!------------------------------------------------density dependency arrays
     ;,DENSITY  ,VISCOSITY   ,DELTAITER
     ;
     ;!------------------------------- Geostatistical inverse problem variables
     ;,PARZ         ,WGT_PAR   ,DERIV  ,IPOS
     ;,IPNT_PAR     ,VARIO_GS          ,ROTMAT_GS
     ;,POSMEAS_GS   ,VMEAS_GS          ,POSDIS_GS
     ;,POSDISAUX_GS ,POSZN_GS          ,DIVZN_GS
     ;,VSTATS_GS    ,TRIM_GS           ,SEARCH_GS
     ;,ESTKRIG_GS   ,ZNWGT_GS          ,CROSSCOV_GS
     ;,DATASC_GS         ,SUPBL_GS
     ;,KRISYS_GS    ,KRISOL_GS         ,CLOSESAM_GS
     ;,KRIGAUX_GS   ,EXDRZN_GS         ,COORDGR_GS
     ;,IVARIO_GS    ,IPOLDRIFT_GS
     ;,IZN_PP_GS    ,IZN_NPP_GS        ,ICROSSCOV_GS
     ;,ISUPBL_GS         ,NUMSB_GS
     ;,MXZONPP_GS   ,ICHECK_GS         ,LDIM_GS
     ;,WGT_UNK      ,PARGOOD           ,PARC_GS
     ;,IFLAG_SIMUL  ,COVPAR_GR
     ;                    ,DTPREVINV
     ;!--------------------------------------------- Statistical analysis arrays
     ;,DEVICESTAT ,EIGENVEC ,RESID     ,RESIDPAR
     ;
     ;!----------------------------------------------------------- String arrays
     ;,DEVNAME
     ;
     ;!================================================================
     ;!--------------- Integer arrays
     ;,IBCOD   ,ISOZ     ,NFTPAR     ,NFNLPAR          
     ;,IVPAR   ,IDMBLCVP ,LXPAREL    ,IXPARNP
     ;
     ;!--------------------------------------------------- Discretization arrays
     ;,KXX      ,LTYPE    ,LNNDEL     ,LDIM
     ;,IBTCO    ,INDPAR
     ;!----------------------------------------------------Sparse storage arrays
     ;,IAD_S    ,IADD_S   ,IADN_S     ,IAD_D
     ;,IADD_D   ,IADN_D   
     ;,iafd_s   ,iafdd_s  ,iafdn_s    ,iafd_d
     ;,iafdd_d  ,iafdn_d   
     ;
     ;!---------------- Integer arrays related to measurements and to dimensions
     ;!---------------- NUMTNOD ,NUMTOBS and NUMTIT
     ;,IODEVICE ,INDEXNOD  ,IOBUTYP   ,IOCALBU
     ;,IOUTYP   ,NBUW      ,NOBUF     ,NOOBSIT
     ;,IOTINT   ,MEASTYP 
     ;
     ;!--------------- Time arrays & Auxiliary array & Non linear integer arrays
     ;,KINT     ,ISOLEQ   , NFNLTIP
     ;,NFNLPRG  ,IPARTNER ,LCOORD
     ; 
     ;!-------------------------------------------- Statistical analysis arrays
     ;
     ;,OBSCLASS,ITYPEPAR)


       IMPLICIT NONE

       INTEGER*4 
     ;IAFLUDSC_COLS         ,IAFLUDSC_ROWS          ,IATRADSC_COLS  
     ;,IATRADSC_ROWS       ,IA_COUPLED_DSC_COLS    ,IA_COUPLED_DSC_ROWS
     ;,IDIMAFLU            ,IDIMDFLU               ,IDIMCFLU
     ;,IDIMATRA            ,IDIMDTRA            ,IDIMDQ     ,IDIMBB                             
     ;,IDIMCOV             ,IDIMFNT                ,IDIMGRAVEL ,IDIMHESS               
     ;,IDIMQ               ,IDIMWORK   ,IDIMDERC    ,IDIMDERH                
     ;,NBAND               ,NBAND1                 ,NBANDCOV                 
     ;,NDEVS               ,NFLAGS                 ,NUMTOBS,NBLCVP
     ;,MAXNB    ,MAXNBF   ,NPBMX        ,MAXPG      ,NTYPEL
     ;,NZPRG    ,NUMTIT   ,NUMTNOD 
     ;,NZTRA
     ;,IDIMDENS ,NMAXF    ,NMAXT
     ;,IOCNSF   ,IOCNST   ,IOCONSRC  ,IOCRITRAP    ,IODENS_INI ,IODIM    
     ;,IOEQT    ,IOFLLI   ,IOFLSAT   ,IOFMLF       ,IOFMLT
     ;,IOINV    ,IOPRHED  ,IOPINITC ,IOPINITH     ,IORTS     ,IOSMFL
     ;,IOSMTP   ,IOSPARSE ,IODIRECT
     ;,IOTRLI   ,IOTRS    ,IOVAR        ,ITPTVAR
     ;

C -------------------------------------------------------------------arrays 

     ;,IFLAGS    ,IOLG_PAR  ,IOPTS     ,IOWRITE  
     ;,IPAR_DIR  ,IPAR_INV  ,LINMET    
     &,NZONE_PAR
     ;!--------------------------------------------------------general variables	 
     ;,IERROR   ,ISOT     ,LMXNDL    ,MAINF    ,NFNL     
     ;,NFNT     ,NINT     ,NOPTS     ,NPAR     ,NPARALG  ,NPAREL   
     ;,NPARF    ,NPARFPRG ,NPARNP    ,NPBFL    ,NPBTP    
     ;,NPPEL    ,NPPNP    ,NSTAT     ,NTDMT    ,NTYPAR   ,NUMEL    
     ;,NUMNP    ,NWRITE   ,NZPAR    
     ;!-------------------------------------------geostatistical problem scalars
     ;,MXROT_GS      ,MXMEASPP_GS
     ;,MXNPP_GS    ,IDIMVAR_GS  ,MXNZON_GS     ,MXDISC_GS   
     ;,MXNVAR_GS   ,MXNPRIM_GS  ,IDIMCROSS_GS  ,MXSAM_GS
     ;,MXKRIG_GS   ,MXCLOSE_GS  ,IDIMIVARIO_GS ,IDIMZONPP_GS 
     ;,MXSB_GS     ,IDIMWGT     ,IOINV_GS
     ;,MXGRPZN     ,NGROUP_ZN     ,IOPT_GS
     ;,IO_KG_GS    ,MXSC_GS     ,IDIMDATASC_GS ,IPNT_PAR     
     ;!-------------------------------------------------------------Matrix types 
     ;,ITYPAFLU  ,ITYPCFLU     ,ITYPDFLU
     &,ITYPATRA  ,ITYPDTRA     ,ITYPFLUDSC
     &,ITYPTRADSC,ITYPACOUPLDSC ,IPOS 
     ;,NPARDET     ,ISIM_GS     ,NFL_SIM   
     ;,NTP_SIM     ,IOVAR_ORIG  ,IORDCH

       REAL*8 

     ; PARC     ,PARM      ,STPAR     ,CFPAREL
     ;,CFPARNP  ,PARNP     ,PAREL     ,COVPAR 
     ;,SOURCE   ,GRAV      ,GRAVEL
     ;,AFLU     ,BFLU      ,HCALAN
     ;,DFLU     ,ALFA      ,HAUX1     ,HAUX2
     ;,CFLU     ,AFLUDSC   ,AFLUDSCF
     ;,DERH     ,BM_ND_FL  ,BM_ZN_FL  ,CNST     
     ;,SOLUTION ,PAR_DIR   ,PAR_INV   ,PAR_WGT
     ;!------------------------------------------ Inverse problem general arrays
     ;,GRAD     ,HESSAUX   ,HESS    ,DLT_PAR
     ;,PARAUX   ,COVINV    ,VJAC
     ;!---------------------------Coordinates and finite element integral arrays
     ;
     ;,COORD    ,BIBI      ,VOLNOD   ,AREA
     ;,GRDFF
     ;!-------------------------------------Auxiliar arrays and Real time arrays
     ;,DAT_VD   ,WORK      ,FNT      ,TIME     
     ;,CAUDAL   ,CONCFLOW
     ;!------------------------------------ Darcy velocity and related variables
     ;,VD       ,QXYZ      ,XNORVD 
     ;,DQDFLU   ,DQDTRA    ,DVDH
     ;,DVDC
     ;!-----------------------------------------------Consistent velocity arrays
     ;,POINTWEIGHT ,GRADLOC, BUOYANCY,DBUOYANCY
     ;,GP_COORD
     ;!---------------------------------------------------------transport arrays
     ;,WATVOL   ,ACTH               
     ;,DERC     ,CAUX1     ,CAUX2     ,DVDP        
     ;,ATRA     ,ATRADSC   ,ATRADSCF  ,BTRA       
     ;,CCALAN   ,DTRA      ,BM_ZN_TT  ,BM_ND_TT  ,DWDH
     ;
     ;!------------------------------------------------------ Observation arrays
     ;,VOBS     ,VOBSC     ,DVOBS   ,TOBS	
     ;,TIT      ,BUDAT     ,EXTNBU  ,WTOBSBU                       
     ;,WTOBSU   ,WTOBSN    ,WTOBST                          
     ;!---------------------------------------------------------- FLOW nonlinear
     ;,HBASE    ,HINI      ,HPREV1    ,DFLUDFLU ,DBFLUDFLU
     ;,HPREV2   ,DNODALRH  ,DPARELDH
     ;!------------------------------------------------------Transport nonlinear
     ;,CBASE    ,CPREV1     ,CPREV2     
     ;,CCALIT   ,DTRADTRA
     ;!-----------------------------------------------Coupled flow and transport
     ;,A_COUPL_DSC ,A_COUPL_DSCF,DFLUDTRA ,DTRADFLU ,DPARELDC
     ;,BCOUPLED    ,DBFLUDTRA,DER_VISC
     ;!------------------------- Common non linear variables (flow or transport
     ;,HCALIT   ,DTMXDS      ,PARACD
     ;
     ;!------------------------------------------------density dependency arrays
     ;,DENSITY  ,VISCOSITY   ,DELTAITER
     ;
     ;!------------------------------- Geostatistical inverse problem variables
     ;,PARZ         ,WGT_PAR           ,DERIV 
     ;,VARIO_GS     ,ROTMAT_GS
     ;,POSMEAS_GS   ,VMEAS_GS          ,POSDIS_GS
     ;,POSDISAUX_GS ,POSZN_GS          ,DIVZN_GS
     ;,VSTATS_GS    ,TRIM_GS           ,SEARCH_GS
     ;,ESTKRIG_GS   ,ZNWGT_GS          ,CROSSCOV_GS
     ;,DATASC_GS         ,SUPBL_GS
     ;,KRISYS_GS    ,KRISOL_GS         ,CLOSESAM_GS
     ;,KRIGAUX_GS   ,EXDRZN_GS         ,COORDGR_GS
     ;,WGT_UNK      ,PARGOOD           ,PARC_GS
     ;,COVPAR_GR
     ;                    ,DTPREVINV   ,FOBJ_WGT 
     ;!--------------------------------------------- Statistical analysis arrays
     ;,DEVICESTAT ,EIGENVEC ,RESID     ,RESIDPAR

       INTEGER*4

     ; IBCOD   ,ISOZ     ,NFTPAR     ,NFNLPAR          
     ;,IVPAR   ,IDMBLCVP ,LXPAREL    ,IXPARNP
     ;
     ;!--------------------------------------------------- Discretization arrays
     ;,KXX      ,LTYPE    ,LNNDEL     ,LDIM
     ;,IBTCO    ,INDPAR
     ;!----------------------------------------------------Sparse storage arrays
     ;,IAD_S    ,IADD_S   ,IADN_S     ,IAD_D
     ;,IADD_D   ,IADN_D   
     ;,IAFD_S   ,IAFDD_S  ,IAFDN_S    ,IAFD_D
     ;,IAFDD_D  ,IAFDN_D   ,INORPAR(NTYPAR)   
     ;
     ;!---------------- Integer arrays related to measurements and to dimensions
     ;!---------------- NUMTNOD ,NUMTOBS and NUMTIT
     ;,IODEVICE ,INDEXNOD  ,IOBUTYP   ,IOCALBU
     ;,IOUTYP   ,NBUW      ,NOBUF     ,NOOBSIT
     ;,IOTINT   ,MEASTYP 
     ;
     ;!--------------- Time arrays & Auxiliary array & Non linear integer arrays
     ;,KINT     ,ISOLEQ   ,NFNLTIP
     ;,NFNLPRG  ,IPARTNER ,LCOORD
     ; 
     ;!------------------------- Geoestatistical inverse problem integer arrays
     ;,IVARIO_GS    ,IPOLDRIFT_GS
     ;,IZN_PP_GS    ,IZN_NPP_GS        ,ICROSSCOV_GS
     ;,ISUPBL_GS         ,NUMSB_GS
     ;,MXZONPP_GS   ,ICHECK_GS         ,LDIM_GS ,IFLAG_SIMUL  
     ;
     ;!-------------------------------------------- Statistical analysis arrays
     ;
     ;,OBSCLASS,ITYPEPAR

      CHARACTER FILENAME(20)*20,DEVNAME(NDEVS)*10

	Character*10,Allocatable::ParName(:)

C     INTERNAL VARIABLES: SCALARS
      INTEGER*4 NROW,iprocess
     &, INTRA,INSTG,INARR,INARRT,INDSP,INDFM ,INPOR ,INFOD,INCRD
     &, INCOE,INCHP ,INCHPT,INQQP ,INQQPT,INALF,INALFT,INCON,INCONT
     &, INDMT,MAXITER_ORIG
     &, i
     &!NO SE UTILIZA
     & ,IOPINVDT
 
      REAL*8 BETAC,CREF,DENSREF,TEMPREF,VISCREF
     &      ,DRELMX_ORIG,DABSMX_ORIG
     &      ,XMARQ_ORIG, VAR_REF,WSPECHEAT
     &      ,WTHERMCON,PRESSURE

       DIMENSION 
     ;      NZONE_PAR(NTYPAR),IPAR_DIR(NPARALG)
     ;     ,LXPAREL(NUMEL,NPAREL)
     ;     ,IOWRITE(NWRITE),IFLAGS(NFLAGS),COORD(NUMNP,3)
     ;     ,IXPARNP(NUMNP,NPARNP),AFLU(NUMNP,IDIMAFLU)
     ;     ,ATRA(NUMNP,2*NBAND+1),PAR_DIR(NPARALG),CNST(6,6,6)
     ;     ,PARZ(NZPAR),IOPTS(NOPTS),HCALIT(NUMNP,NPBFL)

C------------------------- First executable statement

      IPROCESS = 0

      NPBMX = MAX(NPBFL,NPBTP)

      IORDCH = IOPTS(30)

      If(NPAR.GT.0) Then
          Allocate(ParName(NPAR))
      Else
	  Allocate(ParName(1)) !por si...
       End If
       ParName(:)=''

c-parche
C------------------------- This variable is not in the input data.
      IOPINVDT = 0
c-fin-parche
C------------------------- Reads all input data

       CALL ENTDAT 
     ;(PAR_DIR(31)      ,PAR_DIR(32)  ,IDIMCOV    ,IDIMFNT
     ;,IDIMQ            ,IERROR       ,IOCNSF
     ;,IOCNST           ,IODIM        ,IOEQT      ,IOFLLI   
     ;,IOFLSAT          ,IOFMLF       ,IOFMLT     ,IOINV
     ;,IORTS            ,IOTRLI     ,IOTRS
     ;                  ,ISOT         ,LMXNDL     ,MAINF    
     ;,NBAND            ,NDEVS        ,NFNL       ,NFNT     
     ;,NINT             ,NPAR         ,NPAREL     ,NPARF    
     ;,NPARFPRG         ,NPARNP                   ,NPBFL
     ;,NPBMX            ,NPBTP        ,NROW       ,NTDMT    
     ;,NTYPAR           ,NUMEL        ,NUMNP      ,NUMTIT   
     ;,NUMTNOD          ,NUMTOBS      ,NWRITE     ,NZONE_PAR(16)
     ;,NZPAR            ,PAR_DIR(29),PAR_DIR(30)
     ;,ACTH             ,AREA         ,BTRA       ,BUDAT
     ;,CAUDAL           ,CCALIT       ,CFPAREL    
     ;,CFPARNP          ,COORD        ,COVINV     ,COVPAR     
     ;,DTMXDS           ,EXTNBU       ,FNT        ,GRAV       
     ;,HCALIT           ,HCALAN       ,IBCOD      ,IBTCO        
     ;,IFLAGS           ,INDEXNOD     ,INDPAR     ,INORPAR      
     ;,IOBUTYP          ,IOCALBU      ,IODEVICE   ,IOLG_PAR     
     ;,IOTINT           ,MEASTYP      ,IOUTYP     ,IOWRITE ,ISOLEQ       
     ;,ISOZ             ,IVPAR        ,IXPARNP    ,KINT         
     ;,KXX              ,LDIM         ,LNNDEL     ,LTYPE        
     ;,LXPAREL          ,NBUW         ,NFNLPAR    ,NFNLPRG      
     ;,NFNLTIP          ,NFTPAR       ,NOBUF      ,NOOBSIT      
     ;,NZONE_PAR        ,PARACD       ,PARZ       ,PARM         
     ;,QXYZ             ,STPAR        ,TIME       ,TIT          
     ;,TOBS             ,VD           ,VOBS       ,WTOBSBU      
     ;,WTOBSN           ,WTOBST       ,WTOBSU     ,XNORVD       
     ;,DEVNAME          ,FILENAME
     ;,IDMBLCVP         ,NBLCVP       ,NFLAGS
     ;,IO_KG_GS         ,IOPT_GS      ,MXDISC_GS  ,MXNZON_GS
     ;,NGROUP_ZN        ,MXGRPZN      ,DIVZN_GS   ,SUPBL_GS
     ;,POSDIS_GS    ,IPOLDRIFT_GS,SEARCH_GS
     ;,TRIM_GS          ,EXDRZN_GS    ,MXMEASPP_GS,IDIMVAR_GS
     ;,MXNPP_GS         ,MXNVAR_GS    ,POSMEAS_GS ,VMEAS_GS
     ;,VSTATS_GS        ,VARIO_GS
     ;,IVARIO_GS        ,IDIMIVARIO_GS,IDIMWGT    ,IPNT_PAR
     ;,WGT_PAR          ,IOINV_GS     ,POSZN_GS   ,COORDGR_GS
     ;,NPARDET          ,LDIM_GS      ,WGT_UNK    ,PARC
     ;,PAR_WGT          ,IOSMFL     ,IOSMTP  ,PARC_GS
     ;!NUEVOS
     ;,IODENS_INI,ITPTVAR,BETAC,CREF,DENSREF,TEMPREF,VISCREF,WSPECHEAT
     &,WTHERMCON, PARNAME)

C------Sets reference variable depending on what is solved
C------If solving mass fraction, reference variable is temperature.
	 IF (ITPTVAR.EQ.0) THEN
		VAR_REF = TEMPREF
	 ELSE
C------If solving temperature, reference variable is mass fraction.
            VAR_REF = CREF
	    CREF = TEMPREF
	 END IF

       IF (IODIRECT.EQ.0) CALL MAX_CONECT
     &          (KXX      ,LMXNDL  ,LNNDEL    ,MAINF   ,MAXNB   ,NUMEL
     &          ,NUMNP)


       IF (IOEQT.EQ.0) STOP 'END OF DATA READING'

C------------------------- Momentaneamente se asignan los indices

       INTRA =1
       INSTG =2
       INARR =3
       INARRT=4
       INDSP =5
       INDFM =6
       INPOR =7
       INFOD =8
       INCRD =9
       INCOE =10
       INCHP =1
       INCHPT=2
       INQQP =3
       INQQPT=4
       INALF =5
       INALFT=6
       INCON =7
       INCONT=8
       INDMT =9


       CALL PRODAT
     ;(IDIMBB   ,IERROR   ,IODIM    ,IOEQT    ,IOFLLI   
     ;,IOTRLI   ,LDIM     ,LMXNDL   ,MAINF    ,NTDMT    ,NUMEL
     ;,NUMNP    ,ACTH     ,AREA     ,BIBI     ,BTRA     ,CBASE
     ;,CCALAN   ,CNST     ,COORD    ,GRAV     ,GRAVEL   ,GRDFF
     ;,HBASE    ,HCALAN   ,KXX  ,LXPAREL(1,1) ,LNNDEL   ,LTYPE
     ;,VOLNOD   ,COORD(1,1),COORD(1,2),COORD(1,3)
     ;!NUEVOS
     ;,IODENS_INI   ,IPARTNER ,LCOORD   ,POINTWEIGHT,GP_COORD
     &,IDIMGRAVEL)


C____________________________ makes the arrays iad, iadd, iadn that are
C____________________________ used by watsolv
c      IF (IWATSOLV.NE.0) THEN
	   
	   IF (IODIRECT.EQ.0) CALL MAKE_ADJACENCY_ARRAYS
     ;(IODENS_INI ,IPAR_DIR(22),NUMEL      ,NUMNP    
     ;,LMXNDL     ,IPAR_DIR(23),MAXNB    ,MAXNBF    ,NPARALG    
     ;,IPAR_DIR   ,IAD_D      ,IADD_D   ,IADN_D
     ;,IAD_S      ,IADD_S     ,IADN_S   ,IAFD_S   ,IAFDD_S
     ;,IAFDN_S    ,IAFD_D     ,IAFDD_D  ,IAFDN_D  ,KXX
     ;,LNNDEL)


       IF (IFLAGS(2).EQ.1) CALL WRI_ZONE
     ;(   LXPAREL    ,IXPARNP    ,KXX        ,LTYPE      ,LNNDEL
     ;   ,COORD(1,1) ,COORD(1,2) ,COORD(1,3) ,CFPARNP    ,CFPAREL
     ;   ,NPAREL     ,NPARNP     ,NUMNP      ,NUMEL      ,LMXNDL
     ;   ,INTRA      ,INARR      ,INARRT     ,INSTG      ,INDSP
     ;   ,INDFM      ,INCOE      ,INPOR      ,INCRD      ,INFOD
     ;   ,INCHP      ,INCHPT     ,INQQP      ,INQQPT     ,INCON
     ;   ,INALF      ,INALFT     ,INCONT     ,INDMT      ,MAINF 
     ;   ,'ANTES DE INVER' )

C________________ Loop over conditional simulations. 1 if conditional estimation

       DO ISIM_GS=1,MXSC_GS

          NFL_SIM=MAX(1,IOPTS(28)*NPBFL)  ! # of simultaneous flow pb
          NTP_SIM=MAX(1,IOPTS(29)*NPBTP)  ! # of simultaneous tpt pb

         CALL INIT_GEOEST
     ;(DRELMX_ORIG  ,DABSMX_ORIG  ,ISIM_GS    ,IDIMQ    ,IERROR
     ;,IODIM        ,IOEQT        ,IORTS      ,IOTRS    ,MAINF
     ;,MAXITER_ORIG ,NPARALG      ,NUMEL      ,NUMNP    ,NWRITE
     ;,XMARQ_ORIG   ,ACTH         ,CAUDAL
     ;,CCALIT       ,HCALIT       ,IOWRITE    ,IPAR_INV ,LDIM
     ;,LTYPE        ,PAR_DIR      ,PAR_INV    ,QXYZ
     ;,VD           ,XNORVD       ,NFL_SIM    
     ;,NTP_SIM      ,PAREL        ,PARNP      ,NPPEL    ,NPPNP
     ;,IOVAR        ,IOVAR_ORIG)

        CALL INVER
     &           (A_COUPL_DSC        ,A_COUPL_DSCF       ,ACTH
     &           ,AFLU     ,AFLUDSC  ,AFLUDSCF ,ALFA     ,AREA
     &           ,ATRA     ,ATRADSC  ,ATRADSCF ,BCOUPLED ,BETAC    ,BFLU
     &           ,BIBI     ,BM_ND_FL ,BM_ND_TT ,BM_ZN_FL
     &           ,BM_ZN_TT ,BTRA     ,BUOYANCY ,DBUOYANCY
     &           ,CAUDAL   ,CAUX1    ,CAUX2
     &           ,CCALAN   ,CCALIT   ,CFLU
     &           ,CFPAREL  ,CFPARNP  ,CONCFLOW ,COORD    ,COVINV
     &           ,COVPAR   ,CPREV1   ,CPREV2   ,CREF     
     &           ,DAT_VD   ,DBFLUDFLU,DBFLUDTRA
     &           ,DELTAITER,DENSITY  ,DENSREF
     &           ,DER_VISC 
     &           ,DERC     ,DERH
     &           ,DFLU     ,DFLUDFLU ,DFLUDTRA ,DNODALRH ,DPARELDC
     &           ,DPARELDH ,DQDFLU   ,DQDTRA
     &           ,DTMXDS   ,DTPREVINV,DTRA     ,DTRADFLU
     &           ,DTRADTRA ,DVDC     ,DVDH     ,DVDP     ,DVOBS    ,DWDH
     &           ,FILENAME ,FNT
     &           ,FOBJ_WGT ,GP_COORD ,GRAD     ,GRAVEL   
     &           ,GRDFF    ,HAUX1    ,HAUX2    ,HBASE    ,HCALAN
     &           ,HCALIT   ,HESS     ,HESSAUX  ,HPREV1
     &           ,HPREV2
     &           ,IA_COUPLED_DSC_COLS,IA_COUPLED_DSC_ROWS,IAD_D
     &           ,IAD_S    ,IADD_D   ,IADD_S   ,IADN_D
     &           ,IADN_S   ,IAFD_D   ,IAFD_S   ,IAFDD_D  ,IAFDD_S
     &           ,IAFDN_D  ,IAFDN_S  ,IAFLUDSC_COLS      ,IAFLUDSC_ROWS
     &           ,IATRADSC_COLS      ,IATRADSC_ROWS      ,IBCOD
     &           ,IBTCO    ,IDIMAFLU
     &           ,IDIMATRA ,IDIMBB   ,IDIMCFLU ,IDIMCOV  ,IDIMDENS
     &           ,IDIMDERH ,IDIMDERC ,IDIMDFLU ,IDIMDQ   ,IDIMDTRA
     &           ,IDIMFNT  ,IDIMGRAVEL         ,IDIMHESS
     &           ,IDIMQ    ,IDIMWORK
     &           ,IFLAGS   ,INDEXNOD ,INDPAR   ,INORPAR
     &           ,IOCONSRC ,IOCRITRAP,IODENS_INI
     &           ,IODEVICE ,IODIM    ,IODIRECT ,IOEQT    ,IOFLSAT
     &           ,IOFLLI   ,IOFMLF   ,IOFMLT     ,IOINV
     &           ,IOLG_PAR ,IOPINITC ,IOPINITH ,IOPINVDT,IOPTS
     &           ,IORTS    ,IOTRLI   ,IOTRS    ,IOWRITE
     &           ,IPAR_DIR ,IPAR_INV ,IPARTNER ,ISOLEQ
     &           ,ISOT     ,ISOZ     ,IOSPARSE ,IPAR_DIR(25)
     &           ,ITPTVAR  ,ITYPACOUPLDSC      ,ITYPAFLU
     &           ,ITYPATRA ,ITYPCFLU
     &           ,ITYPDFLU ,ITYPDTRA
     &           ,ITYPFLUDSC,ITYPTRADSC        ,IVPAR    ,IXPARNP
     &           ,KINT     ,KXX      ,LDIM
     &           ,LINMET   ,LMXNDL   ,LNNDEL   ,GRADLOC  ,LTYPE
     &           ,LXPAREL  ,MAINF    ,MAXNB    ,MAXNBF
     &           ,2*NUMNP  ,MAXPG    ,MEASTYP  ,NBAND
     &           ,NBAND1   ,NBANDCOV ,NDEVS
     &           ,NFLAGS   ,NFNL
     &           ,NFNLPAR  ,NFNLPRG  ,NFNLTIP  ,NFTPAR
     &           ,NINT     ,NMAXF    ,NMAXT    ,NOOBSIT  ,NOPTS
     &           ,NPAR     ,NPARALG
     &           ,NPAREL   ,NPARF    ,NPARNP
     &           ,NPBFL    ,NPBMX    ,NPBTP    
     &           ,NPPEL    ,NPPNP    ,NROW     ,NSTAT    
     &           ,NTDMT    ,NTYPAR   ,NTYPEL   ,NUMEL    ,NUMNP
     &           ,NUMTIT   ,NUMTNOD  ,NUMTOBS  ,NWRITE
     &           ,NZONE_PAR,NZPAR    ,NZPRG    ,NZTRA    ,DLT_PAR
     &           ,PAR_DIR  ,PAR_INV  ,PAR_WGT  ,PARACD   ,PARAUX
     &           ,PARZ     ,PAREL    ,PARM     ,PARNP    ,POINTWEIGHT
     &           ,PARZ(INORPAR(19)+1)
     &           ,QXYZ     ,SOLUTION
     &           ,SOURCE             ,TIME     ,TIT      ,TOBS
     &           ,VAR_REF  ,VD       ,VISCOSITY,VISCREF  ,VJAC
     &           ,VOBS     ,VOBSC    ,WATVOL   ,WORK
     &           ,WTOBSN   ,WTOBST   ,XNORVD   ,HINI     ,WSPECHEAT
     &           ,WTHERMCON
     ; ,IOINV_GS    ,IDIMDATASC_GS
     ; ,MXGRPZN,NGROUP_ZN,MXMEASPP_GS,IDIMVAR_GS,MXNVAR_GS
     ; ,IDIMCROSS_GS,MXKRIG_GS,IDIMIVARIO_GS,MXNPRIM_GS
     ; ,MXNZON_GS,MXCLOSE_GS,MXNPP_GS,IDIMZONPP_GS,MXSB_GS
     ; ,MXDISC_GS,MXROT_GS,MXSAM_GS,IOPT_GS,IO_KG_GS
     ; ,IZN_PP_GS,IZN_NPP_GS,IPOLDRIFT_GS,IVARIO_GS,NUMSB_GS,ISUPBL_GS
     ; ,ICROSSCOV_GS,POSMEAS_GS,VMEAS_GS,VSTATS_GS,SEARCH_GS
     ; ,TRIM_GS,SUPBL_GS,KRISYS_GS,CLOSESAM_GS,ZNWGT_GS,KRISOL_GS
     ; ,VARIO_GS,ROTMAT_GS,POSZN_GS,KRIGAUX_GS,CROSSCOV_GS
     ; ,ESTKRIG_GS,POSDISAUX_GS,POSDIS_GS,EXDRZN_GS,MXZONPP_GS
     ; ,DATASC_GS,IDIMWGT,IPNT_PAR,WGT_PAR,NPARDET,ICHECK_GS,LDIM_GS
     ; ,COORDGR_GS,PARC,WGT_UNK,IPOS,DERIV,PARGOOD,ISIM_GS
     ; ,PARC_GS,IFLAG_SIMUL,COVPAR_GR, PARNAME)

        DEALLOCATE(PARNAME)
C------------------------- Write results

         CALL WRI 
     ; (IOWRITE(12),IOWRITE(11),IOWRITE(8)
     ; ,IOWRITE(7),IOWRITE(6)
     ; ,IOWRITE (5),IORTS  ,IOPTS(28),IOPTS(29),IOTRS
     ; ,IPROCESS ,LMXNDL   ,NINT     ,NPBFL    ,NPBTP,NUMEL
     ; ,NUMNP    ,NDEVS    ,CCALIT     ,NUMTIT
     ; ,NUMTOBS  ,FILENAME ,HCALIT
     ; ,VOBSC    ,VOBS     ,IBCOD    ,IBTCO
     ; ,IXPARNP(1,INALF) ,IXPARNP(1,INCHP)
     ; ,IXPARNP(1,INCON) ,IXPARNP(1,INCONT),IXPARNP(1,INQQP)
     ; ,KXX    ,LNNDEL
     ; ,LXPAREL(1,INARR) ,LXPAREL(1,INARRT),LXPAREL(1,INCRD)
     ; ,LXPAREL(1,INDFM) ,LXPAREL(1,INDSP) ,LXPAREL(1,INPOR)
     ; ,LXPAREL(1,INSTG) ,LXPAREL(1,INTRA)
     ; ,DEVNAME    ,TIME
     ; ,COORD(1,1),COORD(1,2),COORD(1,3),IOWRITE(2),MAINF,IODEVICE,TIT
     ; ,ISIM_GS  ,ISOLEQ)

C-parche-escribe las presiones en cada nudo en el último tiempo
c         OPEN(UNIT=181,FILE='presiones.dat',STATUS='UNKNOWN')
c           DO I=1,NUMNP
c             PRESSURE= 1000.0*9.81* (HCALIT(I,1)-COORD(I,2))
c             WRITE(181,3743) I,PRESSURE
c 3743        FORMAT(I5,F20.8)
c           END DO
c         CLOSE(181)
c-fin-parche

C------------------------- Write some statistics

         IF (IOVAR.NE.0) CALL STAT_OUTPUT
     ;(IDIMCOV    ,IDIMHESS   ,IDIMWORK  ,IODIM     ,IOVAR     ,IOPRHED  
     ;,ISOT       ,MAINF      ,NBANDCOV  ,NDEVS     ,NFLAGS    ,NPAR      
     ;,NSTAT      ,NTYPAR     ,NUMTOBS   ,NZPAR     ,COVINV    ,COVPAR    
     ;,DEVICESTAT ,DEVNAME    ,GRAD      ,EIGENVEC  ,FOBJ_WGT  ,HESS 
     ;,IFLAGS     ,INORPAR    ,IODEVICE  ,IOLG_PAR  ,ITYPEPAR  ,IVPAR     
     ;,NZONE_PAR  ,OBSCLASS  ,PARC      ,PARM      ,PAR_WGT   
     ;,RESID      ,STPAR      ,TOBS      ,RESIDPAR  ,VJAC      ,VOBS
     ;,VOBSC      ,WORK       ,IOINV     ,MEASTYP
     ;,WGT_UNK    ,IPNT_PAR  ,IDIMWGT   ,IOPT_GS,MXGRPZN
     ;,NPARDET)

       END DO   ! Next conditional simulation



       RETURN 
       END
