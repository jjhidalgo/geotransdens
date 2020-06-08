       SUBROUTINE PROVISIONAL
     ; (NPARALG  ,NSTAT    ,NTYPAR   ,NWRITE
     ; ,FOBJ_WGT ,IOLG_PAR ,IOWRITE  ,IPAR_DIR
     ; ,IPAR_INV ,NZONE_PAR,PAR_DIR  ,PAR_INV  ,PAR_WGT)

*****************************************************************************
* PURPOSE
*     Inicializacion de todas las variables del COMMON famoso
*
* DESCRIPTION
*     En esta subrutina se incluye el COMMON y se inicializan todas las 
*     variables que ahora se han agrupado en varios vectores para 
*     simplificar las llamadas. El objetivo de esta rutina es conseguir 
*     que las subrutinas que siguen teniendo el COMMON sigan funcionando como 
*     hasta ahora. Esto es pasajero, unica,mente se pretende comprobar 
*     que las nuevas modificaciones funcionan adecuadamente. Una vez visto, 
*     esta subrutina de vida efimera deberia desaparecer
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FOBJ_WGT               Array containing all objective function weights for   
*                         state variables (heads, concentrations, etc)          
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  IOLG_PAR               Array containing all logarithmic options of           
*                         estimated parameters                                  
*  IOWRITE                Array containing all output options                   
*  IPAR_DIR               Array containing all integer direct problem           
*                         parameters                                            
*  IPAR_INV               Array containing all integer inverse problem          
*                         parameters                                            
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters                                            
*  PAR_DIR                Array containing all real direct problem              
*                         parameters                                            
*  PAR_INV                Array containing all real inverse problem             
*                         parameters                                            
*  PAR_WGT                Array containing objective function weights for       
*                         all estimated parameters                              
*
* EXTERNAL VARIABLES: SCALARS
*
*  COSMIN                 XMARQ is multiplied by NUMIN if the cosinus of the    
*                         angle between the angle and the parameters increment  
*                         is less then COSMIN during MAXICOS iterations.        
*  DABSMX                 Absolute convergence criterion                        
*  DABSMX1F                                                                     
*  DABSMX1T                                                                     
*  DABSMX2F                                                                     
*  DABSMX2T                                                                     
*  DCITMX                 Testing factor                                        
*  DHITMX                 Testing factor                                        
*  DMINF                  Convergence criterion                                 
*  DRELMX                 Relative convergence criterion                        
*  DRELMX1F                                                                     
*  DRELMX1T                                                                     
*  DRELMX2F                                                                     
*  DRELMX2T                                                                     
*  EPS                    Convergence criterion. Algorithm stops if maximum     
*                         relative change in one parameter is smaller than EPS  
*  EPSFLU                 Time weighting parameter for nonlinear flow problems  
*  EPSTRA                 Time weighting parameter for nonlinear transport      
*                         problems                                              
*  ERRDMS                                                                       
*  FCTDEC                 Time increment decreasing factor                      
*  FCTDVNR                Testing factor                                        
*  FCTINC                 This factor increases the desirable time increment    
*                         of a generic observation interval j, when it is       
*                         smaller than the current time increment               
*  FCTNCV                 Time increment increasing factor                      
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  GMNOR                  Algorithm stops if gradient norm becomes smaller      
*                         than GMNOR                                            
*  GMNOR1                 Convergence criterion                                 
*  INPWR                  Allows writing on MAIN FILE                           
*  IOCMC                  If 1, computed vs. measured values of concent. at all 
*                         observation points (steady state) are written in      
*                         CVM file                                              
*  IOCMH                  If 1, computed vs. measured values of heads at all    
*                         observation points (steady state) are written in      
*                         CVM file                                              
*  IOCRITRAP              Option of treatement on the direct problem            
*                         convergence criteria                                  
*  IOLGALF                Leakage log-scaling index.                            
*  IOLGARR                Areal recharge  log-scaling index.                    
*  IOLGCHP                Prescribed head log-scaling index.                    
*  IOLGCOE                External concentration log-scaling index.             
*  IOLGCRD                Retardation log-scaling index.                        
*  IOLGDFM                Molecular diffusion log-scaling index.                
*  IOLGDSP                Dispersivity log-scaling index.                       
*  IOLGFOD                                                                      
*  IOLGPOR                Porosity log-scaling index.                           
*  IOLGPRG                Generic parameters log-scaling index.                 
*  IOLGQQP                Prescribed flow log-scaling index.                    
*  IOLGSTG                Storage coefficient log-scaling index.                
*  IOLGTRA                Transmissivity log-scaling index                      
*  IOMHC                  If non zero, computed values of concent. at all nodes 
*                         are writen in MCC file every IOMCC observation times  
*  IOMHH                  If non zero, computed values of heads at all nodes    
*                         are writen in MHH file every IOMHH observation times  
*  IOPINITC               Option for the extrapolation of concentrations        
*                         at the next time step in the Newton process.          
*  IOPINITH               Option for the extrapolation of heads or pressures    
*                         at the next time step in the Newton process.          
*  IOPLC                  Controls when computed and/or measured concent. are   
*                         written in PLT file                                   
*  IOPLH                  Controls when computed and/or measured heads are      
*                         written in PLT file                                   
*  IOSEC                  If non zero, time evolution of computed concent.      
*                         at IOSEC 1-D sections is written in SEC file          
*  IOSEH                  If non zero, time evolution of computed heads         
*                         at IOSEH 1-D sections is written in SEC file          
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IOWNR                  Option of printing information about the direct       
*                         problem iterative process evolution                   
*  IOWRC                  In non zero, computed concentrations in all nodes are 
*                         written in main file every IOWRC observation time     
*  IOWRH                  In non zero, computed heads in all nodes are written  
*                         in main file every IOWRH observation time             
*  ITRAPMX                Maximum number of direct problem iterations permited  
*                         in order to look for convergence before time step     
*                         reduction                                             
*  MAXICOS                Number of successive iterations where the previous    
*                         criterion can be violatded                            
*  MAXITER                Maximum number of iterations                          
*  MINCAT                 Minimum number of consecutives convergences,          
*                         achieved after a time step reduction, required        
*                         before to increasing the time increment when it is    
*                         smaller than the desirable time step                  
*  MXNRTF                 Maximum consecutive time steps reduction permitted    
*                         in order to achieve convergence before stopping.      
*  MXNRTT                 Idem for MXNTRF for transport equation                
*  NMTERF1                Maximum number of failed iterations                   
*  NPARALG                Maximum number of algorithm parameters, including     
*                         minimization and simulation (used for dimensioning)   
*  NSTAT                  Maximum number of state variables whose data is used  
*                         for calibration (used for dimensioning)               
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMAX                  Value to multiplcate XMARQ in apropiate iterations    
*  NUMIN                  Value to divide XMARQ in apropiate iterations         
*  NWRITE                 Number of output options (used for dimensioning)      
*  NZALF                  Number of leakage zones                               
*  NZARR                  Number of areal recharge zones                        
*  NZCHP                  Number of prescribed head zones                       
*  NZCOE                  Number of external concentration zones                
*  NZCRD                  Number of retardation Coefficient zones               
*  NZDFM                  Number of molecular difusion zones                    
*  NZDMT                  Number of matrix diffusion zones                      
*  NZDSP                  Number of dispersivity zones                          
*  NZFOD                  Number of zones of first order decay                  
*  NZPOR                  Number of porosity zones                              
*  NZPRG                  Total number of generic parameter zones               
*  NZQQP                  Number of prescribed flow zones                       
*  NZSTG                  Number of storage Coefficient zones                   
*  NZTRA                  Number of transmissivity zones                        
*  OBJCON1                                                                      
*  OBJCON2                                                                      
*  OBJHED1                Value of objective function considered close to that  
*                         corresponding to the minimum of the inverse problem   
*  OBJHED2                Value of the objective function considered notably    
*                         bigger than that expected to find the minimum of the  
*                         inverse problem                                       
*  PERMX1                 Maximum change per iteration of log-transformed       
*                         variables                                             
*  PERMX2                 Maximum relative change per iteration for the rest of 
*                         parameters                                            
*  PHIMAX                 If the above ratio is greather than PHIMAX, it is     
*                         considered a good quadratic aproximation.             
*  PHIMIN                 If the ratio between the actual change on the         
*                         objective function and its quadratic aproximation is  
*                         smaller than PHIMIN, it is considered a poor          
*                         quadratic aproximation.                               
*  RESIDMX1F              Mass balance error per unit time permited over the    
*                         volume corresponding to a node, having guarantee      
*                         of a suitable accuracy                                
*  RESIDMX1T              Mass balance error per unit time permited over the    
*                         volume corresponding to a node, having guarantee      
*                         of a suitable accuracy                                
*  RESIDMX2F              Idem to RESIDMX1 but admitting certaing degree of     
*                         precision error in the solution                       
*  RESIDMX2T              Idem to RESIDMX1 but admitting certaing degree of     
*                         precision error in the solution                       
*  RESIDMXF               Maximum value of mass balance error over the volume   
*                         represented by a node, considered admisible in order  
*                         to assume convergence of the non linear flow equation 
*  RESIDMXT               Same as above with transport                          
*  THETAF                 Time weighting parameter for flow problems            
*  THETAT                 Time weighting parameter for transport problems       
*  XLAMALF                Weighting coefficient of leakage in the               
*                         objective function                                    
*  XLAMARR                Weighting coefficient of areal recharge in the        
*                         objective function                                    
*  XLAMCHP                Weighting coefficient of prescribed head in the       
*                         objective function                                    
*  XLAMCOE                Weighting coefficient of external concent. in the     
*                         objective function                                    
*  XLAMCON                Objective function weight for concentrations          
*  XLAMCRD                Weighting coefficient of retardation in the           
*                         objective function                                    
*  XLAMDFM                Weighting coefficient of  molecular diffusion in the  
*                         objective function                                    
*  XLAMDSP                Weighting coefficient of dispersivity in the          
*                         objective function                                    
*  XLAMLAM                                                                      
*  XLAMPOR                Weighting coefficient of porosity in the              
*                         objective function                                    
*  XLAMPRGF               Weighting coefficient of flow generic param. in the   
*                         objective function                                    
*  XLAMPRGT               Weighting coefficient of transp. generic par. in the  
*                         objective function                                    
*  XLAMQQP                Weighting coefficient of prescribed flow in the       
*                         objective function                                    
*  XLAMSTG                Weighting coefficient of storage in the               
*                         objective function                                    
*  XLAMTRA                Weighting coefficient of transmissivity in the        
*                         objective function                                    
*  XMARQ                  Initial value of Marquardts parameter (0.0)           
*  ZEROF                  Minimum difference between the values of two          
*                         successives computed heads, considered admisible in   
*                         order to apply the convergence criterium given by     
*                         DRELMAX                                               
*  ZEROT                  Idem to ZEROT but in the transport equation           
*
* HISTORY
*
*     AMS      2-1998     First coding (y ultimo, porque se va a usar una 
*                         sola vez)
*****************************************************************************

       IMPLICIT REAL*8(A-H,O-Z)

C______________________________

       INCLUDE 'COMMON.FOR'

       DIMENSION NZONE_PAR(NTYPAR),IOLG_PAR(NTYPAR),IPAR_INV(NPARALG),
     ;   IPAR_DIR(NPARALG),IOWRITE(NWRITE),
     ;   PAR_WGT(NTYPAR),PAR_INV(NPARALG),PAR_DIR(NPARALG),
     ;   FOBJ_WGT(NSTAT)

        NZTRA      =NZONE_PAR (1)
        NZSTG      =NZONE_PAR (2)
        NZARR      =NZONE_PAR (3)
        NZCHP      =NZONE_PAR (4)
        NZQQP      =NZONE_PAR (5)
        NZALF      =NZONE_PAR (6)
        NZDSP      =NZONE_PAR (7)
        NZDSP      =NZONE_PAR (8)
        NZDFM      =NZONE_PAR (9)
        NZPOR      =NZONE_PAR(10)
        NZFOD      =NZONE_PAR(11)
        NZCRD      =NZONE_PAR(12)
        NZCOE      =NZONE_PAR(13)
        NZPRG      =NZONE_PAR(14)
        NZDMT      =NZONE_PAR(16)
        INPWR      =IOWRITE (1)
        IOWAR      =IOWRITE (2)
        IOWRH      =IOWRITE (3)
        IOWRC      =IOWRITE (4)
        IOPLH      =IOWRITE (5)
        IOPLC      =IOWRITE (6)
        IOMHH      =IOWRITE (7)
        IOMHC      =IOWRITE (8)
        IOSEH      =IOWRITE (9)
        IOSEC      =IOWRITE(10)
        IOCMH      =IOWRITE(11)
        IOCMC      =IOWRITE(12)
        IOWNR      =IOWRITE(13)
        XLAMTRA    =PAR_WGT (1)
        XLAMSTG    =PAR_WGT (2)
        XLAMARR    =PAR_WGT (3)
        XLAMCHP    =PAR_WGT (4)
        XLAMQQP    =PAR_WGT (5)
        XLAMALF    =PAR_WGT (6)
        XLAMDSP    =PAR_WGT (7)
        XLAMDSP    =PAR_WGT (8)
        XLAMDFM    =PAR_WGT (9)
        XLAMPOR    =PAR_WGT(10)
        XLAMLAM    =PAR_WGT(11)
        XLAMCRD    =PAR_WGT(12)
        XLAMCOE    =PAR_WGT(13)
        XLAMPRGF   =PAR_WGT(14)
        XLAMPRGT   =PAR_WGT(16)
        XMARQ      =PAR_INV (1)
        PHIMIN     =PAR_INV (2)
        PHIMAX     =PAR_INV (3)
        GMNOR1     =PAR_INV (4)
        GMNOR      =PAR_INV (5)
        DMINF      =PAR_INV (6)
        COSMIN     =PAR_INV (7)
        PERMX1     =PAR_INV (8)
        PERMX2     =PAR_INV (9)
        EPS        =PAR_INV(10)
        NUMIN      =IPAR_INV (1)
        NUMAX      =IPAR_INV (2)
        MAXICOS    =IPAR_INV (3)
        MAXITER    =IPAR_INV (4)
        NMTERF1    =IPAR_INV (5)
        IOLGTRA    =IOLG_PAR (1)
        IOLGSTG    =IOLG_PAR (2)
        IOLGARR    =IOLG_PAR (3)
        IOLGCHP    =IOLG_PAR (4)
        IOLGQQP    =IOLG_PAR (5)
        IOLGALF    =IOLG_PAR (6)
        IOLGDSP    =IOLG_PAR (7)
        IOLGDSP    =IOLG_PAR (8)
        IOLGDFM    =IOLG_PAR (9)
        IOLGPOR    =IOLG_PAR(10)
        IOLGFOD    =IOLG_PAR(11)
        IOLGCRD    =IOLG_PAR(12)
        IOLGCOE    =IOLG_PAR(13)
        IOLGPRG    =IOLG_PAR(14)
        FCTNCV     =PAR_DIR (1)
        FCTDEC     =PAR_DIR (2)
        FCTINC     =PAR_DIR (3)
        DRELMX     =PAR_DIR (4)
        DABSMX     =PAR_DIR (5)
        RESIDMXF   =PAR_DIR (6)
        RESIDMXT   =PAR_DIR (7)
        ZEROF      =PAR_DIR (8)
        ZEROT      =PAR_DIR (9)
        FCTDVNR    =PAR_DIR(10)
        DHITMX     =PAR_DIR(11)
        DCITMX     =PAR_DIR(12)
        OBJHED1    =PAR_DIR(13)
        RESIDMX1F  =PAR_DIR(14)
        DABSMX1F   =PAR_DIR(15)
        DRELMX1F   =PAR_DIR(16)
        OBJHED2    =PAR_DIR(17)
        RESIDMX2F  =PAR_DIR(18)
        DABSMX2F   =PAR_DIR(19)
        DRELMX2F   =PAR_DIR(20)
        OBJCON1    =PAR_DIR(21)
        RESIDMX1T  =PAR_DIR(22)
        DABSMX1T   =PAR_DIR(23)
        DRELMX1T   =PAR_DIR(24)
        OBJCON2    =PAR_DIR(25)
        RESIDMX2T  =PAR_DIR(26)
        DABSMX2T   =PAR_DIR(27)
        DRELMX2T   =PAR_DIR(28)
        THETAF     =PAR_DIR(29)
        THETAT     =PAR_DIR(30)
        EPSFLU     =PAR_DIR(31)
        EPSTRA     =PAR_DIR(32)
        ERRDMS     =PAR_DIR(33)
        MXNRTF     =IPAR_DIR(1)
        MXNRTT     =IPAR_DIR(2)
        MINCAT     =IPAR_DIR(3)
        IOPINITH   =IPAR_DIR(4)
        IOPINITC   =IPAR_DIR(5)
        IOCRITRAP  =IPAR_DIR(6)
        ITRAPMX    =IPAR_DIR(7)
        XLAMCON    =FOBJ_WGT(2)

        IUDIM=10
        IUGRID=11
        IUPAR=12
        IUTIM=13
        IUOBS=14
        IUCAL=15
        IUVAR=16
        MAINF=25
        IUGR1=26
        IUGR2=27
        IUGR3=28
        IUGR4=29
        IUGR5=30
        IUGR6=31
        IUGR7=32
        IUGR8=33
        IUGR9=34
	IUBLH=35
	IUBLC=36
        IAOBS=81
        IAUXH=94
        IAUXC=93      
        IUBALH=83
        IUBALC=84
        IUTIME=85

C---------------  Asignacion de IZ's

       IZTRA=NZTRA
       IF (IZTRA.EQ.0) IZTRA=1
       IZSTG=NZSTG
       IF (IZSTG.EQ.0) IZSTG=1
       IZALF=NZALF
       IF (IZALF.EQ.0) IZALF=1
       IZCHP=NZCHP
       IF (IZCHP.EQ.0) IZCHP=1
       IZQQP=NZQQP
       IF (IZQQP.EQ.0) IZQQP=1
       IZARR=NZARR
       IF (IZARR.EQ.0) IZARR=1
       IZDSP=NZDSP
       IF (IZDSP.EQ.0) IZDSP=1
       IZCRD=NZCRD
       IF (IZCRD.EQ.0) IZCRD=1
       IZDFM=NZDFM
       IF (IZDFM.EQ.0) IZDFM=1
       IZCOE=NZCOE
       IF (IZCOE.EQ.0) IZCOE=1
       IZPOR=NZPOR
       IF (IZPOR.EQ.0) IZPOR=1

       RETURN
       END
