       SUBROUTINE STATE_VARIABLE_INIT
     ;(IDIMFNT  ,INDFLTR  ,INTI       ,IOINV    ,IOPINIT
     ;,NCONVI   ,NINT     ,NPARALG  ,NPARNP     ,NTYPAR   ,NUMITER
     ;,NUMNP    ,NZPAR    ,THETA    ,TICAL      ,TICALAN  ,TIME
     ;,TINC     ,TINCINI  ,TINCLAST ,TINTERVOBS ,BVAR     ,CFPARNP
     ;,FNT      ,IBVAR    ,INORPAR  ,IXPARNP    ,NFTPAR   ,PARC
     ;,VAUX1    ,VCALIT   ,VCALAN   ,SOLUTION   ,VPREV1   ,VPREV2)


********************************************************************************
*
* PURPOSE To initialize the state variable (pressure head, piezometric head or
*         solute concentration at the beginning of the Newton-Raphson process 
*         when solving the nonlinear flow or transport equation.  
*
* DESCRIPTION According to the option selected by user or by the owun code, 
*             the initial value of state variable is computed by extrapolation 
*             or using information gained at previous inverse problem iteration
*             In this case, these information is stored on a direct access files
*             On the other hand, when extrap. criterium is chosen, the info.
*             is contained on several arrays.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  BVAR                   Array containing BFLU or BTRA.
*  CFPARNP                Array containing node coefinient of node j            
*                         corresponding to INpar index parameter zone.          
*  FNT                    Array containing time functions values                
*  IBVAR                  Array containing IBCOD or IBTCO
*  NFTPAR                 Vector containing time function number at every       
*                         parameter zone                                        
*  PARC                   Vector containing calculated values for all           
*                         parameters                                            
*  TIME                   Array containing observation times.
*  VAUX1                  
*  VCALIT
*  VCALAN                                                                       
*  SOLUTION                                                                       
*  VPREV1                                                                       
*  VPREV2                                                                       
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IDIMFNT                First dimension of array FNT, it coincides with       
*                         NFNT if the latter is not zero                        
*  IOINV                  Inverse problem option
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  IOPINIT                Scalar. It is equal to IOPINITH or to IOPINITC
*  IXPARNP                Array containing zone number of each node j,          
*                         corresponding to INpar index parameter zone.          
*  NCONVI                 Number of convergences reached.
*  NINT                   Number of observation times                           
*  NPARNP                 Number of nodal parameters in current problem         
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMITER                Current iteration in inverse problem process          
*  NUMNP                  Number of nodes                                       
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  TICAL                                                                        
*  TICALAN                                                                      
*  TINC                   Current time increment                                
*  TINCINI                Time increment used two time steps ago                
*  TINCLAST                                                                     
*  TINTERVOBS             Length of the current observation time interval       
*
* INTERNAL VARIABLES: SCALARS
*
*  INDFLTR                Indicator variable of flow or trasnport regime
*  INDINIT                Result of COMPARE_INIT. If =0 extrap.
*                                                 If =1 uses inv. p.info
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  COMPARE_INIT                                                                 
*  INIT_EXTRAP                                                                  
*  INIT_INV_PROB_INFO                                                           
*
* HISTORY
*  First coding: German Galarza Nov-1997
*  Revision: Andres Alcolea Oct-1998
*
********************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)


C______________________________Choose the way to initialize state variable
C______________________________as a function of last optimization iteration
C______________________________or as a function of previous solutions. 


      CALL COMPARE_INIT
     ;(INDFLTR  ,INDINIT  ,INTI     ,IOINV    ,IOPINIT  ,NCONVI
     ;,NINT     ,NUMITER  ,NUMNP    ,TICAL    ,TICALAN  ,TINC
     ;,TINCINI  ,TINCLAST ,BVAR     ,IBVAR    ,TIME     ,VAUX1
     ;,VCALIT     ,VCALAN   ,SOLUTION   ,VPREV1   ,VPREV2)  


C______________________________Extrapolation

      IF (INDINIT.EQ.0) THEN

        CALL INIT_EXTRAP
     ;(IDIMFNT  ,INDFLTR  ,INTI       ,NCONVI   ,NINT
     ;,NPARNP   ,NTYPAR   ,NUMNP    ,NZPAR      ,THETA    ,TICALAN
     ;,TINC     ,TINCINI  ,TINCLAST ,TINTERVOBS ,CFPARNP  ,FNT
     ;,IBVAR    ,INORPAR  ,IXPARNP  ,NFTPAR     ,PARC     ,VCALAN
     ;,SOLUTION   ,VPREV1   ,VPREV2)  

      ELSE

C______________________________Uses previous Marquardt iteration information

        CALL INIT_INV_PROB_INFO
     ;(THETA     ,IDIMFNT  ,INDFLTR  ,INTI     ,NINT     ,NPARALG
     ;,NPARNP    ,NUMNP    ,NZPAR    ,TICAL    ,TICALAN  ,TINC
     ;,TINTERVOBS,BVAR     ,CFPARNP  ,FNT      ,IBVAR    ,INORPAR
     ;,IXPARNP   ,NFTPAR   ,PARC     ,TIME     ,VAUX1    ,VCALAN
     ;,SOLUTION)  

      ENDIF

      RETURN
      END
