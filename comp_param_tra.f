       SUBROUTINE COMP_PARAM_TRA
     &           (CCALAN   ,CCALIT   ,CFPAREL  ,CFPARNP  ,DPARELDC
     &           ,DPARELDH ,DTIMET   ,EPSTRA   ,FNT      ,IBTCO
     &           ,IDIMFNT  ,IFLAGS   ,INCLK    ,INCON    ,INDSSTR
     &           ,INORPAR  ,INTI     ,IODENS   ,IOFMLT   ,IOTRLI
     &           ,IPAR_DIR ,IXPARNP  ,KXX      ,LMXNDL   ,LNNDEL
     &           ,LXPAREL  ,MAINF    ,NFLAGS   ,NFNL     ,NFNLPAR
     &           ,NFNLPRG  ,NFNLTIP  ,NFTPAR   ,NINT     ,NPARALG
     &           ,NPAREL   ,NPARNP   ,NPPEL    ,NPPNP    ,NTYPAR
     &           ,NUMEL    ,NUMNP    ,NZONE_PAR,NZPAR    ,NZPRG
     &           ,PARACD   ,PARC     ,PAREL    ,PARNP    ,PRGC
     &           ,XPARAM   ,THETAT   ,TINC     ,TINTERVOBS)

********************************************************************************
*
* PURPOSE
*
*      Computes the values of transport parameters by elements
*
* DESCRIPTION
*
*      Computes the values of transport parameters by elements.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  ACTH                   Aquifer thickness of every element. Cross sectional   
*                         area for 1-D elements, thickness for 2-D elements.    
*  AREA                   Element size (length for 1-D elem, area for 2-D,      
*                         volume for 3-D)                                       
*  CCALAN                 Computed concentrations in the previous time step.    
*  CCALIT                 Computed concentration in last iteration              
*  CFPAREL                Array containing node coefinient of element j         
*                         corresponding to INpar index parameter zone.          
*  FNT                    Array containing time functions values                
*  HCAL                   Computed heads at every node                          
*  HCALAN                 Head level at previous time                           
*  IFLAGS                 Array with different writing options. Used mainly for 
*                         debugging.                                            
*  INORPAR                Array containing the indexes with the location        
*                         of the different paramaters in arrays PARC, PARM,     
*                         IVPAR, NFTPAR, STPAR and FNTPAR                       
*  KXX                    Node numbers of every element (counterclockwise       
*                         order).                                               
*  LNNDEL                 Number of nodes at every element                      
*  LXPAREL                Array containing zone numbers for a given             
*                         element parameter                                     
*  NFNLPAR                Vector containing non-linear function order           
*                         afecting every parameter at each zone.                
*  NFNLPRG                Generic parameter zone number for every nonlinear     
*                         function                                              
*  NFNLTIP                Type of non-linear function                           
*  NFTPAR                 Vector containing time function number at every       
*                         parameter zone                                        
*  NZONE_PAR              Array containing the number of zones of all           
*                         parameters                                            
*  PARACD                 Agreement parameters                                  
*  PARC                   Vector containing calculated values for all           
*                         parameters                                            
*  PAREL                                                                        
*  WATVOL                                                                       
*
* INTERNAL VARIABLES: ARRAYS
*
*  XPARAM                                                                       
*
* EXTERNAL VARIABLES: SCALARS
*
*  DTIMET                 Relative time at which transport right hand side is   
*                         evaluated. If 1, rhs is evaluated at the              
*                         next computing time, if 0, rhs is evaluated           
*                         at the last computing time.                                                                                                 
*  EPSTRA                 Time weighting parameter for nonlinear transport      
*                         problems                                              
*  IDIMFNT                First dimension of array FNT, it coincides with       
*                         NFNT if the latter is not zero                        
*  INDSSTR                Current problem state. If 1, transient state,         
*                         if 0, steady-state                                    
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  IOFMLT                 Transport formulation number                          
*  IOTRS                  Flow regime                                           
*  IPAR_DIR               Array containing all integer direct problem           
*                         parameters                                            
*  LMXNDL                 Maximum number of nodes per element                   
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFLAGS                 Maximum number of allowed flags                       
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NINT                   Number of observation times                           
*  NPARALG                Maximum number of algorithm parameters, including     
*                         minimization and simulation (used for dimensioning)   
*  NPAREL                 Number of element parameters in current problem       
*  NPPEL                                                                        
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMEL                  Number of elements                                    
*  NUMNP                  Number of nodes                                       
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*  PRGC                                                                                                                                                  
*
* INTERNAL VARIABLES: SCALARS
*
*  NNUD                   Number of nodes of the current element                
*  NZPRG                  Total number of generic parameter zones               
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  PARAM_VALUE                                                                  
*
* HISTORY
*
*     AMS      1-1999     First coding
*     AMS      3-1998     Revision
*
*  NOTE This subroutine was formerly called COMP_ELEM_VAR.
*
********************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IDIMFNT  ,INCLK    ,INCON    ,INDSSTR  ,INTI
     &          ,IODENS   ,IOFMLT   ,IOTRLI   ,LMXNDL   ,MAINF
     &          ,NFLAGS   ,NFNL     ,NINT     ,NPARALG  ,NPAREL
     &          ,NPARNP   ,NPPEL    ,NPPNP    ,NTYPAR   ,NUMEL
     &          ,NUMNP    ,NZPAR    ,NZPRG

      INTEGER*4::IBTCO(NUMNP)           ,IFLAGS(NFLAGS)
     &          ,INORPAR(NTYPAR)        ,IPAR_DIR(NPARALG)
     &          ,IXPARNP(NUMNP,NPARNP)  ,KXX(LMXNDL,NUMEL)
     &          ,LXPAREL(NUMEL,NPAREL)  ,LNNDEL(NUMEL)
     &          ,NFNLPAR(NZPAR)         ,NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL)          ,NFTPAR(NZPAR)
     &          ,NZONE_PAR(NTYPAR)

     

      REAL*8::DERSCC,DERSCH,DTIMET,EPSTRA,THETAT,TINC,TINTERVOBS


      REAL*8::CCALAN(NUMNP)        ,CCALIT(NUMNP)
     &       ,CFPAREL(NUMEL,NPAREL),CFPARNP(NUMNP,NPARNP)
     &       ,DPARELDC(NPPEL,NUMEL),DPARELDH(NPPEL,NUMEL)
     &       ,FNT(IDIMFNT,NINT)    ,PARACD(3,NFNL)
     &       ,PARC(NZPAR)          ,PAREL(NUMEL,NPPEL)
     &       ,PARNP(NUMNP,NPPNP)   ,PRGC(NZPRG)
     &       ,XPARAM(8)

C------------------------- Internal

      INTEGER*4::L,NNUD,JJ,J,I,IB
      INTEGER*4::NFNLCLK  ,NFNLCOE  ,NFNLDFM  ,NFNLDSL  ,NFNLDST
     &          ,NFNLFOD  ,NFNLMDF  ,NFNLPOR  ,NFNLRET  ,NFTCLK
     &          ,NFTCOE
     &          ,NFTDFM   ,NFTDSL   ,NFTDST   ,NFTFOD   ,NFTIP
     &          ,NFTMDF,NFTPOR
     &          ,NFTRET

      REAL*8::DTIMETAUX


C------------------------- First executable statement.

      NZPRG = NZONE_PAR(14)

      DO L=1,NUMEL

          NNUD = LNNDEL(L)

C------------------------- Porosity

          JJ = INORPAR(15) + LXPAREL(L,7)
          NFTPOR = NFTPAR(JJ)

          IF (IOTRLI.NE.0) THEN

              NFNLPOR = NFNLPAR(JJ)
              NFTIP = 0
              IF (NFNLPOR.NE.0) NFTIP = NFNLTIP(NFNLPOR)

          ELSE

              NFNLPOR = 0
              NFTIP = 0

          END IF                !IOTRLI.NE.0 
          
          CALL PARAM_VALUE
     &        (CFPAREL(L,7)       ,DERSCC   ,DERSCH   ,DTIMET   ,EPSTRA
     &        ,IDIMFNT  ,INDSSTR  ,INTI     ,IOFMLT   ,L        ,LMXNDL
     &        ,NFLAGS   ,NFNL     ,NFNLPOR  ,NFTIP
     &        ,NFTPOR   ,NINT     ,NNUD     ,NPARALG  ,NUMEL
     &        ,NUMNP    ,NZPRG    ,PARC(JJ) ,PAREL(L,12)        ,PRGC
     &        ,FNT      ,CCALAN   ,CCALIT   ,IFLAGS   ,IPAR_DIR ,KXX
     &        ,NFNLPRG  ,PARACD   ,XPARAM)

C-------------------------  Assigns derivative

          IF (IOTRLI.NE.0) THEN

              DPARELDH(12,L) = DERSCH

              IF (IODENS.EQ.1) THEN

                  DPARELDC(12,L) = DERSCC

              END IF !IODENS.EQ.1

          END IF !IOTRLI.NE.0

          

C------------------------- Molecular diffusion

          IF (NZONE_PAR(9).NE.0) THEN

              JJ = INORPAR(14) + LXPAREL(L,6)
              NFTMDF = NFTPAR(JJ)

              IF (IOTRLI.NE.0) THEN

                  NFNLMDF = NFNLPAR(JJ)
                  NFTIP = 0
                  IF (NFNLMDF.NE.0) NFTIP = NFNLTIP(NFNLMDF)

              ELSE

                  NFNLMDF = 0
                  NFTIP = 0

              END IF            !IOTRLI.NE.0 

              CALL PARAM_VALUE
     &            (CFPAREL(L,6)       ,DERSCC   ,DERSCH   ,DTIMET
     &            ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INTI     ,IOFMLT
     &            ,L        ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLMDF
     &            ,NFTIP    ,NFTMDF   ,NINT
     &            ,NNUD     ,NPARALG  ,NUMEL    ,NUMNP    ,NZPRG
     &            ,PARC(JJ) ,PAREL(L,11)        ,PRGC     ,FNT
     &            ,CCALAN   ,CCALIT   ,IFLAGS   ,IPAR_DIR ,KXX
     &            ,NFNLPRG  ,PARACD   ,XPARAM)

          ELSE

              PAREL(L,11) = 0D0
              DERSCC = 0D0
              DERSCH = 0D0

          END IF !NZONE_PAR(9).NE.0

C-------------------------  Assigns derivative

          IF (IOTRLI.NE.0) THEN  
           
              DPARELDH(11,L) = DERSCH

              IF (IODENS.EQ.1) THEN
                  DPARELDC(11,L) = DERSCC
              END IF !IODENS.EQ.1

          END IF !IOTRLI.NE.0

C------------------------- Longitudinal dispersivity

          IF (NZONE_PAR(7).NE.0) THEN

              JJ = INORPAR(12) + LXPAREL(L,5)
              NFTDSL = NFTPAR(JJ)

              IF (IOTRLI.NE.0) THEN

                  NFNLDSL = NFNLPAR(JJ)
                  NFTIP = 0
                  IF (NFNLDSL.NE.0) NFTIP = NFNLTIP(NFNLDSL)

              ELSE

                  NFNLDSL = 0
                  NFTIP = 0

              END IF            !IOTRLI.NE.0 

              CALL PARAM_VALUE
     &            (CFPAREL(L,5)       ,DERSCC   ,DERSCH   ,DTIMET
     &            ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INTI     ,IOFMLT
     &            ,L        ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLDSL
     &            ,NFTIP    ,NFTDSL   ,NINT
     &            ,NNUD     ,NPARALG  ,NUMEL    ,NUMNP    ,NZPRG
     &            ,PARC(JJ) ,PAREL(L,9)         ,PRGC     ,FNT
     &            ,CCALAN   ,CCALIT   ,IFLAGS   ,IPAR_DIR ,KXX
     &            ,NFNLPRG  ,PARACD   ,XPARAM)

          ELSE

              PAREL(L,9) = 0D0
              DERSCC = 0D0
              DERSCH = 0D0

          END IF !NZONE_PAR(7).NE.0

C-------------------------  Assigns derivative

          IF (IOTRLI.NE.0) THEN   

              DPARELDH(9,L) = DERSCH

              IF (IODENS.EQ.1) THEN

                  DPARELDC(9,L) = DERSCC

              END IF !IODENS.EQ.1

          END IF !IOTRLI.NE.0

C------------------------- Transversal dispersivity

          IF (NZONE_PAR(8).NE.0) THEN

              JJ = INORPAR(13) + LXPAREL(L,5)
              NFTDST = NFTPAR(JJ)

              IF (IOTRLI.NE.0) THEN

                  NFNLDST = NFNLPAR(JJ)
                  NFTIP = 0
                  IF (NFNLDST.NE.0) NFTIP = NFNLTIP(NFNLDST)

              ELSE

                  NFNLDST = 0
                  NFTIP = 0

              END IF            !IOTRLI.NE.0 


              CALL PARAM_VALUE
     &            (CFPAREL(L,5)       ,DERSCC   ,DERSCH   ,DTIMET
     &            ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INTI     ,IOFMLT
     &            ,L        ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLDST
     &            ,NFTIP    ,NFTDST   ,NINT
     &            ,NNUD     ,NPARALG  ,NUMEL    ,NUMNP    ,NZPRG
     &            ,PARC(JJ) ,PAREL(L,10)        ,PRGC     ,FNT
     &            ,CCALAN   ,CCALIT   ,IFLAGS   ,IPAR_DIR ,KXX
     &            ,NFNLPRG  ,PARACD   ,XPARAM)

          ELSE

              PAREL(L,10) = 0D0
              DERSCC = 0D0
              DERSCH = 0D0

          END IF !NZONE_PAR(8).NE.0

C-------------------------  Assigns derivative

          IF (IOTRLI.NE.0) THEN

              DPARELDH(10,L) = DERSCH

              IF (IODENS.EQ.1) THEN

                  DPARELDC(10,L) = DERSCC

              END IF !IODENS.EQ.1

          END IF !IOTRLI.NE.0

C------------------------- Retardation

          IF (NZONE_PAR(12).NE.0) THEN

              JJ = INORPAR(17) + LXPAREL(L,9)

              IF (IOTRLI.NE.0) THEN

                  NFNLRET = NFNLPAR(JJ)
                  NFTIP = 0
                  IF (NFNLRET.NE.0) NFTIP = NFNLTIP(NFNLRET)

              ELSE

                  NFNLRET = 0
                  NFTIP = 0

              END IF !IOTRLI.NE.0 

              CALL PARAM_VALUE
     &            (CFPAREL(L,9)       ,DERSCC   ,DERSCH   ,DTIMET
     &            ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INTI     ,IOFMLT
     &            ,L        ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLPAR(JJ)
     &            ,NFTIP    ,NFTPAR(JJ)         ,NINT
     &            ,NNUD     ,NPARALG  ,NUMEL    ,NUMNP    ,NZPRG
     &            ,PARC(JJ) ,PAREL(L,14)        ,PRGC     ,FNT
     &            ,CCALAN   ,CCALIT   ,IFLAGS   ,IPAR_DIR ,KXX
     &            ,NFNLPRG  ,PARACD   ,XPARAM)

          ELSE

              PAREL(L,14) = 0D0
              DERSCC = 0D0
              DERSCH = 0D0

          END IF !NZONE_PAR(12).NE.0

C-------------------------  Assigns derivative

          IF (IOTRLI.NE.0) THEN   

              DPARELDH(14,L) = DERSCH

              IF (IODENS.EQ.1) THEN

                  DPARELDC(14,L) = DERSCC

              END IF !IODENS.EQ.1

          END IF !IOTRLI.NE.0


C------------------------- First order decay

          IF (NZONE_PAR(11).NE.0) THEN

              JJ = INORPAR(16) + LXPAREL(L,8)

              NFTFOD = NFTPAR(JJ)

              IF (IOTRLI.NE.0) THEN

                  NFNLFOD = NFNLPAR(JJ)
                  NFTIP = 0
                  IF (NFNLFOD.NE.0) NFTIP = NFNLTIP(NFNLFOD)

              ELSE

                  NFNLFOD = 0
                  NFTIP = 0

              END IF            !IOTRLI.NE.0 


              CALL PARAM_VALUE
     &            (CFPAREL(L,8)       ,DERSCC   ,DERSCH   ,DTIMET
     &            ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INTI     ,IOFMLT
     &            ,L        ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLFOD
     &            ,NFTIP    ,NFTFOD   ,NINT
     &            ,NNUD     ,NPARALG  ,NUMEL    ,NUMNP    ,NZPRG
     &            ,PARC(JJ) ,PAREL(L,13)        ,PRGC     ,FNT
     &            ,CCALAN   ,CCALIT   ,IFLAGS   ,IPAR_DIR ,KXX
     &            ,NFNLPRG  ,PARACD   ,XPARAM)
          ELSE

              PAREL(L,13) = 0D0
              DERSCC = 0D0
              DERSCH = 0D0

          END IF !NZONE_PAR(11).NE.0

C-------------------------  Assigns derivative

          IF (IOTRLI.NE.0) THEN   

              DPARELDH(13,L) = DERSCH

              IF (IODENS.EQ.1) THEN

                  DPARELDC(13,L) = DERSCC

              END IF !IODENS.EQ.1

          END IF !IOTRLI.NE.0

C------------------------- External concentration (elements)

          IF (NZONE_PAR(13).NE.0) THEN

              JJ = INORPAR(18) + LXPAREL(L,10)
              NFTCOE = NFTPAR(JJ)

              IF (IOTRLI.NE.0) THEN

                  NFNLCOE = NFNLPAR(JJ)
                  NFTIP = 0
                  IF (NFNLCOE.NE.0) NFTIP = NFNLTIP(NFNLCOE)

              ELSE

                  NFNLCOE = 0
                  NFTIP = 0

              END IF            !IOTRLI.NE.0 


              CALL PARAM_VALUE
     &            (CFPAREL(L,10)      ,DERSCC   ,DERSCH   ,DTIMET
     &            ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INTI     ,IOFMLT
     &            ,L        ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLCOE
     &            ,NFTIP    ,NFTCOE   ,NINT    ,NNUD
     &            ,NPARALG  ,NUMEL    ,NUMNP     ,NZPRG   ,PARC(JJ)
     &            ,PAREL(L,15)        ,PRGC      ,FNT     ,CCALAN
     &            ,CCALIT   ,IFLAGS   ,IPAR_DIR  ,KXX     ,NFNLPRG
     &            ,PARACD   ,XPARAM)

          ELSE

              PAREL(L,15) = 0D0
              DERSCC = 0D0
              DERSCH = 0D0

          END IF !NZONE_PAR(13).NE.0

C-------------------------  Assigns derivative

          IF (IOTRLI.NE.0) THEN

              DPARELDH(15,L) = DERSCH

              IF (IODENS.EQ.1) THEN

                  DPARELDC(15,L) = DERSCC

              END IF !IODENS.EQ.1

          END IF !IOTRLI.NE.0

      END DO !L=1,NUMEL


      IF (IFLAGS(16).EQ.1) THEN

          WRITE(MAINF,2000) 
 2000     FORMAT(10X,'PAREL ARRAY',/,5X,'DSL',10X,'DST',10X,'DFM'
     &           ,10X,'POR',10X,'FOD',10X,'CRD')

          DO L=1,NUMEL

              WRITE(MAINF,2100) (PAREL(L,J),J=9,14)
 2100         FORMAT(6G13.6)

          END DO !L=1,NUMEL

      END IF !IFLAGS(16).EQ.1

C------------------------- External concentration (nodal)
C------------------------- and concentration leakage coefficient (nodal)

      DO I=1,NUMNP

          IB = IBTCO(I)

          IF (IB.NE.0) THEN

              JJ = INORPAR(18) + IXPARNP(I,INCON+INDSSTR)
              NFTCOE = NFTPAR(JJ)

              IF (IOTRLI.NE.0) THEN

                  NFNLCOE = NFNLPAR(JJ)
                  NFTIP = 0
                  IF (NFNLCOE.NE.0) NFTIP = NFNLTIP(NFNLCOE)

              ELSE

                  NFNLCOE = 0
                  NFTIP = 0

              END IF            !IOTRLI.NE.0 

              DTIMETAUX = DTIMET
              IF (IB.EQ.1 .AND. INDSSTR.NE.0) THEN
              
                  DTIMETAUX = DTIMETAUX + (1D0-THETAT)*TINC/TINTERVOBS

	      END IF !IB.EQ.1 .AND. INDSSTR.NE.0

              CALL PARAM_VALUE
     &             (CFPARNP(I,INCON+INDSSTR)   ,DERSCC   ,DERSCH
     &             ,DTIMETAUX
     &             ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INTI     ,IOFMLT
     &             ,I        ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLCOE
     &             ,NFTIP    ,NFTCOE         ,NINT
     &             ,1        ,NPARALG  ,NUMEL    ,NUMNP    ,NZPRG
     &             ,PARC(JJ) ,PARNP(I,4)         ,PRGC     ,FNT
     &             ,CCALAN   ,CCALIT   ,IFLAGS   ,IPAR_DIR ,KXX
     &             ,NFNLPRG  ,PARACD   ,XPARAM)


              IF (IB.EQ.5) THEN

                  JJ = INORPAR(21) + IXPARNP(I,INCLK)
                  NFTCOE = NFTPAR(JJ)

                  IF (IOTRLI.NE.0) THEN

                      NFNLCLK = NFNLPAR(JJ)
                      NFTIP = 0
                      IF (NFNLCLK.NE.0) NFTIP = NFNLTIP(NFNLCLK)

                  ELSE

                      NFNLCLK = 0
                      NFTIP = 0

                  END IF        !IOTRLI.NE.0 

                  CALL PARAM_VALUE
     &                 (CFPARNP(I,INCLK)   ,DERSCC   ,DERSCH   ,DTIMET
     &                 ,EPSTRA   ,IDIMFNT  ,INDSSTR  ,INTI     ,IOFMLT
     &                 ,I        ,LMXNDL   ,NFLAGS   ,NFNL     ,NFNLCLK
     &                 ,NFTIP    ,NFTCLK   ,NINT
     &                 ,1        ,NPARALG  ,NUMEL    ,NUMNP    ,NZPRG
     &                 ,PARC(JJ) ,PARNP(I,6)         ,PRGC     ,FNT
     &                 ,CCALAN   ,CCALIT   ,IFLAGS   ,IPAR_DIR ,KXX
     &                 ,NFNLPRG  ,PARACD   ,XPARAM)


              ELSE

                  CFPARNP(I,INCLK) = 0D0
                  DERSCC = 0D0
                  DERSCH = 0D0

              END IF            !IB.EQ.5

          ELSE

              CFPARNP(I,INCON+INDSSTR) = 0D0
              DERSCC = 0D0
              DERSCH = 0D0

          END IF                !IB.NE.0

      END DO !I=1,NNUD

      IF (IFLAGS(16).EQ.1) THEN

          WRITE(MAINF,3000) 
 3000     FORMAT(10X,'PARNP ARRAY',/,5X,'CHP',10X,'QQP',10X,'ALF'
     &           ,10X,'COE')

          DO I=1,NUMNP

             WRITE(MAINF,2100) (PARNP(I,J),J=1,4)

          END DO !I=1,NUMNP

      END IF !IFLAGS(16).EQ.1


      END SUBROUTINE COMP_PARAM_TRA
