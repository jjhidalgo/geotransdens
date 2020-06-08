       SUBROUTINE WRI_ZONE
     ;(   LXPAREL    ,IXPARNP    ,KXX      ,LTYPE      ,LNNDEL
     ;   ,X          ,Y          ,Z        ,CFPARNP    ,CFPAREL
     ;   ,NPAREL     ,NPARNP     ,NUMNP    ,NUMEL      ,LMXNDL
     ;   ,INTRA      ,INARR      ,INARRT   ,INSTG      ,INDSP
     ;   ,INDFM      ,INCOE      ,INPOR    ,INCRD      ,INFOD
     ;   ,INCHP      ,INCHPT     ,INQQP    ,INQQPT     ,INCON
     ;   ,INALF      ,INALFT     ,INCONT   ,INDMT      ,MAINF 
     ;   ,MSSG    )
***************************************************************
*  PURPOSE: WRITES PARAMETERS ZONES, NODE COORDENATES
*           AND CONECTIVITIES      
***************************************************************
     
       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*20 MSSG
       DIMENSION LXPAREL(NUMEL,NPAREL),IXPARNP(NUMNP,NPARNP)
     ;          ,KXX(LMXNDL,NUMEL),X(NUMNP),Y(NUMNP),Z(NUMNP)
     ;          ,LTYPE(NUMEL),LNNDEL(NUMEL),CFPAREL(NUMEL,NPAREL)
     ;          ,CFPARNP(NUMNP,NPARNP)
       

C_______________________Writes element parameters zones
       WRITE(MAINF,*)
       WRITE(MAINF,*)'___________'//MSSG
       WRITE(MAINF,*)
       WRITE(MAINF,*)'***** LXPAREL VALUE *****'
       WRITE(MAINF,*)
       WRITE(MAINF,3000)
 3000  FORMAT(5X,'ELEM. Z O N E  N U M B E R S',/,5X,44('-'),//,
     ;          '  N.EL.  TRA   ARR   ART   STG   ',
     ;             'DSP   DFM   COE   POR   CRD   FOD'/)      
       DO L=1,NUMEL
         WRITE(MAINF,3100)L,LXPAREL(L,INTRA),
     ;                    LXPAREL(L,INARR),LXPAREL(L,INARRT),
     ;                    LXPAREL(L,INSTG),LXPAREL(L,INDSP),
     ;                    LXPAREL(L,INDFM),LXPAREL(L,INCOE),
     ;                    LXPAREL(L,INPOR),LXPAREL(L,INCRD),
     ;                    LXPAREL(L,INFOD)
       END DO

 3100  FORMAT(13I6)

C_______________________Writes nodal parameters zones

       WRITE(MAINF,*)  
       WRITE(MAINF,*)'***** IXPARNP VALUE *****'
       WRITE(MAINF,*)
       WRITE(MAINF,3300) 
 3300  FORMAT('  N_ND   CHP   CHPT   QQP   QQPT   ALF',
     ;   '   ALFT   CON   CONT   DMT')

       DO L=1,NUMNP
         WRITE(MAINF,3100)L,IXPARNP(L,INCHP),IXPARNP(L,INCHPT),
     ;                    IXPARNP(L,INQQP),IXPARNP(L,INQQPT),
     ;                    IXPARNP(L,INALF),IXPARNP(L,INALFT),
     ;                    IXPARNP(L,INCON),IXPARNP(L,INCONT),
     ;                    IXPARNP(L,INDMT)

       END DO



C_______________________Writes elements connectivities

       
       WRITE(MAINF,*)  
       WRITE(MAINF,*)'***** ELEMENTS CONECTIVITIES *****' 
       WRITE(MAINF,*)
       WRITE(MAINF,3400)
 3400  FORMAT('  N.EL. TYPE N. NOD. N.1  N.2  N.3  N.4  ',
     ;             'N.5  N.6  N.7  N.8   ')

       DO NE=1,NUMEL

         WRITE(MAINF,3500) NE,LTYPE(NE),LNNDEL(NE),
     ;              (KXX(JJ,NE),JJ=1,LNNDEL(NE)),
     ;              (0,JJ=LNNDEL(NE)+1,8)
       END DO

 3500  FORMAT(2I5,3X,9I5,1G11.4)

C_______________________Writes nodes coordinates

         WRITE(MAINF,*)  
         WRITE(MAINF,*)'***** NODE COORDINATES ******'
         WRITE(MAINF,*)
         WRITE(MAINF,3600)
 3600    FORMAT('    N         X         Y         Z      ')          


       DO I=1,NUMNP
         WRITE(MAINF,3700)I,X(I),Y(I),Z(I)
       END DO

 3700  FORMAT(I5,3(1X,F10.3))

C_______________________Writes element coeficients

       WRITE(MAINF,*)
       WRITE(MAINF,*)'***** CFPAREL VALUE *****'
       WRITE(MAINF,*)
       WRITE(MAINF,3800)
 3800  FORMAT(5X,'ELEM.  COEFICIENTS',/,5X,44('-'),//,
     ;          ' N_EL  TRA   ARR   ART   STG   ',
     ;             'DSP   DFM   COE   POR   CRD   FOD'/)
 3900  FORMAT(I5,10F6.3)      
       DO L=1,NUMEL
         WRITE(MAINF,3900)L,CFPAREL(L,INTRA),
     ;                    CFPAREL(L,INARR),CFPAREL(L,INARRT),
     ;                    CFPAREL(L,INSTG),CFPAREL(L,INDSP),
     ;                    CFPAREL(L,INDFM),CFPAREL(L,INCOE),
     ;                    CFPAREL(L,INPOR),CFPAREL(L,INCRD),
     ;                    CFPAREL(L,INFOD)
       END DO

C_______________________Writes nodal coeficients

       WRITE(MAINF,*)  
       WRITE(MAINF,*)'***** CFPARNP VALUE *****'
       WRITE(MAINF,*)
       WRITE(MAINF,4000) 
 4000  FORMAT(' N_ND   CHP  CHPT   QQP  QQPT   ALF',
     ;   '  ALFT   CON  CONT   DMT')
 4100  FORMAT(I5,9F6.3)

       DO L=1,NUMNP
         WRITE(MAINF,4100)L,CFPARNP(L,INCHP),CFPARNP(L,INCHPT),
     ;                    CFPARNP(L,INQQP),CFPARNP(L,INQQPT),
     ;                    CFPARNP(L,INALF),CFPARNP(L,INALFT),
     ;                    CFPARNP(L,INCON),CFPARNP(L,INCONT),
     ;                    CFPARNP(L,INDMT)

       END DO


       RETURN
       END
