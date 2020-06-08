      SUBROUTINE DERQ_NUD
     &          (BETAC    ,CAUDAL   ,CAUX1    ,CCALAN   ,CCALIT
     &          ,CFPARNP  ,CREF     ,DENSREF  ,DERC     ,DTIM
     &          ,EPSFLU   ,FNT      ,HAUX1    ,IBCOD    ,IBTCO
     &          ,IDIMDERC ,IDIMFNT  ,IFLAGS   ,INDSSTR  ,INEW
     &          ,INORPAR  ,INTI     ,IOCAP    ,IODENS   ,IOFLLI
     &          ,IOFMLF   ,IVPAR    ,IXPARNP  ,KXX      ,LMXNDL
     &          ,NFLAGS   ,NFNL     ,NFNLPAR  ,NFNLPRG  ,NFNLTIP
     &          ,NFTPAR   ,NINT     ,NPAR     ,NPARNP   ,NPPNP
     &          ,NPZON    ,NTYPAR   ,NUMEL    ,NUMNP    ,NZONE_PAR
     &          ,NZPAR    ,PARACD   ,PARC     ,PARNP    ,IDIMWGT
     &          ,WGT_PAR  ,IPNT_PAR ,IPOS     ,DERIV)

********************************************************************************
*
* PURPOSE
*
*   Computes the derivative of nodal flow w.r.t. nodal parameters (exp. dep.)
*
* DESCRIPTION
*
*   Computes the derivative of nodal flow w.r.t. nodal parameters, prescribed 
*   head, prescribed flow and leakage (explicit dependence)
*
* EXTERNAL VARIABLES: ARRAYS
*
*  CAUX1                  Array containing concentrations, ponderated by THETAT 
*                         time factor                                           
*  DERC                   Nodal concentration derivatives with respect to       
*                         estimated parameters.                                 
*  HAUX1                  Array containing heads, ponderated by THETAF time     
*                         weight. Is equal to THETAF*HCAL+(1-THETAF)*HCALAN     
*  IBCOD                  Flow boundary condition index                         
*  IFLAGS                 Array with different writing options. Used mainly for
*                         debugging.
*  INORPAR                Array containing the indexes with the location
*                         of the different paramaters in arrays PARC, PARM,
*                         IVPAR, NFTPAR, STPAR and FNTPAR
*  IVPAR                  Vector containing estimation index for all
*                         parameters
*  IXPARNP                Array containing zone number of each node j,
*                         corresponding to INpar index parameter zone.
*  PARC                   Vector containing calculated values for all
*                         parameters
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
*  INDSSTR                If 1, transient flow. If 0, steady state flow
*  NFLAGS                 Maximum number of allowed flags                       
*  NPAR                   Total number of parameters to be estimated            
*  NPARNP                 Number of nodal parameters in current problem         
*  NTYPAR                 Number of different type parameters including         
*                         anisotropy (used for dimensioning)                    
*  NUMNP                  Number of nodes                                       
*  NZPAR                  Total number of zones for all nodal and element       
*                         parameters including each transmissivity tensor       
*                         component.                                            
*
* INTERNAL VARIABLES: SCALARS
*
*  DERQ                   Derivative of nodal flow w.r.t. current parameter
*
*
* HISTORY
*
*     AMS      3-2002     First coding (starting from TRANSIN-II)
*
********************************************************************************

      IMPLICIT NONE

C-------------------- External

      INTEGER*4::IDIMDERC ,IDIMFNT  ,INDSSTR  ,INEW     ,INTI
     &          ,IOCAP    ,IODENS   ,IOFLLI   ,IOFMLF   ,LMXNDL
     &          ,NFLAGS   ,NFNL     ,NINT     ,NPAR     ,NPARNP
     &          ,NPPNP    ,NPZON    ,NTYPAR   ,NUMEL    ,NUMNP
     &          ,NZPAR    ,IDIMWGT  ,IPNT_PAR

      
      REAL*8::DENSREF  ,BETAC    ,CREF     ,DTIM     ,EPSFLU

      REAL*8::DENS
            
      INTEGER*4::IBCOD(NUMNP)     ,IBTCO(NUMNP)         ,INORPAR(NTYPAR)
     &          ,IVPAR(NZPAR)     ,IXPARNP(NUMNP,NPARNP),IFLAGS(NFLAGS)
     &          ,KXX(LMXNDL,NUMEL),NFNLPAR(NZPAR)       ,NFNLPRG(8,NFNL)
     &          ,NFNLTIP(NFNL)    ,NFTPAR(NZPAR)      ,NZONE_PAR(NTYPAR)

      REAL*8::CAUDAL(NUMNP)          ,CAUX1(NUMNP)
     &       ,CCALAN(NUMNP)          ,CCALIT(NUMNP)
     &       ,CFPARNP(NUMNP,NPARNP)  ,DERC(NUMNP,NPAR,IDIMDERC) 
     &       ,FNT(IDIMFNT,NINT)      ,HAUX1(NUMNP)
     &       ,PARC(NZPAR)            ,PARACD(3,NFNL)
     &       ,PARNP(NUMNP,NPPNP)     ,WGT_PAR(IDIMWGT)



C-------------------- Internal

      INTEGER*4::I      ,IP     ,IPAR   ,JJ     ,L     ,NCNF   ,NNUD
     &          ,NPTOT  ,NZ

      REAL*8::CEXT     ,CEXT_CINT,CINT     ,DENSEXT  ,HEAD     ,LEAKCF
     &       ,VAUX

      INTEGER*4::INDEX(12)   ,IPOS(NPAR)

      REAL*8::CFPARAM(12)    ,DERIV(NPAR)  ,XPARAM(8)

C-------------------- First executable statement

      DO I=1,NUMNP

          IF (CAUDAL(I).GT.0
     &       .AND. (IBTCO(I).EQ.2 .OR. IBTCO(I).EQ.3) ) THEN

              CEXT = PARNP(I,4)
              CINT = CAUX1(I)

              IF (IODENS.EQ.1) THEN


                  DENSEXT = DENS(DENSREF,BETAC,CEXT,CREF)

                  CEXT_CINT = DENSEXT*(CEXT - CINT)

              ELSE

                  CEXT_CINT = CEXT - CINT

              END IF !IODENS.EQ.1

              IF (IBCOD(I).EQ.2 .OR. IBCOD(I).EQ.4) THEN

                  NZ = IXPARNP(I,3+INDSSTR)
                  JJ = INORPAR(10)+NZ
                  IP = IVPAR(JJ)

                  IF (IP.GT.0 .OR. IOFLLI.NE.0) THEN


                      IF (IOFLLI.NE.0) THEN

                          NCNF = NFNLPAR(JJ)

                      ELSE

                          NCNF = 0

                      END IF !IOFLLI.NE.0

                      INDEX(1) = JJ

                      CFPARAM(1) = CFPARNP(I,3+INDSSTR)

                      CALL DER_PARAM
     &    (CFPARAM  ,DTIM     ,EPSFLU   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,0        ,IOFMLF   ,IP       ,L
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,NPZON    ,NUMEL
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     &    ,PARACD   ,PARC(INORPAR(19)+1),CCALIT   ,CCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)



                      IF (NPTOT.GT.0) THEN

C---------------------------DERQ = PARNP(I,2)/PARC(INORPAR(10)+NZ)

                          DO IPAR=1,NPTOT

                              DERC(I,IPOS(IPAR),INEW) = 
     &                                           DERC(I,IPOS(IPAR),INEW) 
     &                                         + DERIV(IPAR)*CEXT_CINT


                              IF (IFLAGS(30).EQ.1) THEN
                                  WRITE(77,10) ' >',I,IP,DERIV(IPAR)
   10                             FORMAT (A2,2I6,E20.10)
                              END IF !IFLAGS(30).EQ.1

                          END DO !IPAR=1,NPTOT

                          

                      END IF !NPTOT.GT.0

                  END IF !IP.GT.0 .OR. IOFLLI.NE.0

              ELSE IF (IBCOD(I).GE.3) THEN 

C------------------------- Leakage

C------------------------- Leakage coeficient

                  NZ = IXPARNP(I,5+INDSSTR)
                  JJ = INORPAR(11)+NZ
                  IP = IVPAR(JJ)

                  IF (IP.GT.0 .OR. IOFLLI.NE.0) THEN


                      IF (IOFLLI.NE.0) THEN

                          NCNF = NFNLPAR(JJ)

                      ELSE

                          NCNF = 0

                      END IF !IOFLLI.NE.0

                      INDEX(1) = JJ

                      CFPARAM(1) = CFPARNP(I,5+INDSSTR)

                      CALL DER_PARAM
     &    (CFPARAM  ,DTIM     ,EPSFLU   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,0        ,IOFMLF   ,IP       ,L
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,NPZON    ,NUMEL
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     &    ,PARACD   ,PARC(INORPAR(19)+1),CCALIT   ,CCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)
                      

                      IF (NPTOT.GT.0) THEN

                          HEAD = PARNP(I,1)

                          VAUX = (HEAD-HAUX1(I))*CEXT_CINT

c-------------------------DERQ = PARNP(I,3)/PARC(INORPAR(11)+NZ)*(PARNP(I,1)-HAUX1(I))

                          DO IPAR=1,NPTOT

                              DERC(I,IPOS(IPAR),INEW) = 
     &                                           DERC(I,IPOS(IPAR),INEW)
     &                                           + DERIV(IPAR)*VAUX

                              IF (IFLAGS(30).EQ.1) THEN
                                  WRITE(77,20) I,IP,DERIV(IPAR)
   20                             FORMAT (A2,2I6,E20.10)
                              END IF !IFLAGS(30).EQ.1

                          END DO !IPAR=1,NPTOT

                      END IF !NPTOT.GT.0

                  END IF !IP.GT.0 .OR. IOFLLI.NE.0

C------------------------- Leakage prescribed head

                  NZ = IXPARNP(I,1+INDSSTR)
                  JJ = INORPAR(9)+NZ
                  IP = IVPAR(JJ)

                  IF (IP.GT.0 .OR. IOFLLI.NE.0) THEN


                      IF (IOFLLI.NE.0) THEN

                          NCNF = NFNLPAR(JJ)

                      ELSE

                          NCNF = 0

                      END IF !IOFLLI.NE.0

                      INDEX(1) = JJ

                      CFPARAM(1) = CFPARNP(I,1+INDSSTR)

                      CALL DER_PARAM
     &    (CFPARAM  ,DTIM     ,EPSFLU   ,IDIMFNT  ,INDSSTR  ,INORPAR(19)
     &    ,INTI     ,IOCAP    ,0        ,IOFMLF   ,IP       ,L
     &    ,LMXNDL   ,NFLAGS   ,NFNL     ,NCNF     ,NFNLTIP
     &    ,NFTPAR   ,NINT     ,NNUD     ,NPTOT    ,NPZON    ,NUMEL
     &    ,NUMNP    ,NZPAR    ,NZONE_PAR(14),PARC(JJ),DERIV ,FNT
     &    ,IFLAGS   ,INDEX    ,IPOS     ,IVPAR    ,KXX      ,NFNLPRG
     &    ,PARACD   ,PARC(INORPAR(19)+1),CCALIT   ,CCALAN   ,XPARAM
     ;    ,IDIMWGT  ,NPAR     ,IPNT_PAR ,WGT_PAR)

                      IF (NPTOT.GT.0) THEN

                          LEAKCF = PARNP(I,3)
                          VAUX = LEAKCF*CEXT_CINT

c-------------------------DERQ = PARNP(I,1)/PARC(INORPAR(9)+NZ)*PARNP(I,3)

                          DO IPAR=1,NPTOT

                              DERC(I,IPOS(IPAR),INEW) =
     &                                           DERC(I,IPOS(IPAR),INEW)
     &                                           + DERIV(IPAR)*VAUX
                              IF (IFLAGS(30).EQ.1) THEN
                                  WRITE(77,30) I,IP,DERIV(IPAR)
   30                             FORMAT(2I6,E20.10)
                              END IF !IFLAGS(30).EQ.1

                          END DO !IPAR=1,NPTOT

                      END IF !NPTOT.GT.0

                  END IF !IP.GT.0 .OR. IOFLLI.NE.0

              ENDIF  ! IBCOD(I).GE.3

          ENDIF ! CAUDAL(I)>0 .AND. IBTCO=2,3

      END DO !I=1,NUMNP

      END SUBROUTINE DERQ_NUD
