      SUBROUTINE ENTDAT_MTDZ
     .(NUM_ZON,  IUNIN,  IVPAR,   PARC,INORPAR,NZONE_PAR,NTYPAR,
     .   NZPAR,  IXDMT,  NUMNP,   NPAR,IV_MTDZ,MXIV_MZ,RV_MTDZ,MXRV_MZ, 
     .   IOINV,   IERR,  MAINF ,NPBMX   ,NPARNP)

************************************************************************
*
*   PURPOSE
*     Reads data of the current zone number and initializes everything
*
*   DESCRIPTION
*     This subroutine is for test purposes, it is not meant to be understood,
*     if you understand anything, let me know
*
*   ARGUMENTS
*
*      NUM_ZON             Current matrix diffusion zone number
*      IUNIN               Input unit number
*      MAINF               Output unit number
*      IVPAR (NZPAR)       Global number of all parameters
*      PARC (NZPAR)        Current value of all parameters
*      INORPAR (NTYPAR)    Pointers to each type of parameters
*      NZONE_PAR(NTYPAR)    Number of zones of each type of parameters
*      NTYPAR              Number of types of parameters
*      NZPAR               Total number of parameters
*      NPAR                Number of parameters to be estimated
*      IXDMT (NUMNP)       MZ zone number of every node
*      NUMNP               Number of nodal points
*      IV_MTDZ             Vector of all MTDZ integer variables
*      MXIV_MZ             Dimension of IV_MTDZ
*      RV_MTDZ             Vector of all MTDZ real variables
*      MXRV_MZ             Dimension of RV_MTDZ
*      NPAR                Number of parameters of TRANSIN to be estimated
*      IOINV               Inverse problem option (prob solved if ioInv.gt.0)
*      IERR                Error counter
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER*4 NUM_ZON, NPAR, NUMNP, IV_MTDZ,
     .        MXIV_MZ, MXRV_MZ, I, IERR, IOINV,
     .        IUNIN, MAINF, NDPAR, INORPAR,NZONE_PAR,NTYPAR,IRNEXT,
     .        NZPAR, IVPAR, NPBMX, IXDMT   ,NPARNP

      REAL*8 RV_MTDZ,PARC

      DIMENSION IV_MTDZ(MXIV_MZ) , RV_MTDZ(MXRV_MZ) 

      INCLUDE 'EQUIV_MTDZ.FOR'

* .................................................... Error counter set to zero

      IERR=0

* ..................................................... Write input data heading

      IF (NUM_ZON.EQ.1 .AND. MAINF.GT.0) THEN
         WRITE(MAINF,1001)
 1001    FORMAT(//' MATRIX DIFFUSION ZONES'//
     .          '    ZONE   THICKNESS  DIF.COEF  RETARDAT.',
     .          '   LAMBDA  POROSITY    NUMBER  GEOMETRY'/
     .          '     NO.      Lm      ZONE NO.  ZONE NO.',
     .          '  ZONE NO.  ZONE NO.  OF TERMS    OPTION'/
     .          1x,7('-'),2x,10('-'),6(2x,8('-')))
      END IF

* ................................................................. Read MZ data

      READ (IUNIN,100) MTDZ_NUM, ESP_MZ, (MZ_ZON_NUM(I),I=1,4),
     .              NUM_TER, IOP_AN
  100 FORMAT(I5,F10.0,6I5)
      WRITE (MAINF,1002) MTDZ_NUM, ESP_MZ, (MZ_ZON_NUM(I),I=1,4),
     .              NUM_TER, IOP_AN
 1002 FORMAT(I8,E12.5,6I10)

* ............................................. Check consistency of zone number

      IF (NUM_ZON.NE.MTDZ_NUM) THEN
         IERR=IERR+1
         WRITE(MAINF,1003) NUM_ZON,MTDZ_NUM
 1003    FORMAT(' ERROR: When trying to read zone number:',i5,
     .          ', it has read:',I5)
      END IF

* ................. Initializes parameter values, MZ_PAR_NUM, MZ_EST and I_LM_MZ

      MZ_EST=0
      DO NDPAR=1,4

         CALL INI_PAR_MZ
     .(  NDPAR,  MAINF,  IVPAR,  PARC,INORPAR,NZONE_PAR, NTYPAR, NZPAR,
     .   IOINV,  IERR)

      END DO

* ................................................................. Set pointers

      II_NUD = MXIPRV+MXEQIR+1                                 ! Pointer to NUD_MZ
      NUM_NP_MZ = NUMNP                                    ! Temporary dimension

* ........................ Identify nodal points of the current MZ zone (NUD_MZ)

      CALL ID_NP_MZ ( IV_MTDZ(II_NUD), IXDMT, NUMNP ,NPBMX   ,NPARNP)

      IR_VOL = MXRPRV + 1
      IR_AN = IR_VOL + NUM_NP_MZ
      IR_EXP = IR_AN + 2*NUM_TER
      IR_DEXP = IR_EXP + 2*NUM_TER
      IRNEXT = IR_DEXP

      IF (MZ_EST.NE.0) THEN                     ! DEXP is needed, otherwise, not
         IRNEXT = IRNEXT + 6*NUM_TER
      END IF
         
      IR_RHS = IRNEXT
      IR_REC0 = IR_RHS + NUM_NP_MZ
      IRNEXT = IR_REC0

      IF (I_LM_MZ.NE.0) THEN                    ! REC0 is needed, otherwise, not
         IRNEXT = IRNEXT + NUM_NP_MZ
      END IF

      IR_RECN = IRNEXT
      IR_DREC0 = IR_RECN + NUM_NP_MZ*NUM_TER
      IRNEXT = IR_DREC0

      IF (I_LM_MZ.NE.0) THEN                   ! RV_MTDZ(IR_DREC0) is needed, otherwise, not
         IRNEXT = IRNEXT + NUM_NP_MZ*NPAR
      END IF

      IR_DRECN = IRNEXT

* ................................... Compute actual integer and real dimensions

      IF (IOINV.GT.0) THEN
        MXRV_MZ = IR_DRECN + NUM_NP_MZ*NUM_TER*NPAR
      ELSE
        MXRV_MZ = IR_DRECN
      END IF

      MXIV_MZ = MXIPRV + MXEQIR + NUM_NP_MZ

* .......................................Compute RV_MTDZ(IR_AN) (ALFn**2 and An)

      CALL COMP_AN_MZ ( RV_MTDZ(IR_AN) ,MAINF,IUNIN,IERR)

* .......................................Compute RV_MTDZ(IR_REC0) (SS memory function)

      IF (I_LM_MZ.NE.0) CALL COMP_G0_MZ (RV_MTDZ(IR_AN))

* ..................................... Save MZ variables back to IV and RV_MTDZ

      DO I=1,MXRPRV
         RV_MTDZ(I) = RPRV_MZ(I)
      END DO

      DO I=1,MXIPRV+MXEQIR
         IV_MTDZ(I) = IPRV_MZ(I)
      END DO

      RETURN
      END

      SUBROUTINE COMP_PRE_MTDZ
     .(VOLNOD, NUMNP,IV_MTDZ,MXIV_MZ,RV_MTDZ,MXRV_MZ)

************************************************************************
*
*   PURPOSE
*     Fill VOL_MZ (Volumes of MZ nodes) from VOLNOD
*
*
*   ARGUMENTS
*
*      NUM_ZON             Current matrix diffusion zone number
*      VOLNOD(NUMNP)       Nodal volumes
*      NUMNP               Number of nodal points
*      IV_MTDZ             Vector of all MTDZ integer variables
*      MXIV_MZ             Dimension of IV_MTDZ
*      RV_MTDZ             Vector of all MTDZ real variables
*      MXRV_MZ             Dimension of RV_MTDZ
*      IERR                Error counter
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER NUMNP, IV_MTDZ, MXIV_MZ, MXRV_MZ, I

      REAL*8 RV_MTDZ, VOLNOD

      DIMENSION IV_MTDZ(MXIV_MZ) , RV_MTDZ(MXRV_MZ) , VOLNOD(NUMNP)

      INCLUDE 'EQUIV_MTDZ.FOR'

* ....................................... Store IV and RV_MTDZ onto MZ variables

      DO I=1,MXRPRV
         RPRV_MZ(I) = RV_MTDZ(I) 
      END DO

      DO I=1,MXIPRV+MXEQIR
         IPRV_MZ(I) = IV_MTDZ(I)
      END DO

* ................................ Fill VOL_MZ (Volumes of MZ nodes) from VOLNOD


      CALL COMP_VOL_MZ
     .(  IV_MTDZ(II_NUD), RV_MTDZ(IR_VOL), VOLNOD, NUMNP)


* ..................................... Save MZ variables back to IV and RV_MTDZ

      DO I=1,MXRPRV
         RV_MTDZ(I) = RPRV_MZ(I)
      END DO

      DO I=1,MXIPRV+MXEQIR
         IV_MTDZ(I) = IPRV_MZ(I)
      END DO

      RETURN
      END

      SUBROUTINE UPD_PAR_MTDZ
     .(DELT, PARC, NZPAR,IV_MTDZ,MXIV_MZ,RV_MTDZ,MXRV_MZ)

************************************************************************
* PURPOSE
*    Updates variables as a consequence of a change in model parameters
*
*   ARGUMENTS
*
*      NUM_ZON             Current matrix diffusion zone number
*      DELT                Time increment
*      PARC (NZPAR)        Current value of all parameters
*      NZPAR               Total number of parameters
*      IV_MTDZ             Vector of all MTDZ integer variables
*      MXIV_MZ             Dimension of IV_MTDZ
*      RV_MTDZ             Vector of all MTDZ real variables
*      MXRV_MZ             Dimension of RV_MTDZ
*      IERR                Error counter
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER IV_MTDZ, MXIV_MZ, MXRV_MZ, I, NZPAR

      REAL*8 RV_MTDZ, DELT, PARC

      DIMENSION IV_MTDZ(MXIV_MZ) , RV_MTDZ(MXRV_MZ)

      INCLUDE 'EQUIV_MTDZ.FOR'

* ....................................... Store IV and RV_MTDZ onto MZ variables

      DO I=1,MXRPRV
         RPRV_MZ(I) = RV_MTDZ(I) 
      END DO

      DO I=1,MXIPRV+MXEQIR
         IPRV_MZ(I) = IV_MTDZ(I)
      END DO

* ........... Updates variables as a consequence of a change in model parameters

      IF (MZ_EST.NE.0) THEN
          CALL UPD_PAR_MZ 
     .(  RV_MTDZ(IR_AN), RV_MTDZ(IR_EXP), RV_MTDZ(IR_DEXP), 
     .    DELT, PARC,  NZPAR)
      END IF

* ..................................... Save MZ variables back to IV and RV_MTDZ

      DO I=1,MXRPRV
         RV_MTDZ(I) = RPRV_MZ(I)
      END DO

      DO I=1,MXIPRV+MXEQIR
         IV_MTDZ(I) = IPRV_MZ(I)
      END DO

      RETURN
      END

      SUBROUTINE UPD_DT_MTDZ 
     .(DT_NEW, IV_MTDZ,MXIV_MZ,RV_MTDZ,MXRV_MZ)

************************************************************************
*
*   PURPOSE
*     Updates dimensionless time and related variables in 
*     response to a change in time increment
*
*
*   ARGUMENTS
*
*      NUM_ZON             Current matrix diffusion zone number
*      DT_NEW              New time increment
*      IV_MTDZ             Vector of all MTDZ integer variables
*      MXIV_MZ             Dimension of IV_MTDZ
*      RV_MTDZ             Vector of all MTDZ real variables
*      MXRV_MZ             Dimension of RV_MTDZ
*      IERR                Error counter
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER IV_MTDZ, MXIV_MZ, MXRV_MZ, I

      REAL*8 RV_MTDZ, DT_NEW

      DIMENSION IV_MTDZ(MXIV_MZ) , RV_MTDZ(MXRV_MZ) 

      INCLUDE 'EQUIV_MTDZ.FOR'

* ....................................... Store IV and RV_MTDZ onto MZ variables

      DO I=1,MXRPRV
         IF (I.NE.12 .OR. ICH_IT_MAR.EQ.0) RPRV_MZ(I) = RV_MTDZ(I)
      END DO
      ICH_IT_MAR=0           ! Now we are inside an inverse problem iteration

      DO I=1,MXIPRV+MXEQIR
         IPRV_MZ(I) = IV_MTDZ(I)
      END DO

* ............. Updates time variables in response to a change in time increment

      CALL UPD_DT_MZ 
     .( RV_MTDZ(IR_AN), RV_MTDZ(IR_EXP), RV_MTDZ(IR_DEXP), DT_NEW)

* ..................................... Save MZ variables back to IV and RV_MTDZ

      DO I=1,MXRPRV
         RV_MTDZ(I) = RPRV_MZ(I)
      END DO

      DO I=1,MXIPRV+MXEQIR
         IV_MTDZ(I) = IPRV_MZ(I)
      END DO

      RETURN
      END

      SUBROUTINE COMP_CNTR_MTDZ
     .(NUM_ZON, BTRA, DTRA, NUMNP,
     . IV_MTDZ, MXIV_MZ, RV_MTDZ, MXRV_MZ, IXDMT)

************************************************************************
* PURPOSE
*
*        Adds contribution of the current matrix diffusion zone(MZ) to
*        BTRA and DTRA
*
*   ARGUMENTS
*
*      NUM_ZON             Current matrix diffusion zone number
*      BTRA (NUMNP)        R.H.S. of transport equations
*      DTRA (NUMNP)        Storage matrix of transport equations
*      NUMNP               Number of nodal points
*      IV_MTDZ             Vector of all MTDZ integer variables
*      MXIV_MZ             Dimension of IV_MTDZ
*      RV_MTDZ             Vector of all MTDZ real variables
*      MXRV_MZ             Dimension of RV_MTDZ
*      IERR                Error counter
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER NUM_ZON, NUMNP, IV_MTDZ, MXIV_MZ, MXRV_MZ, I,IXDMT

      REAL*8 RV_MTDZ, BTRA, DTRA

      DIMENSION IV_MTDZ(MXIV_MZ) , RV_MTDZ(MXRV_MZ) ,
     .          BTRA(NUMNP), DTRA(NUMNP)

      INCLUDE 'EQUIV_MTDZ.FOR'

* ....................................... Store IV and RV_MTDZ onto MZ variables

      DO I=1,MXRPRV
         RPRV_MZ(I) = RV_MTDZ(I) 
      END DO

      DO I=1,MXIPRV+MXEQIR
         IPRV_MZ(I) = IV_MTDZ(I)
      END DO

* ...................................... Adds contribution of the current matrix
* .......................................... diffusion zone(MZ) to BTRA and DTRA

      CALL COMP_CNTR_MZ
     .( IV_MTDZ(II_NUD), RV_MTDZ(IR_EXP),  RV_MTDZ(IR_RHS), 
     .  RV_MTDZ(IR_VOL), RV_MTDZ(IR_REC0),   RV_MTDZ(IR_RECN),
     .   BTRA, DTRA, NUMNP,IXDMT,NUM_ZON)

* ..................................... Save MZ variables back to IV and RV_MTDZ

      DO I=1,MXRPRV
         RV_MTDZ(I) = RPRV_MZ(I)
      END DO

      DO I=1,MXIPRV+MXEQIR
         IV_MTDZ(I) = IPRV_MZ(I)
      END DO

      RETURN
      END

      SUBROUTINE COMP_DER_MTDZ
     .(NUM_ZON,    CCAL,CCALAN,   DERC,NUMNP, NPAR,
     . IV_MTDZ, MXIV_MZ, RV_MTDZ, MXRV_MZ, IXDMT)

************************************************************************
* PURPOSE
*    Computes contributions of MZ to DERC (RHS of sensitivity equations)
*
*   ARGUMENTS
*
*      NUM_ZON             Current matrix diffusion zone number
*      CCAL (NUMNP)        C[k+1,i]
*      CCALAN (NUMNP)     C[k,i]
*      DERC (NUMNP,NPAR)   Derivatives of C[k+1,i] with respect to param
*      DERCAN (NUMNP,NPAR)Derivatives of C[k,i] with respect to param
*      NUMNP               Number of nodal points
*      IV_MTDZ             Vector of all MTDZ integer variables
*      MXIV_MZ             Dimension of IV_MTDZ
*      RV_MTDZ             Vector of all MTDZ real variables
*      MXRV_MZ             Dimension of RV_MTDZ
*      IERR                Error counter
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER NUM_ZON, NUMNP, IV_MTDZ, MXIV_MZ, MXRV_MZ, I, NPAR
     ;                        ,IXDMT

      REAL*8 RV_MTDZ, CCAL, CCALAN,  DERC

      DIMENSION IV_MTDZ(MXIV_MZ) , RV_MTDZ(MXRV_MZ)

      INCLUDE 'EQUIV_MTDZ.FOR'

* ....................................... Store IV and RV_MTDZ onto MZ variables

      DO I=1,MXRPRV
         RPRV_MZ(I) = RV_MTDZ(I) 
      END DO

      DO I=1,MXIPRV+MXEQIR
         IPRV_MZ(I) = IV_MTDZ(I)
      END DO

* ................... Contributions of MZ to DERC (RHS of sensitivity equations)

      CALL COMP_DER_MZ
     .( IV_MTDZ(II_NUD), RV_MTDZ(IR_VOL),
     .  RV_MTDZ(IR_EXP), RV_MTDZ(IR_DEXP), RV_MTDZ(IR_RHS),   
     .  RV_MTDZ(IR_RECN), RV_MTDZ(IR_DREC0),
     .  RV_MTDZ(IR_DRECN),   
     . CCAL,CCALAN,   DERC, NUMNP,   NPAR,IXDMT,NUM_ZON)

* ..................................... Save MZ variables back to IV and RV_MTDZ

      DO I=1,MXRPRV
         RV_MTDZ(I) = RPRV_MZ(I)
      END DO

      DO I=1,MXIPRV+MXEQIR
         IV_MTDZ(I) = IPRV_MZ(I)
      END DO

      RETURN
      END

      SUBROUTINE UPD_DREC_MTDZ
     .(NUM_ZON, CCAL,
     . CCALAN,  NUMNP,   DERC,DERCAN,   NPAR,
     . IV_MTDZ, MXIV_MZ, RV_MTDZ, MXRV_MZ, IXDMT)

************************************************************************
* PURPOSE
*   Updates derivatives of vectors RV_MTDZ(IR_REC0) (Io) and RV_MTDZ(IR_RECN) (In)
*
*   ARGUMENTS
*
*      NUM_ZON             Current matrix diffusion zone number
*      CCAL (NUMNP)              C[k+1,i]
*      CCALAN (NUMNP)           C[k,i]
*      DERC (NUMNP,NPAR)         Derivatives of C[k+1,i] w.r.t. parameters
*      DERCAN (NUMNP,NPAR)      Derivatives of C[k,i] w.r.t. parameters
*      NUMNP                     Number of nodal points
*      NPAR                      Number of model parameters
*      IV_MTDZ             Vector of all MTDZ integer variables
*      MXIV_MZ             Dimension of IV_MTDZ
*      RV_MTDZ             Vector of all MTDZ real variables
*      MXRV_MZ             Dimension of RV_MTDZ
*      IERR                Error counter
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER NUM_ZON, NUMNP, IV_MTDZ, MXIV_MZ, MXRV_MZ, I,NPAR
     ;                                  ,IXDMT

      REAL*8 RV_MTDZ, CCAL, CCALAN,  DERC,DERCAN

      DIMENSION IV_MTDZ(MXIV_MZ) , RV_MTDZ(MXRV_MZ)

      INCLUDE 'EQUIV_MTDZ.FOR'

* ....................................... Store IV and RV_MTDZ onto MZ variables

      DO I=1,MXRPRV
         RPRV_MZ(I) = RV_MTDZ(I) 
      END DO

      DO I=1,MXIPRV+MXEQIR
         IPRV_MZ(I) = IV_MTDZ(I)
      END DO

* ....................... Updates derivatives of vectors RV_MTDZ(IR_REC0) (Io) and RV_MTDZ(IR_RECN) (In)

      CALL UPD_DREC_MZ
     .( IV_MTDZ(II_NUD),  RV_MTDZ(IR_EXP),  RV_MTDZ(IR_DEXP),  
     .  RV_MTDZ(IR_RECN), RV_MTDZ(IR_DREC0),
     .  RV_MTDZ(IR_DRECN),   CCAL,
     . CCALAN,  NUMNP,   DERC,DERCAN,   NPAR,IXDMT,NUM_ZON)

* ..................................... Save MZ variables back to IV and RV_MTDZ

      DO I=1,MXRPRV
         RV_MTDZ(I) = RPRV_MZ(I)
      END DO

      DO I=1,MXIPRV+MXEQIR
         IV_MTDZ(I) = IPRV_MZ(I)
      END DO

      RETURN
      END

      SUBROUTINE UPD_REC_MTDZ
     .(NUM_ZON, CCAL,CCALAN,  NUMNP,
     . IV_MTDZ, MXIV_MZ, RV_MTDZ, MXRV_MZ, IXDMT)

************************************************************************
* PURPOSE
*    Updates vectors RV_MTDZ(IR_REC0) (Io) and RV_MTDZ(IR_RECN) (In)
*
* ARGUMENTS
*
*      NUM_ZON             Current matrix diffusion zone number
*      CCAL (NUMNP)        C[k+1,i]
*      CCALAN (NUMNP)     C[k,i]
*      NUMNP               Number of nodal points
*      IV_MTDZ             Vector of all MTDZ integer variables
*      MXIV_MZ             Dimension of IV_MTDZ
*      RV_MTDZ             Vector of all MTDZ real variables
*      MXRV_MZ             Dimension of RV_MTDZ
*      IERR                Error counter
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER NUM_ZON, NUMNP, IV_MTDZ, MXIV_MZ, MXRV_MZ, I,IXDMT

      REAL*8 RV_MTDZ, CCAL, CCALAN

      DIMENSION IV_MTDZ(MXIV_MZ) , RV_MTDZ(MXRV_MZ) ,
     .          CCAL(NUMNP), CCALAN(NUMNP)

      INCLUDE 'EQUIV_MTDZ.FOR'

* ....................................... Store IV and RV_MTDZ onto MZ variables

      DO I=1,MXRPRV
         RPRV_MZ(I) = RV_MTDZ(I) 
      END DO

      DO I=1,MXIPRV+MXEQIR
         IPRV_MZ(I) = IV_MTDZ(I)
      END DO

* ...................................... Updates vectors RV_MTDZ(IR_REC0) (Io) and RV_MTDZ(IR_RECN) (In)

       CALL UPD_REC_MZ
     .( IV_MTDZ(II_NUD), RV_MTDZ(IR_EXP),   RV_MTDZ(IR_REC0),   
     .  RV_MTDZ(IR_RECN),   CCAL,CCALAN,  NUMNP,IXDMT,NUM_ZON)

* ..................................... Save MZ variables back to IV and RV_MTDZ

      DO I=1,MXRPRV
         RV_MTDZ(I) = RPRV_MZ(I)
      END DO

      DO I=1,MXIPRV+MXEQIR
         IV_MTDZ(I) = IPRV_MZ(I)
      END DO

      RETURN
      END

