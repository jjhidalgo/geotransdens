      SUBROUTINE ENTDAT_DMT
     &          (NZDMT    ,IVPAR    ,PARC     ,INORPAR  ,NZONE_PAR
     &          ,NTYPAR   ,NZPAR    ,IXDMT    ,NUMNP    ,NPAR
     &          ,INPWR    ,MAINFI   ,IOINV    ,NPBMX    ,NPARNP)

************************************************************************
*
*   PURPOSE
*     Reads matrix diffusion data 
*
*   DESCRIPTION
*     This subroutine performs the following tasks:
*       1. Read and initialize every matrix diffusion zone
*       2. Reserve space for matrix diffusion
*
*   ARGUMENTS
*      NZDMT               Number of matrix diffusion zones
*      NTDMT               Max. number of terms for matrix diffusion expansion
*      IVPAR (NZPAR)       Global number of all parameters
*      PARC (NZPAR)        Current value of all parameters
*      INORPAR (NTYPAR)    Pointers to each type of parameters
*      NZONE_PAR(NTYPAR)    Number of zones of each type of parameters
*      NTYPAR              Number of types of parameters
*      NZPAR               Total number of parameters
*      NPAR                Number of parameters to be estimated
*      IXDMT (NUMNP)       MZ zone number of every node
*      NUMNP               Number of nodal points
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'MAIN_COMM.FOR'
*      IV             Vector of integer variables
*      LASTII         Last position of IV that has been reserved up to now
*      IIMAX          Dimension of IV
*      RV             Vector real variables
*      LASTII         Last position of IR that has been reserved up to now
*      IRMAX          Dimension of IR

      INCLUDE 'COMMON_DMT.FOR'
*      DELT_LST       Last DELTA_T 
*      IN_IV_DMT      Beginning of matrix diffusion data in IV (pointer)
*      IN_RV_DMT      Beginning of matrix diffusion data in IR (pointer)
*      IUNIN,  MAINF,  ?
*      I_ST_DMT       Index of mat. diff. computation status
*                        0- Nothing has been done yet
*                        1- Data has already been read (ENTDAT_MTDZ)
*                        2- Preliminary calculations done (COMP_PRE_MTDZ)
*                        3- MTDZ variables have been updated in response to 
*                           change in parameters (UPD_PAR_MTDZ)
*                        4- MTDZ variables have been updated in response to 
*                           change in time increment(UPD_DT_MTDZ)
*                        5- Contribution to transport included (COMP_CNTR_MTDZ)
*                        6- Contribution to derivatives included (COMP_DER_MTDZ)
*                        7- Derivatives of Io and In updated (UPD_DREC_MTDZ)
*                        8- Io and In updated (UPD_REC_MTDZ)
      SAVE /DMT/

      INTEGER NZDMT, NUM_ZON, NPAR, NUMNP,
     .        MXIV_MZ, MXRV_MZ, IERR, INPWR, MAINFI,
     .        INORPAR,NZONE_PAR,NTYPAR, IOINV,
     .        NZPAR, IVPAR, IST, IN_IV_MZ, IN_RV_MZ,NPBMX
     ;       , IXDMT   ,NPARNP

      REAL*8 PARC

* .............................................................. IO unit numbers

      IUNIN = 12              ! Input unit number
      MAINF = MAINFI

* .............. Set pointers to beginning of matrix diffusion data in IV and RV

      IN_IV_DMT = LASTII + 1 !Beginning of matrix diffusion data in IV (pointer)
      LASTII = LASTII + 4*NZDMT                ! Reserve four integers per zone 
                             ! (zone pointers to IV and RV and their dimensions)

      IN_RV_DMT = LASTIR + 1 !Beginning of matrix diffusion data in IR (pointer)
*     LASTIR = LASTIR + 0                               ! No real space reserved


* ..................................... Loop to read every matrix diffusion zone

      NUM_MTDZ= NZDMT
      IST=0

      DO NUM_ZON = 1,NZDMT

         IN_IV_MZ = LASTII + 1       ! Beginning of integers at the current zone
         IN_RV_MZ = LASTIR + 1       ! Beginning of reals at the current zone

         IV (IN_IV_DMT + 4*NUM_ZON-4) = IN_IV_MZ ! store the pointer to zone beg
         IV (IN_IV_DMT + 4*NUM_ZON-3) = IN_RV_MZ

         MXIV_MZ = IIMAX - LASTII  ! Max. dimension for integers in current zone
         MXRV_MZ = IRMAX - LASTIR     ! Max. dimension for reals in current zone

         CALL ENTDAT_MTDZ
     .(NUM_ZON,  IUNIN,  IVPAR,  PARC,INORPAR,NZONE_PAR,NTYPAR,
     .   NZPAR,  IXDMT,  NUMNP,   NPAR,   IV(IN_IV_MZ), MXIV_MZ,
     . RV(IN_RV_MZ),   MXRV_MZ,IOINV,   IERR, MAINF ,NPBMX   ,NPARNP)

         IF (IERR.GT.0) THEN           ! Errors when reading mat. diffusion zone

            IST=1
            WRITE (MAINF,1001) IERR,NUM_ZON
 1001       FORMAT (' ERROR:',i5,' ERRORS when reading matrix ',
     .                  'diffusion zone number:',i3)
         END IF

* ..................................................... Update LASTII and LASTIR

         LASTII = LASTII + MXIV_MZ         ! update LASTII with actual int. dim.
         IV (IN_IV_DMT + 4*NUM_ZON-2) = MXIV_MZ

         IF ( LASTII + NUMNP .GT. IIMAX ) THEN                   ! IV size exceeded
            IST=2
            WRITE (MAINF,1010) IIMAX, LASTII, NUM_ZON, LASTII
 1010       FORMAT (' ERROR: the max size of IV is,'I10/
     .              '        it has reached ',I10,' after reading',
     .              ' matrix zone :',i5/
     .              ' You have to increase IIMAX in MAIN_COMMON.FOR',
     .              ' to, at least,',I10)
         END IF

         LASTIR = LASTIR + MXRV_MZ         ! update LASTIR with actual int. dim.
         IV (IN_IV_DMT + 4*NUM_ZON-1) = MXRV_MZ

         IF ( LASTIR +NUMNP .GT. IRMAX ) THEN                    ! RV size exceeded
            IST=2
            WRITE (MAINF,1020) IRMAX, LASTIR, NUM_ZON, LASTIR
 1020       FORMAT (' ERROR: the max size of RV is,'I10/
     .              '        it has reached ',I10,' after reading',
     .              ' matrix zone :',i5/
     .              ' You have to increase IRMAX in MAIN_COMMON.FOR',
     .              ' to, at least,',I10)
         END IF

      END DO ! NUM_ZON

      I_DTRA=LASTIR+1            ! Position of auxiliar array DTRA_DMT
      I_NELBYNUD=LASTII+1        ! Position of auxiliar array NELBYNUD 
                                 ! (number of elem. surrounding each node)
      I_CMP_NELBNUD=0

* ............................................... Stop, if too many errors found

      IF (IST.GT.0) THEN

         WRITE(MAINF,1030)
 1030    FORMAT(///' Too many ERRORS when trying to read ',
     .          'matrix diffusion data'/' SORRY, BUT I STOP')
         STOP

      END IF

* ...................................... Matrix diffusion data read successfully

      I_ST_DMT = 1

      RETURN
      END


      SUBROUTINE COMP_PRE_DMT
     . ( VOLNOD, NUMNP)

************************************************************************
*
*   PURPOSE
*     Preliminary computations
*
*   DESCRIPTION
*     Basically, all this subroutine does is to store volumes associated to 
*     node in every matrix zone
*
*   ARGUMENTS
*      VOLNOD (NUMNP)      Volume of every node
*      NUMNP               Number of nodal points
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'MAIN_COMM.FOR'

      INCLUDE 'COMMON_DMT.FOR'

      INTEGER NUM_ZON, NUMNP, 
     .        MXIV_MZ, MXRV_MZ, IERR, IST,
     .        IN_IV_MZ, IN_RV_MZ

      REAL*8 VOLNOD

* .............................................................. Initialization

      IF (I_ST_DMT.LT.1) THEN

         PRINT 1000
 1000    FORMAT (' ERROR: COMP_PRE_DMT reached ',
     .           ' before reading matrix diffusion data')
         STOP

      END IF

      IST = 0

* ............................................ Loop over matrix diffusion zones

      DO NUM_ZON = 1,NUM_MTDZ

         IERR=0

         IN_IV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-4) 
         IN_RV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-3) 
         MXIV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-2)
         MXRV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-1)


         CALL COMP_PRE_MTDZ
     .(VOLNOD, NUMNP,IV(IN_IV_MZ),MXIV_MZ,RV(IN_RV_MZ),MXRV_MZ)


         IF (IERR.GT.0) THEN           ! Errors when reading mat. diffusion zone

            IST=1
            WRITE(MAINF,1001) IERR,NUM_ZON
 1001       FORMAT (' ERROR:',i5,' ERRORS during',
     .     ' prelim calculations in matrix diffusion zone number:',i3)

         END IF

      END DO ! NUM_ZON

* ............................................... Stop, if too many errors found

      IF (IST.GT.0) THEN

         WRITE(MAINF,1030)
 1030    FORMAT(///' Too many ERRORS during',
     .     ' prelim calculations'/' SORRY, BUT I STOP')
         STOP

      END IF

* ......................................... Preliminary calculations successful

      I_ST_DMT = 2

      RETURN
      END


      SUBROUTINE UPD_PAR_DMT
     .( PARC, NZPAR)

************************************************************************
* PURPOSE
*    Updates variables as a consequence of a change in model parameters
*
*   ARGUMENTS
*
*      PARC (NZPAR)        Current value of all parameters
*      NZPAR               Total number of parameters
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'MAIN_COMM.FOR'

      INCLUDE 'COMMON_DMT.FOR'

      INTEGER NUM_ZON,  NZPAR,
     .        MXIV_MZ, MXRV_MZ, IERR, IST,
     .        IN_IV_MZ, IN_RV_MZ

      REAL*8 DELT, PARC(NZPAR)

* .............................................................. Initialization

      IF (I_ST_DMT.LT.2) THEN

         PRINT 1000
 1000    FORMAT (' ERROR: UPD_PAR_DMT reached ',
     .           ' before preliminary calculations')
         STOP

      ELSE IF (I_ST_DMT .EQ. 2) THEN

         DELT = 1D0

      ELSE

         DELT = DELT_LST

      END IF

      IST = 0

* ............................................ Loop over matrix diffusion zones

      DO NUM_ZON = 1,NUM_MTDZ

         IERR=0

         IN_IV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-4) 
         IN_RV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-3) 
         MXIV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-2)
         MXRV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-1)


         CALL UPD_PAR_MTDZ
     .(DELT, PARC, NZPAR,IV(IN_IV_MZ),MXIV_MZ,RV(IN_RV_MZ),
     . MXRV_MZ)


         IF (IERR.GT.0) THEN                !  Errors during parameters updating

            IST=1
            WRITE(MAINF,1001) IERR,NUM_ZON
 1001       FORMAT (' ERROR from UPD_PAR_DMT:',i5,' ERRORS during',
     .     ' parameters updating in matrix diffusion zone number:',i3)

         END IF

      END DO ! NUM_ZON

* ............................................... Stop, if too many errors found

      IF (IST.GT.0) THEN

         WRITE(MAINF,1030)
 1030    FORMAT(///' Too many ERRORS during',
     .     ' parameters updating'/' SORRY, BUT I STOP')
         STOP

      END IF

* ..............................................  Parameters updating successful

      I_ST_DMT = 3

      RETURN
      END

      SUBROUTINE UPD_DT_DMT
     .( DT_NEW)

************************************************************************
*
*   PURPOSE
*     Updates dimensionless time and related variables in 
*     response to a change in time increment
*
*
*   ARGUMENTS
*
*      DT_NEW              New time increment
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'MAIN_COMM.FOR'

      INCLUDE 'COMMON_DMT.FOR'

      INTEGER NUM_ZON,  
     .        MXIV_MZ, MXRV_MZ, IERR, IST,
     .        IN_IV_MZ, IN_RV_MZ

      REAL*8 DT_NEW

* ............................................................... Initialization

      IF (I_ST_DMT.LT.2) THEN

         PRINT 1000
 1000    FORMAT (' ERROR: UPD_DT_DMT reached ',
     .           ' before preliminary calculations')
         STOP

      END IF

      IST = 0

* ............................................. Loop over matrix diffusion zones

      DO NUM_ZON = 1,NUM_MTDZ

         IERR=0

         IN_IV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-4) 
         IN_RV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-3) 
         MXIV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-2)
         MXRV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-1)

         CALL UPD_DT_MTDZ 
     .(DT_NEW, IV(IN_IV_MZ),MXIV_MZ,RV(IN_RV_MZ),MXRV_MZ)

         IF (IERR.GT.0) THEN            !  Errors during time increment updating

            IST=1
            WRITE(MAINF,1001) IERR,NUM_ZON
 1001       FORMAT (' ERROR from UPD_DT_DMT:',i5,' ERRORS during',
     .   ' time increment updating in matrix diffusion zone number:',i3)

         END IF

      END DO ! NUM_ZON

* ............................................... Stop, if too many errors found

      IF (IST.GT.0) THEN

         WRITE(MAINF,1030)
 1030    FORMAT(///' Too many ERRORS during',
     .     ' time increment updating'/' SORRY, BUT I STOP')
         STOP

      END IF

* ........................................... Time incrment updating  successful

      I_ST_DMT = 4

      RETURN
      END

      SUBROUTINE COMP_CNTR_DMT
     .(BTRA, DTRA, NUMNP,  NUMEL, IDIMDTRA, KXX, LNNDEL, LMXNDL, IXDMT)

************************************************************************
* PURPOSE
*   Adds contribution to BTRA and DTRA
*
* ARGUMENTS
*
*      BTRA (NUMNP)        R.H.S. of transport equations
*      DTRA (NUMEL,IDIMDTRA)        Storage matrix of transport equations
*      NUMNP               Number of nodal points
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'MAIN_COMM.FOR'

      INCLUDE 'COMMON_DMT.FOR'

      INTEGER NUM_ZON,   LMXNDL, 
     .        MXIV_MZ, MXRV_MZ, IERR, IST,
     .        IN_IV_MZ, IN_RV_MZ, 
     .        NUMNP,IXDMT,  NUMEL,  IDIMDTRA, L,  K, NNUD, INUD, 
     ;        KXX(LMXNDL,NUMEL), LNNDEL(NUMEL)

      REAL*8 BTRA, DTRA(NUMEL,IDIMDTRA)

* ............................................................... Initialization

      IF (I_ST_DMT.LT.3) THEN

         PRINT 1000
 1000    FORMAT (' ERROR: COMP_CNTR_DMT reached ',
     .           ' before parameter updating')
         STOP

      END IF

C-------------------------------- Computes the number of elements in each node

      IF (I_CMP_NELBNUD.EQ.0) THEN
         DO L=1,NUMEL
            NNUD=LNNDEL(L)
            DO K=1,NNUD
               INUD=KXX(K,L)
               IV(I_NELBYNUD-1+INUD)=IV(I_NELBYNUD-1+INUD)+1
            ENDDO
         ENDDO
         I_CMP_NELBNUD=1
      ENDIF

C-------------------------------- Second, substracts DTRA by nodes to DTRA by elements
C-------------------------------- (uniformely distributed)

      DO L=1,NUMEL
         NNUD=LNNDEL(L)
         DO K=1,NNUD
            INUD=KXX(K,L)
            DTRA(L,K)=DTRA(L,K)-RV(I_DTRA-1+INUD)/IV(I_NELBYNUD-1+INUD)
         ENDDO
      ENDDO

      IST = 0

* ............................................. Loop over matrix diffusion zones

      DO NUM_ZON = 1,NUM_MTDZ

         IERR=0

         IN_IV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-4) 
         IN_RV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-3) 
         MXIV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-2)
         MXRV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-1)
         
         CALL COMP_CNTR_MTDZ
     .(NUM_ZON, BTRA, RV(I_DTRA),  NUMNP,
     . IV(IN_IV_MZ), MXIV_MZ, RV(IN_RV_MZ), MXRV_MZ, IXDMT)

         IF (IERR.GT.0) THEN   ! Errors computing contributions to BTRA and DTRA

            IST=1
            WRITE(MAINF,1001) IERR,NUM_ZON
 1001       FORMAT (' ERROR from UPD_DT_DMT:',i5,' ERRORS when',
     .   ' computing contributions to BTRA and DTRA',
     .   ' in matrix diffusion zone number:',i3)

         END IF

      END DO ! NUM_ZON

* ............................................... Stop, if too many errors found

      IF (IST.GT.0) THEN

         WRITE(MAINF,1030)
 1030    FORMAT(///' Too many ERRORS when',
     .     ' computing contributions to BTRA and DTRA'/
     .     ' SORRY, BUT I STOP')
         STOP

      END IF

C-------------------------------- Adds DTRA by nodes to DTRA by elements
C-------------------------------- (uniformely distributed)

      DO L=1,NUMEL
         NNUD=LNNDEL(L)
         DO K=1,NNUD
            INUD=KXX(K,L)
            DTRA(L,K)=DTRA(L,K)+RV(I_DTRA-1+INUD)/IV(I_NELBYNUD-1+INUD)
         ENDDO
      ENDDO

* .......................... Computing contributions to BTRA and DTRA successful

      I_ST_DMT = 5

      RETURN
      END




      SUBROUTINE COMP_DER_DMT
     .(   CCAL,CCALAN,   DERC, NUMNP, NPAR,IXDMT)

************************************************************************
* PURPOSE
*    Computes contributions of matrix diffusion to DERC (RHS of sensitivity 
*    equations)
*
* ARGUMENTS
*
*      CCAL (NUMNP)        C[k+1,i]
*      CCALAN (NUMNP)     C[k,i]
*      DERC (NUMNP,NPAR)   Derivatives of C[k+1,i] with respect to param
*      DERCAN (NUMNP,NPAR)Derivatives of C[k,i] with respect to param
*      NUMNP               Number of nodal points
*      NPAR                Number of parameters
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'MAIN_COMM.FOR'

      INCLUDE 'COMMON_DMT.FOR'

      INTEGER NUM_ZON,  
     .        MXIV_MZ, MXRV_MZ, IERR, IST,
     .        IN_IV_MZ, IN_RV_MZ,
     .        NUMNP, NPAR,IXDMT

      REAL*8 CCAL, CCALAN,  DERC

* ............................................................... Initialization

      IF (I_ST_DMT.LT.3) THEN

         PRINT 1000
 1000    FORMAT (' ERROR: COMP_DER_DMT reached ',
     .           ' before parameter updating')
         STOP

      END IF

      IST = 0

* ............................................. Loop over matrix diffusion zones

      DO NUM_ZON = 1,NUM_MTDZ

         IERR=0

         IN_IV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-4) 
         IN_RV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-3) 
         MXIV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-2)
         MXRV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-1)

         CALL COMP_DER_MTDZ
     .(NUM_ZON,    CCAL,CCALAN,   DERC, NUMNP, NPAR,
     . IV(IN_IV_MZ), MXIV_MZ, RV(IN_RV_MZ), MXRV_MZ, IXDMT)

         IF (IERR.GT.0) THEN       ! Errors when computing contributions to DERC

            IST=1
            WRITE(MAINF,1001) IERR,NUM_ZON
 1001       FORMAT (' ERROR from UPD_DT_DMT:',i5,' ERRORS when',
     .   ' computing contributions to DERC',
     .   ' in matrix diffusion zone number:',i3)

         END IF

      END DO ! NUM_ZON

* ............................................... Stop, if too many errors found

      IF (IST.GT.0) THEN

         WRITE(MAINF,1030)
 1030    FORMAT(///' Too many ERRORS when',
     .     ' computing contributions to DERC'/' SORRY, BUT I STOP')
         STOP

      END IF

* ................................... Computing contributions to DERC successful

      I_ST_DMT = 5

      RETURN
      END

      SUBROUTINE UPD_DREC_DMT
     .(   CCAL, CCALAN,   NUMNP,   DERC,DERCAN,   NPAR,   IERR,IXDMT)

************************************************************************
* PURPOSE
*   Updates derivatives of Io and In
*
* ARGUMENTS
*
*      CCAL (NUMNP)        C[k+1,i]
*      CCALAN (NUMNP)     C[k,i]
*      DERC (NUMNP,NPAR)   Derivatives of C[k+1,i] with respect to param
*      DERCAN (NUMNP,NPAR)Derivatives of C[k,i] with respect to param
*      NUMNP               Number of nodal points
*      NPAR                Number of parameters
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'MAIN_COMM.FOR'

      INCLUDE 'COMMON_DMT.FOR'

      INTEGER NUM_ZON,  
     .        MXIV_MZ, MXRV_MZ, IERR, IST,
     .        IN_IV_MZ, IN_RV_MZ,
     .        NUMNP, NPAR,IXDMT

      REAL*8 CCAL, CCALAN,  DERC,DERCAN

* ............................................................... Initialization

      IF (I_ST_DMT.LT.3) THEN

         PRINT 1000
 1000    FORMAT (' ERROR: UPD_DREC_DMT reached ',
     .           ' before parameter updating')
         STOP

      END IF

      IST = 0

* ............................................. Loop over matrix diffusion zones

      DO NUM_ZON = 1,NUM_MTDZ

         IERR=0

         IN_IV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-4) 
         IN_RV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-3) 
         MXIV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-2)
         MXRV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-1)

         CALL UPD_DREC_MTDZ
     .(NUM_ZON, CCAL,
     . CCALAN,  NUMNP,   DERC,DERCAN,   NPAR,
     . IV(IN_IV_MZ), MXIV_MZ, RV(IN_RV_MZ), MXRV_MZ, IXDMT)

         IF (IERR.GT.0) THEN     ! Errors when updating derivatives of Io and In

            IST=1
            WRITE(MAINF,1001) IERR,NUM_ZON
 1001       FORMAT (' ERROR from UPD_DREC_DMT:',i5,' ERRORS when',
     .   ' updating derivatives of Io and In',
     .   ' in matrix diffusion zone number:',i3)

         END IF

      END DO ! NUM_ZON

* ............................................... Stop, if too many errors found

      IF (IST.GT.0) THEN

         WRITE(MAINF,1030)
 1030    FORMAT(///' Too many ERRORS when',
     .     ' updating derivatives of Io and In'/' SORRY, BUT I STOP')
         STOP

      END IF

* ................................. updating derivatives of Io and In successful

      I_ST_DMT = 6

      RETURN
      END

      SUBROUTINE UPD_REC_DMT
     .(   CCAL,CCALAN,  NUMNP,IXDMT)

************************************************************************
* PURPOSE
*    Updates Io and In
*
* ARGUMENTS
*
*      CCAL (NUMNP)        C[k+1,i]
*      CCALAN (NUMNP)     C[k,i]
*      NUMNP               Number of nodal points
*
*   JCR, April,1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'MAIN_COMM.FOR'

      INCLUDE 'COMMON_DMT.FOR'

      INTEGER NUM_ZON,  
     .        MXIV_MZ, MXRV_MZ, IERR, IST,
     .        IN_IV_MZ, IN_RV_MZ,
     .        NUMNP, IXDMT

      REAL*8 CCAL, CCALAN

* ............................................................... Initialization

      IF (I_ST_DMT.LT.3) THEN

         PRINT 1000
 1000    FORMAT (' ERROR: UPD_DREC_DMT reached ',
     .           ' before parameter updating')
         STOP

      END IF

      IST = 0

* ............................................. Loop over matrix diffusion zones

      DO NUM_ZON = 1,NUM_MTDZ

         IERR=0

         IN_IV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-4) 
         IN_RV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-3) 
         MXIV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-2)
         MXRV_MZ = IV (IN_IV_DMT + 4*NUM_ZON-1)

         CALL UPD_REC_MTDZ
     .(NUM_ZON, CCAL,CCALAN,  NUMNP,
     . IV(IN_IV_MZ), MXIV_MZ, RV(IN_RV_MZ), MXRV_MZ, IXDMT)

         IF (IERR.GT.0) THEN                    ! Errors when updating Io and In

            IST=1
            WRITE(MAINF,1001) IERR,NUM_ZON
 1001       FORMAT (' ERROR from UPD_DT_DMT:',i5,' ERRORS when',
     .   ' updating Io and In',
     .   ' in matrix diffusion zone number:',i3)

         END IF

      END DO ! NUM_ZON

* ............................................... Stop, if too many errors found

      IF (IST.GT.0) THEN

         WRITE(MAINF,1030)
 1030    FORMAT(///' Too many ERRORS when',
     .     ' updating Io and In'/' SORRY, BUT I STOP')
         STOP

      END IF

* ................................................ updating Io and In successful

      I_ST_DMT = 7

      RETURN
      END
