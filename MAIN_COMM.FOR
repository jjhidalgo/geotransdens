********************************************************************************
*
* DESCRIPTION Common block containing TRANSIN-IV big arrays, RV for real numbers
*             IV for integers and KV for strings (characters chains).
*
* COMMON VARIABLES: ARRAYS
*
*  IV                     Integer array used to reserve most of the necessary
*                         space for all the integer variables of the problem.
*                         Only some small arrays are not included in it:
*                         IFLAGS, INORPAR, IOLG_PAR, IOWRITE ,IPAR_DIR,
*                         IPAR_INV and NZONE_PAR
*  KV                     Character array used to reserve most of the necessary
*                         space for all the string variables of the problem.
*  RV                     Real array used to reserve most of the necessary
*                         space for all the real variables of the problem.
*                         Only some small arrays are not included in it:
*                         PAR_DIR, PAR_INV, PAR_WGT and FOBJ_WGT
*
* COMMON VARIABLES: SCALARS
*
*  IIMAX                  Used to dimension IV array
*  IKMAX                  Used to dimension KV array
*  IRMAX                  Used to dimension RV array
*  LASTII                 Last memory position allocated in IV
*  LASTIR                 Last memory position allocated in RV
*
* HISTORY: First coding: JCR (April-1999)
*          Revision and header: AAR (March-00)
*
********************************************************************************

      REAL*8 RV
      INTEGER*4 IV, LASTII, LASTIR, IRMAX, IIMAX, IKMAX
      CHARACTER*10 KV
      PARAMETER (IRMAX=60000000, IIMAX=50000000 ,IKMAX=7000)
      COMMON /C_MAIN/ RV(IRMAX), IV(IIMAX), KV(IKMAX), LASTII, LASTIR
