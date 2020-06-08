      REAL*8 FUNCTION GETVALUE2
     ;(KOL  ,MATDIM  ,NBAND  ,ROW  ,SIZE  ,MATRIX)
   
*******************************************************************************
*
* PURPOSE 
* to retrieve an element positioned at a given row (ROW) and column (KOL) of a
* given matrix (MATRIX) of a given size (MATDIM), with largest dimension SIZE 
* in symmetric band storage mode especially for the inverse covariance matrix
*
*******************************************************************************

      IMPLICIT NONE        

                                                   ! External variables: scalars
      INTEGER*4 SIZE,NBAND,ROW,KOL,POS,DUMMY,MATDIM,POSITION
                                                    ! External variables: arrays
      REAL*8 MATRIX(MATDIM)
                                                   ! Internal variables: scalars
      LOGICAL READY, SWAP

      READY =.FALSE.
      SWAP=.FALSE.

      IF (ABS(ROW-KOL) .GT. (NBAND-1)) THEN
         GETVALUE2=0.D0
         READY = .TRUE.
      ENDIF

      IF (ROW .EQ. KOL) THEN
         GETVALUE2 = MATRIX(ROW)
         READY = .TRUE. 
      ENDIF
   
      IF (KOL .GT. ROW .AND. READY .EQV. .FALSE.) THEN
         DUMMY = ROW
         ROW = KOL
         KOL=DUMMY
         SWAP= .TRUE.
      ENDIF

      IF (READY .EQV. .FALSE.) THEN
         POS = POSITION(KOL,SIZE,ROW)
         GETVALUE2 = MATRIX(POS)
      ENDIF

      IF (SWAP .EQV. .TRUE.) THEN
         DUMMY = ROW
         ROW = KOL
         KOL=DUMMY
         SWAP= .TRUE.
      ENDIF  

      RETURN
      END

*******************************************************************************
*******************************************************************************
*******************************************************************************

      REAL FUNCTION GETVALUE3
     ;(KOL  ,NDEVS  ,NPAR  ,NUMTOBS  ,OBSTYPE  ,ROW  ,IODEVICE
     ;,VJAC)

*******************************************************************************
*
* PURPOSE
*
*  This function is used to get values of the jacobian matrix of a certain 
*  measurement type. If the element asked for belongs to the specified obs. 
*  type  the returned value will be an entry in the jacobian matrix, else will 
*  be zero. 
*
*COMMENT 
*
* The function is used because defining a jacobian matrix for each 
* measurement type would ask more of the memory
*
*******************************************************************************

      IMPLICIT NONE
                                                   ! External variables: scalars
      INTEGER*4 NDEVS,NPAR, NUMTOBS,OBSTYPE,ROW,KOL
                                                    ! External variables: arrays
      INTEGER*4 IODEVICE(NDEVS+1,9)
      REAL*8 VJAC(NUMTOBS,NPAR)
                                                   ! Internal variables: scalars      
      INTEGER*4 NOF,NOL,I


      DO I=1,NDEVS                                           ! Loop over devices
        NOF=IODEVICE(I,8)                             ! First obs. within device
        NOL=IODEVICE(I+1,8)-1                          ! Last obs. within device
        IF (ROW .LE. NOL .AND. ROW .GE. NOF) THEN

           IF (OBSTYPE .EQ. IODEVICE(I,1)) THEN                    ! Same class?
               GETVALUE3=VJAC(ROW,KOL)
           ELSE
               GETVALUE3=0.D0
           ENDIF

        ENDIF

      ENDDO                 

      RETURN
      END              
