      SUBROUTINE JVJ
     ;(IDIMCOV  ,IDIMHESS ,NBANDCOV ,NDEVS    ,NPAR     ,NSTAT
     ;,NUMTOBS  ,OBSTYPE  ,COVINV   ,IODEVICE ,OBSCLASS ,OUT
     ;,VJAC)

********************************************************************************
*
* PURPOSE:
*                        t   -1
* calculate the product J * V  *  J    where J(i) is the part of the Jacobian
*                        i   i     i
* matrix related to measurement type i. It has the same dimension as the
* jacobian but all non desired elements are zero. V(i) is the part of the
* covariance matrix of measurement type i and desserves the same commments as J
*
* DESCROPTION
*
* - Step 0: Declaration and initialization of variables
* - Step 1: Product Jt*V-1*J
*
* EXTERNAL VARIABLES: SCALARS
*
*  NUMTOBS                Total number of observations
*  NDEVS                  Number of devices  
*  NPAR                   Total number of parameters to be estimated   
*  IDIMHESS               Used to dimension array HESS. It is equal to          
*                         NPAR*(NPAR+1)/2  
*  IDIMCOV                Used to dimension array COVINV
*  OBSTYPE                scalar between 1 and 10 inicating the measurement type
*  NBANDCOV               Band width of the inverse covariance matrix   
*  NSTAT                  Maximum number of state variables whose data is used  
*                         for calibration (used for dimensioning)    
*
* EXTERNAL VARIABLES: ARRAYS
*
*  COVINV                 Inverse of the covariance matrix 
*  IODEVICE               Column 1: Data type                                   
*                         Column 2: Status for calc. of obs.                    
*                         Column 3: Method of spat. integr.                     
*                         Column 4: Method of temp. integr.                     
*                         Column 5: Number of integr. time                      
*                         Column 6: Covariance matrix type
*                         Column 9: Related problem
*  OBSCLASS               Array containing number of devices where a particular
*                         type of state var. was measured
*                         COLUMN 1 contains number of devices, while COLUMNS 2-?
*                         contain identifiers to those devices
*                         ROW is related to state var. type
*  VJAC                   Jacobian matrix          
*
* INTERNAL VARIABLES: SCALARS
*
*  BAND                   Band nr of the covariance matrix for which calculations
*                         are being carried out
*  ELEM                   the value of the entry in covariance matrix 
*  KOL                    kolumn index of element in matrix
*  I,J                    dummy countersof loop over parameters
*  K1,L1                  dummy counter of loop over devices
*  K2, L2                 dummy counter of loop over measurements
*  KDEV, LDEV             device number
*  JVJIJ                  element(i,j) of the result matrix of the calculation
*  LCOVINVLK              element(l,k) of the covariance matrix
*  LOCALJACKJ             element(k,j) of the jacobian matrix of measurement type
*  LOCALJACLI             element(l,i) of the jacobian matrix of measurement type
*  NOF                    First observation defining a particular device.
*  NOL                    Last observation defining a particular device.
*  POS                    the array index of an element
*  ROW                    row index of element in matrix
*
* EXTERNAL FUNCTIONS REFERENCED
*
*  GETVALUE2              retrieve an element of a given size in symmetric
*                         band storage mode
*  GETVALUE3              get values of the jacobian matrix of a certain 
*                         measurement type.
*
* HISTORY:   LJS (First coding), Dec 2002
*            AAR (Revision and formatting). Jan 2003
*
********************************************************************************

C___________________________ Step 0: Declaration and initialization of variables

      IMPLICIT NONE

                                                   ! EXTERNAL VARIABLES: SCALARS
      INTEGER*4 NUMTOBS, NDEVS, NPAR, IDIMHESS, IDIMCOV,OBSTYPE,
     ;          NBANDCOV,NSTAT
                                                    ! EXTERNAL VARIABLES: ARRAYS
      INTEGER*4 OBSCLASS(NSTAT,NDEVS+1),IODEVICE(NDEVS+1,9)
      REAL*8 VJAC(NUMTOBS,NPAR), COVINV(IDIMCOV),OUT(IDIMHESS)
                                                   ! INTERNAL VARIABLES: SCALARS
      INTEGER*4 I,J,ROW,POS
      INTEGER*4 K1, K2,L1,L2,KNOF,LNOF,KNOL,LNOL,LDEV,KDEV
      REAL*8 GETVALUE2,GETVALUE3,LCOVINVLK,JVJIJ,LOCALJACLJ
     :       ,LOCALJACKI,LCOVINVKK,LOCALJACKJ
  
      CALL ZERO_ARRAY(OUT,IDIMHESS)
      ROW=0

C____________________________________________________ Step 3: multiply Jt*V-1*J

      DO K1=2,OBSCLASS(OBSTYPE,1)+1        !loop over devices 
          KDEV = OBSCLASS(OBSTYPE,K1) 

          IF (IODEVICE(KDEV,6).EQ.1) THEN !diagonal covariance matrix
C=============================DIAGONAL COVARIANCE MATRIX

               KNOF=IODEVICE(KDEV,8)
               KNOL=IODEVICE(KDEV+1,8)-1
               DO K2=KNOF, KNOL                 !loop over measurements 

                  LCOVINVKK=GETVALUE2
     ;            (K2  ,IDIMCOV  ,NBANDCOV  ,K2  ,NUMTOBS  ,COVINV)

                  DO I=1,NPAR       ! Loop over estim. parameters
                      LOCALJACKI=GETVALUE3(I ,NDEVS  ,NPAR  
     ;                ,NUMTOBS  ,OBSTYPE  ,K2  ,IODEVICE,VJAC)  
                       
					 IF (LOCALJACKI .NE. 0d0  )THEN 
                    ! Only to I because the matrix is symmetric
                           DO J=1,I  
                                  
                               LOCALJACKJ=GETVALUE3(J  ,NDEVS  ,NPAR  
     ;                        ,NUMTOBS  ,OBSTYPE ,K2  ,IODEVICE,VJAC)   
                              POS=I*(I+1)/2  - (I-J)    !"out" in symmetric storage mode
                              JVJIJ=LOCALJACKI*LCOVINVKK*LOCALJACKJ 
                              OUT(POS)=OUT(POS)+JVJIJ
                           ENDDO !J=1,I  
                       ENDIF !IF (LOCALJACKI .NE. 0d0  )THEN
                  ENDDO !DO I=1,NPAR
	         ENDDO  !K2

C=============================NONDIAGONAL COVARIANCE MATRIX
	    ELSE
               KNOF=IODEVICE(KDEV,8)
               KNOL=IODEVICE(KDEV+1,8)-1
               DO K2=KNOF, KNOL                 !loop over measurements 
                 
                  DO L1=2,OBSCLASS(OBSTYPE,1)+1 !loop over devices
                     LDEV = OBSCLASS(OBSTYPE,L1) 
                     LNOF=IODEVICE(LDEV,8)
                     LNOL=IODEVICE(LDEV+1,8)-1
                     DO L2=LNOF, LNOL         !loop over measurements
                        
                       LCOVINVLK=GETVALUE2
     ;(K2  ,IDIMCOV  ,NBANDCOV  ,L2  ,NUMTOBS  ,COVINV)
       
        !proceed only if there is covariance between measurements l2 and k2
                       IF(LCOVINVLK .NE. 0D0) THEN
                           
                          DO I=1,NPAR       ! Loop over estim. parameters
                             LOCALJACKI=GETVALUE3(I ,NDEVS ,NPAR  
     ;,NUMTOBS  ,OBSTYPE  ,K2  ,IODEVICE,VJAC)  
                            IF (LOCALJACKI .NE. 0d0  )THEN 
                    ! Only to I because the matrix is symmetric
                               DO J=1,I  
                                  
                                  LOCALJACLJ=GETVALUE3(J  ,NDEVS  
     ; ,NPAR  ,NUMTOBS  ,OBSTYPE ,L2  ,IODEVICE,VJAC)   
                                  POS=I*(I+1)/2  - (I-J)    !"out" in symmetric storage mode
                                  JVJIJ=LOCALJACKI*LCOVINVLK*LOCALJACLJ 
                                  OUT(POS)=OUT(POS)+JVJIJ
                               ENDDO !J=1,I  
                            ENDIF !IF (LOCALJACKI .NE. 0d0  )THEN
                          ENDDO !DO I=1,NPAR         
                       ENDIF ! IF (LCOVINVLK .NE. 0d0  )THEN 
                    ENDDO !DO  L2=LNOF, LNOL
                ENDDO !DO L1=2,OBSCLASS(OBSTYPE,1)+1
             ENDDO !DO K2=KNOF, KNOL
	    ENDIF  !DIAGONAL COVARIANCE MATRIX
      ENDDO  ! DO K1=2,OBSCLASS(OBSTYPE,1)+1
   !  now we are going to get the elements k2,i  and l2,j of the jacobian                      
    
      RETURN 
      END
