      SUBROUTINE TEMP_WEIGHT_OBS
     ;(IOINTTYP ,NO       ,NUMINT   ,NUMTOBS  ,NUMTOBSC ,IOTINT
     ;,TOBS     ,WTOBST)  

***********************************************************************
* PURPOSE
*
* Determines the temporal weights associated with all integration
* times.
*
* DESCRIPTION
*
* The temporal weights are used to calculate the model-output
* corresponding to observations. For each integration time, a weighted
* sum of nodal values is calculated. Next, this value is multiplied
* with the temporal weight. This is done in COMP_OBS. 
*
* Throughout the subroutine, IT is used in comments as abbreviation
* for integration time.
*
* EXTERNAL VARIABLES: ARRAYS
*
*  IOTINT                 Temporal integration (=1) or temporal averaging (=2)  
*  TOBS                   Time of observation                                   
*  WTOBST                 Weight for integration time                           
*
* EXTERNAL VARIABLES: SCALARS
*
*  IOINTTYP               Device time integration type
*  NO                     Last observation of current device
*  NUMINT                 Number of integration times for a given observation
*  NUMTOBS                Total number of observations                          
*  NUMTOBSC               Last observation on previous device
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
* HISTORY
*
*     CK      11-1999     First coding
*     AAR     03-2001     Revision
*
***********************************************************************

      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION WTOBST(NUMTOBS),IOTINT(NUMTOBS),TOBS(NUMTOBS,2)

C______________________________ Loop over observation of current device

      DO I=NUMTOBSC+1,NO

C_______________________One point-in-time

        IF (IOINTTYP.EQ.0) THEN

          WTOBST(I)=1D0

C_______________________ITs defined by simulation

        ELSE IF (IOINTTYP.EQ.1) THEN

          IF (IOTINT(I).EQ.0) THEN                    ! Particular obs.= average
            WTOBST(I)=1D0/(TOBS(I,2)-TOBS(I,1)) 
          ELSE IF (IOTINT(I).EQ.1) THEN                ! Particular obs.= integ.
            WTOBST(I)=1D0                       
          ENDIF

C_______________________Equally distributed simulation times

        ELSE IF (IOINTTYP.EQ.2) THEN

          IF (IOTINT(I).EQ.0) THEN                    ! Particular obs.= average
            WTOBST(I)=1D0/(1D0*NUMINT-1D0)                   
          ELSE IF (IOTINT(I).EQ.1) THEN                ! Particular obs.= integ.
            WTOBST(I)=(TOBS(I,2)-TOBS(I,1))/(NUMINT*1D0-1D0) 
          ENDIF

        ENDIF
      ENDDO

      RETURN

      END
