      SUBROUTINE WRITE_VAR_GS
     ;(IDIMVAR_GS  ,IOINV      ,IOPSC_GS  ,KTYPE_GS   ,MAINF   
     ;,MXMEASPP_GS ,MXNVAR_GS  ,NMEAS_GS  ,NPP_GS     ,NVAR_GS     
     ;,POSMEAS_GS  ,VMEAS_GS   ,VSTATS_GS)

********************************************************************************
*
* PURPOSE  Writes pilot points and sampling locations informations. variables 
*          value and external drift terms (only one group of zones), as well as
*          variables statistics
* 
* EXTERNAL VARIABLES: ARRAYS
*
*  POSMEAS_GS          Array containing coordinates of pilot points and sampling
*                      locations
*  VMEAS_GS            Array containing variables and external drifts values 
*  VSTATS_GS           Array containing variables statistics
*
* EXTERNAL VARIABLES: SCALARS
*  
*  IDIMVAR_GS         Used to dimension array VMEAS_GS
*  IOINV              Inverse problem option
*  IOPSC_GS           0: Conditional estimation; 1: Conditional simulation
*                     2: Unconditional simulation
*  KTYPE_GS           Kriging/cokriging option
*                         0: Simple kriging
*                         1: Residual kriging
*                         2: Kriging with locally varying mean
*                         3: Kriging with external drift (up to four)
*                         4: Simple cokriging
*                         5: Standardized ordinary cokriging
*                         6: Traditional ordinary cokriging
*  MAINF              Main output file number
*  MXMEASPP_GS        Used to dimension array VMEAS_GS
*  MXNVAR_GS          Maximum number of variables to be considered
*  NMEAS_GS           Number of sampling locations
*  NPP_GS             Number of pilot points
*  NVAR_GS            Number of variables (primary+all secondary)
*
* HISTORY: AAR   First coding    (Feb-2002)
*          AAR   Revision        (July-2003)
*
********************************************************************************

C______________________________________________________ Declaration of variables

      IMPLICIT NONE

      INTEGER*4 MAINF,NPP_GS,NMEAS_GS,NVAR_GS,IDIMVAR_GS,MXMEASPP_GS
     ;         ,MXNVAR_GS,KTYPE_GS,IOINV,IOPSC_GS
      REAL*8 POSMEAS_GS(MXMEASPP_GS,3),VMEAS_GS(MXMEASPP_GS,IDIMVAR_GS)
     ;      ,VSTATS_GS(MXNVAR_GS,4)
      INTEGER*4 INDEX,IVAR,J,IEX

C________________ Step 1: Writes pilot points and sampling locations information

      IF (IOINV.GT.0) THEN
        IF (KTYPE_GS.NE.2 .AND. KTYPE_GS.NE.3) THEN

                                                             ! Only Pilot points

          WRITE(MAINF,2200)
 2200     FORMAT(///,30X,' PILOT POINTS LOCATIONS',/
     ;            ,30X,' ===== ====== =========',//
     ;                ,'   ND         X         Y         Z'
     ;              ,/,'===================================',/)

          DO INDEX=1,NPP_GS
            WRITE(MAINF,3000) INDEX,(POSMEAS_GS(INDEX,J),J=1,3)
          END DO

        ELSE

                                              ! Pilot points and assigned drifts
          WRITE(MAINF,2250)
 2250     FORMAT(///,30X,' PILOT POINTS LOCATIONS AND ASSIGNED DRIFTS',/
     ;              ,30X,' ===== ====== ========= === ======== ======',//
     ;                  ,'   ND         X         Y         Z   DRIFT 1'
     ;                  ,'  DRIFT 2    DRIFT 3   DRIFT 4'
     ;                ,/,'============================================='
     ;                   '==============================',/)

          DO INDEX=1,NPP_GS
            WRITE(MAINF,3000) INDEX,(POSMEAS_GS(INDEX,J),J=1,3)
     ;            ,(VMEAS_GS(INDEX,IEX),IEX=MXNVAR_GS+1,IDIMVAR_GS)
          END DO
 
        END IF ! KTYPE_GS.NE.2 .OR. KTYPE_GS.NE.3

      END IF ! IOINV.GT.0
   
      IF (IOPSC_GS.LT.2) THEN
                                        ! Sampling locations and measured values
        WRITE(MAINF,2300)
 2300   FORMAT(///,30X,' SAMPLING LOCATIONS AND MEASURED VALUES',/
     ;            ,30X,' ======== ========= === ======== ======',//
     ;                ,'   ND         X         Y         Z   PRIMARY'
     ;                ,' EXTENS. 1 EXTENS. 2 EXTENS. 3 EXTENS. 4',/,
     ;                 '============================================='
     ;                 '========================================',/)

        DO INDEX=NPP_GS+1,NPP_GS+NMEAS_GS
          WRITE(MAINF,3000) 
     ;     INDEX-NPP_GS,(POSMEAS_GS(INDEX,J),J=1,3)
     ;    ,(VMEAS_GS(INDEX,IVAR),IVAR=1,MXNVAR_GS)
        END DO
                                        ! Sampling locations and assigned drifts
        IF (KTYPE_GS.EQ.2 .OR. KTYPE_GS.EQ.3) THEN

          WRITE(MAINF,2350)
 2350     FORMAT(///,30X,' SAMPLING LOCATIONS AND ASSIGNED DRIFTS',/
     ;              ,30X,' ======== ========= === ======== ======',//
     ;                  ,'   ND         X         Y         Z'
     ;                  ,'   DRIFT 1   DRIFT 2   DRIFT 3   DRIFT 4',/,
     ;                   '==================================='
     ;                   '========================================',/)

          DO INDEX=NPP_GS+1,NPP_GS+NMEAS_GS
            WRITE(MAINF,3000) 
     ;       INDEX-NPP_GS,(POSMEAS_GS(INDEX,J),J=1,3)
     ;      ,(VMEAS_GS(INDEX,IVAR),IVAR=MXNVAR_GS+1,IDIMVAR_GS)
          END DO

        END IF ! KTYPE_GS.NE.2 .OR. KTYPE_GS.NE.3

                                                           ! Variables statistics

        WRITE(MAINF,2003) 
        WRITE(MAINF,2004) INT(VSTATS_GS(1,1)),VSTATS_GS(1,2)
     ;                   ,VSTATS_GS(1,3)*VSTATS_GS(1,3)
        IF (NVAR_GS.GE.2) WRITE(MAINF,2005) INT(VSTATS_GS(2,1))
     ;            ,VSTATS_GS(2,2),VSTATS_GS(2,3)*VSTATS_GS(2,3)
        IF (NVAR_GS.GE.3) WRITE(MAINF,2006) INT(VSTATS_GS(3,1))
     ;            ,VSTATS_GS(3,2),VSTATS_GS(3,3)*VSTATS_GS(3,3)
        IF (NVAR_GS.GE.4) WRITE(MAINF,2007) INT(VSTATS_GS(4,1))
     ;            ,VSTATS_GS(4,2),VSTATS_GS(4,3)*VSTATS_GS(4,3)

 2003   FORMAT(///,25X,' SAMPLING LOCATIONS STATISTICS',/
     ;            ,25X,' ======== ========= ==========',//
     ;            ,'    VARIABLE NUM.   AVERAGE  VARIANCE ',/
     ;            ,'    ======== ====   =======  ======== ',/)

 2004   FORMAT('     PRIMARY',I5,2F10.3)
 2005   FORMAT(' EXTENSIVE 1',I5,2F10.3)
 2006   FORMAT(' EXTENSIVE 2',I5,2F10.3)
 2007   FORMAT(' EXTENSIVE 3',I5,2F10.3)
        
      END IF ! IOPSC_GS.LT.2

 3000 FORMAT(I5,8E10.3)

      RETURN
      END
