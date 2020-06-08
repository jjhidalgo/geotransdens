       SUBROUTINE UPDATE_TIME_INCR
     &    (IOFLLI     ,IOTRLI     ,IODENS     ,IREDTIMH   ,IREDTIMC   
     &    ,IENTRY     ,MINCBI     ,NUCNVAR    ,INDENDDT   ,FCTINC
     &    ,DTINITIAL  ,DTMAXIO    ,TICAL      ,FCTDEC     ,TINTERVOBS
     &    ,DTAVGIO    ,TINC       ,TOLD       ,TINCINI    ,TINCLAST
     &    ,MXNRTV     ,NTRNRV     ,MAINF)

C PURPOSE This routine updates (increases, reduces or mantains) the time
C increment lengt during the time marching solution process. It increases the
C time step when smaller than the maximum permitted during the current
C observation intervale, and reduces it when convergence problems occur. For
C linear problems, it always set the step lengt to a constant value.
C ACRONYM: UPDATE TIME step INCREMENT
C DESCRIPTION Afects the time step by an incremental or decremental factor
C but constraining modifications by mean of controling parameters.
C REFERENCES: author GALARZA G.
C EXTERNALS
C tinc=current time step size. It will be used for trying to solve the
C      equations.
C told=previous time step increment which lead to convergence
C dtmaxio= user defined maximum time increment permitted at the ITOBS
C          observation intervale.
C fctinc=user defined time step incremental factor afecting the time step
C        when it is less than the maximum permitted
c fctdec=user defined time step decremental factor when convergence problems
*  MXNRTV                 Maximum number of time reductions defined by user                                                      
*  NTRNRV                 Number of time reductions                                                      
C mincbi=user defined minimum convergence number of times, obtained
C         after a time step reduction, required for increasing the time
C         increment again
c nucnvar=number of convergences obtined after the last time reduction
C ntrnrf=consequtive time reductions total number solving flow equation
C ntrnrt=consequtive time reductions total number solving transport equation
C iredtimh= indicator of convergence problems in flow equation
C iredtimc= indicator of convergence problems in transport equation
C iownrr=used defined time reduction information printout option
C ioflli=internal indicator of the linear or non-linear narure of the flow equation
C iotrli=internal indicator of the linear or non-linear narure of the trpt equation
C indenddt=indicate if the updated solution time corresponds to the final time of 
C          the current observation time.
C dtavgio=user defined (through variable KINT) uniform time step size in the 
C         current observation intervale, used as criterion by TRANSIN for 
C         defining the time step size.
C tintervobs=lenght of the current observation time intervale
C ticalan=diference between the previous solution time and that corresponding
C         to the initial time in the current observation intervale.
C dtinitial= first time step adopted in this observation intervale.
C ientry=indicator of the total times number that TRANSIN has acceded to
C        the routine update_time during the current observation intervale.
C INTERNALS
C HISTORY      G.G. 5/1997  First coding
C              JHG  12/2003 Control of number of time reductions added.

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C----------------------------------- solving non-linear problems

       IF ((IOFLLI+IOTRLI+IODENS).NE.0) THEN

         IF ((IREDTIMH+IREDTIMC).EQ.0) THEN     !no convergence problems

           TINCLAST=TINC
           TINCINI=TOLD
           TOLD=TINC
           IF (IENTRY.EQ.1) THEN                   !first step of the intervale
             TINC=DTINITIAL
           ELSE                                    !any step of the intervale
             NUCNVAR=NUCNVAR+1


             IF (TINC.LT.DTMAXIO.AND.NUCNVAR.GE.MINCBI)THEN !may increasing
               TINC=TINC*FCTINC
               NUCNVAR=0
               IF (TINC.GT.DTMAXIO) TINC=DTMAXIO            !limits the ste
             ENDIF


             IF (TICAL+TINC.GT.TINTERVOBS) THEN
               TINC=TINTERVOBS-TICAL
               INDENDDT=1
             ENDIF
            
           ENDIF

         ELSE                                          !convergence problems

           NTRNRV = NTRNRV +1
	     WRITE(MAINF,10)
   10     FORMAT(' TIME STEP REDUCTION BECAUSE OF CONVERGENCE PROBLEMS')

C_____________________________Code must stop. Number of time reductions 
C_____________________________greater than user allows

           IF (NTRNRV.GT.MXNRTV) THEN

              WRITE(MAINF,*)
     ;  ' SEVERAL CONVERGENCE PROBLEMS. TRANSIN MUST STOP.',
     ;  ' NUMBER OF TIME REDUCTIONS GREATER THAN ALLOWED.',
     ;  ' SEE THE USER GUIDE CAP.6 FOR HELPING SOLUTION OF THE PROBLEMS'

              WRITE(MAINF,*)' SORRY......'

              STOP

           ENDIF

           NUCNVAR=0
           TINC=TINC*FCTDEC
           INDENDDT=0           ! Unable to reach the forward obs. time

         ENDIF

C --------------------------------- linear problems
       ELSE

         TINCINI=TOLD
         TOLD=TINC
         TINC=DTAVGIO
       ENDIF

       RETURN
       END
