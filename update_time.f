       SUBROUTINE UPDATE_TIME
     ;           (IODENS     ,IOFLLI     ,IOTRLI     ,MINCBI
     ;           ,IREDTIMH   ,IREDTIMC   ,NUCNVAR    ,INDENDDT
     ;           ,IENTRY     ,INTI       ,NINT       ,NUMITER
     ;           ,NPARALG    ,IOPINVDT   ,IOINV      ,KINT
     ;           ,TINC       ,TOLD       ,DTMAXIO    ,FCTINC
     ;           ,FCTDEC     ,DTAVGIO    ,TINTERVOBS ,TICALAN
     ;           ,DTINITIAL  ,TICAL      ,DTIMEF     ,DTIMET
     ;           ,TABSOLUT   ,TIME       ,DTPREVINV  ,DTMXDS
     ;           ,PAR_DIR    ,TINCINI    ,TINCLAST   ,IOEQT
     ;           ,ISOLEQ     ,MXNRTV     ,NTRNRV     ,MAINF)

********************************************************************************
C PURPOSE: This routine performs all computations needed for updating 
C parameters associated with the time step increment. It takes in to
C account if a reduction or increment of the step has to be done.
C For selecting the time step size, it also applies the criteria
C proposed by GALARZA et al., 1997, for taking advantage from the inverse
C problem solution frame.
C
C ACRONYM: UPDATEs all TIME associated parameters
C
C DESCRIPTION
C consists of a serie of conditioned calls to particular routines which
C compute several parameters associated to the variable time in the
C simulation
C
C REFERENCES: author GALARZA G.
C GALARZA G. CARRERA J. AND MEDINA A. Computational aspects in non-linear
C optimization problems. En revision por Int. Jour. of Computational
C methods in engineering
C
C EXTERNALS
C tinc= current time step size. It will be used for trying to solve the
C      equations.
C told= previous time step increment which lead to convergence
C dtmaxio= user defined maximum time increment permitted at the INTI
C          observation intervale.
C fctinc=user defined time step incremental factor afecting the time step
C        when it is less than the maximum permitted
c fctdec=user defined time step decremental factor when convergence problems
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
C tical=diference between the updated solution time and the initial time
C         of the current observation intervale.
C dtimef=realtive (to the observation intervale length) value of (t(k+theta)-tobs(inti)),
C         where the first term is the time where the flow equation will be  solved, 
C         and the second is the initial time of the current observation intervale.
C dtimet= idem for transport
C tabsolut=absolute current solution time
C time=vector containing the observation times
C inti= current observation intervale index.
*  MXNRTV                 Maximum number of time reductions defined by user                                                      
*  NTRNRV                 Number of time reductions                                                      
C nint=total observation times number
C numiter= current inverse problem iteration number
C iopinvdt= user defined option for adopting or not the first time increment
C           in the current INTI observation intervale, equal to
C           that used during the previous inverse iteration at the same
C           INTI intervale, which lead to convergence.
C dtprevinv(nint)= array storing the first time increment used during
C                  the previous inverse iteration, in the INTI intervale,
C                  which lead to convergence.
C ioinv=indicator of the type of inverse problem.
C INTI=current observation intervale number.
C dtmxds(nint)=user defined maximum time step lenght permitted during the current
C             observation intervale.
C kint(nint)=user defined vector containing the total steps number he
C            suggests during the time marching process at the current
C            observation intervale.
C
C INTERNALS
C
C HISTORY  G.G.  5/1997 First coding   
C          JHG  12/2003 Matrix diffusion terms.
C          AMS  2007    Elimination of matrix diffusion terms
********************************************************************************


       IMPLICIT DOUBLE PRECISION (A-H,O-Z)

       DIMENSION PAR_DIR(NPARALG) ,KINT(NINT) ,DTMXDS(NINT)

C----------------------------- Update the indicator of the number of times
C                              that the program has acceeded to this routine
C                              during the current INTI intervale
       IENTRY=IENTRY+1

C----------------------------- Compute time constants associated to the 
C                              current INTI intervale

       IF (IENTRY.EQ.1) CALL COMP_TIME_CONSTANTS
     ;    (KINT       ,IOINV      ,INTI      ,NUMITER    ,IOPINVDT
     ;    ,NINT       ,TIME       ,DTINITIAL ,DTMXDS     ,DTMAXIO
     ;    ,DTPREVINV  ,DTAVGIO    ,TINTERVOBS,IOFLLI     ,IOTRLI)

C----------------------------- Update the time step

       CALL UPDATE_TIME_INCR
     ;    (IOFLLI     ,IOTRLI     ,IODENS     ,IREDTIMH   ,IREDTIMC   
     ;    ,IENTRY     ,MINCBI     ,NUCNVAR    ,INDENDDT   ,FCTINC
     ;    ,DTINITIAL  ,DTMAXIO    ,TICAL      ,FCTDEC     ,TINTERVOBS
     ;    ,DTAVGIO    ,TINC       ,TOLD       ,TINCINI    ,TINCLAST
     ;    ,MXNRTV     ,NTRNRV     ,MAINF)

C----------------------------- Update the solution time

       CALL UPDATE_SOL_TIME
     ;(DTIMEF      ,DTIMET   ,DTINITIAL
     ;,IENTRY   ,INDENDDT  ,INTI     ,IOEQT    ,IREDTIMC
     ;,IREDTIMH    ,NINT     ,TABSOLUT  ,PAR_DIR(29),PAR_DIR(30),TICAL
     ;,TICALAN     ,TINC     ,TINTERVOBS,TIME     ,ISOLEQ)    

       RETURN
       END
