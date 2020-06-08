      SUBROUTINE BALANCE_WRITE_FL
     &          (BM_ND    ,BM_ZN    ,I_PBL    ,I_REC    ,INDBALTYP
     &          ,INTI     ,IOBALH   ,IOBALNP  ,IOCONSRC ,IODENS
     &          ,IUBALH   ,NINT     ,NMAXNP   ,NMAXZP   ,NUMNP
     &          ,NZALF    ,NZARR    ,NZCHP    ,NZCON    ,NZQQP
     &          ,NZSTG    ,TABSOLUT ,TIME)

*******************************************************************************
*
* PURPOSE Writes flow mass balance information.
*
* DESCRIPTION This subroutine considers all posible 
*             cases, mixing these possibilities: 
*               1) Temporal mass balance at: 1.1) Simulation times
*                                            1.2) Observation times
*               2) Global mass balance
*               3) Detailed (nodal) mass balance
*             As a function of entry variables. it writes the appropr.  header
*             and the corresponding information.  
*
* EXTERNAL VARIABLES: ARRAYS
*
*  BM_ZN                  Array containing zonal flow or transport mass balance 
*                         information    
*  BM_ND                  Array containing nodal flow or transport mass balance 
*                         information
*  TIME                   Array containing observation times
*
* INTERNAL VARIABLES: ARRAYS
*
*  TOTAL                  Sum of mass balance contrib. over type of process
*  RELATIVE               Relative mass balance errors                                                      
*
* EXTERNAL VARIABLES: SCALARS
*
*  INDBALTYP              Current type of mass balance information: 
*                             0: steady state mass balance
*                             1: temporal mass balance
*                             2: global mass balance
*  INTI                   Observation time number such that the current         
*                         computation time lies in between observation time     
*                         number INTI and observation time number INTI+1        
*  IOBALH                 Option of computation of zonal flow mass balance 
*  IOBALNP                Option of computation of detailed (nodal) mass bal. It
*                         can be either IOBALDH or IOBALDC
*  I_PBL                  Problem number
*  I_REC                  Number of current simulation                    
*  IUBALH                  Unit number of the main output file (RES.OUT)         
*  NINT                   Number of observation times                           
*  NMAXZP                 Used to dimension BM_ZN (NMAXF or NMAXT)
*  NMAXNP                 Used to dimension BM_ND(NUMNP,NMAXNP). 
*                         Value=6 for flow and =4 for transport,
*  NUMNP                  Number of nodal points
*  NZALF                  Number of leakage zones     
*  NZARR                  Number of areal recharge zones          
*  NZCHP                  Number of prescribed head zones
*  NZQQP                  Number of prescribed flow zones          
*  NZSTG                  Number of storage capacity zones          
*  TABSOLUT               Current time (absolut) 
*
* INTERNAL VARIABLES: SCALARS
*
*  FORMATINI              Dummy index used to allocate position at IUBALH
*  FORMATEND              Dummy index used to allocate position at IUBALH
*  IMAX                   Number of rows of each mass balance info. pack 
*  INOD                   Node number. Dummy variable      
*  IPOINT                 Used to mark positions. See further details at
*                         BALANCE_FL
*  IPRO                   Process counter (dummy)
*  I                      Dummy
*
* INTERNAL VARIABLES: ARRAYS
*
*  RELATIVE               Relative mass balance error
*  TOTAL                  Sum over columns of mass balance information. 
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  BALANCE_WRITE_AUX      Writes flow or transport mass balance information for
*                         a paticular zonal/nodal process.     
*
* HISTORY: JLF: First coding (Dec-2000)
*          AAR: First coding (Dec-2000)
*
********************************************************************************

      IMPLICIT NONE

      INTEGER*4    MAXROW     ,INDBALTYP  ,I_PBL 
     &,INOD       ,INTI       ,IOBALH     ,IOBALNP    
     &,IPOS       ,IPRO       ,I          ,I_REC      ,IUBALH      
     &,NINT       ,NMAXZP     ,NMAXNP     ,NUMNP      ,NZALF      
     &,NZARR      ,NZCHP      ,NZCON
     &,NZQQP      ,NZSTG      ,IROW        ,IINI
     &,IEND  ,IFORMAT    ,IODENS,IOCONSRC

      REAL*8
     & BM_ZN    ,BM_ND    ,ERROR    ,RELATIVE ,TABSOLUT ,TOTAL    
     &,TIME,CONCSTG     

      DIMENSION TIME(NINT),BM_ZN(NMAXZP),BM_ND(NUMNP,NMAXNP)
     &         ,RELATIVE(7),TOTAL(7)

C--------------------  Internal

      INTEGER*4::IPSTG,IPOSC,IPOSH,IPOSQ,IPOSL,IPARR
      CHARACTER*104 LINE
      
C--------------------  Writes titles related to steady state flow mass

      WRITE(IUBALH,100)
      WRITE(IUBALH,1000) I_PBL

      IF (INDBALTYP.EQ.0) THEN ! Steady state mass balance

          WRITE(IUBALH,1010) TIME(INTI+1)  
      
      ELSE IF (INDBALTYP.EQ.1) THEN !Transient

          IF(I_REC.EQ.2) WRITE(IUBALH,1020)

          IF(IOBALH.LT.0) THEN
          
              WRITE(IUBALH,1020) TABSOLUT
              WRITE(IUBALH,1040)

          ELSE IF (IOBALH.GT.0) THEN
          
              WRITE(IUBALH,1020) TIME(INTI+1)
              WRITE(IUBALH,1040)

          END IF

      ELSE IF  (INDBALTYP.EQ.2) THEN !Global
      
          WRITE(IUBALH,1030)
          WRITE(IUBALH,1041)

      END IF !INDBALTYP.EQ.0,1,2,

  100 FORMAT('#')
 1000 FORMAT ('# FLOW PROBLEM #',I2,'.')

 1010 FORMAT('# STEADY-STATE MASS BALANCE. TIME:',G14.6)

 1020 FORMAT('# TRANSIENT MASS BALANCE. TIME',G14.6)
 1030 FORMAT('# GLOBAL MASS BALANCE.')

 1040 FORMAT('#',1X,38('-'))
 1041 FORMAT('#',1X,20('-'))
 1042 FORMAT('#',1X,94('-'))

 1050 FORMAT('# ZN    STORAGE      CONC.STG     RECHARGE     PRS.HEAD
     &     PRS.FLOW     LEAKAGE    SOLUTE MASS')
 
 1060 FORMAT('#',4X,7('[-----------]'))

C-------------------- Writes temporal flow mass balance information at IUBALH
 
      WRITE(IUBALH,100)
      WRITE(IUBALH,1050)
      WRITE(IUBALH,1060)


C-------------------- Initializes pointers to mass balance vector.

      IPSTG = 0
      IPARR = NZSTG*(IODENS + 1)
      IPOSH = IPARR + NZARR
      IPOSQ = IPOSH + NZCHP
      IPOSL = IPOSQ + NZQQP
      IPOSC = IPOSL + NZALF
  
C-------------------- Initialises arrays computation 

      TOTAL(:) = 0D0
      RELATIVE(:) = 0D0

      MAXROW = MAX(NZSTG,NZARR,NZCHP,NZQQP,NZALF,NZCON*IOCONSRC)

      DO IROW=1,MAXROW



          LINE=''

C-------------------- Storage capacity (only transient mass balance)

          IPOS = IPSTG + IROW

          IF(INDBALTYP.NE.0)THEN

              IF (IPOS.LE.NMAXZP) THEN

                  CALL BALANCE_WRITE_AUX 
     &                (BM_ZN(IPOS),IROW,NZSTG,LINE,1,12,TOTAL(1))


                  IF (IODENS.EQ.1) THEN

                      CONCSTG = BM_ZN(IPOS+NZSTG)

                  CALL BALANCE_WRITE_AUX
     &                    (CONCSTG,IROW,NZSTG,LINE,14,25,TOTAL(2))

                  ELSE

                      CALL BALANCE_WRITE_AUX
     &                    (0D0,2,1,LINE,14,25,TOTAL(2))

                  END IF !IODENS.EQ.1

              ELSE

                  CALL BALANCE_WRITE_AUX
     &                (0D0,2,1,LINE,1,12,TOTAL(1))

              END IF !IPOS.LE.NMAXZP

          ELSE

             CALL BALANCE_WRITE_AUX
     &            (0D0,2,1,LINE,1,12,TOTAL(1))

             CALL BALANCE_WRITE_AUX
     &            (0D0,2,1,LINE,14,25,TOTAL(2))

          END IF !INDBALTYP.NE.0


C-------------------- Areal Recharge

          IPOS = IPARR + IROW

          IF (IPOS.LE.NMAXZP) THEN

              CALL BALANCE_WRITE_AUX 
     &            (BM_ZN(IPOS),IROW,NZARR,LINE,27,38,TOTAL(3))

          ELSE

              CALL BALANCE_WRITE_AUX 
     &            (0D0,2,1,LINE,27,38,TOTAL(3))

          END IF !IPOS.LE.NMAXZP


C-------------------- Prescribed head

          IPOS = IPOSH + IROW

          IF (IPOS.LE.NMAXZP) THEN

              CALL BALANCE_WRITE_AUX 
     &            (BM_ZN(IPOS),IROW,NZCHP,LINE,40,51,TOTAL(4))

          ELSE

              CALL BALANCE_WRITE_AUX 
     &            (0D0,2,1,LINE,40,51,TOTAL(4))

          END IF !IPOS.LE.NMAXZP


C-------------------- Presecribed Flow

          IPOS = IPOSQ + IROW

          IF (IPOS.LE.NMAXZP) THEN

              CALL BALANCE_WRITE_AUX 
     &            (BM_ZN(IPOS),IROW,NZQQP,LINE,53,64,TOTAL(5))

          ELSE

              CALL BALANCE_WRITE_AUX 
     &            (0D0,1,1,LINE,53,64,TOTAL(5))


          END IF !IPOS.LE.NMAXZP


C-------------------- Leakage

          IPOS = IPOSL + IROW

          IF (IPOS.LE.NMAXZP) THEN

              CALL BALANCE_WRITE_AUX 
     &            (BM_ZN(IPOS),IROW,NZALF,LINE,66,77,TOTAL(6))

          ELSE

              CALL BALANCE_WRITE_AUX 
     &            (0D0,2,1,LINE,66,77,TOTAL(6))

          END IF !IPOS.LE.NMAXZP


C-------------------- Input mass

          IPOS = IPOSC + IROW

          IF (IODENS.EQ.1 .AND. IOCONSRC.EQ.1 .AND. IPOS.LE.NMAXZP) THEN

              CALL BALANCE_WRITE_AUX 
     &            (BM_ZN(IPOS),IROW,NZCON,LINE,79,90,TOTAL(7))

          ELSE

              CALL BALANCE_WRITE_AUX 
     &            (0D0,2,1,LINE,79,90,TOTAL(7))

          END IF !IPOS.LE.NMAXZP

          WRITE(IUBALH,1500) IROW,LINE


      END DO !IROW=1,MAXROW

 1500 FORMAT(I5,1X,A90)


C-------------------- Writes sum for each column

      WRITE(IUBALH,100)
      WRITE(IUBALH,1600) (TOTAL(I),I=1,7)
 1600 FORMAT('#SUM ',7(1X,G12.6))
 

C-------------------- Computes flow mass balance error

      ERROR = TOTAL(7) + TOTAL(6) + TOTAL(5) + TOTAL(4) + TOTAL(3)
     &      - TOTAL(2) - TOTAL(1)

C-------------------- Writes flow or transport mass balance error and
C-------------------- relative error


      WRITE(IUBALH,100)
      WRITE(IUBALH,1700) ERROR
      WRITE(IUBALH,1800)

      IF(NZSTG.NE.0.AND.INDBALTYP.NE.0.AND.TOTAL(1).NE.0D0) 
     &       RELATIVE(1)=DABS(ERROR*100.D0/TOTAL(1))
      IF(NZSTG.NE.0.AND.INDBALTYP.NE.0.AND.TOTAL(2).NE.0D0) 
     &       RELATIVE(2)=DABS(ERROR*100.D0/TOTAL(2))
      IF(NZARR.NE.0.AND.TOTAL(3).NE.0D0)
     &       RELATIVE(3)=DABS(ERROR*100.D0/TOTAL(3))
      IF(NZCHP.NE.0.AND.TOTAL(4).NE.0D0)
     &       RELATIVE(4)=DABS(ERROR*100.D0/TOTAL(4))
      IF(NZQQP.NE.0.AND.TOTAL(5).NE.0D0)
     &       RELATIVE(5)=DABS(ERROR*100.D0/TOTAL(5))
      IF(NZALF.NE.0.AND.TOTAL(6).NE.0D0)
     &       RELATIVE(6)=DABS(ERROR*100.D0/TOTAL(6))
      IF(NZCON.NE.0.AND.TOTAL(7).NE.0D0)
     &       RELATIVE(7)=DABS(ERROR*100.D0/TOTAL(7))

      WRITE(IUBALH,100)
      WRITE(IUBALH,1900) 
      WRITE(IUBALH,3000) 
      WRITE(IUBALH,1060)
      WRITE(IUBALH,3200) (RELATIVE(I),I=1,7)
      WRITE(IUBALH,100)
      WRITE(IUBALH,1042)

     
 1700 FORMAT('# FLOW MASS BALANCE ERROR:',G14.6)
 1800 FORMAT('#                         --------------')
 1900 FORMAT('# FLOW MASS BALANCE ERROR (%) RELATIVE TO:')
 3000 FORMAT('#       STORAGE      CONC.STG     RECHARGE     PRS.HEAD
     &     PRS.FLOW     LEAKAGE    SOLUTE MASS')

 3200 FORMAT(6X,7(G12.6,1X))

C-------------------- Writes flow and tpt. mass balance for nodes


      IF (IOBALNP.NE.0) THEN                                 ! Balance for nodes

          WRITE(IUBALH,100)
          WRITE(IUBALH,1000) I_PBL

          IF (IOBALNP.EQ.1) THEN

              IF (INDBALTYP.EQ.0) THEN

                  WRITE(IUBALH,2200) TIME(INTI+1)   ! Steady state mass balance
                  WRITE(IUBALH,2210)
      
              ELSE IF (INDBALTYP.EQ.1) THEN

                 WRITE(IUBALH,2300) TIME(INTI+1)     ! Transient mass balance
                 WRITE(IUBALH,2310)

              ELSE IF (INDBALTYP.EQ.2) THEN
              
                  WRITE(IUBALH,2400)  ! Global mass balance
                  WRITE(IUBALH,2410)
           
              END IF

              WRITE(IUBALH,100)
              WRITE(IUBALH,2800)
              WRITE(IUBALH,2900)

          END IF


 2200 FORMAT('# STEADY-STATE NODAL MASS BALANCE. TIME:',G14.6)
 2210 FORMAT('#',1X,38('-'))

 2300 FORMAT('# TRANSIENT NODAL MASS BALANCE AT TIME:',G14.6)
 2310 FORMAT('#',1X,47('-'))

 2400 FORMAT('# GLOBAL NODAL MASS BALANCE')
 2410 FORMAT('#',1X,25('-'))
 

 2800 FORMAT('#       STORAGE      CONC.STG     RECHARGE     PRS.HEAD
     &     PRS.FLOW     LEAKAGE    IN/OUT FLW   SOLUTE MASS')
 
 2900 FORMAT('#',4X,8('[-----------]'))        

          IFORMAT = 13

          IF (IOBALNP.EQ.1) THEN

              DO INOD=1,NUMNP

                  LINE=''

                  IINI = 1
                  IEND = 13

                  DO IPRO=1,NMAXNP  

                      WRITE(LINE(IINI:IEND),3400) BM_ND(INOD,IPRO)
            
            
                      IINI = IINI + IFORMAT
                      IEND = IEND + IFORMAT
            
                  END DO !IPRO=1,NMAXNP

                  WRITE(IUBALH,3500) INOD,LINE

              END DO !I=1,NUMNP

          END IF !IOBALNP.EQ.1

      WRITE(IUBALH,100)
      WRITE(IUBALH,1042)

      END IF

 3400 FORMAT(G12.6)
 3500 FORMAT(I5,1X,A104)

      END SUBROUTINE BALANCE_WRITE_FL
