      SUBROUTINE BALANCE_WRITE_TR
     &          (BM_ND    ,BM_ZN    ,I_PBL    ,I_REC    ,IBTCO
     &          ,INDBALTYP,INTI     ,IOBALC   ,IOBALNP  ,IUBALC
     &          ,NINT     ,NMAXT    ,NUMNP    ,NZCLK    ,NZCOE
     &          ,NZFOD    ,NZMDIF   ,NZPOR    ,NZZOR    ,TABSOLUT
     &          ,TIME)

*******************************************************************************
*
* PURPOSE Writes transport mass balance information.
*
* DESCRIPTION This subroutine is absolutelly general. It considers all posible 
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
*  IOBALC                 Option of computation of zonal transport mass balance 
*  IOBALNP                Option of computation of detailed (nodal) mass bal. It
*                         can be either IOBALDH or IOBALDC
*  I_PBL                  Problem number
*  I_REC                  Number of current simulation                    
*  IUBALC                  Unit number of the main output file (RES.OUT)         
*  NINT                   Number of observation times                           
*  NMAXZP                 Used to dimension BM_ZN (NMAXF or NMAXT)
*  NMAXNP                 Used to dimension BM_ND(NUMNP,NMAXNP). 
*                         Value=6 for flow and =4 for transport,
*  NUMNP                  Number of nodal points
*  NZALF                  Number of leakage zones     
*  NZARR                  Number of areal recharge zones          
*  NZCHP                  Number of prescribed head zones
*  NZCOE                  Number of external concentration zones         
*  NZFOD                  Number of first order decay reaction zones         
*  NZPOR                  Number of porosity zones         
*  NZQQP                  Number of prescribed flow zones          
*  NZSTG                  Number of storage capacity zones          
*  NZZOR                  Number of zero order reactions zones          
*  TABSOLUT               Current time (absolut) 
*
* INTERNAL VARIABLES: SCALARS
*
*  FORMATINI              Dummy index used to allocate position at IUBALC
*  FORMATEND              Dummy index used to allocate position at IUBALC
*  IMAX                   Number of rows of each mass balance info. pack 
*  INOD                   Node number. Dummy variable      
*  IPOINT                 Used to mark positions. See further details at
*                         BALANCE_FL or BALANCE_TR      
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

C-------------------- External

      INTEGER*4::I_PBL    ,I_REC    ,INDBALTYP,INTI     ,IOBALC
     &          ,IOBALNP  ,IUBALC   ,NINT     ,NMAXT    ,NUMNP
     &          ,NZCLK    ,NZCOE    ,NZFOD    ,NZMDIF   ,NZPOR
     &          ,NZZOR

      REAL*8::TABSOLUT

      INTEGER*4::IBTCO(NUMNP)

      REAL*8::TIME(NINT),BM_ZN(NMAXT),BM_ND(NUMNP,12)

C-------------------- Internal
     
      INTEGER*4::I        ,IEND     ,IFORMAT  ,IINI     ,IPCLK
     &          ,IPFLX    ,IPFOD    ,IPMATDIF ,IPPOR    ,IPOS
     &          ,IPOS1    ,IPRO     ,IPZOD    ,IROW     ,IT
     &          ,MAXROW

      REAL*8::ERROR    

      REAL*8::RELATIVE(12),TOTAL(12)

      CHARACTER*156 LINE

            
C-------------------- First executable statement

C-------------------- Writes titles related to steady state tpt. mass 
C-------------------- balance at IUBALC

      WRITE(IUBALC,100)
      WRITE(IUBALC,1000) I_PBL

      IF(INDBALTYP.EQ.0) THEN ! Steady state mass balance
      
          WRITE(IUBALC,1010) TIME(INTI+1) 
          WRITE(IUBALC,1040)
 
      ELSE IF (INDBALTYP.EQ.1) THEN !Transient

          IF (I_REC.EQ.2) WRITE(IUBALC,1020)

          IF (IOBALC.LT.0) THEN
          
              WRITE(IUBALC,1020) TABSOLUT
              WRITE(IUBALC,1040)

          ELSE IF (IOBALC.GT.0) THEN

              WRITE(IUBALC,1020) TIME(INTI+1)
              WRITE(IUBALC,1040)

          END IF !IOBALC.LT.0
      
      ELSE IF (INDBALTYP.EQ.2) THEN ! Global mass balance

         WRITE(IUBALC,1030)
         WRITE(IUBALC,1041)

      END IF !INDBALTYP.EQ.0,1,2

      WRITE(IUBALC,100)
      WRITE(IUBALC,1050)
      WRITE(IUBALC,1060)
      WRITE(IUBALC,1070)
      WRITE(IUBALC,1080)

  100 FORMAT ('#')

 1000 FORMAT ('# TRANSPORT PROBLEM #',I2,'.')

 1010 FORMAT('# STEADY-STATE MASS BALANCE. TIME:',G14.6)

 1020 FORMAT('# TRANSIENT MASS BALANCE. TIME',G14.6)
 1030 FORMAT('# GLOBAL MASS BALANCE.')
 
 1040 FORMAT('# --------------------------------------')
 1041 FORMAT('# ----------------------')
 1042 FORMAT('#',1X,147('-'))

 1050 FORMAT('#'24X,'[MASS FLOW BNDRS. SORTED BY CON. ZONE]')
 1060 FORMAT('#',17X,'[',50('-'),']')
 1070 FORMAT('# ZN    POROSITY     PRS.HEAD     PRS.FLOW     LEAKAGE
     &      RECHARGE     PRS.CONC     MASS INP     CON LEAK      F.O.D.
     &       Z.O.D.      MAT.DIFF')

 1080 FORMAT('#',4X,11('[-----------]'))

C-------------------- Writes temporal/global tpt. mass balance 
C-------------------- information at IUBALC

      MAXROW = MAX(NZPOR,NZFOD,NZZOR,NZCOE,NZCLK,NZMDIF)

C-------------------- Initilizes pointers

      IPPOR = 0
      IPFLX = NZPOR
      IPCLK = IPFLX + NZCOE*6
      IPFOD = IPCLK + NZCLK
      IPZOD = IPFOD + NZFOD
      IPMATDIF = IPZOD + NZZOR

C-------------------- Initialises arrays computation 

      TOTAL(:) = 0D0
      RELATIVE(:) = 0D0


      DO IROW=1,MAXROW
    
          LINE = ''

          IPOS = IPPOR + IROW

C-------------------- Storage capacity (only transient mass balance).

          IF(INDBALTYP.NE.0 .AND .NZPOR.NE.0) THEN

              CALL BALANCE_WRITE_AUX 
     &            (BM_ZN(IPOS),IROW,NZPOR,LINE,1,12,TOTAL(1))

          ELSE

              CALL BALANCE_WRITE_AUX 
     &            (0D0,2,1,LINE,1,12,TOTAL(1))

          END IF !INDBALTYP.NE.0

C-------------------- Ext. Conc. (and mass fluxes).

          IPOS = IPFLX + IROW

          DO I=1,6

              IPOS1 = IPFLX + 6*(IROW-1) + I
              IINI = 14 + 13*(I-1) !13 para dejar un hueco
              IEND = 25 + 13*(I-1)
              IT = I + 1

              IF(IPOS1.LE.NMAXT) THEN

                  CALL BALANCE_WRITE_AUX 
     &                (BM_ZN(IPOS1),IROW,NZCOE,LINE,IINI,IEND,TOTAL(IT))

              ELSE

                  CALL BALANCE_WRITE_AUX 
     &                (0D0,2,1,LINE,IINI,IEND,TOTAL(IT))

              END IF !IPOS.LE.NMAXZP

          END DO !I=1,6


C-------------------- Conc. Leakage.

          IPOS = IPCLK + IROW

          IF(IPOS.LE.NMAXT) THEN 

              CALL BALANCE_WRITE_AUX 
     &            (BM_ZN(IPOS),IROW,NZCLK,LINE,92,103,TOTAL(8))


          ELSE

              CALL BALANCE_WRITE_AUX 
     &            (0D0,2,1,LINE,92,103,TOTAL(8))

          END IF !IPOS.LE.NMAXZP

C-------------------- F.O.D. reactions

          IPOS = IPFOD + IROW

          IF(IPOS.LE.NMAXT) THEN 

              CALL BALANCE_WRITE_AUX 
     &            (BM_ZN(IPOS),IROW,NZFOD,LINE,105,116,TOTAL(9))

          ELSE

              CALL BALANCE_WRITE_AUX 
     &            (0D0,2,1,LINE,105,116,TOTAL(9))

          END IF !IPOS.LE.NMAXZP


C-------------------- Z.O.R. reactions

          IPOS = IPZOD + IROW

          IF(IPOS.LE.NMAXT) THEN

              CALL BALANCE_WRITE_AUX 
     &            (BM_ZN(IPOS),IROW,NZZOR,LINE,118,129,TOTAL(10))

          ELSE

              CALL BALANCE_WRITE_AUX 
     &            (0D0,2,1,LINE,118,129,TOTAL(10))

          END IF !IPOS.LE.NMAXZP

C-------------------- Matrix Diffusion.

          IPOS = IPMATDIF + IROW

          IF(IPOS.LE.NMAXT) THEN

              CALL BALANCE_WRITE_AUX 
     &            (BM_ZN(IPOS),IROW,NZMDIF,LINE,131,142,TOTAL(11))

          ELSE

              CALL BALANCE_WRITE_AUX 
     &            (0D0,2,1,LINE,131,142,TOTAL(11))

          END IF !IPOS.LE.NMAXZP

         WRITE(IUBALC,1500) IROW,LINE
1500     FORMAT(I5,1X,A146)

      END DO !ROW=1,MAXROW

C-------------------- Writes sum for each column

      WRITE(IUBALC,100)
      WRITE(IUBALC,1610) (TOTAL(I),I=1,11)

 1610 FORMAT('#SUM ',11(1X,G12.6))

C-------------------- Computes tpt. mass balance error

      ERROR = SUM(TOTAL(2:11)) - TOTAL(1)

C-------------------- Writes flow or transport mass balance error and
C-------------------- relative error


      WRITE(IUBALC,100)
      WRITE(IUBALC,1710) ERROR
      WRITE(IUBALC,1800)

      IF(NZPOR.NE.0.AND.INDBALTYP.NE.0.AND.TOTAL(1).NE.0D0)
     &         RELATIVE(1)=DABS(ERROR*100.D0/TOTAL(1))

      IF(NZCOE.NE.0) THEN
        
          DO I=2,7
          
              IF(TOTAL(I).NE.0D0)
     &            RELATIVE(I)=DABS(ERROR*100.D0/TOTAL(I))

          END DO !I=1,7

      END IF !NZCOE.NE.0

      IF(NZCLK.NE.0.AND.TOTAL(8).NE.0D0)
     &       RELATIVE(8)=DABS(ERROR*100.D0/TOTAL(8))

      IF(NZFOD.NE.0.AND.TOTAL(9).NE.0D0)
     &       RELATIVE(9)=DABS(ERROR*100.D0/TOTAL(9))

      IF(NZZOR.NE.0.AND.TOTAL(10).NE.0D0)
     &       RELATIVE(10)=DABS(ERROR*100.D0/TOTAL(10))


      WRITE(IUBALC,100)
      WRITE(IUBALC,1910) 
      WRITE(IUBALC,3100)
      WRITE(IUBALC,1080)
      WRITE(IUBALC,3300) (RELATIVE(I),I=1,11)
      WRITE(IUBALC,100)
      WRITE(IUBALC,1042)
      

    
 1710  FORMAT('# TPT. MASS BALANCE ERROR:',G14.6)
 1800  FORMAT('#                         --------------')
 1910  FORMAT('# TPT. MASS BALANCE ERROR (%) RELATIVE TO:')
 3100 FORMAT('#       POROSITY     PRS.HEAD     PRS.FLOW     LEAKAGE
     &      RECHARGE     PRS.CONC     MASS INP     CON LEAK      F.O.D.
     &       Z.O.D.      MAT.DIFF')
 3300  FORMAT(6X,11(G12.6,1X))


C-------------------- Writes tpt. mass balance for nodes

C-------------------- Nodal mass balance header.


      IF (IOBALNP.NE.0) THEN

          WRITE(IUBALC,100)
          WRITE(IUBALC,1000) I_PBL

          IF(INDBALTYP.EQ.0) THEN ! Steady state mass balance

              WRITE(IUBALC,2010) TIME(INTI+1)
              WRITE(IUBALC,2040)

          ELSE IF (INDBALTYP.EQ.1) THEN !Transient mass balance

              WRITE(IUBALC,2020) TIME(INTI+1)
              WRITE(IUBALC,2040)

          ELSE IF (INDBALTYP.EQ.2) THEN ! Global mass balance

              WRITE(IUBALC,2030)
              WRITE(IUBALC,2041)

          END IF !INDBALTYP.EQ.0,1,2

          WRITE(IUBALC,100)
          WRITE(IUBALC,2050)
          WRITE(IUBALC,2060)


 2010 FORMAT('# STEADY-STATE NODAL MASS BALANCE. TIME:',G14.6)

 2020 FORMAT('# TRANSIENT NODAL MASS BALANCE. TIME',G14.6)
 2030 FORMAT('# GLOBAL NODAL MASS BALANCE.')
 
 2040 FORMAT('# --------------------------------------------')
 2041 FORMAT('# --------------------------')
  
 2050 FORMAT('# ND    POROSITY     PRS.HEAD     PRS.FLOW     LEAKAGE
     &      RECHARGE     PRS.CONC     MASS INP     CON LEAK      F.O.D.
     &       Z.O.D.      MAT.DIFF     LAT.FLUX')
 2060 FORMAT('#',4X,12('[-----------]'))

          IFORMAT = 13

          DO I=1,NUMNP

              IF (IOBALNP.EQ.1.OR.(IOBALNP.EQ.2.AND.IBTCO(I).EQ.1)) THEN

                  LINE = ''
                  IINI = 1
                  IEND = 13

                  DO IPRO=1,12

                      WRITE(LINE(IINI:IEND),3410) BM_ND(I,IPRO)

                      IINI = IINI + IFORMAT
                      IEND = IEND + IFORMAT
            
                  END DO !IPRO=1,NMAXNP

                  WRITE(IUBALC,3500) I,LINE
              END IF 
          END DO !I=1,NUMNP                 
    

          WRITE(IUBALC,100)
          WRITE(IUBALC,1042)

      END IF !IOBALNP.NE.0

 3410 FORMAT(G12.6) 
 3500 FORMAT(I5,1X,A156)
      
      END SUBROUTINE BALANCE_WRITE_TR

