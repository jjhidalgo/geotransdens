      SUBROUTINE WRITE_VISUALMESHPLOT
     &          (DIMFILE  ,IAUX     ,IND      ,IOMVV    ,IOSMLT
     &          ,IOTIM    ,IPROB    ,ISOLEQ   ,ISTATEVAR,KXX
     &          ,LMXNDL   ,LNNDEL   ,NINT     ,NUMEL    ,NUMPB
     &          ,NUMNP    ,TIME     ,VCALIT   ,X        ,Y
     &          ,Z)

*****************************************************************************
*
* PURPOSE
*
*     Write Visual Meshplot file
*
* DESCRIPTION
*
*     Write Visual Meshplot file
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  IWRTOPT                1. Flow
*                         2. Flow and transport.
*  IORTS        Transport regime.
*                 0. Steady state transport.
*                 1. Transient transport with prescribed initial conditions.
*                 2. Transient transport with steady-state initial conditions.
*
*  ISOLEQ       Array containing the type of head/concentration
*               solution desired for the user at each obs. time
*                 0. Transient
*                 1. Steady state.
*                 2. Read initial conditions
*                 3. Null initial conditions
*                 4. Not solve.
*
* HISTORY:
*           JHG
*
********************************************************************************

       IMPLICIT NONE

C------------------------- External

      INTEGER*4::IAUX     ,IND      ,IPROB    ,IOMVV    ,IOSMLT   ,IOTIM
     &          ,NUMNP    ,NINT     ,ISTATEVAR,LMXNDL   ,NUMEL
     &          ,NUMPB

      INTEGER*4::ISOLEQ(NINT,4),LNNDEL(NUMEL),KXX(LMXNDL,NUMEL)

      REAL*8::VCALIT(NUMNP),TIME(NINT),X(NUMNP),Y(NUMNP),Z(NUMNP)

      CHARACTER::DIMFILE*20


C------------------------- Internal

      INTEGER*4::I          ,INI        ,INI_IT     ,INODE
     &          ,IP         ,IROOTLEN   ,ISTEP      ,IT
     &          ,IWIDTH     ,JT         ,L          ,LENVARUNITS
     &          ,NNUD       ,NSTEPS     ,NTOT_REG   ,NUM_REG


      CHARACTER::VMSHFILE*20,TYPEL*5,VARUNITS*16,STEPN*5,PB*5
      CHARACTER::strFmt50*40,strFmt60*20,strFmt90*20

C------------------------- Rewinds auxiliar files
      REWIND(IAUX)

C------------------------- Opens file. Takes root from DIM.DAT file.

      IROOTLEN=INDEX(DIMFILE,'DIM',BACK=.TRUE.)-1

       SELECT CASE (IPROB)

          CASE(1)
                PB=''
          CASE(2:9)
              WRITE(PB,1) IPROB
    1         FORMAT('0000',I1)
            CASE(10:99)
              WRITE(PB,2) IPROB
    2         FORMAT('000',I2)
            CASE(100:999)
              WRITE(PB,3) IPROB
    3         FORMAT('00',I3)
            CASE(1000:9999)
              WRITE(PB,4) IPROB
    4         FORMAT('0',I4)
            CASE(10000:99999)
                WRITE(PB,5) IPROB
    5         FORMAT(I5)
        END SELECT !NUMPB

      IF(ISTATEVAR.GT.0) THEN
          VMSHFILE = DIMFILE(1:IROOTLEN) //TRIM(PB)//'VCC.INP'
          VARUNITS = "Mass Fraction, -"
          LENVARUNITS = 16
      ELSE
          VMSHFILE = DIMFILE(1:IROOTLEN) //TRIM(PB)//'VHH.INP'
          VARUNITS = "Head, m"
          LENVARUNITS = 7
      END IF


      OPEN(UNIT=777,FILE=VMSHFILE,STATUS='UNKNOWN')

C------------------------- Counts Steps.

      NSTEPS = 0
      NTOT_REG=NINT
      IF (ISOLEQ(1,1).EQ.1 .AND. IOTIM.NE.0) THEN
         INI=2
      ELSE
         INI=1
      ENDIF

      IF (IOTIM.EQ.0) NTOT_REG=1
      DO IT=INI,NTOT_REG
         JT=MAX(1,IT-1)  ! Should be equal to IT-1 except for
                         ! IT=0 to avoid accessing zero element of ISOLEQ
         IF ( (ISOLEQ(JT,IND).EQ.0 .AND. IOTIM.GT.0 .AND.    !TRANSIENT
     &           ( (MOD(IT,IOMVV).EQ.0 .AND. IT.NE.1) .OR.     !TRANSIENT
     &                                   IT.EQ.NINT) ) .OR.    !TRANSIENT
     &      (IT.EQ.1 .AND. (IOTIM.EQ.0 .OR. IOTIM.EQ.2) ) .OR. ! STEADY
     &       (MOD(IT,IOMVV).EQ.0 .AND. ISOLEQ(JT,IND).EQ.1) ) THEN

            NSTEPS = NSTEPS + 1

            IF (NSTEPS.EQ.1) INI_IT = IT

C------------------------- Next times NUM_REG is the number of problems (if
C------------------------- these are simultaneous)

            NUM_REG=MAX(1,NUMPB*IOSMLT)

          END IF

      END DO !IT=INI,NTOT_REG

C------------------------- Writes INP header
C------------------------- Only data is suppose to vary in time
      WRITE(777,10) NSTEPS
10    FORMAT('# VISUAL MESHPLOT FILE',/,I5,/,'data')

C------------------------- Mesh is included in the first step


C------------------------- Mesh is written

      WRITE(777,20) 'step'//'1',time(INI_IT)
20    FORMAT(A9,' ',G15.5)

C------------------------- Mesh Dimensions
      WRITE(777,30) NUMNP,NUMEL
30    FORMAT(I5,' ',I5)

C------------------------- Mesh Nodes
      DO INODE=1,NUMNP
          WRITE(777,40) INODE,X(INODE),Y(INODE),Z(INODE)
40        FORMAT(I5,' ',G15.5,' ',G15.5,' ',G15.5)
      END DO !INODE

C------------------------- Mesh conectivities

      DO L=1,NUMEL
          NNUD=LNNDEL(L)

          SELECT CASE (NNUD)
              CASE(2)
                  TYPEL='line'
              CASE(3)
                  TYPEL='tri'
              CASE(4)
                  TYPEL='quad'
              CASE(6)
                  TYPEL='prism'
          END SELECT

          strFmt50 = ''
          WRITE (strFmt50,*) NNUD
          strFmt50 = '(I5,1X,I5,1X,A5,1X,'//Trim(AdjustL(strFmt50))
          strFmt50 = Trim(AdjustL(strFmt50)) // '(I5,1X))'

          WRITE(777,strFmt50) L,1,TYPEL,KXX(1:NNUD,L)
c50        FORMAT(I5, ' ',I5,' ',A5,' ',<NNUD>(I5,' '))

      END DO !L


C------------------------- Steps Loop.

      NUM_REG=MAX(1,IPROB*IOSMLT)
      NTOT_REG=NINT
      IF (IOTIM.EQ.0) NTOT_REG=1

      ISTEP = 0

      DO IT=INI,NTOT_REG
         JT=MAX(1,IT-1)  ! Should be equal to IT-1 except for
                         ! IT=0 to avoid accessing zero element of ISOLEQ
         IF ( (ISOLEQ(JT,IND).EQ.0 .AND. IOTIM.GT.0 .AND.    !TRANSIENT
     &           ( (MOD(IT,IOMVV).EQ.0 .AND. IT.NE.1) .OR.     !TRANSIENT
     &                                   IT.EQ.NINT) ) .OR.    !TRANSIENT
     &      (IT.EQ.1 .AND. (IOTIM.EQ.0 .OR. IOTIM.EQ.2) ) .OR. ! STEADY
     &       (MOD(IT,IOMVV).EQ.0 .AND. ISOLEQ(JT,IND).EQ.1) ) THEN

            ISTEP = ISTEP + 1

            DO IP=1,NUM_REG
              READ(IAUX) VCALIT
            ENDDO

C------------------------- Next times NUM_REG is the number of problems (if
C------------------------- these are simultaneous)

            NUM_REG=MAX(1,NUMPB*IOSMLT)

C------------------------- Step header (header for step is already written)
              IF (ISTEP.GT.1) THEN
                  IWIDTH = INT(LOG10(1D0*ISTEP))+1 !1d0*ISTEP para que sea real el argumento de LOG10


                  WRITE(STEPN,'(I5)') ISTEP
                  STEPN = TRIM(ADJUSTL(STEPN))

                  strFmt60 = ''
                  WRITE (strFmt60,*) IWIDTH+4
                  strFmt60 ='(A'//Trim(AdjustL(strFmt60))//',1X,G15.5)'

                  WRITE(777,strFmt60) 'step'//TRIM(stepn),TIME(IT)
c60                FORMAT(A<IWIDTH+4>,' ',G15.5)


              END IF !ISTEP.GT.1
C------------------------- Number of nodal and element variables
              WRITE(777,70) 1,0
70            FORMAT(I5,' ',I5)

C------------------------- Number of components and component size
              WRITE(777,80) 1,1
80            FORMAT(I5,' ',I5)

C------------------------- Kind of variable and units

              strFmt90 = ''
              WRITE (strFmt90,*) LENVARUNITS
              strFmt90 = '(A'//Trim(AdjustL(strFmt90))//')'

              WRITE(777,strFmt90) VARUNITS
c90            FORMAT(A<LENVARUNITS>)

C------------------------- State Variable for current step

              DO I=1,NUMNP

C------------------------- To avoid three digits exponents
C------------------------- (1.0E+123 is otherwise written as 1.0+123).

                  IF (DABS(VCALIT(I)).LE.1.D-100
     &               .AND. VCALIT(I).NE.0.D0) THEN

                      WRITE (777,2250) I,VCALIT(I)
 2250                 FORMAT(I5,' ',G15.7E3)

                  ELSE

                      WRITE (777,2200) I,VCALIT(I)
c 2200                 FORMAT(I5,' ',G15.8)
 2200                 FORMAT(I5,' ',G20.14)

                  END IF !DABS(VCALIT(I)).LE.1.D-100 ...


                END DO !I=1,NUMNP


          END IF
      END DO !IT=INI,NTOT_REG

      CLOSE (777)

      END SUBROUTINE WRITE_VISUALMESHPLOT
