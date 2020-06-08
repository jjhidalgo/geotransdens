      SUBROUTINE WRITE_VELOCITY_VMSH(DIMFILE  ,INTI     ,IODIM
     &                               ,IOFIRST ,KXX      ,LMXNDL
     &                               ,LNNDEL   ,NINT     ,NUMEL
     &                               ,NUMNP    ,TIME     ,VD
     &                               ,X        ,Y)

*****************************************************************************
*
* PURPOSE 
*    
*     Write Visual Meshplot file with the velocities
*
* DESCRIPTION 
*    
*     Write Visual Meshplot file with the velocities
*
*
* EXTERNAL VARIABLES: SCALARS
*
*
* HISTORY: 
*           JHG
*
********************************************************************************

       IMPLICIT NONE

       INTEGER*4::IROOTLEN,NUMNP,NNUD,INODE,L,LENVARUNITS
     &           ,NINT,NUMEL,LMXNDL,I,IWIDTH,INTI,IODIM,IOFIRST

      INTEGER*4::LNNDEL(NUMEL),KXX(LMXNDL,NUMEL)

      REAL*8::VD(IODIM, NUMEL),TIME(NINT),X(NUMNP),Y(NUMNP)

      CHARACTER::VMSHFILE*20,DIMFILE*20,TYPEL*5,VARUNITS*19,STEPN*5
      CHARACTER::BLANK*1
      CHARACTER::strFMt50*40,strFMt60*20,strFMt90*20,strFMt100*20
      DATA BLANK/' '/

C------------------------- Opens file. Takes root from DIM.DAT file.

      IROOTLEN=INDEX(DIMFILE,'DIM',BACK=.TRUE.)-1
      VMSHFILE = DIMFILE(1:IROOTLEN) //'VMSHV.INP'
      VARUNITS = "Darcy's Flow, [L/T]"
      LENVARUNITS = 19

      IF (IOFIRST.EQ.0) THEN

          OPEN(UNIT=778,FILE=VMSHFILE,STATUS='UNKNOWN')

C------------------------- Writes INP header
C------------------------- Only data is suppose to vary in time	 
          WRITE(778,10) NINT-1
   10     FORMAT('# VISUAL MESHPLOT FILE',/,I5,/,'data')

C------------------------- Mesh is included in the first step


C------------------------- Mesh is written

          WRITE(778,20) 'step'//'1',TIME(1)
   20     FORMAT(A9,' ',G15.5)

C------------------------- Mesh Dimensions

          WRITE(778,30) NUMNP,NUMEL
   30     FORMAT(I5,' ',I5)

C------------------------- Mesh Nodes

          DO INODE=1,NUMNP
              WRITE(778,40) INODE,X(INODE),Y(INODE)
   40         FORMAT(I5,' ',G15.5,' ',G15.5)
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
              strFmt50 = Trim(AdjustL(strFmt50))//'(I5,1X))'

              WRITE(778,strFmt50) L,1,TYPEL,KXX(1:NNUD,L)
c  50         FORMAT(I5, ' ',I5,' ',A5,' ',<NNUD>(I5,' '))

          END DO !L

      END IF ! IOFIRST.EQ.0

C------------------------- Steps data.


      IF (INTI.GT.1) THEN
          IWIDTH = INT(LOG10(1D0*INTI))+1 ! 1d0*ISTEP para que sea real el argumento de LOG10

          WRITE(STEPN,'(I5)') INTI
          STEPN = Trim(AdjustL(STEPN))

          strFmt60 = ''
          WRITE (strFmt60,*) IWIDTH+4
          strFmt60 = '(A'//Trim(AdjustL(strFmt60))//',1X,G15.5)'

          WRITE(778,strFmt60) 'step'//TRIM(STEPN),TIME(INTI)
c  60     FORMAT(A<IWIDTH+4>,' ',G15.5)

      END IF
C------------------------- Number of nodal and element variables

      WRITE(778,70) 0,2
   70 FORMAT(I5,' ',I5)

C------------------------- Number of components and component size

      WRITE(778,80) 1,2
   80 FORMAT(I5,' ',I5)

C------------------------- Kind of variable and units
      strFmt90 = ''
      WRITE (strFmt90,*) LENVARUNITS
      strFmt90 = '(A'//Trim(AdjustL(strFmt90))//')'

      WRITE(778,strFmt90) VARUNITS

c  90 FORMAT(A<LENVARUNITS>)

C------------------------- Velocities for each element

      strFmt100 = ''
      WRITE (strFmt100,*) IODIM
      strFmt100 = '(I5,1X,'//Trim(AdjustL(strFmt100))//'(G15.5,A1))'

      DO L=1,NUMEL
          WRITE(778,strFmt100) L,(VD(I,L),BLANK,I=1,IODIM)
c 100    FORMAT(I5,' ',<IODIM>(G15.5,A1))
      END DO

      IF (INTI.EQ.NINT) THEN
          CLOSE(778)
      END IF

      END SUBROUTINE WRITE_VELOCITY_VMSH
