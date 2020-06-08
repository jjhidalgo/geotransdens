       SUBROUTINE WRI 
     ; (IOCMC    ,IOCMH    ,IOMHC    ,IOMHH    ,IOPLC
     ; ,IOPLH    ,IORTS    ,IOSMFL   ,IOSMTP   ,IOTRS
     ; ,IPROCESS ,LMXNDL   ,NINT     ,NPBFL    ,NPBTP    ,NUMEL
     ; ,NUMNP    ,NDEVS    ,CCALIT   ,NUMTIT
     ; ,NUMTOBS  ,FILENAME ,HCALIT
     ; ,VOBSC    ,VOBS     ,IBCOD    ,IBTCO    ,IXALF
     ; ,IXCHP    ,IXCON    ,IXCONT   ,IXQQP
     ; ,KXX      ,LNNDEL   ,LXARR    ,LXARRT   ,LXCRD ,LXDFM
     ; ,LXDSP    ,LXPOR    ,LXSTG    ,LXTRA    ,DEVNAME  ,TIME
     ; ,X        ,Y        ,Z        ,IOWAR    ,MAINF    ,IODEVICE ,TIT
     ; ,ISIM_GS  ,ISOLEQ)


******************************************************************
***   WRITE SOME RESULTS IN DIFFERENT FORMATS FOR GRAPHICS
******************************************************************

        IMPLICIT REAL*8 (A-H,O-Z)

        CHARACTER FILENAME(20)*20,NAMEFIL*20,DEVNAME(NDEVS)*10

        DIMENSION TIME(NINT),HCALIT(NUMNP),CCALIT(NUMNP),
     .            X(NUMNP),Y(NUMNP),Z(NUMNP)
     .            ,KXX(LMXNDL,NUMEL),IBCOD(NUMNP),
     .            IBTCO(NUMNP),IXCHP(NUMNP),
     .            IXQQP(NUMNP),IXALF(NUMNP),LNNDEL(NUMEL),
     .            LXTRA(NUMEL),LXSTG(NUMEL),LXARR(NUMEL),
     .            LXARRT(NUMEL),LXDFM(NUMEL),LXDSP(NUMEL),
     .            LXCRD(NUMEL),LXPOR(NUMEL),IXCON(NUMNP),
     .            IXCONT(NUMNP),
     .            VOBS(NUMTOBS),VOBSC(NUMTOBS+NDEVS),
     .            IODEVICE(NDEVS+1,9)

C------------------------- Format for program CUR

       DO IPROB=1,NPBFL
          IF (IOPLH.NE.0) CALL WR_IOPL 
     ;(81       ,IOPLH    ,1        ,IPROB    ,26      ,NDEVS
     ;,NUMTIT   ,NUMTOBS  ,IODEVICE ,TIT     ,VOBSC(NUMTOBS+1)
     ;,VOBS     ,VOBSC    ,DEVNAME)

       ENDDO

       DO IPROB=1,NPBTP
          IF (IOPLC.NE.0) CALL WR_IOPL 
     ;(81       ,IOPLC    ,2        ,IPROB    ,26      ,NDEVS
     ;,NUMTIT   ,NUMTOBS  ,IODEVICE ,TIT     ,VOBSC(NUMTOBS+1)
     ;,VOBS     ,VOBSC    ,DEVNAME)

       ENDDO

C-------------------------Geophysical observations controlled by IOPLC. 

       DO IPROB=1,NPBTP
          IF (IOPLC.NE.0) CALL WR_IOPL 
     ;(81       ,IOPLC    ,6        ,IPROB    ,26      ,NDEVS
     ;,NUMTIT   ,NUMTOBS  ,IODEVICE ,TIT     ,VOBSC(NUMTOBS+1)
     ;,VOBS     ,VOBSC    ,DEVNAME)

       ENDDO

C------------------------- Format for program MESHPLOT

       NCON=NUMNP+NUMEL-1                !Number of connections

C------------------------- Head

       IF (IOMHH .NE. 0) THEN

          INDEX=1
          NAMEFIL=FILENAME(13)
          IPROC=IPROCESS
          DO IPROB=1,MAX(1,NPBFL*IOSMFL)
             IF (IPROB.NE.1 .OR. ISIM_GS.GT.1) IPROC=1 ! Don't write MSH file
             CALL WR_VAR_MSH
     ;  (94       ,IBCOD    ,INDEX    ,IOMHH    ,1        ,IOSMFL
     ;  ,IOTRS    ,IPROB    ,IPROC    ,27       ,29       ,IXCHP
     ;  ,IXQQP    ,IXALF    ,KXX      ,LXTRA    ,LXSTG    ,LXARR
     ;  ,LXARRT   ,NAMEFIL  ,NCON     ,NINT     ,NUMEL    ,NUMNP
     ;  ,NDEVS    ,NPBFL    ,TIME     ,HCALIT   ,X        ,Y
     ;  ,ISOLEQ)

          END DO !IPROB=1,MAX(1,NPBFL*IOSMFL)

C------------------------ Output for Visualmeshplot
c-parche
c          DO IPROB=1,MAX(1,NPBFL*IOSMFL) !El bucle se repite porque se rebobina el fichero auxiliar
c
c              CALL WRITE_VISUALMESHPLOT
c     &            (FILENAME(1),94       ,INDEX      ,IOMHH    ,IOSMFL
c     &            ,IOTRS      ,IPROB    ,ISOLEQ     ,0        ,KXX
c     &            ,LMXNDL     ,LNNDEL   ,NINT       ,NUMEL    ,NPBFL
c     &            ,NUMNP      ,TIME     ,HCALIT     ,X          ,Y
c     &            ,Z)
c-fin-parche
c
c          END DO !IPROB=1,MAX(1,NPBFL*IOSMFL)


       END IF !IOMHH .NE. 0

C------------------------- Concentrations

       IF (IOMHC .NE. 0) THEN

          INDEX=2
          NAMEFIL=FILENAME(16)
          IPROC=IPROCESS
          DO IPROB=1,MAX(1,NPBTP*IOSMTP)
             IF (IPROB.NE.1 .OR. ISIM_GS.GT.1) IPROC=1 ! Don't write MSH file
             CALL WR_VAR_MSH
     ;  (93       ,IBTCO    ,INDEX    ,IOMHC    ,1        ,IOSMTP
     ;  ,IORTS    ,IPROB    ,IPROC    ,30       ,32       ,IXCON
     ;  ,IXCONT   ,IXCONT   ,KXX      ,LXDFM    ,LXDSP    ,LXPOR
     ;  ,LXCRD    ,NAMEFIL  ,NCON     ,NINT     ,NUMEL    ,NUMNP
     ;  ,NDEVS    ,NPBTP    ,TIME     ,CCALIT   ,X        ,Y
     ;  ,ISOLEQ)

          END DO !IPROB=1,MAX(1,NPBTP*IOSMTP)

C------------------------ Output for Visualmeshplot
c-parche
c          DO IPROB=1,MAX(1,NPBTP*IOSMTP) !El bucle se repite porque se rebobina el fichero auxiliar
c
c              CALL WRITE_VISUALMESHPLOT
c     &            (FILENAME(1),93       ,INDEX      ,IOMHC    ,IOSMTP
c     &            ,IORTS      ,IPROB    ,ISOLEQ     ,1        ,KXX
c     &            ,LMXNDL     ,LNNDEL   ,NINT       ,NUMEL    ,NPBTP
c     &            ,NUMNP      ,TIME     ,CCALIT     ,X          ,Y
c     &            ,Z)
c-fin-parche
c
c          END DO !IPROB=1,MAX(1,NPBTP*IOSMTP)


       END IF !IOMHC .NE. 0

       IF (IOCMH.NE.0) THEN
         DO IPROB=1,NPBFL
        
           CALL WR_MSVSC
     ;(1        ,IOWAR    ,IPROB    ,34       ,MAINF    ,NDEVS
     ;,NUMTOBS  ,IODEVICE ,VOBS     ,VOBSC)

         END DO

       END IF

       IF (IOCMC.NE.0) THEN
         DO IPROB=1,NPBTP
        
           CALL WR_MSVSC
     ;(2        ,IOWAR    ,IPROB    ,34       ,MAINF    ,NDEVS
     ;,NUMTOBS  ,IODEVICE ,VOBS     ,VOBSC)

         END DO

         DO IPROB=1,NPBTP
        
           CALL WR_MSVSC
     ;(6        ,IOWAR    ,IPROB    ,34       ,MAINF    ,NDEVS
     ;,NUMTOBS  ,IODEVICE ,VOBS     ,VOBSC)

           END DO
       END IF

       RETURN
       END

************************************************************************
************************************************************************

       SUBROUTINE WR_VAR_MSH 
     ;  (IAUX     ,IBCC     ,INDEX    ,IOMVV    ,IOP      ,IOSMLT
     ;  ,IOTIM    ,IPROB    ,IPROCESS ,IUGR1    ,IUGR2    ,IX1
     ;  ,IX2      ,IX3      ,KXX      ,LX1      ,LX2      ,LX3
     ;  ,LX4      ,NAMEFIL  ,NCON     ,NINT     ,NUMEL    ,NUMNP
     ;  ,NUOBS    ,NUMPB    ,TIME     ,VCALIT   ,X        ,Y
     ;  ,ISOLEQ)

*****************************************************************************
*
* PURPOSE 
*    
*     Write Meshplot file
*
* DESCRIPTION 
*    
*     Write Meshplot file
*
*
* EXTERNAL VARIABLES: SCALARS
*
*  INDEX                  =1 flow. =2 transport.
*
* HISTORY: 
*           AMS 
*
********************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)

       DIMENSION CUV(100),VCALIT(NUMNP),TIME(NINT),ISOLEQ(NINT,4)
       CHARACTER NAMEVAR(10)*15,NAMEFIL*20
       DATA NAMEVAR(1)/' HEADS         '/
     ;     ,NAMEVAR(2)/' CONCENTRATION '/

C------------------------- Writes MSH file

          IF (IPROCESS.EQ.0)  
     ;       CALL WRIMHF(IUGR1,NAMEFIL,NUMNP,NUMEL,
     ;            NUOBS,NCON,X,Y,IBCC,IX1,IX2,IX3,
     ;            KXX,LX1,LX2,LX3,LX4,IOP)

C------------------------- Starts at the beginning of file IAUX

          REWIND (IAUX)

C------------------------- Reads the number of contour lines and their value

          NCUV=-1  !previous assignation just in case the file DIM is finished
          READ(10,1000,END=9000) NCUV
 1000     FORMAT(I5)
 9000     IF (NCUV.GT.0) READ(10,1100) (CUV(I),I=1,NCUV)
 1100     FORMAT(8F10.0)
           
C------------------------- NUM_REG is the number of records to skip 

          NUM_REG=MAX(1,IPROB*IOSMLT)

          NTOT_REG=NINT
          IF (ISOLEQ(1,1).EQ.1 .AND. IOTIM.NE.0) THEN
             INI=2
          ELSE
             INI=1
          ENDIF
*          IF (IOTIM.EQ.0) NTOT_REG=1
          DO IT=INI,NTOT_REG
             JT=MAX(1,IT-1)  ! Should be equal to IT-1 except for
                             ! IT=0 to avoid accessing zero element of ISOLEQ
             IF ( (ISOLEQ(JT,INDEX).EQ.0 .AND. IOTIM.GT.0 .AND.    !TRANSIENT
     ;               ( (MOD(IT,IOMVV).EQ.0 .AND. IT.NE.1) .OR.     !TRANSIENT
     ;                                       IT.EQ.NINT) ) .OR.    !TRANSIENT
     ;          (IT.EQ.1 .AND. (IOTIM.EQ.0 .OR. IOTIM.EQ.2) ) .OR. ! STEADY
     ;           (MOD(IT,IOMVV).EQ.0 .AND. ISOLEQ(JT,INDEX).EQ.1) ) THEN

                DO IP=1,NUM_REG
                   READ(IAUX) VCALIT
                ENDDO

C------------------------- Next times NUM_REG is the number of problems (if 
C------------------------- these are simultaneous)

                NUM_REG=MAX(1,NUMPB*IOSMLT)

                IF (IT.EQ.1 .OR. ISOLEQ(JT,INDEX).EQ.1) THEN
                   WRITE(IUGR2,2000) NAMEVAR(INDEX)
                ELSE
                   WRITE(IUGR2,2100) NAMEVAR(INDEX),IT,TIME(IT)
                ENDIF
 2000           FORMAT(' STEADY-STATE',A15)
 2100           FORMAT(A15,' AT',I5,'-th OBSERVATION TIME. t=',1P,G12.5)

                DO I=1,NUMNP
                   IF (DABS(VCALIT(I)).LE.1.D-100 .AND. 
     ;                                     VCALIT(I).NE.0.D0) THEN
                      WRITE (IUGR2,2250) I,VCALIT(I)
                   ELSE
                      WRITE (IUGR2,2200) I,VCALIT(I)
                   ENDIF
                ENDDO
                WRITE(IUGR2,1000) NCUV
                IF (NCUV.GT.0) WRITE(IUGR2,2300) (CUV(I),I=1,NCUV)
             ENDIF
c-parche PLT
c 2200        FORMAT(I5,G20.14) !Antes G15.8. Aumentado para diferencias finitas.
 2200        FORMAT(I5,G15.8)
c-fin-parche
 2250        FORMAT(I5,G15.7E3)
 2300        FORMAT(F10.3)
          ENDDO

       RETURN
       END

************************************************************************
************************************************************************

        SUBROUTINE WRIMHF(IU,FILE,NUMNP,NUMEL,NUOBS,NCON,X,Y,
     .                IB,IX1,IX2,IX3,KXX,LX1,LX2,LX3,LX4,IOP)
        IMPLICIT REAL*8 (A-H,O-Z)
        CHARACTER*20 FILE
        DIMENSION X(NUMNP),Y(NUMNP),IB(NUMNP),IX1(NUMNP),IX2(NUMNP),
     .            IX3(NUMNP),KXX(3,NUMEL),LX1(NUMEL),LX2(NUMEL),
     .            LX3(NUMEL),LX4(NUMEL)
        NPMAX=0
        THETA=0
        WRITE(IU,2000) FILE
 2000   FORMAT(A20)
        WRITE(IU,2100)
 2100   FORMAT(/////)
        WRITE(IU,1000) NUMNP,NUMEL,NUOBS,NPMAX,NCON
 1000   FORMAT(5I5)
        WRITE(IU,3000) THETA
 3000   FORMAT(F10.3)
        IF (IOP.EQ.0) THEN
           WRITE(IU,4000) (I,X(I),Y(I),IB(I),IX1(I),0,
     ;           0,I=1,NUMNP)
        ELSE
           WRITE(IU,4000) (I,X(I),Y(I),IB(I),IX1(I),IX2(I),
     ;           IX3(I),I=1,NUMNP)
        END IF
 4000   FORMAT(I5,2g10.4,2I2,8X,I2,8X,I2)
        WRITE(IU,4500) (L,KXX(1,L),KXX(2,L),KXX(3,L),KXX(3,L),
     .                   LX1(L),LX2(L),LX3(L),LX4(L),L=1,NUMEL)
 4500   FORMAT(9I4)
        RETURN
        END

************************************************************************
************************************************************************

      SUBROUTINE WR_IOPL
     ;(IAOB     ,IOPL     ,IOPTION  ,IPROB    ,IUGR     ,NDEVS
     ;,NUMTIT   ,NUMTOBS  ,IODEVICE ,TIT      ,VAUX     ,VOBS     
     ;,VOBSC    ,DEVNAME) 

C______________________________ Variables declaration

      IMPLICIT NONE
      INTEGER*4 IOPTION,ND,NDEVS,NMEAS,NCURV,NOF,NOL,NO,NUMTOBS,
     ;          IUGR,IOPL,NPOINT,NUMTIT,IREAD,IPROB,IAOB,J
      INTEGER*4 IODEVICE(NDEVS+1,9)
      REAL*8 VOBS(NUMTOBS),TIT(NUMTIT),VOBSC(NUMTOBS),
     ;       SIMTIME,VAUX(NDEVS)
      CHARACTER PROB*4,SVAR*6,DEVNAME(NDEVS)*10

C______________________________ Step 1: Identifies problem type (flow/tpt) and 
C______________________________ creates axis.

      IF (IOPTION.EQ.1) THEN                                      ! Flow problem
        PROB='FLOW'        
        SVAR=' HEAD'
      ELSE IF (IOPTION.EQ.2) THEN                                 ! Tpt. problem
        PROB='TPT.'
        SVAR=' CONC.'
      ELSE IF (IOPTION.EQ.3) THEN                                 ! Tpt. problem
        CONTINUE      
      ELSE IF (IOPTION.EQ.4) THEN                                 ! Tpt. problem
        CONTINUE      
      ELSE IF (IOPTION.EQ.6) THEN                                 ! Geophysical data. 
        PROB='TPT.'
        SVAR=' RESIS'  
      END IF

C______________________________ Step 2: Loop over observation devices

      DO ND=1,NDEVS

C______________________________ Step 2.1: Checks if device is taken into account

        IF (IODEVICE(ND,1).EQ.IOPTION .AND.                     ! Device type ok
     ;      IODEVICE(ND,9).EQ.IPROB.AND.     ! Dev. corresponds to current prob.
     ;      IODEVICE(ND,2).GT.0) THEN                  ! Dev. taken into account

C______________________________ Step 2.2: Checks if measurements are available

            NMEAS=0                     ! Number of valid meas. for current dev.
            NCURV=1                               ! Number of curves to be drawn
            NOF=IODEVICE(ND,8)              ! First observation for current dev.
            NOL=IODEVICE(ND+1,8)-1           ! Last observation for current dev.

            DO NO=NOF,NOL                        ! Loop over device observations
              NMEAS=NMEAS+1
            END DO                                            ! Next observation
            IF (NMEAS.GT.0) NCURV=2              

C______________________________ Step 2.3: Writes main header of PLT file (CUR3)

            WRITE(IUGR,1000) NCURV,PROB,IPROB
            
C______________________________ Step 2.4: Writes at measurement times
            
            IF (IOPL.EQ.1) THEN                                    ! Meas. times
              NPOINT=NOL-NOF+1                          ! Number of observations
              
C______________________________ Step 2.4.1: Writes header (CUR3) and axis

              WRITE(IUGR,1100) NPOINT,SVAR,DEVNAME(ND),SVAR

C______________________________ Step 2.4.2: Loop over device observations

              DO NO=NOF,NOL
                 WRITE(IUGR,1200) TIT(NO),VOBSC(NO)
              END DO
                
C______________________________ Step 2.5: Writes at simulation times

            ELSE                                              ! Simulation times

C______________________________ Step 2.5.1: Checks for number of simulations
              
              REWIND (IAOB)
              DO IREAD=1,111111 !Was 11,111 changed to 111,111 ¿?
                                !Changed to 10000 so that NPOINT.LE.99999,i.e,
                                !maximun possible number of time intervals.
                READ(IAOB,END=666) SIMTIME,(VAUX(J),J=1,NDEVS)
              END DO            
 666          NPOINT=IREAD-1

C______________________________ Step 2.5.2: Writes header (CUR3) and axis

              WRITE(IUGR,1300) NPOINT,SVAR,DEVNAME(ND),SVAR
              
C______________________________ Step 2.5.3: Reads again, reordenates and writes

              REWIND(IAOB)
              DO IREAD=1,NPOINT
                READ(IAOB) SIMTIME,(VAUX(J),J=1,NDEVS)
                WRITE(IUGR,1200) SIMTIME,VAUX(ND)
              END DO
                
            END IF                              ! Writting at meas./simul. times

C______________________________ Step 2.6: Writes (if available) measurements

            IF (NMEAS.GT.0) THEN                               ! Meas. available

C______________________________ Step 2.6.1: Writes header (CUR3) and axis

              WRITE(IUGR,1400) NMEAS,SVAR,DEVNAME(ND),SVAR


C______________________________ Step 2.6.2: Loop over device measurements

              DO NO=NOF,NOL
                WRITE(IUGR,1200) TIT(NO),VOBS(NO)
              END DO

            END IF                                         ! Meas. not available

        END IF                                       ! Device must be considered

      END DO                                                       ! Next device

 1000 FORMAT(I5,'   (2F13.0)                             ',
     ;       A4,' PROBLEM: ',I2)
 1100 FORMAT(I5,/,' COMPUTED ',A6,' OBS. DEVICE: ',A10,/,
     ;            ' TIME',/,A6)
c-parche-cambiado formato para PEST (2G20.14) en el PLT
 1200 FORMAT(2G13.6)
c 1200 FORMAT(2G20.14)
C-fin-parche
 1300 FORMAT(I5,/,' COMPUTED ',A6,' OBS. DEVICE: ',A10,/,
     ;            ' TIME',/,A6)
 1400 FORMAT(I5,3X,'-4',/,' MEASURED ',A6,' OBS. DEVICE: ',A10,/,
     ;            ' TIME',/,A6)
      RETURN
      END

************************************************************************
************************************************************************

      SUBROUTINE WR_MSVSC
     ;(IOPTION  ,IOWAR    ,IPROB    ,IUGR     ,MAINF    ,NDEVS
     ;,NUMTOBS  ,IODEVICE ,VOBS     ,VOBSC)

      IMPLICIT NONE
      
C______________________________ Step 0: Declaration of variables

      INTEGER*4 ND,NDEVS,NUMTOBS,IOPTION,IUGR,NO,NOF,NOL,NCURV,IPROB,
     ;          NMEAS,MAINF,IOWAR
      INTEGER*4 IODEVICE(NDEVS+1,9)
      REAL*8 VOBSC(NUMTOBS),VOBS(NUMTOBS)
      CHARACTER PROB*4,SVAR*6

C_____________________ Step 1: Determination of the number of valid observations
C_____________________         corresponding to current problem

      NMEAS=0                     ! Number of valid meas. for current dev.

      DO ND=1,NDEVS                                          ! Loop over devices

C______________________________ Step 1.1: Checks if device is taken into account

        IF (IODEVICE(ND,1).EQ.IOPTION .AND.                     ! Device type ok
     ;      IODEVICE(ND,9).EQ.IPROB.AND.     ! Dev. corresponds to current prob.
     ;      IODEVICE(ND,2).GT.0) THEN                  ! Dev. taken into account

C______________________________ Step 1.2: Checks if measurements are available

            NCURV=1                               ! Number of curves to be drawn
            NOF=IODEVICE(ND,8)              ! First observation for current dev.
            NOL=IODEVICE(ND+1,8)-1           ! Last observation for current dev.

            DO NO=NOF,NOL                        ! Loop over device observations
              NMEAS=NMEAS+1
            END DO                                            ! Next observation

        END IF                                                   ! Device is ok!

      END DO                                                       ! Next device

C______________________________ Step 2: Identifies problem type (flow/tpt) and
C______________________________ creates axis.

      IF (IOPTION.EQ.1) THEN                                      ! Flow problem
        PROB='FLOW'
        SVAR=' HEAD'
      ELSE IF (IOPTION.EQ.2) THEN                                 ! Tpt. problem
        PROB='TPT.'
        SVAR=' CONC.'
      ELSE IF (IOPTION.EQ.3) THEN                                 ! Tpt. problem
        CONTINUE
      ELSE IF (IOPTION.EQ.4) THEN                                 ! Tpt. problem
         CONTINUE
      ELSE IF (IOPTION.EQ.6) THEN                                 ! Geophysical data. 
        PROB='TPT.'
        SVAR=' RESIS'  
      END IF

C______________________________ Step 3: If measurements are available, 
C______________________________         everything is OK; else, a warning is
C______________________________         written and leaves the subroutine

      IF (NMEAS.GT.0) THEN
        NCURV=1
      ELSE

       IF(IOWAR.GT.0) WRITE(MAINF,1000) SVAR,PROB,IPROB
 1000  FORMAT(//,' WARNING: NO MEASUREMENTS AVAILABLE. ',A6, ' - ',A4,
     ;           ' PROBLEM: ',I5,//)
       RETURN

      END IF

C______________________________Step 4: Writes header and axis (CUR3 format)

      WRITE(IUGR,1100) NCURV,PROB,IPROB
 1100 FORMAT(I5,'   (2F13.0)                             ',
     ;       A4,' PROBLEM: ',I2)
      WRITE(IUGR,1200) NMEAS,SVAR,SVAR,PROB,IPROB,SVAR,SVAR
 1200 FORMAT(I5,3X,'-4',/,' COMPUTED',A6,'VS MEASURED',A6,A4,
     ;      ' PROBLEM: ',I2/,' COMPUTED',A6,/,' MEASURED',A6)


C______________________________ Step 5: Loop over observations

      DO ND=1,NDEVS
        IF (IODEVICE(ND,1).EQ.IOPTION .AND.                     ! Device type ok
     ;      IODEVICE(ND,9).EQ.IPROB.AND.     ! Dev. corresponds to current prob.
     ;      IODEVICE(ND,2).GT.0) THEN                  ! Dev. taken into account

          NOF=IODEVICE(ND,8)                ! First observation for current dev.
          NOL=IODEVICE(ND+1,8)-1             ! Last observation for current dev.

          DO NO=NOF,NOL                          ! Loop over device observations
            WRITE(IUGR,1300) VOBSC(NO),VOBS(NO)
          END DO                                              ! Next observation
          
        END IF
      END DO                                                       ! Next device
 1300 FORMAT(2G13.6)
      RETURN
      END

************************************************************************
************************************************************************

       SUBROUTINE WRITE_ARRAY(WW,TITULO,NNN)
       IMPLICIT REAL*8 (A-H,O-Z)
       COMMON /FILENUM/ MAINF,IUDIM,IUGRID,IUPAR,IUTIM,IUOBS,
     .                  IUCAL,IUGR1,IUGR2,IUGR3,IUGR4,IUGR5,IUGR6,
     .                  IUGR7,IAUXH,IAUXC
       CHARACTER*80 TITULO,TITULO1,BLANK
       DIMENSION WW(NNN)
       DATA BLANK/'
     ;                           '/
       LONG=INDEX(TITULO,'@')
       TITULO1=BLANK
       TITULO1(1:LONG-1)=TITULO(1:LONG-1)
       WRITE(MAINF,2000) TITULO1
 2000  FORMAT(//15X,A80)
       K1=NNN/10
       K11=NNN/100
       K21=K1-K11*10
       K2=NNN-K1*10
       WRITE(MAINF,1100) 0
 1100  FORMAT(/,60X,'/',I5,'/')
       DO 100 I=1,K11
          DO 110 J=1,10
 110         WRITE (MAINF,1000) (WW((I-1)*100+(J-1)*10+L),L=1,10)
 1000  FORMAT(10D13.5)
 100      WRITE(MAINF,1100) I*100
       K3=K11*100
       DO 120 J=1,K21
 120      WRITE (MAINF,1000) (WW(K3+(J-1)*10+L),L=1,10)
       K4=K3+K21*10
       WRITE (MAINF,1000) (WW(K4+L),L=1,K2)
       RETURN
       END

************************************************************************
************************************************************************

       SUBROUTINE WRITE_MATRIX(XMM,ND1,ND2,TITULO)
       IMPLICIT REAL*8 (A-H,O-Z)
       COMMON /FILENUM/ MAINF,IUDIM,IUGRID,IUPAR,IUTIM,IUOBS,
     .                  IUCAL,IUGR1,IUGR2,IUGR3,IUGR4,IUGR5,IUGR6,
     .                  IUGR7,IAUXH,IAUXC
       CHARACTER*80 TITULO,TITULO1,BLANK
       DATA BLANK/'
     ;                           '/
       DIMENSION XMM(ND1,ND2)
       LONG=INDEX(TITULO,'@')
       TITULO1=BLANK
       TITULO1(1:LONG-1)=TITULO(1:LONG-1)
       WRITE(MAINF,1000) TITULO1
 1000  FORMAT(//,15X,A80)
       K1=ND2/10
       K3=0
       DO 100 I=1,ND1
        WRITE (MAINF,3000) I
 3000   FORMAT(/,35X,' ROW ',I3,/,' ')
        DO 110 L=1,K1
           K2=(L-1)*10+1
           K3=L*10
 110       WRITE(MAINF,2000) (XMM(I,J),J=K2,K3)
 100    WRITE(MAINF,2000) (XMM(I,J),J=K3+1,ND2)
 2000  FORMAT(10D13.5)
       RETURN
       END

************************************************************************
************************************************************************

       SUBROUTINE WRSYMMAT(A,N,IFILE)
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION A(N*(N+1)/2)
       
       INITROW=1
       IENDROW=1
       LONG=1
       DO I=1,N
          WRITE(IFILE,1000) (A(K),K=INITROW,IENDROW)
          INITROW=IENDROW+1
          IENDROW=INITROW+LONG
          LONG=LONG+1
       END DO
 1000  FORMAT(10E13.3)
       RETURN
       END

************************************************************************
************************************************************************

