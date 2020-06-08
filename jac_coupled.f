      SUBROUTINE JAC_COUPLED
     &          (A_COUPL_DSC        ,A_COUPL_DSCF       ,BCOUPLED
     &          ,DERC     ,DERH     ,IA_COUPLED_DSC_COLS
     &          ,IA_COUPLED_DSC_ROWS,IAD_D    ,IADD_D   ,IADN_D
     &          ,IAFD_D  ,IAFDD_D   ,IAFDN_D  ,IDIMDERC ,IDIMDERH
     &          ,IDIMWORK,IFLAGS    ,INEW     ,INEWT    ,INTI
     &          ,IODIRECT,IPAR_DIR  ,ITERM    ,MAINF    ,MAXNB
     &          ,MAXNBF  ,NBAND1    ,NFLAGS   ,NPAR     ,NPARALG
     &          ,NUMNP   ,PAR_DIR   ,SOLUTION ,WORK)

********************************************************************************
*
* PURPOSE
*
*     Assembles RHS of the flow inverse problem and solve the system when
*     coupled density dependent flow and trasnport is solved.
*
********************************************************************************

      IMPLICIT NONE

C------------------------- External

      INTEGER*4::IA_COUPLED_DSC_COLS,IA_COUPLED_DSC_ROWS,IDIMDERC
     &          ,IDIMDERH ,IDIMWORK ,INEW     ,INEWT    ,INTI
     &          ,IODIRECT ,ITERM    ,MAINF    ,MAXNB    ,MAXNBF
     &          ,NBAND1   ,NFLAGS   ,NPAR     ,NPARALG  ,NUMNP

      
      INTEGER*4::IAD_D(2*MAXNB,2*NUMNP)  ,IADD_D(2*NUMNP)
     &          ,IADN_D(2*NUMNP)         ,IAFD_D(2*MAXNBF,NUMNP*2)
     &          ,IAFDD_D(NUMNP*2)        ,IAFDN_D(NUMNP*2)
     &          ,IFLAGS(NFLAGS)          ,IPAR_DIR(NPARALG)  


      REAL*8::A_COUPL_DSC(IA_COUPLED_DSC_ROWS,IA_COUPLED_DSC_COLS)
     &       ,A_COUPL_DSCF(MAXNB,2*NUMNP)  ,BCOUPLED(2*NUMNP)
     &       ,DERC(NUMNP,NPAR,IDIMDERC)    ,DERH(NUMNP,NPAR,IDIMDERH)
     &       ,PAR_DIR(NPARALG)             ,SOLUTION(2*NUMNP)
     &       ,WORK(IDIMWORK)

C------------------------- Internal

      INTEGER*4::I      ,IPAR   ,J

      CHARACTER::strFmt1*20

C------------------------- First executable statement

      strFmt1 = ''
      WRITE (strFmt1,*) NPAR
      strFmt1 = '(I5,'//Trim(AdjustL(strFmt1))//'F15.10)'

C-------------------- Writes DERC befor solving

      IF (IFLAGS(25).GT.0) THEN

          WRITE(747,*) 'DERH BEFORE INTI = ',INTI
          WRITE(746,*) 'DERC BEFORE INTI = ',INTI
          DO I=1,NUMNP
              WRITE(747,strFmt1) I,(DERH(I,J,INEW),J=1,NPAR)
              WRITE(746,strFmt1) I,(DERC(I,J,INEW),J=1,NPAR)
          END DO !I=1,NUMNP

      END IF !IFLAGS(25).GT.0

      DO IPAR=1,NPAR

          CALL ASSEMBLE_RHS_COUPLED
     &        (BCOUPLED,DERH(1,IPAR,INEW),DERC(1,IPAR,INEWT),NUMNP)


          CALL SOLVE
     &        (IA_COUPLED_DSC_ROWS,IA_COUPLED_DSC_COLS,2*NUMNP   
     &        ,1           ,IPAR_DIR(21)  ,IDIMWORK      ,IODIRECT    ,3
     &        ,INTI        ,IPAR_DIR(18)  ,IPAR_DIR(22)  ,0         
     &        ,ITERM       ,IPAR_DIR(15)  ,MAINF
     &        ,2*NBAND1    ,2*IPAR_DIR(24),IPAR_DIR(20)  ,IPAR_DIR(17)
     &        ,IPAR_DIR(16),PAR_DIR(36)   ,PAR_DIR(37)   ,PAR_DIR(38)
     &        ,A_COUPL_DSC ,A_COUPL_DSCF  ,BCOUPLED      ,IAD_D
     &        ,IADD_D      ,IADN_D        ,IAFD_D        ,IAFDD_D
     &        ,IAFDN_D     ,WORK          ,SOLUTION)


          CALL UNSHUFFLE_RHS
     &        (DERC(1,IPAR,INEWT),DERH(1,IPAR,INEW),NUMNP,SOLUTION)

      END DO !IPAR=1,NPAR

C-------------------- Writes DERC, DERH after solving

      IF (IFLAGS(25).GT.0) THEN

          WRITE(747,*) 'DERH AFTER INTI = ',INTI

          DO I=1,NUMNP
              WRITE(747,strFmt1) I,(DERH(I,J,INEW),J=1,NPAR)
          END DO !I=1,NUMNP

          WRITE(746,*) 'DERC AFTER INTI = ',INTI

          DO I=1,NUMNP
              WRITE(746,strFmt1) I,(DERC(I,J,INEWT),J=1,NPAR)
          END DO !I=1,NUMNP

c    1     FORMAT(I5,<NPAR>F15.10)

      END IF !IFLAGS(25).GT.0

      IF (IFLAGS(4).GT.0) THEN

          WRITE(749,*) 'DERH AFTER INTI = ',INTI
          WRITE(750,*) 'DERC AFTER INTI = ',INTI

          DO I=1,NUMNP

              WRITE(749,strFmt1) I,(DERH(I,J,INEW),J=1,NPAR)
              WRITE(750,strFmt1) I,(DERC(I,J,INEWT),J=1,NPAR)

          END DO !I=1,NUMNP

      END IF !IFLAGS(4).GT.0

      END SUBROUTINE JAC_COUPLED
