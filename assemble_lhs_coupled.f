      SUBROUTINE ASSEMBLE_LHS_COUPLED
     &          (A_COUPLED_DSC  ,AFLU     ,ATRA     ,CFLU     ,DBFLUDFLU
     &          ,DBFLUDTRA      ,DBTRADTRA,DFLU     ,DFLUDFLU 
     &          ,DFLUDTRA       ,DTRA     
     &          ,DTRADFLU 
     &          ,DTRADTRA       ,I_TYPE_A_DSC       ,I_TYPE_C            
     &          ,I_TYPE_D       ,IA_COUPLED_DSC_COLS,IA_COUPLED_DSC_ROWS
     &          ,IAD            ,IADD     ,IADN     ,IDIMAFLU ,IDIMATRA
     &          ,IDIMCFLU       ,IDIMDFLU ,IDIMDTRA ,INDSSTR
     &          ,ISPARSE        ,KXX      ,LMXNDL   ,LNNDEL   ,MAXNB    
     &          ,MAXNN          ,NBAND    ,NUMEL    ,NUMNP    ,TINC     
     &          ,THETAF         ,THETAT   ,EPSFLU   ,EPSTRA)

      IMPLICIT NONE

      INTEGER*4::I_TYPE_A_DSC       ,I_TYPE_C           ,I_TYPE_D
     &          ,IA_COUPLED_DSC_COLS,IA_COUPLED_DSC_ROWS,IDIMAFLU
     &          ,IDIMATRA ,IDIMCFLU ,IDIMDFLU ,IDIMDTRA ,IDSC_COLS
     &          ,IDSC_ROWS,INDSSTR  ,ISPARSE  ,LMXNDL   ,MAXNB
     &          ,MAXNN    ,NBAND    ,NUMEL    ,NUMNP

      INTEGER*4::IAD(MAXNB, MAXNN) ,IADD(MAXNN)   ,IADN(MAXNN)
     &          ,KXX(LMXNDL,NUMEL) ,LNNDEL(NUMEL)

      REAL*8::EPSFLU   ,EPSTRA   ,THETAF   ,THETAT   ,TINC     ,TINCINV

      REAL*8::A_COUPLED_DSC(IA_COUPLED_DSC_ROWS,IA_COUPLED_DSC_COLS)
     &       ,AFLU(NUMEL,IDIMAFLU)    ,ATRA(NUMEL, IDIMATRA)
     &       ,CFLU(NUMEL,IDIMCFLU)    ,DBFLUDFLU(NUMNP)
     &       ,DBFLUDTRA(NUMNP)        ,DBTRADTRA(NUMNP)
     &       ,DFLU(NUMEL,IDIMDFLU)    ,DFLUDFLU(NUMEL,LMXNDL)
     &       ,DFLUDTRA(NUMEL,LMXNDL)  ,DTRA(NUMEL,IDIMDTRA)
     &       ,DTRADFLU(NUMEL,LMXNDL*LMXNDL)
     &       ,DTRADTRA(NUMEL,LMXNDL*LMXNDL)

C------------------------- Initializing some constants multipliers.

      TINCINV = 1D0/(TINC)

C------------------------- Matrix type

      IF (ISPARSE.EQ.0) I_TYPE_A_DSC= 8
      IF (ISPARSE.EQ.1) I_TYPE_A_DSC= 9
      
C------------------------- Row and column number

      IDSC_COLS = IA_COUPLED_DSC_COLS 
      IDSC_ROWS = IA_COUPLED_DSC_ROWS


C---------------------------------------
C------------------------- Block 1 -----
C---------------------------------------

      A_COUPLED_DSC=0D0

C------------------------- Assemble Theta*AFLU

      CALL ASSEMBLE_SHUFFLED
     ;(THETAF     ,6             ,I_TYPE_A_DSC  ,IDIMAFLU,NUMEL  
     ;,IDSC_COLS  ,IDSC_ROWS     ,1             ,LMXNDL    ,NUMEL     
     ;,MAXNB      ,MAXNN         ,AFLU          ,A_COUPLED_DSC,IAD           
     ;,IADD       ,IADN          ,KXX           ,LNNDEL    ,NBAND)


C------------------------- Assemble   tincinv*DFLU

      IF (INDSSTR.EQ.1) THEN
           CALL ASSEMBLE_SHUFFLED 
     ;(TINCINV    ,I_TYPE_D      ,I_TYPE_A_DSC  ,IDIMDFLU  ,NUMEL  
     ;,IDSC_COLS  ,IDSC_ROWS     ,1             ,LMXNDL    ,NUMEL     
     ;,MAXNB      ,MAXNN         ,DFLU          ,A_COUPLED_DSC,IAD           
     ;,IADD       ,IADN          ,KXX           ,LNNDEL    ,NBAND)
	ENDIF
      

C------------------------- Assemble DFLUDFLU

      CALL ASSEMBLE_SHUFFLED
     ;(EPSFLU     ,4             ,I_TYPE_A_DSC  ,LMXNDL*LMXNDL,NUMEL 
     ;,IDSC_COLS  ,IDSC_ROWS     ,1             ,LMXNDL    ,NUMEL     
     ;,MAXNB      ,MAXNN         ,DFLUDFLU      ,A_COUPLED_DSC,IAD           
     ;,IADD       ,IADN          ,KXX           ,LNNDEL    ,NBAND)

C------------------------- Assemble DBFLUDFLU

      CALL ASSEMBLE_SHUFFLED
     ;(-EPSFLU    ,1             ,I_TYPE_A_DSC  ,1         ,NUMNP
     ;,IDSC_COLS  ,IDSC_ROWS     ,1             ,LMXNDL    ,NUMEL     
     ;,MAXNB      ,MAXNN         ,DBFLUDFLU     ,A_COUPLED_DSC,IAD           
     ;,IADD       ,IADN          ,KXX           ,LNNDEL    ,NBAND)


C---------------------------------------
C------------------------- Block 2 -----
C---------------------------------------

C------------------------- Assemble tincinv*CFLU

      IF (INDSSTR.EQ.1) THEN

	   CALL ASSEMBLE_SHUFFLED
     ;(TINCINV    ,I_TYPE_C      ,I_TYPE_A_DSC  ,IDIMCFLU  ,NUMEL    
     ;,IDSC_COLS  ,IDSC_ROWS     ,2             ,LMXNDL    ,NUMEL     
     ;,MAXNB      ,MAXNN         ,CFLU          ,A_COUPLED_DSC,IAD           
     ;,IADD       ,IADN          ,KXX           ,LNNDEL    ,NBAND)

      ENDIF

C------------------------- Assemble DFLUDTRA

      CALL ASSEMBLE_SHUFFLED
     ;(EPSTRA     ,4             ,I_TYPE_A_DSC  ,LMXNDL*LMXNDL,NUMEL  
     ;,IDSC_COLS  ,IDSC_ROWS     ,2             ,LMXNDL    ,NUMEL     
     ;,MAXNB      ,MAXNN         ,DFLUDTRA      ,A_COUPLED_DSC,IAD           
     ;,IADD       ,IADN          ,KXX           ,LNNDEL    ,NBAND)

C------------------------- Assemble DBFLUDTRA

      CALL ASSEMBLE_SHUFFLED
     ;(-EPSTRA    ,1             ,I_TYPE_A_DSC  ,1         ,NUMNP
     ;,IDSC_COLS  ,IDSC_ROWS     ,2             ,LMXNDL    ,NUMEL     
     ;,MAXNB      ,MAXNN         ,DBFLUDTRA     ,A_COUPLED_DSC,IAD           
     ;,IADD       ,IADN          ,KXX           ,LNNDEL    ,NBAND)



C---------------------------------------
C------------------------- Block 3 -----
C---------------------------------------

C------------------------- Assemble DTRADFLU

      CALL ASSEMBLE_SHUFFLED
     ;(EPSFLU     ,4             ,I_TYPE_A_DSC  ,IDIMATRA  ,NUMEL 
     ;,IDSC_COLS  ,IDSC_ROWS     ,3             ,LMXNDL    ,NUMEL     
     ;,MAXNB      ,MAXNN         ,DTRADFLU      ,A_COUPLED_DSC,IAD           
     ;,IADD       ,IADN          ,KXX           ,LNNDEL    ,NBAND)


C---------------------------------------
C------------------------- Block 4 -----
C---------------------------------------


C------------------------- Assemble theta*ATRA

      CALL ASSEMBLE_SHUFFLED
     ;(THETAT     ,4             ,I_TYPE_A_DSC  ,IDIMATRA  ,NUMEL  
     ;,IDSC_COLS  ,IDSC_ROWS     ,4             ,LMXNDL    ,NUMEL     
     ;,MAXNB      ,MAXNN         ,ATRA          ,A_COUPLED_DSC,IAD           
     ;,IADD       ,IADN          ,KXX           ,LNNDEL    ,NBAND)

C------------------------- Assemble   tincinv*DTRA

      IF (INDSSTR.EQ.1) THEN

         CALL ASSEMBLE_SHUFFLED
     ;(TINCINV    ,I_TYPE_D      ,I_TYPE_A_DSC  ,IDIMDTRA  ,NUMEL  
     ;,IDSC_COLS  ,IDSC_ROWS     ,4             ,LMXNDL    ,NUMEL     
     ;,MAXNB      ,MAXNN         ,DTRA          ,A_COUPLED_DSC,IAD           
     ;,IADD       ,IADN          ,KXX           ,LNNDEL    ,NBAND)

      ENDIF

   
C------------------------- Assemble DTRADTRA

      CALL ASSEMBLE_SHUFFLED
     ;(EPSTRA     ,4             ,I_TYPE_A_DSC  ,IDIMATRA  ,NUMEL  
     ;,IDSC_COLS  ,IDSC_ROWS     ,4             ,LMXNDL    ,NUMEL     
     ;,MAXNB      ,MAXNN         ,DTRADTRA      ,A_COUPLED_DSC,IAD           
     ;,IADD       ,IADN          ,KXX           ,LNNDEL    ,NBAND)


C------------------------- Assemble DBTRADTRA

      CALL ASSEMBLE_SHUFFLED
     ;(EPSTRA    ,1             ,I_TYPE_A_DSC  ,1         ,NUMNP
     ;,IDSC_COLS  ,IDSC_ROWS     ,4             ,LMXNDL    ,NUMEL     
     ;,MAXNB      ,MAXNN         ,DBTRADTRA     ,A_COUPLED_DSC,IAD           
     ;,IADD       ,IADN          ,KXX           ,LNNDEL    ,NBAND)

      END SUBROUTINE ASSEMBLE_LHS_COUPLED