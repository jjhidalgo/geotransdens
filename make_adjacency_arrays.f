	SUBROUTINE MAKE_ADJACENCY_ARRAYS
     ;(IODENS     ,IPRECOND   ,NUMEL      ,NUMNP    
     ;,LMXNDL     ,LEVEL      ,MAXNB      ,MAXNBF     ,NPARALG    
     ;,IPAR_DIR   ,IAD_D      ,IADD_D   ,IADN_D
     ;,IAD_S      ,IADD_S     ,IADN_S   ,IAFD_S   ,IAFDD_S
     ;,IAFDN_S    ,IAFD_D     ,IAFDD_D  ,IAFDN_D  ,KXX
     ;,LNNDEL)

*******************************************************************************
C
C PURPOSE
C This subroutine makes the calls needed to generate the arrays that watsolve
C needs to make before the first system can be solved.
C The calls differ when making a system for the flow or transport equation on
C one hand and the ones for making a system for coupled flow and transport 
C with newtons method
C
C
**********************************************************************************
      IMPLICIT NONE
C EXTERNAL VARIABLES: SCALARS
	INTEGER*4 IODENS, NUMEL,NUMNP, LMXNDL, LEVEL, NPARALG
     ;          ,MAXNB,MAXNBF,IPRECOND

C EXTERNAL VARIABLES: ARRAYS
      INTEGER*4 IPAR_DIR(NPARALG),KXX(LMXNDL,NUMEL)
     ;,IAD_D(MAXNB*2,2*NUMNP) ,IADD_D(2*NUMNP)    ,IADN_D(2*NUMNP)
     ;,IAD_S(MAXNB,NUMNP)     ,IADD_S(NUMNP)      ,IADN_S(NUMNP)   
     ;
     ;,IAFD_S(MAXNBF,NUMNP)   ,IAFDD_S(NUMNP)
     ;,IAFDN_S(NUMNP)
     ;,IAFD_D(2*MAXNBF,NUMNP*2),IAFDD_D(NUMNP*2)  
     ;,IAFDN_D(NUMNP*2),LNNDEL(NUMEL)

C INTERNAL VARIABLES: SCALARS
	INTEGER*4 ICAN_CN


	IF (IODENS.EQ.1) THEN  !if we have variable density
!	   IF (NRITCTNFL.NE.0 AND NRITCTNTR.NE.0) ICAN_CN=1
	   IF (IPAR_DIR(10).NE.0 .AND. IPAR_DIR(13).NE.0) ICAN_CN=1
	ENDIF


C_________if watsolv will not be used for coupled system
	

C_____________Create adjecancy arrays for single system 
	 CALL FE_IADM
     ;(KXX,NUMNP,NUMEL,LMXNDL,IADD_S,IAD_S,IADN_S,MAXNB,LNNDEL)
	
C_____________create adjecancy array for preconditioned matrix 
	 IF (IPRECOND.NE.0) CALL SYMFAC
     ;(NUMNP    ,LEVEL     ,MAXNB    ,IAD_S    ,IADD_S   ,IADN_S
     ;,IAFD_S   ,IAFDD_S   ,IAFDN_S    ,MAXNBF)

C_________if watsolv will be used for coupled system
	IF (ICAN_CN.EQ.1) THEN  

      CALL FE_IADM_D
     ;(KXX,NUMNP,NUMEL,LMXNDL,IADD_D,IAD_D,IADN_D,MAXNB,LNNDEL)


C_____________create adjecancy array for preconditioned matrix 
	  IF(IPRECOND.NE.0) CALL SYMFAC     
     ;(2*NUMNP  ,LEVEL     ,2*MAXNB    ,IAD_D    ,IADD_D   ,IADN_D
     ;,IAFD_D   ,IAFDD_D   ,IAFDN_D    ,2*MAXNBF)

	ENDIF

	RETURN
	END SUBROUTINE MAKE_ADJACENCY_ARRAYS

