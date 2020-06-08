      SUBROUTINE MINIMIZATION
     ;(ALF      ,FNEW     ,FOLD     ,GNORM     ,GNORM1    ,IDIMCOV
     ;,IOMIN    ,ISUMFO   ,MAINF    ,MIN_STOP  ,NBANDCOV  ,NDEVS
     ;,NFLAGS   ,NITERF1  ,NITERF2  ,NPAR      ,NPARALG   ,NSTAT
     ;,NUMITER  ,NUMTOBS  ,NWRITE   ,NZPAR     ,OBJCON    ,OBJHED    
     ;,OBJPAR   ,XMAXIM   ,COVINV   ,COVPAR    ,DLT_PAR   ,FOBJ_WGT  
     ;,GRAD     ,HESS     ,HESSAUX  ,IFLAGS    ,IOWRITE   ,IPAR_INV
     ;,PAR_INV  ,PARAUX   ,PARC      ,PARGOOD   ,PARM     ,PARZ
     ;,VJAC     ,VOBS     ,VOBSC     ,WORK      ,WGT_UNK  ,MEASTYP
     ;,PHI      ,OBJGEO)

C_______________________ Step 0: Declaration of variables

      IMPLICIT NONE
                                                              ! Integer external
      INTEGER*4 IOMIN,NPARALG,IDIMCOV,NWRITE,ISUMFO,MAINF,MIN_STOP
     ;         ,NBANDCOV,NFLAGS,NDEVS,NITERF1,NITERF2,NPAR,NSTAT,NUMITER
     ;         ,NUMTOBS,NZPAR
     ;         ,IOWRITE(NWRITE),IPAR_INV(NPARALG),IFLAGS(NFLAGS)
     ;         ,MEASTYP(NUMTOBS)
                                                                 ! Real external
      REAL*8::ALF,FNEW,FOLD,GNORM,GNORM1,OBJCON,OBJGEO,OBJHED,OBJPAR
     ;      ,XMAXIM,PHI
     ;      ,PAR_INV(NPARALG),COVINV(IDIMCOV),COVPAR(NPAR*(NPAR+1)/2)
     ;      ,DLT_PAR(NPAR),FOBJ_WGT(NSTAT),WORK(2*NPAR),WGT_UNK(NPAR)
     ;      ,GRAD(NPAR),HESS(NPAR*(NPAR+1)/2),HESSAUX(NPAR*(NPAR+1)/2)
     ;      ,PARAUX(NPAR),PARC(NPAR),PARM(NPAR),VOBS(NUMTOBS)
     ;      ,VOBSC(NUMTOBS+NDEVS),VJAC(NUMTOBS,NPAR),PARGOOD(NZPAR)
     ;      ,PARZ(NZPAR)
      
C_________________________________Only Marquardt's method is operative

      IF (IOMIN.EQ.1) CALL MARQUARDT
     ;(ALF        ,PAR_INV(7)  ,PAR_INV(6)  ,PAR_INV(10) ,FNEW
     ;,FOLD       ,PAR_INV(5)  ,PAR_INV(4)  ,GNORM       ,GNORM1
     ;,IDIMCOV    ,IOWRITE(15) ,IOWRITE(14) ,ISUMFO      ,MAINF       
     ;,IPAR_INV(3),IPAR_INV(4) ,MIN_STOP    ,NBANDCOV    ,NDEVS       
     ;,NFLAGS     ,NITERF1     ,NITERF2     ,IPAR_INV(5) ,NPAR        
     ;,NSTAT      ,IPAR_INV(2) ,IPAR_INV(1) ,NUMITER     ,NUMTOBS     
     ;,NZPAR      ,OBJCON      ,OBJHED      ,OBJPAR      ,PAR_INV(3)  
     ;,PAR_INV(2) ,PAR_INV(1)  ,XMAXIM      ,COVINV      ,COVPAR      
     ;,DLT_PAR    ,FOBJ_WGT    ,GRAD        ,HESS        ,HESSAUX     
     ;,IFLAGS     ,PARAUX      ,PARC        ,PARGOOD     ,PARM
     ;,PARZ        ,VJAC       ,VOBS        ,VOBSC       ,WGT_UNK
     ;,WORK        ,MEASTYP    ,PHI         ,OBJGEO)

      RETURN 
      END
