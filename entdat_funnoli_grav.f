       SUBROUTINE ENTDAT_FUNNOLI_GRAV
     & (IERROR     ,INPWR      ,IODENS_INI ,IOFLLI     ,IOTRLI
     & ,IOWAR      ,IUPAR      ,MAINF      ,NFNL       ,NROW
     & ,NZPRG      ,FILENAME   ,GRAV       ,NFNLPRG    ,NFNLTIP
     & ,PARACD)

*****************************************************************************
* PURPOSE
*     Reads generic parameter zone number, agreement parameters and/or 
*     gravity direction
*
* DESCRIPTION
*     Reads generic parameter zone number, agreement parameters and/or 
*     gravity direction
*
* EXTERNAL VARIABLES: ARRAYS
*
*  FILENAME               Array containing names for input and output           
*                         data files                                            
*  GRAV                   Gravity array direction 
*  NFNLPRG                Generic parameter zone number for every nonlinear
*                         function
*  NFNLTIP                Type of non-linear function                           
*  PARACD                 Agreement parameters                                  
*
* EXTERNAL VARIABLES: SCALARS
*
*  IERROR                 Current number of errors on input data                
*  INPWR                  Allows writing on MAIN FILE                           
*  IODIM                  Maximum dimension of any element included             
*                         in the problem                                        
*  IOFLLI                 If zero, linear flow problem, otherwise set to 1      
*  IOTRLI                 Idem to IOFLLI, in the case of transport.             
*  IOWAR                  Allows writing warning message in SUBROUTINE ERROR    
*  IUPAR                  Unit number of file PAR                               
*  MAINF                  Unit number of the main output file (RES.OUT)         
*  NFNL                   Total number of non-linear functions required         
*                         in the problem                                        
*  NROW                   Current record number                                 
*  NZPRG                  Total number of generic parameter zones               
*
* INTERNAL VARIABLES: SCALARS
*
*  GRAVNOR                Norm of gravity array
*  LEAUX                  Auxiliar string where the last read line of the       
*                         current input file is stored                          
*
* FUNCTIONS AND SUBROUTINES REFERENCED
*
*  ERROR                  Writes the current error message and error number     
*  LEEL                   Returns a string value containing the current line    
*                         of XXX FILE, if no coment appears.                    
*
* HISTORY
*
*     AMS        1988     First coding
*     SCR      4-1997     Revision and verification
*     AMS      1-1998     Common elimination and addition of header
*****************************************************************************
     
       IMPLICIT REAL*8(A-H,O-Z)
       CHARACTER*20 FILENAME(20),LEEL*100,LEAUX*100
       DIMENSION NFNLTIP(NFNL),NFNLPRG(8,NFNL),PARACD(3,NFNL),GRAV(3)

       IF (IOFLLI.NE.0 .OR. IOTRLI.NE.0) THEN

C------------------------- Writes title in MAIN output file

          IF (INPWR.NE.0) WRITE(MAINF,3000)
 3000     FORMAT(///,'       NON LINEAL FUNCTIONS CHARACTERISTICS',/,
     ;               '       --- ------ --------- ---------------',//,
     ;          '    # TYPE          GROUP PARAMTERS INDEX         ',
     ;          '  AGREEMENT PARAMETERS VALUES',/,
     ;          '    - ----  (1)  (2)  (3)  (4)  (5)  (6)  (7)  (8)',
     ;          '       (1)       (2)       (3)')
 
 3100     FORMAT(2I5,I4,7I5,3E10.4)
           

          DO NFUN=1,NFNL

C------------------------- Reads index and parameters
        
             LEAUX=LEEL(FILENAME,IUPAR,MAINF,NROW,INPWR)
             READ(LEAUX,1000,ERR=9000) NFUN1,NFNLTIP(NFUN1),
     ;           (NFNLPRG(I,NFUN1),I=1,8),(PARACD(I,NFUN1),I=1,3)
 1000        FORMAT(10I5,3F10.0)

C------------------------- Checks data

             IF (NFUN.NE.NFUN1)
     ;          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME, 
     ;               ' NON LINEAL FUNCTION ARE OUT OF SECUENCE',
     ;                 NROW,0,IUPAR,1,6.78)

             IF (NFNLTIP(NFUN1).LE.0)
     ;          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                ' NON LINEAL FUNCTION TYPE ARE OUT OF'//
     ;                ' ORDER',NROW,0,IUPAR,1,6.79)
         
             DO IPG=1,8

                IF  (NFNLPRG(IPG,NFUN1).LT.0 .OR.
     ;                                NFNLPRG(IPG,NFUN1).GT.NZPRG)
     ;              CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                  ' NON LINEAL GROUP PARAMETERS INDEX ARE OUT'
     ;                //' OF ORDER',0,0,IUPAR,1,6.80)
          
             ENDDO

C------------------------- Writes read data in MAIN file
         
             IF (INPWR.NE.0)
     ;          WRITE(MAINF,3100)NFUN1,NFNLTIP(NFUN1),
     ;             (NFNLPRG(I,NFUN1),I=1,8),(PARACD(J,NFUN1),J=1,3)

          ENDDO
       ENDIF    ! IOFLLI.NE.0 .OR. IOTRLI.NE.0


C------------------------- Reads gravity vector when we work with pressure

       IF (IOFLLI.EQ.1 .OR. IODENS_INI.EQ.1) THEN

          LEAUX=LEEL(FILENAME,IUPAR,MAINF,NROW,INPWR)
          READ(LEAUX,1100,ERR=9100) (GRAV(I),I=1,3)
 1100     FORMAT(3F10.0)

          IF (INPWR.NE.0) THEN
              WRITE(MAINF,3200)
 3200         FORMAT(/,'  GRAVITY ARRAY ',/,
     ;                '  ------------- ')

              WRITE(MAINF,3300) (GRAV(I),I=1,3)
 3300         FORMAT(3G10.3)

          END IF !INPWR.NE.0

          IF (IOWAR.NE.0) THEN

C------------------------- Computes gravity norm

             GRAVNOR=DSQRT( GRAV(1)*GRAV(1)+GRAV(2)*GRAV(2)+
     ;                      GRAV(3)*GRAV(3) )

C------------------------- Checks if gravity norm is one

             IF (DABS(GRAVNOR-1.D0).GE.1E-3) 
     ;          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                  'GRAVITY IS NOT A UNITARY ARRAY'
     ;                  ,NROW,0,IUPAR,0,0.00)

C------------------------- Checks if gravity norm is almost zero

             IF (DABS(GRAVNOR).LT.1E-2) 
     ;          CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;                  'GRAVITY NORM IS SMALLER THAN 0.01'
     ;                  ,NROW,0,IUPAR,0,0.00)

C------------------------- Tells if gravity direction follows one of the 
C------------------------- cartesian axes. The code assumes that gravity 
C------------------------- direction follows one axis if its component in 
C------------------------- this axis is at least 50 times larger than the 
C------------------------- rest of its components

C------------------------- Gravity follows X axis

             IF (DABS(GRAV(2)).LT.1E-3 .AND. DABS(GRAV(3)).LT.1E-3
     ;                                 .AND. DABS(GRAV(1)).GT.5E-2) 
     ;          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;                  'YOUR GRAVITY DIRECTION IS X AXIS'
     ;                  ,NROW,0,IUPAR,0,0.00)

C------------------------- Gravity follows Y axis

             IF (DABS(GRAV(1)).LT.1E-3 .AND. DABS(GRAV(3)).LT.1E-3
     ;                                 .AND. DABS(GRAV(2)).GT.5E-2) 
     ;          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;                  'YOUR GRAVITY DIRECTION IS Y AXIS'
     ;                  ,NROW,0,IUPAR,0,0.00)

C------------------------- Gravity follows Z axis

             IF (DABS(GRAV(1)).LT.1E-3 .AND. DABS(GRAV(2)).LT.1E-3
     ;                                 .AND. DABS(GRAV(3)).GT.5E-2) 
     ;          CALL ERROR (IERROR,IOWAR,MAINF,FILENAME,
     ;                  'YOUR GRAVITY DIRECTION IS Z AXIS'
     ;                  ,NROW,0,IUPAR,0,0.00)
           
          ENDIF
       ENDIF

       RETURN
      
9000   CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;           'GENERIC FORTRAN ERROR WHEN READING'
     ;         //'NON-LINEAR FUNCTIONS DATA' 
     ;           ,NROW,0,IUPAR,1,6.81)
       RETURN

9100   CALL ERROR(IERROR,IOWAR,MAINF,FILENAME,
     ;           'GENERIC FORTRAN ERROR WHEN READING'
     ;         //'GRAVITY FLOW PARAMETERS' 
     ;           ,NROW,0,IUPAR,1,6.82)
       RETURN
       END
