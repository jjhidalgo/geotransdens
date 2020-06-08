        SUBROUTINE INTERP_NZNUM 
     ;(   I          ,IN         ,IOLD       ,INALF      ,INALFT
     ;   ,INCHP      ,INCHPT     ,INQQP      ,INQQPT     ,INCON
     ;   ,INCONT     ,INDMT      ,IOEQT      ,IOTRS      ,IORTS
     ;   ,NTDMT      ,NUMNP      ,INPWR      ,NPARNP     ,MAINF
     ;   ,IXPARNP    ,IBTCO      ,IBCOD      ,NZALF      ,NZCOE
     ;   ,NZCHP      ,NZQQP      ,INCLK      ,NZCLK)

*******************************************************************************
* PURPOSE 
*
*       This subroutine interpolates zone numbers for missing nodes    
* 
* DESCRIPTION
*
*       This subroutines starts a zone number interpolaation loop between two
*       non consecutive nodes. Problem type to be solved, and transient or 
*       steady-state, are checked in order to define flow and transport 
*       parameters zone numbers and write them on MAIN file. 
*
* 
* EXTERNAL VARIABLES: ARRAYS 
*
*     IXPARNP(INpar,j)     Array containing zone number of each element j,
*                          corresponding to INpar index parameter zone.
*     IBTCO                Transport boundary condition node
*     IBCOD                Flow boundary condition node
*
*
* EXTERNAL VARIABLES: SCALARS
*
*     I                    Index corresponding to actual node
*     IN                   Index corresponding to last node read
*     INALF                Zone index for leakage in steady-state
*     INALFT               Idem for transient
*     INCHP                Zone index for prescribed head in steady-state
*     INCHPT               Idem for transient
*     INQQP                Zone index for prescribed flow in steady-state
*     INQQPT               Idem for transient state
*     INCON                Zone index for external concentration in steady-state
*     INCONT               Idem for transient
*     INDMT                Zone index for molecular diffusion.
*     IOLD                 Last node read befor IN node
*     IOEQT                Type of problem to be solved
*     IOTRS                Flow regim
*     IORTS                Transport regim
*     NZDMT                Number of matrix diffusion zones
*     NUMNP                Number of nodes
*     INPWR                Printing option
*     NPARNP               Total zone numbers of nodal parameters
*                     
* INTERNAL VARIABLES: SCALARS 
*
*    IAUXIBD               Auxiliary variable containing flow boundary 
*                          condition type to be written in MAIN file
*    IAUXIBT               The same for transport 
*    IAUXCHP               Auxiliary variable conaining prescribed head
*                          zone number or current node to be written in 
*                          in MAIN file
*    IAUXCHPT              The same for transient prescribed head
*    IAUXQQP               The same for prescribed flow
*    IAUXQQPT              The same for transient prescribed flow
*    IAUXALF               The same for leakage
*    IAUXALFT              The same for transient leakage
*    IAUXCON               The same for external concentration
*    IAUXCONT              The same for transient external concentration
*    IAUXDMT               The same for matrix diffusion
*
* HISTORY
*
*     SCR 04-abr-1997     First coding
*     AMS       10-00     Initializes to zero variables IAUX
*
*******************************************************************************

       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION IXPARNP(NUMNP,NPARNP),IBTCO(NUMNP),IBCOD(NUMNP)
  
C------------------------- Initializes to zero some auxiliar variables

       IAUXIBD=  0
       IAUXIBT=  0
       IAUXALF=  0
       IAUXCHP=  0
       IAUXQQP=  0
       IAUXCHPT= 0
       IAUXQQPT= 0
       IAUXALFT= 0
       IAUXCON=  0
       IAUXCONT= 0
       IAUXDMT=  0
       IAUXCLK=  0
     
       NN=IN-IOLD
*_______________________Interpolation loop
             
       DO 100 II=1,NN-1
         IF (IOEQT.NE.2) THEN                                  !Not only flow   

*_______________________Define flow zone number parameters

           IBCOD(I)=IBCOD(IOLD)
           IF (NZALF.NE.0)IXPARNP(I,INALF)=IXPARNP(IOLD,INALF)
           IF (NZCHP.NE.0)IXPARNP(I,INCHP)=IXPARNP(IOLD,INCHP)
           IF (NZQQP.NE.0)IXPARNP(I,INQQP)=IXPARNP(IOLD,INQQP)

*_______________________Define transient flow zone number parameters

            IF (NZCHP.NE.0)IXPARNP(I,INCHPT)=IXPARNP(IOLD,INCHPT)
            IF (NZQQP.NE.0)IXPARNP(I,INQQPT)=IXPARNP(IOLD,INQQPT)
            IF (NZALF.NE.0)IXPARNP(I,INALFT)=IXPARNP(IOLD,INALFT)

         ENDIF                                                  ! IOEQT.NE.2

         IF (IOEQT.NE.1) THEN   !Not only transport

*_______________________Define transport boundary condition type and 
*_______________________external concentration zone number

           IBTCO(I)=IBTCO(IOLD)
           IF (NZCOE.NE.0)IXPARNP(I,INCON)=IXPARNP(IOLD,INCON)
	     IF (NZCLK.NE.0)IXPARNP(I,INCLK)=IXPARNP(IOLD,INCLK)

*_______________________Checks transient and define transient external
*_______________________concentration and matrix diffusion zones numbers

           IF (IORTS.NE.0. AND .NZCOE.NE.0)
     ;       IXPARNP(I,INCONT)=IXPARNP(IOLD,INCONT)
           IF(NTDMT.NE.0)IXPARNP(I,INDMT)=IXPARNP(IOLD,INDMT)
         END IF                                                     ! IOEQT.NE.1

         IF (INPWR.NE.0) THEN

*_______________________Checks type of problem to be solved before writting on 
*_______________________MAIN file. Data from last interpolated node

           IF (IOEQT.NE.2) THEN                      !Not only transport problem

*_______________________Define flow parameters zone number and boundary cond
*_______________________to be written in MAIN file

             IAUXIBD=IBCOD(I)
             IF (NZALF.NE.0)IAUXALF=IXPARNP(I,INALF)
             IF (NZCHP.NE.0)IAUXCHP=IXPARNP(I,INCHP)
             IF (NZQQP.NE.0)IAUXQQP=IXPARNP(I,INQQP)
             IF (IOTRS.NE.0) THEN                               !Check transient

*______________________Define transient parameters zone number
*______________________to be written in MAIN file

 
               IF (NZCHP.NE.0)IAUXCHPT=IXPARNP(I,INCHPT)
               IF (NZQQP.NE.0)IAUXQQPT=IXPARNP(I,INQQPT)
               IF (NZALF.NE.0)IAUXALFT=IXPARNP(I,INALFT)
             END IF                                                 ! IOTRS.NE.0
           END IF                                                   ! IOEQT.NE.2

           IF (IOEQT.NE.1) THEN                           !Not only flow problem

*_______________________Define transport boundary condition type and 
*_______________________external concentration zone nummber
*_______________________to be written in MAIN file

             IAUXIBT=IBTCO(I)
             IF (NZCOE.NE.0)IAUXCON=IXPARNP(I,INCON)
             IF (NZCLK.NE.0)IAUXCLK=IXPARNP(I,INCLK)

*_______________________Define external concentration and matrix diffusion
*_______________________zone number to be written in MAIN file

             IF (IORTS.NE.0) IAUXCONT=IXPARNP(I,INCONT)   !Check transient
             IF(NTDMT.NE.0)IAUXDMT=IXPARNP(I,INDMT)
           END IF                                                   ! IOEQT.NE.1

           WRITE (MAINF,3000) I,IAUXIBD,IAUXIBT,IAUXCHP,IAUXCHPT,
     ;       IAUXQQP,IAUXQQPT,IAUXALF,IAUXALFT,IAUXCON,IAUXCONT,IAUXDMT
     &      ,IAUXCLK
 3000      FORMAT(5X,13I5)
         END IF                                                     ! INPWR.NE.0
         I=I+1
  100    CONTINUE                                             !Next missing node

         RETURN
         END
