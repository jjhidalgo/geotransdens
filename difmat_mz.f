      SUBROUTINE UPD_REC_MZ
     .( NUD_MZ, EXP_MZ,   REC0,   RECN,   CCAL,CCALAN,  NUMNP,IXDMT
     ;                       ,NUM_ZON)

************************************************************************
* PURPOSE
*
*        Updates vectors REC0 (Io) and RECN (In)
*
* DESCRIPTION
*
*        The equations for updating are given in TRANSIN-IV user's guide
*
*        Io[k+1,i] = Io[k,i] + Go*(C[k+1,i]-C[k,i])
*
*        where k = time increment counter
*              i = node number
*              Go= steady-state memory function
*              C = vector of nodal concentrations
*
*        In[k+1,i] = In[k,i]*E1n[1,n] + E2n*(C[k+1,i]-C[k,i])
*
*        where n = term of series expansion counter
*              E1n = EXP( (ALFn**2 + l_D)*DELT_D)
*              E2n = [An*ALFn**2 / (ALFn**2 + l_D)**2]*(1-E1n)/DELT_D
*              DELT_D = dimensionless time increment
*
*        these equations are applied for all terms and nodes of the current 
*        matrix diffusion zone
*
*   ARGUMENTS
*
*        NUD_MZ (NUM_NP_MZ)        List of nodes belonging to the current zone
*        EXP_MZ (2,NUM_NP_MZ)      E1n, E2n
*        REC0 (NUM_NP_MZ)          On output, Io[k+1,i]; on input, Io[k,i]
*        RECN (NUM_TER,NUM_NP_MZ)  On output, In[k+1,i]; on input, In[k,i]
*        CCAL (NUMNP)              C[k+1,i]
*        CCALAN (NUMNP)           C[k,i]
*        NUMNP                     Number of nodal points
*
*
* JCR April, 1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER NUD_MZ, NUMNP , IZ, IG , N,IXDMT,NUM_ZON
      REAL*8  REC0 , RECN , CCAL , CCALAN , EXP_MZ, DELC

      DIMENSION NUD_MZ(NUM_NP_MZ) , REC0(NUM_NP_MZ) ,
     .          RECN(NUM_TER,NUM_NP_MZ) , CCAL(NUMNP) , 
     .          CCALAN(NUMNP),
     .          EXP_MZ (2,NUM_TER),IXDMT(NUMNP)


* ............................................. Loop over all nodes of the zone

      DO IZ = 1, NUM_NP_MZ

         IG = NUD_MZ(IZ)                   ! IG is the global number of node IZ
        if (ixdmt(ig).eq.num_zon) then
         DELC = CCAL(IG) - CCALAN(IG)

         IF (I_LM_MZ.NE.0) THEN                       ! Io needs to be computed
            REC0(IZ) = REC0(IZ) + G0_MZ * DELC
         END IF
                                           ! Loop over all terms(n) to update In
         DO N = 1, NUM_TER
             RECN(N,IZ) = RECN(N,IZ)*EXP_MZ(1,N) + EXP_MZ(2,N)* DELC
         END DO
        endif
      END DO ! IZ

      RETURN
      END


      SUBROUTINE UPD_DREC_MZ
     .( NUD_MZ, EXP_MZ,D_EXP_MZ,  RECN,  DREC0,  DRECN,   CCAL,
     . CCALAN,  NUMNP,   DERC,DERCAN,   NPAR,IXDMT,NUM_ZON)

************************************************************************
* PURPOSE
*
*        Updates derivatives of vectors REC0 (Io) and RECN (In)
*
* DESCRIPTION
*
*        The equations for updating are given in TRANSIN-IV user's guide.
*        These equations are applied for all terms and nodes of the current 
*        matrix diffusion zone
*
*   ARGUMENTS
*
*        NUD_MZ (NUM_NP_MZ)        List of nodes belonging to the current zone
*        EXP_MZ (2,NUM_NP_MZ)      E1n, E2n
*        D_EXP_MZ (2,NUM_NP_MZ,3)  Derivatives of E1n, E2n w.r.t. MZ param
*        REC0 (NUM_NP_MZ)          On output, Io[k+1,i]; on input, Io[k,i]
*        RECN (NUM_TER,NUM_NP_MZ)  On output, In[k+1,i]; on input, In[k,i]
*        DREC0 (NUM_NP_MZ,NPAR)    Derivatives of Io with respect to parameters
*        DRECN (NUM_TER,NUM_NP_MZ,NPAR)    Idem of In
*        CCAL (NUMNP)              C[k+1,i]
*        CCALAN (NUMNP)           C[k,i]
*        DERC (NUMNP,NPAR)         Derivatives of C[k+1,i] with respect to param
*        DERCAN (NUMNP,NPAR)      Derivatives of C[k,i] with respect to param
*        NUMNP                     Number of nodal points
*        NPAR                      Number of model parameters
*
* JCR April, 1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER NUD_MZ , NUMNP , IZ, IG, NPAR , NP, N, ND,IXDMT,NUM_ZON

      REAL*8  RECN, DREC0, DRECN, CCAL, CCALAN, EXP_MZ, D_EXP_MZ,
     .        DELC , DDELC , DERC , DERCAN

      DIMENSION NUD_MZ(NUM_NP_MZ) ,
     .          RECN(NUM_TER,NUM_NP_MZ), DREC0(NUM_NP_MZ,NPAR) ,
     .          DRECN(NUM_TER,NUM_NP_MZ,NPAR) , CCAL(NUMNP) , 
     .          CCALAN(NUMNP), EXP_MZ(2,NUM_TER), 
     .          D_EXP_MZ(2,NUM_TER,3),
     .          DERC(NUMNP,NPAR) , 
     .          DERCAN(NUMNP,NPAR),IXDMT(NUMNP)

* .................. First, compute the generic contributions of all parameters

      DO NP = 1, NPAR                                ! Loop over all parameters

         DO IZ = 1, NUM_NP_MZ         ! Loop over all nodes of the current zone

            IG = NUD_MZ(IZ)
        if (ixdmt(ig).eq.num_zon) then
            DDELC = DERC(IG,NP) - DERCAN(IG,NP)

            IF (I_LM_MZ.NE.0) THEN                             ! d_Io[k+1,i]/d_p
               DREC0(IZ,NP) = DREC0(IZ,NP) + G0_MZ * DDELC
            END IF
                                                         

            DO N = 1, NUM_TER
                                                               ! d_In[k+1,i]/d_p
                DRECN(N,IZ,NP) = DRECN(N,IZ,NP)*EXP_MZ(1,N) + 
     .                           EXP_MZ(2,N)*DDELC
            END DO
        endif
         END DO ! IZ

      END DO ! NP

* ........................... Now, compute the contributions of MZ parameters

      IF (MZ_EST.EQ.0) RETURN           ! No MZ parameter is estimated

      DO ND = 1,3

         NP=MZ_PAR_NUM(ND)

         IF (NP.NE.0) THEN ! This parameter is estimated

            DO IZ = 1, NUM_NP_MZ
               IG = NUD_MZ(IZ)
        if (ixdmt(ig).eq.num_zon) then
               DELC = CCAL(IG) - CCALAN(IG)

               IF (I_LM_MZ.NE.0) THEN                          ! d_Io[k+1,i]/d_p
                  DREC0(IZ,NP) = DREC0(IZ,NP) + D_G0_MZ(ND) * DELC
               END IF

               DO N = 1, NUM_TER
                                                               ! d_In[k+1,i]/d_p
                  DRECN(N,IZ,NP) = DRECN(N,IZ,NP) +
     .                             D_EXP_MZ(1,N,ND)*RECN(N,IZ) + 
     .                             D_EXP_MZ(2,N,ND)*DELC

               END DO  ! N
        endif
            END DO ! IZ

         END IF ! (NP.EQ.0)

      END DO ! ND
            
      RETURN
      END

      SUBROUTINE COMP_DER_MZ
     .( NUD_MZ, VOL_MZ, EXP_MZ,D_EXP_MZ,RHS_MZ,   RECN,
     .   DREC0,  DRECN,   CCAL,CCALAN,   DERC, NUMNP,   NPAR
     ;              ,IXDMT,NUM_ZON)

************************************************************************
* PURPOSE
*
*        Computes contributions of MZ to DERC (RHS of sensitivity equations)
*
* DESCRIPTION
*
*        The equations for updating are given in TRANSIN-IV user's guide.
*        These equations are applied for all terms and nodes of the current 
*        matrix diffusion zone
*
*   ARGUMENTS
*
*        NUD_MZ (NUM_NP_MZ)        List of nodes belonging to the current zone
*        VOL_MZ (NUM_NP_MZ)        Volume of nodes belonging to the current zone
*        AN_MZ  (2,NUM_NP_MZ)      An, ALFn**2
*        EXP_MZ (2,NUM_NP_MZ)      E1n, E2n
*        D_EXP_MZ (2,NUM_NP_MZ,3)  Derivatives of E1n, E2n w.r.t. MZ param
*        RHS_MZ (NUM_NP_MZ)        K'm of nodes belonging to the current zone  
*        REC0 (NUM_NP_MZ)          On output, Io[k+1,i]; on input, Io[k,i]
*        RECN (NUM_TER,NUM_NP_MZ)  On output, In[k+1,i]; on input, In[k,i]
*        DREC0 (NUM_NP_MZ,NPAR)    Derivatives of Io with respect to parameters
*        DRECN (NUM_TER,NUM_NP_MZ,NPAR)    Idem of In
*        CCAL (NUMNP)              C[k+1,i]
*        CCALAN (NUMNP)           C[k,i]
*        DERC (NUMNP,NPAR)         Derivatives of C[k+1,i] with respect to param
*        DERCAN (NUMNP,NPAR)      Derivatives of C[k,i] with respect to param
*        NUMNP                     Number of nodal points
*        NPAR                      Number of model parameters
*
* JCR April, 1999
* AMS derivative temrs corrected May, 2004
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER NUD_MZ, NUMNP , IZ, IG, NPAR, NP, N, ND,IXDMT,NUM_ZON

      REAL*8  VOL_MZ, EXP_MZ, D_EXP_MZ, RHS_MZ,
     .        RECN, DREC0, DRECN,
     .        CCAL, CCALAN, DERC,
     .        DELTA_T, DELC, DLCDT, PDL, PD, PR, SUM

      DIMENSION NUD_MZ(NUM_NP_MZ),
     .          EXP_MZ(2,NUM_TER) , D_EXP_MZ(2,NUM_TER,3) ,
     .          RHS_MZ(NUM_NP_MZ) ,IXDMT(NUMNP), 
     .          VOL_MZ(NUM_NP_MZ) ,
     .          RECN(NUM_TER,NUM_NP_MZ), DREC0(NUM_NP_MZ,NPAR) ,
     .          DRECN(NUM_TER,NUM_NP_MZ,NPAR) , CCAL(NUMNP) , 
     .          CCALAN(NUMNP), DERC(NUMNP,NPAR)

* .......................................................... Auxiliary variables

      DELTA_T = DELT_D*CRD_MZ*ESP_MZ*ESP_MZ/DFM_MZ
      PDL = POR_MZ*DFM_MZ/ESP_MZ/ESP_MZ               ! POR*Dm/(Lm**2)
      PD = POR_MZ*DFM_MZ                              ! POR*Dm
      PR = POR_MZ*CRD_MZ                              ! POR*Rm

* ................... First, compute the generic contributions of all parameters

      DO NP = 1, NPAR

         DO IZ = 1, NUM_NP_MZ

            IG = NUD_MZ(IZ)
                                                                  ! d_K'm[i]/d_p
        if (ixdmt(ig).eq.num_zon) then
            IF (I_LM_MZ.EQ.0) THEN
               SUM=0D0
            ELSE
               SUM = DREC0(IZ,NP)
            END IF

            DO N = 1, NUM_TER
                                                               ! d_In[k+1,i]/d_p
               SUM = SUM - DRECN(N,IZ,NP)*EXP_MZ(1,N)

            END DO

            DERC(IG,NP) = DERC(IG,NP) + VOL_MZ(IZ)*PDL*SUM 
        endif
         END DO  ! IZ

      END DO ! NP


* ........................... Now, compute the contributions of MZ parameters

      IF (MZ_EST.EQ.0) RETURN                  ! No MZ parameter is estimated

      IF (ND_UPD_DKM.EQ.0) THEN        ! Derivatives of Km need to be updated

         CALL UPD_DKM_MZ (D_EXP_MZ)

      END IF


      DO ND = 1,4                                   ! Loop over MZ parameters

         NP=MZ_PAR_NUM(ND)                   ! Global number of NDth MZ parameter

         IF (NP.NE.0) THEN                      ! This parameter is estimated

            DO IZ = 1, NUM_NP_MZ                         ! Loop over MZ nodes

               IG = NUD_MZ(IZ)
        if (ixdmt(ig).eq.num_zon) then
               DELC = CCAL(IG) - CCALAN(IG)
               DLCDT= DELC/DELTA_T

                            ! Contibution of coefficients (por and Dm) of K'm
               IF (ND.EQ.1) THEN  
                  DERC(IG,NP) = DERC(IG,NP) + RHS_MZ(IZ)/DFM_MZ     ! Dm
               ELSE IF (ND.EQ.4) THEN
                  DERC(IG,NP) = DERC(IG,NP) + RHS_MZ(IZ)/POR_MZ    ! por
               END IF

                                      ! Add contribution of derivatives of Km
               DERC(IG,NP) = DERC(IG,NP) - 
     .                               D_KM_MZ(ND)*VOL_MZ(IZ)*DLCDT

C------------------------ Derivatives w.r.t Dm, Rm and lambdam

               IF (ND.LT.4) THEN
                  SUM=0D0                                ! Sum of In*d_E1n/dp
                  DO N = 1, NUM_TER
                      SUM = SUM + RECN(N,IZ)*D_EXP_MZ(1,N,ND)
                  END DO 
                  DERC(IG,NP) = DERC(IG,NP) - VOL_MZ(IZ)*PDL*SUM 

               END IF !(ND.LT.4) 
        endif
            END DO ! IZ

          END IF ! (NP.EQ.0)

      END DO ! ND
 
      RETURN
      END


      SUBROUTINE UPD_DKM_MZ (D_EXP_MZ)

************************************************************************
* PURPOSE
*
*        Computes derivative of Km with respect to matrix diffusion parameters
*
* DESCRIPTION
*
*        The equations for updating are given in TRANSIN-IV user's guide.
*
*   ARGUMENTS
*
*        EXP_MZ (2,NUM_NP_MZ)      E1n, E2n
*        D_EXP_MZ (2,NUM_NP_MZ,3)  Derivatives of E1n, E2n with respect to MZ param
*
* JCR April, 1999
* AMS derivatives corrected May, 2004
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER N, ND

      REAL*8  D_EXP_MZ, DT_P(3), SUM

      DIMENSION D_EXP_MZ(2,NUM_TER,3)

* ......................................... First, compute derivatives of DELT_D
      DT_P(1) = DELT_D/DFM_MZ
      DT_P(2) = - DELT_D/CRD_MZ
      DT_P(3) = 0D0

* ................ Now, compute derivatives of K'm with respect to MZ parameters

      DO ND = 1,4

         D_KM_MZ(ND) = 0D0                                    ! by default, zero

         IF (MZ_PAR_NUM(ND).NE.0) THEN             ! This parameter is estimated

            IF (ND.EQ.4) THEN
               D_KM_MZ(ND) = VAL_KM_MZ/POR_MZ                      ! d(Km)/d(por)

            ELSE                                          ! Remaining parameters

               IF (I_LM_MZ.EQ.0) THEN 
                  SUM = 0D0     ! No contribution from Go (S.S. memory function)
               ELSE
                  SUM = - D_G0_MZ(ND)                         ! Derivative of Go
               END IF

               DO N=1,NUM_TER
                  SUM = SUM + D_EXP_MZ(2,N,ND)           ! Add derivative of E2n
               END DO

               D_KM_MZ(ND) = POR_MZ*CRD_MZ*DELT_D*SUM                ! that's it

                                      ! Add contribution of derivative of DELT_D
               D_KM_MZ(ND) = D_KM_MZ(ND) + VAL_KM_MZ*DT_P(ND)/DELT_D 

               IF (ND.EQ.2) THEN                        ! Add contribution of Rm
                  D_KM_MZ(ND) = D_KM_MZ(ND) + VAL_KM_MZ/CRD_MZ 
               END IF

            END IF  ! (ND.EQ.4)

         END IF   ! (MZ_PAR_NUM(ND).NE.0) 

      END DO ! ND

      ND_UPD_DKM=1                                    ! D_KM_MZ has been updated

      RETURN
      END

      SUBROUTINE UPD_KM_MZ (EXP_MZ)

************************************************************************
* PURPOSE
*
*        Computes Km 
*
* DESCRIPTION
*
*        The equations are given in TRANSIN-IV user's guide (Eq.A.40)
*
*   ARGUMENTS
*
*        EXP_MZ (2,NUM_NP_MZ)      E1n, E2n
*
* JCR April, 1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER N

      REAL*8  EXP_MZ, SUM

      DIMENSION EXP_MZ(2,NUM_TER)


* ............................................ ................Computation of Km 

      IF (I_LM_MZ.EQ.0) THEN 
         SUM = 0D0                                     ! No contribution from Go
      ELSE
         SUM = - G0_MZ               ! Contribution of Go (S.S. memory function)
      END IF

* ................................................................... Sum of E2n

      DO N=1,NUM_TER
         SUM = SUM + EXP_MZ(2,N)                              ! Add value of E2n
      END DO

      VAL_KM_MZ = POR_MZ*CRD_MZ*DELT_D*SUM                ! Km=por*Dm*Delt_D*Sum

      ND_UPD_KM=1                                   ! VAL_KM_MZ has been updated
      ND_UPD_DKM=0                                     ! But its derivatives not

      RETURN
      END

      SUBROUTINE COMP_CNTR_MZ
     .(  NUD_MZ, EXP_MZ,  RHS_MZ, VOL_MZ, REC0,   RECN,
     .   BTRA, DTRA, NUMNP,IXDMT,NUM_ZON)

************************************************************************
* PURPOSE
*
*        Adds contribution of the current matrix diffusion zone(MZ) to
*        BTRA and DTRA
*
* DESCRIPTION
*
*        The equations are given in TRANSIN-IV user's guide
*
*        DTRA[i] = DTRA[i] + V[i]*Km
*        BTRA[i] = BTRA[i] + K'm + V[i]*Km*C[k,i]
*
*        where k = time increment counter
*              i = node number
*              C = vector of nodal concentrations
*
*        these equations are applied for all nodes of the current 
*        matrix diffusion zone
*
*   ARGUMENTS
*
*        NUD_MZ     List of nodes belonging to the current matrix diffusion zone
*        NUM_NP_MZ  Number of nodes belonging to the current matrix diffusion zone
*        REC0         On output, Io[k+1,i]; on input, Io[k,i]
*        RECN         On output, In[k+1,i]; on input, In[k,i]
*        CCAL         C[k+1,i]
*        NUMNP        Number of nodal points
*
* JCR April, 1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER NUD_MZ, NUMNP , IZ, IG , N,IXDMT,NUM_ZON

      REAL*8  EXP_MZ , RHS_MZ, VOL_MZ, REC0 , RECN , 
     .        BTRA, DTRA,
     .        PDL, PR


      DIMENSION NUD_MZ(NUM_NP_MZ),IXDMT(NUMNP), 
     .          EXP_MZ(2,NUM_TER) , RHS_MZ(NUM_NP_MZ) , 
     .          VOL_MZ(NUM_NP_MZ) , REC0(NUM_NP_MZ) ,
     .          RECN(NUM_TER,NUM_NP_MZ), BTRA(NUMNP), DTRA(NUMNP)

* .......................................................... Auxiliary variables

      PDL = POR_MZ*DFM_MZ/ESP_MZ/ESP_MZ               ! POR*Dm/(Lm**2)
      PR = POR_MZ*CRD_MZ                              ! POR*Rm

* ............................................... If needed, updates Km and DTRA

      IF (ND_UPD_KM.EQ.0) THEN                            ! Indeed, it is needed

         DO IZ = 1, NUM_NP_MZ                    ! First, substract Km from DTRA
            IG = NUD_MZ(IZ)
        if (ixdmt(ig).eq.num_zon) then
            DTRA(IG) = DTRA(IG) - VOL_MZ(IZ) * VAL_KM_MZ
        endif
         END DO

         CALL UPD_KM_MZ (EXP_MZ)                                ! Now, update Km

         DO IZ = 1, NUM_NP_MZ                ! And, finally, add Km back to DTRA
            IG = NUD_MZ(IZ)
        if (ixdmt(ig).eq.num_zon) then
            DTRA(IG) = DTRA(IG) + VOL_MZ(IZ) * VAL_KM_MZ 
        endif
         END DO

      END IF ! (ND_UPD_KM.EQ.0)

* ............................................. Computes K'm and adds it to BTRA

      DO IZ = 1, NUM_NP_MZ                       ! for every node of the MZ zone
         IG = NUD_MZ(IZ)

        if (ixdmt(ig).eq.num_zon) then

         IF (I_LM_MZ.EQ.0) THEN
            RHS_MZ(IZ) = 0D0                           ! Initializes K'm to zero
         ELSE

            RHS_MZ(IZ) = REC0(IZ)                 ! Initializes K'm to Io[k+1,i]
         END IF

         DO N = 1, NUM_TER

             RHS_MZ(IZ) = RHS_MZ(IZ) -               ! Adds In[k,i]*exp(alfn...)
     .                    RECN(N,IZ)*EXP_MZ(1,N) 

         END DO

                                  ! Finally, multiplies K'm times POR*Dm/(Lm**2)
         RHS_MZ(IZ) = PDL * VOL_MZ(IZ) * RHS_MZ(IZ)

         BTRA(IG) = BTRA(IG) + RHS_MZ(IZ)                      ! add K'm to BTRA
        endif
      END DO ! IZ

      RETURN
      END


      SUBROUTINE UPD_DT_MZ (  AN_MZ, EXP_MZ, D_EXP_MZ, DT_NEW)

************************************************************************
* PURPOSE
*
*        Updates variables as a consequence of a change in Delta_t
*
* DESCRIPTION
*
*        The equations for updating are given in TRANSIN-IV user's guide.
*        These equations are applied for all terms and nodes of the current 
*        matrix diffusion zone
*
*   ARGUMENTS
*
*        AN_MZ (2,NUM_TER)     AN_MZ(1,n)=ALFn*ALFn;  AN_MZ(2,n)=An
*        EXP_MZ(2,NUM_TER)     E1n = EXP( (ALFn**2 + l_D)*DELT_D)
*                              E2n = [An*ALFn**2 / (ALFn**2 + l_D)**2]*(1-E1n)
*        D_EXP_MZ(2,NUM_TER,3) Derivatives of EXP w.r.t. MZ parameters
*        DT_NEW                New time increment
*
* JCR April, 1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER N
      REAL*8  AN_MZ, EXP_MZ, D_EXP_MZ, DL_P(3), DT_P(3), U, V, 
     .        D_E1N, D_E2N, DT_NEW

      DIMENSION AN_MZ(2,NUM_TER),EXP_MZ(2,NUM_TER),
     .          D_EXP_MZ(2,NUM_TER,3)

* .............................................. Start by upating DELT_D itself

      DELT_D = DT_NEW * DFM_MZ /CRD_MZ /ESP_MZ/ESP_MZ
      
* ............................................ Derivatives of DELTA_D and XLAM_D

      IF (MZ_EST.NE.0) THEN                  ! Derivatives of DELTA_D are needed

         DT_P(1) = DELT_D/DFM_MZ                                     ! w.r.t. Dm
         DT_P(2) = - DELT_D/CRD_MZ                                   ! w.r.t. Rm
         DT_P(3) = 0D0                                           ! w.r.t. XLA_MZ

         IF (I_LM_MZ.NE.0) THEN                                 ! Also of XLAM_D

            DL_P(1) = - XLAM_D/DFM_MZ                                ! w.r.t. Dm
            DL_P(2) = XLAM_D/CRD_MZ                                  ! w.r.t. Rm
            DL_P(3) = XLAM_D/XLA_MZ                              ! w.r.t. XLA_MZ

         END IF

      END IF

* ............................................ Update EXP_MZ and its derivatives

      DO N = 1, NUM_TER

         U = XLAM_D + AN_MZ(1,N)                                ! l_D+(alf_n)**2
         V = AN_MZ(1,N)*AN_MZ(2,N)/U/U                            ! An*ALFn/U**2

         EXP_MZ(1,N) = DEXP (-U*DELT_D)                                ! E1n
         EXP_MZ(2,N) = V*(1d0-EXP_MZ(1,N))/DELT_D                      ! E2n

         IF (MZ_EST.NE.0) THEN        ! ....... Derivatives of EXP_MZ are needed

            D_E1N = - U * EXP_MZ(1,N)            ! First, compute d_E1n/d_DELT_D
                              ! Now, multiply it times the derivatives of DELT_D
            D_EXP_MZ(1,N,1) = D_E1N*DT_P(1)                          ! w.r.t. Dm
            D_EXP_MZ(1,N,2) = D_E1N*DT_P(2)                          ! w.r.t. Rm
            D_EXP_MZ(1,N,3) = D_E1N*DT_P(3)                      ! w.r.t. XLA_MZ

                                                 ! First, compute d_E2n/d_DELT_D
            D_E2N = ( -D_E1N*V/DELT_D) - EXP_MZ(2,N)/DELT_D
                              ! Now, multiply it times the derivatives of DELT_D
            D_EXP_MZ(2,N,1) = D_E2N*DT_P(1)                          ! w.r.t. Dm
            D_EXP_MZ(2,N,2) = D_E2N*DT_P(2)                          ! w.r.t. Rm
            D_EXP_MZ(2,N,3) = D_E2N*DT_P(3)                      ! w.r.t. XLA_MZ

            IF (I_LM_MZ.NE.0) THEN   ! .......................... Also of XLAM_D

               D_E1N = - DELT_D * EXP_MZ(1,N)    ! First, compute d_E1n/d_XLAM_D

                        ! Now, add it multiplied times the derivatives of XLAM_D

               D_EXP_MZ(1,N,1) = D_EXP_MZ(1,N,1) + D_E1N*DL_P(1)     ! w.r.t. Dm
               D_EXP_MZ(1,N,2) = D_EXP_MZ(1,N,2) + D_E1N*DL_P(2)     ! w.r.t. Rm
               D_EXP_MZ(1,N,3) = D_EXP_MZ(1,N,3) + D_E1N*DL_P(3) ! w.r.t. XLA_MZ

                                                 ! First, compute d_E2n/d_XLAM_D
               D_E2N = V*EXP_MZ(1,N) - 2D0*EXP_MZ(2,N)/U

                        ! Now, add it multiplied times the derivatives of XLAM_D

               D_EXP_MZ(2,N,1) = D_EXP_MZ(2,N,1) + D_E2N*DL_P(1)     ! w.r.t. Dm
               D_EXP_MZ(2,N,2) = D_EXP_MZ(2,N,2) + D_E2N*DL_P(2)     ! w.r.t. Rm
               D_EXP_MZ(2,N,3) = D_EXP_MZ(2,N,3) + D_E2N*DL_P(3) ! w.r.t. XLA_MZ

            END IF  ! (I_LM_MZ.NE.0)

         END IF  ! (MZ_EST.NE.0)

      END DO !N

* ................ WARNING; It will have to update Km and its derivatives

      ND_UPD_KM=0
      ND_UPD_DKM=0

      RETURN
      END

      SUBROUTINE UPD_PAR_MZ 
     .(  AN_MZ, EXP_MZ, D_EXP_MZ, DELT, PARC,  NZPAR)

************************************************************************
* PURPOSE
*    Updates variables as a consequence of a change in model parameters
*
* DESCRIPTION
*
*    Start by indeed updating parameter values
*    Then, updates EXP_MZ and its derivatives
*    Finally, it updates XLAM_D, the steady state memory function Go 
*      and its derivatives
*
*   ARGUMENTS
*
*        AN_MZ (2,NUM_TER)        AN_MZ(1,n)=ALFn;  AN_MZ(2,n)=An
*        EXP_MZ(2,NUM_TER)        AUX[1,n] = EXP( (ALFn**2 + l_D)*DELT_D)
*                                 AUX[2,n] = An*ALFn**2 / (ALFn**2 + l_D)**2
*        D_EXP_MZ(2,NUM_TER,3)    Derivatives of EXP w.r.t. MZ parameters
*        DELT                     Time increment
*
* JCR April, 1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER NZPAR

      REAL*8  AN_MZ, EXP_MZ, D_EXP_MZ,
     .        DELT, PARC

      DIMENSION AN_MZ(2,NUM_TER), EXP_MZ(2,NUM_TER),
     .          D_EXP_MZ(2,NUM_TER,3), PARC(NZPAR)

* .................................... Start by indeed updating parameter values

      CALL UPD_PAR_VAL_MZ (   PARC,  NZPAR)

* ................................... XLAM_D and steady state memory function Go

      IF (I_LM_MZ.NE.0) THEN            

         XLAM_D = XLA_MZ * CRD_MZ * ESP_MZ * ESP_MZ / DFM_MZ

         CALL COMP_G0_MZ (AN_MZ)

      END IF

* ....................................... Now, update EXP_MZ and its derivatives

      CALL UPD_DT_MZ (  AN_MZ, EXP_MZ, D_EXP_MZ, DELT)

* ....................... WARNING; It will have to update Km and its derivatives

      ND_UPD_KM=0
      ND_UPD_DKM=0

      RETURN
      END

      SUBROUTINE COMP_G0_MZ (  AN_MZ)

************************************************************************
* PURPOSE
*    Computes steady state memory function Go and its derivatives
*
* DESCRIPTION
*    The equations depend on block geometry
*
* ARGUMENTS
*
*    AN_MZ (2,NUM_TER)   AN_MZ(1,n)=ALFn*ALFn;  AN_MZ(2,n)=An
*
* JCR April, 1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER N

      REAL*8  AN_MZ, U, DER_G0

      DIMENSION AN_MZ(2,NUM_TER)


      IF (IOP_AN.EQ.1) THEN
* .................................................... Slab shaped matrix blocks

           U=DSQRT(XLAM_D)
           G0_MZ = -U*DTANH(U)
           IF (MZ_EST.NE.0) 
     .        DER_G0= G0_MZ/(2D0*XLAM_D) - 1D0/(2D0*(DCOSH(U)**2)) 

      ELSE
* ......................................... Any other geometry for matrix blocks

           G0_MZ = 0D0
           IF (MZ_EST.NE.0) DER_G0 = 0D0

           DO N=1,NUM_TER
              U = AN_MZ(1,N) + XLAM_D
              G0_MZ = G0_MZ - XLAM_D * AN_MZ(2,N) / U
              IF (MZ_EST.NE.0) 
     .          DER_G0= DER_G0 - AN_MZ(2,N)*AN_MZ(1,N)/U/U
           END DO

      END IF
* ....................... Computes derivatives ofsteady state memory function Go

      IF (MZ_EST.NE.0) THEN                             ! Derivatives are needed

                                 ! Multiply DER_G0, (deriv. of Go w.r.t. XLAM_D)
                                               ! times the derivatives of XLAM_D
         D_G0_MZ(1) = DER_G0*(- XLAM_D/DFM_MZ)                       ! w.r.t. Dm
         D_G0_MZ(2) = DER_G0*(XLAM_D/CRD_MZ)                         ! w.r.t. Rm
         D_G0_MZ(3) = DER_G0*(XLAM_D/XLA_MZ)                     ! w.r.t. XLA_MZ

      END IF

      RETURN
      END

      SUBROUTINE COMP_VOL_MZ
     .(  NUD_MZ, VOL_MZ, VOLNOD, NUMNP)

************************************************************************
* PURPOSE
*
*        Identifies the nodes belonging to the current matrix diffusion zone
*
* DESCRIPTION
*
*        For all nodes,
*            first, check if the matrix zone number of the node (IXDMT)
*                   is equal to the current matrix diffusion zone (MTDZ_NUM)
*            second, if so, add the node to the list of nodes of the zone
*
*   ARGUMENTS
*
*        NUD_MZ(NUM_NP_MZ) List of nodes belonging to MTDZ_NUM
*        VOL_MZ(NUM_NP_MZ) Volume associated to nodes belonging to MTDZ_NUM
*        VOLNOD(NUMNP)     Volume associated to all nodes
*        NUMNP             Number of nodal points
*
* JCR April, 1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER NUD_MZ, NUMNP , IZ, IG

      REAL*8 VOL_MZ , VOLNOD

      DIMENSION NUD_MZ(NUM_NP_MZ) , VOL_MZ(NUM_NP_MZ) , VOLNOD(NUMNP)

* ........................................................ Loop over all nodes

      DO IZ = 1, NUM_NP_MZ

         IG=NUD_MZ(IZ)
         VOL_MZ(IZ) = VOLNOD(IG)

      END DO ! IZ

* ................ Define NUM_NP_MZ  Number of nodes belonging to zone MTDZ_NUM

      RETURN
      END

      SUBROUTINE ID_NP_MZ
     .(  NUD_MZ, IXDMT, NUMNP ,NPBMX   ,NPARNP)

************************************************************************
* PURPOSE
*
*        Identifies the nodes belonging to the current matrix diffusion zone
*
* DESCRIPTION
*
*        For all nodes,
*            first, check if the matrix zone number of the node (IXDMT)
*                   is equal to the current matrix diffusion zone (MTDZ_NUM)
*            second, if so, add the node to the list of nodes of the zone
*
*   ARGUMENTS
*
*        NUD_MZ(NUM_NP_MZ) List of nodes belonging to MTDZ_NUM
*        IXDMT(NUMNP)      IXDMT(i) is the matrix diffusion zone of node i
*        NUMNP             Number of nodal points
*
* JCR April, 1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER NUD_MZ, IXDMT, NUMNP , IZ, IG   ,NPARNP
     ;                        ,IPROB,IFL,IOLD,NPBMX

      DIMENSION NUD_MZ(NUM_NP_MZ) , IXDMT(NUMNP,NPARNP,NPBMX) 

* ........................................................ Loop over all nodes
      IZ=0

C----------------------- First loop over problems

      do iprob=1,NPBMX
      DO IG = 1, NUMNP

* ................... First, check if the matrix zone number of the node (IXDMT)
*                       is equal to the current matrix diffusion zone (MTDZ_NUM)

         IF (IXDMT(IG,9,iprob).EQ.MTDZ_NUM) THEN 
                                        !add IG to the list of nodes of the zone

            IZ=IZ+1

C-------------------------- It should add IG to the list of nodes only if it has not
C-------------------------- been accounted yet

            IFL=0
            IF (IZ.GT.1) THEN
               DO IOLD=1,IZ-1
                  IF (NUD_MZ(IOLD).EQ.IG) THEN
                     IFL=1                     ! It has been already accounted
                     GOTO 9000
                  ENDIF
               ENDDO
            ENDIF

 9000       IF (IFL.EQ.0) THEN           ! Not accounted yet
               NUD_MZ(IZ) = IG
            ELSE                         ! Already accounted
               IZ=IZ-1
            ENDIF

         END IF

      END DO ! IG

      enddo  ! iprob

* ................ Define NUM_NP_MZ  Number of nodes belonging to zone MTDZ_NUM

      NUM_NP_MZ = IZ

      RETURN
      END

      SUBROUTINE INI_PAR_MZ
     .(  NDPAR,  MAINF,  IVPAR,  PARC,INORPAR,NZONE_PAR, NTYPAR, NZPAR,
     .   IOINV, IERR)

************************************************************************
* PURPOSE
*
*        Initializes parameter values (PARC), global numbering (MZ_PAR_NUM)...
*
* DESCRIPTION
*
*        1st : Check if zone number is within allowable bounds
*
*        2nd : Assign parameter value from PARC to PAR_MZ
*              [PAR equals DFM (if NDPAR=1),CRD (if NDPAR=2),XLA (if NDPAR=3)
*               and POR (if NDPAR=4)]
*
*        3rd : Assign global zone number from IVPAR to MZ_PAR_NUM
*
*   ARGUMENTS
*
*        NDPAR               NDPAR=1 for DFM, 2 for CRD, 3 for XLA and 4 for POR
*        MAINF               Output unit number
*        IVPAR (NZPAR)       Global number of all parameters
*        PARC (NZPAR)        Current value of all parameters
*        INORPAR (NTYPAR)    Pointers to each type of parameters
*        NZONE_PAR(NTYPAR)    Number of zones of each type of parameters
*        NTYPAR              Number of types of parameters
*        NZPAR               Total number of parameters
*        NPAR                Number of parameters to be estimated
*        IERR                Error counter
*
*
* JCR April, 1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER NDPAR , MAINF , IVPAR , INORPAR , NZONE_PAR ,
     .        NTYPAR , NZPAR , IERR , NTYP(8) , IZN , IZMN, IP,
     .        NZ   ,IOINV
      REAL*8 PARC 
      CHARACTER *10 NOM_PAR(4)

      DIMENSION IVPAR(NZPAR) , PARC(NZPAR) , INORPAR(NTYPAR) , 
     .          NZONE_PAR(NTYPAR)

      DATA NTYP/14,17,16,15,9,12,11,10/, ! Misterious numbers of parameter types
                               ! only known to Agustin Medina and Andres Alcolea
     .     NOM_PAR/'DIFF.COEFF','RETRDATION','DECAY COEF',' POROSITY '/

* ................................... Find pointers to vectors of all parameters

      IP = INORPAR (NTYP(NDPAR))
      NZ = NZONE_PAR (NTYP(NDPAR+4))

* .............................. Check if zone number is within allowable bounds

      IZN = MZ_ZON_NUM(NDPAR)
      IZMN=1                        ! min. zone number for dif. mol. or porosity
      IF ( NDPAR.EQ.2 .OR. NDPAR.EQ.3 ) IZMN = 0 ! min zone num. for retrd,lamda

      IF ( IZN.LT.IZMN .OR. IZN.GT.NZ ) THEN              ! IZN is out of bounds

         WRITE(MAINF,1001) NOM_PAR(NDPAR), IZN
 1001    FORMAT(' ERROR: The ',A10,' zone number is',I5)
         WRITE(MAINF,1002) IZMN,NZ
 1002    FORMAT(8x,'It must be less than',I3,' and greater than',I5)
         IERR=IERR+1

      ELSE

* ................................... Assign parameter value from PARC to PAR_MZ
* .............................. and global zone number from IVPAR to MZ_PAR_NUM

         IF (NDPAR.EQ.1) THEN
                                                    ! Diffusion coefficient, Dm
            DFM_MZ = PARC(IP+IZN)             
            MZ_PAR_NUM(NDPAR) = 0 
            IF (IOINV.GT.0) THEN
               MZ_PAR_NUM(NDPAR) = IVPAR(IP+IZN)
               IF (IVPAR(IP+IZN).GT.0) MZ_EST=MZ_EST+1
            END IF
            MZ_ZON_NUM(NDPAR)=IP+IZN

         ELSE IF (NDPAR.EQ.2) THEN
                                                               ! Retardation, Rm
            IF (IZN.EQ.0) THEN                                   ! by default, 1
               CRD_MZ = 1D0
               MZ_PAR_NUM(NDPAR) = 0
               MZ_ZON_NUM(NDPAR) = 0
            ELSE
               CRD_MZ = PARC (IP+IZN)
               MZ_PAR_NUM(NDPAR) = 0 
               IF (IOINV.GT.0) THEN
                  MZ_PAR_NUM(NDPAR) = IVPAR(IP+IZN)
                  IF (IVPAR(IP+IZN).GT.0) MZ_EST=MZ_EST+1
               END IF
               MZ_ZON_NUM(NDPAR)=IP+IZN
            END IF

         ELSE IF (NDPAR.EQ.3) THEN
                                                        ! Decay coefficient, xla
            IF (NZ.EQ.0 .OR. IZN.EQ.0) THEN                      ! by default, 0
               XLA_MZ = 0D0
               MZ_PAR_NUM(NDPAR) = 0
               I_LM_MZ=0
               MZ_ZON_NUM(NDPAR) = 0
            ELSE
               XLA_MZ = PARC (IP+IZN)
               MZ_PAR_NUM(NDPAR) = 0 
               IF (IOINV.GT.0) THEN
                  MZ_PAR_NUM(NDPAR) = IVPAR(IP+IZN)
                  IF (IVPAR(IP+IZN).GT.0) MZ_EST=MZ_EST+1
               END IF
               I_LM_MZ=1
               MZ_ZON_NUM(NDPAR)=IP+IZN
               XLAM_D = XLA_MZ * CRD_MZ * ESP_MZ * ESP_MZ / DFM_MZ
            END IF

         ELSE IF (NDPAR.EQ.4) THEN
                                                                 ! Porosity, por
            POR_MZ = PARC (IP+IZN)
            MZ_PAR_NUM(NDPAR) = 0 
            IF (IOINV.GT.0) THEN
               MZ_PAR_NUM(NDPAR) = IVPAR(IP+IZN)
               IF (IVPAR(IP+IZN).GT.0) MZ_EST=MZ_EST+1
            END IF
            MZ_ZON_NUM(NDPAR)=IP+IZN

         END IF ! (NDPAR.EQ...)

      END IF ! ( IZN.LT.IZMN .OR. IZN.GT.NZ ) 

      RETURN
      END

      SUBROUTINE UPD_PAR_VAL_MZ
     .(   PARC,  NZPAR)

************************************************************************
* PURPOSE
*
*        Updates parameter values  [typically after an inverse prob. iteration]
*
* DESCRIPTION
*
*        Assign parameter value from PARC to PAR_MZ
*              [PAR equals DFM (if NDPAR=1),CRD (if NDPAR=2),XLA (if NDPAR=3)
*               and POR (if NDPAR=4)]
*        Updating is only needed for parameters that are variable
*
*   ARGUMENTS
*
*        PARC (NZPAR)        Current value of all parameters
*        NZPAR               Total number of parameters
*
*
* JCR April, 1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER NDPAR , NZPAR , IZN

      REAL*8 PARC

      DIMENSION PARC(NZPAR) 

* ........................................................ What am I doing here?

      IF (MZ_EST.EQ.0) THEN                        ! No parameter needs updating

         PRINT 1001
 1001    FORMAT(' WARNING: subroutine UPD_PAR_VAL_MZ called when',
     .          'No parameter needs updating')
         RETURN

      END IF


* ................................... Assign parameter value from PARC to PAR_MZ

      DO NDPAR=1,4

         IZN = MZ_ZON_NUM(NDPAR)

         IF (IZN.GT.0 .AND. MZ_PAR_NUM(NDPAR).GT.0) THEN    ! PAR needs updating

            IF (NDPAR.EQ.1) THEN
                                                    ! Diffusion coefficient, Dm
               DFM_MZ = PARC(IZN)

            ELSE IF (NDPAR.EQ.2) THEN
                                                               ! Retardation, Rm
               CRD_MZ = PARC (IZN)

            ELSE IF (NDPAR.EQ.3) THEN
                                                        ! Decay coefficient, xla
               XLA_MZ = PARC (IZN)

            ELSE IF (NDPAR.EQ.4) THEN
                                                                 ! Porosity, por
               POR_MZ = PARC (IZN)

            END IF ! (NDPAR.EQ...)

         END IF ! ( IZN.GT.0 ) 

      END DO ! NDPAR
      RETURN
      END


      SUBROUTINE COMP_AN_MZ
     .(  AN_MZ, MAINF, IUNIN, IERR)

************************************************************************
* PURPOSE
*    Computes  AN_MZ(1,n)=ALFn*ALFn;  AN_MZ(2,n)=An
*
* DESCRIPTION
*    The equations depend on block geometry
*    Mass balance is kept by ensuring that SUM(An/ALFn/ALFn)=1
*         If needed, the last An is modified
*
* ARGUMENTS
*
*    AN_MZ (2,NUM_TER)   AN_MZ(1,n)=ALFn*ALFn;  AN_MZ(2,n)=An
*    MAINF               Output unit number
*    IUNIN               Input unit number
*    IERR                Error counter
*
*
* JCR April, 1999
************************************************************************

      IMPLICIT NONE

      INCLUDE 'COMMON_MTDZ.FOR'

      INTEGER N, MAINF, IUNIN, IERR

      REAL*8  AN_MZ, PI, SUM

      DIMENSION AN_MZ(2,NUM_TER)


      IF (IOP_AN.EQ.1) THEN
* .................................................... Slab shaped matrix blocks

         PI = 2D0*DASIN(1D0)
         DO N = 1,NUM_TER
            AN_MZ(1,N) = ((2D0*N-1D0)*PI/2D0)*((2D0*N-1D0)*PI/2D0)
            AN_MZ(2,N) = 2D0
         END DO

      ELSE IF (IOP_AN.EQ.2) THEN
* .................................................. Matrix blocks are cylinders

         WRITE(MAINF,1001) IOP_AN
 1001    FORMAT(' ERROR: The option IOP_AN=',I2,' is not operative')
         IERR=IERR+1

      ELSE IF (IOP_AN.EQ.3) THEN
* .................................................... Matrix blocks are spheres

         PI = 2D0*DASIN(1D0)
         DO N = 1,NUM_TER
            AN_MZ(1,N) = (N*PI/3D0)*(N*PI/3D0)
            AN_MZ(2,N) = 2D0/3D0
         END DO
         
      ELSE IF (IOP_AN.EQ.4) THEN
* ...................................................... Matrix blocks are veins

         WRITE(MAINF,1001) IOP_AN
         IERR=IERR+1

      ELSE IF (IOP_AN.EQ.5) THEN
* ............................................. Matrix blocks are hollow spheres

         WRITE(MAINF,1001) IOP_AN
         IERR=IERR+1

      ELSE IF (IOP_AN.EQ.6) THEN
* ............................................. Matrix blocks are hollow spheres

         WRITE(MAINF,1001) IOP_AN
         IERR=IERR+1

      ELSE
* ..................................... Matrix blocks have an arbitrary geometry

         WRITE(MAINF,100)
  100    FORMAT(40X,'        n      ALFn*ALFn           An')
         DO N=1, NUM_TER
            READ (IUNIN,'(2F10.0)') AN_MZ(1,N),AN_MZ(2,N)
            WRITE(MAINF,101) N,AN_MZ(1,N),AN_MZ(2,N)
  101       FORMAT(40X,I10,2E15.5)
         END DO

      END IF

      IF (IERR.GT.0) RETURN

* .................................. Define last An so as to ensure mass balance

      SUM = 0D0
      DO N = 1,NUM_TER-1
         SUM = SUM + AN_MZ(2,N)/AN_MZ(1,N)
      END DO
      AN_MZ(2,NUM_TER) = (1D0 - SUM)*AN_MZ(1,NUM_TER)

      RETURN
      END

