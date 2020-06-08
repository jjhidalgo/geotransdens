      REAL*8  ESP_MZ,            ! Matrix thickness 
     .        DFM_MZ,            ! Dm, Diffusion coefficient of the current zone
     .        CRD_MZ,            ! Rm, Retardation coeff. of the current zone
     .        XLA_MZ,            ! L, First order reaction (decay) coefficient
     .        POR_MZ,            ! Por, Porosity of the current zone
     .        DELT_D,            ! Dimensionless time increment
     .        XLAM_D,            ! Dimensionless lambda
     .        G0_MZ,             ! Go, Steady state memory function
     .        D_G0_MZ,           ! d_Go/d_pD, Derivative of Go w.r.t. MZ param.
     .        VAL_KM_MZ,         ! Km
     .        D_KM_MZ            ! d_Km/d_pD, Derivative of Km w.r.t. MZ param.

      INTEGER MTDZ_NUM,         ! Current matrix diffusion zone number
     .        NUM_TER,          ! Number of terms in the series expansion
     .        NUM_NP_MZ,        ! Number of nodes in the matix diffusion zone
     .        I_LM_MZ,          ! IOLAM for MZ, decay is not active if I_LM_MZ=0
     .        IOP_AN,           ! MAtrix geometry option
     .        MZ_EST,           ! Number of MZ parameters to be estimated
     .        MZ_PAR_NUM,       ! MZ unknown parameters global number 
     .        MZ_ZON_NUM,       ! MZ parameter zone number (1,por;2,Dm;3,Rm;4,lam
     .        ND_UPD_KM,        ! If 0, Km needs to be updated
     .        ND_UPD_DKM,       ! If 0, d_Km/d_p needs to be updated
     ;        ICH_IT_MAR        ! If 1 we have just started a new inverse problem iteration
                                ! 0 otherwise

      COMMON /MTDZ/ ESP_MZ, DFM_MZ, CRD_MZ, XLA_MZ, POR_MZ, DELT_D,
     .              XLAM_D, G0_MZ, D_G0_MZ(3), VAL_KM_MZ, D_KM_MZ(8),
     .              MTDZ_NUM, NUM_TER, NUM_NP_MZ, I_LM_MZ, IOP_AN,
     .              MZ_EST, MZ_PAR_NUM(4), MZ_ZON_NUM(4), ND_UPD_KM,
     .              ND_UPD_DKM,ICH_IT_MAR