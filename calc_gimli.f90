      SUBROUTINE CALC_GIMLI                                              &
     &           (CCAL     ,DCDP     ,DRESDP   ,IFLAGS   ,NDEVGEO        &
     &           ,NFLAGS   ,NPAR     ,NUMEL    ,NUMNP    ,POROSITY       &
     &           ,RESCAL   ,TABSOLUT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     PURPOSE
!     
!     Calls GIMLI to compute resisitivity values and their derivative
!     with respect to nodal concentration at a set of given points.
!     
!     
!     EXTERNAL VARIABLES: ARRAYS
!     
!     CCAL                   Computed concentration at every node.  
!     POROSITY               Porosity at every element.
!     DCDP                   Derivative of the concentration at every node
!                            with respect to the parameters.
!     DRESDP                 Derivative of resisitivity at the observation points with
!                            respect to the parameters.
!     IFLAGS                 Array with different writing options. Used mainly for 
!                            debugging.                                             
!     RESCAL                 Computed resistivities at every node   
!
!     TABSOLUT               Current absolut computation time. Needed to let GIMLI know the ERT 
!                            acquisition scheme. 
!     
!     EXTERNAL VARIABLES: SCALARS
!     
!     NDEVGEO                Number of observation points.
!     NPAR                   Number of parameters.
!     NFLAGS                 Used to dimension IFLAGS.
!     NUMNP                  Number of nodes. 
!     NUMEL                  Number of elements                                      
!     
!
!     INTERNAL VARIABLES: ARRAYS
!
!     NDARRAY_POR           Fortran array with model porosities converted to a Python Numpy Array          
!     NDARRAY_CCAL          Fortran array with model concentrations converted to a Python Numpy Array  
!     NDARRAY_DCDP          Fortran array with DCDP converted to a Python Numpy Array  
!
!     RETVAL                Returned values from Python function
!     
!     RETVAL                Results from RETVAL are stored in this Python Object   (Two values: a Tuple)
!     OBJ_RESCAL            Computed resistivities taken from RETTUPLE and stored in a Generic Object.
!     OBJ_DRESDP            Derivatives taken from RETTUPLE and stored in a Generic Object.
!     NDARRAY_RESCAL        Computed resistivities object converted to a Numpy Array  
!     NDARRAY_DRESDP        DRESDC object converted converted to a Numpy Array
!     PTR_RESCAL            Fortran pointer to computed resistivities inside NDARRAY_RESCAL 
!     PTR_DRESDC            Fortran pointer to derivatives inside NDARRAY_DRESDC
!     
!     IERROR                Integer for error control 
!     GIMLI_ERT             Python module where the function to be called is stored. 
!     ARGS                  Tuple to store the function arguments 
!     PATHS                 Path to the files  
!     
!     HISTORY
!     
!     JHG, APP       6-2020     First coding.
!     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE forpy_mod

      IMPLICIT NONE

!-------------------------External variables

      INTEGER*4::NDEVGEO, NPAR, NFLAGS, NUMNP, NUMEL

      REAL*8:: TABSOLUT
      INTEGER*4::IFLAGS(NFLAGS)

      REAL*8::CCAL(NUMNP), POROSITY(NUMEL), DCDP(NUMNP,NPAR)
      REAL*8::DRESDP(NDEVGEO,NPAR), RESCAL(NDEVGEO)

     
!-------------------------Internal variables
        
      TYPE(NDARRAY) :: NDARRAY_POR, NDARRAY_CCAL, NDARRAY_DCDP, NDARRAY_RESCAL, NDARRAY_DRESDP
    
      REAL*8, POINTER, DIMENSION(:) :: PTR_RESCAL   
      REAL*8, POINTER, DIMENSION(:,:) :: PTR_DRESDP 
	
      TYPE(OBJECT) :: RETVAL, OBJ_RESCAL, OBJ_DRESDCP              
      TYPE(TUPLE) :: ARGS, RETTUPLE                    
               
      INTEGER*4                   :: IERROR
      TYPE(MODULE_PY)             :: GIMLI_ERT
      TYPE(LIST)                  :: PATHS

!-------------------------First executable statment.

      IF (IFLAGS(3).EQ.1) CALL IO_SUB('CALC_GIMLI',0)

!-------------------------CALL TO GIMLI GO HERE

      !INITIALIZE FORPY
      !IERROR = FORPY_INITIALIZE()
        
        
      ! INSTEAD OF SETTING THE ENVIRONMENT VARIABLE PYTHONPATH,
      ! WE CAN ADD THE CURRENT DIRECTORY "." TO SYS.PATH
      IERROR = GET_SYS_PATH(PATHS)
      IERROR = PATHS%APPEND(".")
  
      ! FROM FORTRAN TO PYTHON
      ! CREATE WRAPPERS FOR THE FORTRAN-ARRAYS
      ! WE WRAP THE INPUT FROM TRANSDENS INTO ARRAYS THAT CAN BE USED BY PYTHON. 
      IERROR = NDARRAY_CREATE_NOCOPY(NDARRAY_CCAL, CCAL)
      CALL ERR_PRINT
      IERROR = NDARRAY_CREATE_NOCOPY(NDARRAY_POR, POROSITY)
      CALL ERR_PRINT
      IERROR = NDARRAY_CREATE_NOCOPY(NDARRAY_DCDP, DCDP)
      CALL ERR_PRINT

      ! IMPORT THE PYTHON FUNCTION. THE MAIN ONE THAT FORPY WILL CALL. 
      IERROR = IMPORT_PY(GIMLI_ERT, "ERT_OBS")
      CALL ERR_PRINT
		
      ! CREATE ARGUMENT FOR THE MAIN PYTHON FUNCTION.
      !THE FUNCTION WILL TAKE 'NDARRAY_CCAL', 'NDARRAY_POR', 'NDARRAY_DCDP'
      IERROR = TUPLE_CREATE(ARGS, 4)
      CALL ERR_PRINT
      IERROR = ARGS%SETITEM(0, NDARRAY_CCAL)
      CALL ERR_PRINT
      IERROR = ARGS%SETITEM(1, NDARRAY_POR)
      CALL ERR_PRINT
      IERROR = ARGS%SETITEM(2, NDARRAY_DCDP)
      CALL ERR_PRINT
      IERROR = ARGS%SETITEM(3, TABSOLUT)

      ! CALL THE FUNCTION "RUN_ERT" INSIDE THE "GIMLI_MOD.PY" FILE, WITH THE ARGUMENTS.
      IERROR = CALL_PY(RETVAL, GIMLI_ERT, "RUN_ERT", ARGS)
      CALL ERR_PRINT
        
      ! TAKE THE 2 OUTPUTS OF THE FUNCTION, "RETVAL", AND PUT THEM IN A TUPLE. 
      IERROR = CAST(RETTUPLE, RETVAL)
      CALL ERR_PRINT
        
      ! TAKE THE TWO OUTPUTS FROM THE TUPLE AND PUT THEM INSIDE GENERIC OBJECTS. 
      IERROR = RETTUPLE%GETITEM(OBJ_RESCAL, 0)
      CALL ERR_PRINT
      IERROR = RETTUPLE%GETITEM(OBJ_DRESDCP, 1)
      CALL ERR_PRINT
        
      ! CAST THE GENERIC OBJECTS TO NDARRAYS
      IERROR = CAST(NDARRAY_RESCAL, OBJ_RESCAL)
      CALL ERR_PRINT
      IERROR = CAST(NDARRAY_DRESDP, OBJ_DRESDCP)
      CALL ERR_PRINT
        
      ! USE "GET_DATA" TO ACCESS THE DATA IN THE ARRAYS AND PUT THEM IN FORTRAN ARRAYS "PTR"
      IERROR = NDARRAY_RESCAL%GET_DATA(PTR_RESCAL)
      CALL ERR_PRINT
      IERROR = NDARRAY_DRESDP%GET_DATA(PTR_DRESDP, ORDER='F')
      CALL ERR_PRINT

      DRESDP = PTR_DRESDP
      RESCAL = PTR_RESCAL


      CALL ARGS%DESTROY
      CALL GIMLI_ERT%DESTROY
      CALL RETVAL%DESTROY
      Call PATHS%DESTROY
      !CALL FORPY_FINALIZE

      IF(IFLAGS(3).EQ.1) CALL IO_SUB('COMP_OBS_GEO',1)

      RETURN
      END SUBROUTINE CALC_GIMLI
