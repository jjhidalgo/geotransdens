      SUBROUTINE DEFINENAMES (IOPRHED,TYPENAME)

********************************************************************************
*
* PURPOSE
*
* The names of the types of measurements and parameters are stored 
* in the character array TYPENAME. If you come up with a new type of parameter
* or  measurement, fill in the name here. Do not use more than 10 characters
*
*
* EXTERNAL VARIABLES: ARRAYS
*
*  TYPENAME               Array containing  the names of the state var. and
*                         parameter types in the same order as OBSCLASS and STAT
*
* EXTERNAL VARIABLES: SCALARS
*
*  IOPRHED                Indicates whether the flow state variable state is    
*                         preasure (set to 1) or head (set to 0)                
*
* HISTORY: LJS (Dec. 2002): First coding
*          AAR (Jan. 2003): Revision and formatting
*
********************************************************************************

      IMPLICIT NONE
      CHARACTER*10 TYPENAME(40)
      INTEGER*4 IOPRHED

C______________________________________ State variable types

      IF (IOPRHED .EQ. 1) THEN
          TYPENAME(1) = 'Pressure  '
      ELSE
          TYPENAME(1) = 'Head      '
      ENDIF

      TYPENAME(2) = 'Concent.  '
      TYPENAME(3) = 'Humidity  '
      TYPENAME(4) = 'Flow      '
      TYPENAME(5) = 'Void      '
      TYPENAME(6) = 'Void      '
      TYPENAME(7) = 'Void      '
      TYPENAME(8) = 'Void      '
      TYPENAME(9) = 'Void      '
      TYPENAME(10) = 'Void      '

C______________________________________ Parameter types

      TYPENAME(11) = 'Txx       '     
      TYPENAME(12) = 'Tyy       '     
      TYPENAME(13) = 'Txy       '     
      TYPENAME(14) = 'Tzz       ' 
      TYPENAME(15) = 'Txz       '    
      TYPENAME(16) = 'Tyz       '       
      TYPENAME(17) = 'Storage   '    
      TYPENAME(18) = 'Recharge  '  
      TYPENAME(19) = 'Presc.head'    
      TYPENAME(20) = 'Presc.flow'       
      TYPENAME(21) = 'Leakage   '     
      TYPENAME(22) = 'Long.disp.'    
      TYPENAME(23) = 'Tran.disp.'  
      TYPENAME(24) = 'Mol. diff.'    
      TYPENAME(25) = 'Porosity  '     
      TYPENAME(26) = '1st decay ' 
      TYPENAME(27) = 'Retardatio'    
      TYPENAME(28) = 'Ext. conc.'    
      TYPENAME(29) = 'Gen. par. ' 
      TYPENAME(30) = 'Age coeff '
      TYPENAME(31) = 'Void      '
      TYPENAME(32) = 'Void      '
      TYPENAME(33) = 'Void      '
      TYPENAME(34) = 'Void      '
      TYPENAME(35) = 'Void      '
      TYPENAME(36) = 'Void      '
      TYPENAME(37) = 'Void      '
      TYPENAME(38) = 'Void      '
      TYPENAME(39) = 'Void      '
      TYPENAME(40) = 'Void      '
      
      RETURN
      END
