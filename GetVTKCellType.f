      Integer*4 Function GetVTKCellType(LType) Result (VTKCellType)

C###############################################################################
C
C     Convierte el typo de elemanto de Transdens en el tipo de celda de VTK.
C
C
C     LTYPE(L)=1 : 1-D Element (line in space)         --> VTKCellType = 3
C     LTYPE(L)=2 : 2-D Element (triangular)            --> VTKCellType = 5
C     LTYPE(L)=3 : 2-D Element (cuadrangular)          --> VTKCellType = 9
C     LTYPE(L)=4 : 3-D Element (tetrahedron)           --> VTKCellType = 10
C     LTYPE(L)=5 : 2-D Element (triangle in space)     --> VTKCellType = 5
C     LTYPE(L)=6 : 3-D Element (prismatic)             --> VTKCellType = 13
C
C     Autor: JHG
C     Fecha: Julio, 2007
C
C###############################################################################

      Implicit None

C------------------------- External

      Integer*4::LType !,VTKCellType

C-------------------------
C-------------------------
      Select Case (LType)

          Case(1)

              VTKCellType = 3

          Case(2,5,10)

              VTKCellType = 5

          Case(3)

              VTKCellType = 9

          Case(4)

              VTKCellType = 10

          Case(6,11)

              VTKCellType = 13

      End Select !LType

      End Function GetVTKCellType