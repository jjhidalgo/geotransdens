      Subroutine Write_VTK_Unstructured_Grid
     &          (COORD    ,IOWri1D  ,KXX      ,LMXNDL   ,LNNDEL   ,LTYPE
     &          ,NUMEL    ,NUMNP    ,Title    ,VTKUnit)

C###############################################################################
C
C     Escribe una malla en formto VTK Legacy
C
C
C     Autor: JHG
C     Fecha: Julio, 2007
C            Junio, 2020 - Espacios en blanco elimininados.
C###############################################################################


      Implicit None

C------------------------- External

      Integer*4::IOWri1D  ,LMXNDL   ,NUMEL    ,NUMNP 

      Integer*4::VTKUnit

      Integer*4::GetVTKCellType

      Character*9::Title

      Integer*4::KXX(LMXNDL,NUMEL),LNNDEL(NUMEL),LTYPE(NUMEL)

      Real*8::COORD(NUMNP,3)

C------------------------- Internal

      Integer*4::CellType,I,L,ListSize,NCells,Num1D,NNUD

      Character::strFmt1*20

      Integer*4::Conect(6)

      Logical::Mask1D(NUMEL)

C------------------------- 
C-------------------------


C------------------------- Escribe la cabecera

      NUM1D = 0
      Mask1d = .FALSE.

      Write(VTKUnit,10) Title
      Write(*,10) Title
   10 Format('# vtk DataFile Version 1.0'
     &      ,/,A9
     &      ,/,'ASCII'
     &      ,/,'DATASET UNSTRUCTURED_GRID')

C-------------------------
C------------------------- NUDOS
C------------------------- 

C------------------------- Escribe los nudos


      Write(VTKUnit,20) NUMNP
   20 Format('POINTS',1X,I6,1X,'double')

      Do I=1,NUMNP

          Write(VTKUnit,30) Coord(I,1),Coord(I,2),Coord(I,3)
   30     Format (3(F15.5,1X))

      End Do !I=1,NUMNP

C-------------------------
C------------------------- Celdas
C------------------------- 

C------------------------- Calcula los números que hacen falta para
C------------------------- describir las celdas.

      ListSize = SUM(LNNDEL) + NUMEL

      If (IOWri1D.LT.1) Then

          Mask1D = LTYPE.EQ.1
          Num1D = COUNT(Mask1D)
          ListSize = ListSize - 3*Num1D

      End If

      NCells = NUMEL - (1-IOWri1D)*Num1D

C------------------------- Escribe las celdas (conectividades)
     
      Write(VTKUnit,40) NCells,ListSize
   40 Format ('CELLS',1X,I6,1X,I9)

      

      Do L=1,NUMEL

          NNUD = LNNDEL(L)

          Conect(1:NNUD) = KXX(1:NNUD,L) - 1 !Los nudos se numeran a partir de 0 en vtk

          If (NNUD.NE.2 .OR. (NNUD.EQ.2 .AND. IOWri1D.GT.0) )Then
              strFmt1 = ''
              Write(strFmt1,*) NNUD
              strFmt1 = '(I5,'//Trim(AdjustL(strFmt1))//'(1X,I6))'

              Write(VTKUnit,strFmt1) NNUD,Conect(1:NNUD)
          End If
C   50     Format(I5,<NNUD>(1X,I6))

      End Do !L=1,NUMEL

C------------------------- Escribe el tipo de celdas


      Write(VTKUnit,60) NUMEL - (1-IOWri1D)*Num1D
      
   60 Format(/,'CELL_TYPES',1X,I6)


      Do L=1,NUMEL
      
          If (LTYPE(L).NE.1 .OR. (LTYPE(L).EQ.1 .AND. IOWri1D.GT.0))Then

              CellType = GetVTKCellType(LTYPE(L))
              Write(VTKUnit,70) CellType
   70         Format(I5)

         End If

      End Do !L=1,NUMEL

      End Subroutine Write_VTK_Unstructured_Grid
