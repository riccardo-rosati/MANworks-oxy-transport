# state file generated using paraview version 5.7.0-RC1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# trace generated using paraview version 5.7.0-RC1
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [777, 782]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [1.1, -1.9, 0.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [1.1, -1.9, 12.638425181333009]
renderView1.CameraFocalPoint = [1.1, -1.9, 0.5]
renderView1.CameraParallelScale = 3.1416556144810017
renderView1.Background = [0.32, 0.34, 0.43]
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024L
spreadSheetView1.HiddenColumnLabels = ['Point ID', 'Ct', 'D', 'K', 'LET', 'PO2_t', 'Points', 'Points_Magnitude', 'a1', 'a2', 'a3', 'a4', 'alpha', 'b', 'b1', 'b2', 'beta']
# uncomment following to set a specific view size
# spreadSheetView1.ViewSize = [400, 400]

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.500000)
layout1.AssignView(1, renderView1)
layout1.AssignView(2, spreadSheetView1)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
t_V8_Ctvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Tv8\\Ct_3d1d.vtk'])

# create a new 'Legacy VTK Reader'
t_V13_Ctvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Tv13\\Ct_3d1d.vtk'])

# create a new 'Legacy VTK Reader'
t_V18_Cvvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\v18\\Cv_3d1d.vtk'])

# create a new 'Legacy VTK Reader'
h_V2_Ctvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Hv2\\Ct_3d1d.vtk'])

# create a new 'Calculator'
h_V2_PO2_t = Calculator(Input=h_V2_Ctvtk)
h_V2_PO2_t.ResultArrayName = 'PO2_t'
h_V2_PO2_t.Function = 'Ct/0.0000389'

# create a new 'Calculator'
t_V8_PO2_t = Calculator(Input=t_V8_Ctvtk)
t_V8_PO2_t.ResultArrayName = 'PO2_t'
t_V8_PO2_t.Function = 'Ct/0.0000389'

# create a new 'Legacy VTK Reader'
t_V8_Cvvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Tv8\\Cv_3d1d.vtk'])

# create a new 'Legacy VTK Reader'
h_V2_Cvvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Hv2\\Cv_3d1d.vtk'])

# create a new 'Legacy VTK Reader'
t_V36_Cvvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Tv36\\Cv_3d1d.vtk'])

# create a new 'Calculator'
t_V13_PO2_t = Calculator(Input=t_V13_Ctvtk)
t_V13_PO2_t.ResultArrayName = 'PO2_t'
t_V13_PO2_t.Function = 'Ct/0.0000389'

# create a new 'Clip'
clip7 = Clip(Input=t_V13_PO2_t)
clip7.ClipType = 'Plane'
clip7.Scalars = ['POINTS', 'PO2_t']
clip7.Value = 15.025673550901988

# init the 'Plane' selected for 'ClipType'
clip7.ClipType.Origin = [0.5, 0.5, 0.5]
clip7.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
t_V36_PO2_v = Calculator(Input=t_V36_Cvvtk)
t_V36_PO2_v.ResultArrayName = 'PO2_v'
t_V36_PO2_v.Function = 'Cv/0.00003'

# create a new 'Calculator'
t_V18_PO2_v = Calculator(Input=t_V18_Cvvtk)
t_V18_PO2_v.ResultArrayName = 'PO2_v'
t_V18_PO2_v.Function = 'Cv/0.00003'

# create a new 'Legacy VTK Reader'
t_V2_Cvvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\v2\\Cv_3d1d.vtk'])

# create a new 'Legacy VTK Reader'
h_V8_Ctvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Hv8\\Ct_3d1d.vtk'])

# create a new 'Calculator'
h_V8_PO2_t = Calculator(Input=h_V8_Ctvtk)
h_V8_PO2_t.ResultArrayName = 'PO2_t'
h_V8_PO2_t.Function = 'Ct/0.0000389'

# create a new 'Calculator'
a1 = Calculator(Input=h_V8_PO2_t)
a1.ResultArrayName = 'a1'
a1.Function = '0.5856'

# create a new 'Calculator'
a1_1 = Calculator(Input=t_V8_PO2_t)
a1_1.ResultArrayName = 'a1'
a1_1.Function = '0.179'

# create a new 'Calculator'
h_V2_PO2_v = Calculator(Input=h_V2_Cvvtk)
h_V2_PO2_v.ResultArrayName = 'PO2_v'
h_V2_PO2_v.Function = 'Cv/0.00003'

# create a new 'Calculator'
a1_2 = Calculator(Input=t_V13_PO2_t)
a1_2.ResultArrayName = 'a1'
a1_2.Function = '0.179'

# create a new 'Legacy VTK Reader'
t_V18_Ctvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\v18\\Ct_3d1d.vtk'])

# create a new 'Calculator'
t_V18_PO2_t = Calculator(Input=t_V18_Ctvtk)
t_V18_PO2_t.ResultArrayName = 'PO2_t'
t_V18_PO2_t.Function = 'Ct/0.0000389'

# create a new 'Clip'
clip8 = Clip(Input=t_V18_PO2_t)
clip8.ClipType = 'Plane'
clip8.Scalars = ['POINTS', 'PO2_t']
clip8.Value = 21.633607299326844

# init the 'Plane' selected for 'ClipType'
clip8.ClipType.Origin = [0.5, 0.5, 0.5]
clip8.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Legacy VTK Reader'
h_V8_Cvdvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Hv8\\Cv_3d1d.vtk'])

# create a new 'Calculator'
h_V8_PO2_v = Calculator(Input=h_V8_Cvdvtk)
h_V8_PO2_v.ResultArrayName = 'PO2_v'
h_V8_PO2_v.Function = 'Cv/0.00003'

# create a new 'Calculator'
a1_3 = Calculator(Input=h_V2_PO2_t)
a1_3.ResultArrayName = 'a1'
a1_3.Function = '0.5856'

# create a new 'Legacy VTK Reader'
t_V2_Ctvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\v2\\Ct_3d1d.vtk'])

# create a new 'Calculator'
t_V8_PO2_v = Calculator(Input=t_V8_Cvvtk)
t_V8_PO2_v.ResultArrayName = 'PO2_v'
t_V8_PO2_v.Function = 'Cv/0.00003'

# create a new 'Calculator'
a2 = Calculator(Input=a1_2)
a2.ResultArrayName = 'a2'
a2.Function = '0.0024'

# create a new 'Calculator'
a3 = Calculator(Input=a2)
a3.ResultArrayName = 'a3'
a3.Function = '0.05'

# create a new 'Clip'
clip2 = Clip(Input=h_V8_PO2_t)
clip2.ClipType = 'Plane'
clip2.Scalars = ['POINTS', 'PO2_t']
clip2.Value = 41.18890950461051

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Origin = [0.5, 0.5, 0.5]
clip2.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Legacy VTK Reader'
h_V13_Ctvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Hv13\\Ct_3d1d.vtk'])

# create a new 'Legacy VTK Reader'
t_V36_Ctvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Tv36\\Ct_3d1d.vtk'])

# create a new 'Calculator'
t_V36_PO2_t = Calculator(Input=t_V36_Ctvtk)
t_V36_PO2_t.ResultArrayName = 'PO2_t'
t_V36_PO2_t.Function = 'Ct/0.0000389'

# create a new 'Clip'
clip18 = Clip(Input=t_V36_PO2_t)
clip18.ClipType = 'Plane'
clip18.Scalars = ['POINTS', 'PO2_t']
clip18.Value = 42.57696769308002

# init the 'Plane' selected for 'ClipType'
clip18.ClipType.Origin = [0.5, 0.5, 0.5]
clip18.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
h_V13_PO2_t = Calculator(Input=h_V13_Ctvtk)
h_V13_PO2_t.ResultArrayName = 'PO2_t'
h_V13_PO2_t.Function = 'Ct/0.0000389'

# create a new 'Clip'
clip3 = Clip(Input=h_V13_PO2_t)
clip3.ClipType = 'Plane'
clip3.Scalars = ['POINTS', 'PO2_t']
clip3.Value = 53.25128307570836

# init the 'Plane' selected for 'ClipType'
clip3.ClipType.Origin = [0.5, 0.5, 0.5]
clip3.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
a1_4 = Calculator(Input=h_V13_PO2_t)
a1_4.ResultArrayName = 'a1'
a1_4.Function = '0.5856'

# create a new 'Calculator'
a2_1 = Calculator(Input=a1_4)
a2_1.ResultArrayName = 'a2'
a2_1.Function = '0.0024'

# create a new 'Calculator'
a3_1 = Calculator(Input=a2_1)
a3_1.ResultArrayName = 'a3'
a3_1.Function = '0.05'

# create a new 'Calculator'
a4 = Calculator(Input=a3_1)
a4.ResultArrayName = 'a4'
a4.Function = '0.0031'

# create a new 'Calculator'
b1 = Calculator(Input=a4)
b1.ResultArrayName = 'b1'
b1.Function = '0.2564'

# create a new 'Calculator'
b2 = Calculator(Input=b1)
b2.ResultArrayName = 'b2'
b2.Function = '0.015'

# create a new 'Calculator'
lET = Calculator(Input=b2)
lET.ResultArrayName = 'LET'
lET.Function = '2'

# create a new 'Calculator'
k = Calculator(Input=lET)
k.ResultArrayName = 'K'
k.Function = '2'

# create a new 'Calculator'
a2_2 = Calculator(Input=a1)
a2_2.ResultArrayName = 'a2'
a2_2.Function = '0.0024'

# create a new 'Legacy VTK Reader'
h_V13_Cvvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Hv13\\Cv_3d1d.vtk'])

# create a new 'Calculator'
h_V13_PO2_v = Calculator(Input=h_V13_Cvvtk)
h_V13_PO2_v.ResultArrayName = 'PO2_v'
h_V13_PO2_v.Function = 'Cv/0.00003'

# create a new 'Legacy VTK Reader'
t_V13_Cvvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Tv13\\Cv_3d1d.vtk'])

# create a new 'Calculator'
t_V13_PO2_v = Calculator(Input=t_V13_Cvvtk)
t_V13_PO2_v.ResultArrayName = 'PO2_v'
t_V13_PO2_v.Function = 'Cv/0.00003'

# create a new 'Calculator'
t_V2_PO2_t = Calculator(Input=t_V2_Ctvtk)
t_V2_PO2_t.ResultArrayName = 'PO2_t'
t_V2_PO2_t.Function = 'Ct/0.0000389'

# create a new 'Clip'
clip5 = Clip(Input=t_V2_PO2_t)
clip5.ClipType = 'Plane'
clip5.Scalars = ['POINTS', 'PO2_t']
clip5.Value = 2.8719722553979907

# init the 'Plane' selected for 'ClipType'
clip5.ClipType.Origin = [0.5, 0.5, 0.5]
clip5.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
a1_5 = Calculator(Input=t_V2_PO2_t)
a1_5.ResultArrayName = 'a1'
a1_5.Function = '0.179'

# create a new 'Calculator'
a2_3 = Calculator(Input=a1_5)
a2_3.ResultArrayName = 'a2'
a2_3.Function = '0.0024'

# create a new 'Calculator'
a3_2 = Calculator(Input=a2_3)
a3_2.ResultArrayName = 'a3'
a3_2.Function = '0.05'

# create a new 'Calculator'
a4_1 = Calculator(Input=a3_2)
a4_1.ResultArrayName = 'a4'
a4_1.Function = '0.0031'

# create a new 'Calculator'
b1_1 = Calculator(Input=a4_1)
b1_1.ResultArrayName = 'b1'
b1_1.Function = '0.22'

# create a new 'Calculator'
b2_1 = Calculator(Input=b1_1)
b2_1.ResultArrayName = 'b2'
b2_1.Function = '0.015'

# create a new 'Legacy VTK Reader'
h_V18_Cvvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Hv18\\Cv_3d1d.vtk'])

# create a new 'Calculator'
h_V18_PO2_v = Calculator(Input=h_V18_Cvvtk)
h_V18_PO2_v.ResultArrayName = 'PO2_v'
h_V18_PO2_v.Function = 'Cv/0.00003'

# create a new 'Calculator'
a3_3 = Calculator(Input=a2_2)
a3_3.ResultArrayName = 'a3'
a3_3.Function = '0.05'

# create a new 'Calculator'
a4_2 = Calculator(Input=a3_3)
a4_2.ResultArrayName = 'a4'
a4_2.Function = '0.0031'

# create a new 'Calculator'
b1_2 = Calculator(Input=a4_2)
b1_2.ResultArrayName = 'b1'
b1_2.Function = '0.2564'

# create a new 'Calculator'
b2_2 = Calculator(Input=b1_2)
b2_2.ResultArrayName = 'b2'
b2_2.Function = '0.015'

# create a new 'Calculator'
lET_1 = Calculator(Input=b2_2)
lET_1.ResultArrayName = 'LET'
lET_1.Function = '2'

# create a new 'Calculator'
k_1 = Calculator(Input=lET_1)
k_1.ResultArrayName = 'K'
k_1.Function = '2.5'

# create a new 'Calculator'
aLPHA = Calculator(Input=k_1)
aLPHA.ResultArrayName = 'alpha'
aLPHA.Function = '((a1+a2*LET)*PO2_t + (a3+a4*LET)*K)/(PO2_t+K)'

# create a new 'Calculator'
b = Calculator(Input=aLPHA)
b.ResultArrayName = 'b'
b.Function = '(b1*PO2_t + b2*K)/(PO2_t+K)'

# create a new 'Clip'
clip1 = Clip(Input=h_V2_PO2_t)
clip1.ClipType = 'Plane'
clip1.Scalars = ['POINTS', 'PO2_t']
clip1.Value = 5.7571056356109205

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [0.5, 0.5, 0.5]
clip1.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
bETA = Calculator(Input=b)
bETA.ResultArrayName = 'beta'
bETA.Function = 'b*b'

# create a new 'Calculator'
d = Calculator(Input=bETA)
d.ResultArrayName = 'D'
d.Function = '2'

# create a new 'Calculator'
sf_H_V8 = Calculator(Input=d)
sf_H_V8.ResultArrayName = 'Sf'
sf_H_V8.Function = 'exp(-alpha*D-beta*D*D)'

# create a new 'Calculator'
alphabeta = Calculator(Input=sf_H_V8)
alphabeta.ResultArrayName = 'a_b'
alphabeta.Function = 'alpha/beta'

# create a new 'Clip'
clip10 = Clip(Input=sf_H_V8)
clip10.ClipType = 'Plane'
clip10.Scalars = ['POINTS', 'Sf']
clip10.Value = 0.36887048617417795

# init the 'Plane' selected for 'ClipType'
clip10.ClipType.Origin = [0.5, 0.5, 0.5]
clip10.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Clip'
clip6 = Clip(Input=t_V8_PO2_t)
clip6.ClipType = 'Plane'
clip6.Scalars = ['POINTS', 'PO2_t']
clip6.Value = 6.9507447071822375

# init the 'Plane' selected for 'ClipType'
clip6.ClipType.Origin = [0.5, 0.5, 0.5]
clip6.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
t_V2_PO2_v = Calculator(Input=t_V2_Cvvtk)
t_V2_PO2_v.ResultArrayName = 'PO2_v'
t_V2_PO2_v.Function = 'Cv/0.00003'

# create a new 'Calculator'
a4_3 = Calculator(Input=a3)
a4_3.ResultArrayName = 'a4'
a4_3.Function = '0.0031'

# create a new 'Calculator'
b1_3 = Calculator(Input=a4_3)
b1_3.ResultArrayName = 'b1'
b1_3.Function = '0.22'

# create a new 'Calculator'
b2_3 = Calculator(Input=b1_3)
b2_3.ResultArrayName = 'b2'
b2_3.Function = '0.015'

# create a new 'Calculator'
lET_2 = Calculator(Input=b2_3)
lET_2.ResultArrayName = 'LET'
lET_2.Function = '2'

# create a new 'Calculator'
k_2 = Calculator(Input=lET_2)
k_2.ResultArrayName = 'K'
k_2.Function = '2.5'

# create a new 'Calculator'
aLPHA_1 = Calculator(Input=k_2)
aLPHA_1.ResultArrayName = 'alpha'
aLPHA_1.Function = '((a1+a2*LET)*PO2_t + (a3+a4*LET)*K)/(PO2_t+K)'

# create a new 'Calculator'
b_1 = Calculator(Input=aLPHA_1)
b_1.ResultArrayName = 'b'
b_1.Function = '(b1*PO2_t + b2*K)/(PO2_t+K)'

# create a new 'Calculator'
bETA_1 = Calculator(Input=b_1)
bETA_1.ResultArrayName = 'beta'
bETA_1.Function = 'b*b'

# create a new 'Calculator'
d_1 = Calculator(Input=bETA_1)
d_1.ResultArrayName = 'D'
d_1.Function = '50'

# create a new 'Calculator'
sf_T_V13 = Calculator(Input=d_1)
sf_T_V13.ResultArrayName = 'Sf'
sf_T_V13.Function = 'exp(-alpha*D-beta*D*D)'

# create a new 'Integrate Variables'
integrateVariables2 = IntegrateVariables(Input=sf_T_V13)

# create a new 'Calculator'
calculator3 = Calculator(Input=integrateVariables2)
calculator3.ResultArrayName = 'N'
calculator3.Function = '375'

# create a new 'Calculator'
calculator4 = Calculator(Input=calculator3)
calculator4.ResultArrayName = 'TCP'
calculator4.Function = 'exp(-N*Sf)'

# create a new 'Clip'
clip15 = Clip(Input=sf_T_V13)
clip15.ClipType = 'Plane'
clip15.Scalars = ['POINTS', 'Sf']
clip15.Value = 0.6363680204885862

# init the 'Plane' selected for 'ClipType'
clip15.ClipType.Origin = [0.5, 0.5, 0.5]
clip15.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
alphabeta_1 = Calculator(Input=clip15)
alphabeta_1.ResultArrayName = 'a_b'
alphabeta_1.Function = 'alpha/beta'

# create a new 'Clip'
clip24 = Clip(Input=alphabeta_1)
clip24.ClipType = 'Plane'
clip24.Scalars = ['POINTS', 'a_b']
clip24.Value = 13.851767814230605

# init the 'Plane' selected for 'ClipType'
clip24.ClipType.Origin = [0.5, 0.5, 0.5]
clip24.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Legacy VTK Reader'
h_V18_Ctvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Hv18\\Ct_3d1d.vtk'])

# create a new 'Calculator'
h_V18_PO2_t = Calculator(Input=h_V18_Ctvtk)
h_V18_PO2_t.ResultArrayName = 'PO2_t'
h_V18_PO2_t.Function = 'Ct/0.0000389'

# create a new 'Calculator'
a1_6 = Calculator(Input=h_V18_PO2_t)
a1_6.ResultArrayName = 'a1'
a1_6.Function = '0.5856'

# create a new 'Calculator'
a2_4 = Calculator(Input=a1_6)
a2_4.ResultArrayName = 'a2'
a2_4.Function = '0.0024'

# create a new 'Calculator'
a3_4 = Calculator(Input=a2_4)
a3_4.ResultArrayName = 'a3'
a3_4.Function = '0.05'

# create a new 'Calculator'
a4_4 = Calculator(Input=a3_4)
a4_4.ResultArrayName = 'a4'
a4_4.Function = '0.0031'

# create a new 'Calculator'
b1_4 = Calculator(Input=a4_4)
b1_4.ResultArrayName = 'b1'
b1_4.Function = '0.2564'

# create a new 'Calculator'
b2_4 = Calculator(Input=b1_4)
b2_4.ResultArrayName = 'b2'
b2_4.Function = '0.015'

# create a new 'Calculator'
lET_3 = Calculator(Input=b2_4)
lET_3.ResultArrayName = 'LET'
lET_3.Function = '2'

# create a new 'Calculator'
k_3 = Calculator(Input=lET_3)
k_3.ResultArrayName = 'K'
k_3.Function = '2.5'

# create a new 'Calculator'
aLPHA_2 = Calculator(Input=k_3)
aLPHA_2.ResultArrayName = 'alpha'
aLPHA_2.Function = '((a1+a2*LET)*PO2_t + (a3+a4*LET)*K)/(PO2_t+K)'

# create a new 'Clip'
clip4 = Clip(Input=h_V18_PO2_t)
clip4.ClipType = 'Plane'
clip4.Scalars = ['POINTS', 'PO2_t']
clip4.Value = 57.5796616021879

# init the 'Plane' selected for 'ClipType'
clip4.ClipType.Origin = [0.5, 0.5, 0.5]
clip4.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Clip'
clip20 = Clip(Input=alphabeta)
clip20.ClipType = 'Plane'
clip20.Scalars = ['POINTS', 'a_b']
clip20.Value = 12.584116227780822

# init the 'Plane' selected for 'ClipType'
clip20.ClipType.Origin = [0.5, 0.5, 0.5]
clip20.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
lET_4 = Calculator(Input=b2_1)
lET_4.ResultArrayName = 'LET'
lET_4.Function = '2'

# create a new 'Calculator'
k_4 = Calculator(Input=lET_4)
k_4.ResultArrayName = 'K'
k_4.Function = '2.5'

# create a new 'Calculator'
aLPHA_3 = Calculator(Input=k_4)
aLPHA_3.ResultArrayName = 'alpha'
aLPHA_3.Function = '((a1+a2*LET)*PO2_t + (a3+a4*LET)*K)/(PO2_t+K)'

# create a new 'Calculator'
a2_5 = Calculator(Input=a1_3)
a2_5.ResultArrayName = 'a2'
a2_5.Function = '0.0024'

# create a new 'Calculator'
a3_5 = Calculator(Input=a2_5)
a3_5.ResultArrayName = 'a3'
a3_5.Function = '0.05'

# create a new 'Calculator'
a4_5 = Calculator(Input=a3_5)
a4_5.ResultArrayName = 'a4'
a4_5.Function = '0.0031'

# create a new 'Calculator'
b1_5 = Calculator(Input=a4_5)
b1_5.ResultArrayName = 'b1'
b1_5.Function = '0.2564'

# create a new 'Calculator'
b2_5 = Calculator(Input=b1_5)
b2_5.ResultArrayName = 'b2'
b2_5.Function = '0.051'

# create a new 'Calculator'
lET_5 = Calculator(Input=b2_5)
lET_5.ResultArrayName = 'LET'
lET_5.Function = '2'

# create a new 'Calculator'
k_5 = Calculator(Input=lET_5)
k_5.ResultArrayName = 'K'
k_5.Function = '2.5'

# create a new 'Calculator'
aLPHA_4 = Calculator(Input=k_5)
aLPHA_4.ResultArrayName = 'alpha'
aLPHA_4.Function = '((a1+a2*LET)*PO2_t + (a3+a4*LET)*K)/(PO2_t+K)'

# create a new 'Calculator'
b_2 = Calculator(Input=aLPHA_4)
b_2.ResultArrayName = 'b'
b_2.Function = '(b1*PO2_t + b2*K)/(PO2_t+K)'

# create a new 'Calculator'
bETA_2 = Calculator(Input=b_2)
bETA_2.ResultArrayName = 'beta'
bETA_2.Function = 'b*b'

# create a new 'Calculator'
d_2 = Calculator(Input=bETA_2)
d_2.ResultArrayName = 'D'
d_2.Function = '2'

# create a new 'Calculator'
sf_H_V2 = Calculator(Input=d_2)
sf_H_V2.ResultArrayName = 'Sf'
sf_H_V2.Function = 'exp(-alpha*D-beta*D*D)'

# create a new 'Clip'
clip9 = Clip(Input=sf_H_V2)
clip9.ClipType = 'Plane'
clip9.Scalars = ['POINTS', 'Sf']
clip9.Value = 0.6542898049566793

# init the 'Plane' selected for 'ClipType'
clip9.ClipType.Origin = [0.5, 0.5, 0.5]
clip9.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
alphabeta_2 = Calculator(Input=sf_H_V2)
alphabeta_2.ResultArrayName = 'a_b'
alphabeta_2.Function = 'alpha/beta'

# create a new 'Clip'
clip19 = Clip(Input=alphabeta_2)
clip19.ClipType = 'Plane'
clip19.Scalars = ['POINTS', 'a_b']
clip19.Value = 18.649524100380344

# init the 'Plane' selected for 'ClipType'
clip19.ClipType.Origin = [0.5, 0.5, 0.5]
clip19.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
a1_7 = Calculator(Input=t_V36_PO2_t)
a1_7.ResultArrayName = 'a1'
a1_7.Function = '0.179'

# create a new 'Calculator'
a2_6 = Calculator(Input=a1_7)
a2_6.ResultArrayName = 'a2'
a2_6.Function = '0.0024'

# create a new 'Calculator'
a3_6 = Calculator(Input=a2_6)
a3_6.ResultArrayName = 'a3'
a3_6.Function = '0.05'

# create a new 'Calculator'
a4_6 = Calculator(Input=a3_6)
a4_6.ResultArrayName = 'a4'
a4_6.Function = '0.0031'

# create a new 'Calculator'
b1_6 = Calculator(Input=a4_6)
b1_6.ResultArrayName = 'b1'
b1_6.Function = '0.22'

# create a new 'Calculator'
b2_6 = Calculator(Input=b1_6)
b2_6.ResultArrayName = 'b2'
b2_6.Function = '0.015'

# create a new 'Calculator'
lET_6 = Calculator(Input=b2_6)
lET_6.ResultArrayName = 'LET'
lET_6.Function = '2'

# create a new 'Calculator'
k_6 = Calculator(Input=lET_6)
k_6.ResultArrayName = 'K'
k_6.Function = '2.5'

# create a new 'Calculator'
aLPHA_5 = Calculator(Input=k_6)
aLPHA_5.ResultArrayName = 'alpha'
aLPHA_5.Function = '((a1+a2*LET)*PO2_t + (a3+a4*LET)*K)/(PO2_t+K)'

# create a new 'Calculator'
b_3 = Calculator(Input=aLPHA_5)
b_3.ResultArrayName = 'b'
b_3.Function = '(b1*PO2_t + b2*K)/(PO2_t+K)'

# create a new 'Calculator'
bETA_3 = Calculator(Input=b_3)
bETA_3.ResultArrayName = 'beta'
bETA_3.Function = 'b*b'

# create a new 'Calculator'
d_3 = Calculator(Input=bETA_3)
d_3.ResultArrayName = 'D'
d_3.Function = '20'

# create a new 'Calculator'
sf_T_V36 = Calculator(Input=d_3)
sf_T_V36.ResultArrayName = 'Sf'
sf_T_V36.Function = 'exp(-alpha*D-beta*D*D)'

# create a new 'Clip'
clip17 = Clip(Input=sf_T_V36)
clip17.ClipType = 'Plane'
clip17.Scalars = ['POINTS', 'Sf']
clip17.Value = 0.37064291859081183

# init the 'Plane' selected for 'ClipType'
clip17.ClipType.Origin = [0.5, 0.5, 0.5]
clip17.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
alphabeta_3 = Calculator(Input=sf_T_V36)
alphabeta_3.ResultArrayName = 'a_b'
alphabeta_3.Function = 'alpha/beta'

# create a new 'Clip'
clip27 = Clip(Input=alphabeta_3)
clip27.ClipType = 'Plane'
clip27.Scalars = ['POINTS', 'a_b']
clip27.Value = 10.706793440623288

# init the 'Plane' selected for 'ClipType'
clip27.ClipType.Origin = [0.5, 0.5, 0.5]
clip27.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Integrate Variables'
integrateVariables5 = IntegrateVariables(Input=sf_T_V36)

# create a new 'Calculator'
calculator9 = Calculator(Input=integrateVariables5)
calculator9.ResultArrayName = 'N'
calculator9.Function = '375'

# create a new 'Calculator'
calculator10 = Calculator(Input=calculator9)
calculator10.ResultArrayName = 'TCP'
calculator10.Function = 'exp(-N*Sf)'

# create a new 'Calculator'
a2_7 = Calculator(Input=a1_1)
a2_7.ResultArrayName = 'a2'
a2_7.Function = '0.0024'

# create a new 'Calculator'
a3_7 = Calculator(Input=a2_7)
a3_7.ResultArrayName = 'a3'
a3_7.Function = '0.05'

# create a new 'Calculator'
a4_7 = Calculator(Input=a3_7)
a4_7.ResultArrayName = 'a4'
a4_7.Function = '0.0031'

# create a new 'Calculator'
b1_7 = Calculator(Input=a4_7)
b1_7.ResultArrayName = 'b1'
b1_7.Function = '0.22'

# create a new 'Calculator'
b2_7 = Calculator(Input=b1_7)
b2_7.ResultArrayName = 'b2'
b2_7.Function = '0.015'

# create a new 'Calculator'
lET_7 = Calculator(Input=b2_7)
lET_7.ResultArrayName = 'LET'
lET_7.Function = '2'

# create a new 'Calculator'
k_7 = Calculator(Input=lET_7)
k_7.ResultArrayName = 'K'
k_7.Function = '2.5'

# create a new 'Calculator'
aLPHA_6 = Calculator(Input=k_7)
aLPHA_6.ResultArrayName = 'alpha'
aLPHA_6.Function = '((a1+a2*LET)*PO2_t + (a3+a4*LET)*K)/(PO2_t+K)'

# create a new 'Calculator'
b_4 = Calculator(Input=aLPHA_6)
b_4.ResultArrayName = 'b'
b_4.Function = '(b1*PO2_t + b2*K)/(PO2_t+K)'

# create a new 'Calculator'
bETA_4 = Calculator(Input=b_4)
bETA_4.ResultArrayName = 'beta'
bETA_4.Function = 'b*b'

# create a new 'Calculator'
d_4 = Calculator(Input=bETA_4)
d_4.ResultArrayName = 'D'
d_4.Function = '50'

# create a new 'Calculator'
sf_T_V8 = Calculator(Input=d_4)
sf_T_V8.ResultArrayName = 'Sf'
sf_T_V8.Function = 'exp(-alpha*D-beta*D*D)'

# create a new 'Clip'
clip14 = Clip(Input=sf_T_V8)
clip14.ClipType = 'Plane'
clip14.Scalars = ['POINTS', 'Sf']
clip14.Value = 0.6622941659655256

# init the 'Plane' selected for 'ClipType'
clip14.ClipType.Origin = [0.5, 0.5, 0.5]
clip14.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Integrate Variables'
integrateVariables3 = IntegrateVariables(Input=sf_T_V8)

# create a new 'Calculator'
calculator5 = Calculator(Input=integrateVariables3)
calculator5.ResultArrayName = 'N'
calculator5.Function = '375'

# create a new 'Calculator'
calculator6 = Calculator(Input=calculator5)
calculator6.ResultArrayName = 'TCP'
calculator6.Function = 'exp(-N*Sf)'

# create a new 'Calculator'
alphabeta_4 = Calculator(Input=sf_T_V8)
alphabeta_4.ResultArrayName = 'a_b'
alphabeta_4.Function = 'alpha/beta'

# create a new 'Clip'
clip25 = Clip(Input=alphabeta_4)
clip25.ClipType = 'Plane'
clip25.Scalars = ['POINTS', 'a_b']
clip25.Value = 232.940967044919

# init the 'Plane' selected for 'ClipType'
clip25.ClipType.Origin = [0.5, 0.5, 0.5]
clip25.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
b_5 = Calculator(Input=aLPHA_3)
b_5.ResultArrayName = 'b'
b_5.Function = '(b1*PO2_t + b2*K)/(PO2_t+K)'

# create a new 'Calculator'
bETA_5 = Calculator(Input=b_5)
bETA_5.ResultArrayName = 'beta'
bETA_5.Function = 'b*b'

# create a new 'Calculator'
d_5 = Calculator(Input=bETA_5)
d_5.ResultArrayName = 'D'
d_5.Function = '50'

# create a new 'Calculator'
sf_T_V2 = Calculator(Input=d_5)
sf_T_V2.ResultArrayName = 'Sf'
sf_T_V2.Function = 'exp(-alpha*D-beta*D*D)'

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=sf_T_V2)

# create a new 'Calculator'
calculator1 = Calculator(Input=integrateVariables1)
calculator1.ResultArrayName = 'N'
calculator1.Function = '375'

# create a new 'Calculator'
calculator2 = Calculator(Input=calculator1)
calculator2.ResultArrayName = 'TCP'
calculator2.Function = 'exp(-N*Sf)'

# create a new 'Clip'
clip13 = Clip(Input=sf_T_V2)
clip13.ClipType = 'Plane'
clip13.Scalars = ['POINTS', 'Sf']
clip13.Value = 0.7095083902950026

# init the 'Plane' selected for 'ClipType'
clip13.ClipType.Origin = [0.5, 0.5, 0.5]
clip13.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
alphabeta_5 = Calculator(Input=clip13)
alphabeta_5.ResultArrayName = 'a_b'
alphabeta_5.Function = 'alpha/beta'

# create a new 'Clip'
clip23 = Clip(Input=alphabeta_5)
clip23.ClipType = 'Plane'
clip23.Scalars = ['POINTS', 'a_b']
clip23.Value = 346.5350343892692

# init the 'Plane' selected for 'ClipType'
clip23.ClipType.Origin = [0.5, 0.5, 0.5]
clip23.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
aLPHA_7 = Calculator(Input=k)
aLPHA_7.ResultArrayName = 'alpha'
aLPHA_7.Function = '((a1+a2*LET)*PO2_t + (a3+a4*LET)*K)/(PO2_t+K)'

# create a new 'Calculator'
a1_8 = Calculator(Input=t_V18_PO2_t)
a1_8.ResultArrayName = 'a1'
a1_8.Function = '0.179'

# create a new 'Calculator'
a2_8 = Calculator(Input=a1_8)
a2_8.ResultArrayName = 'a2'
a2_8.Function = '0.0024'

# create a new 'Calculator'
a3_8 = Calculator(Input=a2_8)
a3_8.ResultArrayName = 'a3'
a3_8.Function = '0.05'

# create a new 'Calculator'
a4_8 = Calculator(Input=a3_8)
a4_8.ResultArrayName = 'a4'
a4_8.Function = '0.0031'

# create a new 'Calculator'
b1_8 = Calculator(Input=a4_8)
b1_8.ResultArrayName = 'b1'
b1_8.Function = '0.22'

# create a new 'Calculator'
b2_8 = Calculator(Input=b1_8)
b2_8.ResultArrayName = 'b2'
b2_8.Function = '0.015'

# create a new 'Calculator'
lET_8 = Calculator(Input=b2_8)
lET_8.ResultArrayName = 'LET'
lET_8.Function = '2'

# create a new 'Calculator'
k_8 = Calculator(Input=lET_8)
k_8.ResultArrayName = 'K'
k_8.Function = '2.5'

# create a new 'Calculator'
aLPHA_8 = Calculator(Input=k_8)
aLPHA_8.ResultArrayName = 'alpha'
aLPHA_8.Function = '((a1+a2*LET)*PO2_t + (a3+a4*LET)*K)/(PO2_t+K)'

# create a new 'Calculator'
b_6 = Calculator(Input=aLPHA_8)
b_6.ResultArrayName = 'b'
b_6.Function = '(b1*PO2_t + b2*K)/(PO2_t+K)'

# create a new 'Calculator'
bETA_6 = Calculator(Input=b_6)
bETA_6.ResultArrayName = 'beta'
bETA_6.Function = 'b*b'

# create a new 'Calculator'
d_6 = Calculator(Input=bETA_6)
d_6.ResultArrayName = 'D'
d_6.Function = '20'

# create a new 'Calculator'
sf_T_V18 = Calculator(Input=d_6)
sf_T_V18.ResultArrayName = 'Sf'
sf_T_V18.Function = 'exp(-alpha*D-beta*D*D)'

# create a new 'Integrate Variables'
integrateVariables4 = IntegrateVariables(Input=sf_T_V18)

# create a new 'Calculator'
calculator7 = Calculator(Input=integrateVariables4)
calculator7.ResultArrayName = 'N'
calculator7.Function = '375'

# create a new 'Calculator'
calculator8 = Calculator(Input=calculator7)
calculator8.ResultArrayName = 'TCP'
calculator8.Function = 'exp(-N*Sf)'

# create a new 'Clip'
clip16 = Clip(Input=sf_T_V18)
clip16.ClipType = 'Plane'
clip16.Scalars = ['POINTS', 'Sf']
clip16.Value = 0.5210542768132349

# init the 'Plane' selected for 'ClipType'
clip16.ClipType.Origin = [0.5, 0.5, 0.5]
clip16.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
alphabeta_6 = Calculator(Input=sf_T_V18)
alphabeta_6.ResultArrayName = 'a_b'
alphabeta_6.Function = 'alpha/beta'

# create a new 'Clip'
clip26 = Clip(Input=alphabeta_6)
clip26.ClipType = 'Plane'
clip26.Scalars = ['POINTS', 'a_b']
clip26.Value = 17.20463726860738

# init the 'Plane' selected for 'ClipType'
clip26.ClipType.Origin = [0.5, 0.5, 0.5]
clip26.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
b_7 = Calculator(Input=aLPHA_7)
b_7.ResultArrayName = 'b'
b_7.Function = '(b1*PO2_t + b2*K)/(PO2_t+K)'

# create a new 'Calculator'
bETA_7 = Calculator(Input=b_7)
bETA_7.ResultArrayName = 'beta'
bETA_7.Function = 'b*b'

# create a new 'Calculator'
d_7 = Calculator(Input=bETA_7)
d_7.ResultArrayName = 'D'
d_7.Function = '2'

# create a new 'Calculator'
sf_H_V13 = Calculator(Input=d_7)
sf_H_V13.ResultArrayName = 'Sf'
sf_H_V13.Function = 'exp(-alpha*D-beta*D*D)'

# create a new 'Calculator'
alphabeta_7 = Calculator(Input=sf_H_V13)
alphabeta_7.ResultArrayName = 'a_b'
alphabeta_7.Function = 'alpha/beta'

# create a new 'Clip'
clip21 = Clip(Input=alphabeta_7)
clip21.ClipType = 'Plane'
clip21.Scalars = ['POINTS', 'a_b']
clip21.Value = 12.244588865852243

# init the 'Plane' selected for 'ClipType'
clip21.ClipType.Origin = [0.5, 0.5, 0.5]
clip21.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Clip'
clip11 = Clip(Input=sf_H_V13)
clip11.ClipType = 'Plane'
clip11.Scalars = ['POINTS', 'Sf']
clip11.Value = 0.35611042992083813

# init the 'Plane' selected for 'ClipType'
clip11.ClipType.Origin = [0.5, 0.5, 0.5]
clip11.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
b_8 = Calculator(Input=aLPHA_2)
b_8.ResultArrayName = 'b'
b_8.Function = '(b1*PO2_t + b2*K)/(PO2_t+K)'

# create a new 'Calculator'
bETA_8 = Calculator(Input=b_8)
bETA_8.ResultArrayName = 'beta'
bETA_8.Function = 'b*b'

# create a new 'Calculator'
d_8 = Calculator(Input=bETA_8)
d_8.ResultArrayName = 'D'
d_8.Function = '2'

# create a new 'Calculator'
sf_H_V18 = Calculator(Input=d_8)
sf_H_V18.ResultArrayName = 'Sf'
sf_H_V18.Function = 'exp(-alpha*D-beta*D*D)'

# create a new 'Calculator'
alphabeta_8 = Calculator(Input=sf_H_V18)
alphabeta_8.ResultArrayName = 'a_b'
alphabeta_8.Function = 'alpha/beta'

# create a new 'Clip'
clip22 = Clip(Input=alphabeta_8)
clip22.ClipType = 'Plane'
clip22.Scalars = ['POINTS', 'a_b']
clip22.Value = 12.35753778474169

# init the 'Plane' selected for 'ClipType'
clip22.ClipType.Origin = [0.5, 0.5, 0.5]
clip22.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Clip'
clip12 = Clip(Input=sf_H_V18)
clip12.ClipType = 'Plane'
clip12.Scalars = ['POINTS', 'Sf']
clip12.Value = 0.3590366984772403

# init the 'Plane' selected for 'ClipType'
clip12.ClipType.Origin = [0.5, 0.5, 0.5]
clip12.ClipType.Normal = [0.5, 0.5, 0.5]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from h_V2_Ctvtk
h_V2_CtvtkDisplay = Show(h_V2_Ctvtk, renderView1)

# get color transfer function/color map for 'Ct'
ctLUT = GetColorTransferFunction('Ct')
ctLUT.RGBPoints = [-4.730674845632166e-06, 0.231373, 0.298039, 0.752941, 0.0012622122776519973, 0.865003, 0.865003, 0.865003, 0.0025291552301496267, 0.705882, 0.0156863, 0.14902]
ctLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Ct'
ctPWF = GetOpacityTransferFunction('Ct')
ctPWF.Points = [-4.730674845632166e-06, 0.0, 0.5, 0.0, 0.0025291552301496267, 1.0, 0.5, 0.0]
ctPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
h_V2_CtvtkDisplay.Representation = 'Outline'
h_V2_CtvtkDisplay.ColorArrayName = ['POINTS', 'Ct']
h_V2_CtvtkDisplay.LookupTable = ctLUT
h_V2_CtvtkDisplay.OSPRayScaleArray = 'Ct'
h_V2_CtvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
h_V2_CtvtkDisplay.SelectOrientationVectors = 'Ct'
h_V2_CtvtkDisplay.ScaleFactor = 0.1
h_V2_CtvtkDisplay.SelectScaleArray = 'Ct'
h_V2_CtvtkDisplay.GlyphType = 'Arrow'
h_V2_CtvtkDisplay.GlyphTableIndexArray = 'Ct'
h_V2_CtvtkDisplay.GaussianRadius = 0.005
h_V2_CtvtkDisplay.SetScaleArray = ['POINTS', 'Ct']
h_V2_CtvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
h_V2_CtvtkDisplay.OpacityArray = ['POINTS', 'Ct']
h_V2_CtvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
h_V2_CtvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
h_V2_CtvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
h_V2_CtvtkDisplay.ScalarOpacityFunction = ctPWF
h_V2_CtvtkDisplay.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
h_V2_CtvtkDisplay.ScaleTransferFunction.Points = [6.220206159923691e-06, 0.0, 0.5, 0.0, 0.0004416826122906059, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
h_V2_CtvtkDisplay.OpacityTransferFunction.Points = [6.220206159923691e-06, 0.0, 0.5, 0.0, 0.0004416826122906059, 1.0, 0.5, 0.0]

# show data from h_V8_Ctvtk
h_V8_CtvtkDisplay = Show(h_V8_Ctvtk, renderView1)

# trace defaults for the display properties.
h_V8_CtvtkDisplay.Representation = 'Outline'
h_V8_CtvtkDisplay.ColorArrayName = ['POINTS', 'Ct']
h_V8_CtvtkDisplay.LookupTable = ctLUT
h_V8_CtvtkDisplay.Position = [0.0, -1.2, 0.0]
h_V8_CtvtkDisplay.OSPRayScaleArray = 'Ct'
h_V8_CtvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
h_V8_CtvtkDisplay.SelectOrientationVectors = 'Ct'
h_V8_CtvtkDisplay.ScaleFactor = 0.1
h_V8_CtvtkDisplay.SelectScaleArray = 'Ct'
h_V8_CtvtkDisplay.GlyphType = 'Arrow'
h_V8_CtvtkDisplay.GlyphTableIndexArray = 'Ct'
h_V8_CtvtkDisplay.GaussianRadius = 0.005
h_V8_CtvtkDisplay.SetScaleArray = ['POINTS', 'Ct']
h_V8_CtvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
h_V8_CtvtkDisplay.OpacityArray = ['POINTS', 'Ct']
h_V8_CtvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
h_V8_CtvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
h_V8_CtvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
h_V8_CtvtkDisplay.ScalarOpacityFunction = ctPWF
h_V8_CtvtkDisplay.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
h_V8_CtvtkDisplay.ScaleTransferFunction.Points = [0.0012291321763768792, 0.0, 0.5, 0.0, 0.0019753649830818176, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
h_V8_CtvtkDisplay.OpacityTransferFunction.Points = [0.0012291321763768792, 0.0, 0.5, 0.0, 0.0019753649830818176, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
h_V8_CtvtkDisplay.PolarAxes.Translation = [0.0, -1.2, 0.0]

# show data from h_V13_Ctvtk
h_V13_CtvtkDisplay = Show(h_V13_Ctvtk, renderView1)

# trace defaults for the display properties.
h_V13_CtvtkDisplay.Representation = 'Outline'
h_V13_CtvtkDisplay.ColorArrayName = ['POINTS', 'Ct']
h_V13_CtvtkDisplay.LookupTable = ctLUT
h_V13_CtvtkDisplay.Position = [0.0, -2.4, 0.0]
h_V13_CtvtkDisplay.OSPRayScaleArray = 'Ct'
h_V13_CtvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
h_V13_CtvtkDisplay.SelectOrientationVectors = 'Ct'
h_V13_CtvtkDisplay.ScaleFactor = 0.1
h_V13_CtvtkDisplay.SelectScaleArray = 'Ct'
h_V13_CtvtkDisplay.GlyphType = 'Arrow'
h_V13_CtvtkDisplay.GlyphTableIndexArray = 'Ct'
h_V13_CtvtkDisplay.GaussianRadius = 0.005
h_V13_CtvtkDisplay.SetScaleArray = ['POINTS', 'Ct']
h_V13_CtvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
h_V13_CtvtkDisplay.OpacityArray = ['POINTS', 'Ct']
h_V13_CtvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
h_V13_CtvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
h_V13_CtvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
h_V13_CtvtkDisplay.ScalarOpacityFunction = ctPWF
h_V13_CtvtkDisplay.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
h_V13_CtvtkDisplay.ScaleTransferFunction.Points = [0.001778498524799943, 0.0, 0.5, 0.0, 0.0023644512984901667, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
h_V13_CtvtkDisplay.OpacityTransferFunction.Points = [0.001778498524799943, 0.0, 0.5, 0.0, 0.0023644512984901667, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
h_V13_CtvtkDisplay.PolarAxes.Translation = [0.0, -2.4, 0.0]

# show data from h_V18_Ctvtk
h_V18_CtvtkDisplay = Show(h_V18_Ctvtk, renderView1)

# trace defaults for the display properties.
h_V18_CtvtkDisplay.Representation = 'Outline'
h_V18_CtvtkDisplay.ColorArrayName = ['POINTS', 'Ct']
h_V18_CtvtkDisplay.LookupTable = ctLUT
h_V18_CtvtkDisplay.Position = [0.0, -3.6, 0.0]
h_V18_CtvtkDisplay.OSPRayScaleArray = 'Ct'
h_V18_CtvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
h_V18_CtvtkDisplay.SelectOrientationVectors = 'Ct'
h_V18_CtvtkDisplay.ScaleFactor = 0.1
h_V18_CtvtkDisplay.SelectScaleArray = 'Ct'
h_V18_CtvtkDisplay.GlyphType = 'Arrow'
h_V18_CtvtkDisplay.GlyphTableIndexArray = 'Ct'
h_V18_CtvtkDisplay.GaussianRadius = 0.005
h_V18_CtvtkDisplay.SetScaleArray = ['POINTS', 'Ct']
h_V18_CtvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
h_V18_CtvtkDisplay.OpacityArray = ['POINTS', 'Ct']
h_V18_CtvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
h_V18_CtvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
h_V18_CtvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
h_V18_CtvtkDisplay.ScalarOpacityFunction = ctPWF
h_V18_CtvtkDisplay.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
h_V18_CtvtkDisplay.ScaleTransferFunction.Points = [0.0019505424425005913, 0.0, 0.5, 0.0, 0.0025291552301496267, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
h_V18_CtvtkDisplay.OpacityTransferFunction.Points = [0.0019505424425005913, 0.0, 0.5, 0.0, 0.0025291552301496267, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
h_V18_CtvtkDisplay.PolarAxes.Translation = [0.0, -3.6, 0.0]

# show data from h_V2_PO2_v
h_V2_PO2_vDisplay = Show(h_V2_PO2_v, renderView1)

# trace defaults for the display properties.
h_V2_PO2_vDisplay.Representation = 'Surface'
h_V2_PO2_vDisplay.ColorArrayName = ['POINTS', '']
h_V2_PO2_vDisplay.LineWidth = 2.0
h_V2_PO2_vDisplay.RenderLinesAsTubes = 1
h_V2_PO2_vDisplay.OSPRayScaleArray = 'PO2_v'
h_V2_PO2_vDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
h_V2_PO2_vDisplay.SelectOrientationVectors = 'Cv'
h_V2_PO2_vDisplay.ScaleFactor = 0.1
h_V2_PO2_vDisplay.SelectScaleArray = 'PO2_v'
h_V2_PO2_vDisplay.GlyphType = 'Arrow'
h_V2_PO2_vDisplay.GlyphTableIndexArray = 'PO2_v'
h_V2_PO2_vDisplay.GaussianRadius = 0.005
h_V2_PO2_vDisplay.SetScaleArray = ['POINTS', 'PO2_v']
h_V2_PO2_vDisplay.ScaleTransferFunction = 'PiecewiseFunction'
h_V2_PO2_vDisplay.OpacityArray = ['POINTS', 'PO2_v']
h_V2_PO2_vDisplay.OpacityTransferFunction = 'PiecewiseFunction'
h_V2_PO2_vDisplay.DataAxesGrid = 'GridAxesRepresentation'
h_V2_PO2_vDisplay.PolarAxes = 'PolarAxesRepresentation'
h_V2_PO2_vDisplay.ScalarOpacityUnitDistance = 0.20561578771712258

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
h_V2_PO2_vDisplay.ScaleTransferFunction.Points = [25.006895884871483, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
h_V2_PO2_vDisplay.OpacityTransferFunction.Points = [25.006895884871483, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# show data from h_V8_PO2_v
h_V8_PO2_vDisplay = Show(h_V8_PO2_v, renderView1)

# trace defaults for the display properties.
h_V8_PO2_vDisplay.Representation = 'Surface'
h_V8_PO2_vDisplay.ColorArrayName = ['POINTS', '']
h_V8_PO2_vDisplay.LineWidth = 2.0
h_V8_PO2_vDisplay.RenderLinesAsTubes = 1
h_V8_PO2_vDisplay.Position = [0.0, -1.2, 0.0]
h_V8_PO2_vDisplay.OSPRayScaleArray = 'PO2_v'
h_V8_PO2_vDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
h_V8_PO2_vDisplay.SelectOrientationVectors = 'Cv'
h_V8_PO2_vDisplay.ScaleFactor = 0.1
h_V8_PO2_vDisplay.SelectScaleArray = 'PO2_v'
h_V8_PO2_vDisplay.GlyphType = 'Arrow'
h_V8_PO2_vDisplay.GlyphTableIndexArray = 'PO2_v'
h_V8_PO2_vDisplay.GaussianRadius = 0.005
h_V8_PO2_vDisplay.SetScaleArray = ['POINTS', 'PO2_v']
h_V8_PO2_vDisplay.ScaleTransferFunction = 'PiecewiseFunction'
h_V8_PO2_vDisplay.OpacityArray = ['POINTS', 'PO2_v']
h_V8_PO2_vDisplay.OpacityTransferFunction = 'PiecewiseFunction'
h_V8_PO2_vDisplay.DataAxesGrid = 'GridAxesRepresentation'
h_V8_PO2_vDisplay.PolarAxes = 'PolarAxesRepresentation'
h_V8_PO2_vDisplay.ScalarOpacityUnitDistance = 0.13111097749257294

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
h_V8_PO2_vDisplay.ScaleTransferFunction.Points = [48.47661669676502, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
h_V8_PO2_vDisplay.OpacityTransferFunction.Points = [48.47661669676502, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
h_V8_PO2_vDisplay.PolarAxes.Translation = [0.0, -1.2, 0.0]

# show data from h_V13_PO2_v
h_V13_PO2_vDisplay = Show(h_V13_PO2_v, renderView1)

# trace defaults for the display properties.
h_V13_PO2_vDisplay.Representation = 'Surface'
h_V13_PO2_vDisplay.ColorArrayName = ['POINTS', '']
h_V13_PO2_vDisplay.LineWidth = 2.0
h_V13_PO2_vDisplay.RenderLinesAsTubes = 1
h_V13_PO2_vDisplay.Position = [0.0, -2.4, 0.0]
h_V13_PO2_vDisplay.OSPRayScaleArray = 'PO2_v'
h_V13_PO2_vDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
h_V13_PO2_vDisplay.SelectOrientationVectors = 'Cv'
h_V13_PO2_vDisplay.ScaleFactor = 0.1
h_V13_PO2_vDisplay.SelectScaleArray = 'PO2_v'
h_V13_PO2_vDisplay.GlyphType = 'Arrow'
h_V13_PO2_vDisplay.GlyphTableIndexArray = 'PO2_v'
h_V13_PO2_vDisplay.GaussianRadius = 0.005
h_V13_PO2_vDisplay.SetScaleArray = ['POINTS', 'PO2_v']
h_V13_PO2_vDisplay.ScaleTransferFunction = 'PiecewiseFunction'
h_V13_PO2_vDisplay.OpacityArray = ['POINTS', 'PO2_v']
h_V13_PO2_vDisplay.OpacityTransferFunction = 'PiecewiseFunction'
h_V13_PO2_vDisplay.DataAxesGrid = 'GridAxesRepresentation'
h_V13_PO2_vDisplay.PolarAxes = 'PolarAxesRepresentation'
h_V13_PO2_vDisplay.ScalarOpacityUnitDistance = 0.11224279724452858

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
h_V13_PO2_vDisplay.ScaleTransferFunction.Points = [64.55377442762256, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
h_V13_PO2_vDisplay.OpacityTransferFunction.Points = [64.55377442762256, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
h_V13_PO2_vDisplay.PolarAxes.Translation = [0.0, -2.4, 0.0]

# show data from h_V18_PO2_v
h_V18_PO2_vDisplay = Show(h_V18_PO2_v, renderView1)

# trace defaults for the display properties.
h_V18_PO2_vDisplay.Representation = 'Surface'
h_V18_PO2_vDisplay.ColorArrayName = ['POINTS', '']
h_V18_PO2_vDisplay.LineWidth = 2.0
h_V18_PO2_vDisplay.RenderLinesAsTubes = 1
h_V18_PO2_vDisplay.Position = [0.0, -3.6, 0.0]
h_V18_PO2_vDisplay.OSPRayScaleArray = 'PO2_v'
h_V18_PO2_vDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
h_V18_PO2_vDisplay.SelectOrientationVectors = 'Cv'
h_V18_PO2_vDisplay.ScaleFactor = 0.1
h_V18_PO2_vDisplay.SelectScaleArray = 'PO2_v'
h_V18_PO2_vDisplay.GlyphType = 'Arrow'
h_V18_PO2_vDisplay.GlyphTableIndexArray = 'PO2_v'
h_V18_PO2_vDisplay.GaussianRadius = 0.005
h_V18_PO2_vDisplay.SetScaleArray = ['POINTS', 'PO2_v']
h_V18_PO2_vDisplay.ScaleTransferFunction = 'PiecewiseFunction'
h_V18_PO2_vDisplay.OpacityArray = ['POINTS', 'PO2_v']
h_V18_PO2_vDisplay.OpacityTransferFunction = 'PiecewiseFunction'
h_V18_PO2_vDisplay.DataAxesGrid = 'GridAxesRepresentation'
h_V18_PO2_vDisplay.PolarAxes = 'PolarAxesRepresentation'
h_V18_PO2_vDisplay.ScalarOpacityUnitDistance = 0.1008201900683401

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
h_V18_PO2_vDisplay.ScaleTransferFunction.Points = [73.76881937185924, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
h_V18_PO2_vDisplay.OpacityTransferFunction.Points = [73.76881937185924, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
h_V18_PO2_vDisplay.PolarAxes.Translation = [0.0, -3.6, 0.0]

# show data from t_V2_Ctvtk
t_V2_CtvtkDisplay = Show(t_V2_Ctvtk, renderView1)

# trace defaults for the display properties.
t_V2_CtvtkDisplay.Representation = 'Outline'
t_V2_CtvtkDisplay.ColorArrayName = ['POINTS', 'Ct']
t_V2_CtvtkDisplay.LookupTable = ctLUT
t_V2_CtvtkDisplay.Position = [1.2, 0.0, 0.0]
t_V2_CtvtkDisplay.OSPRayScaleArray = 'Ct'
t_V2_CtvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
t_V2_CtvtkDisplay.SelectOrientationVectors = 'Ct'
t_V2_CtvtkDisplay.ScaleFactor = 0.1
t_V2_CtvtkDisplay.SelectScaleArray = 'Ct'
t_V2_CtvtkDisplay.GlyphType = 'Arrow'
t_V2_CtvtkDisplay.GlyphTableIndexArray = 'Ct'
t_V2_CtvtkDisplay.GaussianRadius = 0.005
t_V2_CtvtkDisplay.SetScaleArray = ['POINTS', 'Ct']
t_V2_CtvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
t_V2_CtvtkDisplay.OpacityArray = ['POINTS', 'Ct']
t_V2_CtvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
t_V2_CtvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
t_V2_CtvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
t_V2_CtvtkDisplay.ScalarOpacityFunction = ctPWF
t_V2_CtvtkDisplay.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
t_V2_CtvtkDisplay.ScaleTransferFunction.Points = [-4.730674845632166e-06, 0.0, 0.5, 0.0, 0.0002281701163155958, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
t_V2_CtvtkDisplay.OpacityTransferFunction.Points = [-4.730674845632166e-06, 0.0, 0.5, 0.0, 0.0002281701163155958, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
t_V2_CtvtkDisplay.PolarAxes.Translation = [1.2, 0.0, 0.0]

# show data from t_V13_Ctvtk
t_V13_CtvtkDisplay = Show(t_V13_Ctvtk, renderView1)

# trace defaults for the display properties.
t_V13_CtvtkDisplay.Representation = 'Outline'
t_V13_CtvtkDisplay.ColorArrayName = ['POINTS', 'Ct']
t_V13_CtvtkDisplay.LookupTable = ctLUT
t_V13_CtvtkDisplay.Position = [1.2, -2.4, 0.0]
t_V13_CtvtkDisplay.OSPRayScaleArray = 'Ct'
t_V13_CtvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
t_V13_CtvtkDisplay.SelectOrientationVectors = 'Ct'
t_V13_CtvtkDisplay.ScaleFactor = 0.1
t_V13_CtvtkDisplay.SelectScaleArray = 'Ct'
t_V13_CtvtkDisplay.GlyphType = 'Arrow'
t_V13_CtvtkDisplay.GlyphTableIndexArray = 'Ct'
t_V13_CtvtkDisplay.GaussianRadius = 0.005
t_V13_CtvtkDisplay.SetScaleArray = ['POINTS', 'Ct']
t_V13_CtvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
t_V13_CtvtkDisplay.OpacityArray = ['POINTS', 'Ct']
t_V13_CtvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
t_V13_CtvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
t_V13_CtvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
t_V13_CtvtkDisplay.ScalarOpacityFunction = ctPWF
t_V13_CtvtkDisplay.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
t_V13_CtvtkDisplay.ScaleTransferFunction.Points = [-5.020732487537316e-07, 0.0, 0.5, 0.0, 0.0011694994755089283, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
t_V13_CtvtkDisplay.OpacityTransferFunction.Points = [-5.020732487537316e-07, 0.0, 0.5, 0.0, 0.0011694994755089283, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
t_V13_CtvtkDisplay.PolarAxes.Translation = [1.2, -2.4, 0.0]

# show data from t_V8_Ctvtk
t_V8_CtvtkDisplay = Show(t_V8_Ctvtk, renderView1)

# trace defaults for the display properties.
t_V8_CtvtkDisplay.Representation = 'Outline'
t_V8_CtvtkDisplay.ColorArrayName = ['POINTS', 'Ct']
t_V8_CtvtkDisplay.LookupTable = ctLUT
t_V8_CtvtkDisplay.Position = [1.2, -1.2, 0.0]
t_V8_CtvtkDisplay.OSPRayScaleArray = 'Ct'
t_V8_CtvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
t_V8_CtvtkDisplay.SelectOrientationVectors = 'Ct'
t_V8_CtvtkDisplay.ScaleFactor = 0.1
t_V8_CtvtkDisplay.SelectScaleArray = 'Ct'
t_V8_CtvtkDisplay.GlyphType = 'Arrow'
t_V8_CtvtkDisplay.GlyphTableIndexArray = 'Ct'
t_V8_CtvtkDisplay.GaussianRadius = 0.005
t_V8_CtvtkDisplay.SetScaleArray = ['POINTS', 'Ct']
t_V8_CtvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
t_V8_CtvtkDisplay.OpacityArray = ['POINTS', 'Ct']
t_V8_CtvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
t_V8_CtvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
t_V8_CtvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
t_V8_CtvtkDisplay.ScalarOpacityFunction = ctPWF
t_V8_CtvtkDisplay.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
t_V8_CtvtkDisplay.ScaleTransferFunction.Points = [-3.0944102036301047e-06, 0.0, 0.5, 0.0, 0.0005438623484224081, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
t_V8_CtvtkDisplay.OpacityTransferFunction.Points = [-3.0944102036301047e-06, 0.0, 0.5, 0.0, 0.0005438623484224081, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
t_V8_CtvtkDisplay.PolarAxes.Translation = [1.2, -1.2, 0.0]

# show data from t_V18_Ctvtk
t_V18_CtvtkDisplay = Show(t_V18_Ctvtk, renderView1)

# trace defaults for the display properties.
t_V18_CtvtkDisplay.Representation = 'Outline'
t_V18_CtvtkDisplay.ColorArrayName = ['POINTS', 'Ct']
t_V18_CtvtkDisplay.LookupTable = ctLUT
t_V18_CtvtkDisplay.Position = [1.2, -3.6, 0.0]
t_V18_CtvtkDisplay.OSPRayScaleArray = 'Ct'
t_V18_CtvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
t_V18_CtvtkDisplay.SelectOrientationVectors = 'Ct'
t_V18_CtvtkDisplay.ScaleFactor = 0.1
t_V18_CtvtkDisplay.SelectScaleArray = 'Ct'
t_V18_CtvtkDisplay.GlyphType = 'Arrow'
t_V18_CtvtkDisplay.GlyphTableIndexArray = 'Ct'
t_V18_CtvtkDisplay.GaussianRadius = 0.005
t_V18_CtvtkDisplay.SetScaleArray = ['POINTS', 'Ct']
t_V18_CtvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
t_V18_CtvtkDisplay.OpacityArray = ['POINTS', 'Ct']
t_V18_CtvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
t_V18_CtvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
t_V18_CtvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
t_V18_CtvtkDisplay.ScalarOpacityFunction = ctPWF
t_V18_CtvtkDisplay.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
t_V18_CtvtkDisplay.ScaleTransferFunction.Points = [7.51473635318689e-05, 0.0, 0.5, 0.0, 0.0016079472843557596, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
t_V18_CtvtkDisplay.OpacityTransferFunction.Points = [7.51473635318689e-05, 0.0, 0.5, 0.0, 0.0016079472843557596, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
t_V18_CtvtkDisplay.PolarAxes.Translation = [1.2, -3.6, 0.0]

# show data from t_V2_PO2_v
t_V2_PO2_vDisplay = Show(t_V2_PO2_v, renderView1)

# trace defaults for the display properties.
t_V2_PO2_vDisplay.Representation = 'Surface'
t_V2_PO2_vDisplay.ColorArrayName = ['POINTS', '']
t_V2_PO2_vDisplay.LineWidth = 2.0
t_V2_PO2_vDisplay.RenderLinesAsTubes = 1
t_V2_PO2_vDisplay.Position = [1.2, 0.0, 0.0]
t_V2_PO2_vDisplay.OSPRayScaleArray = 'PO2_v'
t_V2_PO2_vDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
t_V2_PO2_vDisplay.SelectOrientationVectors = 'Cv'
t_V2_PO2_vDisplay.ScaleFactor = 0.1
t_V2_PO2_vDisplay.SelectScaleArray = 'PO2_v'
t_V2_PO2_vDisplay.GlyphType = 'Arrow'
t_V2_PO2_vDisplay.GlyphTableIndexArray = 'PO2_v'
t_V2_PO2_vDisplay.GaussianRadius = 0.005
t_V2_PO2_vDisplay.SetScaleArray = ['POINTS', 'PO2_v']
t_V2_PO2_vDisplay.ScaleTransferFunction = 'PiecewiseFunction'
t_V2_PO2_vDisplay.OpacityArray = ['POINTS', 'PO2_v']
t_V2_PO2_vDisplay.OpacityTransferFunction = 'PiecewiseFunction'
t_V2_PO2_vDisplay.DataAxesGrid = 'GridAxesRepresentation'
t_V2_PO2_vDisplay.PolarAxes = 'PolarAxesRepresentation'
t_V2_PO2_vDisplay.ScalarOpacityUnitDistance = 0.20561578771712258

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
t_V2_PO2_vDisplay.ScaleTransferFunction.Points = [23.718127825607855, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
t_V2_PO2_vDisplay.OpacityTransferFunction.Points = [23.718127825607855, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
t_V2_PO2_vDisplay.PolarAxes.Translation = [1.2, 0.0, 0.0]

# show data from t_V8_PO2_v
t_V8_PO2_vDisplay = Show(t_V8_PO2_v, renderView1)

# trace defaults for the display properties.
t_V8_PO2_vDisplay.Representation = 'Surface'
t_V8_PO2_vDisplay.ColorArrayName = ['POINTS', '']
t_V8_PO2_vDisplay.LineWidth = 2.0
t_V8_PO2_vDisplay.RenderLinesAsTubes = 1
t_V8_PO2_vDisplay.Position = [1.2, -1.2, 0.0]
t_V8_PO2_vDisplay.OSPRayScaleArray = 'PO2_v'
t_V8_PO2_vDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
t_V8_PO2_vDisplay.SelectOrientationVectors = 'Cv'
t_V8_PO2_vDisplay.ScaleFactor = 0.1
t_V8_PO2_vDisplay.SelectScaleArray = 'PO2_v'
t_V8_PO2_vDisplay.GlyphType = 'Arrow'
t_V8_PO2_vDisplay.GlyphTableIndexArray = 'PO2_v'
t_V8_PO2_vDisplay.GaussianRadius = 0.005
t_V8_PO2_vDisplay.SetScaleArray = ['POINTS', 'PO2_v']
t_V8_PO2_vDisplay.ScaleTransferFunction = 'PiecewiseFunction'
t_V8_PO2_vDisplay.OpacityArray = ['POINTS', 'PO2_v']
t_V8_PO2_vDisplay.OpacityTransferFunction = 'PiecewiseFunction'
t_V8_PO2_vDisplay.DataAxesGrid = 'GridAxesRepresentation'
t_V8_PO2_vDisplay.PolarAxes = 'PolarAxesRepresentation'
t_V8_PO2_vDisplay.ScalarOpacityUnitDistance = 0.13111097749257294

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
t_V8_PO2_vDisplay.ScaleTransferFunction.Points = [13.660525049393375, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
t_V8_PO2_vDisplay.OpacityTransferFunction.Points = [13.660525049393375, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
t_V8_PO2_vDisplay.PolarAxes.Translation = [1.2, -1.2, 0.0]

# show data from t_V13_PO2_v
t_V13_PO2_vDisplay = Show(t_V13_PO2_v, renderView1)

# trace defaults for the display properties.
t_V13_PO2_vDisplay.Representation = 'Surface'
t_V13_PO2_vDisplay.ColorArrayName = ['POINTS', '']
t_V13_PO2_vDisplay.LineWidth = 2.0
t_V13_PO2_vDisplay.RenderLinesAsTubes = 1
t_V13_PO2_vDisplay.Position = [1.2, -2.4, 0.0]
t_V13_PO2_vDisplay.OSPRayScaleArray = 'PO2_v'
t_V13_PO2_vDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
t_V13_PO2_vDisplay.SelectOrientationVectors = 'Cv'
t_V13_PO2_vDisplay.ScaleFactor = 0.1
t_V13_PO2_vDisplay.SelectScaleArray = 'PO2_v'
t_V13_PO2_vDisplay.GlyphType = 'Arrow'
t_V13_PO2_vDisplay.GlyphTableIndexArray = 'PO2_v'
t_V13_PO2_vDisplay.GaussianRadius = 0.005
t_V13_PO2_vDisplay.SetScaleArray = ['POINTS', 'PO2_v']
t_V13_PO2_vDisplay.ScaleTransferFunction = 'PiecewiseFunction'
t_V13_PO2_vDisplay.OpacityArray = ['POINTS', 'PO2_v']
t_V13_PO2_vDisplay.OpacityTransferFunction = 'PiecewiseFunction'
t_V13_PO2_vDisplay.DataAxesGrid = 'GridAxesRepresentation'
t_V13_PO2_vDisplay.PolarAxes = 'PolarAxesRepresentation'
t_V13_PO2_vDisplay.ScalarOpacityUnitDistance = 0.11224279724452858

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
t_V13_PO2_vDisplay.ScaleTransferFunction.Points = [22.098988605042297, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
t_V13_PO2_vDisplay.OpacityTransferFunction.Points = [22.098988605042297, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
t_V13_PO2_vDisplay.PolarAxes.Translation = [1.2, -2.4, 0.0]

# show data from t_V18_PO2_v
t_V18_PO2_vDisplay = Show(t_V18_PO2_v, renderView1)

# trace defaults for the display properties.
t_V18_PO2_vDisplay.Representation = 'Surface'
t_V18_PO2_vDisplay.ColorArrayName = ['POINTS', '']
t_V18_PO2_vDisplay.LineWidth = 2.0
t_V18_PO2_vDisplay.RenderLinesAsTubes = 1
t_V18_PO2_vDisplay.Position = [1.2, -3.6, 0.0]
t_V18_PO2_vDisplay.OSPRayScaleArray = 'PO2_v'
t_V18_PO2_vDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
t_V18_PO2_vDisplay.SelectOrientationVectors = 'Cv'
t_V18_PO2_vDisplay.ScaleFactor = 0.1
t_V18_PO2_vDisplay.SelectScaleArray = 'PO2_v'
t_V18_PO2_vDisplay.GlyphType = 'Arrow'
t_V18_PO2_vDisplay.GlyphTableIndexArray = 'PO2_v'
t_V18_PO2_vDisplay.GaussianRadius = 0.005
t_V18_PO2_vDisplay.SetScaleArray = ['POINTS', 'PO2_v']
t_V18_PO2_vDisplay.ScaleTransferFunction = 'PiecewiseFunction'
t_V18_PO2_vDisplay.OpacityArray = ['POINTS', 'PO2_v']
t_V18_PO2_vDisplay.OpacityTransferFunction = 'PiecewiseFunction'
t_V18_PO2_vDisplay.DataAxesGrid = 'GridAxesRepresentation'
t_V18_PO2_vDisplay.PolarAxes = 'PolarAxesRepresentation'
t_V18_PO2_vDisplay.ScalarOpacityUnitDistance = 0.1008201900683401

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
t_V18_PO2_vDisplay.ScaleTransferFunction.Points = [34.96970748528838, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
t_V18_PO2_vDisplay.OpacityTransferFunction.Points = [34.96970748528838, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
t_V18_PO2_vDisplay.PolarAxes.Translation = [1.2, -3.6, 0.0]

# show data from bETA_8
bETA_8Display = Show(bETA_8, renderView1)

# get color transfer function/color map for 'beta'
betaLUT = GetColorTransferFunction('beta')
betaLUT.RGBPoints = [0.01089074460406303, 0.231373, 0.298039, 0.752941, 0.03606396676330412, 0.865003, 0.865003, 0.865003, 0.06123718892254548, 0.705882, 0.0156863, 0.14902]
betaLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'beta'
betaPWF = GetOpacityTransferFunction('beta')
betaPWF.Points = [0.01089074460406303, 0.0, 0.5, 0.0, 0.06123718892254548, 1.0, 0.5, 0.0]
betaPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
bETA_8Display.Representation = 'Surface'
bETA_8Display.ColorArrayName = ['POINTS', 'beta']
bETA_8Display.LookupTable = betaLUT
bETA_8Display.OSPRayScaleArray = 'beta'
bETA_8Display.OSPRayScaleFunction = 'PiecewiseFunction'
bETA_8Display.SelectOrientationVectors = 'Ct'
bETA_8Display.ScaleFactor = 0.1
bETA_8Display.SelectScaleArray = 'beta'
bETA_8Display.GlyphType = 'Arrow'
bETA_8Display.GlyphTableIndexArray = 'beta'
bETA_8Display.GaussianRadius = 0.005
bETA_8Display.SetScaleArray = ['POINTS', 'beta']
bETA_8Display.ScaleTransferFunction = 'PiecewiseFunction'
bETA_8Display.OpacityArray = ['POINTS', 'beta']
bETA_8Display.OpacityTransferFunction = 'PiecewiseFunction'
bETA_8Display.DataAxesGrid = 'GridAxesRepresentation'
bETA_8Display.PolarAxes = 'PolarAxesRepresentation'
bETA_8Display.ScalarOpacityFunction = betaPWF
bETA_8Display.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
bETA_8Display.ScaleTransferFunction.Points = [0.14570732350093987, 0.0, 0.5, 0.0, 0.1487986635253606, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
bETA_8Display.OpacityTransferFunction.Points = [0.14570732350093987, 0.0, 0.5, 0.0, 0.1487986635253606, 1.0, 0.5, 0.0]

# show data from clip9
clip9Display = Show(clip9, renderView1)

# get color transfer function/color map for 'Sf'
sfLUT = GetColorTransferFunction('Sf')
sfLUT.RGBPoints = [8.3458318405707435e-50, 0.0, 0.0, 1.0, 0.9048033875095829, 1.0, 0.0, 0.0]
sfLUT.ColorSpace = 'HSV'
sfLUT.NanColor = [0.498039215686, 0.498039215686, 0.498039215686]
sfLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Sf'
sfPWF = GetOpacityTransferFunction('Sf')
sfPWF.Points = [8.3458318405707435e-50, 0.0, 0.5, 0.0, 0.9048033875095829, 1.0, 0.5, 0.0]
sfPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
clip9Display.Representation = 'Surface'
clip9Display.ColorArrayName = ['POINTS', 'Sf']
clip9Display.LookupTable = sfLUT
clip9Display.OSPRayScaleArray = 'Sf'
clip9Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip9Display.SelectOrientationVectors = 'Ct'
clip9Display.ScaleFactor = 0.1
clip9Display.SelectScaleArray = 'Sf'
clip9Display.GlyphType = 'Arrow'
clip9Display.GlyphTableIndexArray = 'Sf'
clip9Display.GaussianRadius = 0.005
clip9Display.SetScaleArray = ['POINTS', 'Sf']
clip9Display.ScaleTransferFunction = 'PiecewiseFunction'
clip9Display.OpacityArray = ['POINTS', 'Sf']
clip9Display.OpacityTransferFunction = 'PiecewiseFunction'
clip9Display.DataAxesGrid = 'GridAxesRepresentation'
clip9Display.PolarAxes = 'PolarAxesRepresentation'
clip9Display.ScalarOpacityFunction = sfPWF
clip9Display.ScalarOpacityUnitDistance = 0.10263591818176866

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip9Display.ScaleTransferFunction.Points = [0.43790864199098894, 0.0, 0.5, 0.0, 0.8614589284771672, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip9Display.OpacityTransferFunction.Points = [0.43790864199098894, 0.0, 0.5, 0.0, 0.8614589284771672, 1.0, 0.5, 0.0]

# show data from clip10
clip10Display = Show(clip10, renderView1)

# trace defaults for the display properties.
clip10Display.Representation = 'Surface'
clip10Display.ColorArrayName = ['POINTS', 'Sf']
clip10Display.LookupTable = sfLUT
clip10Display.Position = [0.0, -1.2, 0.0]
clip10Display.OSPRayScaleArray = 'Sf'
clip10Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip10Display.SelectOrientationVectors = 'Ct'
clip10Display.ScaleFactor = 0.1
clip10Display.SelectScaleArray = 'Sf'
clip10Display.GlyphType = 'Arrow'
clip10Display.GlyphTableIndexArray = 'Sf'
clip10Display.GaussianRadius = 0.005
clip10Display.SetScaleArray = ['POINTS', 'Sf']
clip10Display.ScaleTransferFunction = 'PiecewiseFunction'
clip10Display.OpacityArray = ['POINTS', 'Sf']
clip10Display.OpacityTransferFunction = 'PiecewiseFunction'
clip10Display.DataAxesGrid = 'GridAxesRepresentation'
clip10Display.PolarAxes = 'PolarAxesRepresentation'
clip10Display.ScalarOpacityFunction = sfPWF
clip10Display.ScalarOpacityUnitDistance = 0.10263591818176866

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip10Display.ScaleTransferFunction.Points = [0.3615784400971896, 0.0, 0.5, 0.0, 0.3707904092548618, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip10Display.OpacityTransferFunction.Points = [0.3615784400971896, 0.0, 0.5, 0.0, 0.3707904092548618, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip10Display.PolarAxes.Translation = [0.0, -1.2, 0.0]

# show data from clip11
clip11Display = Show(clip11, renderView1)

# trace defaults for the display properties.
clip11Display.Representation = 'Surface'
clip11Display.ColorArrayName = ['POINTS', 'Sf']
clip11Display.LookupTable = sfLUT
clip11Display.Position = [0.0, -2.4, 0.0]
clip11Display.OSPRayScaleArray = 'Sf'
clip11Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip11Display.SelectOrientationVectors = 'Ct'
clip11Display.ScaleFactor = 0.1
clip11Display.SelectScaleArray = 'Sf'
clip11Display.GlyphType = 'Arrow'
clip11Display.GlyphTableIndexArray = 'Sf'
clip11Display.GaussianRadius = 0.005
clip11Display.SetScaleArray = ['POINTS', 'Sf']
clip11Display.ScaleTransferFunction = 'PiecewiseFunction'
clip11Display.OpacityArray = ['POINTS', 'Sf']
clip11Display.OpacityTransferFunction = 'PiecewiseFunction'
clip11Display.DataAxesGrid = 'GridAxesRepresentation'
clip11Display.PolarAxes = 'PolarAxesRepresentation'
clip11Display.ScalarOpacityFunction = sfPWF
clip11Display.ScalarOpacityUnitDistance = 0.10263591818176866

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip11Display.ScaleTransferFunction.Points = [0.35337923923673226, 0.0, 0.5, 0.0, 0.3573571803566705, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip11Display.OpacityTransferFunction.Points = [0.35337923923673226, 0.0, 0.5, 0.0, 0.3573571803566705, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip11Display.PolarAxes.Translation = [0.0, -2.4, 0.0]

# show data from clip12
clip12Display = Show(clip12, renderView1)

# trace defaults for the display properties.
clip12Display.Representation = 'Surface'
clip12Display.ColorArrayName = ['POINTS', 'Sf']
clip12Display.LookupTable = sfLUT
clip12Display.Position = [0.0, -3.6, 0.0]
clip12Display.OSPRayScaleArray = 'Sf'
clip12Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip12Display.SelectOrientationVectors = 'Ct'
clip12Display.ScaleFactor = 0.1
clip12Display.SelectScaleArray = 'Sf'
clip12Display.GlyphType = 'Arrow'
clip12Display.GlyphTableIndexArray = 'Sf'
clip12Display.GaussianRadius = 0.005
clip12Display.SetScaleArray = ['POINTS', 'Sf']
clip12Display.ScaleTransferFunction = 'PiecewiseFunction'
clip12Display.OpacityArray = ['POINTS', 'Sf']
clip12Display.OpacityTransferFunction = 'PiecewiseFunction'
clip12Display.DataAxesGrid = 'GridAxesRepresentation'
clip12Display.PolarAxes = 'PolarAxesRepresentation'
clip12Display.ScalarOpacityFunction = sfPWF
clip12Display.ScalarOpacityUnitDistance = 0.10263591818176866

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip12Display.ScaleTransferFunction.Points = [0.3561836248692617, 0.0, 0.5, 0.0, 0.35963832842885984, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip12Display.OpacityTransferFunction.Points = [0.3561836248692617, 0.0, 0.5, 0.0, 0.35963832842885984, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip12Display.PolarAxes.Translation = [0.0, -3.6, 0.0]

# show data from clip13
clip13Display = Show(clip13, renderView1)

# trace defaults for the display properties.
clip13Display.Representation = 'Surface'
clip13Display.ColorArrayName = ['POINTS', 'Sf']
clip13Display.LookupTable = sfLUT
clip13Display.Position = [1.2, 0.0, 0.0]
clip13Display.OSPRayScaleArray = 'Sf'
clip13Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip13Display.SelectOrientationVectors = 'Ct'
clip13Display.ScaleFactor = 0.1
clip13Display.SelectScaleArray = 'Sf'
clip13Display.GlyphType = 'Arrow'
clip13Display.GlyphTableIndexArray = 'Sf'
clip13Display.GaussianRadius = 0.005
clip13Display.SetScaleArray = ['POINTS', 'Sf']
clip13Display.ScaleTransferFunction = 'PiecewiseFunction'
clip13Display.OpacityArray = ['POINTS', 'Sf']
clip13Display.OpacityTransferFunction = 'PiecewiseFunction'
clip13Display.DataAxesGrid = 'GridAxesRepresentation'
clip13Display.PolarAxes = 'PolarAxesRepresentation'
clip13Display.ScalarOpacityFunction = sfPWF
clip13Display.ScalarOpacityUnitDistance = 0.10263591818176866

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip13Display.ScaleTransferFunction.Points = [0.5098677047985775, 0.0, 0.5, 0.0, 0.9084716761090437, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip13Display.OpacityTransferFunction.Points = [0.5098677047985775, 0.0, 0.5, 0.0, 0.9084716761090437, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip13Display.PolarAxes.Translation = [1.2, 0.0, 0.0]

# show data from clip14
clip14Display = Show(clip14, renderView1)

# trace defaults for the display properties.
clip14Display.Representation = 'Surface'
clip14Display.ColorArrayName = ['POINTS', 'Sf']
clip14Display.LookupTable = sfLUT
clip14Display.Position = [1.2, -1.2, 0.0]
clip14Display.OSPRayScaleArray = 'Sf'
clip14Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip14Display.SelectOrientationVectors = 'Ct'
clip14Display.ScaleFactor = 0.1
clip14Display.SelectScaleArray = 'Sf'
clip14Display.GlyphType = 'Arrow'
clip14Display.GlyphTableIndexArray = 'Sf'
clip14Display.GaussianRadius = 0.005
clip14Display.SetScaleArray = ['POINTS', 'Sf']
clip14Display.ScaleTransferFunction = 'PiecewiseFunction'
clip14Display.OpacityArray = ['POINTS', 'Sf']
clip14Display.OpacityTransferFunction = 'PiecewiseFunction'
clip14Display.DataAxesGrid = 'GridAxesRepresentation'
clip14Display.PolarAxes = 'PolarAxesRepresentation'
clip14Display.ScalarOpacityFunction = sfPWF
clip14Display.ScalarOpacityUnitDistance = 0.10263591818176866

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip14Display.ScaleTransferFunction.Points = [0.4209625873373199, 0.0, 0.5, 0.0, 0.9036257445937314, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip14Display.OpacityTransferFunction.Points = [0.4209625873373199, 0.0, 0.5, 0.0, 0.9036257445937314, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip14Display.PolarAxes.Translation = [1.2, -1.2, 0.0]

# show data from clip15
clip15Display = Show(clip15, renderView1)

# trace defaults for the display properties.
clip15Display.Representation = 'Surface'
clip15Display.ColorArrayName = ['POINTS', 'Sf']
clip15Display.LookupTable = sfLUT
clip15Display.Position = [1.2, -2.4, 0.0]
clip15Display.OSPRayScaleArray = 'Sf'
clip15Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip15Display.SelectOrientationVectors = 'Ct'
clip15Display.ScaleFactor = 0.1
clip15Display.SelectScaleArray = 'Sf'
clip15Display.GlyphType = 'Arrow'
clip15Display.GlyphTableIndexArray = 'Sf'
clip15Display.GaussianRadius = 0.005
clip15Display.SetScaleArray = ['POINTS', 'Sf']
clip15Display.ScaleTransferFunction = 'PiecewiseFunction'
clip15Display.OpacityArray = ['POINTS', 'Sf']
clip15Display.OpacityTransferFunction = 'PiecewiseFunction'
clip15Display.DataAxesGrid = 'GridAxesRepresentation'
clip15Display.PolarAxes = 'PolarAxesRepresentation'
clip15Display.ScalarOpacityFunction = sfPWF
clip15Display.ScalarOpacityUnitDistance = 0.10263591818176866

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip15Display.ScaleTransferFunction.Points = [0.37808927499608425, 0.0, 0.5, 0.0, 0.5742005181144383, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip15Display.OpacityTransferFunction.Points = [0.37808927499608425, 0.0, 0.5, 0.0, 0.5742005181144383, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip15Display.PolarAxes.Translation = [1.2, -2.4, 0.0]

# show data from clip16
clip16Display = Show(clip16, renderView1)

# trace defaults for the display properties.
clip16Display.Representation = 'Surface'
clip16Display.ColorArrayName = ['POINTS', 'Sf']
clip16Display.LookupTable = sfLUT
clip16Display.Position = [1.2, -3.6, 0.0]
clip16Display.OSPRayScaleArray = 'Sf'
clip16Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip16Display.SelectOrientationVectors = 'Ct'
clip16Display.ScaleFactor = 0.1
clip16Display.SelectScaleArray = 'Sf'
clip16Display.GlyphType = 'Arrow'
clip16Display.GlyphTableIndexArray = 'Sf'
clip16Display.GaussianRadius = 0.005
clip16Display.SetScaleArray = ['POINTS', 'Sf']
clip16Display.ScaleTransferFunction = 'PiecewiseFunction'
clip16Display.OpacityArray = ['POINTS', 'Sf']
clip16Display.OpacityTransferFunction = 'PiecewiseFunction'
clip16Display.DataAxesGrid = 'GridAxesRepresentation'
clip16Display.PolarAxes = 'PolarAxesRepresentation'
clip16Display.ScalarOpacityFunction = sfPWF
clip16Display.ScalarOpacityUnitDistance = 0.10263591818176866

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip16Display.ScaleTransferFunction.Points = [0.3671317161151878, 0.0, 0.5, 0.0, 0.4193576082450725, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip16Display.OpacityTransferFunction.Points = [0.3671317161151878, 0.0, 0.5, 0.0, 0.4193576082450725, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip16Display.PolarAxes.Translation = [1.2, -3.6, 0.0]

# show data from t_V36_Ctvtk
t_V36_CtvtkDisplay = Show(t_V36_Ctvtk, renderView1)

# trace defaults for the display properties.
t_V36_CtvtkDisplay.Representation = 'Outline'
t_V36_CtvtkDisplay.ColorArrayName = ['POINTS', 'Ct']
t_V36_CtvtkDisplay.LookupTable = ctLUT
t_V36_CtvtkDisplay.Position = [1.2, -4.8, 0.0]
t_V36_CtvtkDisplay.OSPRayScaleArray = 'Ct'
t_V36_CtvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
t_V36_CtvtkDisplay.SelectOrientationVectors = 'Ct'
t_V36_CtvtkDisplay.ScaleFactor = 0.1
t_V36_CtvtkDisplay.SelectScaleArray = 'Ct'
t_V36_CtvtkDisplay.GlyphType = 'Arrow'
t_V36_CtvtkDisplay.GlyphTableIndexArray = 'Ct'
t_V36_CtvtkDisplay.GaussianRadius = 0.005
t_V36_CtvtkDisplay.SetScaleArray = ['POINTS', 'Ct']
t_V36_CtvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
t_V36_CtvtkDisplay.OpacityArray = ['POINTS', 'Ct']
t_V36_CtvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
t_V36_CtvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
t_V36_CtvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
t_V36_CtvtkDisplay.ScalarOpacityFunction = ctPWF
t_V36_CtvtkDisplay.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
t_V36_CtvtkDisplay.ScaleTransferFunction.Points = [0.0010421141050755978, 0.0, 0.5, 0.0, 0.0022703739814460278, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
t_V36_CtvtkDisplay.OpacityTransferFunction.Points = [0.0010421141050755978, 0.0, 0.5, 0.0, 0.0022703739814460278, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
t_V36_CtvtkDisplay.PolarAxes.Translation = [1.2, -4.8, 0.0]

# show data from t_V36_PO2_v
t_V36_PO2_vDisplay = Show(t_V36_PO2_v, renderView1)

# trace defaults for the display properties.
t_V36_PO2_vDisplay.Representation = 'Surface'
t_V36_PO2_vDisplay.ColorArrayName = ['POINTS', '']
t_V36_PO2_vDisplay.LineWidth = 2.0
t_V36_PO2_vDisplay.RenderLinesAsTubes = 1
t_V36_PO2_vDisplay.Position = [1.2, -4.8, 0.0]
t_V36_PO2_vDisplay.OSPRayScaleArray = 'PO2_v'
t_V36_PO2_vDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
t_V36_PO2_vDisplay.SelectOrientationVectors = 'Cv'
t_V36_PO2_vDisplay.ScaleFactor = 0.1
t_V36_PO2_vDisplay.SelectScaleArray = 'PO2_v'
t_V36_PO2_vDisplay.GlyphType = 'Arrow'
t_V36_PO2_vDisplay.GlyphTableIndexArray = 'PO2_v'
t_V36_PO2_vDisplay.GaussianRadius = 0.005
t_V36_PO2_vDisplay.SetScaleArray = ['POINTS', 'PO2_v']
t_V36_PO2_vDisplay.ScaleTransferFunction = 'PiecewiseFunction'
t_V36_PO2_vDisplay.OpacityArray = ['POINTS', 'PO2_v']
t_V36_PO2_vDisplay.OpacityTransferFunction = 'PiecewiseFunction'
t_V36_PO2_vDisplay.DataAxesGrid = 'GridAxesRepresentation'
t_V36_PO2_vDisplay.PolarAxes = 'PolarAxesRepresentation'
t_V36_PO2_vDisplay.ScalarOpacityUnitDistance = 0.08055650894072297

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
t_V36_PO2_vDisplay.ScaleTransferFunction.Points = [46.94825038313866, 0.0, 0.5, 0.0, 102.61820474018653, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
t_V36_PO2_vDisplay.OpacityTransferFunction.Points = [46.94825038313866, 0.0, 0.5, 0.0, 102.61820474018653, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
t_V36_PO2_vDisplay.PolarAxes.Translation = [1.2, -4.8, 0.0]

# show data from clip17
clip17Display = Show(clip17, renderView1)

# trace defaults for the display properties.
clip17Display.Representation = 'Surface'
clip17Display.ColorArrayName = ['POINTS', 'Sf']
clip17Display.LookupTable = sfLUT
clip17Display.Position = [1.2, -4.8, 0.0]
clip17Display.OSPRayScaleArray = 'Sf'
clip17Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip17Display.SelectOrientationVectors = 'Ct'
clip17Display.ScaleFactor = 0.1
clip17Display.SelectScaleArray = 'Sf'
clip17Display.GlyphType = 'Arrow'
clip17Display.GlyphTableIndexArray = 'Sf'
clip17Display.GaussianRadius = 0.005
clip17Display.SetScaleArray = ['POINTS', 'Sf']
clip17Display.ScaleTransferFunction = 'PiecewiseFunction'
clip17Display.OpacityArray = ['POINTS', 'Sf']
clip17Display.OpacityTransferFunction = 'PiecewiseFunction'
clip17Display.DataAxesGrid = 'GridAxesRepresentation'
clip17Display.PolarAxes = 'PolarAxesRepresentation'
clip17Display.ScalarOpacityFunction = sfPWF
clip17Display.ScalarOpacityUnitDistance = 0.10263591818176866

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip17Display.ScaleTransferFunction.Points = [0.35838569071609133, 0.0, 0.5, 0.0, 0.36913977415231325, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip17Display.OpacityTransferFunction.Points = [0.35838569071609133, 0.0, 0.5, 0.0, 0.36913977415231325, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip17Display.PolarAxes.Translation = [1.2, -4.8, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for betaLUT in view renderView1
betaLUTColorBar = GetScalarBar(betaLUT, renderView1)
betaLUTColorBar.WindowLocation = 'UpperRightCorner'
betaLUTColorBar.Title = 'beta'
betaLUTColorBar.ComponentTitle = ''
betaLUTColorBar.HorizontalTitle = 1
betaLUTColorBar.ScalarBarLength = 0.32999999999999996

# set color bar visibility
betaLUTColorBar.Visibility = 1

# get color legend/bar for sfLUT in view renderView1
sfLUTColorBar = GetScalarBar(sfLUT, renderView1)
sfLUTColorBar.Orientation = 'Horizontal'
sfLUTColorBar.WindowLocation = 'AnyLocation'
sfLUTColorBar.Position = [0.07065459610027855, 0.11636828644501282]
sfLUTColorBar.Title = 'Sf (-)'
sfLUTColorBar.ComponentTitle = ''
sfLUTColorBar.HorizontalTitle = 1
sfLUTColorBar.UseCustomLabels = 1
sfLUTColorBar.CustomLabels = [0.25, 0.38, 0.51, 0.64, 0.77, 0.9]
sfLUTColorBar.AddRangeLabels = 0
sfLUTColorBar.ScalarBarThickness = 20
sfLUTColorBar.ScalarBarLength = 0.2000000000000001

# set color bar visibility
sfLUTColorBar.Visibility = 1

# show color legend
bETA_8Display.SetScalarBarVisibility(renderView1, True)

# show color legend
clip9Display.SetScalarBarVisibility(renderView1, True)

# show color legend
clip10Display.SetScalarBarVisibility(renderView1, True)

# show color legend
clip11Display.SetScalarBarVisibility(renderView1, True)

# show color legend
clip12Display.SetScalarBarVisibility(renderView1, True)

# show color legend
clip13Display.SetScalarBarVisibility(renderView1, True)

# show color legend
clip14Display.SetScalarBarVisibility(renderView1, True)

# show color legend
clip15Display.SetScalarBarVisibility(renderView1, True)

# show color legend
clip16Display.SetScalarBarVisibility(renderView1, True)

# show color legend
clip17Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'spreadSheetView1'
# ----------------------------------------------------------------

# show data from calculator10
calculator10Display = Show(calculator10, spreadSheetView1)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
# ----------------------------------------------------------------