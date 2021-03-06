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
renderView1.ViewSize = [1563, 782]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 0
renderView1.CenterOfRotation = [0.5, 0.5, 0.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [2.037561710312712, 0.6020307575654154, 3.5620389063922078]
renderView1.CameraFocalPoint = [1.2220311386577576, 0.47965438683393297, 0.31918733203831723]
renderView1.CameraViewUp = [-0.014137916063968765, 0.999316510058891, -0.03415599585851173]
renderView1.CameraParallelScale = 0.8660254037844386
renderView1.Background = [0.32, 0.34, 0.43]
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
t_Ctvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\v18\\Ct_3d1d.vtk'])

# create a new 'Calculator'
pO2_t = Calculator(Input=t_Ctvtk)
pO2_t.ResultArrayName = 'PO2_t'
pO2_t.Function = 'Ct/0.0000389'

# create a new 'Calculator'
a1 = Calculator(Input=pO2_t)
a1.ResultArrayName = 'a1'
a1.Function = '0.22'

# create a new 'Calculator'
a2 = Calculator(Input=a1)
a2.ResultArrayName = 'a2'
a2.Function = '0.0024'

# create a new 'Calculator'
a3 = Calculator(Input=a2)
a3.ResultArrayName = 'a3'
a3.Function = '0.05'

# create a new 'Legacy VTK Reader'
h_Ctvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Hv18\\Ct_3d1d.vtk'])

# create a new 'Calculator'
pO2_t_1 = Calculator(Input=h_Ctvtk)
pO2_t_1.ResultArrayName = 'PO2_t'
pO2_t_1.Function = 'Ct/0.0000389'

# create a new 'Calculator'
a1_1 = Calculator(Input=pO2_t_1)
a1_1.ResultArrayName = 'a1'
a1_1.Function = '0.22'

# create a new 'Calculator'
a2_1 = Calculator(Input=a1_1)
a2_1.ResultArrayName = 'a2'
a2_1.Function = '0.0024'

# create a new 'Legacy VTK Reader'
t_Cvvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\v18\\Cv_3d1d.vtk'])

# create a new 'Calculator'
pO2_v = Calculator(Input=t_Cvvtk)
pO2_v.ResultArrayName = 'PO2_v'
pO2_v.Function = 'Cv/0.00003'

# create a new 'Legacy VTK Reader'
h_Cvvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Hv18\\Cv_3d1d.vtk'])

# create a new 'Calculator'
a4 = Calculator(Input=a3)
a4.ResultArrayName = 'a4'
a4.Function = '0.0031'

# create a new 'Calculator'
b1 = Calculator(Input=a4)
b1.ResultArrayName = 'b1'
b1.Function = '0.4'

# create a new 'Calculator'
b2 = Calculator(Input=b1)
b2.ResultArrayName = 'b2'
b2.Function = '0.015'

# create a new 'Calculator'
a3_1 = Calculator(Input=a2_1)
a3_1.ResultArrayName = 'a3'
a3_1.Function = '0.05'

# create a new 'Calculator'
a4_1 = Calculator(Input=a3_1)
a4_1.ResultArrayName = 'a4'
a4_1.Function = '0.0031'

# create a new 'Calculator'
b1_1 = Calculator(Input=a4_1)
b1_1.ResultArrayName = 'b1'
b1_1.Function = '0.4'

# create a new 'Calculator'
b2_1 = Calculator(Input=b1_1)
b2_1.ResultArrayName = 'b2'
b2_1.Function = '0.015'

# create a new 'Calculator'
lET = Calculator(Input=b2_1)
lET.ResultArrayName = 'LET'
lET.Function = '2'

# create a new 'Calculator'
k = Calculator(Input=lET)
k.ResultArrayName = 'K'
k.Function = '2.5'

# create a new 'Calculator'
aLPHA = Calculator(Input=k)
aLPHA.ResultArrayName = 'alpha'
aLPHA.Function = '((a1+a2*LET)*PO2_t + (a3+a4*LET)*K)/(PO2_t+K)'

# create a new 'Calculator'
b = Calculator(Input=aLPHA)
b.ResultArrayName = 'b'
b.Function = '(b1*PO2_t + b2*K)/(PO2_t+K)'

# create a new 'Calculator'
bETA = Calculator(Input=b)
bETA.ResultArrayName = 'beta'
bETA.Function = 'b*b'

# create a new 'Calculator'
dose = Calculator(Input=bETA)
dose.ResultArrayName = 'D'
dose.Function = '2'

# create a new 'Calculator'
pO2_v_1 = Calculator(Input=h_Cvvtk)
pO2_v_1.ResultArrayName = 'PO2_v'
pO2_v_1.Function = 'Cv/0.00003'

# create a new 'Clip'
clip2 = Clip(Input=pO2_t)
clip2.ClipType = 'Plane'
clip2.Scalars = ['POINTS', 'PO2_t']
clip2.Value = 21.633607299326844

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Origin = [0.5, 0.5, 0.5]
clip2.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Clip'
clip1 = Clip(Input=pO2_t_1)
clip1.ClipType = 'Plane'
clip1.Scalars = ['POINTS', 'PO2_t']
clip1.Value = 57.5796616021879

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [0.5, 0.5, 0.5]
clip1.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
lET_1 = Calculator(Input=b2)
lET_1.ResultArrayName = 'LET'
lET_1.Function = '2'

# create a new 'Calculator'
k_1 = Calculator(Input=lET_1)
k_1.ResultArrayName = 'K'
k_1.Function = '2.5'

# create a new 'Calculator'
aLPHA_1 = Calculator(Input=k_1)
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
dose_1 = Calculator(Input=bETA_1)
dose_1.ResultArrayName = 'D'
dose_1.Function = '2'

# create a new 'Calculator'
sf = Calculator(Input=dose_1)
sf.ResultArrayName = 'Sf'
sf.Function = 'exp(-alpha*D-beta*D*D)'

# create a new 'Clip'
clip3 = Clip(Input=sf)
clip3.ClipType = 'Plane'
clip3.Scalars = ['POINTS', 'Sf']
clip3.Value = 0.5210542768132349

# init the 'Plane' selected for 'ClipType'
clip3.ClipType.Origin = [0.5, 0.5, 0.5]
clip3.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
sf_1 = Calculator(Input=dose)
sf_1.ResultArrayName = 'Sf'
sf_1.Function = 'exp(-alpha*D-beta*D*D)'

# create a new 'Clip'
clip4 = Clip(Input=sf_1)
clip4.ClipType = 'Plane'
clip4.Scalars = ['POINTS', 'Sf']
clip4.Value = 0.3590366984772403

# init the 'Plane' selected for 'ClipType'
clip4.ClipType.Origin = [0.5, 0.5, 0.5]
clip4.ClipType.Normal = [0.5, 0.5, 0.5]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from h_Ctvtk
h_CtvtkDisplay = Show(h_Ctvtk, renderView1)

# get color transfer function/color map for 'Ct'
ctLUT = GetColorTransferFunction('Ct')
ctLUT.RGBPoints = [7.51473635318689e-05, 0.231373, 0.298039, 0.752941, 0.0013021512968407478, 0.865003, 0.865003, 0.865003, 0.0025291552301496267, 0.705882, 0.0156863, 0.14902]
ctLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Ct'
ctPWF = GetOpacityTransferFunction('Ct')
ctPWF.Points = [7.51473635318689e-05, 0.0, 0.5, 0.0, 0.0025291552301496267, 1.0, 0.5, 0.0]
ctPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
h_CtvtkDisplay.Representation = 'Outline'
h_CtvtkDisplay.ColorArrayName = ['POINTS', 'Ct']
h_CtvtkDisplay.LookupTable = ctLUT
h_CtvtkDisplay.OSPRayScaleArray = 'Ct'
h_CtvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
h_CtvtkDisplay.SelectOrientationVectors = 'Ct'
h_CtvtkDisplay.ScaleFactor = 0.1
h_CtvtkDisplay.SelectScaleArray = 'Ct'
h_CtvtkDisplay.GlyphType = 'Arrow'
h_CtvtkDisplay.GlyphTableIndexArray = 'Ct'
h_CtvtkDisplay.GaussianRadius = 0.005
h_CtvtkDisplay.SetScaleArray = ['POINTS', 'Ct']
h_CtvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
h_CtvtkDisplay.OpacityArray = ['POINTS', 'Ct']
h_CtvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
h_CtvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
h_CtvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
h_CtvtkDisplay.ScalarOpacityFunction = ctPWF
h_CtvtkDisplay.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
h_CtvtkDisplay.ScaleTransferFunction.Points = [0.0019505424425005913, 0.0, 0.5, 0.0, 0.0025291552301496267, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
h_CtvtkDisplay.OpacityTransferFunction.Points = [0.0019505424425005913, 0.0, 0.5, 0.0, 0.0025291552301496267, 1.0, 0.5, 0.0]

# show data from t_Ctvtk
t_CtvtkDisplay = Show(t_Ctvtk, renderView1)

# trace defaults for the display properties.
t_CtvtkDisplay.Representation = 'Outline'
t_CtvtkDisplay.ColorArrayName = ['POINTS', 'Ct']
t_CtvtkDisplay.LookupTable = ctLUT
t_CtvtkDisplay.Position = [1.5, 0.0, 0.0]
t_CtvtkDisplay.OSPRayScaleArray = 'Ct'
t_CtvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
t_CtvtkDisplay.SelectOrientationVectors = 'Ct'
t_CtvtkDisplay.ScaleFactor = 0.1
t_CtvtkDisplay.SelectScaleArray = 'Ct'
t_CtvtkDisplay.GlyphType = 'Arrow'
t_CtvtkDisplay.GlyphTableIndexArray = 'Ct'
t_CtvtkDisplay.GaussianRadius = 0.005
t_CtvtkDisplay.SetScaleArray = ['POINTS', 'Ct']
t_CtvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
t_CtvtkDisplay.OpacityArray = ['POINTS', 'Ct']
t_CtvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
t_CtvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
t_CtvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
t_CtvtkDisplay.ScalarOpacityFunction = ctPWF
t_CtvtkDisplay.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
t_CtvtkDisplay.ScaleTransferFunction.Points = [7.51473635318689e-05, 0.0, 0.5, 0.0, 0.0016079472843557596, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
t_CtvtkDisplay.OpacityTransferFunction.Points = [7.51473635318689e-05, 0.0, 0.5, 0.0, 0.0016079472843557596, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
t_CtvtkDisplay.PolarAxes.Translation = [1.5, 0.0, 0.0]

# show data from pO2_v_1
pO2_v_1Display = Show(pO2_v_1, renderView1)

# get color transfer function/color map for 'PO2_v'
pO2_vLUT = GetColorTransferFunction('PO2_v')
pO2_vLUT.RGBPoints = [34.96970748528838, 0.0, 0.0, 1.0, 100.0000008692344, 1.0, 0.0, 0.0]
pO2_vLUT.ColorSpace = 'HSV'
pO2_vLUT.NanColor = [0.498039215686, 0.498039215686, 0.498039215686]
pO2_vLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'PO2_v'
pO2_vPWF = GetOpacityTransferFunction('PO2_v')
pO2_vPWF.Points = [34.96970748528838, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]
pO2_vPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
pO2_v_1Display.Representation = 'Surface'
pO2_v_1Display.ColorArrayName = ['POINTS', 'PO2_v']
pO2_v_1Display.LookupTable = pO2_vLUT
pO2_v_1Display.LineWidth = 5.0
pO2_v_1Display.RenderLinesAsTubes = 1
pO2_v_1Display.OSPRayScaleArray = 'PO2_v'
pO2_v_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
pO2_v_1Display.SelectOrientationVectors = 'Cv'
pO2_v_1Display.ScaleFactor = 0.1
pO2_v_1Display.SelectScaleArray = 'PO2_v'
pO2_v_1Display.GlyphType = 'Arrow'
pO2_v_1Display.GlyphTableIndexArray = 'PO2_v'
pO2_v_1Display.GaussianRadius = 0.005
pO2_v_1Display.SetScaleArray = ['POINTS', 'PO2_v']
pO2_v_1Display.ScaleTransferFunction = 'PiecewiseFunction'
pO2_v_1Display.OpacityArray = ['POINTS', 'PO2_v']
pO2_v_1Display.OpacityTransferFunction = 'PiecewiseFunction'
pO2_v_1Display.DataAxesGrid = 'GridAxesRepresentation'
pO2_v_1Display.PolarAxes = 'PolarAxesRepresentation'
pO2_v_1Display.ScalarOpacityFunction = pO2_vPWF
pO2_v_1Display.ScalarOpacityUnitDistance = 0.1008201900683401

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
pO2_v_1Display.ScaleTransferFunction.Points = [73.76881937185924, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
pO2_v_1Display.OpacityTransferFunction.Points = [73.76881937185924, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# show data from pO2_v
pO2_vDisplay = Show(pO2_v, renderView1)

# trace defaults for the display properties.
pO2_vDisplay.Representation = 'Surface'
pO2_vDisplay.ColorArrayName = ['POINTS', 'PO2_v']
pO2_vDisplay.LookupTable = pO2_vLUT
pO2_vDisplay.LineWidth = 5.0
pO2_vDisplay.RenderLinesAsTubes = 1
pO2_vDisplay.Position = [1.5, 0.0, 0.0]
pO2_vDisplay.OSPRayScaleArray = 'PO2_v'
pO2_vDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
pO2_vDisplay.SelectOrientationVectors = 'Cv'
pO2_vDisplay.ScaleFactor = 0.1
pO2_vDisplay.SelectScaleArray = 'PO2_v'
pO2_vDisplay.GlyphType = 'Arrow'
pO2_vDisplay.GlyphTableIndexArray = 'PO2_v'
pO2_vDisplay.GaussianRadius = 0.005
pO2_vDisplay.SetScaleArray = ['POINTS', 'PO2_v']
pO2_vDisplay.ScaleTransferFunction = 'PiecewiseFunction'
pO2_vDisplay.OpacityArray = ['POINTS', 'PO2_v']
pO2_vDisplay.OpacityTransferFunction = 'PiecewiseFunction'
pO2_vDisplay.DataAxesGrid = 'GridAxesRepresentation'
pO2_vDisplay.PolarAxes = 'PolarAxesRepresentation'
pO2_vDisplay.ScalarOpacityFunction = pO2_vPWF
pO2_vDisplay.ScalarOpacityUnitDistance = 0.1008201900683401

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
pO2_vDisplay.ScaleTransferFunction.Points = [34.96970748528838, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
pO2_vDisplay.OpacityTransferFunction.Points = [34.96970748528838, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
pO2_vDisplay.PolarAxes.Translation = [1.5, 0.0, 0.0]

# show data from clip3
clip3Display = Show(clip3, renderView1)

# get color transfer function/color map for 'Sf'
sfLUT = GetColorTransferFunction('Sf')
sfLUT.RGBPoints = [0.35, 0.0, 0.0, 1.0, 0.42, 1.0, 0.0, 0.0]
sfLUT.ColorSpace = 'HSV'
sfLUT.NanColor = [0.498039215686, 0.498039215686, 0.498039215686]
sfLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Sf'
sfPWF = GetOpacityTransferFunction('Sf')
sfPWF.Points = [0.35, 0.0, 0.5, 0.0, 0.42, 1.0, 0.5, 0.0]
sfPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
clip3Display.Representation = 'Surface'
clip3Display.ColorArrayName = ['POINTS', 'Sf']
clip3Display.LookupTable = sfLUT
clip3Display.Position = [1.5, 0.0, 0.0]
clip3Display.OSPRayScaleArray = 'Sf'
clip3Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip3Display.SelectOrientationVectors = 'Ct'
clip3Display.ScaleFactor = 0.1
clip3Display.SelectScaleArray = 'Sf'
clip3Display.GlyphType = 'Arrow'
clip3Display.GlyphTableIndexArray = 'Sf'
clip3Display.GaussianRadius = 0.005
clip3Display.SetScaleArray = ['POINTS', 'Sf']
clip3Display.ScaleTransferFunction = 'PiecewiseFunction'
clip3Display.OpacityArray = ['POINTS', 'Sf']
clip3Display.OpacityTransferFunction = 'PiecewiseFunction'
clip3Display.DataAxesGrid = 'GridAxesRepresentation'
clip3Display.PolarAxes = 'PolarAxesRepresentation'
clip3Display.ScalarOpacityFunction = sfPWF
clip3Display.ScalarOpacityUnitDistance = 0.10263591818176866

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip3Display.ScaleTransferFunction.Points = [0.3671317161151878, 0.0, 0.5, 0.0, 0.4193576082450725, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip3Display.OpacityTransferFunction.Points = [0.3671317161151878, 0.0, 0.5, 0.0, 0.4193576082450725, 1.0, 0.5, 0.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
clip3Display.PolarAxes.Translation = [1.5, 0.0, 0.0]

# show data from clip4
clip4Display = Show(clip4, renderView1)

# trace defaults for the display properties.
clip4Display.Representation = 'Surface'
clip4Display.ColorArrayName = ['POINTS', 'Sf']
clip4Display.LookupTable = sfLUT
clip4Display.OSPRayScaleArray = 'Sf'
clip4Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip4Display.SelectOrientationVectors = 'Ct'
clip4Display.ScaleFactor = 0.1
clip4Display.SelectScaleArray = 'Sf'
clip4Display.GlyphType = 'Arrow'
clip4Display.GlyphTableIndexArray = 'Sf'
clip4Display.GaussianRadius = 0.005
clip4Display.SetScaleArray = ['POINTS', 'Sf']
clip4Display.ScaleTransferFunction = 'PiecewiseFunction'
clip4Display.OpacityArray = ['POINTS', 'Sf']
clip4Display.OpacityTransferFunction = 'PiecewiseFunction'
clip4Display.DataAxesGrid = 'GridAxesRepresentation'
clip4Display.PolarAxes = 'PolarAxesRepresentation'
clip4Display.ScalarOpacityFunction = sfPWF
clip4Display.ScalarOpacityUnitDistance = 0.10263591818176866

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip4Display.ScaleTransferFunction.Points = [0.3561836248692617, 0.0, 0.5, 0.0, 0.35963832842885984, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip4Display.OpacityTransferFunction.Points = [0.3561836248692617, 0.0, 0.5, 0.0, 0.35963832842885984, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for pO2_vLUT in view renderView1
pO2_vLUTColorBar = GetScalarBar(pO2_vLUT, renderView1)
pO2_vLUTColorBar.Orientation = 'Horizontal'
pO2_vLUTColorBar.WindowLocation = 'AnyLocation'
pO2_vLUTColorBar.Position = [0.32, 0.044]
pO2_vLUTColorBar.Title = 'PO2_v (mmHg)'
pO2_vLUTColorBar.ComponentTitle = ''
pO2_vLUTColorBar.HorizontalTitle = 1
pO2_vLUTColorBar.UseCustomLabels = 1
pO2_vLUTColorBar.CustomLabels = [35.0, 48.0, 61.0, 74.0, 87.0, 100.0]
pO2_vLUTColorBar.AddRangeLabels = 0
pO2_vLUTColorBar.ScalarBarLength = 0.3300000000000002

# set color bar visibility
pO2_vLUTColorBar.Visibility = 1

# get color legend/bar for sfLUT in view renderView1
sfLUTColorBar = GetScalarBar(sfLUT, renderView1)
sfLUTColorBar.Orientation = 'Horizontal'
sfLUTColorBar.WindowLocation = 'AnyLocation'
sfLUTColorBar.Position = [0.32, 0.9]
sfLUTColorBar.Title = 'Sf (-)'
sfLUTColorBar.ComponentTitle = ''
sfLUTColorBar.HorizontalTitle = 1
sfLUTColorBar.UseCustomLabels = 1
sfLUTColorBar.CustomLabels = [0.35, 0.37, 0.39, 0.41, 0.42]
sfLUTColorBar.AddRangeLabels = 0
sfLUTColorBar.ScalarBarLength = 0.32999999999999985

# set color bar visibility
sfLUTColorBar.Visibility = 1

# show color legend
pO2_v_1Display.SetScalarBarVisibility(renderView1, True)

# show color legend
pO2_vDisplay.SetScalarBarVisibility(renderView1, True)

# show color legend
clip3Display.SetScalarBarVisibility(renderView1, True)

# show color legend
clip4Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
# ----------------------------------------------------------------