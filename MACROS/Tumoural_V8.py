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
renderView1.CenterOfRotation = [0.5, 0.5, 0.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [1.3155305716549552, 0.6223763707314822, 3.742851574353891]
renderView1.CameraFocalPoint = [0.49999999999999983, 0.49999999999999983, 0.49999999999999983]
renderView1.CameraViewUp = [-0.014137916063968754, 0.999316510058891, -0.03415599585851173]
renderView1.CameraParallelScale = 0.8660254037844386
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
t_Cv_V8vtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Tv8\\Cv_3d1d.vtk'])

# create a new 'Calculator'
pO2_v = Calculator(Input=t_Cv_V8vtk)
pO2_v.ResultArrayName = 'PO2_v'
pO2_v.Function = 'Cv/0.00003'

# create a new 'Legacy VTK Reader'
t_Ct_V8vtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Tv8\\Ct_3d1d.vtk'])

# create a new 'Calculator'
pO2_t = Calculator(Input=t_Ct_V8vtk)
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
lET = Calculator(Input=b2)
lET.ResultArrayName = 'LET'
lET.Function = '2'

# create a new 'Calculator'
k = Calculator(Input=lET)
k.ResultArrayName = 'K'
k.Function = '2.5'

# create a new 'Calculator'
aLPHA = Calculator(Input=k)
aLPHA.ResultArrayName = 'alpha'
aLPHA.Function = """((a1+a2*LET)*PO2_t + (a3+a4*LET)*K)/(PO2_t+K)
"""

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
sf = Calculator(Input=dose)
sf.ResultArrayName = 'Sf'
sf.Function = 'exp(-alpha*D-beta*D*D)'

# create a new 'Clip'
clip1 = Clip(Input=sf)
clip1.ClipType = 'Plane'
clip1.Scalars = ['POINTS', 'Sf']
clip1.Value = 0.6622941659655256

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [0.5, 0.5, 0.5]
clip1.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=sf)

# create a new 'Clip'
clip4 = Clip(Input=pO2_t)
clip4.ClipType = 'Plane'
clip4.Scalars = ['POINTS', 'PO2_t']
clip4.Value = 6.9507447071822375

# init the 'Plane' selected for 'ClipType'
clip4.ClipType.Origin = [0.5, 0.5, 0.5]
clip4.ClipType.Normal = [0.5, 0.5, 0.5]

# create a new 'Calculator'
calculator1 = Calculator(Input=integrateVariables1)
calculator1.ResultArrayName = 'N'
calculator1.Function = '375'

# create a new 'Calculator'
calculator2 = Calculator(Input=calculator1)
calculator2.ResultArrayName = 'TCP'
calculator2.Function = 'exp(-N*Sf)'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from t_Ct_V8vtk
t_Ct_V8vtkDisplay = Show(t_Ct_V8vtk, renderView1)

# get color transfer function/color map for 'Ct'
ctLUT = GetColorTransferFunction('Ct')
ctLUT.RGBPoints = [-3.0944102036301047e-06, 0.231373, 0.298039, 0.752941, 0.0009861352864390938, 0.865003, 0.865003, 0.865003, 0.0019753649830818176, 0.705882, 0.0156863, 0.14902]
ctLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Ct'
ctPWF = GetOpacityTransferFunction('Ct')
ctPWF.Points = [-3.0944102036301047e-06, 0.0, 0.5, 0.0, 0.0019753649830818176, 1.0, 0.5, 0.0]
ctPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
t_Ct_V8vtkDisplay.Representation = 'Outline'
t_Ct_V8vtkDisplay.ColorArrayName = ['POINTS', 'Ct']
t_Ct_V8vtkDisplay.LookupTable = ctLUT
t_Ct_V8vtkDisplay.OSPRayScaleArray = 'Ct'
t_Ct_V8vtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
t_Ct_V8vtkDisplay.SelectOrientationVectors = 'Ct'
t_Ct_V8vtkDisplay.ScaleFactor = 0.1
t_Ct_V8vtkDisplay.SelectScaleArray = 'Ct'
t_Ct_V8vtkDisplay.GlyphType = 'Arrow'
t_Ct_V8vtkDisplay.GlyphTableIndexArray = 'Ct'
t_Ct_V8vtkDisplay.GaussianRadius = 0.005
t_Ct_V8vtkDisplay.SetScaleArray = ['POINTS', 'Ct']
t_Ct_V8vtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
t_Ct_V8vtkDisplay.OpacityArray = ['POINTS', 'Ct']
t_Ct_V8vtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
t_Ct_V8vtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
t_Ct_V8vtkDisplay.PolarAxes = 'PolarAxesRepresentation'
t_Ct_V8vtkDisplay.ScalarOpacityFunction = ctPWF
t_Ct_V8vtkDisplay.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
t_Ct_V8vtkDisplay.ScaleTransferFunction.Points = [-3.0944102036301047e-06, 0.0, 0.5, 0.0, 0.0005438623484224081, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
t_Ct_V8vtkDisplay.OpacityTransferFunction.Points = [-3.0944102036301047e-06, 0.0, 0.5, 0.0, 0.0005438623484224081, 1.0, 0.5, 0.0]

# show data from pO2_v
pO2_vDisplay = Show(pO2_v, renderView1)

# trace defaults for the display properties.
pO2_vDisplay.Representation = 'Surface'
pO2_vDisplay.ColorArrayName = ['POINTS', '']
pO2_vDisplay.LineWidth = 10.0
pO2_vDisplay.RenderLinesAsTubes = 1
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
pO2_vDisplay.ScalarOpacityUnitDistance = 0.13111097749257294

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
pO2_vDisplay.ScaleTransferFunction.Points = [13.660525049393375, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
pO2_vDisplay.OpacityTransferFunction.Points = [13.660525049393375, 0.0, 0.5, 0.0, 100.0000008692344, 1.0, 0.5, 0.0]

# show data from clip1
clip1Display = Show(clip1, renderView1)

# get color transfer function/color map for 'Sf'
sfLUT = GetColorTransferFunction('Sf')
sfLUT.RGBPoints = [0.42, 0.0, 0.0, 1.0, 0.9036257445937314, 1.0, 0.0, 0.0]
sfLUT.ColorSpace = 'HSV'
sfLUT.NanColor = [0.498039215686, 0.498039215686, 0.498039215686]
sfLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Sf'
sfPWF = GetOpacityTransferFunction('Sf')
sfPWF.Points = [0.42, 0.0, 0.5, 0.0, 0.9036257445937314, 1.0, 0.5, 0.0]
sfPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = ['POINTS', 'Sf']
clip1Display.LookupTable = sfLUT
clip1Display.OSPRayScaleArray = 'Sf'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'Ct'
clip1Display.ScaleFactor = 0.1
clip1Display.SelectScaleArray = 'Sf'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'Sf'
clip1Display.GaussianRadius = 0.005
clip1Display.SetScaleArray = ['POINTS', 'Sf']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'Sf']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityFunction = sfPWF
clip1Display.ScalarOpacityUnitDistance = 0.10263591818176866

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [0.4209625873373199, 0.0, 0.5, 0.0, 0.9036257445937314, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [0.4209625873373199, 0.0, 0.5, 0.0, 0.9036257445937314, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for sfLUT in view renderView1
sfLUTColorBar = GetScalarBar(sfLUT, renderView1)
sfLUTColorBar.Orientation = 'Horizontal'
sfLUTColorBar.WindowLocation = 'AnyLocation'
sfLUTColorBar.Position = [0.32, 0.044]
sfLUTColorBar.Title = 'Sf (-)'
sfLUTColorBar.ComponentTitle = ''
sfLUTColorBar.HorizontalTitle = 1
sfLUTColorBar.UseCustomLabels = 1
sfLUTColorBar.CustomLabels = [0.42, 0.5, 0.58, 0.66, 0.74, 0.82, 0.9]
sfLUTColorBar.AddRangeLabels = 0
sfLUTColorBar.ScalarBarThickness = 20
sfLUTColorBar.ScalarBarLength = 0.35

# set color bar visibility
sfLUTColorBar.Visibility = 1

# show color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'spreadSheetView1'
# ----------------------------------------------------------------

# show data from calculator2
calculator2Display = Show(calculator2, spreadSheetView1)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
# ----------------------------------------------------------------