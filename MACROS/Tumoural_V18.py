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
renderView1.CenterOfRotation = [0.5, 0.5, 0.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-3.0499580015642502, 1.836973929106569, -0.9151274367011436]
renderView1.CameraFocalPoint = [0.5000000000000012, 0.4999999999999997, 0.5000000000000006]
renderView1.CameraViewUp = [0.4238884834770319, 0.8739818535700349, -0.23764316360002488]
renderView1.CameraParallelScale = 0.8660254037844386
renderView1.Background = [0.32, 0.34, 0.43]
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024L
spreadSheetView1.HiddenColumnLabels = ['Point ID', 'Ct', 'D', 'K', 'LET', 'PO2', 'Points', 'Points_Magnitude', 'a1', 'a2', 'a3', 'a4', 'alpha', 'b', 'b1', 'b2', 'beta']
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
ct_3d1dvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\v18\\Ct_3d1d.vtk'])

# create a new 'Calculator'
pO2 = Calculator(Input=ct_3d1dvtk)
pO2.ResultArrayName = 'PO2'
pO2.Function = 'Ct/0.0000389'

# create a new 'Calculator'
a1 = Calculator(Input=pO2)
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
aLPHA.Function = '((a1+a2*LET)*PO2 + (a3+a4*LET)*K)/(PO2+K)'

# create a new 'Calculator'
b = Calculator(Input=aLPHA)
b.ResultArrayName = 'b'
b.Function = '(b1*PO2 + b2*K)/(PO2+K)'

# create a new 'Calculator'
bETA = Calculator(Input=b)
bETA.ResultArrayName = 'beta'
bETA.Function = 'b*b'

# create a new 'Calculator'
d = Calculator(Input=bETA)
d.ResultArrayName = 'D'
d.Function = '10'

# create a new 'Calculator'
sf = Calculator(Input=d)
sf.ResultArrayName = 'Sf'
sf.Function = 'exp(-alpha*D-beta*D*D)'

# create a new 'Integrate Variables'
integrateVariables1 = IntegrateVariables(Input=sf)

# create a new 'Calculator'
calculator13 = Calculator(Input=integrateVariables1)
calculator13.Function = 'Sf/0.3333'

# create a new 'Calculator'
calculator14 = Calculator(Input=calculator13)
calculator14.ResultArrayName = 'N'
calculator14.Function = '10000'

# create a new 'Calculator'
calculator15 = Calculator(Input=calculator14)
calculator15.ResultArrayName = 'TCP'
calculator15.Function = 'exp(-N*Result)'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from pO2
pO2Display = Show(pO2, renderView1)

# get color transfer function/color map for 'PO2'
pO2LUT = GetColorTransferFunction('PO2')
pO2LUT.RGBPoints = [1.9318088311534425, 0.231373, 0.298039, 0.752941, 21.633607299326847, 0.865003, 0.865003, 0.865003, 41.33540576750025, 0.705882, 0.0156863, 0.14902]
pO2LUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'PO2'
pO2PWF = GetOpacityTransferFunction('PO2')
pO2PWF.Points = [1.9318088311534425, 0.0, 0.5, 0.0, 41.33540576750025, 1.0, 0.5, 0.0]
pO2PWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
pO2Display.Representation = 'Surface'
pO2Display.ColorArrayName = ['POINTS', 'PO2']
pO2Display.LookupTable = pO2LUT
pO2Display.OSPRayScaleArray = 'PO2'
pO2Display.OSPRayScaleFunction = 'PiecewiseFunction'
pO2Display.SelectOrientationVectors = 'Ct'
pO2Display.ScaleFactor = 0.1
pO2Display.SelectScaleArray = 'PO2'
pO2Display.GlyphType = 'Arrow'
pO2Display.GlyphTableIndexArray = 'PO2'
pO2Display.GaussianRadius = 0.005
pO2Display.SetScaleArray = ['POINTS', 'PO2']
pO2Display.ScaleTransferFunction = 'PiecewiseFunction'
pO2Display.OpacityArray = ['POINTS', 'PO2']
pO2Display.OpacityTransferFunction = 'PiecewiseFunction'
pO2Display.DataAxesGrid = 'GridAxesRepresentation'
pO2Display.PolarAxes = 'PolarAxesRepresentation'
pO2Display.ScalarOpacityFunction = pO2PWF
pO2Display.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
pO2Display.ScaleTransferFunction.Points = [1.9318088311534425, 0.0, 0.5, 0.0, 41.33540576750025, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
pO2Display.OpacityTransferFunction.Points = [1.9318088311534425, 0.0, 0.5, 0.0, 41.33540576750025, 1.0, 0.5, 0.0]

# show data from sf
sfDisplay = Show(sf, renderView1)

# get color transfer function/color map for 'Sf'
sfLUT = GetColorTransferFunction('Sf')
sfLUT.RGBPoints = [7.222457363042114e-08, 0.231373, 0.298039, 0.752941, 0.004832564430341279, 0.865003, 0.865003, 0.865003, 0.009665056636108931, 0.705882, 0.0156863, 0.14902]
sfLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Sf'
sfPWF = GetOpacityTransferFunction('Sf')
sfPWF.Points = [7.222457363042114e-08, 0.0, 0.5, 0.0, 0.009665056636108931, 1.0, 0.5, 0.0]
sfPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
sfDisplay.Representation = 'Surface'
sfDisplay.ColorArrayName = ['POINTS', 'Sf']
sfDisplay.LookupTable = sfLUT
sfDisplay.OSPRayScaleArray = 'Sf'
sfDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
sfDisplay.SelectOrientationVectors = 'Ct'
sfDisplay.ScaleFactor = 0.1
sfDisplay.SelectScaleArray = 'Sf'
sfDisplay.GlyphType = 'Arrow'
sfDisplay.GlyphTableIndexArray = 'Sf'
sfDisplay.GaussianRadius = 0.005
sfDisplay.SetScaleArray = ['POINTS', 'Sf']
sfDisplay.ScaleTransferFunction = 'PiecewiseFunction'
sfDisplay.OpacityArray = ['POINTS', 'Sf']
sfDisplay.OpacityTransferFunction = 'PiecewiseFunction'
sfDisplay.DataAxesGrid = 'GridAxesRepresentation'
sfDisplay.PolarAxes = 'PolarAxesRepresentation'
sfDisplay.ScalarOpacityFunction = sfPWF
sfDisplay.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
sfDisplay.ScaleTransferFunction.Points = [0.3671317161151878, 0.0, 0.5, 0.0, 0.6749768375112821, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
sfDisplay.OpacityTransferFunction.Points = [0.3671317161151878, 0.0, 0.5, 0.0, 0.6749768375112821, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for pO2LUT in view renderView1
pO2LUTColorBar = GetScalarBar(pO2LUT, renderView1)
pO2LUTColorBar.WindowLocation = 'UpperRightCorner'
pO2LUTColorBar.Title = 'PO2'
pO2LUTColorBar.ComponentTitle = ''
pO2LUTColorBar.HorizontalTitle = 1
pO2LUTColorBar.ScalarBarLength = 0.32999999999999996

# set color bar visibility
pO2LUTColorBar.Visibility = 1

# get color legend/bar for sfLUT in view renderView1
sfLUTColorBar = GetScalarBar(sfLUT, renderView1)
sfLUTColorBar.WindowLocation = 'UpperCenter'
sfLUTColorBar.Title = 'Sf'
sfLUTColorBar.ComponentTitle = ''
sfLUTColorBar.HorizontalTitle = 1
sfLUTColorBar.ScalarBarLength = 0.32999999999999996

# set color bar visibility
sfLUTColorBar.Visibility = 1

# show color legend
pO2Display.SetScalarBarVisibility(renderView1, True)

# show color legend
sfDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'spreadSheetView1'
# ----------------------------------------------------------------

# show data from calculator15
calculator15Display = Show(calculator15, spreadSheetView1)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
# ----------------------------------------------------------------