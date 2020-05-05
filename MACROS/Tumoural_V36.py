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
renderView1.CenterOfRotation = [0.5, 0.5, 0.5]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [2.2526169371280464, 3.096966631448539, -0.6748406741092956]
renderView1.CameraFocalPoint = [0.49999999999999994, 0.4999999999999994, 0.500000000000001]
renderView1.CameraViewUp = [-0.10508735842448395, -0.35015414251934207, -0.9307785577546941]
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
t_V36_Ctvtk = LegacyVTKReader(FileNames=['C:\\Users\\Edoardo Rosati\\Desktop\\Tumoural Voronoi\\Tv36\\Ct_3d1d.vtk'])

# create a new 'Calculator'
pO2 = Calculator(Input=t_V36_Ctvtk)
pO2.ResultArrayName = 'PO2_t'
pO2.Function = 'Ct/0.0000389'

# create a new 'Calculator'
a1 = Calculator(Input=pO2)
a1.ResultArrayName = 'a1'
a1.Function = '0.917'

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
b1.Function = '0.36'

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
d = Calculator(Input=bETA)
d.ResultArrayName = 'D'
d.Function = '10'

# create a new 'Calculator'
sf = Calculator(Input=d)
sf.ResultArrayName = 'Sf'
sf.Function = 'exp(-alpha*D-beta*D*D)'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from sf
sfDisplay = Show(sf, renderView1)

# get color transfer function/color map for 'Sf'
sfLUT = GetColorTransferFunction('Sf')
sfLUT.RGBPoints = [9.057703677872791e-10, 0.231373, 0.298039, 0.752941, 2.3200902045646106e-09, 0.865003, 0.865003, 0.865003, 3.734410041341979e-09, 0.705882, 0.0156863, 0.14902]
sfLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Sf'
sfPWF = GetOpacityTransferFunction('Sf')
sfPWF.Points = [9.057703677872791e-10, 0.0, 0.5, 0.0, 3.734410041341979e-09, 1.0, 0.5, 0.0]
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
sfDisplay.ScaleTransferFunction.Points = [0.7998550287874373, 0.0, 0.5, 0.0, 0.8074162605550553, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
sfDisplay.OpacityTransferFunction.Points = [0.7998550287874373, 0.0, 0.5, 0.0, 0.8074162605550553, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

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
sfDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
# ----------------------------------------------------------------