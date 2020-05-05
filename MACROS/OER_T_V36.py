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
renderView1.CameraPosition = [1.3155305716549501, 0.6223763707314816, 3.74285157435387]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
renderView1.CameraViewUp = [-0.014137916063968754, 0.999316510058891, -0.03415599585851173]
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
a1.Function = '0.179'

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
b1.Function = '0.22'

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
aLPHA_A = Calculator(Input=k)
aLPHA_A.ResultArrayName = 'alpha_a'
aLPHA_A.Function = '0.179'

# create a new 'Calculator'
b = Calculator(Input=aLPHA_A)
b.ResultArrayName = 'b'
b.Function = '(b1*PO2_t + b2*K)/(PO2_t+K)'

# create a new 'Calculator'
bETA_A = Calculator(Input=b)
bETA_A.ResultArrayName = 'beta_a'
bETA_A.Function = '0.22'

# create a new 'Calculator'
d = Calculator(Input=bETA_A)
d.ResultArrayName = 'D'
d.Function = '0'

# create a new 'Calculator'
sf = Calculator(Input=d)
sf.ResultArrayName = 'Sf'
sf.Function = 'exp(-alpha_a*D-beta_a*D*D)'

# create a new 'Calculator'
aLPHA_H = Calculator(Input=sf)
aLPHA_H.ResultArrayName = 'alpha_h'
aLPHA_H.Function = """((a1+a2*LET)*PO2_t + (a3+a4*LET)*K)/(PO2_t+K)
"""

# create a new 'Calculator'
bETA_H = Calculator(Input=aLPHA_H)
bETA_H.ResultArrayName = 'beta_h'
bETA_H.Function = 'b*b'

# create a new 'Calculator'
oER_num = Calculator(Input=bETA_H)
oER_num.ResultArrayName = 'num'
oER_num.Function = '2*D*beta_a'

# create a new 'Calculator'
oER_denum = Calculator(Input=oER_num)
oER_denum.ResultArrayName = 'denum'
oER_denum.Function = 'sqrt(alpha_a*alpha_a+4*beta_a*(alpha_h*D+beta_h*D*D)) - alpha_a'

# create a new 'Calculator'
oER = Calculator(Input=oER_denum)
oER.ResultArrayName = 'OER'
oER.Function = 'num/denum'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from oER
oERDisplay = Show(oER, renderView1)

# get color transfer function/color map for 'OER'
oERLUT = GetColorTransferFunction('OER')
oERLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 5.878906683738849e-39, 0.865003, 0.865003, 0.865003, 1.1757813367477812e-38, 0.705882, 0.0156863, 0.14902]
oERLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'OER'
oERPWF = GetOpacityTransferFunction('OER')
oERPWF.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]
oERPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
oERDisplay.Representation = 'Surface'
oERDisplay.ColorArrayName = ['POINTS', 'OER']
oERDisplay.LookupTable = oERLUT
oERDisplay.OSPRayScaleArray = 'OER'
oERDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
oERDisplay.SelectOrientationVectors = 'Ct'
oERDisplay.ScaleFactor = 0.1
oERDisplay.SelectScaleArray = 'OER'
oERDisplay.GlyphType = 'Arrow'
oERDisplay.GlyphTableIndexArray = 'OER'
oERDisplay.GaussianRadius = 0.005
oERDisplay.SetScaleArray = ['POINTS', 'OER']
oERDisplay.ScaleTransferFunction = 'PiecewiseFunction'
oERDisplay.OpacityArray = ['POINTS', 'OER']
oERDisplay.OpacityTransferFunction = 'PiecewiseFunction'
oERDisplay.DataAxesGrid = 'GridAxesRepresentation'
oERDisplay.PolarAxes = 'PolarAxesRepresentation'
oERDisplay.ScalarOpacityFunction = oERPWF
oERDisplay.ScalarOpacityUnitDistance = 0.08665311754517603

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
oERDisplay.ScaleTransferFunction.Points = [1.0717736273780465, 0.0, 0.5, 0.0, 1.1118672278673896, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
oERDisplay.OpacityTransferFunction.Points = [1.0717736273780465, 0.0, 0.5, 0.0, 1.1118672278673896, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for oERLUT in view renderView1
oERLUTColorBar = GetScalarBar(oERLUT, renderView1)
oERLUTColorBar.WindowLocation = 'UpperCenter'
oERLUTColorBar.Title = 'OER'
oERLUTColorBar.ComponentTitle = ''
oERLUTColorBar.HorizontalTitle = 1
oERLUTColorBar.ScalarBarLength = 0.32999999999999996

# set color bar visibility
oERLUTColorBar.Visibility = 1

# show color legend
oERDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(None)
# ----------------------------------------------------------------