# trace generated using paraview version 5.11.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

data_arrays = [
    'x_axis',
    'y_axis',
    'z_axis'
]

scale_factor_input = input("Arrow scale factor")
scale_factor = float(scale_factor_input if scale_factor_input != '' else 0.1)

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
activeSource = GetActiveSource()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=activeSource, GlyphType='Arrow')
glyph1.ScaleFactor = scale_factor
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.OrientationArray = ['POINTS', data_arrays[0]]
glyph1.GlyphMode = 'All Points'
glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')
glyph1Display.Representation = 'Surface'
glyph1Display.AmbientColor = [1.0, 0.0, 0.0]
glyph1Display.DiffuseColor = [1.0, 0.0, 0.0]
glyph1Display.SetScalarColoring(None, 0)

# set active source
SetActiveSource(activeSource)

# create a new 'Glyph'
glyph2 = Glyph(registrationName='Glyph2', Input=activeSource, GlyphType='Arrow')
glyph2.ScaleFactor = scale_factor
glyph2.ScaleArray = ['POINTS', 'No scale array']
glyph2.OrientationArray = ['POINTS', data_arrays[1]]
glyph2.GlyphMode = 'All Points'
glyph2Display = Show(glyph2, renderView1, 'GeometryRepresentation')
glyph2Display.Representation = 'Surface'
glyph2Display.AmbientColor = [1.0, 1.0, 0.0]
glyph2Display.DiffuseColor = [1.0, 1.0, 0.0]
glyph2Display.SetScalarColoring(None, 0)

# set active source
SetActiveSource(activeSource)

# create a new 'Glyph'
glyph3 = Glyph(registrationName='Glyph3', Input=activeSource, GlyphType='Arrow')
glyph3.ScaleFactor = scale_factor
glyph3.ScaleArray = ['POINTS', 'No scale array']
glyph3.OrientationArray = ['POINTS', data_arrays[2]]
glyph3.GlyphMode = 'All Points'
glyph3Display = Show(glyph3, renderView1, 'GeometryRepresentation')
glyph3Display.Representation = 'Surface'
glyph3Display.AmbientColor = [0.3333333333333333, 1.0, 0.0]
glyph3Display.DiffuseColor = [0.3333333333333333, 1.0, 0.0]
glyph3Display.SetScalarColoring(None, 0)

renderView1.Update()
