from paraview.simple import *
import sys
import os
import pandas as pd
import json

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

caseDir = sys.argv[1]
case_path = os.path.join("results", caseDir)
output_dir = "./fracturepaths/"

# create a new 'Xdmf3 Reader S'
outputxdmf = Xdmf3ReaderS(registrationName='output.xdmf', FileName=[f"./results/{caseDir}/output.xdmf"])

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
outputxdmfDisplay = Show(outputxdmf, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
outputxdmfDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera(False, 0.97)

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.0, -0.9499999992549419, 13.4]
renderView1.CameraFocalPoint = [0.0, -0.9499999992549419, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show color bar/color legend
outputxdmfDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'phi'
# phiLUT = GetColorTransferFunction('phi')

# get opacity transfer function/opacity map for 'phi'
# phiPWF = GetOpacityTransferFunction('phi')

# get 2D transfer function for 'phi'
# phiTF2D = GetTransferFunction2D('phi')


# create a new 'Threshold'
threshold1 = Threshold(registrationName='Threshold1', Input=outputxdmf)

# Properties modified on threshold1
threshold1.LowerThreshold = 0.9
threshold1.UpperThreshold = 1.1

# show data in view
threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'

# hide data in view
Hide(outputxdmf, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

animationScene1.GoToLast()

# get layout
layout1 = GetLayout()

# split cell
layout1.SplitHorizontal(0, 0.5)

# set active view
SetActiveView(None)

# Create a new 'SpreadSheet View'
spreadSheetView1 = CreateView('SpreadSheetView')
spreadSheetView1.ColumnToSort = ''
spreadSheetView1.BlockSize = 1024

# show data in view
threshold1Display_1 = Show(threshold1, spreadSheetView1, 'SpreadSheetRepresentation')

# assign view to a particular cell in the layout
AssignViewToLayout(view=spreadSheetView1, layout=layout1, hint=2)

# Properties modified on spreadSheetView1
spreadSheetView1.HiddenColumnLabels = []

# Properties modified on spreadSheetView1
spreadSheetView1.HiddenColumnLabels = ['Point ID', 'displacement', 'displacement_Magnitude', 'phi', 'Points_Magnitude', 'Block Number']

# export view
ExportView(f"./fracturepaths/{caseDir}.csv", view=spreadSheetView1)


# === 5. Exportar propiedades del material a un nuevo CSV ===
props_path = os.path.join(case_path, "simulation_output.txt")
with open(props_path, 'r') as f:
    lines = f.readlines()[2:]
    props = json.loads(''.join(lines))

props_df = pd.DataFrame([props])
props_csv_path = os.path.join(output_dir, f"propiedades_{caseDir}.csv")
props_df.to_csv(props_csv_path, index=False)

print(f"[{caseDir}] Propiedades exportadas: {props_csv_path}")
