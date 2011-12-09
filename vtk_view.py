# Based on http://www.vtk.org/Wiki/VTK/Examples/Python/GeometricObjects/Display/Sphere
import vtk
import ase
 
# create render window, renderer and render window interactor:
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

# create actors
atoms = ase.io.read('CONTCAR', format='vasp')
from ase.visualize.vtk import *
requirevtk()
from ase.visualize.vtk.atoms import *
v = vtkAtoms(atoms)
v.add_cell()

# assign actors to the renderer
for m in v.modules:
    ren.AddActor(v.get_actor(m))
  
# enable user interface interactor
iren.Initialize()
renWin.Render()
iren.Start()

