# cosmo_poly.py

This `cosmo_poly.py` script visualizes cosmo surface from cosmo segment data.

![cosmo surface image](ethanol.png)


## usage

In `vis-poly` folder,
```termnial
python vis_poly.py --inpath ../profiles/GAMESS_TEST/ETHANOL.gout --outpath GAMESS_ETHANOL_COSMO.ply --target surf 
```

### input filetypes (`inpath`):

In the `inpath`  argument, cosmo data must be given. The script supports the following cosmo formats:
- `.cosmo`  from DMol3
- `.cosmo`  from Gaussian09
- `.gout`  from GAMESS

sample files can be found in profiles/DMol3_TEST, profiles/GAMESS_TEST, profiles/GAUSSIAN09_TEST and profiles/UD/cosmo.

### output filetypes (`outpath`):

In the `outpath`  argument, a output file with one of the following extensions must be given. The script supports general mesh formats 
- `ply`: stanford polygon
- `stl`: Standard Triangulated Language
- `obj`: wavefront OBJ
- `off`: object File Format
- `gltf`: GL Transmission Format

and point position data formats
- `.pcd`  point cloud data
- `.xyz`  xyz coordinate data 
- `.xyzn`  xyz coordinate + normal vectors

### types for processing (`types`):

Give `surf`  or `atom` for the `target`  argument. If `surf` is given, segment data is processed. Mesh or point data are saved depending on the extension given in the `outpath` argument.
If `atom`  is given, the position data is saved. In this case, only point data format is accepted.

## Dpendencies

-Open3d: https://pypi.org/project/open3d/


-pyvista: https://pypi.org/project/pyvista/

## What Does This Script Do?

The script reads cosmo data using a parser in to_sigma.py. Then, a surface mesh is constructed from the cosmo segment data using ball pivoting algorithm (open3d functionality).  Pyvista is a handy package for vlizualization of mesh data. 

## How can I vizualize mesh data?

Many software supports mesh data. Some famous packages include
- paraview
- meshlab
- gmsh