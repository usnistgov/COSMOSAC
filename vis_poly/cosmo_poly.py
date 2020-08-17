"""vizualise solvent accessible surface in .cosmo data"""
import os, sys
from io import StringIO
from collections import namedtuple
import pandas
import numpy as np
import open3d as o3d
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../profiles')
import to_sigma as sgm
import csv

class CosmoView:

    def __init__(self, inpath, COSMO_contents):
        self.inpath = inpath
        self.__COSMO_contents = COSMO_contents

    @classmethod
    def read(cls, inpath):
        COSMO_contents = open(inpath).read()
        return cls(inpath, COSMO_contents)

    # private
    def __find_normals(self):
        # calculate normal vectors from point cloud
        normals = []
        for ir, row in self.df.iterrows():
            x, y, z = [_*0.529177249 for _ in (row['x / a.u.'],row['y / a.u.'],row['z / a.u.'])]
            atom = self.df_atom.iloc[int(row['atom'])-1]
            nx, ny, nz = x - atom['x / A'], y - atom['y / A'], z - atom['z / A']
            normals.append([nx,ny,nz])
        return normals

    def __make_polygon(self, COSMO_contents):
        #define
        self.pcd = o3d.geometry.PointCloud()
        self.pcd_atom = o3d.geometry.PointCloud()

        # import .cosmo data using the parser in to_sigma.py
        self.df = sgm.get_seg_DataFrame(COSMO_contents)
        self.df_atom = sgm.get_atom_DataFrame(COSMO_contents)
        self.area_A2, self.volume_A3 = sgm.get_area_volume(COSMO_contents)
        
        xyz = self.df.loc[:, ["x / a.u.", "y / a.u.", "z / a.u."]].values #x, y, z
        val = self.df["charge/area / e/A^2"].values # charge / area
        intensity = np.tile((val - val.min())/(val.max() - val.min()),(3,1)).transpose()
        xyz_atom = self.df_atom.loc[:, ["x / A", "y / A", "z / A"]].values #x, y, z

        # surf
        normals = self.__find_normals()
        self.pcd.points = o3d.utility.Vector3dVector(xyz)
        self.pcd.colors = o3d.utility.Vector3dVector(intensity)
        self.pcd.normals = o3d.utility.Vector3dVector(normals)
        self.pcd.normalize_normals

        # atom
        self.pcd_atom.points = o3d.utility.Vector3dVector(xyz_atom)
      
        # a point cloud is converted to a surface mesh using Ball-Pivoting Algorithm.
        # other possible surface reconstruction algoligthms: Octree, Poission, Powercrust, Delaunay including Alpha shapes, tight cocone etc.
        # some references:
        # https://stackoverflow.com/questions/838761/robust-algorithm-for-surface-reconstruction-from-3d-point-cloud/1990419
        # https://github.com/marcomusy/vedo/blob/master/examples/advanced/recosurface.py
        avg_dist = np.mean(self.pcd.compute_nearest_neighbor_distance())
        radii = [2.0*avg_dist] # change coefficient if the genereted mesh is not watertight.
        voxel_down_pcd = self.pcd.voxel_down_sample(voxel_size=1.2*avg_dist) # this will truncate fine structure.
        self.bpa_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(voxel_down_pcd, o3d.utility.DoubleVector(radii))

    def save(self, outpath, target = 'surf'):
        root, ext = os.path.splitext(outpath)

        # for output filetypes, see http://www.open3d.org/docs/release/tutorial/Basic/file_io.html
        self.__make_polygon(self.__COSMO_contents)
        # for polygon
        if ext == ".ply" or ext == ".stl" or ext == ".obj" or ext == ".off" or ext == ".gltf":
            if target == "surf":
                o3d.io.write_triangle_mesh(outpath, self.bpa_mesh, write_ascii=True)
        elif ext == ".pcd" or ext == ".xyz" or ext == ".xyzn":
            if target == "surf":
                o3d.io.write_point_cloud(outpath, self.pcd, write_ascii=True)
            elif target == "atom":
                o3d.io.write_point_cloud(outpath, self.pcd_atom, write_ascii=True)
        else:
            print("output file type for \"surf\" should be .ply, .stl, .pcd, .obj, .off, .gltf, .xyz or .xyzn ")
            print("output file type for \"atom\" should be .pcd, .xyz or .xyzn ")
        return

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Generate a cosmo surface polygon from a COSMO file')
    parser.add_argument('--inpath', type=str, required=True, nargs=1, help='The path to the cosmo file that is to be processed')
    parser.add_argument('--outpath', type=str, required=True, nargs=1, help='Filetype: ply, .stl, xyz, .pcd, .obj, .off, .gltf. The path to the output of the cosmo surface profile that is to be generated')
    parser.add_argument('--target', type=str, required=True, nargs=1, help='specifysurf or atom. surf: cosmo surface, atom: atom positions')

    # for testing
    arg_str = '--target surf'\
            ' --inpath ../profiles/GAMESS_TEST/ETHANOL.gout'\
            ' --outpath GAMESS_ETHANOL_COSMO.ply'
    args = parser.parse_args(arg_str.split(' '))    

    cosmo = CosmoView.read(args.inpath[0])
    cosmo.save(outpath=args.outpath[0], target = args.target[0])

    # image 
    import pyvista as pv
    if os.path.isfile(args.outpath[0]):
        mesh = pv.read(args.outpath[0])
        mesh = mesh.smooth(n_iter=100)
        cpos = mesh.plot()
