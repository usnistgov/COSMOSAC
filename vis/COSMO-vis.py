import io, pandas, os, sys
import numpy as np
import re
import shutil
import textwrap

sys.path.append('../profiles')
from to_sigma import get_seg_DataFrame, get_atom_DataFrame

def write_and_load(name, ofname = None, overwrite = True, launch = True, bgcolor = '0xffffff'):

    def name_to_inchikey(name):
        df = pandas.read_csv(os.path.join(os.path.dirname(__file__), '../profiles/UD/complist.txt'),sep=' ')
        match = df[(df.NAME==name) | (df.INCHIKEY==name)].INCHIKEY
        if len(match) == 0:
            raise ValueError('Could not match the name:'+name)
        key = match.iloc[0]
        return key
    inchikey = name_to_inchikey(name)
    path = os.path.join(os.path.dirname(__file__), '../profiles/UD/cosmo/'+inchikey+'.cosmo')
    contents = open(path).read()
    atoms = get_atom_DataFrame(contents)
    segments = get_seg_DataFrame(contents)

    out_atoms = []
    for ir, row in atoms.iterrows():
        # See http://www.ch.ic.ac.uk/rzepa/mim/domestic/html/atomcols.htm
        colors = dict(H='white',O='red',Cl='green',N='blue',C='grey',S='yellow',P='orange')
        out_atoms.append([colors.get(row['atom'],'yellow'),row['x / A'],row['y / A'],row['z / A']])

    out_segments = []
    for ir, row in segments.iterrows():
        x, y, z = [_*0.529177249 for _ in (row['x / a.u.'],row['y / a.u.'],row['z / a.u.'])]
        atom = atoms.iloc[int(row['atom'])-1]
        nx, ny, nz = x-atom['x / A'], y-atom['y / A'], z-atom['z / A']
        rn = (row['area / A^2']/np.pi)**(0.5)
        charge = row['charge / e']
        out_segments.append([rn,charge,x,y,z,nx,ny,nz])

    # Open the base file
    with open(os.path.join(os.path.dirname(__file__),'vis.html')) as fp:
        contents = fp.read()

    # Replace the includes inside the HTML file
    for file in ['Lut.js', 'OrbitControls.js', "three.min.js"]:
        the_js = '<script>\n' + open(os.path.join(os.path.dirname(__file__), file), 'r').read() + '\n</script>\n'
        contents = re.sub(r'<script src="'+file+'"></script>', lambda x: the_js, contents)

    # Inject the data from the dependent JS files
    for file in ["Data.js"]:
        contents = re.sub(r'<script src="'+file+'"></script>', 
                          r'\n<script>var spheres = {0:s};\nvar segments = {1:s};\n</script>'.format(str(out_atoms), str(out_segments)), 
                          contents)

    # Set the background color
    contents = contents.replace('%BGCOLOR%',bgcolor)

    if os.path.exists(ofname) and not overwrite:
        raise ValueError('Could not write file '+ofname+'; file already exists.  Pass --overwrite flag to force overwrite')

    with open(ofname,'w') as fp:
        fp.write(contents)

    if launch:
        os.startfile(ofname)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Visualize the patches and atoms in a COSMO file')
    parser.add_argument('--name', type=str, required=True, nargs=1, help='The name or InChI key for the molecule')
    parser.add_argument('--ofname', type=str, required=True, nargs=1, help='The path of the output file')
    parser.add_argument('--bgcolor', type=str, required=False, nargs=1, default='0xffffff', help='The background color for the scene, specified as a string; 0xffffff would be white')
    parser.add_argument('--overwrite', action='store_const',
                        const=True, default=False,
                        help='overwrite the output file without asking')
    parser.add_argument('--no-launch', action='store_const',
                        const=True, default=False,
                        help="Don't launch the file")

    args = parser.parse_args()
    # args = parser.parse_args('--name NEOPENTANE --ofname SAC-vis.html --overwrite'.split(' '))  # For testing
    write_and_load(name = args.name[0], 
                   ofname = args.ofname[0], 
                   overwrite = args.overwrite, 
                   launch = not args.no_launch, 
                   bgcolor = args.bgcolor[0])