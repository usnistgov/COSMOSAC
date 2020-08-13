"""
Script to inject the name into the generated profiles from complist.txt and replace NaN values with null values (what JSON standard requires)
"""
import os
import glob
import pandas
import json

df = pandas.read_csv('complist.txt', sep=r'\s+', engine='python', index_col='INCHIKEY')
for f in glob.glob('sigma3/*.sigma'):
    with open(f) as fp:
        lines = fp.readlines()
        meta = json.loads(lines[0].replace('# meta: ',''))
        stdikey = os.path.split(f)[1].split('.')[0]
        if 'CAS' in meta:
            del meta['CAS']
        meta['name'] = df.loc[stdikey, 'NAME']
        meta['standard_INCHIKEY'] = stdikey
        lines[0] = '# meta: ' + json.dumps(meta) + '\n'
    with open(f, 'w') as fp:
        fp.write(''.join(lines))

for f in glob.glob('sigma3/*.sigma'):
    with open(f) as fp:
        contents = fp.read().replace(': NaN',': null')
    with open(f, 'w') as fp:
        fp.write(contents)