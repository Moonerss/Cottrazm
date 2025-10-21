import argparse
import os


def ME_normalize():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--INDIR', type=str, help='Input directory of cell ranger output', nargs='?', default=None)
    parser.add_argument('-o', '--OUTDIR', type=str, help='Outpur directory for analysis result', nargs='?', default=None)
    parser.add_argument('-n', '--NAME', type=str, help='The name of input sample', nargs='?', default=None)

    # Collect args
    args = parser.parse_args()
    argdict = vars(args)

    import stlearn as st
    import scanpy as sc
    import numpy as np
    from numpy import random,mat
    from pathlib import Path
    import pandas as pd
    from scipy import io,sparse

    print (argdict['NAME'], "start SME normalize")

    #read data
    data=st.Read10X(path = argdict['INDIR'])
    data.var_names_make_unique()
    data.layers['raw_count']=data.X
    #tile data
    TILE_PATH=Path(os.path.join(argdict['INDIR'],'{0}_tile'.format(argdict['NAME'])))
    TILE_PATH.mkdir(parents=True,exist_ok=True)
    
    #tile morphology
    st.pp.tiling(data,TILE_PATH,crop_size=40)
    st.pp.extract_feature(data)

    ###process data
    st.pp.normalize_total(data)
    st.pp.log1p(data)

    #gene pca dimention reduction
    st.em.run_pca(data,n_comps=50,random_state=0)
    
    #stSME to normalise log transformed data
    st.spatial.SME.SME_normalize(data, use_data="raw",weights =  "weights_matrix_gd_md")

    #convert SME_norm data to sparesmatrix
    raw_SME_normalized = mat(data.obsm['raw_SME_normalized'])
    raw_SME_normalizedA = sparse.csr_matrix(raw_SME_normalized)
    print ("matrix convert ok!")
    
    io.mmwrite(os.path.join(argdict['OUTDIR'],'{0}_raw_SME_normalizeA.mtx'.format(argdict['NAME'])),raw_SME_normalizedA)
    print("Morphology adjusted ok!")

if __name__ == '__main__':
    ME_normalize()
