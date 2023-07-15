import numpy as np
from numba import jit
from scipy.sparse import csr_matrix

@jit(nopython=True)
def filtergenes(sitesAndgenes_matrix,
                matrix_rows,
                Pvalues,
                useful_keys,
                P_cutoff=0.05):
  
    shape = sitesAndgenes_matrix.shape
    useful_matrix = np.zeros(shape=(len(useful_keys),shape[1]))
    for ik in range(len(useful_keys)):
        k = useful_keys[ik]
        useful_matrix[ik]=sitesAndgenes_matrix[matrix_rows.index(k)]

    Pvalues = Pvalues.astype(np.single)
    t2 = (useful_matrix.transpose()*Pvalues)
    P_matrix = t2.transpose()
    t2 = useful_matrix.sum(axis=0)
    t4 = np.where(np.all([P_matrix<P_cutoff,useful_matrix>0],axis=0),1,0)
    t4 = t4.sum(axis=0)
    testable_index= np.where(np.all([t2>=2,t4.astype(np.single)/t2>0.5],axis=0))[0]
    control_index= np.where(np.all([t2>=3,t4==0],axis=0))[0]

    return (testable_index,control_index)

def filtergenes_sparse(sitesAndgenes_sparse_matrix,
                matrix_rows,
                matrix_columns,
                Pvalues,
                useful_keys,
                exclude_intronic,
                P_cutoff=0.05):
    
    shape = sitesAndgenes_sparse_matrix.shape
    assert shape[1] == len(matrix_columns) 

    selecting_rows = [ matrix_rows.index(k) for k in useful_keys]
    print ("Rows selected")
    useful_sparse_matrix = sitesAndgenes_sparse_matrix[selecting_rows]
   
    t = useful_sparse_matrix.toarray()
    if exclude_intronic:
        useful_sparse_matrix = np.where(np.all([t>0,t<10],axis=0),1,0)
    else:
        useful_sparse_matrix = np.where(t>0,1,0)
    t = None

    t1 = useful_sparse_matrix.transpose()*Pvalues
    P_values_array = t1.transpose()
    t1 = None 
    t2 = useful_sparse_matrix.sum(axis=0)
    t4 = np.where(np.all([P_values_array<P_cutoff, useful_sparse_matrix>0],axis=0),1,0)
    t4 = t4.sum(axis=0)
    testable_index= np.where(np.all([t2>=2,t4.astype(np.single)/t2>0.5],axis=0))[0]
    name_test_genes = np.asarray(matrix_columns)[testable_index]
 
    control_index= np.where(np.all([t2>=3,t4==0],axis=0))[0]
    name_control_genes = np.asarray(matrix_columns)[control_index]

    return name_test_genes,csr_matrix(useful_sparse_matrix[:,testable_index]), name_control_genes, csr_matrix(useful_sparse_matrix[:,control_index])

def filtergenes_sparse2(sitesAndgenes_sparse_matrix,
                matrix_columns,
                Pvalues,
                P_cutoff=0.05):
    
    shape = sitesAndgenes_sparse_matrix.shape
    assert shape[1] == len(matrix_columns) 
    useful_sparse_matrix = np.where(sitesAndgenes_sparse_matrix.toarray()>0,np.int8(1),np.int8(0))
    t1 = useful_sparse_matrix.transpose()*np.asarray(Pvalues,dtype=np.single)
    P_values_array = t1.transpose()
    t1 = None 
    t2 = useful_sparse_matrix.sum(axis=0)
    t4 = np.where(np.all([P_values_array<np.single(P_cutoff), useful_sparse_matrix>np.int8(0)],axis=0),np.int8(1),np.int8(0))
    t4 = t4.sum(axis=0)
    testable_index= np.where(np.all([t2>=np.int8(2),t4.astype(np.single)/t2 >np.single(0.5)],axis=0))[0]
    name_test_genes = np.asarray(matrix_columns)[testable_index]
    powerful_index= np.where(np.all([t2>=np.int8(2),t4.astype(np.single)/t2<=np.single(0.5)],axis=0))[0]
    ase_index = np.where(np.all([t2>=np.int8(2),t4/t2 ==1],axis=0))[0]
    part_ase_index=np.where(np.all([t2>=np.int8(2),t4.astype(np.single)/t2>np.single(0.5),t4.astype(np.single)/t2<np.single(1.0)],axis=0))[0]
    control_index= np.where(np.all([t2>=np.int8(3),t4==np.int8(0)],axis=0))[0]
    name_control_genes = np.asarray(matrix_columns)[control_index]

    return {'index_test_genes': testable_index,
            'index_control_genes':control_index,
            'index_powerful_genes':powerful_index,
            'index_ase_genes':ase_index,
            'index_part_ase_genes':part_ase_index}
