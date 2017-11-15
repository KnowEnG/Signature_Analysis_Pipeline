"""
@author: The KnowEnG dev team
"""
import os
import numpy as np
import pandas as pd
import knpackage.toolbox as kn
import knpackage.distributed_computing_utils as dstutil

from   sklearn.metrics import silhouette_score
from   sklearn.metrics.pairwise import cosine_similarity
from   scipy.stats              import spearmanr


def run_similarity(run_parameters):
    """ Performs similarity analysis and saves the similarity matrix.

    Args:
        run_parameters: parameter set dictionary.
    """

    expression_name     = run_parameters["spreadsheet_name_full_path"]
    signature_name      = run_parameters["signature_name_full_path"  ]
    similarity_measure  = run_parameters["similarity_measure"        ]

    expression_df       = kn.get_spreadsheet_df(expression_name)
    signature_df        = kn.get_spreadsheet_df(signature_name )
    
    samples_names       = expression_df.columns
    signatures_names    =  signature_df.columns
    signatures_names    = [i.split('.')[0] for i in signatures_names]
    signature_df.columns= signatures_names

    similarity_mat = generate_similarity_mat(expression_df, signature_df,similarity_measure)
    # similarity_mat = map_similarity_range(similarity_mat, 0)
    similarity_df  = pd.DataFrame(similarity_mat, index=samples_names, columns=signatures_names)
    if(similarity_measure == "spearman"):
        accuracy = calculate_accuracy(similarity_df)
        print(accuracy)
    save_final_samples_signature(similarity_df, run_parameters)


def run_cc_similarity(run_parameters):
    """ Performs similarity analysis with bootstraps and saves the similarity matrix.

    Args:
        run_parameters: parameter set dictionary.
    """
    tmp_dir        = 'tmp_cc_similarity'
    run_parameters = update_tmp_directory(run_parameters, tmp_dir)

    expression_name      = run_parameters["spreadsheet_name_full_path"]
    signature_name       = run_parameters["signature_name_full_path"  ]
    similarity_measure  = run_parameters["similarity_measure"         ]
    number_of_bootstraps = run_parameters['number_of_bootstraps'      ]
    processing_method    = run_parameters['processing_method'         ]

    expression_df        = kn.get_spreadsheet_df(expression_name)
    signature_df         = kn.get_spreadsheet_df(signature_name )

    samples_names       = expression_df.columns
    signatures_names    =  signature_df.columns
    signatures_names    = [i.split('.')[0] for i in signatures_names]
    signature_df.columns= signatures_names

    expression_mat      = expression_df.as_matrix()
    signature_mat       =  signature_df.as_matrix()
    if   processing_method == 'serial':
         for sample in range(0, number_of_bootstraps):
           run_cc_similarity_signature_worker(expression_df, signature_df, run_parameters, sample)

    elif processing_method == 'parallel':
         find_and_save_cc_similarity_parallel(expression_df, signature_df, run_parameters, number_of_bootstraps)

    else:
        raise ValueError('processing_method contains bad value.')

    # consensus_df = form_consensus_df(run_parameters, expression_df, signature_df)
    similarity_df = assemble_similarity_df(expression_df, signature_df, run_parameters)

    similarity_df  = pd.DataFrame(similarity_df.values, index=samples_names, columns=signatures_names)
    if(similarity_measure == "spearman"):
        accuracy = calculate_accuracy(similarity_df)
        print(accuracy)
    save_final_samples_signature(similarity_df, run_parameters)

    kn.remove_dir(run_parameters["tmp_directory"])


def run_net_similarity(run_parameters):
    """ Run random walk first to smooth expression and signature 
    then perform similarity analysis and save the similarity matrix.

    Args:
        run_parameters: parameter set dictionary.
    """
    expression_name     = run_parameters["spreadsheet_name_full_path"]
    signature_name      = run_parameters["signature_name_full_path"  ]
    gg_network_name     = run_parameters['gg_network_name_full_path' ]
    similarity_measure  = run_parameters["similarity_measure"        ]

    expression_df       = kn.get_spreadsheet_df(expression_name)
    signature_df        = kn.get_spreadsheet_df(signature_name )

    samples_names       = expression_df.columns
    signatures_names    =  signature_df.columns
    signatures_names    = [i.split('.')[0] for i in signatures_names]
    signature_df.columns= signatures_names

    network_mat, unique_gene_names = kn.get_sparse_network_matrix(gg_network_name)
    # network_mat                    = kn.normalize_sparse_mat_by_diagonal(network_mat)
    
    expression_df                  = kn.update_spreadsheet_df(expression_df, unique_gene_names)
    signature_df                   = kn.update_spreadsheet_df(signature_df, unique_gene_names)

    expression_mat                 = expression_df.as_matrix()
    signature_mat                  = signature_df.as_matrix()

    expression_mat, iterations = kn.smooth_matrix_with_rwr(expression_mat, network_mat, run_parameters)
    signature_mat,  iterations = kn.smooth_matrix_with_rwr(signature_mat,  network_mat, run_parameters)

    expression_df.iloc[:] = expression_mat
    signature_df.iloc[:]  = signature_mat

    similarity_mat = generate_similarity_mat(expression_df, signature_df,similarity_measure)
    # similarity_mat = map_similarity_range(similarity_mat, 0)
    similarity_df  = pd.DataFrame(similarity_mat, index=samples_names, columns=signatures_names)
    if(similarity_measure == "spearman"):
        accuracy = calculate_accuracy(similarity_df)
        print(accuracy)

    save_final_samples_signature(similarity_df, run_parameters)


def run_cc_net_similarity(run_parameters):
    """ wrapper: call sequence to perform signature analysis with
        random walk smoothing and bootstrapped similarity and save results.

    Args:
        run_parameters: parameter set dictionary.
    """
    tmp_dir = 'tmp_cc_similarity_'
    run_parameters = update_tmp_directory(run_parameters, tmp_dir)

    expression_name      = run_parameters["spreadsheet_name_full_path"]
    signature_name       = run_parameters["signature_name_full_path"  ]
    gg_network_name      = run_parameters['gg_network_name_full_path' ]
    similarity_measure   = run_parameters["similarity_measure"        ]
    number_of_bootstraps = run_parameters['number_of_bootstraps'      ]
    processing_method    = run_parameters['processing_method'         ]

    expression_df        = kn.get_spreadsheet_df(expression_name)
    signature_df         = kn.get_spreadsheet_df(signature_name )

    samples_names        = expression_df.columns
    signatures_names     =  signature_df.columns
    signatures_names     = [i.split('.')[0] for i in signatures_names]
    signature_df.columns = signatures_names

    network_mat, unique_gene_names = kn.get_sparse_network_matrix(gg_network_name)
    # network_mat                    = kn.normalize_sparse_mat_by_diagonal(network_mat)
    
    expression_df                  = kn.update_spreadsheet_df(expression_df, unique_gene_names)
    signature_df                   = kn.update_spreadsheet_df(signature_df, unique_gene_names)

    expression_mat                 = expression_df.as_matrix()
    signature_mat                  = signature_df.as_matrix()

    expression_mat, iterations = kn.smooth_matrix_with_rwr(expression_mat, network_mat, run_parameters)
    signature_mat,  iterations = kn.smooth_matrix_with_rwr(signature_mat,  network_mat, run_parameters)

    expression_df.iloc[:] = expression_mat
    signature_df.iloc[:]  = signature_mat

    if   processing_method == 'serial':
         for sample in range(0, number_of_bootstraps):
            run_cc_similarity_signature_worker(expression_df, signature_df, run_parameters, sample)

    elif processing_method == 'parallel':
         find_and_save_cc_similarity_parallel(expression_df, signature_df, run_parameters, number_of_bootstraps)

    else:
        raise ValueError('processing_method contains bad value.')

    # consensus_df = form_consensus_df(run_parameters, expression_df, signature_df)
    similarity_df = assemble_similarity_df(expression_df, signature_df, run_parameters)
    similarity_df  = pd.DataFrame(similarity_df.values, index=samples_names, columns=signatures_names)
    if(similarity_measure == "spearman"):
        accuracy = calculate_accuracy(similarity_df)
        print(accuracy)
    save_final_samples_signature(similarity_df, run_parameters)
    kn.remove_dir(run_parameters["tmp_directory"])


def find_and_save_cc_similarity_parallel(expression_df, signature_df, run_parameters, local_parallelism):
    """ central loop: compute components for the similarity matrix by

    Args:
        expression_df    : genes x samples
        signature_df     : genes x samples
        run_parameters   : dictionary of run-time parameters
        local_parallelism: parallelism option
    """
    import knpackage.distributed_computing_utils as dstutil

    jobs_id          = range(0, local_parallelism)
    zipped_arguments = dstutil.zip_parameters(expression_df, signature_df, run_parameters, jobs_id)

    if 'parallelism' in run_parameters:
        parallelism = dstutil.determine_parallelism_locally(local_parallelism, run_parameters['parallelism'])

    else:
        parallelism = dstutil.determine_parallelism_locally(local_parallelism)

    dstutil.parallelize_processes_locally(run_cc_similarity_signature_worker, zipped_arguments, parallelism)


def run_cc_similarity_signature_worker(expression_df, signature_df, run_parameters, sample):
    """Worker to execute nmf_clusters in a single process

    Args:
        expression_mat: genes x samples matrix.
        signature_mat: genes x samples matrix.
        run_parameters: dictionary of run-time parameters.
        sample: each loops.

    Returns:
        None

    """
    import knpackage.toolbox as kn
    import numpy as np

    rows_sampling_fraction = run_parameters["rows_sampling_fraction"]
    similarity_measure     = run_parameters['similarity_measure'   ]
    sampled_expression_df  = expression_df.sample(frac = rows_sampling_fraction,random_state=sample)
    sampled_signature_df   =  signature_df.sample(frac = rows_sampling_fraction,random_state=sample)
    sampled_similarity_mat = generate_similarity_mat(sampled_expression_df, sampled_signature_df, similarity_measure)

    save_a_signature_to_tmp(sampled_similarity_mat, run_parameters, sample)
    
def save_a_signature_to_tmp(sampled_similarity_mat, run_parameters, sequence_number):
    """ save a sampled_similarity_mat in temorary files with sequence_number appended names.
    Args:
        run_parameters: parmaeters including the "tmp_directory" name.
        sequence_number: temporary file name suffix.
    """
    tmp_dir = run_parameters["tmp_directory"]

    os.makedirs(tmp_dir, mode=0o755, exist_ok=True)

    hname_e = os.path.join(tmp_dir, 'tmp_h_e_%d'%(sequence_number))

    with open(hname_e, 'wb') as fh0:
        sampled_similarity_mat.dump(fh0)

def assemble_similarity_df(expression_df, signature_df, run_parameters):
    """ compute the similarity df from the express dataframe and signature dataframe
        formed by the bootstrap "temp_*" files.
    Args:
        run_parameters: parameter set dictionary with "tmp_directory" key.
        expression_df: dataframe of expression data.
        signature_df: dataframe of signature data.
    Returns:
        similarity_df: similarity_df with the value to be similarity matrix
    """

    processing_method     = run_parameters['processing_method'    ]
    tmp_directory         = run_parameters['tmp_directory'        ]
    cluster_shared_volumn = run_parameters['cluster_shared_volumn']
    number_of_bootstraps  = run_parameters['number_of_bootstraps' ]

    if processing_method == 'distribute':
        tmp_dir = os.path.join(cluster_shared_volumn,
                               os.path.basename(os.path.normpath(tmp_directory)))
    else:
        tmp_dir = tmp_directory
        
    dir_list         = os.listdir(tmp_dir)
    expression_names = expression_df.columns
    signatures_names =  signature_df.columns
    similarity_mat   = np.zeros((expression_names.shape[0], signatures_names.shape[0]))

    for tmp_f in dir_list:
        if tmp_f[0:8] == 'tmp_h_e_':
            hname_e = os.path.join(tmp_dir, 'tmp_h_e_' + tmp_f[8:len(tmp_f)])

            sampled_similarity_mat = np.load(hname_e)

            similarity_mat += sampled_similarity_mat

    similarity_mat /= number_of_bootstraps

    # similarity_mat  = map_similarity_range(similarity_mat, 0)
    similarity_df   = pd.DataFrame(similarity_mat, index=expression_names, columns=signatures_names)

    return similarity_df



def generate_similarity_mat(expression_df, signature_df,similarity_measure):
    """generate matrix which save the similarity value of input dataframes

    Args:
        expression_df: genes x samples dataframe.
        signature_df:  genes x samples dataframe.
        
    Returns:
        similarity_mat: matrix with similarity values
    """

    genes_in_expression =  expression_df.index
    genes_in_signature  =   signature_df.index

    common_genes        = kn.find_common_node_names(genes_in_expression, genes_in_signature)
    expression_mat      = expression_df.loc[common_genes, :].values
    signature_mat       =  signature_df.loc[common_genes, :].values
    nx                  = expression_mat.shape[1]

    if   (similarity_measure == "cosine" ):
          similarity_mat      = cosine_similarity(expression_mat.T, signature_mat.T)
          # print(similarity_mat.shape)
    elif (similarity_measure == "spearman"):
          similarity_mat      = spearmanr(expression_mat, signature_mat)[0]
          # print(expression_mat)
          similarity_mat      = np.abs(similarity_mat[0:nx,nx:] )
          # print(similarity_mat.shape)
    return similarity_mat


def save_final_samples_signature(result_df, run_parameters):
    """ wtite .tsv file that assings a cluster number label to the sample_names.

    Args:
        result_df: result dataframe
        run_parameters: write path (run_parameters["results_directory"]).
    """
    fn_result = get_output_file_name(run_parameters, 'result', 'viz')
    result_df.to_csv(fn_result, sep='\t', float_format='%g')


def get_output_file_name(run_parameters, prefix_string, suffix_string='', type_suffix='tsv'):
    """ get the full directory / filename for writing
    Args:
        run_parameters: dictionary with keys: "results_directory", "method" and "correlation_measure"
        prefix_string:  the first letters of the ouput file name
        suffix_string:  the last letters of the output file name before '.tsv'

    Returns:
        output_file_name:   full file and directory name suitable for file writing
    """
    results_directory  = run_parameters["results_directory" ]
    method             = run_parameters['method'            ]
    similarity_measure = run_parameters['similarity_measure']

    output_file_name   = os.path.join(results_directory, prefix_string + '_' + method + '_' + similarity_measure)
    output_file_name   = kn.create_timestamped_filename(output_file_name) + '_' + suffix_string + '.' + type_suffix

    return output_file_name


def update_tmp_directory(run_parameters, tmp_dir):
    ''' Update tmp_directory value in rum_parameters dictionary

    Args:
        run_parameters: run_parameters as the dictionary config
        tmp_dir: temporary directory prefix subjected to different functions

    Returns:
        run_parameters: an updated run_parameters

    '''
    processing_method     = run_parameters['processing_method'    ]
    cluster_shared_volumn = run_parameters['cluster_shared_volumn']
    run_directory         = run_parameters["run_directory"        ]

    if processing_method == 'distribute':
        run_parameters["tmp_directory"] = kn.create_dir(cluster_shared_volumn, tmp_dir)
    else:
        run_parameters["tmp_directory"] = kn.create_dir(run_directory, tmp_dir)

    return run_parameters

def calculate_accuracy(similarity_df):
    """ Calculate accuracy given similarity dataframe and benchmark result

    Args:
        similarity_df: result dataframe from run_similarity methods
    """
    result = similarity_df.idxmax(axis=1, skipna=True)
    benchmark = pd.read_csv('../data/spreadsheets/label_validation.txt', index_col=None, header=None, sep='\t')
    ret_li = result.values
    ben_li = benchmark.values.reshape((1,-1))[0]

    common = ret_li==ben_li
    common[common==True] = 1
    common[common==False] = 0
    accuracy = sum(common)/len(ret_li)
    return accuracy
