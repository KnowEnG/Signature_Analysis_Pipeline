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

def run_cos(run_parameters):
    """ Performs cosine similarity analysis and saves  results.

    Args:
        run_parameters: parameter set dictionary.
    """

    expression_name = run_parameters["spreadsheet_name_full_path"]
    signature_name  = run_parameters["signature_name_full_path"  ]

    expression_df   = kn.get_spreadsheet_df(expression_name)
    signature_df    = kn.get_spreadsheet_df(signature_name )

    samples_names    = expression_df.columns
    signatures_names =  signature_df.columns

    cos_mat = generate_cos_mat(expression_df, signature_df)
    cos_mat = map_cos_range(cos_mat, 0)
    cos_df  = pd.DataFrame(cos_mat, index=samples_names, columns=signatures_names)

    save_final_samples_signature(cos_df, run_parameters)


def run_cc_cos(run_parameters):
    """ wrapper: call sequence to perform signature analysis with
        consensus clustering and write results.

    Args:
        run_parameters: parameter set dictionary.
    """
    tmp_dir = 'tmp_cc_cos'
    run_parameters = update_tmp_directory(run_parameters, tmp_dir)

    expression_name = run_parameters["spreadsheet_name_full_path"]
    signature_name  = run_parameters["signature_name_full_path"  ]

    expression_df   = kn.get_spreadsheet_df(expression_name)
    signature_df    = kn.get_spreadsheet_df(signature_name )

    expression_mat = expression_df.as_matrix()
    signature_mat   =  signature_df.as_matrix()

    number_of_bootstraps = run_parameters['number_of_bootstraps']
    processing_method = run_parameters['processing_method']

    if processing_method == 'serial':
        for sample in range(0, number_of_bootstraps):
            run_cc_cos_signature_worker(expression_mat, signature_mat, run_parameters, sample)

    elif processing_method == 'parallel':
        find_and_save_cc_cos_signature_parallel(spreadsheet_mat, run_parameters, number_of_bootstraps)

    else:
        raise ValueError('processing_method contains bad value.')

    consensus_df = form_consensus_df(
        run_parameters, expression_df, signature_df)
    
    save_final_samples_signature(consensus_df, run_parameters)
    kn.remove_dir(run_parameters["tmp_directory"])

def run_net_cos(run_parameters):
    """ Run random walk first to smooth spreadsheet
    and perform cosine similarity analysis and saves results.

    Args:
        run_parameters: parameter set dictionary.
    """

    expression_name = run_parameters["spreadsheet_name_full_path"]
    signature_name  = run_parameters["signature_name_full_path"  ]
    gg_network_name = run_parameters['gg_network_name_full_path']

    expression_df   = kn.get_spreadsheet_df(expression_name)
    signature_df    = kn.get_spreadsheet_df(signature_name )

    samples_names    = expression_df.columns
    signatures_names =  signature_df.columns

    network_mat, unique_gene_names = kn.get_sparse_network_matrix(gg_network_name)
    network_mat = kn.normalize_sparse_mat_by_diagonal(network_mat)
    
    expression_df = kn.update_spreadsheet_df(expression_df, unique_gene_names)
    signature_df = kn.update_spreadsheet_df(signature_df, unique_gene_names)

    expression_mat = expression_df.as_matrix()
    signature_mat = signature_df.as_matrix()

    expression_mat, iterations = kn.smooth_matrix_with_rwr(expression_mat, network_mat, run_parameters)
    signature_mat, iterations = kn.smooth_matrix_with_rwr(signature_mat, network_mat, run_parameters)

    expression_df.iloc[:] = expression_mat
    signature_df.iloc[:] = signature_mat

    cos_mat = generate_cos_mat(expression_df, signature_df)
    cos_mat = map_cos_range(cos_mat, 0)
    cos_df  = pd.DataFrame(cos_mat, index=samples_names, columns=signatures_names)

    save_final_samples_signature(cos_df, run_parameters)

def run_cc_net_cos(run_parameters):
    """ wrapper: call sequence to perform signature analysis with
        random walk smoothing and consensus clustering and write results.

    Args:
        run_parameters: parameter set dictionary.
    """
    tmp_dir = 'tmp_cc_cos'
    run_parameters = update_tmp_directory(run_parameters, tmp_dir)

    expression_name = run_parameters["spreadsheet_name_full_path"]
    signature_name  = run_parameters["signature_name_full_path"  ]
    gg_network_name = run_parameters['gg_network_name_full_path']

    expression_df   = kn.get_spreadsheet_df(expression_name)
    signature_df    = kn.get_spreadsheet_df(signature_name )

    samples_names    = expression_df.columns
    signatures_names =  signature_df.columns

    network_mat, unique_gene_names = kn.get_sparse_network_matrix(gg_network_name)
    network_mat = kn.normalize_sparse_mat_by_diagonal(network_mat)
    
    expression_df = kn.update_spreadsheet_df(expression_df, unique_gene_names)
    signature_df = kn.update_spreadsheet_df(signature_df, unique_gene_names)

    expression_mat = expression_df.as_matrix()
    signature_mat = signature_df.as_matrix()

    expression_mat, iterations = kn.smooth_matrix_with_rwr(expression_mat, network_mat, run_parameters)
    signature_mat, iterations = kn.smooth_matrix_with_rwr(signature_mat, network_mat, run_parameters)

    expression_df.iloc[:] = expression_mat
    signature_df.iloc[:] = signature_mat

    number_of_bootstraps = run_parameters['number_of_bootstraps']
    processing_method = run_parameters['processing_method']

    if processing_method == 'serial':
        for sample in range(0, number_of_bootstraps):
            run_cc_cos_signature_worker(expression_mat, signature_mat, run_parameters, sample)

    elif processing_method == 'parallel':
        find_and_save_cc_cos_signature_parallel(spreadsheet_mat, run_parameters, number_of_bootstraps)

    else:
        raise ValueError('processing_method contains bad value.')

    consensus_df = form_consensus_df(
        run_parameters, expression_df, signature_df)
    
    save_final_samples_signature(consensus_df, run_parameters)
    kn.remove_dir(run_parameters["tmp_directory"])


def find_and_save_cc_cos_signature_parallel(spreadsheet_mat, run_parameters, local_parallelism):
    """ central loop: compute components for the consensus matrix by
        non-negative matrix factorization.

    Args:
        spreadsheet_mat: genes x samples matrix.
        run_parameters: dictionary of run-time parameters.
        local_parallelism: parallelism option
    """
    import knpackage.distributed_computing_utils as dstutil

    jobs_id = range(0, local_parallelism)
    zipped_arguments = dstutil.zip_parameters(spreadsheet_mat, run_parameters, jobs_id)
    if 'parallelism' in run_parameters:
        parallelism = dstutil.determine_parallelism_locally(local_parallelism, run_parameters['parallelism'])
    else:
        parallelism = dstutil.determine_parallelism_locally(local_parallelism)
    dstutil.parallelize_processes_locally(run_cc_cos_signature_worker, zipped_arguments, parallelism)


def run_cc_cos_signature_worker(expression_mat, signature_mat, run_parameters, sample):
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

    np.random.seed(sample)
    rows_sampling_fraction = run_parameters["rows_sampling_fraction"]
    cols_sampling_fraction = run_parameters["cols_sampling_fraction"]
    expression_mat_T, sample_permutation_e = kn.sample_a_matrix(expression_mat.T,
                                                             rows_sampling_fraction, cols_sampling_fraction)
    
    signature_mat_T, sample_permutation_s = kn.sample_a_matrix(signature_mat.T,
                                                             rows_sampling_fraction, cols_sampling_fraction)

    save_a_signature_to_tmp(expression_mat_T.T, sample_permutation_e, signature_mat_T.T, sample_permutation_s, run_parameters, sample)


def save_a_signature_to_tmp(expression_mat, sample_permutation_e, signature_mat, sample_permutation_s, run_parameters, sequence_number):
    """ save one h_matrix and one permutation in temorary files with sequence_number appended names.

    Args:
        expression_mat: permutation x k size matrix.
        sample_permutation_e: indices of expression_mat rows permutation.
        signature_mat: permutation x k size matrix.
        sample_permutation_s: indices of signature_mat rows permutation.
        run_parameters: parmaeters including the "tmp_directory" name.
        sequence_number: temporary file name suffix.
    """
    import os
    import numpy as np

    tmp_dir = run_parameters["tmp_directory"]

    os.makedirs(tmp_dir, mode=0o755, exist_ok=True)

    hname_e = os.path.join(tmp_dir, 'tmp_h_e_%d'%(sequence_number))
    pname_e = os.path.join(tmp_dir, 'tmp_p_e_%d'%(sequence_number))

    hname_s = os.path.join(tmp_dir, 'tmp_h_s_%d'%(sequence_number))
    pname_s = os.path.join(tmp_dir, 'tmp_p_s_%d'%(sequence_number))

    with open(hname_e, 'wb') as fh0:
        expression_mat.dump(fh0)
    with open(pname_e, 'wb') as fh1:
        sample_permutation_e.dump(fh1)
    with open(hname_s, 'wb') as fh2:
        signature_mat.dump(fh2)
    with open(pname_s, 'wb') as fh3:
        sample_permutation_s.dump(fh3)


def form_consensus_df(run_parameters, expression_df_orig, signature_df_orig):
    """ compute the consensus df from the express dataframe and signature dataframe
        formed by the bootstrap "temp_*" files.

    Args:
        run_parameters: parameter set dictionary with "tmp_directory" key.
        expression_df_orig: dataframe of expression data.
        signature_df_orig: dataframe of signature data.

    Returns:
        cos_df: cos_df with the value to be consensus matrix
    """

    if run_parameters['processing_method'] == 'distribute':
        tmp_dir = os.path.join(run_parameters['cluster_shared_volumn'],
                               os.path.basename(os.path.normpath(run_parameters['tmp_directory'])))
    else:
        tmp_dir = run_parameters["tmp_directory"]
        
    dir_list = os.listdir(tmp_dir)
    samples_names    = expression_df_orig.columns
    signatures_names =  signature_df_orig.columns
    cos_array = np.zeros((samples_names.shape[0], signatures_names.shape[0]))

    for tmp_f in dir_list:
        if tmp_f[0:8] == 'tmp_p_e_':
            pname_e = os.path.join(tmp_dir, tmp_f)
            hname_e = os.path.join(tmp_dir, 'tmp_h_e_' + tmp_f[8:len(tmp_f)])
            pname_s = os.path.join(tmp_dir, 'tmp_p_s_' + tmp_f[8:len(tmp_f)])
            hname_s = os.path.join(tmp_dir, 'tmp_h_s_' + tmp_f[8:len(tmp_f)])

            expression_mat = np.load(hname_e)
            signature_mat  = np.load(hname_s)
            sample_permutation_e = np.load(pname_e)
            sample_permutation_s = np.load(pname_s)

            expression_df = pd.DataFrame(expression_mat)
            expression_df.index = expression_df_orig.index[sample_permutation_e]
            expression_df.columns = expression_df_orig.columns
            
            signature_df = pd.DataFrame(signature_mat)
            signature_df.index = signature_df_orig.index[sample_permutation_s]
            signature_df.columns = signature_df_orig.columns

            cos_mat = generate_cos_mat(expression_df, signature_df)
            cos_array += cos_mat

    cos_array /= run_parameters['number_of_bootstraps']
    cos_array = map_cos_range(cos_array, 0)
    cos_df = pd.DataFrame(cos_array, index=samples_names, columns=signatures_names)
    return cos_df

def generate_cos_mat(expression_df, signature_df):
    """generate matrix which save the cosine value of input dataframes

    Args:
        expression_df: genes x samples dataframe.
        signature_df: genes x samples dataframe.
        
    Returns:
        cos_mat: matrix with cosine value
    """
    genes_in_expression =  expression_df.index
    genes_in_signature  =   signature_df.index

    common_genes        = kn.find_common_node_names(genes_in_expression, genes_in_signature)
    expression_mat      = expression_df.loc[common_genes, :].values
    signature_mat       =  signature_df.loc[common_genes, :].values
    cos_mat             = cosine_similarity(expression_mat.T, signature_mat.T)
    return cos_mat

def map_cos_range(cos_mat, axis_val):
    """Normalize cosine matrix via given axis

    Args:
        cos_mat: sample1 x sample2 matrix.
        axis_val: given axis.
        
    Returns:
        cos_mat: normalized cosine matrix with 0 and 1.
    """
    max_value_row_index = np.argmax(cos_mat, axis=axis_val)
    num_of_cols =  len(cos_mat[0])

    cos_mat[max_value_row_index, range(num_of_cols)] = 1
    cos_mat[cos_mat!=1] = 0

    return cos_mat

def save_final_samples_signature(result_df, run_parameters):
    """ wtite .tsv file that assings a cluster number label to the sample_names.

    Args:
        result_df: result dataframe
        run_parameters: write path (run_parameters["results_directory"]).
    """
    result_df.to_csv(get_output_file_name(run_parameters, 'result', 'viz'), sep='\t')

def get_output_file_name(run_parameters, prefix_string, suffix_string='', type_suffix='tsv'):
    """ get the full directory / filename for writing
    Args:
        run_parameters: dictionary with keys: "results_directory", "method" and "correlation_measure"
        prefix_string:  the first letters of the ouput file name
        suffix_string:  the last letters of the output file name before '.tsv'

    Returns:
        output_file_name:   full file and directory name suitable for file writing
    """
    output_file_name = os.path.join(run_parameters["results_directory"], prefix_string + '_' + run_parameters['method'])
    output_file_name = kn.create_timestamped_filename(output_file_name) + '_' + suffix_string + '.' + type_suffix

    return output_file_name


def update_tmp_directory(run_parameters, tmp_dir):
    ''' Update tmp_directory value in rum_parameters dictionary

    Args:
        run_parameters: run_parameters as the dictionary config
        tmp_dir: temporary directory prefix subjected to different functions

    Returns:
        run_parameters: an updated run_parameters

    '''
    if (run_parameters['processing_method'] == 'distribute'):
        run_parameters["tmp_directory"] = kn.create_dir(run_parameters['cluster_shared_volumn'], tmp_dir)
    else:
        run_parameters["tmp_directory"] = kn.create_dir(run_parameters["run_directory"], tmp_dir)

    return run_parameters

