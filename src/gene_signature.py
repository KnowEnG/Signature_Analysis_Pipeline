"""
Created on Fri Jul 14 11:43:45 2017
@author: The KnowEnG dev team
"""

def similarity(run_parameters):
    '''similarity signature'''
    from gene_signature_toolbox import run_similarity
    run_similarity(run_parameters) 

def cc_similarity(run_parameters):
    '''kmeans consensus signature of the similarity-based similarity'''
    from gene_signature_toolbox import run_cc_similarity
    run_cc_similarity(run_parameters)

def net_similarity(run_parameters):
    '''net-similarity signature "'''
    from gene_signature_toolbox import run_net_similarity
    run_net_similarity(run_parameters)

def cc_net_similarity(run_parameters):
    '''kmeans consensus signature of the net-similarity-based similarity'''
    from gene_signature_toolbox import run_cc_net_similarity
    run_cc_net_similarity(run_parameters)

SELECT = {
    "similarity":similarity,
    "cc_similarity":cc_similarity,
    "net_similarity":net_similarity,
    "cc_net_similarity":cc_net_similarity}

def main():
    """
    This is the main function to perform sample similarity
    """
    import sys
    from knpackage.toolbox import get_run_directory_and_file
    from knpackage.toolbox import get_run_parameters
    
    run_directory, run_file = get_run_directory_and_file(sys.argv)
    run_parameters          = get_run_parameters(run_directory, run_file)

    SELECT[run_parameters["method"]](run_parameters)

if __name__ == "__main__":
    main()
