"""
Created on Fri Jul 14 11:43:45 2017
@author: The KnowEnG dev team
"""

def cos(run_parameters):
    '''cos signature'''
    from gene_signature_toolbox import run_cos
    run_cos(run_parameters) 

def cc_cos(run_parameters):
    '''kmeans consensus signature of the cos-based similarity'''
    from gene_signature_toolbox import run_cc_cos
    run_cc_cos(run_parameters)

def net_cos(run_parameters):
    '''net-cos signature "'''
    from gene_signature_toolbox import run_net_cos
    run_net_cos(run_parameters)

def cc_net_cos(run_parameters):
    '''kmeans consensus signature of the net-cos-based similarity'''
    from gene_signature_toolbox import run_cc_net_cos
    run_cc_net_cos(run_parameters)

SELECT = {
    "cos":cos,
    "cc_cos":cc_cos,
    "net_cos":net_cos,
    "cc_net_cos":cc_net_cos}

def main():
    """
    This is the main function to perform sample similarity
    """
    import sys
    from knpackage.toolbox import get_run_directory_and_file
    from knpackage.toolbox import get_run_parameters
    
    run_directory, run_file = get_run_directory_and_file(sys.argv)
    run_parameters = get_run_parameters(run_directory, run_file)
    SELECT[run_parameters["method"]](run_parameters)

if __name__ == "__main__":
    main()
