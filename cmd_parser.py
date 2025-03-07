import argparse
import sys

# Format partition
def parse_partition(value):
    parts = value.split(',')

    if len(parts) != 3:
        first_part = [int(x) for x in parts[0].split('+')]
        second_part = int(parts[1])
        return [first_part, second_part]

    else:
        return [int(x) for x in parts]


#"skygrid_populationsize": '32', #integer
         #"skygrid_gridpoints": '31.0',
         #'skygrid_cutoff':'8.0',
         #'clock_model': 'ucld',
         #'ucld_mean': '0.003',
         #'ucld_stdev': '0.05',
        # 'ucld_meaninrealspace': 'true',
        # 'substitution_model': 'hky',
         #'hky_kappa': '2.0',
         #'gamma_categories': '4',
         #'gamma_alpha': '0.5',
         #'population_model': 'constant',
         #'partitions': [[1,2],3],
        #'empirical_tree_distribution': None,
       # 'chain_length':'10000000',
        #'log_every':'1000',
        #"gtr_rates_value": '1.0',
        #"gtr_rates_dimension": '6',
        #'continuous_phylogeo': 'true',
       # 'continuous_phylogeo_jitter': '0.001',
        #'trait_names': ['lat', 'lon'],

def parse_args():
    parser = argparse.ArgumentParser(description='make a BEAST v1.10.4 XML')

    #Sequence parameters
    seq_group = parser.add_argument_group('Sequence options')
    seq_group.add_argument("--fasta",
                           dest="fasta",
                           required=True,
                           help="The fasta file for analysis")

    seq_group.add_argument("--codon-partitioning",
                           dest="partitions",
                           type=parse_partition,
                           #action="store_false",
                           help="comma separated list describing codon partitioning. For example, the codon partitioning required for SRD06 would be 1+2,3")

    # Substitution model parameters
    site_group = parser.add_argument_group('Substitution model options')

    site_group.add_argument("--substitution-model",
                             dest="substitution_model",
                             choices=['hky', 'gtr'],
                             help="Select substitution model. Currently limited to HKY and GTR")

    site_group.add_argument("--gamma",
                            dest="use_gamma",
                            action="store_true",
                            help='Select the number of discretised gamma categories to model rate variation between sites.')

    site_group.add_argument("--gamma-categories",
                            dest="gamma_categories",
                            default='4',
                            help='Select the number of discretised gamma categories to model rate variation between sites.')

    # Clock parameters
    clock_group = parser.add_argument_group('Molecular clock model options')

    clock_group.add_argument("--clock-model",
                            dest="clock_model",
                            choices=['ucld', 'strict'],
                            help="Select clock model. Currently limited to strict and uncorrelated relaxed lognormal")

    # Tree parameters
    tree_group = parser.add_argument_group('Tree model options')
    tree_group.add_argument("--tree-model",
                            dest="tree_model",
                            choices=['constant', 'skygrid'],
                            help="Select tree model. Currently limited to constant and skygrid")

    tree_group.add_argument("--skygrid-cutoff",
                            help='Time at last transition point',
                            dest="skygrid_cutoff")

    tree_group.add_argument("--skygrid-grids",
                            default='50',
                            help="Number of parameters",
                            dest="skygrid_grids")

    tree_group.add_argument("--empirical-tree-model",
                            help="Flag to run an empirical tree distribution model",
                            default=False,
                            action="store_true",
                            dest='empirical_tree_model')

    tree_group.add_argument("--empirical-tree-distribution",
                            dest="empirical_tree_distribution",
                            required='--empirical-tree-model' in sys.argv,
                            help="Treefile containing posterior tree distribution")


    # Relaxed clock priors
    prior_group = parser.add_argument_group("Prior Options")
    prior_group.add_argument("--ucld-mean",  # change to start values (check this is correct)!!!
                             dest="ucld_mean",
                             default='lognormal,0.001,0.001,true',
                             help="Specify the prior distribution for uncorrelated lognormal relaxed clock mean. See ReadMe for details of how to declare prior distributions.")

    prior_group.add_argument("--ucld-stdev",
                             dest="ucld_stdev",
                             default='exponential,0.3333333333333333',
                             help="Specify the prior distribution for uncorrelated lognormal relaxed clock standard deviation")

    # Strict clock priors
    prior_group.add_argument("--clock-rate",  # change to start values (check this is correct)!!!
                             dest="clock_rate",
                             default='lognormal,0.001,0.001,true',
                             help="Specify the prior distribution for strict clock substitution rate")

    # HKY model priors
    prior_group.add_argument("--hky-kappa",
                             dest="hky_kappa",
                             default='lognormal,1.0,1.25,false',
                             help='Specify the prior distribution for HKY transition-transversion parameter')

    # GTR model priors
    site_group.add_argument("--gtr-ac",
                            dest="gtr_ac",
                            default= 'gamma,0.05,20.0',
                            help='Specify the prior distribution for GTR A-C substitution parameter')

    site_group.add_argument("--gtr-ag",
                            dest="gtr_ag",
                            default='gamma,0.05,20.0',
                            help='Specify the prior distribution for GTR A-G substitution parameter')

    site_group.add_argument("--gtr-at",
                            dest="gtr_at",
                            default='gamma,0.05,20.0',
                            help='Specify the prior distribution for GTR A-T substitution parameter')

    site_group.add_argument("--gtr-cg",
                            dest="gtr_cg",
                            default='gamma,0.05,20.0',
                            help='Specify the prior distribution for GTR C-G substitution parameter')

    site_group.add_argument("--gtr-gt",
                            dest="gtr_gt",
                            default='gamma,0.05,20.0',
                            help='Specify the prior distribution for GTR G-T substitution parameter')

    # Gamma
    prior_group.add_argument("--gamma-alpha",
                             dest="gamma_alpha",
                             default='exponential,0.5',
                             help='Select alpha parameter of gamma distribution to model rate variation between sites.')

    # Relative rate prior (HKY)
    prior_group.add_argument("--allMus",
                             dest="all_mus",
                             default='uniform,0,100',
                             help='Specify the prior distribution for relative rates amongst partitions parameter' )

    # Constant Population prior
    prior_group.add_argument("--constant-pop",
                             dest="constant_population",
                             default='lognormal,10,100,true',
                             help='Specify the prior distribution for coalescent population size parameter')

    # Traits
    trait_group = parser.add_argument_group("Trait Options")
    trait_group.add_argument("--continuous-phylogeo",
                             action="store_true",
                             dest="continuous_phylogeo",
                             help="Flag to run a continuous phylogeographic analysis")

    trait_group.add_argument("--continuous-phylogeo-coords",
                             dest="continuous_trait_file",
                             required='--continuous-phylogeo' in sys.argv,
                             help="Comma delimited file with headers 'taxon,lat,lon' for continuous phylogeographic analysis")

    trait_group.add_argument("--continuous-phylogeo-jitter",
                             dest="continuous_phylogeo_jitter",
                             default="0.001",
                             help="Jitter for duplicate points in continuous phylogeographic analysis")

    # General settings for MCMC run
    gen_group = parser.add_argument_group('General options')
    gen_group.add_argument("--chain-length",
                           dest="chain_length",
                           default="100000000",
                           help="Number of states to run the MCMC chain for. Default=100m")

    gen_group.add_argument("--log-every",
                           dest="log_every",
                           default="10000",
                           help="How often to log tree and log files. Default=10,000")

    gen_group.add_argument("--file-stem",
                           dest="file_stem",
                           help="File stem for analysis")

    args = parser.parse_args()
    return args