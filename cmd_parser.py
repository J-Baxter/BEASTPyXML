from WriteEddieShell.BEASTPyXML.block_functions import *
from WriteEddieShell.BEASTPyXML.operators import *
from WriteEddieShell.BEASTPyXML.mcmc import *


filename = 'USUV_2025Feb10_alldata_aligned_formatted_noFLI_NFLG.fasta_subsample1.fasta'

args = {"skygrid_populationsize": '32', #integer
         "skygrid_gridpoints": '31.0',
         'skygrid_cutoff':'8.0',
         'clock_model': 'ucld',
         'ucld_mean': '0.003',
         'ucld_stdev': '0.05',
         'ucld_meaninrealspace': 'true',
         'substitution_model': 'hky',
         'hky_kappa': '2.0',
         'gamma_categories': '4',
         'gamma_alpha': '0.5',
         'population_model': 'constant',
         'partitions': [[1,2],3],
        'empirical_tree_distribution': None,
        'chain_length':'250000000',
        'log_every':'25000',
        'file_stem':'USUV_2025Feb10_alldata_aligned_formatted_noFLI_NFLG'}



# Make base root of XML
root = etree.Element('beast', version='1.0')

# Format dates and sequences
taxa, seq_list = read_fasta(filename)
dates = parse_dates(taxa)
date_decimal, date_precision = format_dates(dates)

# Write taxa and alignment blocks
test = write_taxa_block(root, taxa, date_decimal, date_precision)
test = write_alignment_block(test, taxa, seq_list)

# write patterns block
if args['partitions']:
    test = write_patterns_block(test, args['partitions'])


# population
test = write_treeprior_block(test, args['population_model'])

# tree model
test = write_treemodel_block(test, taxa, date_precision)

# tree & tmrca statistics
test = write_treelengthstatistic_block(test)
test = write_tmrcastatistic_block(test)

# write population prior:
if re.search(args['population_model'], 'constant'):
    test = write_coalescentlikelihood_block(test)

elif re.search(args['population_model'], 'skygrid'):
    test = write_skygridlikelihood_block(test)


# write block for clock models
if re.search(args['clock_model'], 'ucld'):
    test = write_relaxedclock_block(test, args)

#elif re.search(args['clock_model'], 'strict'):
    #test

# write HKY model
if args['partitions']:
    for partition in args['partitions']:
        test = write_hky_block(test, partition, args)
else:
    test = write_hky_block(test, '', args)

if args['partitions']:
    for partition in args['partitions']:
     test = write_site_block(test, partition, args)
else:
    test = write_site_block(test, '', args)

# Compound parameter only called if partitions are present
if args['partitions']:
    test = write_compound_block(test, args)

# Likelihood of tree given sequence data
test = write_treedatalikelihood_block(test, args)

# Operators
test = write_operator_block(test, args, date_precision, taxa)

# MCMC
test = write_mcmc(test, args, date_precision, taxa)

# Report
test = write_report(test)

# Save to file
xml_string = etree.tostring(test, pretty_print=True, encoding="utf-8").decode()
with open("output.xml", "w", encoding="utf-8") as f:
    f.write(xml_string)