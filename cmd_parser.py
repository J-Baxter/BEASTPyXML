import re

from WriteEddieShell.BEASTPyXML.block_functions import write_taxa_block, write_alignment_block, write_patterns_block, \
    write_treeprior_block, write_treemodel_block, write_treelengthstatistic_block, write_tmrcastatistic_block, \
    write_coalescentlikelihood_block, write_hky_block, write_site_block, write_compound_block, \
    write_skygridlikelihood_block, write_treedatalikelihood_block

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
         'partitions': [[1,2],3]}


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

test = write_treedatalikelihood_block(test, args)

xml_string = etree.tostring(test, pretty_print=True, encoding="utf-8").decode()

b_xml = etree.tostring(test, pretty_print=True, encoding="utf-8")
with open("test.xml", "wb") as f:
    f.write(b_xml)

xml_string = etree.tostring(root, pretty_print=True, encoding="utf-8").decode()

# Save to a file
with open("output.xml", "w", encoding="utf-8") as f:
    f.write(xml_string)