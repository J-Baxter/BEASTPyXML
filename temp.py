from csv import DictReader

# from WriteEddieShell.BEASTPyXML.empirical_tree_model import write_empirical_tree_model
from block_functions import *
from operators import *
from mcmc import *
from cmd_parser import *
from continuous_phylogeo import *
from empirical_tree_model import *


def main():
    args = parse_args()
    taxa, seq_list = read_fasta(args.fasta)

    if args.continuous_phylogeo:
        with open(args.continuous_trait_file, mode="r", encoding="utf-8") as file:
            reader = DictReader(file)  # Automatically maps headers to values
            latlon = {row["taxon"]: [str(row["lat"]), str(row["long"])] for row in reader}

    # Make base root of XML
    root = etree.Element('beast', version='1.0')

    # Format dates and sequences
    dates = parse_dates(taxa)
    date_decimal, date_precision = format_dates(dates)

    # Write taxa and alignment blocks
    test = write_taxa_block(root, taxa, date_decimal, date_precision)

    if not args.empirical_tree_model:
        test = write_alignment_block(test, taxa, seq_list)

    if args.continuous_phylogeo:
        test = write_taxon_traits(test, latlon)

    # write patterns block
    if not args.empirical_tree_model:
        test = write_patterns_block(test, args.partitions)

        # population
        test = write_treeprior_block(test, args.tree_model)

        # tree model
        test = write_treemodel_block(test, taxa, date_precision)

    else:
        # Empirical tree model block
        test = write_empirical_tree_model(test, args)

    # tree & tmrca statistics
    if not args.empirical_tree_model:
        test = write_treelengthstatistic_block(test)
        test = write_tmrcastatistic_block(test)

        # write population prior:
        if args.tree_model == 'constant':
            test = write_coalescentlikelihood_block(test)

        elif args.tree_model == 'skygrid':
            test = write_skygridlikelihood_block(test, args)

        # write block for clock models
        if args.clock_model == 'ucld':
            test = write_relaxedclock_block(test, args)

        elif args.clock_model == 'strict':
            test = write_strictclock_block(test)

        # write HKY model
        if args.substitution_model == 'hky':
            if args.partitions:
                for partition in args.partitions:
                    test = write_hky_block(test, partition)
            else:
                test = write_hky_block(test, '')

        # write GTR model
        elif args.substitution_model == 'gtr':
            if args.partitions:
                for partition in args.partitions:
                    test = write_gtr_block(test, partition)
            else:
                test = write_gtr_block(test, '')

        if args.partitions:
            for partition in args.partitions:
                test = write_site_block(test, partition, args)
        else:
            test = write_site_block(test, '', args)

        # Compound parameter only called if partitions are present
        if args.partitions:
            test = write_compound_block(test, args)

    if args.continuous_phylogeo:
        test = write_multivariatemodel_block(test)

    if not args.empirical_tree_model:
        # Likelihood of tree given sequence data
        test = write_treedatalikelihood_block(test, args)

    if args.continuous_phylogeo:
        test = write_gammarrw_block(test)
        test = write_gammarrwlikelihood_block(test, args)
        test = write_multivariatestats_block(test)

    # Operators
    test = write_operator_block(test, args, date_precision, taxa)

    # MCMC
    test = write_mcmc(test, args, date_precision, taxa)

    # Report
    test = write_report(test)

    # Save to file
    xml_string = etree.tostring(test, pretty_print=True, encoding="utf-8", xml_declaration=True, method="xml",
                                standalone="yes").decode()
    with open(args.file_stem + ".xml", "w", encoding="utf-8") as f:
        f.write(xml_string)
        f.close()


if __name__ == '__main__':
    main()
