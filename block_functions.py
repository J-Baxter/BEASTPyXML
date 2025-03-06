import re
from datetime import datetime
from lxml import etree
from Bio import AlignIO

# Import fasta
def read_fasta(filename):
    with open(filename) as handle:
        alignment = AlignIO.read(handle, "fasta")

    # Extract taxon names
    taxa = []
    for record in alignment:
        taxa.append(record.id)

    # Extract sequences
    sequence = []
    for record in alignment:
        sequence.append(record.seq)

    return taxa, sequence


# If dates are within fasta name, extract dates
def parse_dates(taxa):
    dates = []
    for taxon in taxa:
        match = re.search(r'(\d{4}-\d{2}-\d{2}$)', taxon)
        if match:
            dates.append(match.group(0))
        else:
            match = re.search(r'(\d{4}-\d{2}$)', taxon)
            if match:
                dates.append(match.group(0))
            else:
                match = re.search(r'(\d{4}$)', taxon)
                if match:
                    dates.append(match.group(0))
                else:
                    dates.append(None)
    return dates


# convert dates to decimals
def decimal_date(date):
    dt = datetime.strptime(date, "%Y-%m-%d")
    start_of_year = datetime(dt.year, 1, 1)
    start_of_next_year = datetime(dt.year + 1, 1, 1)

    # Compute the fraction of the year
    year_fraction = (dt - start_of_year).total_seconds() / (start_of_next_year - start_of_year).total_seconds()

    # Compute decimal year
    return dt.year + year_fraction


# format parsed dates (also calculates uncertainty)
def format_dates(dates):
    date_decimal = []
    precision = []

    for date in dates:
        if date.count('-') == 2:
            dt = decimal_date(date)
            date_decimal.append(round(dt, 5))
            precision.append(0)

        elif date.count('-') == 1:
            tmp = date + '-15'
            dt = decimal_date(tmp)
            date_decimal.append(round(dt, 4))
            precision.append(0.08333333333333333)

        elif date.count('-') == 0:
            tmp = date + '-07-02'
            dt = decimal_date(tmp)
            date_decimal.append(round(dt, 0))
            precision.append(1.0)

    return date_decimal, precision


# Write basic taxa block
def write_taxa_block(x, taxa, dates_formatted, precision):
    taxa_block = etree.SubElement(x, 'taxa', id='taxa')
    n = len(taxa)

    for i in range(n):
        date_string =  str(dates_formatted[i])
        tmp = etree.SubElement(taxa_block, "taxon", id=taxa[i])

        if precision[i] > 0:
            precision_string = str(precision[i])
            etree.SubElement(tmp,
                             "date",
                             value=date_string,
                             direction='forwards',
                             units='years',
                             uncertainty=precision_string)
        else:
            etree.SubElement(tmp,
                             "date",
                             value=date_string,
                             direction='forwards',
                             units='years')

    return x


# Write sequence alignments
def write_alignment_block(x, taxa, sequences):
    aln_block = etree.SubElement(x, 'alignment', id='alignment', dataType='nucleotide')

    n = len(taxa)
    for i in range(n):
        seq = etree.SubElement(aln_block, "sequence")
        etree.SubElement(seq,"taxon",idref=taxa[i])
        seq.text = str(sequences[i])

    return x


# Write initial patterns block
def write_patterns_block(x, partition):
    if not partition:
        tmp2 = etree.SubElement(x, 'patterns', id='patterns', attrib={"from": '1'}, strip="false")
        etree.SubElement(tmp2, 'alignment', idref='alignment')

    else:
        for p in partition:
            if isinstance(p, list):
                name = 'CP' + str(p[0]) + '+' + str(p[1])
                tmp = etree.SubElement(x, 'mergePatterns', id=name+'.patterns')

                n_patterns = len(p)
                for i in range(n_patterns):
                    tmp2 = etree.SubElement(tmp, 'patterns', attrib={"from":str(p[i])}, every="3", strip="false")
                    etree.SubElement(tmp2, 'alignment', idref='alignment')

            else:
                name = 'CP' + str(p)
                tmp = etree.SubElement(x, 'patterns', id=name+'.patterns', attrib={"from":str(p)}, every="3", strip="false")
                etree.SubElement(tmp, 'alignment', idref='alignment')

    return x


# Write tree prior block (where treeprior is a string - currently either 'constant' or 'skygrid')
def write_treeprior_block(x, treeprior):
    tmp = etree.SubElement(x, 'constantSize', units="years")
    tmp2 = etree.SubElement(tmp, "populationSize")
    tmp3 = etree.SubElement(tmp2, "parameter", id="constant.popSize")

    tmp4 = etree.SubElement(x, 'coalescentSimulator', id="startingTree")
    etree.SubElement(tmp4, "taxa", idref="taxa")
    tmp5 = etree.SubElement(tmp4, "constantSize")

    if re.match(treeprior, 'constant'):
        tmp.set('id', "constant")
        tmp3.set('id', "constant.popSize")
        tmp3.set('value', "1.0")
        tmp3.set('lower', "0.0")
        tmp5.set('idref', "constant")

    elif re.match(treeprior, 'skygrid'):
        tmp.set('id', "initialDemo")
        tmp3.set('id', "initialDemo.popSize")
        tmp3.set('value', "100")
        tmp5.set('idref', "initialDemo")

    return x


# write basic tree model block
def write_treemodel_block(x, taxa, precision):
    treemodel_block = etree.SubElement(x, 'treeModel', id='treeModel')
    etree.SubElement(treemodel_block, "coalescentTree", idref="startingTree")

    tmp = etree.SubElement(treemodel_block, "rootHeight")
    etree.SubElement(tmp, "parameter", id="treeModel.rootHeight")

    tmp = etree.SubElement(treemodel_block, "nodeHeights", internalNodes='true')
    etree.SubElement(tmp, "parameter", id="treeModel.internalNodeHeights")

    tmp = etree.SubElement(treemodel_block, "nodeHeights", internalNodes='true', rootNode="true")
    etree.SubElement(tmp, "parameter", id="treeModel.allInternalNodeHeights")

    n = len(taxa)
    for i in range(n):
        if precision[i] > 0:
            tmp = etree.SubElement(treemodel_block, "leafHeight", taxon=taxa[i])
            etree.SubElement(tmp, "parameter", id='age(' + taxa[i] + ')')

    return x


def write_treelengthstatistic_block(x):
    tmp = etree.SubElement(x, 'treeLengthStatistic', id='treeLength')
    etree.SubElement(tmp, "treeModel", idref="treeModel")

    return x


def write_tmrcastatistic_block(x):
    tmp = etree.SubElement(x, 'tmrcaStatistic', id="age(root)", absolute="true")
    etree.SubElement(tmp, "treeModel", idref="treeModel")

    return x


# Coalescent Likelihood Block
def write_coalescentlikelihood_block(x):
    tmp = etree.SubElement(x, 'coalescentLikelihood', id="coalescent")
    tmp2 = etree.SubElement(tmp, "model")
    etree.SubElement(tmp2, "constantSize", idref="constant")
    tmp2 = etree.SubElement(tmp, "populationTree")
    etree.SubElement(tmp2, "treeModel", idref="treeModel")

    return x


# Skygrid Likelihood Block
def write_skygridlikelihood_block(x, parameters):

    # Extract parameters
    sg_popsize = parameters.skygrid_grids
    sg_gridpoints = str(float(parameters.skygrid_grids)+1)
    sg_cutoff = parameters.skygrid_cutoff

    # Write XML
    tmp = etree.SubElement(x, 'gmrfSkyGridLikelihood', id="skygrid")
    tmp2 = etree.SubElement(tmp, "populationSizes")
    etree.SubElement(tmp2, "parameter", id="skygrid.logPopSize", dimension=sg_popsize, value='1.0')

    tmp2 = etree.SubElement(tmp, "precisionParameter")
    etree.SubElement(tmp2, "parameter", id="skygrid.precision", value='0.1', lower = '0.0')

    tmp2 = etree.SubElement(tmp, "numGridPoints")
    etree.SubElement(tmp2, "parameter", id="skygrid.numGridPoints", value=sg_gridpoints)

    tmp2 = etree.SubElement(tmp, "cutOff")
    etree.SubElement(tmp2, "parameter", id="skygrid.cutOff", value=sg_cutoff)

    tmp2 = etree.SubElement(tmp, "populationTree")
    etree.SubElement(tmp2, "treeModel", idref="treeModel")

    return x


def write_relaxedclock_block(x, parameters):
    ucld_mean = re.split(',', parameters.ucld_mean)[1]
    ucld_stdev = re.split(',', parameters.ucld_stdev)[1]

    # Write XML
    # discretised branch rates
    tmp = etree.SubElement(x, 'discretizedBranchRates', id="branchRates")
    etree.SubElement(tmp, "treeModel", idref="treeModel")
    # distribution
    tmp2 = etree.SubElement(tmp, "distribution")
    tmp3 = etree.SubElement(tmp2, "logNormalDistributionModel", meanInRealSpace='true')
    tmp4 = etree.SubElement(tmp3, "mean")
    etree.SubElement(tmp4, "parameter", id="ucld.mean", value=ucld_mean, lower="0.0")
    tmp4 = etree.SubElement(tmp3, "stdev")
    etree.SubElement(tmp4, "parameter", id="ucld.stdev", value=ucld_stdev, lower="0.0")
    # rate categories
    tmp2 = etree.SubElement(tmp, "rateCategories")
    etree.SubElement(tmp2, "parameter", id="branchRates.categories")

    # rate statistics
    block_names = ['meanRate', 'coefficientOfVariation']

    for block_name in block_names:
        tmp = etree.SubElement(x, 'rateStatistic', id=block_name, name=block_name, mode='mean', internal="true", external="true")
        etree.SubElement(tmp, "treeModel", idref="treeModel")
        etree.SubElement(tmp, "branchRates", idref="branchRates")

    #covariance statistic
    tmp = etree.SubElement(x, 'rateCovarianceStatistic', id='covariance', name='covariance')
    etree.SubElement(tmp, "treeModel", idref="treeModel")
    etree.SubElement(tmp, "branchRates", idref="branchRates")

    return x


def write_strictclock_block(x):
    tmp = etree.SubElement(x, 'strictClockBranchRates', id="branchRates")
    tmp2 = etree.SubElement(tmp, 'rate')
    etree.SubElement(tmp2, "parameter", id='clock.rate', value='1.0')

    tmp = etree.SubElement(x, 'rateStatistic', id='meanRate', name='meanRate', mode='mean', internal="true",
                           external="true")
    etree.SubElement(tmp, "treeModel", idref="treeModel")
    etree.SubElement(tmp, "branchRates", idref="branchRates")

    return x


# HKY -remember we're passing partition outside the function:


def write_hky_block(x, partition):
    if not partition:
        name = ''

    else:
        if isinstance(partition, list):
            name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
        else:
            name = 'CP' + str(partition) + '.'

    tmp = etree.SubElement(x, 'HKYModel', id=name+'hky')

    # Frequency block
    if not x.findall('.//frequencies'):
        tmp2 = etree.SubElement(tmp, 'frequencies')
        tmp3 = etree.SubElement(tmp2, 'frequencyModel', dataType='nucleotide')
        tmp4 = etree.SubElement(tmp3, 'frequencies')
        etree.SubElement(tmp4, 'parameter', id='frequencies', value='0.25 0.25 0.25 0.25')
    else:
        tmp2 = etree.SubElement(tmp, 'frequencies')
        tmp3 = etree.SubElement(tmp2, 'frequencyModel', dataType='nucleotide')
        tmp4 = etree.SubElement(tmp3, 'frequencies')
        etree.SubElement(tmp4, 'parameter', idref='frequencies')


    # Kappa block
    tmp2 = etree.SubElement(tmp, 'kappa')
    etree.SubElement(tmp2, 'parameter', id= name+'kappa', value='1.0', lower='0.0')

    return x


def write_gtr_block(x, partition):
    # Extract parameters
    #gtr_rates_value = parameters["gtr_rates_value"]
    #gtr_rates_dimension = parameters["gtr_rates_dimension"]

    #gtr_rate_ac = parameters.gtr_ac
    #gtr_rate_ag = parameters.gtr_ag
    #gtr_rate_at = parameters.gtr_at
    #gtr_rate_cg = parameters.gtr_cg
    #gtr_rate_gt = parameters.gtr_gt

    #gtr_rates_value = [gtr_rate_ac, gtr_rate_ag, gtr_rate_at, gtr_rate_cg, gtr_rate_gt]

    if not partition:
        name = ''

    else:
        if isinstance(partition, list):
            name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
        else:
            name = 'CP' + str(partition) + '.'

    tmp = etree.SubElement(x, 'gtrModel', id=name + 'gtr')

    # Frequency block
    if not x.findall('.//frequencies'):
        tmp2 = etree.SubElement(tmp, 'frequencies')
        tmp3 = etree.SubElement(tmp2, 'frequencyModel', dataType='nucleotide')
        tmp4 = etree.SubElement(tmp3, 'frequencies')
        etree.SubElement(tmp4, 'parameter', id='frequencies', value='0.25 0.25 0.25 0.25')
    else:
        tmp2 = etree.SubElement(tmp, 'frequencies')
        tmp3 = etree.SubElement(tmp2, 'frequencyModel', dataType='nucleotide')
        tmp4 = etree.SubElement(tmp3, 'frequencies')
        etree.SubElement(tmp4, 'parameter', idref='frequencies')

    # rates block

    #if all(i == gtr_rates_value[0] for i in gtr_rates_value):
        #tmp2 = etree.SubElement(tmp, 'rates')
        #etree.SubElement(tmp2, 'parameter', id=name + 'gtr.rates', value=gtr_rates_value[0], dimension='6', lower='0.0')

   # else:
    rates = ['AC', 'AG', 'AT', 'CG', 'GT']
    for i in range(len(rates)):
        tmp2 = etree.SubElement(tmp, 'rate'+rates[i])
        etree.SubElement(tmp2, 'parameter', id=name+'gtr.'+rates[i], value='1.0', lower='0.0')

    return x


# Site model
def write_site_block(x, partition, parameters):
    # Extract parameters

    substitution_model = parameters.substitution_model
    use_gamma = parameters.use_gamma
    gamma_categories = parameters.gamma_categories
    #gamma_alpha = parameters["gamma_alpha"]

    if not partition:
        name = ''

    else:
        if isinstance(partition, list):
            name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
        else:
            name = 'CP' + str(partition) + '.'

    # Link to substitution model
    tmp = etree.SubElement(x, 'siteModel', id=name+'siteModel')

    if re.search('hky', substitution_model):
        tmp2 = etree.SubElement(tmp, 'substitutionModel')
        etree.SubElement(tmp2, substitution_model.upper()+'Model', idref=name+substitution_model)

        tmp3 = etree.SubElement(tmp, 'relativeRate')
        etree.SubElement(tmp3, 'parameter', id=name + 'mu', value='1.0', lower='0.0')

    elif re.search('gtr', substitution_model):
        tmp2 = etree.SubElement(tmp, 'substitutionModel')
        etree.SubElement(tmp2, substitution_model + 'Model', idref=name + substitution_model)

        tmp3 = etree.SubElement(tmp, 'relativeRate')
        etree.SubElement(tmp3, 'parameter', id=name + 'mu', value='1.0', lower='0.0')

    if use_gamma:
        tmp4 = etree.SubElement(tmp, 'gammaShape', gammaCategories=gamma_categories)
        etree.SubElement(tmp4, 'parameter',  id=name+'alpha', value='0.5', lower='0.0')

    return x


# Define compound parameter
# Only called if partition is present
def write_compound_block(x, parameters):
    tmp = etree.SubElement(x, 'compoundParameter', id='allMus')
    partitions = parameters.partitions

    for partition in partitions:
        if isinstance(partition, list):
            name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
        else:
            name = 'CP' + str(partition) + '.'

        etree.SubElement(tmp, 'parameter', idref=name+'mu')

    return x


# Write block describing likelihood for tree given sequence data
# Must loop through partitions within function
def write_treedatalikelihood_block(x, parameters):
    partitions = parameters.partitions
    clock_model = parameters.clock_model

    tmp = etree.SubElement(x, 'treeDataLikelihood', id='treeLikelihood', useAmbiguities="false")

    if not partitions:
        tmp2 = etree.SubElement(tmp, 'partition')
        etree.SubElement(tmp2, 'patterns', idref='patterns')
        etree.SubElement(tmp2, 'siteModel', idref='siteModel')

    else:
        for partition in partitions:
            if isinstance(partition, list):
                name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
            else:
                name = 'CP' + str(partition) + '.'

            tmp2 = etree.SubElement(tmp, 'partition')
            etree.SubElement(tmp2, 'patterns', idref=name+'patterns')
            etree.SubElement(tmp2, 'siteModel', idref=name+'siteModel')

    etree.SubElement(tmp, 'treeModel', idref="treeModel")

    if re.match(clock_model, 'strict'):
        etree.SubElement(tmp, 'strictClockBranchRates', idref="branchRates")

    elif re.match(clock_model, 'ucld'):
        etree.SubElement(tmp, 'discretizedBranchRates', idref="branchRates")

    return x
