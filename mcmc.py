import re

from louis import logLevels
from lxml import etree

# Lognormal Prior
def write_lognormal_prior(x, parameter_name, **kwargs):
    mu = kwargs.get('mu', None)
    sigma = kwargs.get('sigma', None)
    offset = kwargs.get('offset', None)
    mean = kwargs.get('mean', None)
    stdev = kwargs.get('stdev', None)
    mean_in_real_space = kwargs.get('mean_in_real_space', None)

    attributes = {"mu": mu,
                  "sigma": sigma,
                  "offset": offset,
                  "meanInRealSpace": mean_in_real_space,
                  "mean": mean,
                  "stdev": stdev}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write Attribute
    tmp = etree.SubElement(x, 'logNormalPrior', clean_attributes)
    etree.SubElement(tmp, 'parameter', idref=parameter_name)

    return x


# Uniform prior
def write_uniform_prior(x, parameter_name, **kwargs):
    lower = kwargs.get('lower', None)
    upper = kwargs.get('upper', None)

    attributes = {"lower": lower,
                  "upper": upper}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write Attribute
    tmp = etree.SubElement(x, 'uniformPrior', clean_attributes)
    etree.SubElement(tmp, 'parameter', idref=parameter_name)

    return x


# Exponential prior
def write_exponential_prior(x, parameter_name, **kwargs):
    offset = kwargs.get('offset', None)
    scale = kwargs.get('scale', None)
    mean = kwargs.get('mean', None)


    attributes = {"offset": offset,
                  "scale": scale,
                  "mean":mean}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write Attribute
    tmp = etree.SubElement(x, 'exponentialPrior', clean_attributes)
    etree.SubElement(tmp, 'parameter', idref=parameter_name)

    return x


# Exponential prior
def write_oneonx_prior(x, parameter_name):

    # Write Attribute
    tmp = etree.SubElement(x, 'oneOnXPrior')
    etree.SubElement(tmp, 'parameter', idref=parameter_name)

    return x

def write_gamma_prior(x, parameter_name, **kwargs):
    offset = kwargs.get('offset', None)
    scale = kwargs.get('scale', None)
    shape = kwargs.get('shape', None)

    attributes = {"offset": offset,
                  "scale": scale,
                  "shape": shape}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write Attribute
    tmp = etree.SubElement(x, 'gammaPrior', clean_attributes)
    etree.SubElement(tmp, 'parameter', idref=parameter_name)

    return x


def write_prior(x, parameters):
    tmp = etree.SubElement(x, 'prior', id='prior')

    # HKY substitution model
    if re.search(parameters['substitution_model'], 'hky'):
        if not parameters['partitions']:
            write_lognormal_prior(tmp, 'kappa', mu="1.0", sigma="1.25", offset="0.0")
        else:
            for partition in parameters['partitions']:
                if isinstance(partition, list):
                    name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
                else:
                    name = 'CP' + str(partition) + '.'

                write_lognormal_prior(tmp, name+'kappa', mu="1.0", sigma="1.25", offset="0.0")
        write_uniform_prior(tmp, 'allMus', lower='0.0', upper='100.0')

    # GTR substitution model
    # if re.search(parameters['substitution_model']. 'gtr'): etc etc

    # Base frequencies
    write_uniform_prior(tmp, 'frequencies', lower='0.0', upper='1.0')

    # Gamma heterogeneity across sites
    if parameters['gamma_alpha']:
        if not parameters['partitions']:
            write_exponential_prior(tmp, 'alpha', mean="0.5", offset="0.0")
        else:
            for partition in parameters['partitions']:
                if isinstance(partition, list):
                    name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
                else:
                    name = 'CP' + str(partition) + '.'

                write_exponential_prior(tmp, name+'alpha', mean="0.5", offset="0.0")

    # Uncorrelated lognormal relaxed clock
    if re.search(parameters['clock_model'], 'ucld'):
        write_lognormal_prior(tmp, 'ucld.mean', mean="0.005", stdev="0.005", offset="0.0")
        write_exponential_prior(tmp, 'ucld.stdev', mean="0.3333333333333333", offset="0.0")
        etree.SubElement(tmp, 'discretizedBranchRates', idref="branchRates")

    # Strict clock
    if re.search(parameters['clock_model'], 'strict'):
        etree.SubElement(tmp, 'strictClockBranchRates', idref="branchRates")


    # Tree model operators
    if not parameters['empirical_tree_distribution']:
        etree.SubElement(tmp, 'coalescentLikelihood', idref="coalescent")

        # Constant population
        if re.search(parameters['population_model'], 'constant'):
            write_oneonx_prior(tmp, "constant.popSize")

        # Non parameteric skygrid
        if re.search(parameters['population_model'], 'skygrid'):
            write_gamma_prior(tmp, 'skygrid.precision', shape="0.001", scale="1000.0", offset="0.0")
            etree.SubElement(tmp, 'gmrfSkyGridLikelihood', idref="skygrid")


    # Traits go here
    return x

def write_likelihood(x):
    tmp = etree.SubElement(x, 'likelihood', id='likelihood')
    etree.SubElement(tmp, 'treeDataLikelihood', idref='treeLikelihood')

    return x

def write_column(x, label, element):
    tmp = etree.SubElement(x, 'column', label=label, dp='4', width='12')
    etree.SubElement(tmp, element, idref=label.lower())
    return x

def write_screenlog(x, parameters):
    tmp = etree.SubElement(x, 'log', id='screenlog', logEvery=parameters['log_every'])
    write_column(tmp,'Joint', 'joint')
    write_column(tmp, 'Prior', 'prior')
    write_column(tmp, 'Likelihood', 'likelihood')
    write_column(tmp, 'age(root)', 'tmrcaStatistic')

    if re.search(parameters['clock_model'], 'ucld'):
        write_column(tmp, 'ucld.mean', 'parameter')

    elif re.search(parameters['clock_model'], 'strict'):
        write_column(tmp, 'clock.rate', 'parameter')

    return x


def write_filelog(x, parameters, precision, taxa):
    tmp = etree.SubElement(x, 'log', id='fileLog', logEvery=parameters['log_every'], fileName=parameters['file_stem']+'.log', overwrite='false')
    etree.SubElement(tmp, 'joint', idref='joint')
    etree.SubElement(tmp, 'prior', idref='prior')
    etree.SubElement(tmp, 'likelihood', idref='treeModel.rootHeight')
    etree.SubElement(tmp, 'tmrcaStatistic', idref='age(root)')
    etree.SubElement(tmp, 'treeLengthStatistic', idref='treeLength')

    if not parameters['empirical_tree_distribution']:

        # Constant population
        if re.search(parameters['population_model'], 'constant'):
            etree.SubElement(tmp, 'parameter', idref="constant.popSize")

        # Non parameteric skygrid
        if re.search(parameters['population_model'], 'skygrid'):
            etree.SubElement(tmp, 'parameter', idref="skygrid.precision")
            etree.SubElement(tmp, 'parameter', idref="skygrid.logPopsize")
            etree.SubElement(tmp, 'parameter', idref="skygrid.cutOff")

    # HKY substitution model
    if re.search(parameters['substitution_model'], 'hky'):
        if not parameters['partitions']:
            etree.SubElement(tmp, 'parameter', idref="kappa")
        else:
            for partition in parameters['partitions']:
                if isinstance(partition, list):
                    name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
                else:
                    name = 'CP' + str(partition) + '.'

                etree.SubElement(tmp, 'parameter', idref=name+"kappa")
        etree.SubElement(tmp, 'compoundParameter', idref="allMus")

    # GTR substitution model
    # if re.search(parameters['substitution_model']. 'gtr'): etc etc

    # Base frequencies
    etree.SubElement(tmp, 'parameter', idref="frequencies")

    # Gamma heterogeneity across sites
    if parameters['gamma_alpha']:
        if not parameters['partitions']:
            etree.SubElement(tmp, 'parameter', idref="alpha")
        else:
            for partition in parameters['partitions']:
                if isinstance(partition, list):
                    name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
                else:
                    name = 'CP' + str(partition) + '.'

                etree.SubElement(tmp, 'parameter', idref=name+"alpha")

    # Uncorrelated lognormal relaxed clock
    if re.search(parameters['clock_model'], 'ucld'):
        etree.SubElement(tmp, 'parameter', idref="ucld.mean")
        etree.SubElement(tmp, 'parameter', idref="ucld.stdev")
        etree.SubElement(tmp, 'discretizedBranchRates', idref="branchRates")
        etree.SubElement(tmp, 'rateStatistic', idref="meanRate")
        etree.SubElement(tmp, 'rateStatistic', idref="coefficientOfVariation")
        etree.SubElement(tmp, 'rateCovarianceStatistic', idref="covariance")

    # Strict clock
    if re.search(parameters['clock_model'], 'strict'):
        etree.SubElement(tmp, 'strictClockBranchRates', idref="branchRates")

    # Precision sampling for tip-dates
    if precision:
        [etree.SubElement(tmp, 'parameter', idref='age(' + taxa[i] + ')') for i, z in enumerate(precision) if z > 0]

    # Tree model operators
    if not parameters['empirical_tree_distribution']:
        etree.SubElement(tmp, 'treeDataLikelihood', idref='treeLikelihood')

    # Constant population
        if re.search(parameters['population_model'], 'constant'):
            etree.SubElement(tmp, 'coalescentLikelihood', idref='coalescent')
    # Non parameteric skygrid
        if re.search(parameters['population_model'], 'skygrid'):
            etree.SubElement(tmp, 'gmrfSkyGridLikelihood', idref='skygrid')

    return x


def write_treelog(x, parameters):
    tmp = etree.SubElement(x, 'logTree', id="treeFileLog", logEvery=parameters['log_every'] ,nexusFormat="true", fileName=parameters['file_stem'] + '.trees', sortTranslationTable="true")
    etree.SubElement(tmp, 'treeModel', idref='treeModel')
    tmp2=etree.SubElement(tmp, 'trait', name='rate', tag='rate')

    if re.search(parameters['clock_model'], 'ucld'):
        etree.SubElement(tmp2, 'discretizedBranchRates', idref='branchRates')
    if re.search(parameters['clock_model'], 'strict'):
        etree.SubElement(tmp2, 'strictClockBranchRates', idref="branchRates")

    etree.SubElement(tmp, 'joint', idref='joint')

    #traits go here

    return x



def write_mcmc(x, parameters, precision, taxa):
    tmp = etree.SubElement(x, 'mcmc', id='mcmc',  chainLength=parameters['chain_length'], autoOptimize="true", operatorAnalysis= parameters['file_stem']+'.ops')

    # joint block including prior and likelihood
    tmp2 = etree.SubElement(tmp, 'joint', id='joint')
    write_prior(tmp2, parameters = parameters)
    write_likelihood(tmp2)

    # operator
    etree.SubElement(tmp, 'operators', idref='operators')

    # screen log
    write_screenlog(tmp, parameters = parameters)

    # file log
    write_filelog(tmp, parameters = parameters, precision = precision, taxa = taxa)

    #tree log
    write_treelog(tmp, parameters = parameters)

    return x

def write_report(x):
     tmp = etree.SubElement(x, 'report')
     tmp2 = etree.SubElement(tmp, 'property', name='timer')
     etree.SubElement(tmp2, 'mcmc', idref='mcmc')

     return x


