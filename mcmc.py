import re
from lxml import etree

def write_prior(x, parameter_name, prior_string):

    dist, params = re.split(',', prior_string,1)
    params = params.split(',')

    if dist == 'lognormal':
            if params[2]:
                write_lognormal_prior(x,
                                      parameter_name,
                                      mean_in_real_space = 'true',
                                      mean = params[0],
                                      stdev = params[1])
            else:
                write_lognormal_prior(x,
                                      parameter_name,
                                      mean_in_real_space = 'false',
                                      mu = params[0],
                                      sigma = params[1])

    #if dist == 'normal':
        #write_normal_prior(x,
                           #parameter_name,
                           #mean=params[0],
                           #stdev=params[1])

    if dist == 'uniform':
        write_uniform_prior(x,
                           parameter_name,
                           lower=params[0],
                           upper=params[1])

    if dist == 'gamma':
        write_gamma_prior(x,
                          parameter_name,
                          # offset=params[0],
                          scale=params[0],
                          shape=params[1])

    if dist == 'exponential':
        write_exponential_prior(x,
                                parameter_name,
                               # offset=params[0],
                               # scale=params[0],
                                mean=params[0])

    return x



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


# 1/X prior - this is an improper prior so not included
#def write_oneonx_prior(x, parameter_name):

    # Write Attribute
    #tmp = etree.SubElement(x, 'oneOnXPrior')
    #etree.SubElement(tmp, 'parameter', idref=parameter_name)

    #return x

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


def write_dirichlet_prior(x, parameter_name, **kwargs):
    sums_to = kwargs.get('sums_to', None)
    alpha = kwargs.get('alpha', None)
    counts = kwargs.get('counts', None)

    attributes = {"sumsTo": sums_to,
                  "alpha": alpha,
                  "counts": counts}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write Attribute
    tmp = etree.SubElement(x, 'dirichletPrior', clean_attributes)
    etree.SubElement(tmp, 'parameter', idref=parameter_name)

    return x

def write_asr_block(x, parameters):

    possible_traits = {key: vars(parameters).get(key) for key in ['continuous_phylogeo', 'asr_sequence', 'dta']}

    traits_to_reconstruct = {key: value for (key, value) in possible_traits.items() if value is not None}
    #traits_to_reconstruct = {key: value for key, value in vars(parameters).items() if value is not None and value is not False}

    for key,value in traits_to_reconstruct.items():
        if key == 'continuous_phylogeo':
            name = 'location'
            tag = 'location'
            idref = 'location.traitLikelihood'
            element = 'multivariateTraitLikelihood'

        #elif key == 'asr_sequence':
            #pass

        else:
            name = key+'.states'
            tag = key
            idref = key+'.traitLikelihood'
            element = 'ancestralTreeLikelihood'

        tmp = etree.SubElement(x, 'trait', name=name, tag=tag)
        etree.SubElement(tmp, element, idref=idref)

    return x


def write_prior_block(x, parameters):
    tmp = etree.SubElement(x, 'prior', id='prior')

    # HKY substitution model
    if parameters.substitution_model == "hky":
        if not parameters.partitions:
            write_prior(tmp, 'kappa', parameters.hky_kappa)
            #write_lognormal_prior(tmp, 'kappa', mu="1.0", sigma="1.25", offset="0.0")
        else:
            for partition in parameters.partitions:
                if isinstance(partition, list):
                    name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
                else:
                    name = 'CP' + str(partition) + '.'

                #write_lognormal_prior(tmp, name+'kappa', mu="1.0", sigma="1.25", offset="0.0")
                write_prior(tmp, name+'kappa', parameters.hky_kappa)
            #write_uniform_prior(tmp, 'allMus', lower='0.0', upper='100.0')
            write_prior(tmp, 'allMus', parameters.all_mus)
        #write_uniform_prior(tmp, 'frequencies', lower='0.0', upper='1.0')


    # GTR substitution model
    if parameters.substitution_model == "gtr":
        rates = ['AC', 'AG', 'AT', 'CG', 'GT']
        if not parameters.partitions:
            for i in range(len(rates)):
                param = f"gtr_{rates[i].lower()}"
                write_prior(tmp, 'gtr.'+rates[i], getattr(parameters, param))
                #write_dirichlet_prior(tmp, "gtr.rates", alpha="1.0", sums_to="6.0")
        else:
            for partition in parameters.partitions:
                if isinstance(partition, list):
                    name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
                else:
                    name = 'CP' + str(partition) + '.'

                for i in range(len(rates)):
                    param = f"gtr_{rates[i].lower()}"
                    write_prior(tmp, name+  'gtr.' + rates[i], getattr(parameters, param))

                    #write_dirichlet_prior(tmp, name + "gtr.rates", alpha="1.0", sums_to="6.0")
            write_prior(tmp, 'allMus', parameters.all_mus)
        #write_dirichlet_prior(tmp, 'frequencies', alpha="1.0", sums_to="1.0")

    # Base frequencies
    #write_uniform_prior(tmp, 'frequencies', lower='0.0', upper='1.0')
    write_dirichlet_prior(tmp, 'frequencies', alpha="1.0", sums_to="1.0")

    # Gamma heterogeneity across sites
    if parameters.use_gamma:
        if not parameters.partitions:
            write_prior(tmp, 'alpha', parameters.gamma_alpha)
            #write_exponential_prior(tmp, 'alpha', mean="0.5", offset="0.0")
        else:
            for partition in parameters.partitions:
                if isinstance(partition, list):
                    name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
                else:
                    name = 'CP' + str(partition) + '.'
                write_prior(tmp,  name+'alpha', parameters.gamma_alpha)
                #write_exponential_prior(tmp, name+'alpha', mean="0.5", offset="0.0")

    # Uncorrelated lognormal relaxed clock
    if parameters.clock_model == 'ucld':
        write_prior(tmp, 'ucld.mean', parameters.ucld_mean)
        write_prior(tmp, 'ucld.stdev', parameters.ucld_stdev)


        #write_lognormal_prior(tmp, 'ucld.mean', mean= parameters['ucld_mean_mean'], stdev= parameters['ucld_mean_stdev'], offset="0.0")
        #write_exponential_prior(tmp, 'ucld.stdev', mean="0.3333333333333333", offset="0.0")
        etree.SubElement(tmp, 'discretizedBranchRates', idref="branchRates")

    # Strict clock
    if parameters.clock_model == 'strict':
        write_prior(tmp, 'clock.rate', parameters.clock_rate)
        etree.SubElement(tmp, 'strictClockBranchRates', idref="branchRates")


    # Tree model operators
    if not parameters.empirical_tree_distribution:

        # Constant population
        if parameters.tree_model == "constant":
            etree.SubElement(tmp, 'coalescentLikelihood', idref="coalescent")
            #write_oneonx_prior(tmp, "constant.popSize")
            write_prior(tmp, "constant.popSize", parameters.constant_population)

        # Non parameteric skygrid
        if parameters.tree_model == "skygrid":
            write_gamma_prior(tmp, 'skygrid.precision', shape="0.001", scale="1000.0", offset="0.0")
            etree.SubElement(tmp, 'gmrfSkyGridLikelihood', idref="skygrid")

    # Traits go here
    if parameters.continuous_phylogeo:
        etree.SubElement(tmp, 'distributionLikelihood', idref="location.diffusion.prior")
        etree.SubElement(tmp, 'multivariateWishartPrior', idref="location.precisionPrior")

    return x

def write_likelihood(x, parameters):
    tmp = etree.SubElement(x, 'likelihood', id='likelihood')
    etree.SubElement(tmp, 'treeDataLikelihood', idref='treeLikelihood')

    if parameters.continuous_phylogeo:
        etree.SubElement(tmp, 'multivariateTraitLikelihood', idref="location.traitLikelihood")

    return x

def write_column(x, label, element):
    tmp = etree.SubElement(x, 'column', label=label, dp='4', width='12')
    etree.SubElement(tmp, element, idref=label.lower())
    return x

def write_screenlog(x, parameters):
    tmp = etree.SubElement(x, 'log', id='screenlog', logEvery=parameters.log_every)
    write_column(tmp,'Joint', 'joint')
    write_column(tmp, 'Prior', 'prior')
    write_column(tmp, 'Likelihood', 'likelihood')
    write_column(tmp, 'age(root)', 'tmrcaStatistic')

    if parameters.clock_model == 'ucld':
        write_column(tmp, 'ucld.mean', 'parameter')

    elif parameters.clock_model == 'strict':
        write_column(tmp, 'clock.rate', 'parameter')

    return x


def write_filelog(x, parameters, precision, taxa):
    tmp = etree.SubElement(x, 'log', id='fileLog', logEvery=parameters.log_every, fileName=parameters.file_stem+'.log', overwrite='false')
    etree.SubElement(tmp, 'joint', idref='joint')
    etree.SubElement(tmp, 'prior', idref='prior')
    etree.SubElement(tmp, 'likelihood', idref='likelihood')
    etree.SubElement(tmp, 'parameter', idref='treeModel.rootHeight')
    etree.SubElement(tmp, 'tmrcaStatistic', idref='age(root)')
    etree.SubElement(tmp, 'treeLengthStatistic', idref='treeLength')

    if not parameters.empirical_tree_distribution:

        # Constant population
        if parameters.tree_model == "constant":
            etree.SubElement(tmp, 'parameter', idref="constant.popSize")

        # Non parameteric skygrid
        if parameters.tree_model == "skygrid":
            etree.SubElement(tmp, 'parameter', idref="skygrid.precision")
            etree.SubElement(tmp, 'parameter', idref="skygrid.logPopSize")
            etree.SubElement(tmp, 'parameter', idref="skygrid.cutOff")

    # HKY substitution model
    if parameters.substitution_model == 'hky':
        if not parameters.partitions:
            etree.SubElement(tmp, 'parameter', idref="kappa")
        else:
            for partition in parameters.partitions:
                if isinstance(partition, list):
                    name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
                else:
                    name = 'CP' + str(partition) + '.'

                etree.SubElement(tmp, 'parameter', idref=name+"kappa")


    # GTR substitution model
    if parameters.substitution_model == 'gtr':
        rates = ['AC', 'AG', 'AT', 'CG', 'GT']

        if not parameters.partitions:
            for i in range(len(rates)):
                etree.SubElement(tmp, 'parameter', idref="gtr."+rates[i])
        else:
            for partition in parameters.partitions:
                if isinstance(partition, list):
                    name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
                else:
                    name = 'CP' + str(partition) + '.'

                for i in range(len(rates)):
                    etree.SubElement(tmp, 'parameter', idref=name + "gtr." + rates[i])


    # Base frequencies
    etree.SubElement(tmp, 'parameter', idref="frequencies")

    if parameters.partitions:
        etree.SubElement(tmp, 'compoundParameter', idref="allMus")

    # Gamma heterogeneity across sites
    if parameters.use_gamma:
        if not parameters.partitions:
            etree.SubElement(tmp, 'parameter', idref="alpha")
        else:
            for partition in parameters.partitions:
                if isinstance(partition, list):
                    name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
                else:
                    name = 'CP' + str(partition) + '.'

                etree.SubElement(tmp, 'parameter', idref=name+"alpha")

    # Uncorrelated lognormal relaxed clock
    if parameters.clock_model == 'ucld':
        etree.SubElement(tmp, 'parameter', idref="ucld.mean")
        etree.SubElement(tmp, 'parameter', idref="ucld.stdev")
        etree.SubElement(tmp, 'discretizedBranchRates', idref="branchRates")
        etree.SubElement(tmp, 'rateStatistic', idref="meanRate")
        etree.SubElement(tmp, 'rateStatistic', idref="coefficientOfVariation")
        etree.SubElement(tmp, 'rateCovarianceStatistic', idref="covariance")

    # Strict clock
    if parameters.clock_model == 'strict':
        etree.SubElement(tmp, 'rateStatistic', idref="meanRate")
        etree.SubElement(tmp, 'parameter', idref="clock.rate")
        etree.SubElement(tmp, 'strictClockBranchRates', idref="branchRates")

    # Precision sampling for tip-dates
    if precision:
        [etree.SubElement(tmp, 'parameter', idref='age(' + taxa[i] + ')') for i, z in enumerate(precision) if z > 0]

    # Continuous Phylogeography
    if parameters.continuous_phylogeo:
        etree.SubElement(tmp, 'matrixParameter', idref="location.precision")
        etree.SubElement(tmp, 'correlation', idref="location.correlation")
        etree.SubElement(tmp, 'matrixInverse', idref="location.varCovar")
        etree.SubElement(tmp, 'continuousDiffusionStatistic', idref="location.diffusionRate")
        etree.SubElement(tmp, 'multivariateTraitLikelihood', idref="location.traitLikelihood")

    # Tree model operators
    if not parameters.empirical_tree_distribution:
        etree.SubElement(tmp, 'treeDataLikelihood', idref='treeLikelihood')

    # Constant population
        if parameters.tree_model == "constant":
            etree.SubElement(tmp, 'coalescentLikelihood', idref='coalescent')

    # Non parameteric skygrid
        if parameters.tree_model == "skygrid":
            etree.SubElement(tmp, 'gmrfSkyGridLikelihood', idref='skygrid')

    return x


def write_treelog(x, parameters):
    tmp = etree.SubElement(x, 'logTree', id="treeFileLog", logEvery=parameters.log_every ,nexusFormat="true", fileName=parameters.file_stem + '.trees', sortTranslationTable="true")
    etree.SubElement(tmp, 'treeModel', idref='treeModel')

    if not parameters.empirical_tree_distribution:
        tmp2 = etree.SubElement(tmp, 'trait', name='rate', tag='rate')
        if parameters.clock_model == 'ucld':
            etree.SubElement(tmp2, 'discretizedBranchRates', idref='branchRates')
        if parameters.clock_model == 'strict':
            etree.SubElement(tmp2, 'strictClockBranchRates', idref="branchRates")

    etree.SubElement(tmp, 'joint', idref='joint')

    #traits go here
    if parameters.continuous_phylogeo:
        write_asr_block(x, parameters)

        etree.SubElement(tmp, 'multivariateDiffusionModel', idref="location.diffusionModel")
        etree.SubElement(tmp, 'multivariateTraitLikelihood', idref="location.traitLikelihood")
        tmp2 = etree.SubElement(tmp, 'trait', tag="location.rate")
        etree.SubElement(tmp2, 'arbitraryBranchRates', idref="location.diffusion.branchRates")

    return x



def write_mcmc(x, parameters, precision, taxa):
    tmp = etree.SubElement(x, 'mcmc', id='mcmc',  chainLength=parameters.chain_length, autoOptimize="true", operatorAnalysis= parameters.file_stem+'.ops')

    # joint block including prior and likelihood
    tmp2 = etree.SubElement(tmp, 'joint', id='joint')
    write_prior_block(tmp2, parameters = parameters)
    write_likelihood(tmp2, parameters = parameters)

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


