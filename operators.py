import re
from lxml import etree


def write_scaleoperator_block(x, parameter_name, **kwargs):
    scale_factor = kwargs.get('scale_factor', None)
    scale_all = kwargs.get('scale_all', None)
    scale_all_independently = kwargs.get('scale_all_independently', None)
    weight = kwargs.get('weight', None)
    auto_optimize = kwargs.get('auto_optimize', None)
    df = kwargs.get('df', None)

    attributes = {"scaleFactor": scale_factor,
                  "scaleAll": scale_all,
                  "scaleAllIndependently": scale_all_independently,
                  "weight": weight,
                  "autoOptimize": auto_optimize,
                  "df": df}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write Attribute
    tmp = etree.SubElement(x, 'scaleOperator', clean_attributes)
    etree.SubElement(tmp, 'parameter', idref=parameter_name)

    return x


def write_uniformoperator_block(x, parameter_name, **kwargs):
    weight = kwargs.get('weight', None)
    lower = kwargs.get('lower', None)
    upper = kwargs.get('upper', None)

    attributes = {"weight": weight,
                  "lower": lower,
                  "upper": upper}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write attribute
    tmp = etree.SubElement(x, 'uniformOperator', clean_attributes)
    etree.SubElement(tmp, 'parameter', idref=parameter_name)

    return x


def write_uniformintegeroperator_block(x, parameter_name, **kwargs):
    weight = kwargs.get('weight', None)
    lower = kwargs.get('lower', None)
    upper = kwargs.get('upper', None)
    count = kwargs.get('count', None)

    attributes = {"weight": weight,
                  "lower": lower,
                  "upper": upper,
                  "count": count}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write attribute
    tmp = etree.SubElement(x, 'uniformIntegerOperator', clean_attributes)
    etree.SubElement(tmp, 'parameter', idref=parameter_name)

    return x


def write_deltaexchange_block(x, parameter_name, **kwargs):
    delta = kwargs.get('delta', None)
    parameter_weights = kwargs.get('parameter_weights', None)
    weight = kwargs.get('weight', None)
    auto_optimize = kwargs.get('auto_optimize', None)
    integer = kwargs.get('integer', None)

    attributes = {"delta": delta,
                  "parameterWeights": parameter_weights,
                  "weight": weight,
                  "autoOptimize": auto_optimize,
                  "integer": integer}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write Attribute
    tmp = etree.SubElement(x, 'deltaExchange', clean_attributes)
    etree.SubElement(tmp, 'parameter', idref=parameter_name)

    return x


def write_narrowexchange_block(x, **kwargs):
    weight = kwargs.get('weight', None)

    attributes = {"weight": weight}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write Attribute
    tmp = etree.SubElement(x, 'narrowExchange', clean_attributes)
    etree.SubElement(tmp, 'treeModel', idref='treeModel')

    return x


def write_wideexchange_block(x, **kwargs):
    weight = kwargs.get('weight', None)

    attributes = {"weight": weight}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write Attribute
    tmp = etree.SubElement(x, 'wideExchange', clean_attributes)
    etree.SubElement(tmp, 'treeModel', idref='treeModel')

    return x


def write_updownoperator_block(x, up_parameter, down_parameter, **kwargs):
    scale_factor = kwargs.get('scale_factor', None)
    weight = kwargs.get('weight', None)
    auto_optimize = kwargs.get('auto_optimize', None)

    attributes = {"scaleFactor": scale_factor,
                  "weight": weight,
                  "autoOptimize": auto_optimize}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write Attribute
    tmp = etree.SubElement(x, 'upDownOperator', clean_attributes)

    up = etree.SubElement(tmp, 'up')
    etree.SubElement(up, 'parameter', idref=up_parameter)

    down = etree.SubElement(tmp, 'down')
    etree.SubElement(down, 'parameter', idref=down_parameter)

    return x


def write_subtreeslide_block(x, **kwargs):
    size = kwargs.get('size', None)
    weight = kwargs.get('weight', None)
    auto_optimize = kwargs.get('auto_optimize', None)
    target_acceptance = kwargs.get('target_acceptance', None)
    swap_in_random_rate = kwargs.get('swap_in_random_rate', None)
    swap_in_random_trait = kwargs.get('swap_in_random_trait', None)
    gaussian = kwargs.get('gaussian', None)

    attributes = {"size": size,
                  "weight": weight,
                  "autoOptimize": auto_optimize,
                  "targetAcceptance": target_acceptance,
                  "swapInRandomRate": swap_in_random_rate,
                  "swapInRandomTrait": swap_in_random_trait,
                  "gaussian": gaussian}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write Attribute
    tmp = etree.SubElement(x, 'subtreeSlide', clean_attributes)
    etree.SubElement(tmp, 'treeModel', idref='treeModel')

    return x


def write_wilsonbalding_block(x, **kwargs):
    weight = kwargs.get('weight', None)
    attributes = {"weight": weight}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write Attribute
    tmp = etree.SubElement(x, 'wilsonBalding', clean_attributes)
    etree.SubElement(tmp, 'treeModel', idref='treeModel')

    return x


def write_swapoperator_block(x, parameter_name, **kwargs):
    size = kwargs.get('size', None)
    weight = kwargs.get('weight', None)
    auto_optimize = kwargs.get('auto_optimize', None)

    attributes = {"size": size,
                  "weight": weight,
                  "autoOptimize": auto_optimize}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write Attribute
    tmp = etree.SubElement(x, 'swapOperator', clean_attributes)
    etree.SubElement(tmp, 'parameter', idref=parameter_name)

    return x


def write_gmrfupdateroperator_block(x, **kwargs):
    size = kwargs.get('size', None)
    weight = kwargs.get('weight', None)
    scale_factor = kwargs.get('scale_factor', None)

    attributes = {"size": size,
                  "weight": weight,
                  "scaleFactor": scale_factor}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write Attribute
    tmp = etree.SubElement(x, 'gmrfGridBlockUpdateOperator', clean_attributes)
    etree.SubElement(tmp, 'gmrfSkyrideLikelihood', idref='skygrid')

    return x


def write_precisiongibbs_block(x, **kwargs):
    weight = kwargs.get('weight', None)
    attributes = {"weight": weight}

    # Filter out missing arguments
    clean_attributes = {key: value for (key, value) in attributes.items() if value is not None}

    # Write Attribute
    tmp = etree.SubElement(x, 'precisionGibbsOperator', clean_attributes)
    etree.SubElement(tmp, 'multivariateTraitLikelihood', idref='location.traitLikelihood')
    etree.SubElement(tmp, 'multivariateWishartPrior', idref="location.precisionPrior")

    return x


def write_operator_block(x, parameters, precision, taxa):
    tmp = etree.SubElement(x, 'operators', id='operators', optimizationSchedule='default')

    if not parameters.empirical_tree_model:

        # HKY substitution model
        if parameters.substitution_model == 'hky':
            if not parameters.partitions:
                write_scaleoperator_block(tmp, 'kappa', scale_factor='0.75', weight='1')
            else:
                for partition in parameters.partitions:
                    if isinstance(partition, list):
                        name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
                    else:
                        name = 'CP' + str(partition) + '.'

                    write_scaleoperator_block(tmp, name + 'kappa', scale_factor='0.75', weight='1')
                write_deltaexchange_block(tmp, 'allMus', delta='0.01', parameter_weights="6870 3435", weight="3")
            write_deltaexchange_block(tmp, parameter_name='frequencies', delta='0.01', weight='1')

        # GTR substitution model
        if parameters.substitution_model == 'gtr':
            rates = ['AC', 'AG', 'AT', 'CG', 'GT']
            if not parameters.partitions:
                for i in range(len(rates)):
                    write_deltaexchange_block(tmp, "gtr." + rates[i], delta="0.01", weight="1")
            else:
                for partition in parameters.partitions:
                    if isinstance(partition, list):
                        name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
                    else:
                        name = 'CP' + str(partition) + '.'

                    for i in range(len(rates)):
                        write_deltaexchange_block(tmp, name + "gtr." + rates[i], delta="0.01", weight="1")
                write_deltaexchange_block(tmp, 'allMus', delta='0.01', parameter_weights="568 568 568", weight="3")
            write_deltaexchange_block(tmp, parameter_name='frequencies', delta='0.01', weight='1')

        # Gamma heterogeneity across sites
        if parameters.use_gamma:
            if not parameters.partitions:
                write_scaleoperator_block(tmp, 'alpha', scale_factor='0.75', weight='1')
            else:
                for partition in parameters.partitions:
                    if isinstance(partition, list):
                        name = 'CP' + str(partition[0]) + '+' + str(partition[1]) + '.'
                    else:
                        name = 'CP' + str(partition) + '.'

                    write_scaleoperator_block(tmp, name + 'alpha', scale_factor='0.75', weight='1')

        # Uncorrelated lognormal relaxed clock
        if parameters.clock_model == 'ucld':
            write_scaleoperator_block(tmp, "ucld.mean", scale_factor='0.75', weight='3')
            write_scaleoperator_block(tmp, "ucld.stdev", scale_factor='0.75', weight='3')
            write_updownoperator_block(tmp, "treeModel.allInternalNodeHeights", "ucld.mean", scale_factor='0.75',
                                       weight='3')
            write_swapoperator_block(tmp, "branchRates.categories", size='1', weight='10', auto_optimize='false')
            write_uniformintegeroperator_block(tmp, "branchRates.categories", weight='10')

        # Strict clock
        if parameters.clock_model == 'strict':
            pass

        # Tree model operators
        if not parameters.empirical_tree_distribution:
            write_subtreeslide_block(tmp, size='1.0', weight='30', gaussian='true')
            write_narrowexchange_block(tmp, weight='30')
            write_wideexchange_block(tmp, weight='3')
            write_wilsonbalding_block(tmp, weight='3')
            write_scaleoperator_block(tmp, "treeModel.rootHeight", scale_factor='0.75', weight='3')
            write_uniformoperator_block(tmp, "treeModel.internalNodeHeights", weight='30')

        # Constant population
        if parameters.tree_model == "constant":
            write_scaleoperator_block(tmp, 'constant.popSize', scale_factor='0.75', weight='3')

        # Non parameteric skygrid
        if parameters.tree_model == "skygrid":
            write_gmrfupdateroperator_block(tmp, scale_factor="1.0", weight="2")
            write_scaleoperator_block(tmp, 'skygrid.precision', scale_factor='0.75', weight='1')

        # Precision sampling for tip-dates
        if precision:
            [write_uniformoperator_block(tmp, 'age(' + taxa[i] + ')', weight='1') for i, z in enumerate(precision) if
             z > 0]

    # Traits go here
    if parameters.continuous_phylogeo:
        write_scaleoperator_block(tmp, 'location.diffusion.rates', scale_factor='0.75', weight='30')
        write_precisiongibbs_block(tmp, weight='2')

    if parameters.empirical_tree_model:
        tmp = etree.SubElement(x, 'empiricalTreeDistributionOperator', weight="3")
        etree.SubElement(tmp, 'empiricalTreeDistributionModel', idref='treeModel')

    return x
