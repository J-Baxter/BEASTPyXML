import re
from lxml import etree

def write_taxon_traits(x, parameters, traits):
    trait_names = parameters['trait_names']
    n_traits = len(trait_names)

    for n in range(n_traits):
        for key, value in traits.items():
            taxon = x.find(f".//taxon[@id='{key}']")

            if taxon is not None:
                attribute = etree.SubElement(taxon, 'attr', name=trait_names[n])

                if isinstance(value, list) and len(value) > n:
                    attribute.text = value[n]
                elif isinstance(value, str):
                    attribute.text = value
                else:
                    attribute.text = "NA"

    coords = [i for i, name in enumerate(trait_names) if re.search('lat|lon', name)]

    if len(coords) == 2:
        lat = next((i for i, name in enumerate(trait_names) if re.search('lat', name)), None)
        lon = next((i for i, name in enumerate(trait_names) if re.search('lon', name)), None)

        for key, value in traits.items():
            taxon = x.find(f".//taxon[@id='{key}']")

            if all(v is not None for v in [taxon, value[lat], value[lon]]):
                joint_attribute = etree.SubElement(taxon, 'attr', name='location')
                joint_attribute.text = str(value[lat] + ' ' + value[lon])

            else:
                joint_attribute = etree.SubElement(taxon, 'attr', name='location')
                joint_attribute.text = 'NA NA'

    return x


def write_multivariatemodel_block(x):
    tmp = etree.SubElement(x, 'multivariateDiffusionModel', id='location.diffusionModel')
    tmp2 = etree.SubElement(tmp, 'precisionMatrix')
    tmp3 = etree.SubElement(tmp2, 'matrixParameter', id="location.precision")
    etree.SubElement(tmp3, "parameter", id="location.precision.col1", value="0.05 0.002")
    etree.SubElement(tmp3, "parameter", id="location.precision.col2", value="0.002 0.05")

    tmp = etree.SubElement(x, 'multivariateWishartPrior', id='location.precisionPrior', df='2')
    tmp2 = etree.SubElement(tmp, 'scaleMatrix')
    tmp3 = etree.SubElement(tmp2, 'matrixParameter')
    etree.SubElement(tmp3, "parameter", value="1.0 0.0")
    etree.SubElement(tmp3, "parameter", value="0.0 1.0")

    tmp2 = etree.SubElement(tmp, 'data')
    etree.SubElement(tmp2, 'parameter', idref='location.precision')

    return x


def write_cauchyrrw_block(x):
    tmp = etree.SubElement(x, 'arbitraryBranchRates', id='location.diffusion.branchRates')
    etree.SubElement(tmp, "treeModel", idref="treeModel")
    tmp2 = etree.SubElement(tmp, 'rates')
    etree.SubElement(tmp2, 'parameter', id= "location.diffusion.rates" , lower='0.0')

    tmp = etree.SubElement(x, 'distributionLikelihood', id="location.diffusion.prior")
    tmp2 = etree.SubElement(tmp, 'data')
    etree.SubElement(tmp2, 'parameter', idref='location.diffusion.rates')

    tmp2 = etree.SubElement(tmp, 'distribution')
    tmp3 = etree.SubElement(tmp2, 'onePGammaDistributionModel')
    tmp4 = etree.SubElement(tmp3, 'shape')
    etree.SubElement(tmp4, 'parameter', value='0.5')

    return x


def write_cauchyrrwlikelihood_block(x, parameter):
    tmp = etree.SubElement(x, 'multivariateTraitLikelihood',
                           id='location.traitLikelihood',
                           traitName='location',
                           useTreeLength='true',
                           scaleByTime='true',
                           reportAsMultivariate='true',
                           reciprocalRates='true',
                           integrateInternalTraits='true')

    etree.SubElement(tmp, 'multivariateDiffusionModel', idref="location.diffusionModel")
    etree.SubElement(tmp, 'treeModel', idref="treeModel")

    tmp2 = etree.SubElement(tmp, 'traitParameter')
    etree.SubElement(tmp2, 'parameter', id='leaf.location')

    if parameter['continuous_phylogeo_jitter']:
        jitter = str(parameter['continuous_phylogeo_jitter'])
        tmp2 = etree.SubElement(tmp, 'jitter', window=jitter+' '+jitter, duplicatesOnly='true')
        etree.SubElement(tmp2, 'parameter', idref='leaf.location')

    tmp2 = etree.SubElement(tmp, 'conjugateRootPrior')
    tmp3 = etree.SubElement(tmp2, 'meanParameter')
    etree.SubElement(tmp3, 'parameter', value="0.0 0.0")

    tmp3 = etree.SubElement(tmp2, 'priorSampleSize')
    etree.SubElement(tmp3, 'parameter', value="0.000001")

    etree.SubElement(tmp, 'arbitraryBranchRates', idref="location.diffusion.branchRates")

    return x


def write_multivariatestats_block(x):
    tmp = etree.SubElement(x, 'correlation', id='location.correlation', dimension1='1', dimension2='2')
    etree.SubElement(tmp, "matrixParameter", idref="location.precision")

    tmp = etree.SubElement(x, 'matrixInverse', id='location.varCovar')
    etree.SubElement(tmp, "matrixParameter", idref="location.precision")

    tmp = etree.SubElement(x, 'continuousDiffusionStatistic', id="location.diffusionRate", greatCircleDistance="true")
    etree.SubElement(tmp, "matrixParameter", idref="location.traitLikelihood")

    return x
