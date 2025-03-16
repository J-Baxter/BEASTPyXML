# General datatype block
from lxml import etree

def write_generaldatatype_block(x, trait, trait_name):
    # Comment
    comment = etree.Comment('General data type for discrete trait model,' + trait_name)
    x.insert(1, comment)

    tmp = etree.SubElement(x,'generalDataType', id=trait_name+'.dataType')

    # one line for each state to set up CTM
    unique_states = set(trait.values())
    comment = etree.Comment('Number of States = ' + len(unique_states))
    tmp.insert(1, comment)

    for unique_state in unique_states:
        etree.SubElement(tmp, 'state', code=unique_state)

    return x

def write_attributepatterns_block(x, trait_name):
    # Comment
    comment = etree.Comment('Data pattern for discrete trait,' + trait_name)
    x.insert(1, comment)

    tmp = etree.SubElement(x, 'attributePatterns', id=trait_name+".pattern", attribute=trait_name)
    etree.SubElement(tmp, 'taxa', idref="taxa")
    etree.SubElement(tmp, 'generalDataType', idref=trait_name+".dataType")

    return x

    # take the strict clock (+ rate statistic) from other blocks


# Currently written for symmetric model. Need to understand how dimension is calculated
# and changes required for asymmentric model.
def write_generalsubstitutionmodel_block(x, trait, trait_name, model):
    comment = etree.Comment('symmetric CTMC model for discrete state reconstructions')
    x.insert(1, comment)
    tmp = etree.SubElement(x, 'generalSubstitutionModel', id=trait_name+".model", randomizeIndicator="false")
    etree.SubElement(tmp, 'generalDataType', idref=trait_name+".dataType")

    # frequency sub-block
    tmp2 = etree.SubElement(tmp, 'frequencies')
    tmp3 = etree.SubElement(tmp2, 'frequencyModel',  id=trait_name+".frequencyModel", normalize="true")
    etree.SubElement(tmp3, 'generalDataType', idref=trait_name+'.dataType')
    tmp4 = etree.SubElement(tmp3, 'frequencies')

    unique_states = set(trait.values())
    etree.SubElement(tmp4, 'parameter', id=trait_name+".frequencies", dimension=len(unique_states))

    # rates sub-block
    comment = etree.Comment('rates and indicators')
    tmp.insert(10000000, comment)

    tmp2 = etree.SubElement(tmp, 'rates')
    etree.SubElement(tmp2, "parameter", id=trait_name+".rates", dimension=trait_dimension, value="1.0", lower="0.0")

    # site sub-block
    tmp = etree.SubElement(x, 'siteModel', id=trait_name+".siteModel")
    tmp2 = etree.SubElement(tmp, 'substitutionModel')
    etree.SubElement(tmp2, 'generalSubstitutionModel', idref=trait_name+".model")

    return x



#{x + '.dataType' for x in set(thisdict.values())}

