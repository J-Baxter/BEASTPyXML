from lxml import etree


def write_empirical_tree_model(x, parameters):
    comment = etree.Comment('Create Empirical Tree Model')
    x.insert(1000, comment)
    tmp = etree.SubElement(x,
                           'empiricalTreeDistributionModel', id='treeModel',
                           fileName=parameters.empirical_tree_distribution)
    etree.SubElement(tmp, 'taxa', idref='taxa')

    tmp = etree.SubElement(x, 'statistic', id="treeModel.currentTree", name="Current Tree")

    etree.SubElement(tmp, 'empiricalTreeDistributionModel', idref='treeModel')

    return x


def write_empiricaltree_operator(x):
    comment = etree.Comment('Empirical Tree Operator')
    x.insert(1000, comment)
    tmp = etree.SubElement(x, 'empiricalTreeDistributionOperator', weight="3")
    etree.SubElement(tmp, 'empiricalTreeDistributionModel', idref='treeModel')

    return x
