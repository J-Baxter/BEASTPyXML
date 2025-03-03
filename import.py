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

def write_alignment_block(x, taxa, sequences):
    aln_block = etree.SubElement(x, 'alignment', id='alignment', dataType='nucleotide')

    n = len(taxa)
    for i in range(n):
        seq = etree.SubElement(aln_block, "sequence")
        etree.SubElement(seq,"taxon",idref=taxa[i])
        seq.text = str(sequences[i])

    return x


<treeModel id="treeModel">
		<coalescentTree idref="startingTree"/>
		<rootHeight>
			<parameter id="treeModel.rootHeight"/>
		</rootHeight>
		<nodeHeights internalNodes="true">
			<parameter id="treeModel.internalNodeHeights"/>
		</nodeHeights>
		<nodeHeights internalNodes="true" rootNode="true">
			<parameter id="treeModel.allInternalNodeHeights"/>
		</nodeHeights>
	</treeModel>


def write_treemodel_block(x, taxa, precision, tree_model):
    treemodel_block = etree.SubElement(x, 'treeModel', id='treeModel')
    etree.SubElement(treemodel_block, "rootHeight", idref="startingTree")
    etree.SubElement(treemodel_block, "nodeHeights", idref="startingTree")
    etree.SubElement(treemodel_block, "nodeHeights", idref="startingTree")


    return x

taxa, seq_list = read_fasta('USUV_2025Feb10_alldata_aligned_formatted_noFLI_NFLG.fasta_subsample1.fasta')
dates = parse_dates(taxa)
date_decimal, precision = format_dates(dates)

root = etree.Element('beast', version='1.0')
test_3 = write_taxa_block(root, taxa, date_decimal, precision)
test_4 = write_alignment_block(test_3, taxa, seq_list)
xml_string = etree.tostring(test_4, pretty_print=True, encoding="utf-8").decode()

print(xml_string)