import argparse
import textwrap

import pysam

# =============================================================================

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
Analysis of MHC prediction results in combination with SNPs in a particular protein sequence.
------------------------------------------
'''))
# Parse command line arguments
# -i FILE -t TASK -h HELP -o OUTPUT -th Threshold

parser.add_argument("-v", "--vcf", metavar='file.vcf.gz', dest="vcf", help="VCF file", type=str)
parser.add_argument("-o", "--output", metavar='file.pdf', dest="output", help="Output file to save", type=str)

args = parser.parse_args()


# =============================================================================

def main():
    vcf_reader = pysam.VariantFile(str(args.vcf), 'r')

    data = []
    for record in vcf_reader:
        data.append([record.chrom, record.pos, record.start, record.stop, record.id, record.ref, list(record.alts),
                     ])
