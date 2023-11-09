"""
Extract from vcf-file columns
REF POS ALT FORMAT:FREQ
to tsv-file
"""

from os import path
import sys


def vcf_parse(input_path: str, output_path: str = None) -> None:
    """
    Extract from vcf-file columns
    REF POS ALT FORMAT:FREQ
    to tsv-file
    Default output_filename: input_filename_selected.tsv
    """
    vcf_results = []
    if output_path is None:
        file_name = '.'.join(path.basename(input_path).split('.')[:-1]) + '_selected.tsv'
        output_path = path.join(path.dirname(input_path), file_name)
    with open(input_path) as file:
        line = file.readline()
        while line.startswith('##'):
            line = file.readline()
        for line in file:
            info_all = line.strip().split('\t')
            info_last_col = info_all[9].split(':')
            info_select = (info_all[3], int(info_all[1]), info_all[4], float(info_last_col[6].strip('%')))
            vcf_results.append(info_select)
    with open(output_path, mode='w') as file:
        header = 'REF\tCOORD\tALT\tFREQ,%'
        file.write(header + '\n')
        for vcf_info in vcf_results:
            info_to_string = [str(i) for i in vcf_info]
            line_to_write = '\t'.join(info_to_string)
            file.write(line_to_write + '\n')


script_input = sys.argv[1]

script_output = None

if len(sys.argv) > 2:
    script_output = sys.argv[2]

vcf_parse(script_input, script_output)
