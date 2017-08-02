import os.path
import functools

import utils


# {1: {'gene_identity': 0, 'gene_mean_read_coverage': 0.0, 'gene_number_positions_multiple_alleles': 0, 'header': 'fyuA', 'gene_coverage': 0.0, 'gene_low_coverage': 100.0}}


def simplify_data_by_gene(data_by_gene):
    cleaned_data_by_gene = {}
    for counter, data in data_by_gene.items():
        cleaned_data_by_gene[data['header']] = {'gene_identity': data['gene_identity'], 'gene_coverage': data['gene_coverage'], 'gene_depth': data['gene_mean_read_coverage']}
    return cleaned_data_by_gene


def possible_types(data_by_gene, typing_rules_file, min_gene_coverage, min_gene_identity, min_gene_depth):
    data_by_gene = simplify_data_by_gene(data_by_gene)

    possible_pathotypes = []
    with open(typing_rules_file, 'rtU') as reader:
        genes = []
        for line in reader:
            line = line.splitlines()[0]
            if len(line) > 0:
                line = line.split('\t')
                if line[0].startswith('#'):
                    genes = map(str.lower, line[1:])
                else:
                    profile = line[1:]
                    congruence = []
                    for x, gene_requirement in enumerate(profile):
                        gene_requirement = True if gene_requirement == '1' else False if gene_requirement == '0' else None
                        if gene_requirement is None:
                            congruence.append(True)
                        else:
                            if data_by_gene[genes[x]]['gene_coverage'] >= min_gene_coverage and data_by_gene[genes[x]]['gene_identity'] >= min_gene_identity and data_by_gene[genes[x]]['gene_depth'] >= min_gene_depth:
                                gene_present = True
                            else:
                                gene_present = False

                            if gene_present == gene_requirement:
                                congruence.append(True)
                            else:
                                congruence.append(False)
                    if all(congruence):
                        possible_pathotypes.append(line[0])
    return possible_pathotypes


module_timer = functools.partial(utils.timer, name='Module Typing')


@module_timer
def typing(data_by_gene, typing_rules_file, min_gene_coverage, min_gene_identity, min_gene_depth, outdir):
    possible_pathotypes = possible_types(data_by_gene, typing_rules_file, min_gene_coverage, min_gene_identity, min_gene_depth)
    with open(os.path.join(outdir, 'patho_typing.report.txt'), 'wt') as writer:
        if len(possible_pathotypes) > 0:
            writer.write('\n'.join(possible_pathotypes) + '\n')
            print '\n' + 'Pathotypes found:' + '\n'
            print '\n'.join(possible_pathotypes) + '\n'
        else:
            writer.write('NA' + '\n')
            print '\n' + 'It was not possible to identify any possible pathotype match' + '\n'

    return None, None
