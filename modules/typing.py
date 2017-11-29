import os.path
import functools

import utils


# {1: {'gene_identity': 0, 'gene_mean_read_coverage': 0.0, 'gene_number_positions_multiple_alleles': 0, 'header': 'fyuA', 'gene_coverage': 0.0, 'gene_low_coverage': 100.0}}


def simplify_data_by_gene(data_by_gene):
    cleaned_data_by_gene = {}
    for counter, data in data_by_gene.items():
        cleaned_data_by_gene[data['header'].lower()] = {'gene_identity': data['gene_identity'], 'gene_coverage': data['gene_coverage'], 'gene_depth': data['gene_mean_read_coverage']}
    return cleaned_data_by_gene


def possible_types(data_by_gene, typing_rules_file, min_gene_coverage, min_gene_identity, min_gene_depth):
    data_by_gene = simplify_data_by_gene(data_by_gene)

    possible_pathotypes = []
    genes_present = []
    with open(typing_rules_file, 'rtU') as reader:
        genes = []
        for line in reader:
            line = line.splitlines()[0]
            if len(line) > 0:
                line = line.split('\t')
                if line[0].startswith('#'):
                    genes = line[1:]
                else:
                    profile = line[1:]
                    congruence = []
                    for x, gene_requirement in enumerate(profile):
                        if data_by_gene[genes[x].lower()]['gene_coverage'] >= min_gene_coverage and data_by_gene[genes[x].lower()]['gene_identity'] >= min_gene_identity and data_by_gene[genes[x].lower()]['gene_depth'] >= min_gene_depth:
                            gene_present = True
                            genes_present.append(genes[x])
                        else:
                            gene_present = False

                        gene_requirement = True if gene_requirement == '1' else False if gene_requirement == '0' else None
                        if gene_requirement is not None:
                            if gene_present == gene_requirement:
                                congruence.append(True)
                            else:
                                congruence.append(False)
                        else:
                            congruence.append(True)
                    if all(congruence):
                        possible_pathotypes.append(line[0])
    return possible_pathotypes, list(set(genes_present))


module_timer = functools.partial(utils.timer, name='Module Typing')


@module_timer
def typing(data_by_gene, typing_rules_file, min_gene_coverage, min_gene_identity, min_gene_depth, outdir):
    possible_pathotypes, genes_present = possible_types(data_by_gene, typing_rules_file, min_gene_coverage, min_gene_identity, min_gene_depth)
    with open(os.path.join(outdir, 'patho_typing.report.txt'), 'wt') as writer_report:
        with open(os.path.join(outdir, 'patho_typing.extended_report.txt'), 'wt') as writer_extended_report:
            writer_extended_report.write('#Pathotypes_found' + '\n')
            if len(possible_pathotypes) > 0:
                writer_report.write('\n'.join(possible_pathotypes) + '\n')
                writer_extended_report.write('\n'.join(possible_pathotypes) + '\n')
                print '\n' + 'Pathotypes found:' + '\n'
                print '\n'.join(possible_pathotypes) + '\n'
            else:
                writer_report.write('NA' + '\n')
                writer_extended_report.write('NA' + '\n')
                print '\n' + 'It was not possible to identify any possible pathotype match' + '\n'
            writer_extended_report.write('\n' + '#Genes_present' + '\n')
            writer_extended_report.write('\n'.join(genes_present) + '\n')

    return None, None
