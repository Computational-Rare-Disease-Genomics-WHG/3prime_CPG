"""
filter_hail.py

This script filters a Hail context file to only include CpG variants on MANE transcripts.

Running this script:
    python3 filter_hail.py --hail-context <input_file> --output <output_file> --mane-summary <mane_summary_file>
"""

from typing import Union
import os
import sys
import argparse
import hail as hl


def flip_base(base: hl.expr.StringExpression) -> hl.expr.StringExpression:
    """
    Flip the single base
    """
    return hl.switch(base).when('A', 'T').when('T', 'A').when('G', 'C').when('C', 'G').default(base)


def reverse_complement_bases(bases: hl.expr.StringExpression) -> hl.expr.StringExpression:
    """
    Find the reverse complement bases
    """
    return hl.delimit(hl.range(bases.length() - 1, -1, -1).map(lambda i: flip_base(bases[i])), '')


def collapse_strand(ht: Union[hl.Table, hl.MatrixTable]) -> hl.Table:
    """
    Collapse a strand ht
    """
    collapse_expr = {
        'ref': hl.if_else(((ht.ref == 'G') | (ht.ref == 'T')), reverse_complement_bases(ht.ref), ht.ref),
        'alt': hl.if_else(((ht.ref == 'G') | (ht.ref == 'T')), reverse_complement_bases(ht.alt), ht.alt),
        'context': hl.if_else(
            ((ht.ref == 'G') | (ht.ref == 'T')), reverse_complement_bases(ht.context), ht.context
        ),
        'was_flipped': (ht.ref == 'G') | (ht.ref == 'T'),
    }
    return ht.annotate(**collapse_expr)

def annotate_variant_types(
    t: Union[hl.MatrixTable, hl.Table], heptamers: bool = False
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Adds cpg, transition, and variant_type, variant_type_model columns
    """
    mid_index = 3 if heptamers else 1
    transition_expr = (
        ((t.ref == 'A') & (t.alt == 'G'))
        | ((t.ref == 'G') & (t.alt == 'A'))
        | ((t.ref == 'T') & (t.alt == 'C'))
        | ((t.ref == 'C') & (t.alt == 'T'))
    )
    cpg_expr = (
        (t.ref == 'G') & (t.alt == 'A') & (t.context[mid_index - 1 : mid_index] == 'C')
    ) | ((t.ref == 'C') & (t.alt == 'T') & (t.context[mid_index + 1 : mid_index + 2] == 'G'))
    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(transition=transition_expr, cpg=cpg_expr)
    else:
        t = t.annotate(transition=transition_expr, cpg=cpg_expr)
    variant_type_expr = hl.case().when(t.cpg, 'CpG').when(t.transition, 'non-CpG transition').default('transversion')
    variant_type_model_expr = hl.if_else(t.cpg, t.context, 'non-CpG')

    return t.annotate(
            variant_type=variant_type_expr,
            variant_type_model=variant_type_model_expr
    )


def trimer_from_heptamer(t: Union[hl.MatrixTable, hl.Table]) -> hl.Table:
    """
    Context is heptameric, function provides the trimer (i.e. filter to a 1 bp context)
    """
    trimer_expr = hl.if_else(hl.len(t.context) == 7, t.context[2:5], t.context)
    return t.annotate(context=trimer_expr)


def annotate_cpgs(ht: hl.Table) -> hl.Table:
    """
    Annotate CpG sites in a Hail table
    @param ht: Hail table
    @return: Hail table with CpG annotations
    """

    # Set context length
    str_len = 3

    # Collapse the strand
    ht = collapse_strand(ht)

    ht = ht.filter((hl.len(ht.ref) == 1) & (hl.len(ht.alt) == 1) & ht.context.matches(f'[ATCG]{{{str_len}}}'))

    # Annotate variant types
    ht = annotate_variant_types(ht, heptamers=False)

    # Filter to CpG sites
    ht = ht.filter(ht.cpg)

    # Return the Hail table
    return ht


def filter_to_mane(
    mt: hl.Table, mane_ht: hl.Table, vep_root: str = 'vep'
) -> hl.Table:
    """
    Filters the vep output to the MANE select transcript.
    Similar to gnomad.utils.vep.filter_vep_to_canonical_transcripts
    from the gnomAD package.

    @param ht : gnomad hail table to filter
    @param mane_ht : Mane summary data
    @param vep_root : prefix for the vep columns in the gnomad file
    @returns mt : a hail table where the transcript consequences
                are only on Mane Select for each variant
    """
    mane_list = hl.literal(mane_ht.ensembl_transcript_stable_id.collect())
    mane_select = mt[vep_root].transcript_consequences.filter(lambda csq: mane_list.contains(csq.transcript_id))
    vep_data = mt[vep_root].annotate(transcript_consequences=mane_select)
    return mt.annotate_rows(**{vep_root: vep_data})


def main(args):
    """Main entry point """

    # Check if hail file exists
    if not os.path.exists(args.hail_context):
        print("Hail context file not found")
        sys.exit(1)

    # Initialize hail
    hl.init(log='/tmp/hail.log', default_reference='GRCh38')

    # Import hail context table
    ht = hl.read_table(args.hail_context)

    # Annotate CpG sites
    ht = annotate_cpgs(ht)

    # Filter items not on MANE transcripts
    mane_ht = hl.import_table(args.mane_summary)
    # Remove version number
    mane_ht = mane_ht.annotate(ensembl_transcript_stable_id=mane_ht.ensembl_transcript_stable_id.split(".")[0])
    ht = filter_to_mane(ht, mane_ht)

    # Write out hail context table as a TSV file
    ht.export(args.output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This script filters Hail")
    parser.add_argument('--hail-context', help='Hail context file', required=True)
    parser.add_argument('--output', help='Output file', required=True)
    parser.add_argument('--mane-summary', help='MANE summary file', required=True)
    args = parser.parse_args()
    main(args)
