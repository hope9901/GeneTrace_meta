"""Conservative interpretation of observed annotation membership."""


def scope_aware_presence(observed_copies, compartments, assembly_scope):
    """Keep positive observations; expose missing assay scope as N/A, not loss.

    A returned zero is still only unobserved annotation membership. It is never
    a validated biological-loss call. Mixed/unknown compartments are not guessed.
    """
    if observed_copies > 0:
        return 1
    if compartments == {'nuclear_or_unclassified'} and assembly_scope.get('nuclear_coding_loci') == 0:
        return 'N/A'
    if compartments == {'organelle'} and assembly_scope.get('explicit_organelle_bp') == 0:
        return 'N/A'
    if compartments == {'plasmid'} and assembly_scope.get('explicit_plasmid_bp') == 0:
        return 'N/A'
    return 0
