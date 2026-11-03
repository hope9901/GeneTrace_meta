"""Conservative guards for interpreting regional raw-read evidence."""

def quality_encoding_assessment(character_counts):
    """Observed FASTQ encoding alone cannot diagnose sequencing quality failure."""
    chars={c for c,n in character_counts.items() if n>0}
    zero_only=chars=={'!'}
    return dict(all_observed_qualities_encode_phred_zero=zero_only,
                meaning='quality_information_unresolved' if zero_only else 'encoded_quality_values_present' if chars else 'no_quality_observations',
                bq_filter_can_remove_all_observed_bases=zero_only,
                sequencing_failure_inferred=False,
                warning='All observed qualities are !. Check conversion/original files before using base-quality thresholds; these values alone do not prove poor sequencing.' if zero_only else '')

def regional_evidence_assessment(covered_fraction, spanning_mq20_records, partial_input, full_file_checksum_verified=False):
    if not 0<=covered_fraction<=1 or spanning_mq20_records<0:
        raise ValueError('Invalid evidence metrics')
    supported=covered_fraction==1 and spanning_mq20_records>0
    return dict(classification='read_supported_sequence_presence_candidate' if supported else 'unresolved_read_support',
                partial_input=partial_input,full_file_checksum_verified=full_file_checksum_verified,
                low_coverage_diagnosis_allowed=not partial_input and full_file_checksum_verified,
                biological_loss_inferred=False,functional_status='not_assessed',
                independent_molecule_support='not_established')
