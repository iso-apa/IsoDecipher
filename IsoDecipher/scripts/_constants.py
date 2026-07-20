CODING_BIOTYPES = {
    "protein_coding",
    "IG_C_gene", "IG_V_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene",
    "TR_C_gene", "TR_V_gene", "TR_D_gene", "TR_J_gene",
}

IG_BIOTYPES = {
    "IG_C_gene", "IG_V_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene",
}

TR_BIOTYPES = {
    "TR_C_gene", "TR_V_gene", "TR_D_gene", "TR_J_gene",
}

IG_PSEUDO_BIOTYPES = {
    "IG_C_pseudogene", "IG_V_pseudogene",
    "IG_D_pseudogene", "IG_J_pseudogene",
}

TR_PSEUDO_BIOTYPES = {
    "TR_C_pseudogene", "TR_V_pseudogene",
    "TR_D_pseudogene", "TR_J_pseudogene",
}

PAS_STRONG   = {'AATAAA', 'ATTAAA'}
PAS_MODERATE = {
    'AGTAAA', 'TATAAA', 'CATAAA', 'GATAAA',
    'AATATA', 'AATACA', 'AATAGA', 'AAGAAA', 'ACTAAA',
}


def get_pas_motif_strength(motif_str, panel_source):
    """
    Returns: 'strong' | 'moderate' | 'weak' | 'unknown_gtf' | 'unknown_pa2'
    motif_str: raw PAS_motif field (e.g. 'AATAAA@5' or 'AATAAA@5;ATTAAA@12')
    """
    if not isinstance(motif_str, str):
        return 'unknown_gtf' if panel_source == 'GTF' else 'unknown_pa2'
    first = motif_str.split(';')[0]
    motif = first.split('@')[0]
    if motif in PAS_STRONG:
        return 'strong'
    elif motif in PAS_MODERATE:
        return 'moderate'
    else:
        return 'weak'
