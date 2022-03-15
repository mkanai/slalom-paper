import argparse
import hail as hl

from gnomad.utils.vep import (
    process_consequences,
    filter_vep_to_canonical_transcripts,
    get_most_severe_consequence_for_summary,
    CSQ_CODING_HIGH_IMPACT,
    CSQ_CODING_MEDIUM_IMPACT,
    CSQ_CODING_LOW_IMPACT,
    CSQ_NON_CODING,
)
from hail.genetics import reference_genome
from fm_insights.utils import register_log, annotate_bed

coding_high = hl.set(CSQ_CODING_HIGH_IMPACT)
coding_medium = hl.set(CSQ_CODING_MEDIUM_IMPACT)
coding_low = hl.set(CSQ_CODING_LOW_IMPACT)
non_coding = hl.set(CSQ_NON_CODING)

bed_files = {
    "GRCH37": [
        "gs://finemapping-insights/annotations/baselineLD_v2.2/Promoter_UCSC.bed",
        "gs://finemapping-insights/annotations/Ulirsch_v1.0/DHSmerged_Ulirsch.bed",
        "gs://finemapping-insights/annotations/Ulirsch_v1.0/Roadmap_H3K27ac_Ulirsch.bed",
        "gs://finemapping-insights/annotations/Ulirsch_v1.0/CA_H3K27ac_Ulirsch.bed",
    ],
    "GRCh38": [
        "gs://meta-finemapping-simulation/annotations_hg38/Promoter_UCSC.bed",
        "gs://meta-finemapping-simulation/annotations_hg38/DHSmerged_Ulirsch.bed",
        "gs://meta-finemapping-simulation/annotations_hg38/Roadmap_H3K27ac_Ulirsch.bed",
        "gs://meta-finemapping-simulation/annotations_hg38/CA_H3K27ac_Ulirsch.bed",
    ],
}


gnomad_latest_versions = {"GRCh37": "2.1.1", "GRCh38": "3.1.2"}
gnomad_v2_pops = ["afr", "amr", "asj", "eas", "fin", "nfe", "nfe_est", "nfe_nwe", "nfe_onf", "nfe_seu"]
gnomad_v3_pops = ["afr", "ami", "amr", "asj", "eas", "mid", "fin", "nfe", "oth", "sas"]


def annotate_consequence_category(csq_expr, annot_location="consequence_category"):
    annot_expr = {
        annot_location: hl.case()
        .when(coding_high.contains(csq_expr), "coding_high")
        .when(coding_medium.contains(csq_expr), "coding_medium")
        .when(coding_low.contains(csq_expr), "coding_low")
        .when(non_coding.contains(csq_expr), "non_coding")
        .or_missing()
    }
    return annot_expr


def main(args):
    reference_genome = args.reference_genome

    if reference_genome == "GRCh37":
        from gnomad.resources.grch37.gnomad import public_release

        ht = public_release("genomes").versions[gnomad_latest_versions[reference_genome]].ht()
        freq_index_dict = ht.freq_index_dict.collect()[0]
        freq_expr = {pop: ht.freq[freq_index_dict[f"gnomad_{pop}"]] for pop in gnomad_v2_pops}
        freq_expr.update({"all": ht.freq[freq_index_dict[f"gnomad"]]})
    elif reference_genome == "GRCh38":
        from gnomad.resources.grch38.gnomad import public_release

        ht = public_release("genomes").versions[gnomad_latest_versions[reference_genome]].ht()
        freq_index_dict = ht.freq_index_dict.collect()[0]
        freq_expr = {pop: ht.freq[freq_index_dict[f"{pop}-adj"]] for pop in gnomad_v3_pops}
        freq_expr.update({"all": ht.freq[freq_index_dict[f"adj"]]})
    else:
        raise ValueError("Invalid --reference-genome")

    ht = ht.annotate(freq=hl.struct(**freq_expr))
    ht = filter_vep_to_canonical_transcripts(ht)
    ht = process_consequences(ht)
    ht = get_most_severe_consequence_for_summary(ht)

    # extract most severe
    ht = ht.select(
        freq=ht.freq,
        most_severe=hl.if_else(hl.is_defined(ht.most_severe_csq), ht.most_severe_csq, "intergenic_variant"),
        gene_most_severe=ht.vep.worst_csq_for_variant_canonical.gene_symbol,
        lof=ht.vep.worst_csq_for_variant_canonical.lof,
        hgnc_id=ht.vep.worst_csq_for_variant_canonical.hgnc_id,
        hgvsp=ht.vep.worst_csq_for_variant_canonical.hgvsp,
        transcript_id=ht.vep.worst_csq_for_variant_canonical.transcript_id,
        polyphen_prediction=ht.vep.worst_csq_for_variant_canonical.polyphen_prediction,
        polyphen_score=ht.vep.worst_csq_for_variant_canonical.polyphen_score,
        sift_prediction=ht.vep.worst_csq_for_variant_canonical.sift_prediction,
        sift_score=ht.vep.worst_csq_for_variant_canonical.sift_score,
        protein_coding=ht.protein_coding,
    )

    ht = ht.select_globals()
    ht = ht.annotate(**annotate_consequence_category(ht.most_severe))
    ht = annotate_bed(ht, bed_files=bed_files[reference_genome], reference_genome=reference_genome)
    ht = ht.annotate(
        consequence=(
            hl.case(missing_false=True)
            .when(hl.is_defined(ht.lof) & (ht.lof != "LC"), "pLoF")
            .when(
                (ht.lof == "LC")
                | (ht.consequence_category == "coding_high")
                | (ht.consequence_category == "coding_medium"),
                "Missense",
            )
            .when(ht.consequence_category == "coding_low", "Synonymous")
            .when(ht.most_severe == "3_prime_UTR_variant", "UTR3")
            .when(ht.most_severe == "5_prime_UTR_variant", "UTR5")
            .when(ht.Promoter_UCSC == 1, "Promoter")
            .when(
                (ht.DHSmerged_Ulirsch == 1) & ((ht.Roadmap_H3K27ac_Ulirsch == 1) | (ht.CA_H3K27ac_Ulirsch == 1)), "CRE"
            )
            .default("Non-genic")
        )
    )
    ht.describe()
    ht = ht.checkpoint(
        f"gs://meta-finemapping-simulation/gnomad/gnomad.genomes.r{gnomad_latest_versions[args.reference_genome]}.sites.most_severe.ht",
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-genome", type=str, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    register_log()

    main(args)
