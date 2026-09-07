// IMPORT MODULES
include { SAVANA_RUN      } from '../../../modules/nf-core/savana/run/main'
include { SAVANA_CLASSIFY } from '../../../modules/nf-core/savana/classify/main'
include { SAVANA_CNA      } from '../../../modules/nf-core/savana/cna/main'

workflow PAIRED_SAVANA {

    take:
    tumor_normal_input // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_germline_vcf, phased_germline_tbi]
    fasta              // [[:], fasta]
    fai                // [[:], fai]
    contigs             // [[:], contigs]  -- one contig per line, restricts SAVANA to canonical chromosomes

    main:
    // .first() re-establishes a value channel: fasta/fai are singletons reused for every sample,
    // but .map() on fai strips that property, so plain .combine() would only fire once for the run.
    ch_fasta_fai = fasta.combine(fai.map { _meta, index -> index }).first()
    // ch_fasta_fai: [[:], fasta, fai]

    //
    // MODULE: SAVANA_RUN (label: process_high)
    // Input:  [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
    //         [[:], fasta, fai]
    // Output: .sv_breakpoints_vcf -- [meta, vcf]  -- raw (unclassified) SV breakpoints
    //
    tumor_normal_input
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, _phased_vcf, _phased_tbi ->
            return [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
        }
        .set { savana_run_input }
    // savana_run_input: [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]

    SAVANA_RUN (
        savana_run_input,
        ch_fasta_fai,
        contigs
    )

    //
    // MODULE: SAVANA_CLASSIFY (label: process_medium)
    // Input:  [meta, sv_breakpoints_vcf]  -- SAVANA_RUN's raw breakpoints
    // Output: .somatic_vcf -- [meta, vcf]  -- classified somatic SV VCF, optional (empty if no somatic SVs)
    //
    SAVANA_CLASSIFY (
        SAVANA_RUN.out.sv_breakpoints_vcf
    )

    //
    // MODULE: SAVANA_CNA (label: process_high)
    // SAVANA's own combined `savana to`/`savana_main` workflow feeds CNA the *classified* somatic
    // breakpoints (args.breakpoints = args.somatic_output), not the raw ones -- classification
    // filters out germline/artefact junctions before they inform CNA binning.
    // Input:  [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, snp_vcf, [], breakpoints]
    //         [[:], fasta, fai]
    //         contigs / blacklist / g1000_vcf ([]/[]/[]: snp_vcf takes priority in paired mode)
    //
    tumor_normal_input
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_vcf, _phased_tbi ->
            def allele_counts_het_snps = []
            return [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_vcf, allele_counts_het_snps]
        }
        // somatic_vcf is optional -- fall back to [] (no breakpoints) rather than dropping the sample
        .join(SAVANA_CLASSIFY.out.somatic_vcf, remainder: true)
        .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_vcf, allele_counts_het_snps, somatic_vcf ->
            return [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_vcf, allele_counts_het_snps, somatic_vcf ?: []]
        }
        .set { savana_cna_input }
    // savana_cna_input: [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, snp_vcf, [], breakpoints]

    SAVANA_CNA (
        savana_cna_input,
        ch_fasta_fai,
        [],  // contigs -- SAVANA_CNA's own module doesn't take a genome-aware contigs file; test profile scopes it via ext.args instead
        [],  // blacklist (unused)
        []   // g1000_vcf (unused -- snp_vcf takes priority in paired mode)
    )

    emit:
    somatic_vcf = SAVANA_CLASSIFY.out.somatic_vcf  // [meta, vcf]  -- classified somatic SV VCF
    cn_calls    = SAVANA_CNA.out.cna               // [meta, tsv]  -- segmented absolute copy number
}
