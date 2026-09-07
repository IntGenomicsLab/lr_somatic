// IMPORT MODULES
include { SAVANA_TO } from '../../../modules/nf-core/savana/to/main'

workflow TUMORONLY_SAVANA {

    take:
    tumor_input // [meta, tumor_bam, tumor_bai, phased_germline_vcf, phased_germline_tbi]  -- tumor-only, no matched normal
    fasta       // [[:], fasta]
    fai         // [[:], fai]
    contigs     // [[:], contigs]  -- one contig per line, restricts SAVANA to canonical chromosomes
    g1000_vcf   // [[:], '1000g_hg38' | '1000g_t2t']  -- bundled 1000 Genomes population SNP set (genome-specific)

    main:
    // .first() re-establishes a value channel: fasta/fai are singletons reused for every sample,
    // but .map() on fai strips that property, so plain .combine() would only fire once for the run.
    ch_fasta_fai = fasta.combine(fai.map { _meta, index -> index }).first()
    // ch_fasta_fai: [[:], fasta, fai]

    //
    // MODULE: SAVANA_TO (label: process_high)
    // Chains breakpoint-calling + classification + CNA internally. Tumor-only has no matched
    // germline control, so allele counting uses the bundled 1000g population SNP set
    // (--g1000_vcf) rather than the tumour's own phased VCF, per SAVANA's README.
    // Input:  [meta, tumor_bam, tumor_bai, [], [], []]  -- snp_vcf/allele_counts_het_snps/breakpoints unused
    //         [[:], fasta, fai] / contigs / [[:],[]] blacklist / g1000_vcf
    // Output: .somatic_vcf -- [meta, vcf]  -- classified somatic SV VCF (savana classify, run internally)
    //         .cna         -- [meta, tsv]  -- segmented absolute copy number (savana cna, run internally)
    //
    tumor_input
        .map { meta, tumor_bam, tumor_bai, _phased_vcf, _phased_tbi ->
            def snp_vcf = []
            def allele_counts_het_snps = []
            def breakpoints = []
            return [meta, tumor_bam, tumor_bai, snp_vcf, allele_counts_het_snps, breakpoints]
        }
        .set { savana_to_input }
    // savana_to_input: [meta, tumor_bam, tumor_bai, [], [], []]

    SAVANA_TO (
        savana_to_input,
        ch_fasta_fai,
        contigs,
        [[:], []],  // blacklist (unused)
        g1000_vcf
    )

    emit:
    somatic_vcf = SAVANA_TO.out.somatic_vcf  // [meta, vcf]  -- classified somatic SV VCF
    cn_calls    = SAVANA_TO.out.cna          // [meta, tsv]  -- segmented absolute copy number
}
