//
// Reshape the VEP plugin releases that ship as zip archives (REVEL and EVE)
//

include { UNZIP as UNZIP_REVEL } from '../../modules/nf-core/unzip/main.nf'
include { UNZIP as UNZIP_EVE   } from '../../modules/nf-core/unzip/main.nf'
include { WGET as WGET_REVEL   } from '../../modules/nf-core/wget/main'
include { WGET as WGET_EVE     } from '../../modules/nf-core/wget/main'
include { VEPPLUGIN_REVEL      } from '../../modules/local/vepplugin/revel/main.nf'
include { VEPPLUGIN_EVE        } from '../../modules/local/vepplugin/eve/main.nf'

include { resolveVepPlugins    } from './utils_nfcore_lrsomatic_pipeline'

workflow PREPARE_VEP_PLUGINS {

    main:

    // Resolved here rather than taken as an input: a workflow input would arrive in a channel
    def plugins = resolveVepPlugins()
    def prepare = plugins.prepare

    ch_versions = channel.empty()
    def staged = [ channel.fromList(plugins.ready_files) ]

    //
    // MODULES: WGET_REVEL -> UNZIP_REVEL -> VEPPLUGIN_REVEL (labels: process_single, process_single, process_medium)
    // Input:  the REVEL release, as a URL or as a local zip
    // Output: .files -- revel_grch38.tsv.gz and its index
    // A remote release goes through wget: REVEL's host answers 403 to a request carrying no
    // User-Agent and EVE's redirects HTTPS to HTTP, so neither can be staged as a path input
    //
    if (prepare.containsKey('vep_revel')) {
        if (prepare['vep_revel'].toString().contains('://')) {
            WGET_REVEL (
                channel.value([ [ id: 'revel' ], prepare['vep_revel'] ])
            )

            ch_revel_zip = WGET_REVEL.out.outfile
            ch_versions = ch_versions.mix(WGET_REVEL.out.versions)
        }
        else {
            ch_revel_zip = channel.value([ [ id: 'revel' ], file(prepare['vep_revel'], checkIfExists: true) ])
        }

        UNZIP_REVEL (
            ch_revel_zip
        )

        VEPPLUGIN_REVEL (
            UNZIP_REVEL.out.unzipped_archive.map { _meta, dir -> dir }
        )

        staged << VEPPLUGIN_REVEL.out.files
        ch_versions = ch_versions.mix(UNZIP_REVEL.out.versions)
    }

    //
    // MODULES: WGET_EVE -> UNZIP_EVE -> VEPPLUGIN_EVE (labels: process_single, process_single, process_medium)
    // Input:  the EVE release, as a URL or as a local zip -- one VCF per protein
    // Output: .files -- eve_merged.vcf.gz and its index
    //
    if (prepare.containsKey('vep_eve')) {
        if (prepare['vep_eve'].toString().contains('://')) {
            WGET_EVE (
                channel.value([ [ id: 'eve' ], prepare['vep_eve'] ])
            )

            ch_eve_zip = WGET_EVE.out.outfile
            ch_versions = ch_versions.mix(WGET_EVE.out.versions)
        }
        else {
            ch_eve_zip = channel.value([ [ id: 'eve' ], file(prepare['vep_eve'], checkIfExists: true) ])
        }

        UNZIP_EVE (
            ch_eve_zip
        )

        VEPPLUGIN_EVE (
            UNZIP_EVE.out.unzipped_archive.map { _meta, dir -> dir }
        )

        staged << VEPPLUGIN_EVE.out.files
        ch_versions = ch_versions.mix(UNZIP_EVE.out.versions)
    }

    // A value channel, since both the germline and the somatic VEP task read it. collect() emits
    // nothing on an empty upstream, so ifEmpty is what carries the no-plugins case.
    ch_extra_files = staged
        .inject(channel.empty()) { acc, ch -> acc.mix(ch) }
        .flatten()
        .collect()
        .ifEmpty([])

    emit:
    extra_files = ch_extra_files // channel: value list of plugin .pm and data files
    versions    = ch_versions    // channel: versions.yml files
}
