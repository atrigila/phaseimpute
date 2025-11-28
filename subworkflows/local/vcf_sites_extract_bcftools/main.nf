include { BCFTOOLS_CONVERT              } from '../../../modules/nf-core/bcftools/convert'
include { BCFTOOLS_VIEW                 } from '../../../modules/nf-core/bcftools/view'
include { GAWK as GAWK_COMMA            } from '../../../modules/nf-core/gawk'
include { TABIX_BGZIP as BGZIP_COMMA    } from '../../../modules/nf-core/tabix/bgzip'
include { GAWK as GAWK_NOCOMMA          } from '../../../modules/nf-core/gawk'
include { TABIX_BGZIP as BGZIP_NOCOMMA  } from '../../../modules/nf-core/tabix/bgzip'

workflow VCF_SITES_EXTRACT_BCFTOOLS {
    take:
    ch_vcf          // channel: [ [id, chr], vcf, index ]
    ch_fasta        // channel: [ [genome], fasta, fai ]

    main:

    ch_versions = Channel.empty()
    ch_fasta = ch_fasta.map { meta, fasta, _fai -> [meta, fasta] }

    // Convert VCF to Hap and Legend files
    BCFTOOLS_CONVERT(ch_vcf, ch_fasta, [])
    ch_versions = ch_versions.mix(BCFTOOLS_CONVERT.out.versions.first())

    // Extract sites positions
    BCFTOOLS_VIEW(ch_vcf, [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())

    // Transform posfile to TSV
    GAWK_NOCOMMA(BCFTOOLS_CONVERT.out.legend, [], false)
    ch_versions = ch_versions.mix(GAWK_NOCOMMA.out.versions.first())

    // Transform posfile to TSV with ','
    GAWK_COMMA(BCFTOOLS_CONVERT.out.legend, [], false)
    ch_versions = ch_versions.mix(GAWK_COMMA.out.versions.first())

    // Compress TSV
    BGZIP_COMMA(GAWK_COMMA.out.output)
    ch_versions = ch_versions.mix(BGZIP_COMMA.out.versions.first())

    BGZIP_NOCOMMA(GAWK_NOCOMMA.out.output)
    ch_versions = ch_versions.mix(BGZIP_NOCOMMA.out.versions.first())

    // Join extracted sites and index
    ch_posfile = BCFTOOLS_VIEW.out.vcf
        .join(BCFTOOLS_VIEW.out.tbi)
        .join(BCFTOOLS_CONVERT.out.hap)
        .join(BCFTOOLS_CONVERT.out.legend)
        .join(BGZIP_COMMA.out.output)
        .join(BGZIP_NOCOMMA.out.output)

    emit:
    posfile       = ch_posfile          // channel: [ [id, chr], vcf, csi, hap, legend, posfile_comma, posfile_nocomma ]
    versions      = ch_versions         // channel: [ versions.yml ]
}
