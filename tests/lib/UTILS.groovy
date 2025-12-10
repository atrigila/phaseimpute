// Helper functions for pipeline tests
class UTILS {
    public static def getRecursiveFileNames(fileOrDir, outputDir, variants_md5=true, txt_md5=true) {
        def f = new File(fileOrDir.toString())
        if ( f.isDirectory() ) {
            return f.listFiles()
                .collect { getRecursiveFileNames(it, outputDir, variants_md5, txt_md5) }
                .flatten()
                .sort { it.file }
        }
        def check = null
        def file_name = f.toString().replace("${outputDir}/","")
        if ( (f.name.endsWith('.txt') || f.name.endsWith('.log')) && txt_md5 ) {
            check = path(f.toString()).md5
            return ["file": file_name, "md5": check ]
        } else if ( f.name.endsWith('.vcf.gz') && variants_md5 ) {
            check = path(f.toString()).vcf.variantsMD5
            return ["file": file_name, "variants": check ]
        } else{
            return ["file": file_name]
        }
    }

    public static def vcfDetails(filePath) {
        def summary = path(filePath).vcf.summary.replaceAll(", phasedAutodetect=(false|true)", "")
        def samples = path(filePath).vcf.header.getGenotypeSamples().sort()
        return [summary: summary, samples: samples]
    }
}
