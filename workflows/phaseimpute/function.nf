def chunkPrepareChannel(ch_chunks, tool) {
    def ch_chunks_branched = ch_chunks.branch{ meta, txt ->
            txt: txt != []
                return [meta, txt]
            empty: true
                error "ERROR: Empty chunks file provided for ${tool}"
        }
    if(tool == "glimpse1"){
        return ch_chunks_branched.txt.map { chr, txt -> [chr, file(txt)]}
            .splitCsv(header: ['ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'], sep: "\t", skip: 0)
            .map { meta, it -> [meta, it["RegionIn"], it["RegionOut"]]}
    } else if(tool == "quilt") {
        def ch_chunks_txt = ch_chunks_branched.txt.map { metaC, txt -> [metaC, file(txt)]}
            .splitCsv(header: ['ID', 'Chr', 'RegionIn', 'RegionOut', 'Size1', 'Size2'], sep: "\t", skip: 0)
            .map { metaC, it ->
                def startEnd = it["RegionIn"].split(':')[1].split('-')
                [ metaC, metaC.chr, startEnd[0], startEnd[1] ]
            }
        return ch_chunks_txt
    } else {
        error "ERROR: Only 'glimpse1' and 'quilt' output format are supported. Got ${tool}"
    }
}
