require(doMC)
registerDoMC(4)
require(data.table)

basedir="/srv/gsfs0/projects/montgomery/bballiu/PIVUS/DATA/RNAseq/Filtered_bam/"

suffix=".Aligned.out_mapq255_sorted"

files=list.files(basedir, glob2rx(paste0("*",suffix,".bam")))

samples=sapply(strsplit(files, ".", fixed=T), function(g) g[1]) 

all_hets=foreach(sample=samples, .combine = c, .errorhandling = "stop") %do% {
    
    outfile=paste0(basedir,"allelic_counts/",sample,suffix,".txt.gz")
            cat(outfile)
            if (!file.exists(outfile)) {
                        cat("...doesn't exist\n")
                                return(NULL)
                    }
            cat("\n")
            b=fread(paste0("zcat < ",outfile))
            b=as.data.frame(b)
            colnames(b)=c("chr","pos","r","y","ref","alt")
            paste0(b$chr,":",b$pos)
    }

ta=table(all_hets)
to_keep=names(ta)[ta>1]

all_samp=rbindlist( foreach(sample=samples, .errorhandling = "remove") %dopar% {
    outfile=paste0(basedir,"allelic_counts/",sample,suffix,".txt.gz")
            cat(outfile)
            if (!file.exists(outfile)) {
                        cat("...doesn't exist\n")
                                return(NULL)
                    }
    b=fread(paste0("zcat < ",outfile))
    b=as.data.frame(b)
    colnames(b)=c("chr","pos","r","y","ref","alt")
    b$sample=sample
    b[ paste0(b$chr,":",b$pos) %in% to_keep, ]
} )

all_samp=as.data.frame(all_samp)

gzf=gzfile(paste0(basedir,"ase",suffix,".txt.gz"), "w")
write.table(all_samp, file=gzf, quote=F, row.names=F, col.names = T, sep="\t")
close(gzf)

# Get gene SNP mapping
exons=fread("zcat < ~/pivus/exons.saf.gz")
exons=as.data.frame(exons)

chroms=unique(exons$Chr)
exons_by_chrom=setNames(foreach(chrom=chroms) %dopar% exons[exons$Chr==chrom,],chroms)

genes=exons[,1:2]
genes=genes[!duplicated(genes$GeneID),]
rownames(genes)=genes$GeneID

snpmeta=all_samp[,c("chr","pos")]
snpmeta=snpmeta[!duplicated(snpmeta),]

snpmeta$variantID=with(snpmeta, paste(chr,pos,sep=":"))

snps_by_chrom=setNames(foreach(chrom=chroms) %dopar% snpmeta[snpmeta$chr==chrom,],chroms)

gene_snps=rbindlist( foreach(gene=genes$GeneID) %dopar% {
      chrom=genes[gene,"Chr"]
        exons_for_gene=exons_by_chrom[[chrom]][exons_by_chrom[[chrom]]$GeneID==gene,]
        snps=unique( foreach(i=seq_len(nrow(exons_for_gene)), .combine = c) %do% {
                snps_by_chrom[[chrom]]$variantID[ (exons_for_gene[i,"Start"] <= snps_by_chrom[[chrom]]$pos) & (exons_for_gene[i,"End"] >= snps_by_chrom[[chrom]]$pos) ]
            } )
        if (length(snps)==0) NULL else data.frame(gene=gene,snp=snps)
  })

gene_snps=as.data.frame(gene_snps)

gene_snps$gene=as.character(gene_snps$gene)
gene_snps$snp=as.character(gene_snps$snp)

gzf=gzfile(paste0(basedir, "gene_snp_mapping",suffix,".txt.gz"),"w")
write.table(gene_snps, file=gzf, col.names = T, row.names = F, quote = F, sep="\t")
close(gzf)

