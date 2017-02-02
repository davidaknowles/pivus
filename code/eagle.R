#!/srv/gsfs0/software/R/3.3.1/bin/Rscript

library(abind)
require(data.table)
require(doMC)

#registerDoMC(detectCores()-1)
registerDoMC(10)

source("/home/dak33/Dropbox/eagle/eagle/beta_binomial_models/bb_glm_flips_rep_prior_gene.R")

#dat=fread("zcat < link2data/RNAseq/STAR_G-T/allelic_counts/ase.txt.gz")
dat=fread("zcat < link2data/DATA/RNAseq/Filtered_bam/ase.Aligned.out_mapq255_sorted.txt.gz")
setDF(dat)

gene_snps=fread("zcat < link2data/DATA/RNAseq/Filtered_bam/gene_snp_mapping.Aligned.out_mapq255_sorted.txt.gz")
setDF(gene_snps)

dat$snp=with(dat, paste(chr,pos,sep=":"))

require(reshape2)

ref=dcast(dat, snp ~ sample, value.var = "r")
rownames(ref)=ref$snp
ref$snp=NULL
ref=as.matrix(ref)

alt=dcast(dat, snp ~ sample, value.var = "y")
rownames(alt)=alt$snp
alt$snp=NULL
alt=as.matrix(alt)

stopifnot(all(colnames(alt)==colnames(ref)))
stopifnot(all(rownames(alt)==rownames(ref)))

alt[is.na(alt)]=0
ref[is.na(ref)]=0

n=alt+ref

samp_meta=read.table("~/pivus/link2data/DATA/Pheno_Covariates/expCovariates_eQTL.txt", stringsAsFactors = F)

samp_meta$samp=(seq_len(126)+1) %/% 2

samp_meta=samp_meta[ colnames(alt), c("Age","samp") ]

colnames(samp_meta)[1]="timepoint"

snpmeta=rownames(alt)

is_het=stan_model("~/Dropbox/eagle/eagle/beta_binomial_models/is_het.stan")

filter_data=function(gene) {

    snps=gene_snps[gene_snps$gene==gene,"snp"]
    i=snpmeta %in% snps

    cast_me=function(mat) do.call(abind,c( apply(mat, 1, function(g) {
        temp=cbind(samp_meta, a=g)
        #d=reshape2::dcast(temp, samp ~ timepoint , value.var="a")
        #rownames(d)=d$samp
        #d$samp=NULL
        #d
        unstack( temp, samp ~ timepoint )
    }), along=3 ))
        
    a=cast_me( alt[i,,drop=F] )
    a[is.na(a)]=0

    nh=cast_me( n[i,,drop=F] )
    nh[is.na(nh)]=0 # need this to account for missing samples

    dummy=foreach(snp_index=seq_len(dim(nh)[3])) %do% {
        as=a[,,snp_index]
        nhh=nh[,,snp_index]
        ind_to_keep=rowSums(nhh)>0
        treat_to_keep=colSums(nhh)>0

        as=as[ind_to_keep,treat_to_keep,drop=F]
        nhh=nhh[ind_to_keep,treat_to_keep,drop=F]

        o=  optimizing(is_het, dat=list(N=nrow(nhh), T=ncol(nhh), errorRate=0.01, concShape=1.001, concRate=0.001, ys=as, ns=nhh), as_vector=F)
        eo=exp(o$par$probs)
        pr=sweep(eo, 1, rowSums(eo), "/")

        homo=pr[,1]<0.95

        cat("Removing",sum(homo),"individual(s) from SNP",dimnames(a)[[3]][snp_index],"\n")

        nh[which(ind_to_keep)[homo],treat_to_keep,snp_index]=0
        a[which(ind_to_keep)[homo],treat_to_keep,snp_index]=0
    }

    snp_to_keep=apply(nh>0, 3, any)
    a=a[,,snp_to_keep,drop=F]
    nh=nh[,,snp_to_keep,drop=F]

    list(a=a,nh=nh)
}

genes=unique(gene_snps$gene)

p=foreach(gene=genes, .combine=c) %dopar% {

    filtered_data=filter_data(gene)
    if (is.null(filtered_data)) return(NA)

#    tryCatch({
    res <- betabinomial_glm_flips_rep_gene(filtered_data$a,filtered_data$nh,verbose=F, iterations = 5000)
    res$lrtp[1]
#} , error=function(g) NA )
}

gzf=gzfile("eagle_gene_new.txt.gz","w")
write.table(data.frame(gene=genes,p=p), file=gzf, row.names = F, quote = F, sep="\t")
close(gzf)
cat("Done!\n")
stop()

require(ggplot2)

#res=read.table("../data/eagle_gene.txt.gz",sep="\t",header=T, stringsAsFactors = F)

res=read.table("../data/eagle_gene_1conc.txt.gz",sep="\t",header=T, stringsAsFactors = F)


rownames(snpmeta)=snpmeta$variantID

require(dplyr)
gene_id=gene_snps %>% group_by(gene) %>% summarize(id=do.call(paste,as.list(snp)))
class(gene_id)="data.frame"
back_track=gene_id %>% group_by(id) %>% summarise(gene=do.call(paste,as.list(gene)))
class(back_track)="data.frame"
rownames(back_track)=back_track$id
gene_id$ambiguous=back_track[gene_id$id, "gene"]

unique_gene=gene_id[!duplicated(gene_id$id),]
rownames(res)=res$gene
res=res[unique_gene$gene,]
res$ambiguous=unique_gene$ambiguous

res$q=bh(res$p)
nsig=sum(res$q<.1, na.rm=T)
multiqq(list(gene=res$p))



pdf("../figures/eagle_gene1.pdf",width=12,height=10)
foreach(i=order(res$p)[1:nsig]) %do% {
      gene=res$gene[i]
        dat=filter_data(gene)
        melted=melt(dat$a/dat$nh)
        melted_n=melt(dat$nh)
        tokeep=!is.na(melted$value)
        melted=cbind(melted[tokeep,], n=melted_n[tokeep,"value"])
        colnames(melted)=c("Individual","TimePoint","SNP","AllelicRatio","Coverage")
        melted$TimePoint=as.numeric(melted$TimePoint)

        the_title=paste0(res$ambiguous[i]," p=",format(res$p[i],digits=1)," q=",format(res$q[i],digits=1))

        #levels(melted$SNP)=snpmeta[levels(melted$SNP),"position"]
        melted$inter=interaction( melted$TimePoint, melted$SNP )
        melted$ii=as.numeric(melted$inter)
        print( ggplot(melted, aes( inter, AllelicRatio, label=Individual, col=SNP, shape=Individual)) + geom_point( aes( size=Coverage)) + ylim(0,1) + theme_bw(base_size = 16) + geom_line(aes(ii,AllelicRatio),alpha=.5)  + scale_shape_manual(values=seq_along(levels(melted$Individual))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Treatment.SNP") + ggtitle(the_title) )

        NULL
  }
dev.off()
