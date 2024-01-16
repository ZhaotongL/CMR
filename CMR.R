library(readr)
library(data.table)
library(dplyr)
library(TwoSampleMR)
library(MendelianRandomization)
library(optparse)

option_list = list(
  make_option("--exp_gwas", action="store", default=NA, type='character',
              help="Path to exposure GWAS summary statistics (COJO format + chr + pos) [required]"),
  make_option("--out_gwas", action="store", default=NA, type='character',
              help="Path to outcome GWAS summary statistics (COJO format + chr + pos) [required]"),
  make_option("--bfile", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam) for LD reference panel splitted into 22 chromosomes [required]"),
  make_option("--pfile", action="store", default=NA, type='character',
              help="Path to PLINK2 pgen input file prefix (minus pgen/psam/pvar)for LD reference panel splitted into 22 chromosomes [required if --bfile is not provided]"),
  make_option("--LDref_SNP", action="store", default=NA, type='character',
              help="Path to a text file containing all SNPs present in the LD reference panel [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--tmp", action="store", default=NA, type='character',
              help="Path to temporary files [required]"),
  make_option("--p_threshold", action="store", default=5e-8, type='double',
              help="P-value threshold with the exposure to select IVs [default: %default]"),
  make_option("--COJOp_threshold", action="store", default=5e-6, type='double',
              help="Conditional/joint P-value threshold to select SNPs in COJO [default: %default]"),
  make_option("--PLINK2", action="store", default="plink2", type='character',
              help="Path to plink2 executable [%default]"),
  make_option("--GCTA", action="store", default="gcta_nr_robust", type='character',
              help="Path to gcta executable [%default]"),
  make_option("--ldblock", action="store", default="fourier_ls-all.bed", type='character',
              help="Path to LD Block file [%default]"),
  make_option("--MR_methods", action="store", default="IVW,Median,Mode", type='character',
              help="MR methods to be run, seperated by ',' [%default]"),	
  make_option("--standardMR", action="store_true", default=TRUE,
              help="Run standard MR analysis or not. [default: %default]"),		  
  make_option("--save.data", action="store_true", default=TRUE,
              help="Save the IVs data for running MR and CMR or not. [default: %default]")          
)

opt = parse_args(OptionParser(option_list=option_list))
GCTA = opt$GCTA
PLINK2 = opt$PLINK2

LDref_SNP = fread(opt$LDref_SNP)


temp_dir <- opt$tmp

exp_df_cojo = fread(opt$exp_gwas)
out_df_cojo = fread(opt$out_gwas)
colnames(exp_df_cojo) = c('SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N','chr','pos')
colnames(out_df_cojo) = c('SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N','chr','pos')
common_SNP = intersect(exp_df_cojo$SNP,out_df_cojo$SNP)
common_SNP = intersect(common_SNP,LDref_SNP)
exp_df_cojo = subset(exp_df_cojo, exp_df_cojo$SNP %in% common_SNP)
out_df_cojo = subset(out_df_cojo, out_df_cojo$SNP %in% exp_df_cojo$SNP)
gwasX = sprintf("%s/exposure.cojo", temp_dir)
fwrite(subset(exp_df_cojo, select=c('SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N')),gwasX,sep='\t')

exp_dat = exp_df_cojo %>% select(rsid = SNP, pval=p) %>% filter(pval < opt$p_threshold) 
clump_dat = NULL
for(chr in 1:22){
    clump_dat_tmp = ld_clump_local(dat=exp_dat,bfile=sprintf('%s_chr%s',opt$bfile,chr)) 
    clump_dat = rbind(clump_dat,clump_dat_tmp)
}
exp_mr_dat = subset(exp_df_cojo,exp_df_cojo[['SNP']] %in% clump_dat[['rsid']])


    
gwasY = sprintf("%s/outcome.cojo", temp_dir)
fwrite(subset(out_df_cojo, select=c('SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N')),gwasY,sep='\t')

out_mr_dat = subset(out_df_cojo,out_df_cojo[['SNP']] %in% clump_dat[['rsid']])
out_mr_dat = out_mr_dat[match(exp_mr_dat$SNP,out_mr_dat$SNP),]


out_mr_dat_2smr = out_mr_dat %>% select(SNP,beta.outcome=b,se.outcome=se,effect_allele.outcome=A1,other_allele.outcome=A2,eaf.outcome=freq,pval.outcome=p,samplesize.outcome=N,chr,pos)
exp_mr_dat_2smr = exp_mr_dat %>% select(SNP,beta.exposure=b,se.exposure=se,effect_allele.exposure=A1,other_allele.exposure=A2,eaf.exposure=freq,pval.exposure=p,samplesize.exposure=N,chr,pos)
exp_mr_dat_2smr$id.exposure=1
exp_mr_dat_2smr$exposure='exposure'
out_mr_dat_2smr$id.outcome=2
out_mr_dat_2smr$outcome='outcome'
harmo_dat = harmonise_data(exp_mr_dat_2smr,out_mr_dat_2smr,2)
harmo_dat = harmo_dat %>% filter(mr_keep==TRUE)
exp_mr_dat = harmo_dat %>% select(SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure,samplesize.exposure,chr,pos)
out_mr_dat = harmo_dat %>% select(SNP,effect_allele.outcome,other_allele.outcome,eaf.outcome,beta.outcome,se.outcome,pval.outcome,samplesize.outcome,chr,pos)
colnames(exp_mr_dat) = colnames(out_mr_dat) = c('SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N','chr','pos')
N_x = median(exp_mr_dat$N)
N_y = median(out_mr_dat$N)

exp_mr_dat$chr = as.numeric(exp_mr_dat$chr)
exp_mr_dat$pos = as.numeric(exp_mr_dat$pos)


## Marginal statistics - standard MR
all_SNP = exp_mr_dat %>% arrange(chr,pos)
b_exp = exp_mr_dat$b[match(all_SNP$SNP,exp_mr_dat$SNP)]
se_exp = exp_mr_dat$se[match(all_SNP$SNP,exp_mr_dat$SNP)]
pval_exp = exp_mr_dat$p[match(all_SNP$SNP,exp_mr_dat$SNP)]
b_out = out_mr_dat$b[match(all_SNP$SNP,out_mr_dat$SNP)]
se_out = out_mr_dat$se[match(all_SNP$SNP,out_mr_dat$SNP)]
pval_out = out_mr_dat$p[match(all_SNP$SNP,out_mr_dat$SNP)]

## Conditional statistics - CMR Stage 1
b_exp_cond = se_exp_cond = pval_exp_cond = b_out_cond = se_out_cond = pval_out_cond = NULL
ldBlocks = read.table(opt$ldblock,header=TRUE)
ldBlocks$chr = gsub(pattern='chr',replacement='',x=ldBlocks$chr)
ldBlocks_res = NULL
    
for(i in 1:nrow(ldBlocks)){
    curChr = ldBlocks[i,1]
    curMin = ldBlocks[i,2]
    curMax = ldBlocks[i,3]
    curSnps = all_SNP[all_SNP$pos > curMin & all_SNP$pos <= curMax & all_SNP$chr == curChr,]
    if(nrow(curSnps) > 0){
        ldBlocks_res = rbind(ldBlocks_res,c(ldBlocks[i,],nrow(curSnps)))
        curLDFILE = sprintf('%s_chr%s', opt$bfile, curChr)
        bfile = sprintf('%s/block%s', temp_dir, i)
        if(!file.exists(paste0(bfile,'.bed'))){
                plink_command = sprintf("%s --pfile %s --chr %s --from-bp %s --to-bp %s --make-bed --out %s",
                                PLINK2, curLDFILE, curChr, curMin, curMax, bfile)
                plk_msg=system(plink_command,intern=verbose)
        }
        cojo_command = sprintf("%s --bfile %s --cojo-file %s --cojo-slct --cojo-p %s --out %s/out_block%s", GCTA, bfile, gwasY, opt$COJOp_threshold, temp_dir, i )
        cojo_msg=system(cojo_command,intern=verbose)
        exp_snp = curSnps$SNP
        if(file.exists(sprintf("%s/out_block%s.jma.cojo",temp_dir, i))){
            out_cojo = fread(sprintf("%s/out_block%s.jma.cojo",temp_dir,i))
            out_snp = out_cojo$SNP

            snp = sprintf("%s/block%s_snp.txt",temp_dir,i)
            write.table(union(exp_snp, out_snp), snp, quote=F,row.names=F,col.names=F)
            
            joint_cm = sprintf("%s --bfile %s --cojo-file %s --extract %s --cojo-joint --out %s/Xblock%s", GCTA, bfile, gwasX, snp, temp_dir,i)
            cojo_msg=system(joint_cm,intern=verbose)
            
            
            if(!file.exists(sprintf("%s/Xblock%s.jma.cojo",temp_dir,i))){
                Ysnp = sprintf("%s/block%s_Ysnp.txt",temp_dir,i)
                write.table(out_snp, Ysnp, quote=F,row.names=F,col.names=F)
                cond_cm = sprintf("%s --bfile %s --cojo-file %s --extract %s --cojo-cond %s --out %s/Xblock%s", GCTA, bfile, gwasX, snp, Ysnp,temp_dir,i)
                cojo_msg=system(cond_cm,intern=verbose)
                cond_cm = sprintf("%s --bfile %s --cojo-file %s --extract %s --cojo-cond %s --out %s/Yblock%s", GCTA, bfile, gwasY, snp, Ysnp,temp_dir,i)
                cojo_msg=system(cond_cm,intern=verbose)
                cond_cojo_X = fread(sprintf("%s/Xblock%s.cma.cojo",temp_dir,i))
                b_exp_cond = c(b_exp_cond, cond_cojo_X$bC[match(exp_snp,cond_cojo_X$SNP)])
                se_exp_cond = c(se_exp_cond, cond_cojo_X$bC_se[match(exp_snp,cond_cojo_X$SNP)])
                pval_exp_cond = c(pval_exp_cond, cond_cojo_X$pC[match(exp_snp,cond_cojo_X$SNP)])

                cond_cojo_Y = fread(sprintf("%s/Yblock%s.cma.cojo",temp_dir,i))
                b_out_cond = c(b_out_cond, cond_cojo_Y$bC[match(exp_snp,cond_cojo_Y$SNP)])
                se_out_cond = c(se_out_cond, cond_cojo_Y$bC_se[match(exp_snp,cond_cojo_Y$SNP)])
                pval_out_cond = c(pval_out_cond, cond_cojo_Y$pC[match(exp_snp,cond_cojo_Y$SNP)])

            }else{
                joint_cm = sprintf("%s --bfile %s --cojo-file %s --extract %s --cojo-joint --out %s/Yblock%s", GCTA, bfile, gwasY, snp, temp_dir,i)
                cojo_msg=system(joint_cm,intern=verbose)

                joint_cojo_X = fread(sprintf("%s/Xblock%s.jma.cojo",temp_dir,i))
                b_exp_cond = c(b_exp_cond, joint_cojo_X$bJ[match(exp_snp,joint_cojo_X$SNP)])
                se_exp_cond = c(se_exp_cond, joint_cojo_X$bJ_se[match(exp_snp,joint_cojo_X$SNP)])
                pval_exp_cond = c(pval_exp_cond, joint_cojo_X$pJ[match(exp_snp,joint_cojo_X$SNP)])

                joint_cojo_Y = fread(sprintf("%s/Yblock%s.jma.cojo",temp_dir,i))
                b_out_cond = c(b_out_cond, joint_cojo_Y$bJ[match(exp_snp,joint_cojo_Y$SNP)])
                se_out_cond = c(se_out_cond, joint_cojo_Y$bJ_se[match(exp_snp,joint_cojo_Y$SNP)])
                pval_out_cond = c(pval_out_cond, joint_cojo_Y$pJ[match(exp_snp,joint_cojo_Y$SNP)])
            }
        }else{
            b_exp_cond = c(b_exp_cond, harmo_dat$beta.exposure[match(exp_snp,harmo_dat$SNP)])
            se_exp_cond = c(se_exp_cond, harmo_dat$se.exposure[match(exp_snp,harmo_dat$SNP)])
            b_out_cond = c(b_out_cond, harmo_dat$beta.outcome[match(exp_snp,harmo_dat$SNP)])
            se_out_cond = c(se_out_cond, harmo_dat$se.outcome[match(exp_snp,harmo_dat$SNP)])
            pval_exp_cond = c(pval_exp_cond, harmo_dat$pval.exposure[match(exp_snp,harmo_dat$SNP)])
            pval_out_cond = c(pval_out_cond, harmo_dat$pval.outcome[match(exp_snp,harmo_dat$SNP)])
        }
        rm_cmd = sprintf('rm %s.*', bfile)
        system(rm_cmd)
        }
}

 
    if(opt$save.data){
        dat = list(SNP=all_SNP$SNP,
             b_exp = b_exp, b_out = b_out, se_exp = se_exp, se_out = se_out, pval_exp = pval_exp, pval_out = pval_out,
             b_exp_cond = b_exp_cond, b_out_cond = b_out_cond, se_exp_cond = se_exp_cond, se_out_cond = se_out_cond,
             pval_exp_cond = pval_exp_cond, pval_out_cond = pval_out_cond)
   }else{
        dat = list()
   }
    rm_cmd = sprintf('rm %s/*cojo*', temp_dir)
    system(rm_cmd)
}
unlink(temp_dir, recursive = TRUE)

MR = unique( c(unlist(strsplit(opt$MR_methods,',')),"top1") )
M = length(MR)

res_MR = NULL
## standard MR analysis
if(opt$standardMR){
    Marg_mr_obj = mr_input(bx= dat$b_exp, bxse = dat$se_exp, by = dat$b_out, byse = dat$se_out, snp = dat$SNP)
    for(i in 1:M){
        if(MR[i]=='IVW'){
            mr_res = mr_ivw(Marg_mr_obj)
            res_MR = rbind(res_MR, c('IVW',mr_res@Estimate,mr_res@StdError,mr_res@Pvalue))
        }
        if(MR[i]=='Median'){
            mr_res = mr_median(Marg_mr_obj)
            res_MR = rbind(res_MR, c('Median',mr_res@Estimate,mr_res@StdError,mr_res@Pvalue))
        }
        if(MR[i]=='Mode'){
            mr_res = mr_mbe(Marg_mr_obj)
            res_MR = rbind(res_MR, c('Mode',mr_res@Estimate,mr_res@StdError,mr_res@Pvalue))
        }
        if(MR[i]=='cML'){
            mr_res = mr_cML(Marg_mr_obj, n = min(N_x,N_y))
            res_MR = rbind(res_MR, c('cML',mr_res@Estimate,mr_res@StdError,mr_res@Pvalue))
        }
    }
    res_MR = as.data.frame(res_MR)
    colnames(res_MR) = c('Method','b','se','pval')
}

## MR analysis - CMR Stage 2
keep_ind_cond = which(!is.na(dat$b_exp_cond)) 
nsnp_cond = length(keep_ind_cond)
Joint_mr_obj = mr_input(bx= dat$b_exp_cond[keep_ind_cond], bxse = dat$se_exp_cond[keep_ind_cond], by = dat$b_out_cond[keep_ind_cond], byse = dat$se_out_cond[keep_ind_cond], snp = dat$SNP[keep_ind_cond])

res_CMR = NULL
for(i in 1:M){
    if(MR[i]=='IVW'){
        mr_res = mr_ivw(Joint_mr_obj)
        res_CMR = rbind(res_MR, c('IVW',mr_res@Estimate,mr_res@StdError,mr_res@Pvalue))
    }
    if(MR[i]=='Median'){
        mr_res = mr_median(Joint_mr_obj)
        res_CMR = rbind(res_MR, c('Median',mr_res@Estimate,mr_res@StdError,mr_res@Pvalue))
    }
    if(MR[i]=='Mode'){
        mr_res = mr_mbe(Joint_mr_obj)
        res_CMR = rbind(res_MR, c('Mode',mr_res@Estimate,mr_res@StdError,mr_res@Pvalue))
    }
    if(MR[i]=='cML'){
        mr_res = mr_cML(Joint_mr_obj, n = min(N_x,N_y))
        res_CMR = rbind(res_MR, c('cML',mr_res@Estimate,mr_res@StdError,mr_res@Pvalue))
    }
}
res_CMR = as.data.frame(res_MR)
colnames(res_CMR) = c('Method','b','se','pval')

dat$MR_res = res_MR
dat$CMR_res = res_CMR
saveRDS(dat,sprintf('%sCMR.rds',opt$out))