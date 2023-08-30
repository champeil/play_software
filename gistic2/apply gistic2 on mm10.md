# create gistic2 reference mat file for mm10
## preview
- the gistic2 only got mat files for human (hg16-hg38)
- when applying the co-occured cytobands analysis to another species, we need to create the reference file
  - so we use matlab to create the reference files

## prepare and start
- matlab (in gistic2 the version of the matlab is R2014a)
- gistic2
- hg19.mat file
- python and R to modify the necessary files

## first: analysis hg19.mat
- use matlab to load hg19.mat `load('hg19.mat')`
  - three variables in the workspace
    - cyto: struct format to record the cytobands info
      - chr: chromosome name without 'chr' prefix (1:22,X,Y) in char class
      - chrn: the order of the chromosome name (1:22,23,24) in double class
      - name: the names of cytobands (1q36.33) in char class
      - start: the start of cytobands in double class
      - end: the end of cytobands in double class
      - stain: the result of Giemsa stain in char class
    - rg: to record the detail transcript's or ncrna's annotations
      - refseq: ref_id (only use NM_ and NR_) in char class
      - gene: gene name or discription in char class
      - symb: gene symbol in char class
      - locus_id: unique gene record, there use the entrez_id in int_32 class
      - chr: chromosome name of gene with 'chr' prefix in char class
      - strand: strand of gene in char class
      - start: the start of the ref_id in int_32 class
      - end: the end of the ref_id in int_32 class
      - cds_start: the cds start of ref_id in int_32 class
      - cds_end: the cds end of ref_id in int_32 class
      - status: annotation status (for example: reviewed) in char class
      - chrn: the order of the chromosome name (1:22,23,24) in double class
        - the field must in double class or else will cause error of different class (other type vs double)
    - rg_info: metadata of the reference file, just edit it causually
      - two column: Field and Value
        - Field: assembly, source, path, url, script, makeversion, maker, data
        - Value: the value of corresponding Field

## then create mm10 reference mat file and recompile
### mm10 mat file creation 
#### cyto table
- source: ucsc, just douwnload and use as cyto table
- [http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/cytoBand.txt.gz](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/cytoBand.txt.gz)
#### rg table
- method1: gtf files from refGene, Ensembl and GENCODE (have no try)
- method2: [Table Browser (ucsc.edu) | UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1682280320_PyXA25KMEfQIMXupGWI6kZ0Os7bu&clade=mammal&org=Mouse&db=mm10&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=refGene&hgta_regionType=genome&position=chr12%3A56%2C694%2C976-56%2C714%2C605&hgta_outputType=primaryTable&hgta_outFileName=result2.txt)
  - parameter
    - clase: mammal
    - genome: mouse
    - assembly: GRCm38/mm10
    - group: genes and gene predictions
    - track: NCBI refseq
    - table: refseq all
  - result: obtain the table with refseq, cds_start, cds_end, gene info of all mrna and ncrna
  - then use biomaRt package to annotate and get the rest information
#### rg_info
- just do it in the specified format
#### combine as mat file 
### code
```R
# this script is for creating files for RG and RG_INFO
# author: laojp
# time: 2023.08.28
# position: SYSUCC bioinformatic platform

# first read the necessary packages
packages <- c("dplyr","stringr","tidyverse","biomaRt","rtracklayer")
for (i in packages) {
  if (!require(i, character.only = TRUE)) {
    stop(paste("package: ", i, " is not exits"))
  }
}

# biomart for converting the refid
listMarts(host="https://nov2020.archive.ensembl.org/")
mm10_mart <- useMart("ensembl",host="https://nov2020.archive.ensembl.org/")
listDatasets(mm10_mart)
mm10_mart <- useDataset("mmusculus_gene_ensembl", mm10_mart)
listFilters(mm10_mart)

# get NM and NR, convert and combine
refseq <- read.table("C:/Users/Administrator/Desktop/result (1).txt",header=TRUE) %>%
  dplyr::select(name,cdsStart,cdsEnd) %>%
  dplyr::filter(stringr::str_detect(name,pattern="^NM")) %>%
  dplyr::mutate(name=stringr::str_remove_all(name,pattern="\\..*")) %>%
  dplyr::distinct() %>%
  dplyr::left_join(biomaRt::getBM(attributes=c("refseq_mrna","description","external_gene_name","chromosome_name",
                                               "transcript_start","transcript_end","entrezgene_id",
                                               "strand"),filter="refseq_mrna",
                                  values=.$name[grep("^NM",.$name)] %>% unlist() %>% as.character(),mart=mm10_mart),
                   by=c("name"="refseq_mrna")) %>%
  dplyr::bind_rows(read.table("C:/Users/Administrator/Desktop/result (1).txt",header=TRUE) %>%
                     dplyr::select(name,cdsStart,cdsEnd) %>%
                     dplyr::filter(stringr::str_detect(name,pattern="^NR")) %>%
                     dplyr::mutate(name=stringr::str_remove_all(name,pattern="\\..*")) %>%
                     dplyr::distinct() %>%
                     dplyr::left_join(biomaRt::getBM(attributes=c("refseq_ncrna","description","external_gene_name","chromosome_name",
                                                                  "transcript_start","transcript_end","entrezgene_id",
                                                                  "strand"),filter="refseq_ncrna",
                                                     values=.$name[grep("^NR",.$name)] %>% unlist() %>% as.character(),mart=mm10_mart),
                                      by=c("name"="refseq_ncrna"))) %>%
  dplyr::filter(!is.na(description)) %>%
  dplyr::filter(!is.na(external_gene_name)) %>%
  dplyr::filter(chromosome_name %in% c(1:19,"X","Y")) %>%
  dplyr::filter(!is.na(entrezgene_id)) %>%
  dplyr::rename(refseq=name,
                gene=description,
                symb=external_gene_name,
                locus_id=entrezgene_id,
                chr=chromosome_name,
                strand=strand,
                start=transcript_start,
                end=transcript_end,
                cds_start=cdsStart,
                cds_end=cdsEnd) %>%
  dplyr::mutate(chrn=as.numeric(unlist(lapply(chr,function(x){ifelse(x=="X",20,ifelse(x=="Y",21,x))}))),
                chr=paste("chr",chr,sep = ""),
                status="Reviewed",
                strand=unlist(lapply(strand,function(x){ifelse(x==-1,0,1)}))) %>%
  dplyr::select(refseq,gene,symb,locus_id,chr,strand,start,end,cds_start,cds_end,status,chrn)
                            
write.table(refseq,"C:/Users/Administrator/Desktop/rg.txt",sep="\t",quote = FALSE,row.names = FALSE)
```
```matlab
% first read cyto table
% ps: because chr column is a mixture of numeric and character, so we need to set options
opt=detectImportOptions("C:\Users\Administrator\Desktop\cytoBand.csv")
opt=setvartype(opt,{'chr'},'char')
opt=setvartype(opt,{'start'},'int32')
opt=setvartype(opt,{'xEnd'},'int32')
cyto=readtable("C:\Users\Administrator\Desktop\cytoBand.csv",opt)
cyto = sortrows(cyto,"chr","ascend")
cyto = table2struct(cyto)

% then read eg table
opt=detectImportOptions("C:\Users\Administrator\Desktop\rg.txt")
opt=setvartype(opt,{'chr'},'char')
opt=setvartype(opt,{'start'},'int32')
opt=setvartype(opt,{'xEnd'},'int32')
opt=setvartype(opt,{'cds_start'},'int32')
opt=setvartype(opt,{'cds_end'},'int32')
opt=setvartype(opt,{'chrn'},'double')
rg=readtable("C:\Users\Administrator\Desktop\rg.txt",opt) 
rg = table2struct(rg)

% then read rg_info table
rg_info=struct(assembly=char("GRCm38/mm10"),source=char("UCSC"),path=char("https://genome.ucsc.edu/"),url=char(""),script=char(""),make_version=char(""),maker=char("champeil"),date=char("21-Aug-2023"))

% save mat file
save C:\Users\Administrator\Desktop\mm10.mat cyto rg rg_info
```
- ps: because the 'end' name is conflict with the parameter of readtable in matlab R2023a, so the 'end' name will change to 'X_End' automatically and we need to change it manually with python
```python
import scipy.io
import numpy as np
mat_data=scipy.io.loadmat("/home/laojp/software/gistic_mm10/refgenefiles/mm10.mat")
mat_data['cyto'].dtype.names=('chr', 'chrn', 'name', 'start', 'end', 'stain')
mat_data['rg'].dtype.names=('refseq','gene','symb','locus_id','chr','strand','start','end','cds_start','cds_end','status','chrn')
scipy.io.savemat("/home/laojp/software/gistic_mm10/refgenefiles/mm10_2.mat",mat_data)
```
```matlab
% then i notice the dimention num of hg19.mat is reversed and the class of rg.chrn of mm10_2.mat is unit8, so we need to modified again in matlab
cyto=transpose(cyto)
rg=transpose(rg)
for i = 1:size(rg,1) rg(i).chrn=double(rg(i).chrn); end
save C:\Users\Administrator\Desktop\mm10.mat cyto rg rg_info
```

# recompile gp_gistic2_from_seg execute file
## preview
- as the reference said, when using mat file of other species, we need to re-compile the source dir to get gp_gistic2_from_seg execute file
  - for some scripts is coded with human chromosome as reference (hard-coding chromosome?)
- the organization form of mouse chromosome is different with human's
  - because the centrosomes of mouse chromosomes near the edge, so only q band in mouse cytoband files
- so we need to modify the source code and recompile with matlab

## process
### modify the source code file
#### modify source/RefGeneInfo.m 
```matlab
% line 15: change the total number of chromosome
nchr = 21;
% line 29: change the chromosome name
RGI.chr.symb = {'1','2','3','4','5','6','7','8','9','10',...'11','12','13','14','15','16','17','18',...'19','X','Y'};
% line 40-42: change the autosomal region
RGI.txt2num('20') = 20;   
RGI.txt2num('21') = 21;    
RGI.chr.autosomal = (1:21)<20; 
```
#### modify source/normalize_by_arm_length.m
```matlab
% line 85, 86: because each chromosome has two bands (p,q) in human and only one in mouse, so cancel the '2*' to ignore the "index exceed dimensions error"
isq = find(D.pos(Q(:,2))>=armstart(Q(:,1))')
spans_cent = intersect(find(D.pos(Q(:,2))<armstart(Q(:,1))'), ...
  find(D.pos(Q(:,3))>=armstart(Q(:,1))'));
% line 114, 115, 120, 126, 127, 132: the same reasons to ignore the "index exceed dimensions error"
norms(isq) = armlengths_by_snp(Q(isq,1)');
norms(isp) = armlengths_by_snp((Q(isp,1)-1)');
norms(spans_cent) = armlengths_by_snp(QQ(spans_cent,1)');
norms(isq) = armlengths(Q(isq,1)');
norms(isp) = armlengths((Q(isp,1)-1)');
norms(spans_cent) = armlengths(QQ(spans_cent,1)');
```
#### modify source/arm_medians.m and source/add_cyto.m (linked by D.armn variant)
```matlab
% add_cyto.m line 52-56
if ~isempty(find(cyto(p).name=='p'))
  C.armn(snp)=2;
else
  C.armn(snp)=1;
end
% arm_medians.m line 30-33、37
arm_medians = nan(Nchr,size(D.dat,2));
arm_marks = zeros(Nchr,1);                                                                                                   
arm_names = repmat({''},Nchr,1);                                                                                             
pq = 'q';
for a=1
% ps: when use other type of cytoband, then modify the code according to the organization of the cytoband
```
### recompile with matlab to get new module
- move to the dir of source/ dir
- enter matlab and input the compile code (according to the github of broad gistic2)
  - if use linux in gistic2, then need to recompile with matlab in linux system
  - I try to recompile in windows but not work
- replace the new gp_gistic2_from_seg with the old one in [gistic2_dir]/gp_gistic2_from_seg  
```matlab 
% gp_gistic2_from_seg cannot be change for the name of source/gp_gistic2_from_seg.m
mcc -v -m -w enable -I [absolute_path_to_source]/source gp_gistic2_from_seg
```
## about MEX-file invalid error
- error detail
  - Invalid MEX-file '/home/laojp/.mcrCache8.3/gp_gis1/toolbox/optim/optim/private/barrierConvexQPmex.mexa64': libmwactiveset.so: cannot open shared object file: No such file or directory
  - Error in ipqpcommon>ipConvexQP (line 98)
- reasons: without corresponding lib file in [gistic2_dir]/MATLAB_RUN_TIME/bin/glnxa64
- solve: copy corresponding lib files from the installation dir of matlab (bin/glnxa64) to [gistic2_dir]/MATLAB_RUN_TIME/bin/glnxa64

## about "arm 1: 1q 19236 markers" output
- cause gistic2 is designed for human, D.armn variant only got {1,2} for {p,q} respectively
  - if not modify,then will output "arm 1: 1p 0 markers" and it is wired obviously
  - when apply to other species which has {p,q,r} band, it also only output {p,q} but not {r}
  - so we need to move back to gistic_broad_analysis.m code and move back to add_cyto.m and arm_medians.m and modify the codes
  - the details are in the "modify the source code file" parts

# fulture (a big pie in the sky)
- modify all code and create a GISTIC2 which can apply to all species

# ps
- the Chinese version of this tutorial is in my obsidian md file, tell me if you want
- the related files are uploaded to the play_software/gistic2/

# citation
- Zhang M, Wen H, Liang M, Qin Y, Zeng Q, Luo D, Zhong X, Li S. Diagnostic Value of Sylvian Fissure Hyperechogenicity in Fetal SAH. AJNR Am J Neuroradiol. 2022 Apr;43(4):627-632. doi: 10.3174/ajnr.A7449. Epub 2022 Mar 10. PMID: 35272984; PMCID: PMC8993207.

# reference
- [how to build refgene.mat file · Issue #5 · sbamin/canine_gistic2 (github.com)](https://github.com/sbamin/canine_gistic2/issues/5)
- [sbamin/canine_gistic2: Canine CanFam3.1 GISTIC2 copy number analysis. For debugging GISTIC2 .mat file, please contact developers at https://www.genepattern.org/modules/docs/GISTIC_2.0) or https://groups.google.com/a/broadinstitute.org/g/gistic-forum (github.com)](https://github.com/sbamin/canine_gistic2)
- [Cytoband | Integrative Genomics Viewer (broadinstitute.org)](https://software.broadinstitute.org/software/igv/Cytoband)
- [Custom cytoBand GenomicRanges (for mouse genome) · Issue #10 · dozmorovlab/multiHiCcompare (github.com)](https://github.com/dozmorovlab/multiHiCcompare/issues/10)
- [sbamin/canine_gistic2: Canine CanFam3.1 GISTIC2 copy number analysis. For debugging GISTIC2 .mat file, please contact developers at https://www.genepattern.org/modules/docs/GISTIC_2.0) or https://groups.google.com/a/broadinstitute.org/g/gistic-forum (github.com)](https://github.com/sbamin/canine_gistic2)
- [gistic2/support/release_gistic_source.sh at master · broadinstitute/gistic2 (github.com)](https://github.com/broadinstitute/gistic2/blob/master/support/release_gistic_source.sh#L36)
- [Where can I download the length of short and long arms for each chromosome (biostars.org)](https://www.biostars.org/p/383786/)





       








  
