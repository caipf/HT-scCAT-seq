workflow main{
	#ATAC
	Map[String,Array[String]] refConfig_a = {
		"hg19" : ["hg19","hs","/hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scATAC/database/hg19/bwa_mem2_index/genome.fa","/hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scATAC/database/hg19/regions/tss.bed","/hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scATAC/database/hg19/regions/promoter_scATAC.bed","chrM",""],
		"mm10" : ["mm10","mm","/hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scATAC/database/mm10/bwa_mem2_index/genome.fa","/hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scATAC/database/mm10/regions/tss.bed","/hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scATAC/database/mm10/regions/promoter_scATAC.bed","chrM","/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/mawen/pipeline/01.scATAC/01.scATAC_v3.0/database/chromap_index/mm10.index"],
		"hg38" : ["hg38","hs","/hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scATAC/database/hg38/bwa_mem2_index/genome.fa","/hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scATAC/database/hg38/regions/tss.bed","/hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scATAC/database/hg38/regions/promoter_scATAC.bed","chrM","/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/test/4_scATAC_chromap/index/genome.index"],
		"malu" : ["malu","2403862610","/hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scATAC/database/malu/bwa_mem2_index/malu_new.fa","/hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scATAC/database/malu/regions/TSS_malu.bed","/hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scATAC/database/malu/regions/Promoter_4k_malu.bed","chrMT",""],
		"macaque" : ["macaque","2871842584","/hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scATAC/database/macaque/bwa_mem2_index/Macaca_fascicularis_5.0.dna_sm.toplevel.fa","/hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scATAC/database/macaque/regions/TSS_Macaca_fascicularis.bed","/hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scATAC/database/macaque/regions/Promoter_4k_Macaca_fascicularis.bed","chrMT","/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/mawen/07.d2c_test/03.chromap_diff_library_test/index/macaca.genome.index"] 
	}
	String Fastq1_atac
	String Fastq2_atac
	String SampleID
	String ProjectID
	String reference
  	String? readStructure
  	String ID=SampleID
  	String PID = "sentieon.q -P "+ProjectID+"_sentieon"
	String refdir=refConfig_a[reference][2]
	String Rscript="/jdfssz1/ST_SUPERCELLS/PUB/scCAT/pipeline/V1/bin/Rscript"
    String python="/jdfssz1/ST_SUPERCELLS/PUB/scCAT/software/python"
	String tss=refConfig_a[reference][3]
	String promo=refConfig_a[reference][4]
	String chrmt = refConfig_a[reference][5]
	String ref_index = refConfig_a[reference][6]
	String runID= ID
	String? defaultConfig="/jdfssz1/ST_SUPERCELLS/PUB/scCAT/scripts/C4scATAClib_seqT1_R1_70_R2_50.json"
	String config=select_first([readStructure,defaultConfig])
	String species=refConfig_a[reference][0]
	String macs2sp=refConfig_a[reference][1]
	String whitelist="/jdfssz1/ST_SUPERCELLS/PUB/scCAT/scripts/whitelist.txt"
    String chrome_size="/jdfssz1/ST_SUPERCELLS/PUB/scCAT/software/bin/d2c_v1.4.4/bin/anno/bedtools/chrom_mm10.sizes"
	Int ?ForceFrag
	#RNA
	Map[String,Array[String]] refConfig_r = {
		"mm10" : ["mm10","mouse","/jdfssz1/ST_SUPERCELLS/PUB/scRNA/pipeline/v3.1.5/database/Mouse/mm10_gex_2020A","/jdfssz1/ST_SUPERCELLS/PUB/scRNA/pipeline/v3.1.5/database/Mouse/genes.gtf"],
		"hg38" : ["hg38","human","/jdfssz1/ST_SUPERCELLS/PUB/scRNA/pipeline/v3.1.5/database/Human/Human_2020A","/jdfssz1/ST_SUPERCELLS/PUB/scRNA/pipeline/v3.1.5/database/Human/genes.gtf"]
	}
	String Outdir
	String Species=refConfig_r[reference][1]
	String cDNA_Fastq1
	String cDNA_Fastq2
	String BeadsBarcode
	Int? expectCellNum
	Int? forceCellNum
	String Refdir=refConfig_r[reference][2]
	String Gtf=refConfig_r[reference][3]
	String Fai
	String SampleName
	Boolean? includeIntron = true
    String root="/jdfssz1/ST_SUPERCELLS/PUB/scCAT/pipeline/V1"
	String lib_RNA="/jdfssz1/ST_SUPERCELLS/PUB/scCAT/scripts/lib_RNA.sh"
	String lib_ATAC="/jdfssz1/ST_SUPERCELLS/PUB/scCAT/scripts/lib_ATAC.sh"
	String Soupx


	call makedir {
		input:
		outdir=Outdir
	}

	#RNA
	call cDNAanno{
		input:
		outdir=makedir.dir,
		fastq1=cDNA_Fastq1,
		fastq2=cDNA_Fastq2,
		barcode=BeadsBarcode,
		barcode_counts_raw="cDNA_barcode_counts_raw.txt",
		report="cDNA_sequencing_report.csv",
		refdir=Refdir,
		gtf=Gtf,
        root=root,
		intron=includeIntron
	}


	#ATAC
	call mapping {
		input:
		fastq1=Fastq1_atac,
		fastq2=Fastq2_atac,
		outdir=makedir.dir,
		ref_index=ref_index,
		ref=refdir,
		whitelist=whitelist,
        root=root
	}
	call deconvolution { 
		input: 
		lib=lib_ATAC,
		alnbed=mapping.alnbed,
		outdir=makedir.dir,
		ID=ID,
		species = species,
		ForceFrag=ForceFrag, 
		chrmt=chrmt ,
        root=root
	}

	#RNA
	call getM280UMI{
		input:
		outdir=makedir.dir,
		lib=lib_RNA,
		sampleName=SampleName,
		finalbam=cDNAanno.finalsortbam,
		expectCell=expectCellNum,
		forceCell=forceCellNum,
		rawmatrix=cDNAanno.rawmatrix,
		beads_stat=cDNAanno.beads_stat,
		ATAC_barcodeTran=deconvolution.barcodeTran,
        Rscript=Rscript,
        root=root
	}
	call statSaturation{
		input:
		outdir=makedir.dir,
		anno_decon=getM280UMI.anno_decon,
		sampleName=SampleName,
        python=python,
        root=root
	}

	#ATAC
	call peakcount {
		input:
		lib=lib_ATAC,
		listtxt=deconvolution.listtxt,
		macs2sp= macs2sp,
		promo= promo,
		outdir=makedir.dir,
		chrmt=chrmt,
		FragmentFile=deconvolution.FragmentFile,
		ID=ID,
        root=root,
		Rscript=Rscript
	}
	call bdg2bw {
		input:
		ID=ID,
		outdir=makedir.dir,
		pileup_bdg=peakcount.pileup_bdg,
		lambda_bdg=peakcount.lambda_bdg,
		chrome_size=chrome_size,
        root=root
	}

	#RNA
	call count{
		input:
		outdir=makedir.dir,
		anno_decon=getM280UMI.anno_decon,
        cellid=getM280UMI.cellid,
        root=root
	}
	if (includeIntron==true) {
		call splice_matrix{
			input:
			outdir=makedir.dir,
			anno_decon=getM280UMI.anno_decon,
			countmatrix=count.matrixdir,
            cellid=getM280UMI.cellid,
            root=root
		}
	}

    call overlap_matrix{
        input:
        outdir=makedir.dir,
        rna=count.matrixdir,
        atac=peakcount.Peak,
        Rscript=Rscript,
        FragmentFile=deconvolution.FragmentFile,
        lib_RNA=lib_RNA,
        lib_ATAC=lib_ATAC,
		tss=tss,
		ID=ID,
		sampleName=SampleName,
        species=Species,
        root=root,
		RawMatrix=cDNAanno.rawmatrix,
		barcodeTran=getM280UMI.barcodeTran,
		Soupx=Soupx,
		chrmt=chrmt
    }

    call report{
        input:
		outdir=makedir.dir, 
		root=root,
		fastq1_rna=cDNA_Fastq1,
		fastq2_rna=cDNA_Fastq2, 
		fastq1_atac=Fastq1_atac,
		fastq2_atac=Fastq2_atac,
		ID=SampleName,
		refdir_rna=refdir,
		refdir=ref_index,
		species=species,
		cell_report_2=overlap_matrix.cell_report_2,
		mergestat=getM280UMI.mergestat,
		matrixdir=count.matrixdir,
		DNAsequencing=cDNAanno.DNAsequencing,
		mapping=cDNAanno.mapping,
		annostat=cDNAanno.annostat,
		cluster=overlap_matrix.cluster,
		marker=overlap_matrix.marker,
		QCfig_RNA=overlap_matrix.QCfig_RNA,
		saturation=statSaturation.saturation
    }
	
}


## task.1 Create directory ##
task makedir{
	String outdir
	command<<<
		# mkdir -p ${outdir}
		# mkdir -p ${outdir}/RNA/01.cDNAAnno
		# mkdir -p ${outdir}/RNA/02.M280UMI_stat
		# mkdir -p ${outdir}/RNA/03.Matrix

        # mkdir -p ${outdir}/ATAC
		# mkdir -p ${outdir}/ATAC/02.d2cfile
		# mkdir -p ${outdir}/ATAC/01.out
		# mkdir -p ${outdir}/ATAC/01.out/Peak
		# mkdir -p ${outdir}/ATAC/01.out/Promoter

        # mkdir -p ${outdir}/Joint
        # mkdir -p ${outdir}/Joint/plot
        # mkdir -p ${outdir}/Joint/report
		# mkdir -p ${outdir}/Joint/report/div
        # mkdir -p ${outdir}/Joint/report/base64
        # mkdir -p ${outdir}/Joint/report/table
		# mkdir -p ${outdir}/Joint/report/ATAC
		# mkdir -p ${outdir}/Joint/report/RNA
	>>>
	output{
		String dir="${outdir}"
	}
}

task cDNAanno{
	String fastq1
	String fastq2
	String barcode
	String barcode_counts_raw
	String report
	String outdir
	String refdir
	String gtf
    String root
	Boolean? intron
	command<<<
		# export LD_LIBRARY_PATH=/ldfssz1/ST_BI/USER/software/lib/:$LD_LIBRARY_PATH
		# export PATH=/jdfssz1/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/bin:$PATH
		# export LD_LIBRARY_PATH="/jdfssz1/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/lib:/jdfssz1/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/lib64:$LD_LIBRARY_PATH"

		# echo "${fastq1}" |sed 's/,/\n/g' > ${outdir}/RNA/01.cDNAAnno/in1.list
		# echo "${fastq2}" |sed 's/,/\n/g' > ${outdir}/RNA/01.cDNAAnno/in2.list

		# bcPara=${outdir}/RNA/01.cDNAAnno/bc.Para
		# echo "in1=${outdir}/RNA/01.cDNAAnno/in1.list" > $bcPara
		# echo "in2=${outdir}/RNA/01.cDNAAnno/in2.list" >> $bcPara
		# echo "config=${barcode}" >> $bcPara
		# echo "cbdis=${outdir}/RNA/01.cDNAAnno/${barcode_counts_raw}" >> $bcPara
		# echo "report=${outdir}/RNA/01.cDNAAnno/${report}" >> $bcPara
		# echo "adapter=/jdfssz1/ST_SUPERCELLS/PUB/scRNA/pipeline/v3.1.5/common/file/adapter.txt" >> $bcPara
		
		# mkdir -p ${outdir}/RNA/01.cDNAAnno/RawMatrix

		# /usr/bin/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" ${root}/bin/scStar --outSAMattributes singleCell --outSAMtype BAM Unsorted --genomeDir ${refdir} --outFileNamePrefix ${outdir}/RNA/01.cDNAAnno/ --stParaFile $bcPara --outSAMmode NoQS --runThreadN 10 --limitOutSJcollapsed 10000000 --limitIObufferSize 350000000
		
		# /usr/bin/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" ${root}/bin/Anno -I ${outdir}/RNA/01.cDNAAnno/Aligned.out.bam -a ${gtf} -L ${outdir}/RNA/01.cDNAAnno/${barcode_counts_raw} -o ${outdir}/RNA/01.cDNAAnno -c 10 -m chrM -B ${barcode} ${true="--intron" false="" intron } --anno 1
		# rm -rf ${outdir}/01.cDNAAnno/tmp.*.bam
		# /usr/bin/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" ${root}/bin/samtools sort -@ 10 ${outdir}/RNA/01.cDNAAnno/final.bam -o ${outdir}/RNA/01.cDNAAnno/final.sorted.bam
		# rm -rf ${outdir}/01.cDNAAnno/tmp.*.bam
		# /usr/bin/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" ${root}/bin/samtools index ${outdir}/RNA/01.cDNAAnno/final.sorted.bam
		# /usr/bin/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" ${root}/bin/PISA count -one-hit -@ 10 -cb CB -anno-tag GN -umi UB -outdir ${outdir}/RNA/01.cDNAAnno/RawMatrix ${outdir}/RNA/01.cDNAAnno/final.sorted.bam

		# echo 'cDNAanno Done!'
	>>>
	output{
		String beads_stat="${outdir}/RNA/01.cDNAAnno/beads_stat.txt"
		File mapping="${outdir}/RNA/01.cDNAAnno/alignment_report.csv"
		File annostat="${outdir}/RNA/01.cDNAAnno/anno_report.csv"
		File finalsortbam="${outdir}/RNA/01.cDNAAnno/final.sorted.bam"
		File sequencing="${outdir}/RNA/01.cDNAAnno/${report}"
		File barcodeFile="${outdir}/RNA/01.cDNAAnno/${barcode_counts_raw}"
		File rawmatrix="${outdir}/RNA/01.cDNAAnno/RawMatrix"
		File DNAsequencing="${outdir}/RNA/01.cDNAAnno/cDNA_sequencing_report.csv"
	}
}


task  mapping{
	String fastq1
	String fastq2
	String outdir
	String ref_index
	String ref
	String ?lib_ATAC
	String whitelist
    String root
	Int cpp=4
	Int mem=20
	command <<<
		# export PATH=/jdfssz1/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/bin:$PATH
		# export LD_LIBRARY_PATH="/jdfssz1/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/lib:/jdfssz1/ST_BIGDATA/USER/zhaofuxiang/lib/gcc-9.1.0/lib64:$LD_LIBRARY_PATH"
		# /usr/bin/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" ${root}/bin/chromap --preset atac --bc-error-threshold 0 --trim-adapters -x ${ref_index} -r ${ref} -1 ${sep=',' fastq1} -2 ${sep=',' fastq2} -o ${outdir}/ATAC/01.out/aln.bed --barcode ${sep=',' fastq1} --barcode-whitelist ${whitelist} --read-format bc:0:19,r1:20:-1 -t ${cpp} 2> ${outdir}/ATAC/01.out/alignment_report.tsv
		
		# echo 'mapping Done!'
	>>>
	runtime{
		backend:"Local"
		cpu:cpp
		memory:"${mem} GB"
	}
	output {
		String alnbed="${outdir}/ATAC/01.out/aln.bed"
	}
}


## task.4 d2c ## 
task deconvolution {
	String alnbed
	String outdir
	String ?lib
	String species
	String chrmt
	String ID
    String root
	Int ?ForceFrag
	Int ?mapq
	Int cpp=1
	Int mem=7
	command <<<	
	
		export LD_LIBRARY_PATH="/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/pipeline/scATAC/software/v1.3.7/gcclib/lib:/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/pipeline/scATAC/software/v1.3.7/gcclib/lib64:$LD_LIBRARY_PATH" && export PATH="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/Python-3/bin:$PATH"
		${root}/bin/d2c_v1.4.4/bin/d2c merge -i ${alnbed} --mapq 30 --bf ${default=0 ForceFrag} -o ${outdir}/ATAC/02.d2cfile -c 10 -n ${ID} --mc ${chrmt} -r ${species} --sat --bt1 CB
		cp ${outdir}/ATAC/02.d2cfile/${ID}.CorrelationBarcodes.tsv.gz ${outdir}/Joint/report/ATAC/plot_input2_Jaccard_Overlap_Knee.csv.gz
		cp ${outdir}/ATAC/02.d2cfile/${ID}.barcodeCount.tsv ${outdir}/Joint/report/ATAC/plot_input1_Bead_Barcode_Knee.csv
		tail -n +2 ${outdir}/ATAC/02.d2cfile/${ID}.Metadata.tsv |awk '{print $1}' > ${outdir}/ATAC/02.d2cfile/list.${ID}.txt
		
		echo 'deconvolution Done!'
	>>>
	runtime{
		backend:"Local"
		cpu:cpp
		memory:"${mem} GB"
	}
	output {
	String FragmentFile="${outdir}/ATAC/02.d2cfile/${ID}.fragments.tsv.gz"
	String listtxt= "${outdir}/ATAC/02.d2cfile/list.${ID}.txt"
	File FragmentFile1="${outdir}/ATAC/02.d2cfile/${ID}.fragments.tsv.gz"
	File listtxt1= "${outdir}/ATAC/02.d2cfile/list.${ID}.txt"
	File barcodeTran="${outdir}/ATAC/02.d2cfile/${ID}.barcodeMerge.tsv"
	}
}

task getM280UMI{
	String outdir
	String root
	String sampleName
	String finalbam
	String beads_stat
	String lib
	Int? expectCell
	Int? forceCell
	String rawmatrix
	String ATAC_barcodeTran
    String Rscript
	command<<<
		# if [ -f ${default=abjdbashj lib} ]; then
		# source ${lib}
		# fi
		
		# perl ${root}/scripts/to16.pl ${ATAC_barcodeTran} > ${outdir}/RNA/02.M280UMI_stat/${sampleName}_barcodeTranslate_16.txt
        # awk '{print $2}' ${outdir}/RNA/02.M280UMI_stat/${sampleName}_barcodeTranslate_16.txt > ${outdir}/RNA/02.M280UMI_stat/cellid.txt
		# cp ${ATAC_barcodeTran} ${outdir}/RNA/02.M280UMI_stat/${sampleName}_barcodeTranslate.txt
		# /usr/bin/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" ${root}/bin/tagAdd -n 4  -bam ${finalbam} -file ${outdir}/RNA/02.M280UMI_stat/${sampleName}_barcodeTranslate_16.txt -out ${outdir}/RNA/02.M280UMI_stat/anno_decon.bam -tag_check CB:Z: -tag_add DB:Z: &&
		# perl ${root}/scripts/cell_stat.pl -c ${beads_stat} -m ${ATAC_barcodeTran} -o ${outdir}/RNA/02.M280UMI_stat/ &&
		# ${Rscript} ${root}/scripts/CellMergeStat.R -I ${outdir}/RNA/02.M280UMI_stat/${sampleName}_barcodeTranslate_16.txt -O ${outdir}/RNA/02.M280UMI_stat -n ${sampleName}

		# echo 'getM280UMI Done!'
	>>>
	output{
		File anno_decon="${outdir}/RNA/02.M280UMI_stat/anno_decon.bam"
		File mergestat="${outdir}/RNA/02.M280UMI_stat/merge_cell.stat"
		File barcodeTran="${outdir}/RNA/02.M280UMI_stat/${sampleName}_barcodeTranslate_16.txt"
        File cellid="${outdir}/RNA/02.M280UMI_stat/cellid.txt"
	}

}
task statSaturation{
	String outdir
	String sampleName
	String anno_decon
    String python
    String root
	command<<<
		# /usr/bin/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" ${root}/bin/samtools index ${anno_decon}
		# /usr/bin/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" ${python} ${root}/scripts/saturation_sp.py -i ${anno_decon} -o ${outdir}/RNA/02.M280UMI_stat --quality 20 --threads 10 -n ${sampleName}

		# echo 'statSaturation Done!'
	>>>
	output{
		File saturation = "${outdir}/RNA/02.M280UMI_stat/${sampleName}_Saturation.tsv"
	}

}

## task.5 get Raw matrix ##
task peakcount {
	String outdir
	String ?lib
	String macs2sp
	String promo
	String ID
	String FragmentFile
	String listtxt
	String chrmt
	Int cpp=1
	Int mem=14
    String root
	String Rscript
	command <<<
		# if [ -f ${default=abjdbashj lib} ]; then
		# source ${lib}
		# fi

		# /usr/bin/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" ${root}/bin/macs2 callpeak -t ${FragmentFile} -f BED -g ${macs2sp} -n ${ID} -B -q 0.001 --nomodel --outdir ${outdir}/ATAC/01.out
		# /usr/bin/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)"  ${Rscript} ${root}/scripts/C4scATAC_Cluster_AnnotationV3.R -I ${outdir}/ATAC/01.out/${ID}_peaks.narrowPeak -F ${FragmentFile} -C ${listtxt} -G ${promo} -MT ${chrmt} -Q ${outdir}/ATAC/02.d2cfile/${ID}.Metadata.tsv -O ${outdir}/ATAC

		# echo 'peakcount Done!'
	>>>
	runtime{
		backend:"Local"
		cpu:cpp
		memory:"${mem} GB"
	}
	output {
		String Qcfie2 = "${outdir}/ATAC/01.out/${ID}_peaks.narrowPeak"
		File Qcfie21 = "${outdir}/ATAC/01.out/${ID}_peaks.narrowPeak"
		File pileup_bdg = "${outdir}/ATAC/01.out/${ID}_treat_pileup.bdg"
		File lambda_bdg = "${outdir}/ATAC/01.out/${ID}_control_lambda.bdg"
        String Peak = "${outdir}/ATAC/01.out/Peak"
	}
}

task bdg2bw {
    String ID
    String outdir
    String pileup_bdg
    String lambda_bdg
    String chrome_size
    String root
    command <<<
        # sed -i 1d ${pileup_bdg}
        # sed -i 1d ${lambda_bdg}
        # ${root}/bin/macs2 bdgcmp -t ${pileup_bdg} -c ${lambda_bdg} -o ${outdir}/ATAC/01.out/${ID}_FE.bdg -m FE

        # ${root}/bin/bedClip ${outdir}/ATAC/01.out/${ID}_FE.bdg ${chrome_size} ${outdir}/ATAC/01.out/${ID}_FE.clip
        # LC_COLLATE=C sort -k1,1 -k2,2n ${outdir}/ATAC/01.out/${ID}_FE.clip  > ${outdir}/ATAC/01.out/${ID}_FE.sort.clip

        # ${root}/bin/bedGraphToBigWig ${outdir}/ATAC/01.out/${ID}_FE.sort.clip ${chrome_size} ${outdir}/ATAC/01.out/${ID}_FE.bw

		# echo 'bdg2bw Done!'
    >>>
}


#RNA
task count{
	String outdir
	String anno_decon
    String cellid
    String root
	command<<<
		# mkdir -p ${outdir}/RNA/03.Matrix/FilterMatrix
		# /usr/bin/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" ${root}/bin/PISA count -one-hit -@ 10 -cb DB -anno-tag GN -umi UB -outdir ${outdir}/RNA/03.Matrix/FilterMatrix ${anno_decon} -list ${cellid}
		
		# echo 'count Done!'
	>>>
	output{
		String matrixdir = "${outdir}/RNA/03.Matrix/FilterMatrix"
	}
}
task splice_matrix{
	String outdir
	String anno_decon
	String countmatrix
    String cellid
    String root
	command<<<
		# mkdir -p ${outdir}/RNA/03.Matrix/SpliceMatrix
		# mkdir -p ${outdir}/RNA/03.Matrix/RNAVelocityMatrix
		# /usr/bin/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" ${root}/bin/PISA count -one-hit -@ 10 -cb DB -ttype E,S -anno-tag GN -umi UB -outdir ${outdir}/RNA/03.Matrix/SpliceMatrix ${anno_decon} -list ${cellid}
		# /usr/bin/time -f "%e(elapsed) %U(system) %P(cpu) %M(max RAM KB)" ${root}/bin/PISA count -one-hit -velo -@ 10 -cb DB -anno-tag GN -umi UB -outdir ${outdir}/RNA/03.Matrix/RNAVelocityMatrix ${anno_decon} -list ${cellid}
		
		# echo 'splice matrix Done!'
	>>>
	output{
		String matrixdir = "${outdir}/RNA/03.Matrix/RNAVelocityMatrix"
	}
}

task overlap_matrix{
    String rna
    String atac
    String Rscript
    String outdir
    String FragmentFile
    String ?dim
    Float ?doublepercentage
    String ?mitpercentage
    String ?minfeatures
	String ?PCusage
    Float ?resolution=0.5
    String sampleName
    String lib_RNA
    String lib_ATAC
    String tss
    String species
    String ID
    String root
	String RawMatrix
	String barcodeTran
	String Soupx
	String chrmt
    command<<<
        # ${Rscript} ${root}/scripts/wnn.R -R ${rna} -A ${atac} -SN ${sampleName} -F ${FragmentFile} -SP ${species} -RM ${RawMatrix} -BT ${barcodeTran} -D ${default=20 dim} -MP ${default=10 mitpercentage} -MF ${default=200 minfeatures} -PC ${default=50 PCusage} -RES ${resolution} -O ${outdir}/Joint -MT ${chrmt} -IS ${Soupx}

		
        # # ATAC
        # ${Rscript} ${root}/scripts/plot_TSSEnrichment_FragSize.R -T ${tss} -F ${FragmentFile} -G ${ID} -O ${outdir}/Joint/report/ATAC

		# nf=$(gunzip -c ${outdir}/ATAC/02.d2cfile/${ID}.fragments.tsv.gz | awk '{if($3-$2<147) print $0}' | wc -l)
		# mn=$(gunzip -c ${outdir}/ATAC/02.d2cfile/${ID}.fragments.tsv.gz | awk '{if($3-$2>147 && $3-$2<294) print $0}' | wc -l)
		# total=$(gunzip -c ${outdir}/ATAC/02.d2cfile/${ID}.fragments.tsv.gz | wc -l)
		# echo "$nf $total" | awk '{print "Fraction of nucleosome-free regions:"$1/$2*100"%"}' >> ${outdir}/Joint/report/ATAC/5_2.library.QC.csv
		# echo "$mn $total" | awk '{print "Fraction of fragments mono-nucleosome regions:"$1/$2*100"%"}' >> ${outdir}/Joint/report/ATAC/5_2.library.QC.csv

		# echo 'overlap matrix Done!'
    >>>
    output{
		File QCfig_RNA="${outdir}/Joint/plot/${sampleName}_raw_QCplot.svg"
        String Qcfie_ATAC = "${outdir}/Joint/report/ATAC/5_2.library.QC.csv"
		File Qcfie1_ATAC = "${outdir}/Joint/report/ATAC/5_2.library.QC.csv"
		File clusterRDS="${outdir}/Joint/WNN_annotation.RDS"
		File cell_report_2="${outdir}/Joint/cell_report_2.csv"
		File cluster="${outdir}/Joint/cluster_rna.csv"
		File marker="${outdir}/Joint/${sampleName}_rna_marker.csv"
    }
}

task report{
    String outdir
	String root
	String fastq1_rna
	String fastq2_rna
	String fastq1_atac
	String fastq2_atac
	String ID
	String refdir_rna
	String refdir
	String species
	String cell_report_2
	String mergestat
	String matrixdir
	String DNAsequencing
	String mapping
	String annostat
	String cluster
	String marker
	String QCfig_RNA
	String saturation
	
	command <<<
        # atac
		echo "Sample ID:${ID}" > ${outdir}/Joint/report/ATAC/2.sample.csv
		# echo "FASTQ path:${sep=',' fastq1_atac},${sep=',' fastq2_atac}" >>  ${outdir}/Joint/report/ATAC/2.sample.csv
		echo "FASTQ path:${fastq1_atac},${fastq2_atac}" >>  ${outdir}/Joint/report/ATAC/2.sample.csv

		echo "Pipeline version:v1.0.0" >> ${outdir}/Joint/report/ATAC/2.sample.csv
		echo "Reference path:${refdir}" >> ${outdir}/Joint/report/ATAC/2.sample.csv

		grep "Number of reads:" ${outdir}/ATAC/01.out/alignment_report.tsv | sed 's/\.//g' | awk -F ":" '{printf "Total number of reads :%'"'"'0.0f\n",$2/2}' > ${outdir}/Joint/report/ATAC/3.sequencing.csv
		grep "Number of barcodes in whitelist\|Number of reads:" ${outdir}/ATAC/01.out/alignment_report.tsv |sed 's/\.//g'|awk -F ':' '{print $2}'|awk 'NR==1{tmp=$1}NR>1{printf "Reads pairs with a valid barcode:%'"'"'0.2f%\n",$1*2/tmp*100}'  >> ${outdir}/Joint/report/ATAC/3.sequencing.csv

		grep "Number of mapped reads:\|Number of barcodes in whitelist" ${outdir}/ATAC/01.out/alignment_report.tsv |sed 's/\.//g'|awk -F ':' '{print $2}' |awk 'NR==1{tmp=$1}NR>1{printf "Reads pairs with a valid barcode:%'"'"'0.2f%\n",tmp/($1*2)*100}' >> ${outdir}/Joint/report/ATAC/3.mapping.csv
		grep "Number of uniquely mapped reads:" ${outdir}/ATAC/01.out/alignment_report.tsv | sed 's/\.//g'| awk -F ":" '{printf "Uniquely mapped reads:%'"'"'0.0f\n",$2}' >> ${outdir}/Joint/report/ATAC/3.mapping.csv
		# echo "Mitochondria ratio:-" >> ${outdir}/Joint/report/ATAC/3.mapping.csv

		grep "bead_cutoff" ${outdir}/ATAC/02.d2cfile/${ID}.d2cCutoff.tsv | awk '{printf "Bead threshold:%'"'"'0.0f\n",$2}' > ${outdir}/Joint/report/ATAC/4.cells.csv
		wc -l ${outdir}/ATAC/02.d2cfile/${ID}.barcodeMerge.tsv > ${outdir}/Joint/report/ATAC/wc_barcode.txt
		awk -F "," '{printf "Bead number:%'"'"'0.0f\n",$1}' ${outdir}/Joint/report/ATAC/wc_barcode.txt >> ${outdir}/Joint/report/ATAC/4.cells.csv
		grep "cor_cutoff" ${outdir}/ATAC/02.d2cfile/${ID}.d2cCutoff.tsv | awk '{printf "cor threshold:%0.5f\n",$2}' >> ${outdir}/Joint/report/ATAC/4.cells.csv
		rm ${outdir}/Joint/report/ATAC/wc_barcode.txt



		${root}/bin/python /jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/pipeline/scCAT/pipeline/V1.1/html/barcode1.py --outPath ${outdir}/Joint
		${root}/bin/python /jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/pipeline/scCAT/pipeline/V1.1/html/jaccard.py --outPath ${outdir}/Joint
		${root}/bin/python /jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/pipeline/scCAT/pipeline/V1.1/html/svg_to_base64_string.py --outPath ${outdir}
		${root}/bin/python /jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/pipeline/scCAT/pipeline/V1.1/html/st.py --outPath ${outdir} --ID ${ID}




		# rna
		 "SampleName,${ID}" >${outdir}/Joint/report/RNA/1.cell_report.csv
		echo "Species,${species}" >>${outdir}/Joint/report/RNA/1.cell_report.csv
		`${root}/bin/Rscript /jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/pipeline/scCAT/pipeline/V1.1/html/cell_report_1.R -I ${mergestat} -M ${matrixdir} -O ${outdir}/Joint/report/RNA/`
		`cat ${outdir}/Joint/report/RNA/cell_report_1.csv ${cell_report_2} >>${outdir}/Joint/report/RNA/1.cell_report.csv`
		`rm ${outdir}/Joint/report/RNA/cell_report_1.csv`

		echo "FASTQ path:${fastq1_rna},${fastq2_rna}" >  ${outdir}/Joint/report/RNA/2.sample.csv
		echo "Reference path:${refdir_rna}" >> ${outdir}/Joint/report/RNA/2.sample.csv

		/jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/mawen/src/tools/Rscript /jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/pipeline/scCAT/pipeline/V1.1/html/get_new_cutoff.R --matrix ${outdir}/RNA/01.cDNAAnno/RawMatrix --tran ${outdir}/ATAC/02.d2cfile/${ID}.barcodeMerge.tsv --hex /hwfssz5/ST_SUPERCELLS/P21Z10200N0090/songyumo/pipelines/scRNA/version/V3.1.5/common/file/hex_barcode_map.tsv --output ${outdir}/Joint/report/RNA

		perl ${root}/scripts/count2per_for_report.pl ${DNAsequencing} > ${outdir}/Joint/report/RNA/3.cDNA_sequencing_report.csv
		num=`cat ${DNAsequencing} |sed -n 2p | awk -F ',' '{print $2}'` && map=`cat ${mapping} |sed -n 1p|awk -F ',' '{print $2}' ` && percent_2=`awk 'BEGIN{printf "%.2f%%\n",('$map'/'$num')*100}'` && printf "Raw reads,$num\nMapped reads,$percent_2%\n" > ${outdir}/Joint/report/RNA/4.alignment_reporttmp.csv
		`cat ${mapping} |sed -n "3,6p"  >> ${outdir}/Joint/report/RNA/4.alignment_reporttmp.csv`
		perl ${root}/scripts/count2per_for_report.pl ${outdir}/Joint/report/RNA/4.alignment_reporttmp.csv > ${outdir}/Joint/report/RNA/4.alignment_report.csv
		rm -rf ${outdir}/Joint/report/RNA/4.alignment_reporttmp.csv

		rna=$(tail -n 1 ${outdir}/Joint/report/RNA/11.saturation.tsv |awk '{print $NF}') && atac=$(tail -n 1 ${outdir}/ATAC/02.d2cfile/${ID}.sequenceSaturation.tsv |awk '{print $3}') && rnap=$(printf "%.2f" `echo "scale=2;($rna/100)"|bc`) && atacp=$(printf "%.2f" `echo "scale=2;($atac)"|bc`) && printf "RNA saturation ratio:$rnap\nATAC saturation ratio:$atacp\n" > ${outdir}/Joint/report/saturation.csv

		`cp ${annostat} ${outdir}/Joint/report/RNA/5.anno_report.csv`
		`cp ${cluster} ${outdir}/Joint/report/RNA/9.cluster.csv`
		`cp ${marker} ${outdir}/Joint/report/RNA/10.marker.csv`
		`cp ${QCfig_RNA} ${outdir}/Joint/report/RNA/7.${ID}_raw_QCplot.svg`
		`cp ${saturation} ${outdir}/Joint/report/RNA/11.saturation.tsv`
		cp ${outdir}/RNA/02.M280UMI_stat/${ID}_CellNumber_merge.svg ${outdir}/Joint/report/RNA/6.${ID}_CellNumber_merge.svg



		# all
		${root}/bin/python /jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/pipeline/scCAT/pipeline/V1.1/html/pre_process.py --outPath ${outdir}/Joint

		${root}/bin/python /jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/pipeline/scCAT/pipeline/V1.1/html/generate.py --outPath ${outdir}/Joint --ID ${ID} --intron true --htmlTemplate /jdfssz1/ST_SUPERCELLS/P21Z10200N0090/Automated/USER/songyumo/pipeline/scCAT/pipeline/V1.1/html/template-joint-atac.html 

		rm -rf ${outdir}/Joint/report/saturation.csv
	>>>
}


