#!/usr/bin/env nextflow
nextflow.enable.dsl=2



// make csv file with headers from the given input

process make_csv {
	publishDir "${params.out_dir}"
	input:
	path(fastq_input)
	output:
	path("samplelist.csv")
	
	script:
	"""
	makecsv.sh ${fastq_input}

	"""

}

//merge fastq files for each SampleName and create a merged file for each SampleNames
process merge_fastq {
	publishDir "${params.out_dir}/merged"
	label "low"
	input:
	tuple val(SampleName),path(SamplePath)
	output:
	tuple val(SampleName),path("${SampleName}.{fastq,fastq.gz}"),emit:reads

	shell:
	"""
	count=\$(ls -1 ${SamplePath}/*.gz 2>/dev/null | wc -l)
	
	
		if [[ "\${count}" != "0" ]]
		then
			cat ${SamplePath}/*.fastq.gz > ${SampleName}.fastq.gz
					
		else
			count=\$(ls -1 ${SamplePath}/*.fastq 2>/dev/null | wc -l)
			if [[ "\${count}" != "0" ]]
			then
				cat ${SamplePath}/*.fastq > ${SampleName}.fastq
				
			fi
		fi
	"""
}

//trim barcodes and adapter using porechop

process porechop {
	label "high"
	publishDir "${params.out_dir}/trimmed"
	input:
	tuple val(SampleName),path(SamplePath)
	output:
	tuple val(SampleName),path ("${SampleName}_trimmed.fastq")
	script:
	"""
	porechop -i ${SamplePath} -o ${SampleName}_trimmed.fastq
	"""
}

process dragonflye {
    label "high"
    publishDir "${params.out_dir}/Assembly",mode:"copy"
    input:
    tuple val(SampleName),path(SamplePath)
	val(medaka_model)
    output:
	tuple val(SampleName),path("${SampleName}_flye.fasta")
	path("${SampleName}_flye-info.txt"),emit:flyeinfo
    script:
    """
    dragonflye --reads ${SamplePath} --outdir ${SampleName}_assembly --gsize 2.4M --nanohq 
    # rename fasta file with samplename
    mv "${SampleName}_assembly"/flye.fasta "${SampleName}"_flye.fasta
    # rename fasta header with samplename
    sed -i 's/contig/${SampleName}_contig/g' "${SampleName}_flye.fasta"
     # rename flyeinfo file and contnents
    mv "${SampleName}_assembly"/flye-info.txt "${SampleName}"_flye-info.txt
    sed -i 's/contig/${SampleName}_contig/g' "${SampleName}_flye-info.txt"
    """
}


process medaka {
	publishDir "${params.out_dir}/medaka",mode:"copy"
	label "high"
	input:
	tuple val(SampleName),path(SamplePath)
	tuple val(SampleName),path(draft_assembly)
	path("${SampleName}_flye-info.txt")
	val (medaka_model)
	output:
	tuple val(SampleName),emit:sample
	tuple val(SampleName),path("${SampleName}_assembly.fasta"),emit:assembly
	path("${SampleName}_flye-info.txt"),emit:flyeinfo
	
	script:
	"""
	
	medaka_consensus -i ${SamplePath} -d ${draft_assembly} -o ${SampleName}_medaka_assembly --bacteria
		
	mv ${SampleName}_medaka_assembly/consensus.fasta ${SampleName}_assembly.fasta
	
	"""

}


process busco {
    label "low"
    publishDir "${params.out_dir}/busco",mode:"copy"
    input:
    tuple val(SampleName),path(cons)
    output:
    path ("${SampleName}_busco.txt")
    script:

    """
    busco -i ${cons} -m genome -l bacteria_odb10 -o ${SampleName}_busco_results
	mv ${SampleName}_busco_results/*.txt ${SampleName}_busco.txt

    """
}

process mlst {
	publishDir "${params.out_dir}/mlst/",mode:"copy"
	label "low"
	input:
	tuple val(SampleName),path(assembly)
	output:
	path("${SampleName}_MLST.csv")
	script:
	"""
	mlst ${assembly} > ${SampleName}_MLST.csv
	sed -i 's,_assembly.fasta,,g' ${SampleName}_MLST.csv
	"""
}

process abricate{
	publishDir "${params.out_dir}/abricate/",mode:"copy"
	label "low"
	input:
	tuple val(SampleName),path(consensus)
	path(sero_db)
	path(vcf_db)
	output:
	path("${SampleName}_sero.csv"),emit:sero
	path("${SampleName}_vf.csv"),emit:vif
	path("${SampleName}_AMR.csv"),emit:AMR
	
	script:
	"""
	abricate --datadir ${sero_db} --db Ssuis_serotype -minid 60  -mincov 60 --quiet ${consensus} 1> ${SampleName}_sero.csv
	sed -i 's,_assembly.fasta,,g' ${SampleName}_sero.csv
	abricate -datadir ${vcf_db} --db Ssuis_vfdb ${consensus} 1> ${SampleName}_vf.csv
	sed -i 's,_assembly.fasta,,g' ${SampleName}_vf.csv
	abricate --db card ${consensus} 1> ${SampleName}_AMR.csv
	sed -i 's,_assembly.fasta,,g' ${SampleName}_AMR.csv
	cat *vf.csv > sample_vf.csv
	
	"""
	
}


process minimap2 {
        publishDir "${params.out_dir}/minimap2/",mode:"copy"
		label "low"
        input:
        tuple val(SampleName),path(SamplePath)
		path(reference)
		output:
        tuple val(SampleName),path ("${SampleName}.sam")
        script:
        """
        minimap2 -ax map-ont ${reference} ${SamplePath} > ${SampleName}.sam
        """
}
process samtools {
		 publishDir "${params.out_dir}/samtools/",mode:"copy"
		label "low"
        input:
        tuple val(SampleName),path(samfile)
		path (reference)
        output:
		tuple val(SampleName),path("${SampleName}_sorted.bam")
		path ("${SampleName}_sorted.bai")
		path (reference)
        script:
        """
        samtools view -b -S ${samfile} > ${SampleName}.bam
		samtools sort ${SampleName}.bam > ${SampleName}_sorted.bam
		samtools faidx ${reference} 
		samtools index ${SampleName}_sorted.bam > ${SampleName}_sorted.bai
		
		"""
}
process bcftools {
	 publishDir "${params.out_dir}/vcf/",mode:"copy"
		label "low"
        input:
		tuple val(SampleName),path(sorted_bam)
		path (sorted_bai)
		path(reference)
		
        output:
       	path("${SampleName}_vcf.csv")
		script:
        """
		bcftools mpileup -Ob -f ${reference} ${sorted_bam} > ${SampleName}.bcf
		bcftools call -mv -Ob "${SampleName}.bcf" > ${SampleName}.vcf
		bcftools view "${SampleName}.vcf" > ${SampleName}_vcf.csv
		"""
}

process serotype {
	publishDir "${params.out_dir}/serotype/",mode:"copy"
	label "low"
	input:
	path(abricatefile)
	path(vcf_file)
	path(csvfile)
	output:
	path ("*serotype.csv")
	script:
	"""
	
	sero.sh ${csvfile}
	"""
}


process make_limsfile {
	label "low"
	publishDir "${params.out_dir}/LIMS",mode:"copy"
	input:
	path (serotyping_results)
	path (vf_results)
	path (amr_results)
	path (mlst_results)
	path (software_version)
	output:
	path("*_LIMS_file.csv")
	path("sero_file.csv"),emit:sero
	path("MLST_file_*.csv"),emit:mlst
	
	script:
	"""
	LIMS_file.sh
	
	date=\$(date '+%Y-%m-%d_%H-%M-%S')
	awk 'FNR==1 && NR!=1 { while (/^#F/) getline; } 1 {print}' ${mlst_results} > MLST_file.csv
	# add header to mlst file
	sed -i \$'1 i\\\nSAMPLE\tSCHEME\tST\taroA\tcpn60\tdpr\tgki\tmutS\trecA\tthrA' MLST_file_\${date}.csv
	

	"""
}

process make_report {
	label "low"
	publishDir "${params.out_dir}/",mode:"copy"
	input:
	path(rmdfile)
	path(serofile)
	path (mlstfile)
	path(busco)
	path (samplelist)
	path (vffiles)
	path (amrfiles)
	output:
	path("Ssuis_report.html")

	script:

	"""
	
	cp ${rmdfile} rmdfile_copy.Rmd
	cp ${samplelist} samples.csv
	cp ${mlstfile} mlstfile.csv
	cp ${serofile} serofile.csv

	Rscript -e 'rmarkdown::render(input="rmdfile_copy.Rmd",params=list(sero="serofile.csv",mlst="mlstfile.csv",csv="samples.csv"),output_file="Ssuis_report.html")'
	"""

}







workflow {
    data=Channel
	.fromPath(params.input)
	merge_fastq(make_csv(data).splitCsv(header:true).map { row-> tuple(row.SampleName,row.SamplePath)})
	
	// Merge fastq files for each sample

	// based on the optional argument trim barcodes using porechop, assemble using dragonflye and polish using medaka
    if (params.trim_barcodes){
		porechop(merge_fastq.out)
		dragonflye(porechop.out,params.medaka_model) 
		medaka(porechop.out,dragonflye.out,params.medaka_model)
	} else {
        dragonflye(merge_fastq.out,params.medaka_model) 
		medaka(merge_fastq.out,dragonflye.out,params.medaka_model)          
    }
	versionfile=file("${baseDir}/software_version.csv")
	//checking completeness of assembly
    busco(medaka.out.assembly)
	//mlst
    mlst (medaka.out.assembly)
	//abricate AMR,serotyping and virulence factors
	sero_db=("${baseDir}/Ssuis_serotype_db")
	vcf_db=("${baseDir}/Ssuis_VF_db")
    abricate (medaka.out.assembly,sero_db,vcf_db)
	reference=file("${baseDir}/Ssuis_cps2K.fasta")
	//variant calling to resolve serotype
	minimap2 (medaka.out.assembly,reference)
	samtools(minimap2.out,reference)
	bcftools(samtools.out)
	serotype(abricate.out.sero.collect(),bcftools.out.collect(),make_csv.out)
    
	versionfile=file("${baseDir}/software_version.csv")
	//make lims file
    make_limsfile (serotype.out.collect(),abricate.out.vif.collect(),abricate.out.AMR.collect(),mlst.out.collect(),versionfile)
	
	//report generation

	rmd_file=file("${baseDir}/Ssuis_report.Rmd")
	make_report (rmd_file,make_limsfile.out.sero,make_limsfile.out.mlst,busco.out.collect(),make_csv.out,abricate.out.vif.collect(),abricate.out.AMR.collect())


}

