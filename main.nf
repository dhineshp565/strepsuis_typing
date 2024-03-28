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
    output:
    val(SampleName),emit:sample
	tuple val(SampleName),path("${SampleName}_flye.fasta"),emit:assembly
	path("${SampleName}_flye-info.txt"),emit:flyeinfo
    script:
    """
    dragonflye --reads ${SamplePath} --outdir ${SampleName}_assembly --model r941_e81_hac_g514 --gsize 2.4M --nanohq --medaka 1
    # rename fasta file with samplename
    mv "${SampleName}_assembly"/flye.fasta "${SampleName}"_flye.fasta
    # rename fasta header with samplename
    sed -i 's/contig/${SampleName}_contig/g' "${SampleName}_flye.fasta"
     # rename flyeinfo file and contnents
    mv "${SampleName}_assembly"/flye-info.txt "${SampleName}"_flye-info.txt
    sed -i 's/contig/${SampleName}_contig/g' "${SampleName}_flye-info.txt"
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
	tuple val(SampleName),path(consensus)
	output:
	path("${SampleName}_MLST.csv")
	script:
	"""
	mlst --legacy --scheme ssuis ${consensus} > ${SampleName}_MLST.csv
	sed -i 's,_flye.fasta,,g' ${SampleName}_MLST.csv
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
	abricate -datadir ${vcf_db} --db Ssuis_vfdb ${consensus} 1> ${SampleName}_vf.csv
	sed -i 's,_flye.fasta,,g' ${SampleName}_vf.csv
	abricate --db card ${consensus} 1> ${SampleName}_AMR.csv
	sed -i 's,_flye.fasta,,g' ${SampleName}_AMR.csv
	
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
	path ("*_serotype.csv")
	script:
	"""
	
	sero.sh ${csvfile}
	"""
}

    
process resolve_sero {
		publishDir "${params.out_dir}/res_serotype/",mode:"copy"
		label "medium"
        input:
        tuple val(SampleName),path(assembly)
		path(reference)
		tuple val(SampleName),path(serotypefile)
		output:
		path ("${SampleName}_serotype.csv"),emit:sero
		path ("${SampleName}_vcf.csv"),emit:vcf
		script:
        """
		if grep -wq "cps-1" ${serotypefile} || grep -wq "cps-2" ${serotypefile};then

			minimap2 -ax map-ont ${reference} ${assembly} > ${SampleName}.sam
			samtools view -b -S ${SampleName}.sam > ${SampleName}.bam
			samtools sort ${SampleName}.bam > ${SampleName}_sorted.bam
			samtools faidx ${reference} 
			samtools index ${SampleName}_sorted.bam > ${SampleName}_sorted.bai
			bcftools mpileup -Ob -f ${reference} ${SampleName}_sorted.bam > ${SampleName}.bcf
			bcftools call -mv -Ob ${SampleName}.bcf > ${SampleName}.vcf
			bcftools view ${SampleName}.vcf > ${SampleName}_vcf.csv
			if grep -q -e "cps-2" ${serotypefile} && grep -q -e "483" ${SampleName}_vcf.csv ;then
				sed -i 's,serotype-2,serotype-1/2,g' ${serotypefile}
				sed -i 's/_flye.fasta//g' ${serotypefile}
				cut -f 1,2,6,10,11,15 ${serotypefile} > ${SampleName}_serotype.csv
			fi
			if grep -q -e "cps-1" ${serotypefile} && grep -q -e "483" ${SampleName}_vcf.csv ;then
				sed -i 's,serotype-1,serotype-14,g' ${serotypefile}
				sed -i 's,_flye.fasta,,g' ${serotypefile}
				cut -f 1,2,6,10,11,15 ${serotypefile} > ${SampleName}_serotype.csv
			fi
		else
				sed -i 's,_flye.fasta,,g' ${serotypefile}
				cut -f 1,2,6,10,11,15 ${serotypefile} > ${SampleName}_serotype.csv
				echo "no vcf" > ${SampleName}_vcf.csv
		fi


        """

}

process make_limsfile {
	label "low"
	publishDir "${params.out_dir}/LIMS",mode:"copy"
	input:
	path (serotyping_results)
	path (mlst_results)
	output:
	path("Ssuis_sero_file.csv")
	path("Ssuis_MLST_file.csv")
	script:
	"""
	awk 'FNR==1 && NR!=1 { while (/^#F/) getline; } 1 {print}' ${serotyping_results} > Ssuis_sero_file.csv
	sed -i 's/#FILE/id/g' Ssuis_sero_file.csv
	awk 'FNR==1 && NR!=1 { while (/^FILE/) getline; } 1 {print}' ${mlst_results} > Ssuis_MLST_file.csv
	sed -i 's/FILE/id/g' Ssuis_MLST_file.csv
	"""
}

process make_report {
	label "low"
	publishDir "${params.out_dir}/reports",mode:"copy"
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
	cp ${serofile} serofile.csv
	cp ${mlstfile} mlstfile.csv

	Rscript -e 'rmarkdown::render(input="rmdfile_copy.Rmd",params=list(sero="serofile.csv",mlst="mlstfile.csv",csv="samples.csv"),output_file="Ssuis_report.html")'
	"""

}





workflow {
    data=Channel
	.fromPath(params.input)
	merge_fastq(make_csv(data).splitCsv(header:true).map { row-> tuple(row.SampleName,row.SamplePath)})
	// Merge fastq files for each sample

	// based on the optional argument trim barcodes using porechop and assemble using dragonflye
    if (params.trim_barcodes){
		porechop(merge_fastq.out)
		dragonflye(porechop.out) 
	} else {
        dragonflye(merge_fastq.out)           
    }
	versionfile=file("${baseDir}/software_version.csv")
    busco(dragonflye.out.assembly)
    mlst (dragonflye.out.assembly)
	sero_db=("${baseDir}/Ssuis_serotype_db")
	vcf_db=("${baseDir}/Ssuis_VF_db")
    abricate (dragonflye.out.assembly,sero_db,vcf_db)
	reference=file("${baseDir}/Ssuis_cps2K.fasta")
	minimap2 (dragonflye.out.assembly,reference)
	samtools(minimap2.out,reference)
	bcftools(samtools.out)
	serotype(abricate.out.sero.collect(),bcftools.out.collect(),make_csv.out)
    //resolve_sero(dragonflye.out.assembly,reference,abricate.out)
    make_limsfile (serotype.out.collect(),mlst.out.collect())
	

	rmd_file=file("${baseDir}/Ssuis_report.Rmd")
	make_report (rmd_file,make_limsfile.out,busco.out.collect(),make_csv.out,abricate.out.vif.collect(),abricate.out.AMR.collect())


}

