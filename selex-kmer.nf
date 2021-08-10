#!/usr/bin/env nextflow
"""
========================================================
Groovy Helper Functions
========================================================    
"""
def filter_selex_rounds(String round_name) {
    if (params.round_order == null || params.round_order == "" || params.round_order.size() == 0) return true;
    // // building regex to check files if they match the rounds specified in YOUR_SELEX
    // round_regex = "(" + params.round_order.join('|') + ")" + params.trim_delimiter + ".*"; // TODO replace trim_delimiter with autodetect
    // return s.matches(round_regex);
    
    return params.round_order.contains(round_name);
}

def get_round_id(String round_name) {
    if (params.round_order == null || params.round_order == "" || params.round_order.size() == 0) return 0;
    else return params.round_order.indexOf(round_name);
}

"""
========================================================
Make output directory if it doesn't exist
========================================================    
"""
dir_output = file(params.output_dir)
if (!dir_output.exists()) {
    if (!dir_output.mkdir()) {
        println("Couldn't create output directory.")
    }
}


"""
========================================================
Read FASTA-files from specified input directory and order by round_order (if present)
========================================================    
"""
fasta_files = Channel.fromPath(params.input_dir + "/" + params.fasta_pattern, checkIfExists:true, type: "file")
fasta_files
    .map { it -> [get_round_id(it.getSimpleName()), it.getSimpleName(), it].flatten() }
    .view()
    .filter { filter_selex_rounds(it[1]) }
    .set{ fasta_files_filtered }
    
Channel.fromPath(params.derep_fasta)
	.into { derep_fasta_blast_db; derep_fasta_blastn_query; derep_fasta_pcr_removal; derep_fasta_scoring }
   

"""
========================================================
Removal of PCR Duplicates (Mutations, Indels)
========================================================
"""


process pcr_dup_removal_blast_db {
//	conda 'bioconda::blast'

    input:
		file(derep_fasta) from derep_fasta_blast_db
    output:
        file("blast.db") into cleaned_derep_fasta
    script:
    """
    	makeblastdb -in $derep_fasta -dbtype nucl -out blast.db/db
    """
}

derep_fasta_blastn_query
	.splitFasta(by: 1000, file:true)
	.set { derep_1000 }

process pcr_dup_removal_blastn {
//	conda 'bioconda::blast'
	publishDir "${params.output_dir}/",
	mode: "copy"

    input:
		file("blast.db") from cleaned_derep_fasta
		each file(derep_fasta_1000) from derep_1000
    output:
        file("blast.${derep_fasta_1000.baseName}.csv") into pcr_dup_removal_blastn_results
    script:
    """
    	blastn -task blastn-short -query $derep_fasta_1000 -num_threads 1 -evalue 0.05  -strand plus -db blast.db/db -outfmt 6 -gapopen 0 -gapextend 4 -penalty -3 -reward 2 > blast.${derep_fasta_1000.baseName}.csv 
    	#-word_size 5 -gapopen 0 -gapextend 4 -penalty -3 -reward 2
    """
}

process pcr_dup_removal_blastn_concat_results {
//	conda 'bioconda::blast'
	publishDir "${params.output_dir}/",
	mode: "copy"

    input:
		file(blast_csv) from pcr_dup_removal_blastn_results.collect()
    output:
        file("blastn_results.csv") into pcr_dup_removal_blastn_results_concat
    script:
    """
		find -L . -wholename './blast.*.csv' | sort | xargs cat > blastn_results.csv
    """
}




process pcr_dup_removal_all {
//	conda 'pandas networkx'
	publishDir "${params.output_dir}/",
	mode: "copy"

    echo true
    input:
    	file(derep) from derep_fasta_pcr_removal
		file(blastn_csv) from pcr_dup_removal_blastn_results_concat
    output:
         file("aptamers.clean.fasta") into derep_pcr_removed
    script:
    """
		remove_pcr_duplicates.py -f $derep -b $blastn_csv -o aptamers.clean.fasta
    """
}

"""
========================================================
Analysing SELEX Enrichment
========================================================
"""


fasta_files_filtered
	.combine(derep_pcr_removed)
	.set { fasta_files_pcr_removed }


process count_kmers {
//    conda 'numpy'
    publishDir "${params.output_dir}/",
        pattern: '*{csv}',
        mode: "copy"
                
    input:
        tuple val(round_id), val(round_name), file(fasta), file(cleaned_derep) from fasta_files_pcr_removed
    output:
        tuple val(round_id), val(round_name), file(fasta), file("kmers_${round_name}.csv") into kmer_counts
        
    """
        kmer_count.py -i $fasta --cleaned $cleaned_derep -k ${params.k} -u > kmers_${round_name}.csv
    """
}

kmer_counts
    .into { kmer_counts1; kmer_counts2 }

kmer_counts1
    .combine(kmer_counts2)
    .filter { (it[0] < it[4]) }
    .set{ fisher_kmer_combinations }



process kmer_fisher_test {
//	conda 'scipy'
    publishDir "${params.output_dir}/",
	mode: "copy"

    //echo true
    input:
		tuple val(id1), val(name1), file(fasta1), file(kmers1), val(id2), val(name2), file(fasta2), file(kmers2) from fisher_kmer_combinations
    output:
		tuple val(name1), val(name2), file("fisher_kmer_${name1}_${name2}.csv") into fisher_kmer_channel
    script:
    """
		fisher_testing.py -i1 $kmers1 -i2 $kmers2 -k ${params.k} > fisher_kmer_${name1}_${name2}.csv
    """
}


derep_csv = Channel.fromPath(params.derep_csv)
fisher_kmer_channel
	.combine(derep_fasta_scoring)
	.combine(derep_csv)
	.into { fisher_aptamer_channel; fisher_aptamer_channel2 }




process aptamer_scoring {
//	conda 'anaconda::scipy pandas'
    publishDir "${params.output_dir}/",
	mode: "copy"

    //echo true
    input:
		tuple val(name1), val(name2), file(fisher_kmers), file(derep_fasta), file(derep_csv) from fisher_aptamer_channel
    output:
		tuple val(name1), val(name2), file("scores_${name1}_${name2}.csv") into scoring_channel
    script:
    """
		score_aptamers.py -i $fisher_kmers -k ${params.k} -f $derep_fasta -n1 ${name1} -n2 ${name2} > scores_${name1}_${name2}.csv
    """
}

/*
process aptamer_scoring_2 {
conda 'scipy pandas'
    publishDir "${params.output_dir}/",
	mode: "copy"

    //echo true
    input:
		tuple val(name1), val(name2), file(fisher_kmers), file(derep_fasta), file(derep_csv) from fisher_aptamer_channel2
    output:
		tuple val(name1), val(name2), file("aptamer_fisher_${name1}_${name2}.csv") into scoring_channel2
    script:
    """
		fisher_aptamer_testing.py -i $derep_csv -r1 $name1 -r2 $name2 > aptamer_fisher_${name1}_${name2}.csv
    """
}*/


