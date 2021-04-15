#!/usr/bin/env nextflow
"""
========================================================
Groovy Helper Functions
========================================================    
"""
def reverse_complement(String s) {
    complement(s.reverse());
}

def complement(String s) {
    def acgt_map = [
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "a": "t",
        "c": "g",
        "g": "c",
        "t": "a"
    ];

    char[] sc = new char[s.length()];
    for (int i = 0; i < s.length(); i++) {
        sc[i] = acgt_map[s[i]];
    }
    new String(sc);
}

def remove_all_extensions(String s) {
    s.substring(0, s.indexOf("."));
}

def trim_fn(String s) {    
    if (params.trim_filenames) {
        return s.substring(0, s.indexOf(params.trim_delimiter));
    } else {
        return s;
    }
}

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
   

"""
========================================================
Analysing SELEX Enrichment
========================================================
"""

process count_kmers {
    publishDir "${params.output_dir}/",
        pattern: '*{csv}',
        mode: "copy"
                
    input:
        tuple val(round_id), val(round_name), file(fasta) from fasta_files_filtered
    output:
        tuple val(round_id), val(round_name), file(fasta), file("kmers_${round_name}.csv") into kmer_counts
        
    """
        kmer_count.py -i $fasta -k ${params.k} -u > kmers_${round_name}.csv
    """
}

kmer_counts
    .into { kmer_counts1; kmer_counts2 }

kmer_counts1
    .combine(kmer_counts2)
//    .map { it -> [it[0][0], it[0][1], it[0][2], it[0][3],
//		  it[1][0], it[1][1], it[1][2], it[1][3]] } 
//    .flatten()
    .filter { (it[0] < it[4]) }
    .view()
    .set{ fisher_kmer_combinations }

derep_fasta = Channel.fromPath(params.derep_fasta)
fisher_kmer_combinations_derep = fisher_kmer_combinations
	.combine(derep_fasta)

process kmer_fisher_test {
    publishDir "${params.output_dir}/",
	mode: "copy"

    echo true
    input:
	tuple val(id1), val(name1), file(fasta1), file(kmers1), val(id2), val(name2), file(fasta2), file(kmers2), file(derep) from fisher_kmer_combinations_derep
    output:
	tuple val(name1), val(name2), file("fisher_${name1}_${name2}.csv"), file("scores_${name1}_${name2}.csv") into fisher_channel
    script:
    """
	fisher_testing.py -i1 $kmers1 -i2 $kmers2 -k ${params.k} > fisher_${name1}_${name2}.csv
	score_aptamers.py -i fisher_${name1}_${name2}.csv -k ${params.k} -f $derep > scores_${name1}_${name2}.csv
    """
}


