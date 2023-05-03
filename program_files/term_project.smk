import json

rule all:
    input:
        config["target"]

rule ncbi_download:
    output:
        "assembly_files/{id}.fna"
    params:
        accession = lambda wildcards: config["seqs_to_use"][wildcards.id]["accession"],
        id = "{id}"
    conda:
        "snakemake/conda/ncbi_datasets.yaml"
    shell:
        """
        mkdir ncbi_temp/{params.accession}
        cd ncbi_temp/{params.accession}
        datasets download genome accession {params.accession}
        unzip ncbi_dataset.zip
        cd ../..
        mv ncbi_temp/{params.accession}/ncbi_dataset/data/{params.accession}/{params.accession}_{params.id}_genomic.fna assembly_files/{params.id}.fna
        rm -r ncbi_temp/{params.accession}
        """

rule fraggenescan:
    input:
        "assembly_files/{id}.fna"
    output:
        "fraggenescan_output/{id}.fragscan.faa"
    params:
        output = "{id}.fragscan"
    threads: 10
    conda:
        "snakemake/conda/fraggenescan.yaml"
    shell:
        """
        mkdir fragscan_{wildcards.id}
        FragGeneScan -s {input} -o fragscan_{wildcards.id}/{params.output} -w 1 -t complete -p {threads}
        mv fragscan_{wildcards.id}/{params.output}.faa fraggenescan_output/
        rm -r fragscan_{wildcards.id}
        """

rule diamond:
    input:
        "fraggenescan_output/{id}.fragscan.faa"
    output:
        "diamond_output/{id}_fragscan-diamond.out"
    threads: 10
    conda:
        "snakemake/conda/diamond.yaml"
    shell:
        "diamond blastp -p {threads} -q {input} --db db_files/RepeatPeps.lib.dmnd > {output}"

rule reindex:
    input:
        "diamond_output/{id}_fragscan-diamond.out"
    output:
        "diamond_output/{id}_diamond-hits.bed"
    shell:
        "python3 program_files/fragscan_index.py {input} {output}"

rule get_fasta_windows:
    input:
        diamond_hits = "diamond_output/{id}_diamond-hits.bed",
        genomic = "assembly_files/{id}.fna"
    output:
        out = "te_windows/{id}_windows.fasta",
        report = "logs_reports/{id}_windows.out"
    params:
        window_size = 10000
    conda:
        "snakemake/conda/python.yaml"
    shell:
        "python3 program_files/ltr-window.py diamond_window -i {input.diamond_hits} -f {input.genomic} -o {output.out} -r {output.report} -w {params.window_size}"

#For non-chromosome level assemblies, if there is time:
'''
└ $ bwa mem ref_calbicans_SC5314.fasta ASM2576641v1_windows.fasta -t 32 > bwa_test.sam
└ $ samtools view -S -b bwa_test.sam > bwa_test.bam
└ $ bedtools bamtobed -i bwa_test.bam > bwa_test.bed
'''

rule group:
    input:
        expand("logs_reports/{ids}_windows.out", ids=config["seqs_to_use"])
    output:
        "align.sh"
    conda:
        "snakemake/conda/biopython.yaml"
    shell:
        "python3 program_files/te_group.py {input}"

rule alignment:
    input:
        "align.sh"
    output:
        "align_out.txt"
    threads: 15
    conda:
        "snakemake/conda/clustal_omega.yaml"
    shell:
        "bash align.sh > {output}"

rule ltr_finder:
    input:
        "te_windows/{id}_windows.fasta"
    output:
        "ltr_finder/{id}_windows.ltr"
    conda:
        "snakemake/conda/LTR_finder.yaml"
    shell:
        "ltr_finder {input} > {output}"

rule te_seqs:
    input:
        ltr_finder = "ltr_finder/{id}_windows.ltr",
        gdna = "assembly_files/{id}.fna"
    output:
        "te_seqs/{id}_TEs.fasta"
    conda:
        "snakemake/conda/python.yaml"
    shell:
        "python3 program_files/ltr-window.py te_seqs -i {input.ltr_finder} -g {input.gdna} -o {output}"

rule cat_TEs:
    input:
        expand("te_seqs/{ids}_TEs.fasta", ids=config["seqs_to_use"])
    output:
        "te_seqs/all_TEs.fasta"
    shell:
        "cat {input} > {output}"

rule blast_ref:
    input:
        "te_seqs/all_TEs.fasta"
    output:
        "blast/all_TEs-blastvref.out"
    conda:
        "snakemake/conda/blast.yaml"
    threads: 10
    shell:
        "blastn -query {input} -db db_files/goodwin_lib.fasta -out {output} -num_threads {threads}"

rule self_db:
    input:
        "te_seqs/all_TEs.fasta"
    output:
        "blast/all_TEs.fasta",
        "blast/all_TEs.fasta.nhr",
        "blast/all_TEs.fasta.nsq",
        "blast/all_TEs.fasta.nin"
    conda:
        "snakemake/conda/blast.yaml"
    shell:
        """
        cp {input} blast/all_TEs.fasta
        makeblastdb -in blast/all_TEs.fasta -dbtype nucl
        """

rule blast_self:
    input:
        seqs = "te_seqs/all_TEs.fasta",
        db_file = "blast/all_TEs.fasta.nhr"
    output:
        "blast/all_TEs-blastvself.out"
    conda:
        "snakemake/conda/blast.yaml"
    threads: 30
    shell:
        "blastn -query {input.seqs} -db {input.seqs} -out {output} -num_threads {threads}"

rule cross_check:
    input:
        vs_library = "blast/all_TEs-blastvref.out",
        vs_self = "blast/all_TEs-blastvself.out"
    output:
        "logs_reports/TEs_found.txt"
    conda:
        "snakemake/conda/python.yaml"
    shell:
        "python3 program_files/ltr-window.py blast_output --vs_library {input.vs_library} --vs_self {input.vs_self} -c program_files/snakemake/config.json > {output}"

rule download_reference_tes:
    output:
        "reference_tes/{ref}.fasta"
    params:
        url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id={ref}&extrafeat=null&fmt_mask=0&retmode=html&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=1000000"
    shell:
        '''
        wget -O {output} "{params.url}"
        '''

rule download_reference_genbank:
    output:
        "reference_tes/{ref}.gbff"
    shell:
        '''
        curl -o {output} -s "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&sendto=on&log$=seqview&db=nuccore&sort=&id={ref}&from=begin&to=end&maxplex=3"
        '''

rule te_family_fasta:
    input:
        reference = "reference_tes/{ref}.fasta",
        all_seqs = "te_seqs/all_TEs.fasta"
    output:
        "phylogeny/{ref}_family.fasta"
    params:
        config = lambda wildcards: f"program_files/snakemake/config.json"
    conda:
        "snakemake/conda/python.yaml"
    shell:
        "python3 program_files/ltr-window.py fasta_family -i {wildcards.ref} --ref {input.reference} -c {params.config} -o {output}"

rule unknown_family_fasta:
    input:
        "te_seqs/all_TEs.fasta"
    output:
        "phylogeny/family_{unknown}.fasta"
    conda:
        "snakemake/conda/python.yaml"
    shell:
        "python3 program_files/ltr-window.py fasta_family -i {wildcards.unknown} -c program_files/snakemake/config.json -o {output}"

rule align_family:
    input:
        "phylogeny/{family_fasta}.fasta"
    output:
        "phylogeny/{family_fasta}.aln"
    threads: 2
    conda:
        "snakemake/conda/clustal_omega.yaml"
    shell:
        "clustalo -i {input} -o {output} --threads={threads}"

rule iqtree:
    input:
        "phylogeny/{family_fasta}.aln"
    output:
        "phylogeny/{family_fasta}.aln.contree"
    conda:
        "snakemake/conda/iqtree.yaml"
    threads: 4
    shell:
        """
        iqtree -s {input} -m TEST -nt {threads} -redo >> phylogeny/iq_tree_logs
        iqtree -s {input} -m $(awk '$1 ~ /^Best-fit/ {{print $6}}' {input}.iqtree) -bb 1000 -nt {threads} -redo >> phylogeny/iq_tree_logs
        """

rule related_ltrs:
    output:
        "phylogeny/{family_fasta}_group.fasta"
    shell:
        """
        group=$(awk -F. '{{print $1}}' <<< {output} | awk -F/ '{{print $2}}')
        output=$(echo phylogeny/${{group}}.fasta)
        touch $output
        count=$(awk -F_ '{{print NF}}' <<<$group)
        for ((i=1;i<=$count;i++))
        do 
            family=$(awk -F_ '{{print $'$i'}}' <<< $group)
            accession=$(awk '$0 ~ / '$family'/ {{print substr($1,2)}}' reference_tes/*)
            familyfile=$(echo phylogeny/${{accession}}_family.fasta)
            if [ ${{familyfile}} != "phylogeny/_family.fasta" ]; then 
                cat $familyfile >> $output
            fi
        done
        """

"""See all trees:
awk '$1 ~ /^\+/ || $1 ~ /^\|/' phylogeny/*.iqtree

touch phylogeny/whio_tara_titi.fasta
for tefamily in whio tara titi
do
    cat phylogeny/$(awk '$0 ~ / tefamily/ {print substr($1,2)}' reference_tes/*)_family.fasta >> whio_tara_titi.fasta
"""
