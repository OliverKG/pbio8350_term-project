import sys
import os
import random
import subprocess
import json

params = {
    "-t":40,#threads
    "-q":"Candida albicans",#query
    "--assembly-level":"contig"#minimum assembly level
}
for x in range(len(sys.argv)):
    if sys.argv[x] in params: pass
    elif sys.argv[x-1] in params:
        params[sys.argv[x-1]] = sys.argv[x]

run_id = random.randint(10000,99999)

summary_file = f"summary_{run_id}.json"

subprocess.run(f"datasets summary genome taxon '{params['-q']}' > {summary_file}", shell=True)

with open(summary_file) as json_file:
    summary = json.load(json_file)

config = {"run_id":run_id,"assembly_level":params["--assembly-level"],"seqs_to_use":{}}
for assembly in summary["reports"]:
    accesion = assembly["accession"]
    name = assembly["assembly_info"]["assembly_name"]
    assembly_level = assembly["assembly_info"]["assembly_level"]
    if not assembly_level in config:
        config[assembly_level] = {}
    config[assembly_level][name] = {
        "accession":accesion
    }
    if "infraspecific_names" in assembly["organism"]:
        for identifier in ["strain","isolate"]:
            if identifier in assembly["organism"]["infraspecific_names"]:
                config[assembly_level][name][identifier] = assembly["organism"]["infraspecific_names"][identifier]

for assembly_level in ["Chromosome","Scaffold","Contig"]:
    config["seqs_to_use"].update(config[assembly_level])
    if params["--assembly-level"].lower() == assembly_level.lower():
        break

os.remove(summary_file)

config_file = f"program_files/snakemake/config.json"

def dump(config, config_file, target="same"):
    if not target == "same":
        config["target"] = target
    with open(config_file,"w") as output:
        json.dump(config, output)

def update_from(config_file):
    with open(config_file) as new_config:
        config = json.load(new_config)
    return config

os.remove("logs_reports/TEs_found.txt")
dump(config, config_file, target="logs_reports/TEs_found.txt")
subprocess.run(f"snakemake -s program_files/term_project.smk --configfile={config_file} --use-conda --cores {params['-t']}", shell=True)
config = update_from(config_file)
te_families = {f"{te_family}_family":config["known_families"][te_family] for te_family in config["known_families"]}
te_families.update({f"family_{unknown_family}":config["unknown_families"][unknown_family] for unknown_family in config["unknown_families"]})
fasta_family_output = [f"phylogeny/{te_family}.fasta" for te_family in te_families]
phylogeny_outputs = [f"phylogeny/{te_family}.aln.contree" for te_family in te_families if len(te_families[te_family]) > 2]
related_ltr_groups = [
    #"alpha_san_xi",-problem was just duplicated sequences, could go back and fix
    "lambda_zeta",
    #"omega_mu",
    "nu_iota",
    #"omicron_pi_rho",
    "chi_phi",
    "eta_epsilon",
    "whio_tara_titi",
    #"moa_episemon",
    "sampi_theta"]
related_ltr_outputs = [f"phylogeny/{group}_group.aln.contree" for group in related_ltr_groups]
dump(config, config_file, target = fasta_family_output + phylogeny_outputs + related_ltr_outputs)
subprocess.run(f"snakemake -s program_files/term_project.smk --configfile={config_file} --use-conda --cores {params['-t']}", shell=True)
