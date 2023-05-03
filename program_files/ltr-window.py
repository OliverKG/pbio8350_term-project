import sys
import subprocess
import json

def diamond_window(args):
    params = {
        "-i":"",#input (bed format)
        "-f":"",#genomic sequence (fasta format)
        "-o":False,#output path
        "-r":False,#report path
        "-w":6000 #window size
        } 
    params = get_params(params,args)

    diamond_out = params["-i"]
    genome_assembly = params["-f"]
    window_size = int(params["-w"])
    window_radius = int(window_size/2)

    f = open(diamond_out)
    diamond_out = f.read().strip().split("\n")
    f.close()

    chr = {}
    for line in diamond_out:
        line = line.split("\t")
        start = int(line[1])
        end = int(line[2])
        if line[0] not in chr:
            chr[line[0]] = [[start-window_radius, end+window_radius, [",".join([str(start),str(end),line[3],line[4],line[5]])]]]
        else:
            novel = True
            for x in range(len(chr[line[0]])):
                if(start > chr[line[0]][x][0] and end < chr[line[0]][x][1]):
                    novel = False
                    if(start-window_radius < chr[line[0]][x][0]): chr[line[0]][x][0] = start-window_radius
                    if(end+window_radius > chr[line[0]][x][1]): chr[line[0]][x][1] = end+window_radius
                    chr[line[0]][x][2].append(",".join([str(start),str(end),line[3],line[4],line[5]]))
            if novel == True:
                chr[line[0]].append([start-window_radius, end+window_radius, [",".join([str(start),str(end),line[3],line[4],line[5]])]])

    def strain_name(namelist):
        namelist = namelist.split("\\n")
        name = namelist[0]
        name = name.split(",")[0]
        name = name.split("bctg")[0].strip()
        name = " ".join(name.split(" ")[:-1])
        patternlist = {
            "albicans":1,
            "strain":1,
            "Chromosome":0,
            "chromosome":0,
            "ctg":0,
            "Contig":0,
            "contig":0
        }
        for pattern in patternlist:
            name = name.split(pattern)
            name = name[min(len(name)-1,patternlist[pattern])]
        return name.strip().replace(" ","-")

    def contig(line):
        line = line.split("\\n")[0]
        contig_type = "contig"
        if("hromosome" in line): 
            contig_type = "chromosome"
            line = line.split("hromosome")[1].strip()
            name = "chr" + line.split(" ")[0].replace(",","")
        else: name = line.split(",")[0].split(" ")[-1]
        if "bctg" in name:
            name = " ".join(line.split(",")[0].split(" ")[-2:])
        return contig_type, name

    with get_output(params) as output, get_report(params) as report:
        strain = strain_name(awk("^>",genome_assembly))
        report.write("strain\tcontig_accession\tcontig_type\tcontig_name\tstart\tend\n")
        for chromosome in chr:
            sequence = seq(chromosome,genome_assembly)
            for line in chr[chromosome]:
                contig_type, contig_name = contig(awk(f"^>{chromosome}",genome_assembly))
                start = max(0,line[0])
                end = min(len(sequence)-1,line[1])
                seq_name = f"{chromosome}_{start}_{end}"
                all_data = ';'.join(line[2])
                if("LTR" in all_data):
                    output.write(f">{seq_name}\n{sequence[start:end+1]}\n")
                    report.write(f"{strain}\t{chromosome}\t{contig_type}\t{contig_name}\t{start}\t{end}\n")

def te_seqs(args):
    params = {
        "-i":"",#input (ltr)
        "-g":"",#genomic fasta
        "-o":False#output path
        } 
    params = get_params(params,args)

    input_file = open(params["-i"],"r")
    ltr_string = input_file.read()
    input_file.close()

    chromosomes = {}

    class Prediction:
        def __init__(self, ltr_data):
            self.data = ltr_data.split("\n")
            self.window = self.data[1][self.data[1].index(":")+1:self.data[1].index("Len")].strip()
            self.chromosome = self.window[:self.window.rindex("_")]
            self.chromosome = self.chromosome[:self.chromosome.rindex("_")]
            self.get_info()
            self.locate()
            if self.chromosome not in chromosomes:
                chromosomes[self.chromosome] = set()
            chromosomes[self.chromosome].add(self)
        
        def get_info(self):
            info = {
                "location":"Location",
                #"score":"Score",
                #"status":"Status",
                #"5p_ltr":"5'-LTR",
                #"3p_ltr":"3'-LTR",
                #"5p_tg":"5'-TG",
                #"3p_ca":"3'-CA",
                #"tsr":"TSR",
                #"sharpness":"Sharpness",
                #"strand":"Strand +",
                #"ppt":"PPT"
                }
            for line in self.data:
                for info_sought in info:
                    if info[info_sought] in line:
                        info[info_sought] = line[line.index(":")+1:].strip()
            self.__dict__.update(info)
        
        def locate(self):
            seq_start = int(self.window.split("_")[-2])
            subseq_start = int(self.location[:self.location.index("-")].strip())
            subseq_end = int(self.location[self.location.index("-")+1:self.location.index("Len")].strip())
            self.start = seq_start+subseq_start
            self.end= seq_start+subseq_end
            self.name = f"{self.chromosome}_{self.start}_{self.end}"

    predictions = [Prediction(sequence) for sequence in ltr_string.split("Predict protein Domains") if "Score" in sequence]

    for chromosome in chromosomes:
        sequence = seq(chromosome, params["-g"])
        for prediction_obj in chromosomes[chromosome]:
            prediction_obj.sequence = sequence[prediction_obj.start:prediction_obj.end+1]
            prediction_obj.fasta = f">{prediction_obj.name}\n{prediction_obj.sequence}\n"
    
    with get_output(params) as output:
        for prediction_obj in predictions:
            output.write(prediction_obj.fasta)

def blast_output(args):
    params = {
        "--vs_library":"",#blast vs library
        "--vs_self":"",#blast vs self
        "-c":""#config file
        } 
    params = get_params(params,args)

    input_file = open(params["--vs_library"],"r")
    vs_library = input_file.read()
    input_file.close()
    input_file = open(params["--vs_self"],"r")
    vs_self = input_file.read()
    input_file.close()

    library_seqs = {}

    class te_sequence:
        def __init__(self, library_blast):
            self.raw_out = library_blast.split("\n")
            self.name = self.raw_out[0].strip()
            hits = library_blast[library_blast.index("Sequences producing significant alignments"):]
            hits = hits[hits.index("\n"):hits.index(">")].split("\n")
            self.hits = {}
            self.cross_hits = {}
            for row in hits:
                hit = row.strip().split()
                if len(hit) > 0:
                    hit = hit[0]
                    self.lib_hit(hit, library_blast)

        def lib_hit(self, hit, info):
            self.hits[hit] = self.hit(hit, info)
            if not hit in library_seqs: library_seqs[hit] = set()
            library_seqs[hit].add(self)

        def hit(self, hit, info_dump):
            hit_info = info_dump[info_dump.index(">")+1:]
            hit_info = hit_info[hit_info.index(hit):]
            if ">" in hit_info: hit_info = hit_info[:hit_info.index(">")]
            hit_array = {
                "score":"Score =",
                "e":"Expect =",
                "identity":"Identities = ",
                "gaps":"Gaps ="
                }
            for info in hit_array:
                info_list = []
                raw_info = hit_info.split(hit_array[info])[1:]
                for entry in raw_info:
                    entry = entry.strip().split(",")[0]
                    info_list.append(entry.split("\n")[0].strip())
                hit_array[info] = info_list
            hit_array["description"] = hit_info[len(hit):hit_info.index("\n")].strip()
            return hit_array

        def cross_hit(self, cross_ref):
            hits = cross_ref[cross_ref.index("Sequences producing significant alignments"):]
            hits = hits[hits.index("\n"):hits.index(">")].split("\n")
            for row in hits:
                hit = row.strip().split()
                if len(hit) > 0:
                    hit = hit[0]
                    if not hit == self.name: self.cross_hits[hit] = self.hit(hit, cross_ref)
        
        def print(self):
            print(f" - {self.name}")
            if len(self.cross_hits) > 0:
                for hit in self.cross_hits:
                    print(f"    - {hit}")

    te_seqs = [te_sequence(result) for result in vs_library.split("Query=") if "significant alignments" in result]
    te_seqs = {te.name:te for te in te_seqs}
    cross_ref = [result for result in vs_self.split("Query=") if "significant alignments" in result]
    cross_ref = {ref.strip().split()[0]:ref for ref in cross_ref}
    unidentified_TEs = {ref:cross_ref[ref] for ref in cross_ref if not ref in te_seqs}
    unidentified_TE_groups = []
    unid_keys = list(unidentified_TEs.keys())
    while len(unid_keys) > 0:
        search_TErm = unid_keys.pop(0)
        content = unidentified_TEs[search_TErm]
        group = {search_TErm}
        for x in range(len(unid_keys)):
            te = unid_keys[x+1-len(group)]
            if te in content:
                group.add(te)
                unid_keys.remove(te)
        unidentified_TE_groups.append(group)
    unidentified_TE_groups = sorted(unidentified_TE_groups, key=len, reverse=True)

    #for seq in te_seqs:
        #if seq in cross_ref:
        #    te_seqs[seq].cross_hit(cross_ref[seq])
    output = {"known_families":{},"unknown_families":{}}
    for ref in library_seqs:
        output["known_families"][ref] = []
        print(ref)
        for te in library_seqs[ref]:
            output["known_families"][ref].append(te.name)
            te.print()
    
    ungrouped = False
    for x in range(len(unidentified_TE_groups)):
        if(len(unidentified_TE_groups[x]) > 1):
            output["unknown_families"][f"unknown_{x+1}"] = list(unidentified_TE_groups[x])
            print(f"Unidentified TE group {x+1}")
            for te in unidentified_TE_groups[x]:
                print(f" - {te}")
        else:
            if not ungrouped:
                print(f"Unidentified TEs with no matches:")
                ungrouped = True
            print(f" - {unidentified_TE_groups[x].pop()}")
    
    with open(params["-c"]) as json_file:
        old_config = json.load(json_file)
    output.update(old_config)

    output["chromosome_at"] = {}
    for genome in output["seqs_to_use"]:
        chromosome_list = awk("^>",f"assembly_files/{genome}.fna").replace(">","").split("\\n")
        chromosome_list = [chromosome.split()[0].replace("b'","") for chromosome in chromosome_list]
        for chromosome in chromosome_list:
            output["chromosome_at"][chromosome] = genome

    with open(params["-c"],"w") as output_file:
        json.dump(output, output_file)

def fasta_family(args):
    params = {
        "--ref":"",#reference
        "-i":"",#identity
        "-c":"",#config file
        "-o":""#output
        }
    params = get_params(params,args)

    with open(params["-c"]) as json_file:
        config = json.load(json_file)

    if params["--ref"] == "":
        sequences = config["unknown_families"]
    else:
        sequences = config["known_families"]
    sequences = sequences[params["-i"]]
    chromosomes = config["chromosome_at"]

    with open(params["-o"],"w") as output:
        if not params["--ref"] == "":
            output.write(str(subprocess.run(f"cat {params['--ref']}", shell=True, capture_output=True).stdout).replace("b'","").replace("\\n'","").replace("\\n","\n")+"\n")
        for sequence in sequences:
            chromosome = sequence[:sequence.rindex("_")]
            chromosome = chromosome[:chromosome.rindex("_")]
            output.write(f">{chromosomes[chromosome]}_{sequence}\n{seq(sequence, 'te_seqs/all_TEs.fasta')}\n")

def window(args):
    with open("program_files/snakemake/config.json") as new_config:
        config = json.load(new_config)

    for genome in config["Chromosome"]:
        file = f"assembly_files/{genome}.fna"
        seq_name = awk(f"hromosome R", file,0)
        seq_name = seq_name[seq_name.index(">")+1:]
        seq_name = seq_name[:seq_name.index(" ")]
        print(f">{genome}_chr5")
        print(seq(seq_name, file)[1000:40000])


def awk(string,file, column=1):
    awk = subprocess.run(f"awk '${column} ~ /{string}/' {file}", shell=True, capture_output=True)
    return str(awk.stdout)

def seq(name, fasta_file):
    sed = subprocess.run(f"sed -n '/{name}/, />/{{ /{name}/! {{ />/! p }} }}' {fasta_file}", shell=True, capture_output=True)
    return str(sed.stdout).replace("\\n","").replace("b'","").replace("'","")

class get_output:
    def __init__(self, params):
        self.params = params
        pass
    def __enter__(self):
        if self.params["-o"] == False:
            return console
        else: return open(self.params["-o"],"w")
    def __exit__(self, *args):
        pass
    
class get_report:
    def __init__(self, params):
        self.params = params
        pass
    def __enter__(self):
        if self.params["-r"] == False:
            if self.params["-o"] == False:
                return void
            else: return console
        else: return open(self.params["-r"],"w")
    def __exit__(self, *args):
        pass

class console:
    def __init__(self):
        pass
    def write(string):
        print(string.strip())

class void:
    def __init__(self):
        pass
    def write(string):
        pass

def get_params(param_dict,args):
    for x in range(len(args)):
        if args[x] in param_dict: pass
        elif args[x-1] in param_dict:
            param_dict[args[x-1]] = args[x]
    return param_dict

if __name__ == '__main__':
    globals()[sys.argv[1]](sys.argv[2:])
