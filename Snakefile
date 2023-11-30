from pathlib import Path
import os

# Setup config
configfile: "config/config.yaml"

# Setup paths
benchmark_dir_path=Path(config["benchmark_dir"])
cluster_log_dir_path=Path(config["cluster_log_dir"])
log_dir_path = Path(config["log_dir"])

genomes_dir_path=Path(config["genomes_dir"])
results_dir_path=Path(config["results_dir"])
quartet_dir_path=results_dir_path/config["quartet_dir"]
chromosomes_dir_path=results_dir_path/config["chromosomes_dir"]

if "genome_assembly_ids" not in config:
  config["genome_assembly_ids"] = [f.stem for f in genomes_dir_path.iterdir() if f.is_dir()]

# Создание словаря, где ключи - это идентификаторы сборок геномов, а их значения - это список их хромосомных скаффолдов по файлу .whitelist 
dictionary={}
if "assembly_scaffold_ids" not in config:
  for i in range(len(config["genome_assembly_ids"])):
    assembly_id=config["genome_assembly_ids"][i]
    whitelist_path=genomes_dir_path/f"{assembly_id}/{assembly_id}.whitelist"
    with open(whitelist_path, 'r') as file:
      scaffold_ids=[line.strip() for line in file]
    dictionary[assembly_id]=scaffold_ids
config["assembly_scaffold_ids"] = dictionary
#print(config["assembly_scaffold_ids"])
#print(config["assembly_scaffold_ids"].keys())
#print(config["assembly_scaffold_ids"].values())

def expand_result_template_list(dictionary, template):
    template_list=[]
    for assembly_id, scaffold_ids in dictionary.items():
        for scaffold_id in scaffold_ids:
            template_list+=expand(template, assembly_id=assembly_id, scaffold_id=scaffold_id)
    #print(template_list)
    return template_list

output_files = [
  expand(quartet_dir_path/"{assembly_id}/tmp/{assembly_id}.centro.chr.txt",assembly_id=config["assembly_scaffold_ids"].keys()),
  expand(genomes_dir_path/"{assembly_id}/{assembly_id}.filtered.fasta",assembly_id=config["assembly_scaffold_ids"].keys()),
  expand_result_template_list(dictionary, chromosomes_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.fasta.2.7.7.80.10.500.750.dat.gff.gz"),
  expand_result_template_list(dictionary, chromosomes_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.fasta.2.7.7.80.10.500.750.dat.gff.tab")
  ]

localrules: all

rule all:
  input:
    output_files

include: "workflow/rules/quartet.smk"