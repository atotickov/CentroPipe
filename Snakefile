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
gff_tab_dir_path=results_dir_path/config["gff_tab_dir"]
fasta_dir_path=results_dir_path/config["fasta_dir"]


def create_wildcards_dictionary():
  '''
  Creating a dictionary where the keys are identifiers of genome assemblies,
  and the values of those keys are lists of identifiers for only chromosome scaffolds 
  corresponding to the respective genome assembly. Only scaffolds listed in 
  the input .whitelist files for each genome assembly are taken into account.
  '''
  if "genome_assembly_ids" not in config:
    config["genome_assembly_ids"] = [f.stem for f in genomes_dir_path.iterdir() if f.is_dir()]
  dictionary={} # wildcards_dictionary: 'assembly_id': ['scaffold_ids']
  if "assembly_scaffold_ids" not in config:
    for i in range(len(config["genome_assembly_ids"])):
      assembly_id=config["genome_assembly_ids"][i]
      whitelist_path=genomes_dir_path/f"{assembly_id}/{assembly_id}.whitelist"
      with open(whitelist_path, 'r') as file:
        scaffold_ids=[line.strip() for line in file]
      dictionary[assembly_id]=scaffold_ids
  config["assembly_scaffold_ids"] = dictionary
  return dictionary

wildcards_dictionary = create_wildcards_dictionary()

def get_quartet_params(quartet_params, global_variable=False) -> str:
  '''
  Obtaining all quarTeT program parameters in two formats: 
  either for shell usage or for substitution into file names.
  '''
  params_order = ['Match', 'Mismatch', 'Delta', 'PM', 'PI', 'Minscore', 'l', 'MinPeriod', 'MaxPeriod']
  params = {p: '' for p in params_order}
  for params_dictionary in quartet_params: # yaml dictionary format print: [{'Match': 2}, {'Mismatch': 7}, {'Delta': 7}, {'PM': 80}, {'PI': 10}, {'Minscore': 500}, {'MaxPeriod': 750}, {'MinPeriod': 650}, {'l': 10}]
    for key, value in params_dictionary.items():
      params[key] = str(value)
  if global_variable:
    all_params = [f"{params['Match']}.{params['Mismatch']}.{params['Delta']}.{params['PM']}.{params['PI']}.{params['Minscore']}.{params['l']}.{params['MinPeriod']}.{params['MaxPeriod']}"]
    return all_params
  else:
    quartet_params_str = [f"--trf {params['Match']} {params['Mismatch']} {params['Delta']} {params['PM']} {params['PI']} {params['Minscore']} -r {params['l']} -n {params['MinPeriod']} -m {params['MaxPeriod']}"]
    return ' '.join(quartet_params_str)

all_quartet_params = get_quartet_params(quartet_params=config["quartet_params"], global_variable=True)

def expand_result_template_list(wildcards_dictionary, template):
  '''
  Function that utilizes the provided template with curly braces for substituting available wildcards. 
  It returns a list of strings containing all combinations of the given template.
  '''
  template_list=[]
  for assembly_id, scaffold_ids in wildcards_dictionary.items():
    for scaffold_id in scaffold_ids:
      template_list+=expand(template, assembly_id=assembly_id, scaffold_id=scaffold_id, all_quartet_params=all_quartet_params)
  return template_list

localrules: all

rule all:
  input:
    expand_result_template_list(wildcards_dictionary, genomes_dir_path/"{assembly_id}/{assembly_id}.filtered.fasta"),
    expand_result_template_list(wildcards_dictionary, quartet_dir_path/"{assembly_id}/tmp/{assembly_id}.centro.chr.txt"),
    expand_result_template_list(wildcards_dictionary, quartet_dir_path/"{assembly_id}/tmp/{assembly_id}.centro.chr.bed"),
    expand_result_template_list(wildcards_dictionary, gff_tab_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.raw.dat.gff.gz"),
    expand_result_template_list(wildcards_dictionary, gff_tab_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.raw.dat.gff.tab"),
    expand_result_template_list(wildcards_dictionary, gff_tab_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.filtered.dat.gff.tab"),
    expand_result_template_list(wildcards_dictionary, gff_tab_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.raw.dat.gff.stats.tab"),
    expand_result_template_list(wildcards_dictionary, gff_tab_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.filtered.dat.gff.stats.tab"),
    expand_result_template_list(wildcards_dictionary, fasta_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.filtered.fasta"),
    expand_result_template_list(wildcards_dictionary, fasta_dir_path/"{assembly_id}/{scaffold_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.rotated.fasta"),
    #expand(genomes_dir_path/"{assembly_id}/{assembly_id}.filtered.fasta",assembly_id=config["assembly_scaffold_ids"].keys()),
    # expand_result_template_list(wildcards_dictionary, gff_tab_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.YES"),
    # [[expand(gff_tab_dir_path / "{assembly_id}/{assembly_id}.{scaffold_id}.fasta.2.7.7.80.10.500.750.dat.gff.gz", assembly_id=[assembly_id], scaffold_id=[scaffold_id]) for scaffold_id in dictionary[assembly_id]] for assembly_id in dictionary]

include: "workflow/rules/quartet.smk"
