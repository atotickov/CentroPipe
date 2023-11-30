ruleorder: scaffold_filtering_by_length > run_quartet > dat_files_conversion

# Правило для фильтрации всех скаффолдов по определенной минимальной длине
rule scaffold_filtering_by_length:
        input:
                fasta=genomes_dir_path/"{assembly_id}"/"{assembly_id}.fasta",
                lenfile=genomes_dir_path/"{assembly_id}"/"{assembly_id}.len"
        output:
                genomes_dir_path/"{assembly_id}/{assembly_id}.filtered.fasta"
        conda:
                "../../%s" % config["data_processing_conda_config"]
        params:
                minlen=config["scaffold_filtering_by_length.params"]
        log:
                std=log_dir_path/"scaffold_filtering_by_length.{assembly_id}.log",
                cluster_log=cluster_log_dir_path/"scaffold_filtering_by_length.{assembly_id}.cluster.log",
                cluster_err=cluster_log_dir_path/"scaffold_filtering_by_length.{assembly_id}.cluster.err"
        benchmark:
                benchmark_dir_path/"scaffold_filtering_by_length.{assembly_id}.benchmark.txt"
        resources:
                threads=config["scaffold_filtering_by_length.threads"],
                time=config["scaffold_filtering_by_length.time"],
                mem_mb=config["scaffold_filtering_by_length.mem_mb"]
        shell:
                "awk '$2 >= {params.minlen} {{print $1}}' {input.lenfile} | xargs samtools faidx {input.fasta} > {output} 2> {log.std} "

# Правило для запуска квартета на сборках геномов с фильтрованными скаффолдами по минимальной длине скаффолдов
rule run_quartet:
        input:
                genomes_dir_path/"{assembly_id}"/"{assembly_id}.filtered.fasta"
        output:
                quartet_outdir=directory(quartet_dir_path/"{assembly_id}"),
                dat_files_outdir=directory(quartet_dir_path/"{assembly_id}/tmp/trfdat"),
                centromere_coordinate_file=quartet_dir_path/"{assembly_id}/tmp/{assembly_id}.centro.chr.txt"
        conda:
                "../../%s" % config["quartet_conda_config"]
        params:
                options=config["quartet.params"],
        log:
                std=log_dir_path/"quartet.{assembly_id}.log",
                cluster_log=cluster_log_dir_path/"quartet.{assembly_id}.cluster.log",
                cluster_err=cluster_log_dir_path/"quartet.{assembly_id}.cluster.err"
        benchmark:
                benchmark_dir_path/"quartet.{assembly_id}.benchmark.txt"
        resources:
                threads=config["quartet.threads"],
                time=config["quartet.time"],
                mem_mb=config["quartet.mem_mb"]
        shell:
                "MYPWD=$(pwd); "
                "cd {output.quartet_outdir}; "
                "quartet.py CentroMiner -i $MYPWD/{input} -t {resources.threads} -p {wildcards.assembly_id} {params.options} > $MYPWD/{log.std} 2>&1 "

# Правило для конвертации .dat файлов в .gff и .gff.tab. Конвертируются лишь скаффолды, перечисленные в whitelist сборки
rule dat_files_conversion:
        input:
                dat_files_outdir=quartet_dir_path/"{assembly_id}/tmp/trfdat"
        output:
                gff_files=chromosomes_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.fasta.2.7.7.80.10.500.750.dat.gff.gz",
                tab_files=chromosomes_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.fasta.2.7.7.80.10.500.750.dat.gff.tab"
        conda:
                "../../%s" % config["data_processing_conda_config"]
        params:
                options=config["chromosomes.params"],
                chromosomes_outdir=directory(chromosomes_dir_path),
        log:
                std=log_dir_path/"chromosomes.{assembly_id}.{scaffold_id}.log",
                cluster_log=cluster_log_dir_path/"chromosomes.{assembly_id}.{scaffold_id}.cluster.log",
                cluster_err=cluster_log_dir_path/"chromosomes.{assembly_id}.{scaffold_id}.cluster.err"
        benchmark:
                benchmark_dir_path/"chromosomes.{assembly_id}.{scaffold_id}.benchmark.txt"
        resources:
                threads=config["chromosomes.threads"],
                time=config["chromosomes.time"],
                mem_mb=config["chromosomes.mem_mb"]
        shell:
                "mkdir -p {params.chromosomes_outdir}; "
                "for i in {input}/{wildcards.assembly_id}.{wildcards.scaffold_id}.fasta.2.7.7.80.10.500.750.dat; do TRF.py -i $i -o {output.gff_files}; done > {log.std} 2>&1; "
                "GFF_to_TAB.py {output.gff_files} {output.tab_files} > {log.std} 2>&1 "













# rule processed:
#         input:
#                 dat_files_outdir=quartet_dir_path/"{assembly_id}/tmp/trfdat"
#         output:
#                 processed_files=quartet_dir_path/"{assembly_id}/tmp/trfdat/{assembly_id}.{scaffold_id}.fasta.2.7.7.80.10.500.750.PROCESSED",
#         shell:
#                 "[-f "'{input}'" ] && touch "'{output.OK}'" || {echo "'ERROR! {input} not found.'"; exit 1; }"


# rule processed:
#         input:
#                 quartet_dir_path/"{assembly_id}/tmp/trfdat"
#         output:
#                 quartet_dir_path/"{assembly_id}/tmp/trfdat/{assembly_id}.{scaffold_id}.fasta.2.7.7.80.10.500.750.OK"
#         shell:
#                 "mkdir -p {params}; "
#                 "[-f "'{input}'" ] && touch "'{output.ok}'" || {echo "'ERROR! {input} not found.'"; exit 1; }"


# rule dat_files_conversion:
#         input:
#                 input_file=quartet_dir_path/"{assembly_id}/tmp/trfdat/{assembly_id}.{scaffold_id}.fasta.2.7.7.80.10.500.750.dat"
#         output:
#                 gff_files=chromosomes_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.fasta.2.7.7.80.10.500.750.dat.gff.gz",
#                 tab_files=chromosomes_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.fasta.2.7.7.80.10.500.750.dat.tab"
#         params:
#                 options=config["chromosomes.params"],
#                 chromosomes_outdir=directory(chromosomes_dir_path/"{assembly_id}"),
#                 assembly_id=lambda wildcards: wildcards.assembly_id,
#                 scaffold_id=lambda wildcards: wildcards.scaffold_id
#         wildcard_constraints:
#                 assembly_id="|".join(config['assembly_scaffold_ids'].keys()),
#                 scaffold_id=lambda wildcards: "|".join(config['assembly_scaffold_ids'][wildcards.assembly_id])
#         log:
#                 std=log_dir_path/"chromosomes.{assembly_id}.{scaffold_id}.log",
#                 cluster_log=cluster_log_dir_path/"chromosomes.{assembly_id}.{scaffold_id}.cluster.log",
#                 cluster_err=cluster_log_dir_path/"chromosomes.{assembly_id}.{scaffold_id}.cluster.err"
#         benchmark:
#                 benchmark_dir_path/"chromosomes.{assembly_id}.{scaffold_id}.benchmark.txt"
#         resources:
#                 threads=config["chromosomes.threads"],
#                 time=config["chromosomes.time"],
#                 mem_mb=config["chromosomes.mem_mb"]
#         shell:
#                 "mkdir -p {params.chromosomes_outdir}; "
#                 "TRF.py -i {input} -o {output.gff_file} > {log.std} 2>&1; "
#                 "GFF_to_TAB.py {output.gff_file} {output.tab_file} > {log.std} 2>&1 "



# rule dat_files_conversion:
#         input:
#                 expand_result_template_list(dictionary, 'quartet_dir_path/{assembly_id}/tmp/trfdat/{assembly_id}.{scaffold_id}.fasta.2.7.7.80.10.500.750.dat')
#         output:
#                 gff_files=expand_result_template_list(dictionary, 'chromosomes_dir_path/{assembly_id}/{assembly_id}.{scaffold_id}.dat.gff.gz'),
#                 tab_files=expand_result_template_list(dictionary, 'chromosomes_dir_path/{assembly_id}/{assembly_id}.{scaffold_id}.fasta.2.7.7.80.10.500.750.dat.tab')
#         conda:
#                 "../../%s" % config["data_processing_conda_config"]
#         params:
#                 options=config["chromosomes.params"],
#                 chromosomes_outdir=directory(chromosomes_dir_path/"{assembly_id}"),
#         log:
#                 std=log_dir_path/"chromosomes.{assembly_id}.{scaffold_id}.log",
#                 cluster_log=cluster_log_dir_path/"chromosomes.{assembly_id}.{scaffold_id}.cluster.log",
#                 cluster_err=cluster_log_dir_path/"chromosomes.{assembly_id}.{scaffold_id}.cluster.err"
#         benchmark:
#                 benchmark_dir_path/"chromosomes.{assembly_id}.{scaffold_id}.benchmark.txt"
#         resources:
#                 threads=config["chromosomes.threads"],
#                 time=config["chromosomes.time"],
#                 mem_mb=config["chromosomes.mem_mb"]
#         shell:
#                 "mkdir -p {params.chromosomes_outdir}; "
#                 "TRF.py -i {input} -o {output.gff_file} > {log.std} 2>&1; "
#                 "GFF_to_TAB.py {output.gff_file} {output.tab_file} > {log.std} 2>&1 "

# rule processed:
#     input:
#         expand_dat_files_from_quartet_trfdat
#     output:
#         chromosomes_dir_path / "{assembly_id}.processed.txt"
#     shell:
#         "echo {input} 'hello' > {output}"