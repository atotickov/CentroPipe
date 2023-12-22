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
                minlen=config["scaffold_filtering_by_length_params"]
        log:
                std=log_dir_path/"scaffold_filtering_by_length.{assembly_id}.log",
                cluster_log=cluster_log_dir_path/"scaffold_filtering_by_length.{assembly_id}.cluster.log",
                cluster_err=cluster_log_dir_path/"scaffold_filtering_by_length.{assembly_id}.cluster.err"
        benchmark:
                benchmark_dir_path/"scaffold_filtering_by_length.{assembly_id}.benchmark.txt"
        resources:
                threads=config["scaffold_filtering_by_length_threads"],
                time=config["scaffold_filtering_by_length_time"],
                mem_mb=config["scaffold_filtering_by_length_mem_mb"]
        shell:
                " awk '$2 >= {params.minlen} {{print $1}}' {input.lenfile} | xargs samtools faidx {input.fasta} > {output} 2> {log.std} "

# Правило для запуска квартета на сборках геномов с фильтрованными скаффолдами по минимальной длине скаффолдов
rule run_quartet:
        input:
                fasta=genomes_dir_path/"{assembly_id}"/"{assembly_id}.filtered.fasta",
                len_file=genomes_dir_path/"{assembly_id}"/"{assembly_id}.len"
        output:
                quartet_outdir=directory(quartet_dir_path/"{assembly_id}"),
                tmp_outdir=directory(quartet_dir_path/"{assembly_id}/tmp"),
                dat_files_outdir=directory(quartet_dir_path/"{assembly_id}/tmp/trfdat"),
                centromere_coordinate_txt_file=quartet_dir_path/"{assembly_id}/tmp/{assembly_id}.centro.chr.txt",
                centromere_coordinate_bed_file=quartet_dir_path/"{assembly_id}/tmp/{assembly_id}.centro.chr.bed"
        conda:
                "../../%s" % config["quartet_conda_config"]
        params:
                options=get_quartet_params(quartet_params=config["quartet_params"], global_variable=False) # результат функции - строка вида: --trf 2 7 7 80 10 500 -r 10 -n 650 -m 750
        log:
                std=log_dir_path/"run_quartet.{assembly_id}.log",
                cluster_log=cluster_log_dir_path/"run_quartet.{assembly_id}.cluster.log",
                cluster_err=cluster_log_dir_path/"run_quartet.{assembly_id}.cluster.err"
        benchmark:
                benchmark_dir_path/"run_quartet.{assembly_id}.benchmark.txt"
        resources:
                threads=config["quartet_threads"],
                time=config["quartet_time"],
                mem_mb=config["quartet_mem_mb"]
        shell:
                " MYPWD=$(pwd); "
                " cd {output.quartet_outdir}; "
                " quartet.py CentroMiner -i $MYPWD/{input.fasta} -t {resources.threads} -p {wildcards.assembly_id} {params.options} > $MYPWD/{log.std} 2>&1; "
                " cd $MYPWD/{output.dat_files_outdir}/; "
                " for datfile in *.dat; do mv $datfile ${{datfile%.*.*.*.*.*.*.*.*.*}}.{all_quartet_params}.raw.dat; done; >> $MYPWD/{log.std} 2>&1; "
                " cd $MYPWD/{output.tmp_outdir}/; "
                " $MYPWD/workflow/scripts/create_bed_file.py -q $MYPWD/{output.centromere_coordinate_txt_file} -l $MYPWD/{input.len_file} -o $MYPWD/{output.centromere_coordinate_bed_file} >> $MYPWD/{log.std} 2>&1;"

# Правило для конвертации .dat файлов в .gff и .gff.tab. Конвертируются лишь скаффолды, перечисленные в whitelist сборки
rule dat_files_conversion:
        input:
                dat_files_outdir=quartet_dir_path/"{assembly_id}/tmp/trfdat"
        output:
                gff_files=gff_tab_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.raw.dat.gff.gz",
                tab_files=gff_tab_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.raw.dat.gff.tab"
        conda:
                "../../%s" % config["data_processing_conda_config"]
        params:
                options=config["dat_files_conversion_params"],
                gff_tab_outdir=directory(gff_tab_dir_path),
        log:
                std=log_dir_path/"dat_files_conversion.{assembly_id}.{scaffold_id}.{all_quartet_params}.log",
                cluster_log=cluster_log_dir_path/"dat_files_conversion.{assembly_id}.{scaffold_id}.{all_quartet_params}.cluster.log",
                cluster_err=cluster_log_dir_path/"dat_files_conversion.{assembly_id}.{scaffold_id}.{all_quartet_params}.cluster.err"
        benchmark:
                benchmark_dir_path/"dat_files_conversion.{assembly_id}.{scaffold_id}.{all_quartet_params}.benchmark.txt"
        resources:
                threads=config["dat_files_conversion_threads"],
                time=config["dat_files_conversion_time"],
                mem_mb=config["dat_files_conversion_mem_mb"]
        shell:
                " mkdir -p {params.gff_tab_outdir}; "
                " TRF.py -i {input}/{wildcards.assembly_id}.{wildcards.scaffold_id}.{wildcards.all_quartet_params}.raw.dat -o {output.gff_files} > {log.std} 2>&1; "
                " GFF_to_TAB.py {output.gff_files} {output.tab_files} >> {log.std} 2>&1; "

rule tab_files_filtration:
        input:
                raw_tab_file=gff_tab_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.raw.dat.gff.tab"
        output:
                filtered_tab_file=gff_tab_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.filtered.dat.gff.tab"
        conda:
                "../../%s" % config["data_processing_conda_config"]
        params:
                min_monomer_copies=config["tab_files_filtration_params"][0]["MinMonomerCopies"],
                max_monomer_size=config["tab_files_filtration_params"][1]["MaxMonomerSize"],
                min_monomer_size=config["tab_files_filtration_params"][2]["MinMonomerSize"]
        log:
                std=log_dir_path/"tab_files_filtration.{assembly_id}.{scaffold_id}.{all_quartet_params}.log",
                cluster_log=cluster_log_dir_path/"tab_files_filtration.{assembly_id}.{scaffold_id}.{all_quartet_params}.cluster.log",
                cluster_err=cluster_log_dir_path/"tab_files_filtration.{assembly_id}.{scaffold_id}.{all_quartet_params}.cluster.err"
        benchmark:
                benchmark_dir_path/"tab_files_filtration.{assembly_id}.{scaffold_id}.{all_quartet_params}.benchmark.txt"
        resources:
                threads=config["tab_files_filtration_threads"],
                time=config["tab_files_filtration_time"],
                mem_mb=config["tab_files_filtration_mem_mb"]
        shell:
                " MYPWD=$(pwd); "
                " cd {gff_tab_dir_path}; "
                " awk 'NR==1 {{print}} $11>={params.min_monomer_copies} && ($12>={params.min_monomer_size} && $12<={params.max_monomer_size}) {{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$7\"\t\"$8\"\t\"$9\"\t\"$10\"\t\"$11\"\t\"$12\"\t\"$13\"\t\"$14\"\t\"$15\"\t\"$16\"\t\"$17\"\t\"$18\"\t\"$19\"\t\"$20\"\t\"$21\"\t\"$22}}' $MYPWD/{input} > $MYPWD/{output} 2> $MYPWD/{log.std}; "

rule get_stats:
        input:
                raw_tab_files=gff_tab_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.raw.dat.gff.tab",
                filtered_tab_files=gff_tab_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.filtered.dat.gff.tab"
        output:
                stats_raw_files=gff_tab_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.raw.dat.gff.stats.tab",
                stats_filtered_files=gff_tab_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.filtered.dat.gff.stats.tab"
        conda:
                "../../%s" % config["data_processing_conda_config"]
        params:
                options=config["get_stats_params"],
        log:
                std=log_dir_path/"get_stats.{assembly_id}.{scaffold_id}.{all_quartet_params}.log",
                cluster_log=cluster_log_dir_path/"get_stats.{assembly_id}.{scaffold_id}.{all_quartet_params}.cluster.log",
                cluster_err=cluster_log_dir_path/"get_stats.{assembly_id}.{scaffold_id}.{all_quartet_params}.cluster.err"
        benchmark:
                benchmark_dir_path/"get_stats.{assembly_id}.{scaffold_id}.{all_quartet_params}.benchmark.txt"
        resources:
                threads=config["get_stats_threads"],
                time=config["get_stats_time"],
                mem_mb=config["get_stats_mem_mb"]
        shell:
                " MYPWD=$(pwd); "
                " cd {gff_tab_dir_path}; "
                " $MYPWD/workflow/scripts/stats.number_of_monomers.py -i $MYPWD/{input.raw_tab_files} -o $MYPWD/{output.stats_raw_files} > $MYPWD/{log.std} 2>&1; "
                " $MYPWD/workflow/scripts/stats.number_of_monomers.py -i $MYPWD/{input.filtered_tab_files} -o $MYPWD/{output.stats_filtered_files} >> $MYPWD/{log.std} 2>&1; "

rule get_fasta:
        input:
                filtered_tab_files=gff_tab_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.filtered.dat.gff.tab"
        output:
                filtered_fasta_files=fasta_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.filtered.fasta"
        conda:
                "../../%s" % config["data_processing_conda_config"]
        params:
                options=config["get_fasta_params"],
        log:
                std=log_dir_path/"get_fasta.{assembly_id}.{scaffold_id}.{all_quartet_params}.log",
                cluster_log=cluster_log_dir_path/"get_fasta.{assembly_id}.{scaffold_id}.{all_quartet_params}.cluster.log",
                cluster_err=cluster_log_dir_path/"get_fasta.{assembly_id}.{scaffold_id}.{all_quartet_params}.cluster.err"
        benchmark:
                benchmark_dir_path/"get_fasta.{assembly_id}.{scaffold_id}.{all_quartet_params}.benchmark.txt"
        resources:
                threads=config["get_fasta_threads"],
                time=config["get_fasta_time"],
                mem_mb=config["get_fasta_mem_mb"]
        shell:
                " MYPWD=$(pwd); "
                " cd {gff_tab_dir_path}; "
                " awk -F'\t' 'NR>1 {{print \">\"$1\".\"$9, \"period: \"$4\"-\"$5\" | copies: \"$11\" | consensus size: \"$12\"\\n\"$21}}' $MYPWD/{input} > $MYPWD/{output} 2> $MYPWD/{log.std}; "

rule rotation:
        input:
                filtered_fasta_files=fasta_dir_path/"{assembly_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.filtered.fasta"
        output:
                rotated_fasta_files=fasta_dir_path/"{assembly_id}/{scaffold_id}/{assembly_id}.{scaffold_id}.{all_quartet_params}.rotated.fasta"
        conda:
                "../../%s" % config["general_conda_config"]
        params:
                max_difference=config["rotation_params"][0]["MaxDifference"],
                e_value=config["rotation_params"][1]["EValue"]
        log:
                std=log_dir_path/"rotation.{assembly_id}.{scaffold_id}.{all_quartet_params}.log",
                cluster_log=cluster_log_dir_path/"rotation.{assembly_id}.{scaffold_id}.{all_quartet_params}.cluster.log",
                cluster_err=cluster_log_dir_path/"rotation.{assembly_id}.{scaffold_id}.{all_quartet_params}.cluster.err"
        benchmark:
                benchmark_dir_path/"rotation.{assembly_id}.{scaffold_id}.{all_quartet_params}.benchmark.txt"
        resources:
                threads=config["rotation_threads"],
                time=config["rotation_time"],
                mem_mb=config["rotation_mem_mb"]
        shell:
                " MYPWD=$(pwd); "
                " cd {fasta_dir_path}/{wildcards.assembly_id}/{wildcards.scaffold_id}; "
                " rotate_sequences.py -i $MYPWD/{input} -o $MYPWD/{fasta_dir_path}/{wildcards.assembly_id}/{wildcards.scaffold_id}/$(basename {output} .rotated.fasta) -m {params.max_difference} -e {params.e_value} > $MYPWD/{log.std} 2>&1; "
