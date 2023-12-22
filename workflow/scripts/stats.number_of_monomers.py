#!/usr/bin/env python3

import argparse
import pandas as pd

# Функция для подсчета уникальных значений в заданном столбце DataFrame
def count_unique_values(df, column_name):
    unique_counts = df[column_name].value_counts().reset_index()  # Получение количества уникальных значений и их подсчет
    unique_counts.columns = ['monomer_size', 'monomer_count']     # Переименование столбцов
    return unique_counts

# Функция для вычисления статистик по группам в DataFrame
def calculate_stats_by_group(df, group_by_column, columns_to_calculate):
    grouped = df.groupby(group_by_column)                      # Группировка данных по заданному столбцу
    stats = grouped.agg(**columns_to_calculate).reset_index()  # Вычисление статистик для каждой группы
    return stats
    

# Функция для вычисления содержания GC-пар в последовательности
def calculate_gc_content(df, sequence_column, group_by_column):
    df[sequence_column] = df[sequence_column].astype(str)         # Преобразование столбца в строковый тип данных
    # Вычисление процента GC-пар и создание статистики по группам
    df['GC_content'] = df[sequence_column].apply(lambda x: (x.lower().count('g') + x.lower().count('c')) / len(x) * 100)
    gc_stats = calculate_stats_by_group(df, group_by_column, {'mean_GC': ('GC_content', 'mean'), 'median_GC': ('GC_content', 'median')})
    return gc_stats

# Функция для обработки файла .gff.tab
def process_tab_file(input_file, output_file):
    # Считывание данных из файла .gff.tab в DataFrame
    data = pd.read_csv(input_file, delimiter='\t')
    
    # Получение статистики по уникальным значениям в столбце 'consensus_size'
    unique_values_count = count_unique_values(data, 'consensus_size')
    #print(unique_values_count)

    # Вычисление статистик по длине последовательности в зависимости от 'consensus_size'
    sequence_lengths_stats = calculate_stats_by_group(data, 'consensus_size', {
        'mean_seq_length': ('repeat_seq', lambda x: x.str.len().mean()),
        'median_seq_length': ('repeat_seq', lambda x: x.str.len().median()),
        'min_seq_length': ('repeat_seq', lambda x: x.str.len().min()),
        'max_seq_length': ('repeat_seq', lambda x: x.str.len().max())
    })
    #print(sequence_lengths_stats)
    
    # Вычисление статистик по количеству копий в зависимости от 'consensus_size'
    copies_statistics = calculate_stats_by_group(data, 'consensus_size', {
        'mean_copies': ('copies', 'mean'),
        'median_copies': ('copies', 'median')
    })
    #print(copies_statistics)

    # Вычисление статистик содержания GC-пар для каждого 'consensus_size'
    gc_content_stats = calculate_gc_content(data, 'cons_seq', 'consensus_size')
    #print(gc_content_stats)

    # Объединение всех статистик в один DataFrame
    result = pd.merge(unique_values_count, sequence_lengths_stats, left_on='monomer_size', right_on='consensus_size')
    result = pd.merge(result, copies_statistics, on='consensus_size')
    result = pd.merge(result, gc_content_stats, on='consensus_size')
    
    # Удаление лишнего столбца 'consensus_size', оставив только 'monomer_size'
    result.drop(columns=['consensus_size'], inplace=True)

    # Переупорядочивание столбцов в DataFrame
    result = result[['monomer_size', 'monomer_count', 'mean_copies', 'median_copies', 'mean_GC', 'median_GC',
                     'mean_seq_length', 'median_seq_length', 'min_seq_length', 'max_seq_length']]

    # Запись полученных статистик в выходной файл
    result.to_csv(output_file, sep='\t', index=False, float_format='%.1f')

if __name__ == "__main__":
    # Создание парсера для аргументов командной строки
    parser = argparse.ArgumentParser(description='Obtaining statistical data from a .gff.tab file containing information for just ONE (!) scaffold.')
    parser.add_argument('-i', '--input', dest='input_file', type=str, help='Input .gff.tab file.', required=True)
    parser.add_argument('-o', '--output', dest='output_file', type=str, help='Output .tab file.', required=True)
    args = parser.parse_args()

    # Получение имен входного и выходного файлов из аргументов командной строки
    input_file = args.input_file
    output_file = args.output_file

    # Обработка входного файла для вычисления статистик и сохранения результатов в выходной файл
    process_tab_file(input_file, output_file)