#!/usr/bin/env python3

import argparse

def read_file(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            data.append(line.strip().split())
    return data

def write_bed_file(data, len_data, output_file):
    centro_dict = {}
    for entry in data:
        end_coord = entry[2]
        for len_entry in len_data:
            if end_coord == len_entry[1]:
                chrom = len_entry[0]
                if chrom not in centro_dict:
                    centro_dict[chrom] = set()
                centro_dict[chrom].add((entry[3], entry[4]))

    with open(output_file, 'w') as file:
        for chrom, coords in centro_dict.items():
            for coord_pair in coords:
                file.write(f"{chrom}\t{coord_pair[0]}\t{coord_pair[1]}\n")

def main():
    parser = argparse.ArgumentParser(description='Process centro file and length file to generate bed file')
    parser.add_argument('-q', '--centro_file', type=str, help='Path to quartet .centro.chr.txt file')
    parser.add_argument('-l', '--len_file', type=str, help='Path to .len file')
    parser.add_argument('-o', '--output_file', type=str, help='Output bed file name')

    args = parser.parse_args()
    centro_data = read_file(args.centro_file)
    len_data = read_file(args.len_file)
    output_file = args.output_file

    write_bed_file(centro_data, len_data, output_file)

if __name__ == "__main__":
    main()
