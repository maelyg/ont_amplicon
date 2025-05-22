import os
import argparse
import csv

def find_fastq_dirs(base_dir, single=False, extension=".fastq.gz"):
    fastq_entries = []
    if single:
        for root, _, files in os.walk(base_dir):
            for file in files:
                if file.endswith(extension):
                    full_path = os.path.join(root, file)
                    #rel_path = os.path.relpath(file, base_dir)
                    #subdir = rel_path.split(os.sep)[0]
                    file = file.replace('.fastq.gz', '') 
                    fastq_entries.append((file, full_path))
                    break  # No need to check more files in this directory
    else:
        for root, dirs, files in os.walk(base_dir):
            for file in files:
                if file.endswith(extension):
                    rel_path = os.path.relpath(root, base_dir)
                    subdir = rel_path.split(os.sep)[0]
                    
                    fastq_entries.append((subdir, root))
                    break  # No need to check more files in this directory

    return sorted(fastq_entries)

def main():
    parser = argparse.ArgumentParser(description="Find directories or subdirectories containing FASTQ files.")
    parser.add_argument("-d",  "--directory", help="Path to base directory to search")
    parser.add_argument("-o", "--output", default="index.csv", help="Output csv file to save the directory paths")
    parser.add_argument("-s", "--single", help="Specify whether there is only a single fastq.gz file per sample", action="store_true")

    args = parser.parse_args()

    base_dir = os.path.abspath(args.directory)
    output_path = os.path.abspath(args.output)


    if not os.path.isdir(base_dir):
        print(f"Error: '{base_dir}' is not a valid directory.")
        return

    fastq_entries = find_fastq_dirs(base_dir, args.single)

    if fastq_entries:
        with open(output_path, "w", newline="") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["sampleid","fastq_path","target_organism","target_gene","target_size","fwd_primer","rev_primer","test","method"])  # header
            for filename, dir_path in fastq_entries:
                if args.single:
                    writer.writerow([filename, dir_path, "", "", "", "", "", "", ""])
                else:
                    writer.writerow([filename, f"{dir_path}/*fastq.gz", "", "", "", "", "", "", ""])  # empty columns for other fields
        print(f"Found {len(fastq_entries)} directories. Output written to: {output_path}")
    else:
        print("No .fastq.gz files found in any subdirectories.")

if __name__ == "__main__":
    main()