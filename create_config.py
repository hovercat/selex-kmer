#!/usr/bin/env python
import os


######################################################
# Argument parsing functions
######################################################
def parse_string(s, default=None):
    s = s.strip()
    if s is not None and s != "":
        return s
    else:
        return default


def parse_list(s, delimiter=' ', default=None):
    s = s.strip()
    return s.split(delimiter)


def parse_boolean(s, default=None):
    s = s.strip()
    if s.isnumeric():
        return int(s) == 1
    else:
        return default


def parse_float(s, default=None):
    s = s.strip()
    if s.isnumeric():
        return float(s)
    else:
        return default


def parse_int(s, default=None):
    s = s.strip()
    if s.isnumeric():
        return int(s)
    else:
        return default


def make_sane_string(s):
    return ''.join([char for char in s.strip() if char.isalpha() or char.isdigit() or char == ' ']).rstrip()


######################################################
# Starting config creation
######################################################
arguments = dict()

print("Creating a config file for selex-kmer workflow.")
print(
    "Please answer the following questions. If you don't know an answer, in some cases you can press enter without an answer to opt for the default value.")

# Experiment Name
while not arguments.get("experiment"):
    arguments["experiment"] = parse_string(input("What's your SELEX experiment's unique identifying name?\n"))
    if not arguments["experiment"]:
        print("There is no default value for the experiment's name.")

# Input Directory
arguments["input_dir"] = parse_string(input("In which directory are the preprocessed FASTA-files stored? (default: "
                                            "./input_dir)\n"), "./input_dir")

# Output Directory
arguments["output_dir"] = parse_string(
    input("In which directory to write output files to? (default: './output_{}/)\n".format(make_sane_string(arguments["experiment"]))),
    './output_{}'.format(make_sane_string(arguments["experiment"]))
)

while not arguments.get("derep_csv"):
    arguments["derep_csv"] = parse_string(input("In which directory is the dereplicated csv file stored as provided by selex-assess?\n"))
while not arguments.get("derep_fasta"):
    arguments["derep_fasta"] = parse_string(input("In which directory is the dereplicated fasta file stored as provided by selex-assess?\n"))

# Check input directory for FASTA-Files
if os.path.exists(arguments["input_dir"]):
    fasta_files = [f for f in os.listdir(arguments["input_dir"]) if
                   os.path.isfile(os.path.join(arguments["input_dir"], f))
                   and (f.lower().endswith(".fasta") or f.lower().endswith(".fa"))]
    fasta_files.sort()

    print("Files found in the input directory '{}'".format(arguments["input_dir"]))
    print('\n'.join(fasta_files))
else:
    print("Input directory '{}' does not exist.".format(arguments["input_dir"]))

# Ask for round names
while not arguments.get("round_order"):
    arguments["round_order"] = parse_list(
        input(
            "Please provide the list of SELEX rounds in sequential order, separated by a space. (e.g. R0 R2 R4 R6)\n"),
        delimiter=" ")
    if not arguments.get("round_order"):
        print("Please specify the round order of your FASTA-files!")

# Ask for FASTA-Pattern
while not arguments.get("fasta_pattern"):
    arguments["fasta_pattern"] = parse_string(input(
        "Please provide a FASTA-search pattern, including a wild-card. For example: \nIf you had R0.fasta your search pattern would be: \'*1.fasta\'\n"))
    if not arguments["fasta_pattern"]:
        print("Please specify a FASTA-search pattern!")

while not arguments.get("k"):
    arguments["k"] = parse_int(input("What's your kmer size? default: k=6.\n"), 6)
while not arguments.get("cpus"):
    arguments["cpus"] = parse_int(input("What's your max number of available CPUs?\n"))

# Finished reading config parameters.
print("")
config_file_name = "./{}.config".format(make_sane_string(arguments["experiment"]))
with open(config_file_name, "w") as config_file:
    config_file.write("params {\n")
    for key, value in arguments.items():
        config_file.write("\t")
        if isinstance(value, list):
            config_file.write('{key} = ["'.format(key=key))
            config_file.write('","'.join(value))
            config_file.write('"]\n')
        elif isinstance(value, int):
            config_file.write('{key} = {value}\n'.format(key=key, value=value))
        elif isinstance(value, str):
            config_file.write('{key} = "{value}"\n'.format(key=key, value=value))
        elif isinstance(value, bool):
            if value:
                config_file.write('{key} = true\n')
            else:
                config_file.write('{key} = false\n')
    config_file.write("}\n")

print("Done.")
print("Your config file was saved to '{}'.".format(config_file_name))
