#!/usr/bin/env python3

## SeA-SnaP wrapper to run pipelines and helper functions
## version: 0.9.8
## author: J.P.Pett (patrick.pett@bihealth.de)


import os, sys, time, shutil, argparse, json, yaml, re
import pandas as pd
from pathlib import Path
from tools.pipeline_tools import CovariateFileTool, SampleInfoTool

SCRIPT_DIR = Path(sys.path[0])
CONFIGS = dict(
	DE="DE_config.yaml",
	mapping="mapping_config.yaml",
	sc="sc_config.yaml",
)

CLUSTER_CONFIG = "cluster_config.json"

CLUSTER_START = dict(
	sge="export /opt/sge/lib/lx-amd64/libdrmaa.so; qsub -cwd -V -pe smp 1 -l h_vmem=4G -l h_rt=100:00:00 -P control -j y -o pipeline_log.out -e pipeline_log.err run_pipeline.sh",
	slurm="unset DRMAA_LIBRARY_PATH; unset DISPLAY; sbatch -c 1 --mem-per-cpu=4G -t 60:00:00 -p medium -o pipeline_log.out -e pipeline_log.err run_pipeline.sh",
)


############################## HELPER FUNCTIONS

def setup_working_directory(args):
	"""
	setup a working directory for running the pipeline
	"""
	choice = args.configs
	
	# paths of config files
	config_files = [SCRIPT_DIR / val for key, val in CONFIGS.items() if key in choice]

	# create working directory
	working_dir = Path(time.strftime(args.dirname))
	try:
		working_dir.mkdir(parents=True)
		print("working directory {} created...".format(str(working_dir)))
	except FileExistsError:
		print("Error: directory {} already exists!".format(str(working_dir)))
		raise

	# copy config files
	for configf in config_files:
		shutil.copy(str(configf), str(working_dir / configf.name))
		
	cl_config = SCRIPT_DIR / CLUSTER_CONFIG
	shutil.copy(str(cl_config), str(working_dir / cl_config.name))

	# symlink to wrapper
	(working_dir / "sea-snap").symlink_to(SCRIPT_DIR / "sea-snap.py")

def generate_sample_info(args):
	"""
	generate a sample info file before running the mapping pipeline
	"""
	config_files = args.config_files
	print(f"\nloading config files: {', '.join(config_files)}\n")
	
	sit = SampleInfoTool(*config_files)

	if args.add_ext:
		sit.allowed_read_extensions += args.add_ext
	
	# fill info
	if args.get_from=="parse_dir":
		sit.update_sample_info(library_default=args.library_default)
	elif args.get_from=="yaml":
		sit.read_yaml(args.input_file)
	elif args.get_from=="tsv":
		sit.read_table(args.input_file, sep=args.sep)
	elif args.get_from=="sodar":
		sit.parse_isatab(args.input_file)
		print(f"...ISA-tab parsed. Sample IDs are: {list(sit.sample_info)}")
		sit.update_sample_info(library_default=args.library_default, add=True)
	
	# write to file
	if args.write_to=="yaml":
		sit.write_yaml(args.output + ".yaml")
	elif args.write_to=="tsv":
		sit.write_table(args.output + ".tsv", sep=args.sep)
	print("\nsample info {} auto-generated. EDIT BEFORE RUNNING PIPELINE!\n".format(args.output))

def generate_covariate_file(args):
	"""
	generate a covariate file before running the DE pipeline
	"""
	steps        = [args.step] if args.step else ["star", "salmon", "feature_counts"]
	extension    = args.extension
	config_files = args.config_files

	for step in steps:
		if not extension:
			default_ext = dict(star="gene_counts.tab", salmon="sf", feature_counts="feature_counts")
			if step in default_ext:
				extension = default_ext[step]
			else:
				raise ValueError(
					f"Cannot choose a default extension for step '{step}', "
					"please set the extension on the command line."
					)
		print(f"\nSearch files with step '{step}' and extension '{extension}'...")

		cft = CovariateFileTool(*config_files)

		# fill 5 mandatory columns
		if args.tpm:
			cft.update_covariate_data(step, extension, {"tpm": ("tpm_calculator", "tsv")})
		else:
			cft.update_covariate_data(step, extension)

		# add custom columns
		if args.add_cols:
			for col in args.add_cols:
				col_name = col[0]
				col_dict = {item.split(":")[0] : "".join(item.split(":")[1:]).split(",") for item in col[1:]}
				cft.add_column(col_name, col_dict)

		# write to file
		cft.write_covariate_file(args.output)
		print("\ncovariate file {} auto-generated. EDIT BEFORE RUNNING PIPELINE!".format(args.output))

def show_matrix(args):
	"""
	print a model matrix to console
	"""
	with open(args.config_file, 'r') as stream:
		try:
			config_dict = yaml.safe_load(stream)
		except yaml.YAMLError as exc:
			print(exc)
	design = config_dict["experiment"]["design_formula"]
	print(f"design formula: {design}")
	col_names = re.findall("[^~()/:*+\s0-9\-]+", design)
	print(f"column names: {str(col_names)}")
	cov_data = pd.read_csv(args.covariate_file, sep="\t", header=0, dtype=str)
	expressions = ["-e \"{} <- c('{}')\"".format(col_name, "','".join(cov_data[col_name])) for col_name in col_names]
	cmd = "Rscript --vanilla {} -e \"model.matrix({})\"".format(" ".join(expressions), design)
	print(cmd)
	os.system(cmd)

def select_contrast(args):
	"""
	run a GUI (Shiny) to help setting contrasts for DE pipeline
	"""
	tool_path = SCRIPT_DIR / "tools" / "ContrastSelector.R"
	os.system('R --vanilla -e "source(\'{}\'); app"'.format(str(tool_path)))
	
def cleanup_cluster_log(args):
	"""
	delete log files created during cluster execution
	"""
	def rmdir(dir_name):
		dir_path = Path(dir_name)
		if dir_path.is_dir():
			for item in dir_path.iterdir():
				if item.is_dir():
					rmdir(item)
				else:
					item.unlink()
			dir_path.rmdir()

	rmdir("cluster_log")
	for p in Path(".").glob("temp_snakemake*.sh"): p.unlink()
	for p in Path(".").glob("run_pipeline.sh"): p.unlink()
	for p in Path(".").glob("core.*"): p.unlink()
	for p in Path(".").glob("*log*"): 
		if p.is_file(): 
			p.unlink()


############################## RUN PIPELINES

def run_pipeline(snakefile, args):
	"""
	run a pipeline locally or on cluster
	"""
	# local command
	if args.mode in ["local", "l"]:
		command = "snakemake --snakefile {sfile}".format(sfile = str(SCRIPT_DIR / snakefile))
		if args.snake_options:
			command += " " + " ".join(args.snake_options)
	# cluster command
	elif args.mode in ["cluster", "c"]:
		with open(CLUSTER_CONFIG, "r") as json_file:
			data = json.load(json_file)
		Path("cluster_log").mkdir(exist_ok=True)
		run_script = Path("run_pipeline.sh")
		s_command  = "#!/bin/bash\nsnakemake --snakefile {sfile}".format(sfile = str(SCRIPT_DIR / snakefile))
		if args.snake_options: s_command += " " + " ".join(args.snake_options)
		s_command += " " + data["__set_run_command__"]["snake_opt"]
		s_command += " " + "--cluster-config " + CLUSTER_CONFIG
		if args.slurm:
			s_command += " " + data["__set_run_command__"]["run_command_slurm"]
		else:
			s_command += " " + data["__set_run_command__"]["run_command_sge"]
		run_script.write_text(s_command)
		command = "set -e;" + CLUSTER_START["slurm" if args.slurm else "sge"]
		print(command)
	# run
	os.system(command)


def run_mapping_pipeline(args):
	"""
	run mapping pipeline
	"""
	run_pipeline("mapping_pipeline.snake", args)


def run_DE_pipeline(args):
	"""
	run differential expression (DE) pipeline
	"""
	run_pipeline("DE_pipeline.snake", args)


def run_sc_pipeline(args):
	"""
	run single cell pipeline
	"""
	run_pipeline("sc_pipeline.snake", args)


############################## DEFINE PARSER

parser = argparse.ArgumentParser(description="run SeA-SnaP pipelines and helpers")
subparsers = parser.add_subparsers(title="subcommands",   metavar="COMMAND", help='use COMMAND -h for more information')

### HELPERS

#--- parser for setup_working_directory
parser_working_dir = subparsers.add_parser('working_dir', help="setup a working directory for running the pipeline")
parser_working_dir.add_argument('--dirname', '-d',   default="results_%Y_%m_%d/", help="name of directory")
parser_working_dir.add_argument('--configs', '-c', nargs='+', default=["mapping", "DE", "sc"], choices=["mapping", "DE", "sc"], help="configs to be imported")
parser_working_dir.set_defaults(func=setup_working_directory)

#--- parser for generate_sample_info
parser_sample_info = subparsers.add_parser('sample_info', help="generate sample info for mapping pipeline", description=
"""Generate sample info (yaml) file for the mapping pipeline.""")
parser_sample_info.add_argument('--library_default', '-l', default="unstranded", choices=["unstranded","forward","reverse"], help="default strandedness for all samples")
parser_sample_info.add_argument('--config_files',    '-c', nargs='+', default=["mapping_config.yaml"], help="config files to be loaded")
parser_sample_info.add_argument('--add_ext',         '-a', nargs='+', default=[], help="custom extension for fastq files")
parser_sample_info.add_argument('--output',          '-o', default="sample_info", help="name of sample info file")
parser_sample_info.add_argument('--input',           '-i', default="sample_info.tsv", dest="input_file", help="import from this file; only needed if --from is used")
parser_sample_info.add_argument('--sep',             '-s', default="\t", help='separator for importing or exporting tables with --from tsv or --to tsv (default "\\t")')
parser_sample_info.add_argument('--from',            '-f', default="parse_dir", dest="get_from", choices=["parse_dir","yaml","tsv","sodar"], help="import sample info file type")
parser_sample_info.add_argument('--to',              '-t', default="yaml", dest="write_to", choices=["yaml","tsv"], help="export sample info file type")
parser_sample_info.set_defaults(func=generate_sample_info)

#--- parser for generate_covariate_file
parser_covariate_file = subparsers.add_parser('covariate_file', help="generate a covariate file for DE pipeline", description=
"""Generate a covariate file for the DE pipeline.
Five mandatory columns are automatically generated.
Additional columns can be added with --col NAME LEVELS, 
where NAME is the column name and LEVELS can be specified in two ways:
1) by group:level pairs, e.g. gr1:lvl1 gr2:lvl1 gr3:lvl2
2) by level:groups list, e.g. lvl1:gr1,gr2 lvl2:gr3""")
parser_covariate_file.add_argument('step', nargs='?', help="name of the rule to collect outputs from (e.g. 'salmon', 'star' or 'feature_counts')")
parser_covariate_file.add_argument('extension', nargs='?', help="name of the output file extension (e.g. 'sf', 'gene_counts.tab' or 'feature_counts')")
parser_covariate_file.add_argument('--config_files', nargs='+', default=["DE_config.yaml"], help="config files to be loaded")
parser_covariate_file.add_argument('--output', default="covariate_file.txt", help="name of covariate file")
parser_covariate_file.add_argument('--col', nargs='+',   action='append', dest='add_cols', help="add a column, use e.g.: --col NAME gr1:lvl1 gr2:lvl1 gr3:lvl2 ...")
parser_covariate_file.add_argument('--tpm', action='store_true', help="attach TPM column with output of TPMcalculator for display with DE pipeline")
parser_covariate_file.set_defaults(func=generate_covariate_file)

#--- parser for select_contrast
parser_select_contrast = subparsers.add_parser('select_contrast', help="display information to help choosing contrast", description=
"""Runs a Shiny App to help choosing contrasts (Open the displayed link in a browser to view).
With --console set, instead prints a model matrix based on config- and covariate file to console.""")
parser_select_contrast.add_argument('--console',        '-c',    action="store_const", const=show_matrix, dest="func", help="print model matrix to console")
parser_select_contrast.add_argument('--config_file',    '-conf', default="DE_config.yaml",     help="with --console: config file to be loaded")
parser_select_contrast.add_argument('--covariate_file', '-cov' , default="covariate_file.txt", help="with --console: name of covariate file")
parser_select_contrast.set_defaults(func=select_contrast)

### PIPELINES

#--- parser for mapping pipeline
parser_mapping = subparsers.add_parser('mapping', help="run mapping pipeline")
parser_mapping.add_argument('mode', choices=["local","l","cluster","c"], help="run locally or on cluster?")
parser_mapping.add_argument('--slurm', action='store_true', help="run using SLURM; default is SGE; only used in cluster mode")
parser_mapping.add_argument('snake_options', nargs=argparse.REMAINDER, help="pass options to snakemake (...)")
parser_mapping.set_defaults(func=run_mapping_pipeline)

#--- parser for DE pipeline
parser_DE = subparsers.add_parser('DE', help="run DE pipeline")
parser_DE.add_argument('mode', choices=["local","l","cluster","c"], help="run locally or on cluster?")
parser_DE.add_argument('--slurm', action='store_true', help="run using SLURM; default is SGE; only used in cluster mode")
parser_DE.add_argument('snake_options', nargs=argparse.REMAINDER, help="pass options to snakemake (...)")
parser_DE.set_defaults(func=run_DE_pipeline)

#--- parser for sc pipeline
parser_DE = subparsers.add_parser('sc', help="run single cell pipeline")
parser_DE.add_argument('mode', choices=["local","l","cluster","c"], help="run locally or on cluster?")
parser_DE.add_argument('--slurm', action='store_true', help="run using SLURM; default is SGE; only used in cluster mode")
parser_DE.add_argument('snake_options', nargs=argparse.REMAINDER, help="pass options to snakemake (...)")
parser_DE.set_defaults(func=run_sc_pipeline)

#--- parser for cleanup_cluster_log
parser_cleanup_log = subparsers.add_parser('cleanup_log', help="delete log files from cluster execution")
parser_cleanup_log.set_defaults(func=cleanup_cluster_log)


############################## PARSE ARGUMENTS

if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)
args = parser.parse_args()
#print(args)

############################## EXECUTE CHOSEN FUNCTION

args.func(args)



