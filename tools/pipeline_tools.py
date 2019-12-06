## SeA-SnaP tools used in pipelines
## version: 0.9.5
## author: J.P.Pett (patrick.pett@bihealth.de)

import sys, os, re, shutil, hashlib, itertools, yaml, pandas as pd
from collections import namedtuple, Mapping, OrderedDict
from copy import deepcopy
from time import strftime
from warnings import warn
from pathlib import Path
from glob import iglob, glob
#from snakemake.io import glob_wildcards

yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_dict(dict(data)))

##################################################################################################################################
#---------------------------------------------- base class for handling file paths ----------------------------------------------#
##################################################################################################################################

class PipelinePathHandler:
	"""
	Generates file paths from patterns defined in config.
	
	Loads and checks paths from config file upon initialization and provides
	methods to fill and expand the wildcards.
	"""
	
	allowed_wildcards          = ["step", "extension"]
	required_wildcards_out_log = ["step", "extension"]
	required_wildcards_in      = []
	
	def __init__(self, workflow, test_allowed_wildcards=True):
		# load test_config for validity check of config values
		#if type(test_config) is dict:
		#	self.test_config = test_config
		#elif isinstance(test_config, str):
		#	with open(test_config, "r") as stream:
		#		try:
		#			self.test_config = yaml.safe_load(stream)
		#		except yaml.YAMLError as exc:
		#			print(exc)
		#else:
		#	raise TypeError("Wrong type of argument test_config: must be dict or str.")
	
		self.snakemake_workflow = workflow
		
		self.out_path_pattern = self.snakemake_workflow.config["pipeline_param"]["out_path_pattern"]
		self.log_path_pattern = self.snakemake_workflow.config["pipeline_param"]["log_path_pattern"]
		self.in_path_pattern  = self.snakemake_workflow.config["pipeline_param"]["in_path_pattern"]
		
		self.wildcard_constraints = self._prepare_inpathpattern()
		
		self.out_dir_pattern = "/".join(self.out_path_pattern.split("/")[:-1])
		
		self.out_path_wildcards    = self._get_wildcard_list(self.out_path_pattern)
		self.log_path_wildcards    = self._get_wildcard_list(self.log_path_pattern)
		self.in_path_wildcards     = self._get_wildcard_list(self.in_path_pattern)
		self.outdir_path_wildcards = self._get_wildcard_list(self.out_dir_pattern)
		
		self.input_choice = self._set_input_choice(self.snakemake_workflow.config)
		
		# check consistency of loaded input and output patterns
		self._test_config_input(test_allowed_wildcards)
		
		# wildcard values
		self.wildcard_values = self._get_wildcard_values_from_input(self.in_path_pattern)
		
		# neutral replacement for optional wildcards
		self.opt_wildcard_placeholders = {w: "{{{}}}".format(w) for w in set(self.in_path_wildcards)-set(self.required_wildcards_in)}
		
		# path pattern of out directory may miss some wildcards
		placeholders_in_outdir_path  = set(self.outdir_path_wildcards) & set(self.opt_wildcard_placeholders.keys())
		self.opt_wildc_placeh_outdir = {key: self.opt_wildcard_placeholders[key] for key in placeholders_in_outdir_path}
	
	
	#---------------------------------------------------- helper methods ----------------------------------------------------#
	
	def _test_config_input(self, test_allowed_wildcards):
		""" test whether wildcards in path patterns are valid """
		if test_allowed_wildcards:
			# test if all wildcards allowed
			if not all(x in self.allowed_wildcards for x in self.out_path_wildcards+self.log_path_wildcards+self.in_path_wildcards):
				raise ValueError("Error in config file: unknown wildcards. allowed wildcards: {}".format(self.allowed_wildcards))
		
		# test for required wildcards in all patterns
		if not all(x in self.out_path_wildcards and x in self.log_path_wildcards for x in self.required_wildcards_out_log):
			raise ValueError("Error in config file: 'step', 'extension', 'sample' and 'name' must be wildcards "
					"in out_path_pattern and log_path_pattern")
		if not all(x in self.in_path_wildcards for x in self.required_wildcards_in):
			raise ValueError("Error in config file: 'step', 'extension', 'sample' and 'name' must be wildcards "
					"in out_path_pattern and log_path_pattern")

	def _test_config_general(self, base_dict, check_dict):
		""" test whether values set in the config are valid """
		for key, val in check_dict.items():
		
			if key not in base_dict and key not in ["__num__", "__opt__", "__any__", "__any_other__"]:
				raise IndexError("Error in config file: The required key {} was not defined!".format(key))
			
			# if number
			if key == "__num__":
				if not isinstance(base_dict[key], (int, float)):
					raise ValueError("Error in config file: value of {} must be a number!".format(key))
				assert len(val)==2
				if val[0] and not base_dict[key] > val[0]: 
					raise ValueError("Error in config file: value of {} must be >{}!".format(key, val[0]))
				if val[1] and not base_dict[key] < val[1]: 
					raise ValueError("Error in config file: value of {} must be <{}!".format(key, val[1]))
			
			# if string
			if isinstance(val, str):
				if isinstance(base_dict[key], str):
					if not re.fullmatch(base_dict[key], val):
						raise ValueError("Error in config file: value of {} does not match {}!".format(key, val))
				else:
					raise ValueError("Error in config file: value of {} should be a string! got: {}".format(key, val))
			# if dict
			elif isinstance(val, dict):
				pass
			
			if key in base_dict and isinstance(base_dict[key], dict) and isinstance(check_dict[key], Mapping):
				base_dict[key] = dict_merge(base_dict[key], check_dict[key])
			else:
				base_dict[key] = check_dict[key]
		return base_dict
		
	def _set_input_choice(self, config):
		""" create input choice dictionary from config ensuring a standardized structure """
		input_choice = {}
		if "input_choice" in config["pipeline_param"] and isinstance(config["pipeline_param"]["input_choice"], dict): 
			for key, val in config["pipeline_param"]["input_choice"].items():
				if val:
					if   isinstance(val, str):
						input_choice[key] = [val]
					elif isinstance(val, list):
						input_choice[key] = val
					else:
						raise TypeError("Error in config: input_choice elements must have dict values "
								"str or list, not {}!".format(type(val)))
		return input_choice
		
	def _get_wildcard_list(self, pattern):
		return re.findall("{([^}]+)}", pattern)
		
	def _make_name(self, description):
		return re.sub("\.|/|-|\s", "_", description)
		
	@staticmethod
	def _md5(fname):
		hash_md5 = hashlib.md5()
		with open(fname, "rb") as f:
			for chunk in iter(lambda: f.read(4096), b""):
				hash_md5.update(chunk)
		return hash_md5.hexdigest()
		
	def _get_wildcard_combinations(self, wildcard_values):
		""" go through wildcard values and get combinations """
		
		combinations = []
		WildcardComb = namedtuple("WildcardComb", [s for s in wildcard_values])
		wildcard_comb_num = len(wildcard_values["sample"])
		assert all([len(val)==wildcard_comb_num for val in wildcard_values.values()])
		for index in range(wildcard_comb_num):
			combinations.append(WildcardComb(**{key: wildcard_values[key][index] for key in wildcard_values}))
		return combinations

	def _prepare_inpathpattern(self):
		""" read and remove wildcard constraints from in_path_pattern """
		wildcards = re.findall("{([^{}]+)}", self.in_path_pattern)
		wildcard_constraints = {}
		for wildcard in wildcards:
			comp = wildcard.split(",")
			if len(comp)>1:
				wildcard_constraints[comp[0]] = comp[1]
				self.in_path_pattern = self.in_path_pattern.replace(wildcard, comp[0])
		return wildcard_constraints
		
	def _get_wildcard_values_from_input(self, input_pattern, unix_style=True):
		""" go through files in input path and get values matching the wildcards """
		
		glob_pattern =       re.sub("{[^}./]+}",   "*", input_pattern)
		wildcards    =   re.findall("{([^}./]+)}",      input_pattern)
		input_files  = iglob(glob_pattern + ("*" if glob_pattern[-1]!="*" else ""), recursive=True)
			
		wildcard_values = {w:[] for w in wildcards}
		for inp in input_files:
			self._get_wildcard_values_from_file_path(inp, input_pattern, wildc_val=wildcard_values, unix_style=unix_style)
		return wildcard_values
		
	def _wildc_replace(self, matchobj):
		""" method used with re.sub to generate match pattern from path pattern """
		wildc_name = matchobj.group(1)
		if wildc_name in self.wildcard_constraints:
			return "({})".format(self.wildcard_constraints[wildc_name].replace("//","/"))
		elif wildc_name == "extension":
			return "([^}/]+)"
		else:
			return "([^}./]+)"
		
	def _get_wildcard_values_from_file_path(self, filename, input_pattern, wildc_val={}, unix_style=True):
		""" get values matching wildcards from given file path """
				
		match_pattern   =       re.sub("\\\\{([^}./]+)\\\\}", self._wildc_replace, re.escape(input_pattern))
		wildcards       =   re.findall("{([^}./]+)}",                                        input_pattern)
		if unix_style:
			match_pattern = re.sub(r"\\\*\\\*", "[^{}]*",   match_pattern)
			match_pattern = re.sub(r"(?<!\[\^{}\]\*)\\\*",      "[^{}./]*", match_pattern)
			
		wildcard_values = wildc_val if wildc_val else {w:[] for w in wildcards}
		
		matches = re.match(match_pattern, filename).groups()
		assert len(matches)==len(wildcards)
		seen = set()
		for index,wildc in enumerate(wildcards):
			if not wildc in seen:
				wildcard_values[wildc].append(matches[index])
				seen.add(wildc)
				
		return wildcard_values
		
	def _choose_input(self, wildcards, choice_name, options):
		""" called by wrapper choose_input() """
		if wildcards and hasattr(wildcards, choice_name):
			input_from = getattr(wildcards, choice_name)
		elif choice_name in self.input_choice:
			inp_choice = self.input_choice[choice_name]
			if   isinstance(inp_choice, str):
				input_from = inp_choice
			elif isinstance(inp_choice, list):
				input_from = inp_choice[0]
			else:
				raise TypeError("Error choosing input: {} is not a valid type for input choice! "
						"(input_choice in config for {})".format(type(inp_choice), choice_name))
		else:
			raise KeyError("Error choosing input: no wildcard '{}' passed and no input choice for '{}' "
					"specified in config!".format(choice_name, choice_name))
					
		for inp in options:
			if input_from == inp["step"]:
				return self.file_path(**inp)
		raise ValueError("Error choosing input: no valid mapping input type specified! (got: {})".format(input_from))
		
	#-------------------------------------------- methods used in snakemake file --------------------------------------------#
		
	def choose_input(self, choice_name, options):
		""" 
		choose rule to import from based on wildcard or config
		
		One option is that a wildcard with name choice_name is passed, giving the name of a rule that should be imported from.
		This is useful, when the pipeline should compute results for both inputs in one run. Then the wildcard in the ouput file
		paths distinguishes these runs.
		The other option is that the choice of input is specified in the config file (under pipeline_params.input_choice).
		This is useful, if only one of the possible inputs should be used in one run.
		All parameters for both options should be provided and the behaviour of the pipeline will depend on, whether a corresponding
		wildcard was provided in the output path pattern, computing results for all or just one input.
		
		In the config file under pipeline_params.input_choice for each choice_name a list of input rules can be specified. In the wildcard
		mode (when a respective wildcard was provided in the output path pattern) expansion will happen over the whole list to compute
		results. In the config mode (when the wildcard was not included) only the first element of the list will be used.
		
		:param options:  wildcard arguments to generate different inputs (list of dicts)
		:param choice_name:  name of the wildcard (string) used for the decision
		:param wildcards:  wildcards object from snakemake
		:returns: a function that is used by snakemake to obtain the input file path
		"""
		return lambda wildcards: self._choose_input(wildcards, choice_name, options)
		
	def wildcard_values_from(self, filepath, in_path_pattern=True):
		"""
		Parse a file path and return a dict of wildcard names to values.
		
		:param filepath: a file path (string)
		:param in_path_pattern: use in_path_pattern to extract wildcards, if False use out_path_pattern
		:returns: a dict of wildcard names to values
		"""
		if in_path_pattern:
			return self._get_wildcard_values_from_file_path(filepath, self.in_path_pattern)
		else:
			return self._get_wildcard_values_from_file_path(filepath, self.out_path_pattern)
		
	def log(self, out_log, script, step, extension, **kwargs):
		"""
		create various log files and return a path to the script file for execution
		
		:param out_log: the main log file of rule 'step' that will contain output and error messages
		:param script: the script with wildcards filled by snakemake
		:param step: the rule for which logs shall be generated
		:param extension: the extension of the script file, e.g. '.R' or '.sh'
		:returns: path to a file containing the script
		"""
		script_file = self.file_path(step, extension,        log=True, **kwargs)
		config_yaml = self.file_path(step, "config.yaml",    log=True, **kwargs)
		conda_list  = self.file_path(step, "conda_list.txt", log=True, **kwargs)
		conda_info  = self.file_path(step, "conda_info.txt", log=True, **kwargs)
		conda_env   = self.file_path(step, "conda_env.yaml", log=True, **kwargs)
		
		# write script to file
		with open(script_file, "w") as f: f.write(script)
		
		# dump merged config
		with open(config_yaml, "w") as f: yaml.dump(self.snakemake_workflow.config, f, default_flow_style=False)
		
		# write conda logs
		os.system("conda list > {}".format(conda_list))
		os.system("conda info > {}".format(conda_info))
		os.system("conda env export > {}".format(conda_env))
		
		# git and snakefile info to out_log
		git_dir     = self.snakemake_workflow.basedir
		out_log_abs = str(Path(out_log).resolve())
		os.system('echo "--------------- git info ---------------" > {}'.format(out_log))
		os.system('cd {}; echo "name:" $(git rev-parse --show-toplevel) >> {}'.format(git_dir, out_log_abs))
		os.system('cd {}; echo "branch:" $(git symbolic-ref HEAD) >> {}'.format(git_dir, out_log_abs))
		os.system('cd {}; echo "hash:" $(git rev-parse HEAD) >> {}'.format(git_dir, out_log_abs))
		os.system('cd {}; {{ echo "status:"; git status -sb; }} >> {}'.format(git_dir, out_log_abs))
		os.system('echo "--------------- snakefile --------------" >> {}'.format(out_log))
		os.system('echo "snakefile name: {}" >> {}'.format(self.snakemake_workflow.snakefile, out_log))
		os.system('echo "snakefile md5 : {}" >> {}'.format(self._md5(self.snakemake_workflow.snakefile), out_log))
		os.system('echo "----------------------------------------" >> {l}; echo >> {l}'.format(l=out_log))
		
		return script_file

	def export(self):
		"""
		export selected results into a separate folder structure (as configured in config file)
		"""
		export_spec = self.snakemake_workflow.config["export"]
		for pattern in export_spec["path_pattern"]:
			pat = strftime(pattern).replace("{GENOME}", self.snakemake_workflow.config["organism"]["genome_version"])
			wildcards = self._get_wildcard_list(pat)
			key_wcs   = [wc for wc in wildcards if wc[:6] == "files:"]
			# no nested keys for now:
			assert len(key_wcs)==1
			key = key_wcs[0].split(":")[1:]
			assert len(key)>0
			pat = pat.replace("{{{}}}".format(key_wcs[0]), key[0])
			for arg_dct in export_spec["_".join(key)]:
				mode = "file"
				if "dir" in arg_dct:
					mode = "dir"
					arg_dct = arg_dct["dir"]
				extra_wcs = set(wildcards) - set(key_wcs) - (set(arg_dct) & set(wildcards))
				# only either 'sample' (mapping) or 'contrast' (DE) for now:
				assert len(extra_wcs)<=1
				if extra_wcs:
					extra_wc   = list(extra_wcs)[0]
					wc_in_dct  = {k:v for k,v in arg_dct.items() if k in wildcards}
					search_pat = self.out_path_pattern if not "log" in arg_dct else self.log_path_pattern
					
					if mode == "file":
						source = self.expand_path(**arg_dct)
					else:
						source = self.file_path(**{**arg_dct, "extension": "{extension}"})
						source = str(Path(source).parent / "**")
						source = glob(source.replace("{{{}}}".format(extra_wc), "*"), recursive=True)
						search_pat = str(Path(search_pat).parent)
						
					get_wc = self._get_wildcard_values_from_file_path
					target = [pat.format(**{**wc_in_dct, extra_wc: get_wc(src, search_pat)[extra_wc][0]}) for src in source]
				else:
					source = [self.file_path(**arg_dct)]
					target = [pat.format(**{k:v for k,v in arg_dct.items() if k in wildcards})]
				assert len(source)==len(target)
				for i in range(len(source)):
					if mode == "dir": target[i] = str(Path(target[i]) / Path(source[i]).name)
					print("copy {} to {}...".format(source[i], target[i]))
					Path(target[i]).parent.mkdir(exist_ok = True, parents = True)
					shutil.copy2(source[i], target[i])
					Path(target[i] + ".md5").write_text(self._md5(target[i]))
		
		
##################################################################################################################################
#----------------------------------------------- child class for mapping pipline ------------------------------------------------#
##################################################################################################################################

class MappingPipelinePathHandler(PipelinePathHandler):
	""" path handler for mapping pipeline """

	allowed_wildcards          = ["step", "extension", "sample", "mate", "batch", "flowcell", "lane", "library"]
	required_wildcards_out_log = ["step", "extension", "sample"]
	required_wildcards_in      = ["sample"]
	
	def __init__(self, workflow, **kwargs):
		super().__init__(workflow, **kwargs)
		
		self.samples    = self.snakemake_workflow.config["sample_info"]
		self.sample_ids = list(self.samples.keys())
		
		# wildcard value-combinations parsed from input directory
		self.wildcard_combs  = self._get_wildcard_combinations_per_sample(self.wildcard_values)
		
		# paths to static files
		self.data_paths = self.snakemake_workflow.config["organism"]
		
		# test if all sample IDs can be found in input path
		if not set(self.sample_ids).issubset(set(self.wildcard_values["sample"])):
			raise ValueError("Error in config file: not all specified samples could be found using given input path pattern")
		
		# write log with parsed values to file for debugging
		self._write_log(sample="all")
		
		
	#---------------------------------------------------- helper methods ----------------------------------------------------#
	
	def _test_config_input(self, test_allowed_wildcards):
		super()._test_config_input(test_allowed_wildcards)
		
		# test whether wildcards used consistently across patterns
		if not set(self.out_path_wildcards) == set(self.log_path_wildcards) == set(self.in_path_wildcards) | set(self.required_wildcards_out_log) - set(self.required_wildcards_in):
			raise ValueError("Error in config file: out_path_pattern, log_path_pattern and in_path_pattern "
			"do not contain the same wildcards. out: {}, log: {}, in: {}".format(set(self.out_path_wildcards), set(self.log_path_wildcards), set(self.in_path_wildcards)))
			
		# test if all wildcards used in outdir path pattern
		if not set(self.out_path_wildcards) == set(self.outdir_path_wildcards) | set(["extension"]):
			raise ValueError("Error in config file: all wildcards of out and log dir should be used in folder names (exception: {{extension}}), otherwise "
			"different rules might compute output in the same folder, which can lead to mixed or deleted intermediate files. " 
			"in folder name: {}, all: {}".format(set(self.outdir_path_wildcards), set(self.out_path_wildcards)))
	
	def _get_wildcard_combinations_per_sample(self, wildcard_values):
		""" go through wildcard values and get combinations per sample """
		
		wildcard_comb_num = len(wildcard_values["sample"])
		assert all([len(val)==wildcard_comb_num for val in wildcard_values.values()])
		WildcardComb = namedtuple("WildcardComb", [s for s in wildcard_values])
		per_sample_comb = {s:[] for s in wildcard_values["sample"]}
		for index in range(wildcard_comb_num):
			sample = wildcard_values["sample"][index]
			per_sample_comb[sample].append(WildcardComb(**{key: wildcard_values[key][index] for key in wildcard_values}))
		return per_sample_comb
		
	def _write_log(self, **kwargs):
		filename = self.file_path(step="MappingPipelinePathHandler", extension="log", log=True, batch="all_batches", 
		                          flowcell="all_flowcells", lane="all_lanes", library="all_libraries", **kwargs)
		os.makedirs(os.path.dirname(filename), exist_ok=True)
		with open(filename, "w") as f:
			f.write("pattern input:\n   in path pattern: {}\n   out path pattern: {}\n   log path pattern: {}\n\n"
			"parsed wildcards:\n   in path wildcards: {}\n   out/log path wildcards: {}\n   out dir wildcards: {}\n\n"
			"samples:\n   sample IDs: {}\n   sample info: {}\n\n"
			"wildcard values:\n   per wildcard: {}\n   combinations: {}\n\n"
			"neutral replacement of optional wildcards:\n   out/log: {}\n   out dir: {}".format(self.in_path_pattern, self.out_path_pattern, self.log_path_pattern, 
			self.in_path_wildcards, self.out_path_wildcards, self.outdir_path_wildcards, self.sample_ids, self.samples, self.wildcard_values, self.wildcard_combs, 
			self.opt_wildcard_placeholders, self.opt_wildc_placeh_outdir))
			
			
	#-------------------------------------------- methods used in snakemake file --------------------------------------------#
	
	def get_fastq_pairs(self, wildcards, mate=0, mate_key=""):
		"""
		Generate paths to one or more .fastq input files for a given name extension (e.g. paired end extension).
		
		:param wildcards: rule wildcards
		:param mate: index of extension type to return (based on mate_key)
		:param mate_key: key to config dict where different extensions to a sample ID are defined (e.g. paired end extensions)
		:returns: a list with paths to specified input files, [] if extension at given index (mate) does not exist
		"""
		kwargs_out = {key: getattr(wildcards, key, val) for key, val in self.opt_wildcard_placeholders.items()}
		pattern_list = []
		seen = set()
		for comb in self.wildcard_combs[wildcards.sample]:
			#TODO: case of ignored wildcards? meant for e.g. allFlowcell
			kwargs_filled = {key: getattr(comb, key) if ("{" in val or val not in self.wildcard_values[key]) and key != "mate" else val for key,val in kwargs_out.items()}
			kwargs_id_tup = tuple(kwargs_filled[key] for key in sorted(kwargs_filled))
			if kwargs_id_tup not in seen:
				seen.add(kwargs_id_tup)
				if mate_key:
					mate_lst = self.samples[wildcards.sample][mate_key]
					if len(mate_lst)>mate: kwargs_filled["mate"] = mate_lst[mate]
					pattern  = self.in_path_pattern.format(sample = wildcards.sample, **kwargs_filled)
					pattern += self.samples[wildcards.sample]["read_extension"]
				else:
					if "mate" not in kwargs_filled: kwargs_filled["mate"] = "*"
					pattern = self.in_path_pattern.format(sample=wildcards.sample, **kwargs_filled) + self.samples[wildcards.sample]["read_extension"]
				pattern_list.append(pattern)
		return [path for pat in pattern_list for path in iglob(pat, recursive=True)]
		
	def file_path(self, step, extension, sample="{sample}", log=False, **kwargs):
		"""
		Generate single path for intermediate and output or log files.
		
		:param step:  Snakemake rule for which the paths are generated
		:param extension: file extension for the generated file path
		:param sample: sample ID to be included in the file path
		:param log: generate path to logfile if log==True otherwise generate path to output/intermediate file
		:param **kwargs: if used specify replacement for {batch}, {flowcell}, {lane}, etc. ...
		"""
		path_pattern = self.log_path_pattern if log else self.out_path_pattern
		kwargs_out = {key: kwargs[key] if key in kwargs else val for key, val in self.opt_wildcard_placeholders.items()}
		return path_pattern.format(step=step, extension=extension, sample=sample, **kwargs_out)
		
	def out_dir_name(self, step, sample="{sample}", **kwargs):
		"""
		Generate single path to intermediate and output file directory.
		
		:param step:  Snakemake rule for which the paths are generated
		:param sample: sample ID to be included in the file path
		"""
		kwargs_out = {key: kwargs[key] if key in kwargs else val for key, val in self.opt_wildc_placeh_outdir.items()}
		return self.out_dir_pattern.format(step=step, sample=sample, **kwargs_out)
		
	def expand_path(self, step, extension="", **kwargs):
		"""
		Generate multiple paths for intermediate and output files corresponding to different sample IDs (and e.g. paired end extensions).

		Always generates paths over all sample IDs and their combinations with keyword arguments (kwargs) that were not
		defined like, e.g. flowcell="all_flowcells", where "all_flowcells" is an arbitrary string that will be filled in the generated paths
		instead of the real flowcell IDs. Only keyword arguments that correspond to wildcards in the config path pattern (with the same name) are used.
		
		:param step:  Snakemake rule name (string) for which the paths are generated
		:param extension: file extension for the generated file path; '' to generate paths to directories
		:param **kwargs: replacement strings for optional wildcards, e.g. batch, flowcell, lane (see description)
		:returns: list of paths
		"""
		kwargs_out = {key: kwargs[key] if key in kwargs else val for key, val in self.opt_wildcard_placeholders.items()}
		paths = []
		for sample in self.sample_ids:
			seen = set()
			for comb in self.wildcard_combs[sample]:
				kwargs_filled = {key: getattr(comb, key) if "{" in val else val for key,val in kwargs_out.items()}
				kwargs_id_tup = tuple(kwargs_filled[key] for key in sorted(kwargs_filled))
				if kwargs_id_tup not in seen:
					seen.add(kwargs_id_tup)
					if extension:
						out_path_pattern = self.file_path(step, extension, sample=sample, **kwargs_filled)
					else:
						out_path_pattern = self.out_dir_name(step, sample=sample, **kwargs_filled)
					paths.append(out_path_pattern)
		return paths
	
	def link_index(self, step, sample="{sample}", entry="nopath", **kwargs):
		"""
		Generate symbolic link to index folder (if provided in config).
		
		A key-value pair is assumed to be in the organism config file with the name of the step (e.g. "star_index") as key 
		and a path to a corresponding index/input file. If no index should be linked, keep the key, but leave the value empty.
		
		:param step:  Snakemake rule for which the paths are generated
		:param sample: sample ID to be included in the file path
		:param entry: set the path to be linked explicitly
		:param **kwargs: if used specify replacement for {batch}, {flowcell}, {lane}, etc. ...
		"""
		loc = Path(self.out_dir_name(step, sample, **kwargs))
		if entry != "nopath":
			index = entry
		elif step in self.data_paths:
			index = self.data_paths[step]
		else:
			raise ValueError("Error linking index: no keyword provided in config to set index for {}!".format(step))
		if index:
			ind = Path(index)
			loc.parent.mkdir(exist_ok=True, parents=True)
			if not loc.is_dir():
				# create
				loc.symlink_to(ind.resolve(), target_is_directory=True)
			elif not loc.samefile(index):
				# update
				loc.unlink()
				loc.symlink_to(ind.resolve(), target_is_directory=True)


##################################################################################################################################
#------------------------------------------------- child class for DE pipline ---------------------------------------------------#
##################################################################################################################################

class DEPipelinePathHandler(PipelinePathHandler):
	""" path handler for differential expression pipeline """

	allowed_wildcards          = ["step", "extension", "sample", "mate", "batch", "flowcell", "lane", "contrast", "mapping", "library"]
	required_wildcards_out_log = ["step", "extension", "contrast"]
	required_wildcards_in      = ["step", "extension", "sample"]
	
	def __init__(self, workflow, **kwargs):
		super().__init__(workflow, **kwargs)
		
		self.contrasts         = self.snakemake_workflow.config["contrasts"]["contrast_list"]
		self.contrast_defaults = self.snakemake_workflow.config["contrasts"]["defaults"]
		self.contrast_ids = [self._make_name(contr["title"]) + "_ID{}".format(index) for index, contr in enumerate(self.contrasts)] #used in pipeline to fill in wildcards
		
		# write log with parsed values to file for debugging
		self._write_log(contrast="all")
		
		
	#---------------------------------------------------- helper methods ----------------------------------------------------#
		
	def _get_filtered_wildcard_values(self, step, extension):
		""" return wildcard values filtered by step and extension """
		
		wildcard_values_filtered = {key:[] for key in self.wildcard_values}
		wildcard_comb_num = len(self.wildcard_values["sample"])
		assert all([len(val)==wildcard_comb_num for val in self.wildcard_values.values()])
		for index in range(wildcard_comb_num):
			if self.wildcard_values["step"][index]==step and self.wildcard_values["extension"][index]==extension:
				for key in self.wildcard_values: wildcard_values_filtered[key].append(self.wildcard_values[key][index])
		return wildcard_values_filtered
		
	def _write_log(self, **kwargs):
		filename = self.file_path(step="DEPipelinePathHandler", extension="log", log=True, **kwargs)
		os.makedirs(os.path.dirname(filename), exist_ok=True)
		with open(filename, "w") as f:
			f.write("pattern input:\n   in path pattern: {}\n   out path pattern: {}\n   log path pattern: {}\n\n"
			"parsed wildcards:\n   in path wildcards: {}\n   out/log path wildcards: {}\n   out dir wildcards: {}\n\n"
			"wildcard values:\n   per wildcard: {}\n\n"
			"neutral replacement of optional wildcards:\n   out/log: {}\n   out dir: {}".format(self.in_path_pattern, self.out_path_pattern, self.log_path_pattern, 
			self.in_path_wildcards, self.out_path_wildcards, self.outdir_path_wildcards, self.wildcard_values, 
			self.opt_wildcard_placeholders, self.opt_wildc_placeh_outdir))
			
			
	#-------------------------------------------- methods used in snakemake file --------------------------------------------#
	
	def get_contrast(self, contrast_ID):
		""" 
		get contrast dict of parameters from contrast ID, as e.g. used in wildcard
		(merge defaults and contrast-specific parameters that are defined in config)
		"""
		def dict_merge(base_dict, add_dict):
			""" merge recursively into dict (and replace)"""
			for key, val in add_dict.items():
				if key in base_dict and isinstance(base_dict[key], dict) and isinstance(add_dict[key], Mapping):
					base_dict[key] = dict_merge(base_dict[key], add_dict[key])
				else:
					base_dict[key] = add_dict[key]
			return base_dict
			
		b_dct = deepcopy(self.contrast_defaults)
		a_dct = self.contrasts[int(contrast_ID.split("_ID")[-1])]
		
		return dict_merge(b_dct, a_dct)
		
	def get_contrast_id_dict(self):
		""" return dict of contrast titles to IDs """
		titles = [contr["title"] for contr in self.contrasts]
		if len(set(titles)) < len(titles):
			warn("titles of contrasts defined in config are not unique! using last ID for each title")
		return {contr["title"]: self.contrast_ids[i] for i,contr in enumerate(self.contrasts)}
	
	def file_path(self, step, extension, contrast="{contrast}", log=False, **kwargs):
		"""
		Generate single path for intermediate and output or log files.
		
		:param step:  Snakemake rule for which the paths are generated
		:param extension: file extension for the generated file path
		:param contrast: contrast ID to be included in the file path
		:param log: generate path to logfile if log==True otherwise generate path to output/intermediate file
		:param **kwargs: if used specify replacement for {batch}, {flowcell}, {lane}, etc. ...
		"""
		path_pattern = self.log_path_pattern if log else self.out_path_pattern
		kwargs_out = {key: kwargs[key] if key in kwargs else val for key, val in self.opt_wildcard_placeholders.items()}
		return path_pattern.format(step=step, extension=extension, contrast=contrast, **kwargs_out)
	
	def expand_path(self, step, extension, if_set=False, **kwargs):
		"""
		Generate multiple paths for intermediate and output files corresponding to different contrast IDs
		
		:param step:  Snakemake rule name (string) for which the paths are generated
		:param extension: file extension for the generated file path
		:param if_set: dict entry of a contrast used to filter for which to expand over; expand over all contrast if False
		:param **kwargs: replacement strings for optional wildcards, e.g. batch, flowcell, lane (see description)
		:returns: list of paths
		"""
		input_choices = set(self.input_choice)
		kwargs_out = {key: kwargs[key] if key in kwargs else val for key, val in self.opt_wildcard_placeholders.items()}
		expand_input_choices = list(input_choices & set(kwargs_out))
		choice_kwargs = [{expand_input_choices[i]: choice for i,choice in enumerate(choice_comb)} 
				for choice_comb in itertools.product(*[self.input_choice[eic] for eic in expand_input_choices])]
		
		paths = []
		for contr in (c_id for i, c_id in enumerate(self.contrast_ids) if any({k:v}==if_set for k,v in self.contrasts[i].items()) or not if_set):
			if not expand_input_choices:
				paths.append(self.file_path(step, extension, contrast=contr, **kwargs_out))
			else:
				kwargs_out = {key: val for key, val in kwargs_out.items() if key not in expand_input_choices}
				for kwargs_ch in choice_kwargs:
					paths.append(self.file_path(step, extension, contrast=contr, **kwargs_out, **kwargs_ch))
		return paths
		
	def log_generated_files(self, save_to="", **kwargs):
		"""
		Scan all generated output/intermediate files and list them in a file with their corresponding wildcard values.
		This file can be used e.g. by the pipeline Rmd report to identify paths to files without the handler.
		"""
		input_files = iglob(re.sub("{[^}./]+}",   "*", self.out_path_pattern))
		wildcards   =   re.findall("{([^}./]+)}",      self.out_path_pattern)
			
		table_cols = {w:[] for w in wildcards}
		table_cols["filename"] = []
		for inp in input_files:
			table_cols["filename"].append(inp)
			self._get_wildcard_values_from_file_path(inp, self.out_path_pattern, wildc_val=table_cols)
		assert all([len(val)==len(table_cols["filename"]) for val in table_cols.values()])
		
		table = pd.DataFrame(table_cols)
		if not save_to:
			save_to = self.file_path(step="DEPipelinePathHandler", extension="tsv", contrast="all", **kwargs)
		table.to_csv(save_to, sep="\t", index=False)


##################################################################################################################################
#----------------------------------------------------- covariate file tool ------------------------------------------------------#
##################################################################################################################################

class CovariateFileTool(PipelinePathHandler):
	""" Tool to generate a covariate file before running the pipeline; use same config as for DEPipelinePathHandler """
	
	allowed_wildcards          = DEPipelinePathHandler.allowed_wildcards
	required_wildcards_out_log = DEPipelinePathHandler.required_wildcards_out_log
	required_wildcards_in      = DEPipelinePathHandler.required_wildcards_in
	
	def __init__(self, config_yaml, *add_yaml):
		with open(config_yaml, 'r') as stream:
			try:
				config_dict = yaml.safe_load(stream)
			except yaml.YAMLError as exc:
				print(exc)
		for yml in add_yaml:
			with open(yml, 'r') as stream:
				try:
					config_dict = {**config_dict, **yaml.safe_load(stream)}
				except yaml.YAMLError as exc:
					print(exc)
					
		self.in_path_pattern = config_dict["pipeline_param"]["in_path_pattern"]
		
		self.wildcard_constraints = self._prepare_inpathpattern()
		
		self.wildcard_values = self._get_wildcard_values_from_input(self.in_path_pattern)

		self.opt_wildcard_placeholders = {w: "{{{}}}".format(w) for w in set(self.wildcard_values)-set(self.required_wildcards_in)}
		
		self.covariate_data = pd.DataFrame({"filename": [], "md5": [], "group":[], "replicate":[], "label":[]})
		
		
	#---------------------------------------------------- helper methods ----------------------------------------------------#
		
	def _get_wildcard_combinations(self, wildcard_values, step, extension):
		""" go through wildcard values and get combinations """
		
		combinations = []
		WildcardComb = namedtuple("WildcardComb", [s for s in wildcard_values])
		wildcard_comb_num = len(wildcard_values["sample"])
		assert all([len(val)==wildcard_comb_num for val in wildcard_values.values()])
		for index in range(wildcard_comb_num):
			if wildcard_values["step"][index]==step and wildcard_values["extension"][index]==extension:
				combinations.append(WildcardComb(**{key: wildcard_values[key][index] for key in wildcard_values}))
		return combinations
		
	def _get_mapping_input(self, step, extension, wildcards):
		"""
		Generate paths to input files from read mapping, e.g. from STAR or Salmon.
		
		:param wildcards: rule wildcards
		:returns: a list with paths to specified input files
		"""
		wildcard_combs = self._get_wildcard_combinations(self.wildcard_values, step, extension)
		wildcard_placeholders = {"sample":"{sample}", **self.opt_wildcard_placeholders}
		kwargs_out = {key: getattr(wildcards, key, val) for key, val in wildcard_placeholders.items()}
		pattern_list = []
		seen = set()
		for comb in wildcard_combs:
			kwargs_filled = {key: getattr(comb, key) if "{" in val or val not in self.wildcard_values[key] else val for key,val in kwargs_out.items()}
			kwargs_id_tup = tuple(kwargs_filled[key] for key in sorted(kwargs_filled))
			if kwargs_id_tup not in seen:
				seen.add(kwargs_id_tup)
				pattern_list.append(self.in_path_pattern.format(step=step, extension=extension, **kwargs_filled))
		return [path for pat in pattern_list for path in iglob(pat, recursive=True)]
		
		
	#---------------------------------------------------- access methods ----------------------------------------------------#
		
	def update_covariate_data(self, step, extension):
		""" 
		fill mandatory columns of the covariate data frame by searching the input path specified in the config file
		
		:param step:  Snakemake rule name for which the files are searched
		:param extension: file extension of the searched files
		"""
		files = self._get_mapping_input(step, extension, wildcards={})
		md5   = [self._md5(f) for f in files]
		group = [self._get_wildcard_values_from_file_path(f,self.in_path_pattern)["sample"][0] for f in files]
		
		replicate, num_g = [], {}
		for g in group:
			if g not in num_g: num_g[g] = 1
			else: num_g[g] += 1
			replicate.append(num_g[g])
		label = ["{}_{}".format(a,b) for a,b in zip(group, replicate)]
		
		self.covariate_data = pd.DataFrame({"filename":files, "md5":md5, "group":group, "replicate":replicate, "label":label})
		
	def add_column(self, name, levels):
		""" 
		add a custom column to the covariate data frame
		
		levels can be either a list (order important!) or a dictionary. If levels is a dictionary it can be of two forms:
		- {<group(string)>: <level(string)>, ...}
		- {<level(string)>: [<group(string)>, <group(string)>, ...], ...}
		
		:param name:  name of the column
		:param levels: levels of the column
		"""
		if type(levels)==dict:
			if any(l not in self.covariate_data.group for l in levels):
				self.covariate_data[name] = [l for g in self.covariate_data.group for l,gs in levels.items() if g in gs]
			else:
				self.covariate_data[name] = [levels[g] for g in self.covariate_data.group]
		elif type(levels)==list:
			self.covariate_data[name] = levels
		
	def write_covariate_file(self, filename):
		""" write covariate data frame to file """
		self.covariate_data.to_csv(filename, sep="\t", index=False)
		
	
##################################################################################################################################
#------------------------------------------------------ sample info tool --------------------------------------------------------#
##################################################################################################################################


class SampleInfoTool(PipelinePathHandler):
	""" Tool to generate a sample info file before running the pipeline; use same config as for MappingPipelinePathHandler """
	
	allowed_wildcards          = MappingPipelinePathHandler.allowed_wildcards
	required_wildcards_out_log = MappingPipelinePathHandler.required_wildcards_out_log
	required_wildcards_in      = MappingPipelinePathHandler.required_wildcards_in
	
	allowed_read_extensions    = [".fastq", ".fastq.gz"]
	
	def __init__(self, config_yaml, *add_yaml):
		with open(config_yaml, "r") as stream:
			try:
				config_dict = yaml.safe_load(stream)
			except yaml.YAMLError as exc:
				print(exc)
		for yml in add_yaml:
			with open(yml, "r") as stream:
				try:
					config_dict = {**config_dict, **yaml.safe_load(stream)}
				except yaml.YAMLError as exc:
					print(exc)
					
		self.in_path_pattern  = config_dict["pipeline_param"]["in_path_pattern"]
		
		self.wildcard_constraints = self._prepare_inpathpattern()
		
		self.sample_info = {}
		
		
	#---------------------------------------------------- helper methods ----------------------------------------------------#
		
	def _get_wildcard_values_from_read_input(self, unix_style=True):
		""" go through files in input path and get values matching the wildcards """
		
		glob_pattern  =      re.sub("{[^}./]+}",           "*",                           self.in_path_pattern)
		wildcards     =  re.findall("{([^}./]+)}",                                        self.in_path_pattern)
		match_pattern =      re.sub("\\\\{([^}./]+)\\\\}", self._wildc_replace, re.escape(self.in_path_pattern))
		input_files   = glob(glob_pattern + ("*" if glob_pattern[-1]!="*" else ""), recursive=True)
		if unix_style:
			match_pattern = re.sub(r"\\\*\\\*",            "[^{}]*",   match_pattern)
			match_pattern = re.sub(r"(?<!\[\^{}\]\*)\\\*", "[^{}./]*", match_pattern)
		
		wildcard_values = {w:[] for w in wildcards}
		for inp in input_files:
			self._get_wildcard_values_from_file_path(inp, self.in_path_pattern, wildc_val=wildcard_values, unix_style=unix_style)
			
		print("\ninput files:\n{}".format("\n".join(input_files)))
		print(f"\nmatch pattern:\n{match_pattern}")
		return {**wildcard_values, "read_extension": [f.replace(re.match(match_pattern,f).group(0), "") for f in input_files]}
		
	def _convert_str_entries_to_lists(self, key="paired_end_extensions"):
		""" for importing lists from table entries """
		for smpl_info in self.sample_info.values():
			smpl_info[key] = [s.replace("'","").replace('"',"") for s in re.findall("[^\[\]\s,]+", smpl_info[key])]
			
	def _add_info_fields(self, add_dict):
		""" add fields from add_dict to self.sample_info if they are not already present """
		for sample, fields in add_dict.items():
			if sample in self.sample_info:
				s_info = self.sample_info[sample]
				for f_key,f_val in fields.items():
					if f_key not in s_info: s_info[f_key] = f_val
			
		
	#---------------------------------------------------- access methods ----------------------------------------------------#
		
	def update_sample_info(self, library_default="unstranded", add=False):
		""" 
		fill mandatory info about sample by searching the input path specified in the config file.
		
		attention: stranded is initially set to library_default for all samples! 
		This information has to be edited manually in table or yaml, if libraries were prepared differently.
		
		:param library_default:  options: ["unstranded", "forward", "reverse"]
		"""
		wildcard_values = self._get_wildcard_values_from_read_input()
		wildcard_combs  = [comb for comb in self._get_wildcard_combinations(wildcard_values) if comb.read_extension in self.allowed_read_extensions]
		print("\nextracted combinations:\n{}".format("\n".join("\t".join(i) for i in [wildcard_combs[0]._fields] + wildcard_combs)))
		
		sample_info = {}
		for comb in wildcard_combs:
			if comb.sample not in sample_info:
				sample_info[comb.sample] = {"stranded":library_default, "read_extension":comb.read_extension}
				sample_info[comb.sample]["paired_end_extensions"] = [getattr(comb, "mate", "")]
			elif hasattr(comb, "mate"):
				paired_end_ext     = getattr(comb, "mate", "")
				paired_end_ext_lst = sample_info[comb.sample]["paired_end_extensions"]
				if paired_end_ext_lst == [""]:
					raise ValueError("Error compiling sample information: sample {} has names with and without paired end extensions".format(comb.sample))
				if paired_end_ext not in paired_end_ext_lst:
					paired_end_ext_lst.append(paired_end_ext)
		if add:
			# add missing fields
			self._add_info_fields(sample_info)
		else:
			# overwrite
			self.sample_info = sample_info
		
	def write_table(self, filename, sep="\t"):
		"""
		write sample info to table
		"""
		tab = pd.DataFrame(self.sample_info).transpose()
		tab.to_csv(filename, sep=sep)
		
	def read_table(self, filename, sep="\t"):
		"""
		read sample info from table
		"""
		tab = pd.read_csv(filename, sep=sep, index_col=0).transpose()
		self.sample_info = tab.to_dict()
		self._convert_str_entries_to_lists("paired_end_extensions")

	def write_yaml(self, filename):
		"""
		write sample info to yaml
		"""
		with open(filename, "w") as f:
			yaml.dump({"sample_info": self.sample_info}, f, default_flow_style=False)
			
	def read_yaml(self, filename):
		"""
		read sample info from yaml
		"""
		with open(filename, "r") as f:
			try:
				self.sample_info = yaml.safe_load(f)["sample_info"]
			except yaml.YAMLError as exc:
				print(exc)
				
	def parse_isatab(self, filename):
		"""
		parse sample info from ISA-tab table
		"""
		tab = pd.read_csv(filename, sep="\t", index_col=False)
		tab.dropna(axis="columns", how="all", inplace=True)
		parse_conf = Path(os.path.realpath(__file__)).parent / Path("ISAtab_parse_conf.yaml")
		with open(str(parse_conf), "r") as stream:
			try:
				parse_conf = yaml.safe_load(stream)
			except yaml.YAMLError as exc:
				print(exc)
		
		def find_column(key):
			if not parse_conf[key]["columns"]: return None
			for col in tab.columns:
				for col_regex in parse_conf[key]["columns"]:
					if re.match(col_regex, col): return col
		def map_value(key, value):
			if "as_is" in parse_conf[key] and parse_conf[key]["as_is"]: return value
			value = str(value).lower()
			for val_regex, val_repl in parse_conf[key]["value"].items():
				if re.fullmatch(val_regex, value): 
					return re.sub(val_regex, val_repl, value) if isinstance(val_repl, str) else val_repl
		
		sample_info_cols = {}
		for key in parse_conf:
			column = find_column(key)
			print(f"key: {key} | column: {column}")
			if column is not None: 
				sample_info_cols[key] = [map_value(key, val) for val in tab[column]]
		
		self.sample_info = {sample_info_cols["id"][i]: 
		                   	{key: val[i] for key,val in sample_info_cols.items() if key != "id"} 
		                   for i in range(len(sample_info_cols["id"]))}
			
		
##################################################################################################################################
#------------------------------------------------------ class for report --------------------------------------------------------#
##################################################################################################################################

class ReportTool(PipelinePathHandler):

	entry_heading_pattern = "#HSTART.*\n(.*{{ENTRY_NAME}}.*)\n.*#HEND"
	insert_pattern        = "#>.*INSERT.*<#"
	entry_name_wildcard   = "{{ENTRY_NAME}}", "{{ENTRY_ID}}"

	def __init__(self, pph, id_dicts={}):
		#if type(config) is dict:
		#	config_dict = config
		#elif isinstance(config, str):
		#	with open(config, "r") as stream:
		#		try:
		#			config_dict = yaml.safe_load(stream)
		#		except yaml.YAMLError as exc:
		#			print(exc)
		#else:
		#	raise TypeError("Wrong type of argument config: must be dict or str.")
		self.path_handler = pph
		config_dict       = self.path_handler.snakemake_workflow.config
			
		if config_dict["pipeline_param"]["report_snippets"]:
			self.report_snippet_base_dir = Path(config_dict["pipeline_param"]["report_snippets"])
		else:
			self.report_snippet_base_dir = Path(sys.path[0]) / "report"
		
		self.report_snippet_building_plan = config_dict["report"]["report_snippets"]
		
		self.report_snippet_defaults      = config_dict["report"]["defaults"]
		
		self.contrast_list = [contr["title"] for contr in config_dict["contrasts"]["contrast_list"]] if "contrasts" in config_dict else ""
		
		# define substitutions for report generation
		self.substitutions = {"__contrasts__": self.contrast_list}
		
		self.id_dicts = id_dicts
	
	
	#---------------------------------------------------- helper methods ----------------------------------------------------#
		
	def _split_template(self, template_text):
		before, after = re.split(self.insert_pattern, template_text, maxsplit=1, flags=re.MULTILINE)
		return (before, after)
		
	def _make_id(self, name):
		entry_name, templ_name = name
		if templ_name in self.id_dicts:
			return self.id_dicts[templ_name][entry_name]
		else:
			return self._make_name(entry_name)
		
	def _insert_entry_name(self, text, name):
		entry_name = name[0]
		return text.replace(self.entry_name_wildcard[0], entry_name).replace(self.entry_name_wildcard[1], self._make_id(name))
		
	def _get_entry_heading_code(self, template_text):
		return re.search(self.entry_heading_pattern, template_text).group(1)
		
	def _get_entry_heading(self, heading_code, name):
		return heading_code.replace(self.entry_name_wildcard[0], name)
		
	def _rem_entry_heading_code(self, template_text):
		return re.sub(self.entry_heading_pattern, "", template_text)
		
	def _assemble_entries(self, entries, path, snippet_name, entry_heading_code):
		""" assemble entry list (e.g. list of contrasts) """
	
		entry_text = []
		
		# add subsection entries
		if type(entries) is str: entries = self.substitutions[entries] if entries in self.substitutions else [entries]
		for entry in entries:
			if type(entry) is str:
			
				entry_name = entry
				if snippet_name not in self.report_snippet_defaults:
					raise KeyError("Error compiling report snippets for {snip} {entr}! (no snippets provided and "
							"key {snip} not found in config defaults)".format(snip=snippet_name, entr=entry))
				sub_snippet_list = self.report_snippet_defaults[snippet_name]
				
			elif isinstance(entry, dict):
			
				assert len(entry)==1
				entry_name = list(entry.keys())[0]
				
				sub_snippet_list  = entry[entry_name]
				
			else:
				raise TypeError("Error in report snippet building plan! (expected str or dict, got {})".format(type(entry)))
			
			entry_text += [self._get_entry_heading(entry_heading_code, entry_name), "\n\n" , 
					self._assemble_template(sub_snippet_list, path, snippet_name, entry_heading_code, (entry_name, snippet_name))]
			
		return "".join(entry_text)
		
	def _assemble_template(self, snippet_list, path, snippet_name, entry_heading_code, entry=("","")):
		""" assemble snippet list """
		
		snippet_text = []
		
		# add snippets
		for snippet in snippet_list:
			if type(snippet) is str:
				
				snippet_file = path / snippet
				snippet_cont = self._insert_entry_name(snippet_file.read_text(), entry)
				
				requirements = re.findall("(?<=#REQUIRE)\s+{{(\S+)}}", snippet_cont)
				snippet_cont = re.sub(    "#REQUIRE\s+{{\S+}}\n+", "", snippet_cont)
				req_fields   = ["step","extension","contrast"]
				
				if all(Path(self.path_handler.file_path( **dict(zip( req_fields, req.split("-") )) )).exists() for req in requirements):
					snippet_text.append(snippet_cont)
					
			elif isinstance(snippet, dict):
				assert len(snippet)==1
				snippet_key   = list(snippet.keys())[0]
				snippet_value = list(snippet.values())[0]
				
				if snippet_key == "__list__":
					# add list entries
					snippet_text.append( self._assemble_entries(snippet_value, path, snippet_name, entry_heading_code) )
				else:
					# load template
					sub_template_path  = path / snippet_key / (snippet_key + "_main_template.Rmd")
					sub_template_text  = sub_template_path.read_text()
					# extract heading
					entry_heading_code = self._get_entry_heading_code(sub_template_text)
					sub_template_text  = self._rem_entry_heading_code(sub_template_text)
					# prepare for insert
					temp_begin, temp_end = self._split_template(sub_template_text)
					
					# add sub-section
					sub_section_text = self._assemble_template(snippet_value, path/snippet_key, snippet_key, entry_heading_code, entry)
					snippet_text.append(temp_begin + sub_section_text + temp_end)
				
			else:
				raise TypeError("Error in report snippet building plan! (expected str or dict, got {})".format(type(snippet)))
			
		return  "\n".join(snippet_text)
		
		
	#---------------------------------------------------- access methods ----------------------------------------------------#
	
	def generate_report(self):
		"""
		Generate a report by assembling Rmd snippets as defined in the config.
		
		Starts in report directory and recursively alternates between concatenating lists of snippets 
		and creating sub-lists of entries (e.g. contrasts).
		"""
		template_path = self.report_snippet_base_dir / "report_main_template.Rmd"
		template_text = template_path.read_text()
		
		# generate report
		temp_begin, temp_end = self._split_template(template_text)
		report_text = self._assemble_template(self.report_snippet_building_plan, self.report_snippet_base_dir, 
							snippet_name="other", entry_heading_code="# {{ENTRY_NAME}}")
		
		return temp_begin + report_text + temp_end








