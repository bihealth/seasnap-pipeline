## SeA-SnaP tools used in pipelines
## version: 0.9.5
## author: J.P.Pett (patrick.pett@bihealth.de)

import sys, os, re, shutil, hashlib, itertools, yaml, pandas as pd
from builtins import isinstance, TypeError
from collections.abc import Mapping, Iterable
from collections import namedtuple, OrderedDict
from types import SimpleNamespace
from contextlib import contextmanager
from copy import deepcopy
from time import strftime
from warnings import warn
import warnings
from pathlib import Path
from glob import iglob, glob

# from snakemake.io import glob_wildcards

yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_dict(dict(data)))

warnings.simplefilter("always")


##################################################################################################################################
# ---------------------------------------------- base class for handling file paths ----------------------------------------------#
##################################################################################################################################

class PipelinePathHandler:
	"""
	Generates file paths from patterns defined in config.
	
	Loads and checks paths from config file upon initialization and provides
	methods to fill and expand the wildcards.
	"""

	allowed_wildcards = ["step", "extension"]
	required_wildcards_out_log = ["step", "extension"]
	required_wildcards_in = []
	wildcard_fix_values = {}

	def __init__(self, workflow, test_config=False, test_allowed_wildcards=True):
		self.test_config = self._load_test_config(test_config)
		if self.test_config:
			self._test_config_general(workflow.config, self.test_config)

		self.snakemake_workflow = workflow

		self.out_path_pattern = self.snakemake_workflow.config["pipeline_param"]["out_path_pattern"]
		self.log_path_pattern = self.snakemake_workflow.config["pipeline_param"]["log_path_pattern"]
		self.in_path_pattern = self.snakemake_workflow.config["pipeline_param"]["in_path_pattern"]

		self.wildcard_constraints = self._prepare_inpathpattern()

		self.out_dir_pattern = "/".join(self.out_path_pattern.split("/")[:-1])

		self.out_path_wildcards = self._get_wildcard_list(self.out_path_pattern)
		self.log_path_wildcards = self._get_wildcard_list(self.log_path_pattern)
		self.in_path_wildcards = self._get_wildcard_list(self.in_path_pattern)
		self.outdir_path_wildcards = self._get_wildcard_list(self.out_dir_pattern)

		self.input_choice = self._set_input_choice(self.snakemake_workflow.config)

		# check consistency of loaded input and output patterns
		if self.test_config:
			self._test_config_input(test_allowed_wildcards)

		# wildcard values
		self.wildcard_values = self._get_wildcard_values_from_input(self.in_path_pattern)

		# neutral replacement for optional wildcards
		self.opt_wildcard_placeholders = {w: "{{{}}}".format(w) for w in
										  set(self.in_path_wildcards) - set(self.required_wildcards_in)}

		# path pattern of out directory may miss some wildcards
		placeholders_in_outdir_path = set(self.outdir_path_wildcards) & set(self.opt_wildcard_placeholders.keys())
		self.opt_wildc_placeh_outdir = {key: self.opt_wildcard_placeholders[key] for key in placeholders_in_outdir_path}

	# ---------------------------------------------------- helper methods ----------------------------------------------------#

	def _load_test_config(self, test_config):
		""" load test_config for validity check of config values """
		if test_config:
			if type(test_config) is dict:
				return test_config
			elif isinstance(test_config, str):
				with open(test_config, "r") as stream:
					try:
						return yaml.safe_load(stream)
					except yaml.YAMLError as exc:
						print(exc)
			else:
				raise TypeError("Wrong type of argument test_config: must be dict or str.")
		else:
			return None

	def _test_config_input(self, test_allowed_wildcards):
		""" test whether wildcards in path patterns are valid """
		if test_allowed_wildcards:
			# test if all wildcards allowed
			if not all(x in self.allowed_wildcards for x in
					   self.out_path_wildcards + self.log_path_wildcards + self.in_path_wildcards):
				raise ValueError(
					"Error in config file: unknown wildcards. allowed wildcards: {}".format(self.allowed_wildcards))

		# test for required wildcards in all patterns
		if not all(
				x in self.out_path_wildcards and x in self.log_path_wildcards for x in self.required_wildcards_out_log):
			raise ValueError("Error in config file: 'step', 'extension', 'sample' and 'name' must be wildcards "
							 "in out_path_pattern and log_path_pattern")
		if not all(x in self.in_path_wildcards for x in self.required_wildcards_in):
			raise ValueError("Error in config file: 'step', 'extension', 'sample' and 'name' must be wildcards "
							 "in out_path_pattern and log_path_pattern")

	def _test_config_general(self, base_dict, check_dict):
		""" test whether values set in the config are valid """
		for key, val in check_dict.items():

			## key test
			# if wildcard
			if key == "__any__":
				for glob_key in base_dict:
					if glob_key not in check_dict:
						self._test_config_general(base_dict, {glob_key: val})
			# if reference
			elif key == "__any_other__":
				if isinstance(val, str):
					ref_dict = self.test_config
					for r_key in val.split(":"): ref_dict = ref_dict[r_key]
				elif isinstance(val, dict):
					ref_dict = val
				for opt_key in ref_dict:
					if opt_key in base_dict:
						self._test_config_general(base_dict, {opt_key: ref_dict[opt_key]})
			# if missing
			elif key not in base_dict and key not in ["__num__", "__opt__"]:
				raise KeyError("Error in config file: The required key '{}' was not defined!".format(key))

			## value test
			# if string
			elif isinstance(val, str):
				if isinstance(base_dict[key], str):
					if not re.fullmatch(val, base_dict[key], re.DOTALL):
						raise ValueError("Error in config file: value of '{}' ('{}') does not match '{}'!".format(key,
																												  base_dict[
																													  key],
																												  val))
				else:
					raise TypeError(
						"Error in config file: value of '{}' should be a string! got: {}".format(key, base_dict[key]))
			# if list
			elif isinstance(val, list):
				if isinstance(base_dict[key], list):
					for item in base_dict[key]:
						self._test_config_general({key: item}, {key: val[0]})
				else:
					raise TypeError(
						"Error in config file: value of '{}' should be a list! got: {}".format(key, base_dict[key]))
			# if None
			elif val is None:
				if base_dict[key] is not None:
					raise TypeError(
						"Error in config file: value of '{}' should be None! got: {}".format(key, base_dict[key]))
			# if bool
			elif isinstance(val, bool):
				if not isinstance(base_dict[key], bool):
					raise TypeError(
						"Error in config file: value of '{}' should be bool! got: {}".format(key, base_dict[key]))
			# if dict
			elif isinstance(val, dict):
				# if number range
				if len(val) == 1 and list(val)[0] == "__num__":
					if not isinstance(base_dict[key], (int, float)):
						raise TypeError("Error in config file: value of '{}' must be a number!".format(key))
					num_range = list(val.values())[0]
					assert len(num_range) == 2
					if num_range[0] and not base_dict[key] >= num_range[0]:
						raise ValueError("Error in config file: value of '{}' must be >={}!".format(key, num_range[0]))
					if num_range[1] and not base_dict[key] <= num_range[1]:
						raise ValueError("Error in config file: value of '{}' must be <={}!".format(key, num_range[1]))
				else:
					others = {}
					for v_key, v_val in val.items():
						# if option list (OR)
						if v_key == "__opt__":
							options, error_num = v_val, 0
							for option in options:
								try:
									self._test_config_general(base_dict, {key: option})
								except (KeyError, ValueError, TypeError) as err:
									if not re.match("\"?Error in config file", str(err)):
										raise
									error_num += 1
								except:
									raise
							if error_num == len(options):
								raise ValueError("Error in config file: no valid option for '{}' ('{}')! "
												 "must be one of: {}".format(key, base_dict[key], options))
						else:
							others[v_key] = v_val
					# others (AND)
					if others:
						if isinstance(base_dict[key], dict):
							self._test_config_general(base_dict[key], others)
						else:
							raise TypeError("Error in config file: value of '{}' must be a dict!".format(key))

	def _set_input_choice(self, config):
		""" create input choice dictionary from config ensuring a standardized structure """
		input_choice = {}
		if "input_choice" in config["pipeline_param"] and isinstance(config["pipeline_param"]["input_choice"], dict):
			for key, val in config["pipeline_param"]["input_choice"].items():
				if val:
					if isinstance(val, str):
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

	@staticmethod
	def _is_older(file1, file2):
		file1, file2 = Path(file1), Path(file2)
		assert file1.is_file()
		if file2.is_file():
			return file1.stat().st_mtime < file2.stat().st_mtime
		elif file2.is_dir():
			newest_f = float("inf")
			for p in file2.rglob("*"):
				if p.stat().st_mtime < newest_f:
					newest_f = p.stat().st_mtime
			return file1.stat().st_mtime < newest_f
		else:
			raise TypeError(f"Path {file2} is not a file or directry!")

	@staticmethod
	def _get_names_values(data):
		if isinstance(data, Mapping):
			data_values = list(data.values())
			data_keys = list(data.keys())
		elif isinstance(data, Iterable):
			data_values = list(data)
			data_keys = None
		else:
			raise TypeError(f"Cannot get names and values for type {type(data)}!")
		return data_keys, data_values

	def _get_wildcard_combinations(self, wildcard_values):
		""" go through wildcard values and get combinations """

		combinations = []
		WildcardComb = namedtuple("WildcardComb", [s for s in wildcard_values])
		wildcard_comb_num = len(wildcard_values["sample"])
		assert all([len(val) == wildcard_comb_num for val in wildcard_values.values()])
		for index in range(wildcard_comb_num):
			combinations.append(WildcardComb(**{key: wildcard_values[key][index] for key in wildcard_values}))
		return combinations

	def _prepare_inpathpattern(self):
		""" read and remove wildcard constraints from in_path_pattern """
		wildcards = re.findall("{([^{}]+)}", self.in_path_pattern)
		wildcard_constraints = {}
		for wildcard in wildcards:
			comp = wildcard.split(",")
			if len(comp) > 1:
				wildcard_constraints[comp[0]] = comp[1] + "|" + self.wildcard_fix_values[comp[0]]
				self.in_path_pattern = self.in_path_pattern.replace(wildcard, comp[0])
		return wildcard_constraints

	def _get_wildcard_fix_values(self, fix):
		inv = ["!", "-"]  # invert selection if these characters precede wildcard name
		if isinstance(fix, str):
			if fix == "all":
				return self.wildcard_fix_values
			else:
				fix = [fix]
		if isinstance(fix, list):
			include = set([w for w in fix if w[0] not in inv])
			exclude = set([w[1:] for w in fix if w[0] in inv])
			if exclude and not include:
				include = set(self.wildcard_fix_values)  # all wildcards
			for wcd in include | exclude:
				if wcd not in self.allowed_wildcards:
					raise ValueError(f"Unknown wildcard: '{wcd}'!")
			use_wcd = include - exclude
			return {w: self.wildcard_fix_values[w] for w in use_wcd}
		else:
			raise ValueError(f"Wrong argument for 'fix': {fix}.")

	def _get_wildcard_values_from_input(self, input_pattern, unix_style=True, verbose=False):
		""" go through files in input path and get values matching the wildcards """

		glob_pattern = re.sub("{[^}./]+}", "*", input_pattern)
		wildcards = re.findall("{([^}./]+)}", input_pattern)
		input_files = glob(glob_pattern + ("*" if glob_pattern[-1] != "*" else ""), recursive=True)
		match_tup = self._get_match_pattern_and_wildcards(input_pattern, unix_style)

		if verbose:
			print("\ninput files:\n{}".format("\n".join(input_files)))
			print(f"\nmatch pattern:\n{match_tup[0]}")

		wildcard_values = {w: [] for w in wildcards}
		for inp in input_files:
			self._get_wildcard_values_from_file_path(inp, match_tup=match_tup, wildc_val=wildcard_values)
		return wildcard_values

	def _wildc_replace(self, matchobj):
		""" method used with re.sub to generate match pattern from path pattern """
		wildc_name = matchobj.group(1)
		if wildc_name in self.wildcard_constraints:
			return "({})".format(self.wildcard_constraints[wildc_name].replace("//", "/"))
		elif wildc_name == "extension":
			return "([^}/]+)"
		else:
			return "([^}./]+)"

	def _get_wildcard_values_from_file_path(
			self, filename, input_pattern=None, match_tup=None, wildc_val=None, unix_style=True, verbose=False
		):
		"""
		get values matching wildcards from given file path
		:param match_tup: tuple of (match_pattern, wildcards)
		"""
		if not input_pattern and not match_tup:
			raise ValueError("Need either 'input_pattern' or 'match_tup' as argument.")

		if match_tup:
			match_pattern, wildcards = match_tup
		else:
			match_pattern, wildcards = self._get_match_pattern_and_wildcards(input_pattern, unix_style)

		if verbose:
			print(f"\nmatch pattern:\n{match_pattern}")

		wildcard_values = wildc_val if wildc_val else {w: [] for w in wildcards}

		match_obj = re.match(match_pattern, filename)
		if not match_obj:
			warn(
				f"Skipping file in working dir: '{filename}', because not matching match pattern: '{match_pattern}' ..will not be tracked.")
			found = False
		else:
			matches = match_obj.groups()
			assert len(matches) == len(wildcards)
			seen = set()
			for index, wildc in enumerate(wildcards):
				if not wildc in seen:
					wildcard_values[wildc].append(matches[index])
					seen.add(wildc)
			found = True

		return found, wildcard_values

	def _get_match_pattern_and_wildcards(self, input_pattern, unix_style):
		match_pattern = re.sub("\\\\{([^}./]+)\\\\}", self._wildc_replace, re.escape(input_pattern))
		if unix_style:
			match_pattern = re.sub(r"\\\*\\\*", "[^{}]*", match_pattern)
			match_pattern = re.sub(r"(?<!\[\^{}\]\*)\\\*", "[^{}./]*", match_pattern)
		wildcards = re.findall("{([^}./]+)}", input_pattern)
		return match_pattern, wildcards

	def _collect_generated_files(self, path_pattern=None):
		if not path_pattern: path_pattern = self.out_path_pattern

		input_files = iglob(re.sub("{[^}./]+}", "*", path_pattern))
		wildcards = re.findall("{([^}./]+)}", path_pattern)

		table_cols = {w: [] for w in wildcards}
		table_cols["filename"] = []
		for inp in input_files:
			found, _ = self._get_wildcard_values_from_file_path(inp, input_pattern=path_pattern, wildc_val=table_cols)
			if found:
				table_cols["filename"].append(inp)
		assert all([len(val) == len(table_cols["filename"]) for val in table_cols.values()])

		return table_cols

	def _choose_input(self, wildcards, choice_name, options, func):
		""" called by wrapper choose_input() """
		if wildcards and hasattr(wildcards, choice_name):
			input_from = getattr(wildcards, choice_name)
		elif choice_name in self.input_choice:
			inp_choice = self.input_choice[choice_name]
			if isinstance(inp_choice, str):
				input_from = inp_choice
			elif isinstance(inp_choice, list):
				input_from = inp_choice[0]
			else:
				raise TypeError(
					"Error choosing input: {} is not a valid type for input choice! "
					"(input_choice in config for {})".format(type(inp_choice), choice_name)
					)
		else:
			raise KeyError(
				"Error choosing input: no wildcard '{}' passed and no input choice for '{}' "
				"specified in config!".format(choice_name, choice_name)
				)

		for inp in options:
			inp = {**dict(wildcards.items()), **inp}
			if input_from == inp["step"]:
				if func == "file_path":
					return self.file_path(**inp)
				elif func == "out_dir_name":
					return self.out_dir_name(**inp)
				elif callable(func):
					return func(**inp)
				else:
					raise ValueError("Error choosing input: wrong func argument: {}".format(str(func)))
		raise ValueError(
			"Error choosing input: no valid mapping input type specified! (got: {})".format(input_from)
			)

	def load_config_from_path(self, path, path_handler=None):
		if not path_handler:
			path_handler = self
		config_file = path_handler.file_path("pipeline_report", "yaml", path_pattern=path, fix="all")
		with open(config_file, "r") as stream:
			try:
				config_dict = yaml.safe_load(stream)
			except yaml.YAMLError as exc:
				print(exc)
		return config_dict

	#-------------------------------------------- methods used in snakemake file --------------------------------------------#
		
	def choose_input(self, choice_name, options, func="file_path"):
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
		return lambda wildcards: self._choose_input(wildcards, choice_name, options, func=func)

	def wildcard_values_from(self, filepath, in_path_pattern=True):
		"""
		Parse a file path and return a dict of wildcard names to values.
		
		:param filepath: a file path (string)
		:param in_path_pattern: use in_path_pattern to extract wildcards, if False use out_path_pattern
		:returns: a dict of wildcard names to values
		"""
		if in_path_pattern:
			return self._get_wildcard_values_from_file_path(filepath, input_pattern=self.in_path_pattern)[1]
		else:
			return self._get_wildcard_values_from_file_path(filepath, input_pattern=self.out_path_pattern)[1]

	@staticmethod
	def get_r_repr(data, to_type=None, round_float=None):
		"""
		Translate python data structure into R representation for vector (i.e. c(...)) or list (i.e. list(...)).
		If to_type==None, choose c() or list() automatically. Named vector/list is created if data is of type
		Mapping (e.g. a dictionary).

		:param data: data to convert into R representation
		:param to_type: "vector" or "list"; if None, choose automatically
		:param round_float: round floats to N digits; if None do not round
		:returns: a string with R representation of data
		"""
		if not isinstance(data, Iterable):
			# single value
			if isinstance(data, bool):
				return "TRUE" if data else "FALSE"
			elif isinstance(data, (int, float)):
				return str(round(data, round_float) if round_float else data)
			elif data is None:
				return "NULL"
		elif isinstance(data, str):
			# string
			return '"' + data + '"'
		else:
			# other Iterable
			data_keys, data_values = PipelinePathHandler._get_names_values(data)
			# auto-determine target type to_type
			same_type = all(isinstance(i, type(data_values[0])) for i in data_values)
			no_iter_inside = not isinstance(data_values[0], Iterable) or isinstance(data_values[0],
																					str) if data_values else True
			use_type = "vector" if same_type and no_iter_inside else "list"
			if to_type:
				if to_type == "vector" and to_type != use_type:
					warn("Cannot use 'vector' for all Iterables, using 'list'!")
				else:
					use_type = to_type
			# set substitute string
			if use_type == "vector":
				subs = "c({})"
			elif use_type == "list":
				subs = "list({})"
			else:
				raise ValueError(f"Wrong value {use_type} for to_type!")
			# fill values
			data_values = (
				PipelinePathHandler.get_r_repr(value, round_float=round_float, to_type=to_type)
				for value in data_values
				)
			if data_keys:
				return subs.format(", ".join(f'"{k}"={v}' for k, v in zip(data_keys, data_values)))
			else:
				return subs.format(", ".join(data_values))
		raise TypeError(f"No R representation set for {type(data)}!")

	def file_path(self, *args, **kwargs):
		pass

	def out_dir_name(self, *args, **kwargs):
		pass

	def expand_path(self, *args, **kwargs):
		pass

	def log(self, out_log, script, step, extension, fix=None, **kwargs):
		"""
		create various log files and return a path to the script file for execution
		
		:param out_log: the main log file of rule 'step' that will contain output and error messages
		:param script: the script with wildcards filled by snakemake
		:param step: the rule for which logs shall be generated
		:param extension: the extension of the script file, e.g. '.R' or '.sh'
		:param fix: fix some wildcards; e.g. fix=["sample"] has the same effect as passing sample="all_samples", fix="all" fixes all wildcards
		:returns: path to a file containing the script
		"""
		script_file = self.file_path(step, extension, log=True, fix=fix, **kwargs)
		config_yaml = self.file_path(step, "config.yaml", log=True, fix=fix, **kwargs)
		conda_list = self.file_path(step, "conda_list.txt", log=True, fix=fix, **kwargs)
		conda_info = self.file_path(step, "conda_info.txt", log=True, fix=fix, **kwargs)
		conda_env = self.file_path(step, "conda_env.yaml", log=True, fix=fix, **kwargs)

		# write script to file
		with open(script_file, "w") as f: f.write(script)

		# dump merged config
		with open(config_yaml, "w") as f: yaml.dump(self.snakemake_workflow.config, f, default_flow_style=False)

		# write conda logs
		os.system("conda list > {}".format(conda_list))
		os.system("conda info > {}".format(conda_info))
		os.system("conda env export > {}".format(conda_env))

		# git and snakefile info to out_log
		git_dir = self.snakemake_workflow.basedir
		out_log_abs = str(Path(out_log).resolve())
		os.system('echo "--------------- git info ---------------" > {}'.format(out_log))
		os.system('cd {}; echo "name:" $(git rev-parse --show-toplevel) >> {}'.format(git_dir, out_log_abs))
		os.system('cd {}; echo "branch:" $(git symbolic-ref HEAD) >> {}'.format(git_dir, out_log_abs))
		os.system('cd {}; echo "hash:" $(git rev-parse HEAD) >> {}'.format(git_dir, out_log_abs))
		os.system('cd {}; {{ echo "status:"; git status -sb; }} >> {}'.format(git_dir, out_log_abs))
		os.system('echo "--------------- snakefile --------------" >> {}'.format(out_log))
		os.system('echo "snakefile name: {}" >> {}'.format(self.snakemake_workflow.snakefile, out_log))
		os.system('echo "snakefile md5 : {}" >> {}'.format(self._md5(self.snakemake_workflow.snakefile), out_log))
		os.system('echo "----------------------------------------" >> "{l}"; echo >> "{l}"'.format(l=out_log))

		return script_file

	def export(self, config_key="export"):
		"""
		export selected results into a separate folder structure (as configured in config file)
		"""
		export_spec = self.snakemake_workflow.config[config_key]
		blueprint = export_spec["blueprint"]
		none_context = contextmanager(lambda: iter([None]))()
		with (open(blueprint["file"], "w") if blueprint and blueprint["file"] else none_context) as bp_out:
			for pattern in export_spec["path_pattern"]:
				# --- go through path patterns (that specify which location to copy to)
				pat = strftime(pattern).replace(
					"{GENOME}",
					self.snakemake_workflow.config["organism"]["genome_version"]
					)
				wildcards = self._get_wildcard_list(pat)
				key_wcs = [wc for wc in wildcards if wc[:6] == "files:"]
				assert len(key_wcs) == 1  # no nested keys for now
				key = key_wcs[0].split(":")[1:]
				assert len(key) > 0
				pat = pat.replace("{{{}}}".format(key_wcs[0]), key[0])
				# --- read specification for fetching files
				for opt_dct in export_spec["_".join(key)]:
					# --- read options
					if "files" in opt_dct:
						mode = "files"
					elif "dir" in opt_dct:
						mode = "dir"
					else:
						raise ValueError("Error in export: no mode (files or dir) specfied in config!")
					compress = opt_dct["compress"] if "compress" in opt_dct else None
					compress_list = opt_dct["compress_list"] if "compress_list" in opt_dct else []
					exclude = opt_dct["exclude"] if "exclude" in opt_dct else []
					tar_excl = " ".join(f'--exclude="{i}"' for i in exclude)
					arg_dct = opt_dct[mode]
					suffix = opt_dct["suffix"] if "suffix" in opt_dct else None

					# --- extract files
					extra_wcs = set(wildcards) - set(key_wcs) - (set(arg_dct) & set(wildcards))
					assert len(extra_wcs) <= 1  # only either 'sample' (mapping) or 'contrast' (DE) for now
					if extra_wcs:
						extra_wc = list(extra_wcs)[0]
						wc_in_dct = {k: v for k, v in arg_dct.items() if k in wildcards}
						search_pat = self.out_path_pattern if "log" not in arg_dct else self.log_path_pattern

						if mode == "files":
							sourcef = self.expand_path(**arg_dct)
						else:
							sourcef = self.file_path(**{**arg_dct, "extension": "{extension}"})
							if compress and not compress_list:
								list_dir = Path(sourcef).parent / ""
								if suffix: list_dir = list_dir / suffix
								sourcef = glob(str(list_dir).replace(f"{{{extra_wc}}}", "*"))
							else:
								source_tmp = []
								list_dir = Path(sourcef).parent
								if suffix: list_dir = list_dir / suffix
								dirs = glob(str(list_dir).replace(f"{{{extra_wc}}}", "*"))
								for d in dirs:
									if Path(d).is_file():
										source_tmp.append(str(d))
									else:
										for fp in Path(d).iterdir():
											f = str(fp)
											if os.path.isdir(f) and str(fp.name) not in compress_list:
												source_tmp.extend(
													glob(str(Path(f) / "**"), recursive=True)
												)
											# elif fp.with_suffix("").stem not in compress_list or os.path.isdir(f):
											else:
												source_tmp.append(f)
								sourcef = source_tmp
							search_pat = str(Path(search_pat).parent)

						get_wc = self._get_wildcard_values_from_file_path
						target = [
							pat.format(**{**wc_in_dct, extra_wc: get_wc(src, input_pattern=search_pat)[1][extra_wc][0]})
							for src in sourcef
							]
					else:
						sourcef = [self.file_path(**arg_dct)]
						target = [pat.format(**{k: v for k, v in arg_dct.items() if k in wildcards})]
					# --- copy files or write blueprint
					assert len(sourcef) == len(target)
					for i in range(len(sourcef)):
						if mode == "dir" and (not compress or compress_list):
							target[i] = str(Path(target[i]) / Path(sourcef[i]).name)
						f_src, f_trg = Path(sourcef[i]).resolve(), Path(target[i])
						if f_src.exists() and f_src.suffix != ".md5":
							print("\n...copy {} to {} ...\n".format(sourcef[i], target[i]))
							if not compress_list or f_src.name in compress_list:
								if compress == "zip":
									to_zip = str(f_src.with_suffix('.zip'))
									if not Path(to_zip).exists() or self._is_older(to_zip, f_src):
										print(f"compression: create {to_zip} ...")
										os.system(f"cd {str(f_src.parent)}; zip -r {to_zip} {str(f_src.name)}")
									sourcef[i], f_src = to_zip, Path(to_zip)
									if mode == "dir" and compress_list:
										zip_trg = f_trg.with_suffix(f_trg.suffix + '.zip')
										target[i], f_trg = str(zip_trg), zip_trg
								elif compress == "tar":
									to_tar = str(f_src.with_suffix('.tar.gz'))
									if not Path(to_tar).exists() or self._is_older(to_tar, f_src):
										print(f"compression: create {to_tar} ...")
										if extra_wcs:
											print(extra_wc)
											extra_wc_val = self._get_wildcard_values_from_file_path(
												sourcef[i],
												input_pattern=os.path.join(self.out_dir_pattern, "**"),
												unix_style=True
											)[1]
											print(extra_wc_val)
											tar_excl_upd = tar_excl.replace(f"{{{extra_wc}}}", extra_wc_val[extra_wc][0])
										else:
											tar_excl_upd = tar_excl
										os.system(f"cd {str(f_src.parent)}; tar -czf {to_tar} {tar_excl_upd} {str(f_src.name)}")
									sourcef[i], f_src = to_tar, Path(to_tar)
									if mode == "dir" and compress_list:
										tar_trg = f_trg.with_suffix(f_trg.suffix + '.tar.gz')
										target[i], f_trg = str(tar_trg), tar_trg
							if bp_out:
								print(blueprint["command"].format(
									src=f_src,
									dest=target[i]
									), file=bp_out)
								# -- create md5 sum
								md5_path = Path(sourcef[i] + ".md5")
								if not md5_path.exists() or self._is_older(md5_path, sourcef[i]):
									md5_path.write_text(self._md5(sourcef[i]))
								print(blueprint["command"].format(
									src=md5_path.resolve(),
									dest=f_trg.with_suffix(f_trg.suffix + ".md5")
									), file=bp_out)
							else:
								Path(target[i]).parent.mkdir(exist_ok=True, parents=True)
								shutil.copy2(sourcef[i], target[i])
								Path(target[i] + ".md5").write_text(self._md5(target[i]))
						elif not f_src.exists():
							warn(f"Source file {str(f_src)} does not exist!")


##################################################################################################################################
# ----------------------------------------------- child class for mapping pipline ------------------------------------------------#
##################################################################################################################################

class MappingPipelinePathHandler(PipelinePathHandler):
	""" path handler for mapping pipeline """

	allowed_wildcards = ["step", "extension", "sample", "mate", "batch", "flowcell", "lane", "library", "lib_type"]
	required_wildcards_out_log = ["step", "extension", "sample"]
	required_wildcards_in = ["sample"]
	wildcard_fix_values = dict(sample="all_samples", mate="all_mates", batch="all_batches", flowcell="all_flowcells",
							   lane="all_lanes", library="all_libraries", lib_type="all_libtypes")

	def __init__(self, workflow, test_config=False, **kwargs):
		super().__init__(workflow, test_config, **kwargs)

		self.samples = self.snakemake_workflow.config["sample_info"]
		self.sample_ids = list(self.samples.keys())
		self._test_config_samples()

		# wildcard value-combinations parsed from input directory
		self.wildcard_combs = self._get_wildcard_combinations_per_sample(self.wildcard_values)

		# paths to static files
		self.data_paths = self.snakemake_workflow.config["organism"]

		# test if all sample IDs can be found in input path
		if test_config and not set(self.sample_ids).issubset(set(self.wildcard_values["sample"])):
			raise ValueError("Error in config file: not all specified samples could be found using given input path pattern")

		# write log with parsed values to file for debugging
		self._write_log()


	#---------------------------------------------------- helper methods ----------------------------------------------------#

	def _test_config_input(self, test_allowed_wildcards):
		super()._test_config_input(test_allowed_wildcards)

		# test whether wildcards used consistently across patterns
		if not set(self.out_path_wildcards) == set(self.log_path_wildcards) == set(self.in_path_wildcards) | set(
				self.required_wildcards_out_log) - set(self.required_wildcards_in):
			raise ValueError(
				"Error in config file: out_path_pattern, log_path_pattern and in_path_pattern "
				"do not contain the same wildcards. out: {}, log: {}, in: {}".format(
					set(self.out_path_wildcards),
					set(self.log_path_wildcards),
					set(self.in_path_wildcards)
				)
			)

		# test if all wildcards used in outdir path pattern
		if not set(self.out_path_wildcards) == set(self.outdir_path_wildcards) | {"extension"}:
			raise ValueError(
				"Error in config file: all wildcards of out and log dir should be used in folder names (exception: {{extension}}), otherwise "
				"different rules might compute output in the same folder, which can lead to mixed or deleted intermediate files. "
				"in folder name: {}, all: {}".format(set(self.outdir_path_wildcards), set(self.out_path_wildcards))
			)

	def _test_config_samples(self):
		# test if mate wildcard used if multiple read extensions are used
		if all([
				any(len(val["paired_end_extensions"]) > 1 for _, val in self.samples.items()),
				"mate" not in self.in_path_wildcards,
			]
		):
			raise ValueError(
				"Error in config file: some samples have more than one 'paired end extension' "
				"set in the sample info file, but {mate} wildcard in the in_path_pattern is not set. "
			)

	def _get_wildcard_combinations_per_sample(self, wildcard_values):
		""" go through wildcard values and get combinations per sample """

		wildcard_comb_num = len(wildcard_values["sample"])
		assert all([len(val) == wildcard_comb_num for val in wildcard_values.values()])
		WildcardComb = namedtuple("WildcardComb", [s for s in wildcard_values])
		per_sample_comb = {s: [] for s in wildcard_values["sample"]}
		for index in range(wildcard_comb_num):
			sample = wildcard_values["sample"][index]
			per_sample_comb[sample].append(
				WildcardComb(**{key: wildcard_values[key][index] for key in wildcard_values}))
		return per_sample_comb

	def _write_log(self, **kwargs):
		filename = self.file_path(step="MappingPipelinePathHandler", extension="log", log=True, fix="all", **kwargs)
		os.makedirs(os.path.dirname(filename), exist_ok=True)
		with open(filename, "w") as f:
			f.write(
				"pattern input:\n   in path pattern: {}\n   out path pattern: {}\n   log path pattern: {}\n\n"
				"parsed wildcards:\n   in path wildcards: {}\n   out/log path wildcards: {}\n   out dir wildcards: {}\n\n"
				"samples:\n   sample IDs: {}\n   sample info: {}\n\n"
				"wildcard values:\n   per wildcard: {}\n   combinations: {}\n\n"
				"neutral replacement of optional wildcards:\n   out/log: {}\n   out dir: {}".format(
					self.in_path_pattern, self.out_path_pattern, self.log_path_pattern,
					self.in_path_wildcards, self.out_path_wildcards, self.outdir_path_wildcards, self.sample_ids,
					self.samples, self.wildcard_values, self.wildcard_combs, self.opt_wildcard_placeholders,
					self.opt_wildc_placeh_outdir
				)
			)


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
			# TODO: case of ignored wildcards? meant for e.g. allFlowcell
			kwargs_filled = {key: getattr(comb, key) if ("{" in val or val not in self.wildcard_values[
				key]) and key != "mate" else val for key, val in kwargs_out.items()}
			kwargs_id_tup = tuple(kwargs_filled[key] for key in sorted(kwargs_filled))
			if kwargs_id_tup not in seen:
				seen.add(kwargs_id_tup)
				if mate_key:
					mate_lst = self.samples[wildcards.sample][mate_key]
					kwargs_filled["mate"] = mate_lst[mate] if len(mate_lst) > mate else mate_lst[0]
					pattern  = self.in_path_pattern.format(sample = wildcards.sample, **kwargs_filled)
					pattern += self.samples[wildcards.sample]["read_extension"]
				else:
					if "mate" not in kwargs_filled or mate == "*": kwargs_filled["mate"] = "*"
					pattern = self.in_path_pattern.format(sample=wildcards.sample, **kwargs_filled) + \
							  self.samples[wildcards.sample]["read_extension"]
				pattern_list.append(pattern)
		paths = [path for pat in pattern_list for path in iglob(pat, recursive=True)]
		paths.sort()
		return paths

	def file_path(self, step, extension, sample="{sample}", log=False, fix=None, path_pattern=None, **kwargs):
		"""
		Generate single path for intermediate and output or log files.
		
		:param step:  Snakemake rule for which the paths are generated
		:param extension: file extension for the generated file path
		:param sample: sample ID to be included in the file path
		:param log: generate path to logfile if log==True otherwise generate path to output/intermediate file
		:param fix: fix some wildcards; e.g. fix=["sample"] has the same effect as passing sample="all_samples", fix="all" fixes all wildcards
		:param path_pattern: explicitly set the path pattern (string)
		:param **kwargs: if used specify replacement for {batch}, {flowcell}, {lane}, etc. ...
		"""
		if fix:
			fixed_wildcards = self._get_wildcard_fix_values(fix)
			if sample == "{sample}" and "sample" in fixed_wildcards:
				sample = fixed_wildcards["sample"]
			kwargs = {**{k: v for k, v in fixed_wildcards.items() if k != "sample"}, **kwargs}
		if not path_pattern: path_pattern = self.log_path_pattern if log else self.out_path_pattern
		kwargs_out = {key: kwargs[key] if key in kwargs else val for key, val in self.opt_wildcard_placeholders.items()}
		return path_pattern.format(step=step, extension=extension, sample=sample, **kwargs_out)

	def out_dir_name(self, step, sample="{sample}", fix=None, **kwargs):
		"""
		Generate single path to intermediate and output file directory.
		
		:param step:  Snakemake rule for which the paths are generated
		:param sample: sample ID to be included in the file path
		:param fix: fix some wildcards; e.g. fix=["sample"] has the same effect as passing sample="all_samples", fix="all" fixes all wildcards
		"""
		if fix:
			fixed_wildcards = self._get_wildcard_fix_values(fix)
			if sample == "{sample}" and "sample" in fixed_wildcards:
				sample = fixed_wildcards["sample"]
			kwargs = {**{k: v for k, v in fixed_wildcards.items() if k != "sample"}, **kwargs}
		kwargs_out = {key: kwargs[key] if key in kwargs else val for key, val in self.opt_wildc_placeh_outdir.items()}
		return self.out_dir_pattern.format(step=step, sample=sample, **kwargs_out)

	def expand_path(self, step, extension="", fix=None, **kwargs):
		"""
		Generate multiple paths for intermediate and output files corresponding to different sample IDs (and e.g. paired end extensions).

		Always generates paths over all sample IDs and their combinations with keyword arguments (kwargs) that were not
		defined like, e.g. flowcell="all_flowcells", where "all_flowcells" is an arbitrary string that will be filled in the generated paths
		instead of the real flowcell IDs. Only keyword arguments that correspond to wildcards in the config path pattern (with the same name) are used.
		
		:param step:  Snakemake rule name (string) for which the paths are generated
		:param extension: file extension for the generated file path; '' to generate paths to directories
		:param fix: fix some wildcards; e.g. fix=["sample"] has the same effect as passing sample="all_samples", fix="all" fixes all wildcards
		:param **kwargs: replacement strings for optional wildcards, e.g. batch, flowcell, lane (see description)
		:returns: list of paths
		"""
		if fix:
			fixed_wildcards = self._get_wildcard_fix_values(fix)
			kwargs = {**{k: v for k, v in fixed_wildcards.items() if k != "sample"}, **kwargs}
		kwargs_out = {key: kwargs[key] if key in kwargs else val for key, val in self.opt_wildcard_placeholders.items()}
		paths = []
		for sample in self.sample_ids:
			seen = set()
			for comb in self.wildcard_combs[sample]:
				kwargs_filled = {key: getattr(comb, key) if "{" in val else val for key, val in kwargs_out.items()}
				kwargs_id_tup = tuple(kwargs_filled[key] for key in sorted(kwargs_filled))
				if kwargs_id_tup not in seen:
					seen.add(kwargs_id_tup)
					if extension:
						out_path_pattern = self.file_path(step, extension, sample=sample, **kwargs_filled)
					else:
						out_path_pattern = self.out_dir_name(step, sample=sample, **kwargs_filled)
					paths.append(out_path_pattern)
		return paths

	def link_index(self, step, sample="{sample}", entry="nopath", subdir="", add_done=False, fix=None, **kwargs):
		"""
		Generate symbolic link to index folder (if provided in config).
		
		A key-value pair is assumed to be in the organism config file with the name of the step (e.g. "star_index") as key 
		and a path to a corresponding index/input file. If no index should be linked, keep the key in the config file,
		but leave the value empty.
		
		:param step:  Snakemake rule for which the paths are generated
		:param sample: sample ID to be included in the file path
		:param entry: set the path to be linked explicitly
		:param fix: fix some wildcards; e.g. fix=["sample"] has the same effect as passing sample="all_samples", fix="all" fixes all wildcards
		:param **kwargs: if used specify replacement for {batch}, {flowcell}, {lane}, etc. ...
		"""
		if fix:
			fixed_wildcards = self._get_wildcard_fix_values(fix)
			if sample == "{sample}" and "sample" in fixed_wildcards:
				sample = fixed_wildcards["sample"]
			kwargs = {**{k: v for k, v in fixed_wildcards.items() if k != "sample"}, **kwargs}
		loc = Path(self.out_dir_name(step, sample, **kwargs)) / subdir
		if entry != "nopath":
			index = entry
		elif step in self.data_paths:
			index = self.data_paths[step]
		else:
			raise ValueError("Error linking index: no keyword provided in config to set index for {}!".format(step))
		if index:
			ind = Path(index)
			loc.parent.mkdir(exist_ok=True, parents=True)
			if not loc.exists():
				# create
				loc.symlink_to(ind.resolve())
				if add_done:
					Path(self.file_path(step, extension="done", sample=sample, **kwargs)).touch()
			elif not loc.samefile(index):
				# update
				loc.unlink()
				loc.symlink_to(ind.resolve())
				if add_done:
					Path(self.file_path(step, extension="done", sample=sample, **kwargs)).touch()

	def log_generated_files(self, save_to="", path_pattern=None, **kwargs):
		"""
		Scan all generated output/intermediate files and list them in a file with their corresponding wildcard values.
		This file can be used e.g. by the pipeline Rmd report to identify paths to files without the handler.
		"""
		table_cols = self._collect_generated_files(path_pattern=path_pattern)
		table = pd.DataFrame(table_cols)
		if not save_to:
			save_to = self.file_path(
				step="MappingPipelinePathHandler", extension="tsv", fix="all", path_pattern=path_pattern, **kwargs
			)
		table.to_csv(save_to, sep="\t", index=False)


##################################################################################################################################
# ------------------------------------------------- child class for DE pipline ---------------------------------------------------#
##################################################################################################################################

class DEPipelinePathHandler(PipelinePathHandler):
	""" path handler for differential expression pipeline """

	allowed_wildcards = ["step", "extension", "sample", "mate", "batch", "flowcell", "lane", "contrast", "mapping",
						 "library"]
	required_wildcards_out_log = ["step", "extension", "contrast"]
	required_wildcards_in = ["step", "extension", "sample"]
	wildcard_fix_values = dict(contrast="all")

	def __init__(self, workflow, test_config=False, **kwargs):
		super().__init__(workflow, test_config, **kwargs)

		self.contrasts = self.snakemake_workflow.config["contrasts"]["contrast_list"]
		self.contrast_defaults = self.snakemake_workflow.config["contrasts"]["defaults"]
		self.contrast_ids = [self._make_contrast_id(contr["title"], index) for index, contr in
							 enumerate(self.contrasts)]  # used in pipeline to fill in wildcards

		# edit the config dictionary
		for i, (c_dict, c_id) in enumerate(zip(self.contrasts, self.contrast_ids)):
			c_dict["ID"] = c_id
			c_merged = self.get_contrast(c_id)
			self.contrasts[i] = c_merged

		# write log with parsed values to file for debugging
		self._write_log(contrast="all")

	# ---------------------------------------------------- helper methods ----------------------------------------------------#

	def _get_filtered_wildcard_values(self, step, extension):
		""" return wildcard values filtered by step and extension """

		wildcard_values_filtered = {key: [] for key in self.wildcard_values}
		wildcard_comb_num = len(self.wildcard_values["sample"])
		assert all([len(val) == wildcard_comb_num for val in self.wildcard_values.values()])
		for index in range(wildcard_comb_num):
			if self.wildcard_values["step"][index] == step and self.wildcard_values["extension"][index] == extension:
				for key in self.wildcard_values: wildcard_values_filtered[key].append(self.wildcard_values[key][index])
		return wildcard_values_filtered

	def _write_log(self, **kwargs):
		filename = self.file_path(step="DEPipelinePathHandler", extension="log", log=True, **kwargs)
		os.makedirs(os.path.dirname(filename), exist_ok=True)
		with open(filename, "w") as f:
			f.write(
				"pattern input:\n   in path pattern: {}\n   out path pattern: {}\n   log path pattern: {}\n\n"
				"parsed wildcards:\n   in path wildcards: {}\n   out/log path wildcards: {}\n   out dir wildcards: {}\n\n"
				"wildcard values:\n   per wildcard: {}\n\n"
				"neutral replacement of optional wildcards:\n   out/log: {}\n   out dir: {}".format(
					self.in_path_pattern, self.out_path_pattern, self.log_path_pattern,
					self.in_path_wildcards, self.out_path_wildcards, self.outdir_path_wildcards, self.wildcard_values,
					self.opt_wildcard_placeholders, self.opt_wildc_placeh_outdir
				)
			)

	def _make_contrast_id(self, contrast_title, index):
		return self._make_name(contrast_title) + "_ID{}".format(index)

	@staticmethod
	def _check_dict_set(check_dict, set_val):
		try:
			if isinstance(set_val, dict):
				tmp_d = check_dict
				all_set = True
				for k, v in set_val.items():
					tmp_d = tmp_d[k]
					all_set = all_set and DEPipelinePathHandler._check_dict_set(tmp_d, v)
				return all_set
			else:
				return check_dict == set_val
		except KeyError:
			return False
		except:
			raise

	# -------------------------------------------- methods used in snakemake file --------------------------------------------#

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

	def get_contrast_id_dict(self, contrasts):
		""" return dict of contrast titles to IDs """
		titles = [contr["title"] for contr in contrasts]
		ids = []
		if len(set(titles)) < len(titles):
			warn("titles of contrasts defined in config are not unique! using last ID for each title")
		return {title: self._make_contrast_id(title, i) for i, title in enumerate(titles)}

	def file_path(self, step, extension, contrast="{contrast}", log=False, path_pattern=None, fix=None, **kwargs):
		"""
		Generate single path for intermediate and output or log files.
		
		:param step:  Snakemake rule for which the paths are generated
		:param extension: file extension for the generated file path
		:param contrast: contrast ID to be included in the file path
		:param log: generate path to logfile if log==True otherwise generate path to output/intermediate file
		:param path_pattern: explicitly set the path pattern (string)
		:param fix: fix some wildcards; e.g. fix=["sample"] has the same effect as passing sample="all_samples", fix="all" fixes all wildcards
		:param **kwargs: if used specify replacement for {batch}, {flowcell}, {lane}, etc. ...
		"""
		if fix:
			fixed_wildcards = self._get_wildcard_fix_values(fix)
			if contrast == "{contrast}" and "contrast" in fixed_wildcards:
				contrast = fixed_wildcards["contrast"]
			kwargs = {**{k: v for k, v in fixed_wildcards.items() if k != "contrast"}, **kwargs}
		if not path_pattern: path_pattern = self.log_path_pattern if log else self.out_path_pattern
		kwargs_out = {key: kwargs[key] if key in kwargs else val for key, val in self.opt_wildcard_placeholders.items()}
		return path_pattern.format(step=step, extension=extension, contrast=contrast, **kwargs_out)

	def expand_path(self, step, extension, if_set=False, fix=None, **kwargs):
		"""
		Generate multiple paths for intermediate and output files corresponding to different contrast IDs
		
		:param step:  Snakemake rule name (string) for which the paths are generated
		:param extension: file extension for the generated file path
		:param if_set: dict entry of a contrast used to filter for which to expand over; expand over all contrast if False
		:param fix: fix some wildcards; e.g. fix=["contrast"] has the same effect as passing contrast="all", fix="all" fixes all wildcards
		:param **kwargs: replacement strings for optional wildcards, e.g. batch, flowcell, lane (see description)
		:returns: list of paths
		"""
		if fix:
			fixed_wildcards = self._get_wildcard_fix_values(fix)
			kwargs = {**{k: v for k, v in fixed_wildcards.items() if k != "contrast"}, **kwargs}
		input_choices = set(self.input_choice)
		kwargs_out = {key: kwargs[key] if key in kwargs else val for key, val in self.opt_wildcard_placeholders.items()}
		expand_input_choices = list(input_choices & set(kwargs_out))
		choice_kwargs = [{expand_input_choices[i]: choice for i, choice in enumerate(choice_comb)}
						 for choice_comb in
						 itertools.product(*[self.input_choice[eic] for eic in expand_input_choices])]

		paths = []
		for contr in (c_id for i, c_id in enumerate(self.contrast_ids) if
					  any(self._check_dict_set({k: v}, if_set) for k, v in self.contrasts[i].items()) or not if_set):
			if not expand_input_choices:
				paths.append(self.file_path(step, extension, contrast=contr, **kwargs_out))
			else:
				kwargs_out = {key: val for key, val in kwargs_out.items() if key not in expand_input_choices}
				for kwargs_ch in choice_kwargs:
					paths.append(self.file_path(step, extension, contrast=contr, **kwargs_out, **kwargs_ch))
		return paths

	def log_generated_files(self, save_to="", path_pattern=None, **kwargs):
		"""
		Scan all generated output/intermediate files and list them in a file with their corresponding wildcard values.
		This file can be used e.g. by the pipeline Rmd report to identify paths to files without the handler.
		"""
		table_cols = self._collect_generated_files(path_pattern=path_pattern)
		table = pd.DataFrame(table_cols)
		if not save_to:
			save_to = self.file_path(step="DEPipelinePathHandler", extension="tsv", contrast="all",
									 path_pattern=path_pattern, **kwargs)
		table.to_csv(save_to, sep="\t", index=False)


##################################################################################################################################
# ----------------------------------------------------- covariate file tool ------------------------------------------------------#
##################################################################################################################################

class CovariateFileTool(PipelinePathHandler):
	""" Tool to generate a covariate file before running the pipeline; use same config as for DEPipelinePathHandler """

	allowed_wildcards = DEPipelinePathHandler.allowed_wildcards
	required_wildcards_out_log = DEPipelinePathHandler.required_wildcards_out_log
	required_wildcards_in = DEPipelinePathHandler.required_wildcards_in

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

		if "in_path_pattern" not in config_dict["pipeline_param"] or not config_dict["pipeline_param"]["in_path_pattern"]:
			raise ValueError("in_path_pattern needs to be set in the config file.")

		self.in_path_pattern = config_dict["pipeline_param"]["in_path_pattern"]

		self.wildcard_constraints = self._prepare_inpathpattern()

		self.wildcard_values = self._get_wildcard_values_from_input(self.in_path_pattern, verbose=True)

		if not self.wildcard_values["step"]:
			raise ValueError(
				"Error extracting wildcards: no wildcards found, because in_path_pattern did not match any files!\n"
				"in_path_pattern: {}".format(self.in_path_pattern)
				)

		self.opt_wildcard_placeholders = {w: "{{{}}}".format(w) for w in
										  set(self.wildcard_values) - set(self.required_wildcards_in)}

		self.covariate_data = pd.DataFrame({"filename": [], "md5": [], "group": [], "replicate": [], "label": []})

	# ---------------------------------------------------- helper methods ----------------------------------------------------#

	def _get_wildcard_combinations(self, wildcard_values, step, extension):
		""" go through wildcard values and get combinations """

		combinations = []
		WildcardComb = namedtuple("WildcardComb", [s for s in wildcard_values])
		wildcard_comb_num = len(wildcard_values["sample"])
		assert all([len(val) == wildcard_comb_num for val in wildcard_values.values()])
		for index in range(wildcard_comb_num):
			if wildcard_values["step"][index] == step and wildcard_values["extension"][index] == extension:
				combinations.append(WildcardComb(**{key: wildcard_values[key][index] for key in wildcard_values}))
		return combinations

	def _get_mapping_input(self, step, extension, wildcards, verbose=False):
		"""
		Generate paths to input files from read mapping, e.g. from STAR or Salmon.
		
		:param wildcards: rule wildcards
		:returns: a list with paths to specified input files
		"""
		wildcard_combs = self._get_wildcard_combinations(self.wildcard_values, step, extension)

		if not wildcard_combs:
			raise ValueError(f"No files found with combination of wildcards 'step': '{step}', 'extension': '{extension}'")

		if verbose:
			print(
				"\nextracted combinations:\n{}".format(
					"\n".join("\t".join(i) for i in [wildcard_combs[0]._fields] + wildcard_combs)
				)
			)

		wildcard_placeholders = {"sample": "{sample}", **self.opt_wildcard_placeholders}

		kwargs_out = {key: getattr(wildcards, key, val) for key, val in wildcard_placeholders.items()}
		pattern_list = []
		seen = set()
		for comb in wildcard_combs:
			kwargs_filled = {key: getattr(comb, key) if "{" in val or val not in self.wildcard_values[key] else val for
							 key, val in kwargs_out.items()}
			kwargs_id_tup = tuple(kwargs_filled[key] for key in sorted(kwargs_filled))
			if kwargs_id_tup not in seen:
				seen.add(kwargs_id_tup)
				pattern_list.append(self.in_path_pattern.format(step=step, extension=extension, **kwargs_filled))
		return [path for pat in pattern_list for path in iglob(pat, recursive=True)]


	#---------------------------------------------------- access methods ----------------------------------------------------#

	def update_covariate_data(self, step, extension, other = None):
		""" 
		fill mandatory columns of the covariate data frame by searching the input path specified in the config file
		
		:param step:  Snakemake rule name for which the files are searched
		:param extension: file extension of the searched files
		"""
		files = sorted(
			self._get_mapping_input(step, extension, SimpleNamespace(), verbose=True)
		)
		files_wildcards = {
			f: self._get_wildcard_values_from_file_path(f, input_pattern=self.in_path_pattern)[1]
			for f in files
		}

		# add md5, group, replicate, label per file
		md5   = [self._md5(f) for f in files]
		group = [files_wildcards[f]["sample"][0] for f in files]

		replicate, num_g = [], {}
		for g in group:
			if g not in num_g:
				num_g[g] = 1
			else:
				num_g[g] += 1
			replicate.append(num_g[g])
		label = ["{}_{}".format(a,b) for a,b in zip(group, replicate)]

		# add any extra files if given
		extra_files = {}
		if other:
			for col, (a_step, a_ext) in other.items():
				self._get_mapping_input(a_step, a_ext, SimpleNamespace(), verbose=True)  # for printing
				add_paths = []
				for f in files:
					wildcards = {
						k: v[0]
						for k, v in files_wildcards[f].items()
						if k not in ["step", "extension"]
					}
					add_p = self._get_mapping_input(a_step, a_ext, SimpleNamespace(**wildcards))
					if not len(add_p) == 1:
						raise ValueError(
							f"There should be only one file with these wildcards: \n{str(wildcards)}\n"
							f"found: \n{str(add_p)}"
						)
					add_paths.append(add_p[0])
				extra_files[col] = add_paths

		# create data frame
		self.covariate_data = pd.DataFrame(
			{"filename":files, "md5":md5, "group":group, "replicate":replicate, "label":label, **extra_files}
		)

	def add_column(self, name, levels):
		"""
		add a custom column to the covariate data frame
		
		levels can be either a list (order important!) or a dictionary. If levels is a dictionary it can be of two forms:
		- {<group(string)>: <level(string)>, ...}
		- {<level(string)>: [<group(string)>, <group(string)>, ...], ...}
		
		:param name:  name of the column
		:param levels: levels of the column
		"""
		if type(levels) == dict:
			if any(l not in self.covariate_data.group for l in levels):
				self.covariate_data[name] = [l for g in self.covariate_data.group for l, gs in levels.items() if
											 g in gs]
			else:
				self.covariate_data[name] = [levels[g] for g in self.covariate_data.group]
		elif type(levels) == list:
			self.covariate_data[name] = levels

	def write_covariate_file(self, filename):
		""" write covariate data frame to file """
		self.covariate_data.to_csv(filename, sep="\t", index=False)


##################################################################################################################################
# ------------------------------------------------------ sample info tool --------------------------------------------------------#
##################################################################################################################################


class SampleInfoTool(PipelinePathHandler):
	""" Tool to generate a sample info file before running the pipeline; use same config as for MappingPipelinePathHandler """
	allowed_wildcards          = MappingPipelinePathHandler.allowed_wildcards + ["lib_type"]
	required_wildcards_out_log = MappingPipelinePathHandler.required_wildcards_out_log
	required_wildcards_in      = MappingPipelinePathHandler.required_wildcards_in
	wildcard_fix_values        = MappingPipelinePathHandler.wildcard_fix_values

	allowed_read_extensions    = [".fastq", ".fastq.gz", ".fq", ".fq.gz"]

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

		self.in_path_pattern = config_dict["pipeline_param"]["in_path_pattern"]

		self.wildcard_constraints = self._prepare_inpathpattern()

		self.sample_info = {}

	# ---------------------------------------------------- helper methods ----------------------------------------------------#

	def _get_wildcard_values_from_read_input(self, unix_style=True):
		""" go through files in input path and get values matching the wildcards """

		glob_pattern = re.sub("{[^}./]+}", "*", self.in_path_pattern)
		wildcards = re.findall("{([^}./]+)}", self.in_path_pattern)
		match_pattern = re.sub("\\\\{([^}./]+)\\\\}", self._wildc_replace, re.escape(self.in_path_pattern))
		input_files = glob(glob_pattern + ("*" if glob_pattern[-1] != "*" else ""), recursive=True)
		if unix_style:
			match_pattern = re.sub(r"\\\*\\\*", "[^{}]*", match_pattern)
			match_pattern = re.sub(r"(?<!\[\^{}\]\*)\\\*", "[^{}./]*", match_pattern)

		print("\ninput files:\n{}".format("\n".join(input_files)))
		print(f"\nmatch pattern:\n{match_pattern}")

		wildcard_values = {w: [] for w in wildcards}
		for inp in input_files:
			self._get_wildcard_values_from_file_path(
				inp, input_pattern=self.in_path_pattern, wildc_val=wildcard_values, unix_style=unix_style
			)

		return {
			**wildcard_values,
			"read_extension": [
				f.replace(re.match(match_pattern, f).group(0), "")
				for f in input_files
				if re.match(match_pattern, f)
			]
		}

	def _convert_str_entries_to_lists(self, key="paired_end_extensions"):
		""" for importing lists from table entries """
		for smpl_info in self.sample_info.values():
			smpl_info[key] = [s.replace("'", "").replace('"', "") for s in re.findall("[^\[\]\s,]+", smpl_info[key])]

	def _add_info_fields(self, add_dict):
		""" add fields from add_dict to self.sample_info if they are not already present """
		for sample, fields in add_dict.items():
			if sample in self.sample_info:
				s_info = self.sample_info[sample]
				for f_key, f_val in fields.items():
					if f_key not in s_info: s_info[f_key] = f_val

	# ---------------------------------------------------- access methods ----------------------------------------------------#

	def update_sample_info(self, library_default="unstranded", add=False):
		"""
		fill mandatory info about sample by searching the input path specified in the config file.
		
		attention: stranded is initially set to library_default for all samples! 
		This information has to be edited manually in table or yaml, if libraries were prepared differently.
		
		:param library_default:  options: ["unstranded", "forward", "reverse"]
		"""
		wildcard_values = self._get_wildcard_values_from_read_input()
		if not wildcard_values["sample"]:
			raise ValueError(
				"Error extracting wildcards: in_path_pattern did not match any file path!\n"
				f"in_path_pattern: {self.in_path_pattern}\n"
			)
		wildcard_combs = [comb for comb in self._get_wildcard_combinations(wildcard_values) if comb.read_extension in self.allowed_read_extensions]
		if not wildcard_combs:
			raise ValueError(
				"Error extracting wildcards: read extension not matched!\n"
				f"list of allowed extensions: {self.allowed_read_extensions}"
			)
		print(
			"\nextracted combinations:\n{}".format(
				"\n".join(
					"\t".join(i)
					for i in [wildcard_combs[0]._fields] + wildcard_combs
				)
			)
		)

		sample_info = {}
		for comb in wildcard_combs:
			if comb.sample not in sample_info:
				sample_info[comb.sample] = {"stranded": library_default, "read_extension": comb.read_extension}
				sample_info[comb.sample]["paired_end_extensions"] = [getattr(comb, "mate", "")]
				sample_info[comb.sample]["lib_types"] = {comb.lib_type: None} if hasattr(comb, "lib_type") else {}
			else:
				if hasattr(comb, "mate"):
					paired_end_ext = getattr(comb, "mate", "")
					paired_end_ext_lst = sample_info[comb.sample]["paired_end_extensions"]
					if paired_end_ext_lst == [""]:
						raise ValueError(
							"Error compiling sample information: sample {} has names with and without paired end extensions".format(
								comb.sample))
					if paired_end_ext not in paired_end_ext_lst:
						paired_end_ext_lst.append(paired_end_ext)
						paired_end_ext_lst.sort()
				if hasattr(comb, "lib_type"):
					lib_types = sample_info[comb.sample]["lib_types"]
					if lib_types == {}:
						raise ValueError(
							"Error compiling sample information: sample {} has names with and without library type information".format(
								comb.lib_type))
					if comb.lib_type not in lib_types:
						lib_types[comb.lib_type] = None
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
		print(f"isatab parse config file: {str(parse_conf)}")
		with open(str(parse_conf), "r") as stream:
			try:
				parse_conf = yaml.safe_load(stream)
			except yaml.YAMLError as exc:
				print(exc)

		def find_columns(key):
			if not parse_conf[key]["columns"]: return None
			cols = []
			for col_regex in parse_conf[key]["columns"]:
				for col in tab.columns:
					if re.match(col_regex, col):
						cols.append(col)
			return cols

		def map_value(key, value, col):
			if "as_is" in parse_conf[key] and parse_conf[key]["as_is"]:
				return value
			value = str(value)
			for val_regex, val_repl in parse_conf[key]["value"].items():
				if re.fullmatch(val_regex, value):
					if isinstance(val_repl, str):
						return re.sub(val_regex, val_repl, value)
					elif isinstance(val_repl, dict):
						repl_dict = {}
						for k_repl, v_repl in val_repl.items():
							add_key = re.sub(val_regex, k_repl, value)
							if add_key == add_key and add_key != 'nan':
								col_mode = False
								if r"\col_mode: " in v_repl:
									v_repl = v_repl.replace(r"\col_mode: ", "")
									col_mode = True

								if col_mode:
									opts = {k: v for pair in v_repl.split(",") for k, v in [pair.split(":")]}
									repl_dict.update({add_key: opts[col]})
								else:
									repl_dict.update({add_key: re.sub(val_regex, v_repl, value)})
						return list(repl_dict.items())
					else:
						return val_repl

		sample_info_cols = {}
		for key in parse_conf:
			add = parse_conf[key]["add"] if "add" in parse_conf[key] else False
			columns = find_columns(key)
			print(f"key: {key} | column: {str(columns)}")
			if columns is not None:
				for col in columns:
					if key not in sample_info_cols:
						sample_info_cols[key] = [map_value(key, val, col) for val in tab[col]]
					elif add:
						for sic, val in zip(sample_info_cols[key], tab[col]):
							if isinstance(sic, list):
								sic.extend(map_value(key, val, col))
							else:
								raise ValueError(
									f"Error in ISAtab_parse_conf.yaml: cannot 'add', because {sic} is not a list!"
									)

		self.sample_info = {
			sample_id: {
				key: val[i]
				for key, val in sample_info_cols.items()
				if key != "id"
				}
			for i, sample_id in enumerate(sample_info_cols["id"])
		}


##################################################################################################################################
# ------------------------------------------------------ class for report --------------------------------------------------------#
##################################################################################################################################

class ReportTool(PipelinePathHandler):

	entry_heading_pattern = "#HSTART.*\n([\s\S]*{{ENTRY_NAME}}[\s\S]*)\n.*#HEND"
	insert_pattern        = "#>.*INSERT.*<#"
	entry_name_wildcard   = "{{ENTRY_NAME}}", "{{ENTRY_ID}}"

	def __init__(self, pph, profile="DE"):
		self.path_handler = pph
		config_dict = self.path_handler.snakemake_workflow.config

		if config_dict["pipeline_param"]["report_snippets"]:
			self.report_snippet_base_dir = Path(config_dict["pipeline_param"]["report_snippets"])
		else:
			self.report_snippet_base_dir = Path(sys.path[0]) / "report"

		self.use_results = self._make_use_results_dict(config_dict)
		self.merge_mode  = bool(config_dict["report"]["merge"]) if "merge" in config_dict["report"] else False
		self.config = config_dict

		# define substitutions for report generation
		# self.substitutions = {"__contrasts__": self.path_handler.get_contrast_id_dict}

		# define profiles (different kinds of reports)
		if profile == "DE":
			# main report of DE pipeline
			self.start_template = "report_main_template.Rmd"
			self.report_snippet_base_dir = Path(sys.path[0]) / "report" / "Rmd" / "DE_report"
			self.report_snippet_building_plan = config_dict["report"]["report_snippets"]
			self.report_snippet_defaults = config_dict["report"]["defaults"]
			self.id_dict_for_analysis = self._DE_id_dict_from_path
			self.req_fields = ["step", "extension", "contrast"]
		elif profile == "circRNA":
			# report for circRNA analysis
			self.start_template = "circRNA_report_main_template.Rmd"
			self.report_snippet_base_dir = Path(sys.path[0]) / "report" / "Rmd" / "circRNA_analysis"
			self.report_snippet_building_plan = config_dict["circRNA_report"]["report_snippets"]
			self.report_snippet_defaults = config_dict["circRNA_report"]["defaults"]
			self.id_dict_for_analysis = self._mapping_id_dict_from_path
			self.req_fields = ["step", "extension", "sample"]
		elif profile == "sc_analysis":
			# report for circRNA analysis
			self.start_template = "sc_analysis_main_template.ipynb"
			self.report_snippet_base_dir = Path(sys.path[0]) / "report" / "ipynb" / "sc_analysis"
			self.report_snippet_building_plan = config_dict["jupyter_notebook"]["report_snippets"]
			self.report_snippet_defaults = config_dict["jupyter_notebook"]["defaults"]
			self.id_dict_for_analysis = self._mapping_id_dict_from_path
			self.req_fields = ["step", "extension", "sample"]

		if config_dict["pipeline_param"]["report_snippets"]:
			self.report_snippet_base_dir = Path(config_dict["pipeline_param"]["report_snippets"])

		self.snippet_path = [ self.report_snippet_base_dir ]
		if config_dict["report"]["path"]:
			self.snippet_path=[ Path(p) for p in config_dict["report"]["path"].split(os.pathsep) ] + self.snippet_path 


		self._id_cache = {}

	# ---------------------------------------------------- helper methods ----------------------------------------------------#

	def _make_use_results_dict(self, config_dict):
		""" unify different definitions of the merge specification; it can be str, list, dict """
		merge_spec = config_dict["report"]["merge"] if "merge" in config_dict["report"] else None
		if isinstance(merge_spec, str): merge_spec = [merge_spec]
		if isinstance(merge_spec, list): merge_spec = {"analysis": [path for path in merge_spec]}
		if isinstance(merge_spec, dict):
			for k, v in merge_spec.items():
				if isinstance(v, str): merge_spec[k] = [v]
			if "analysis" not in merge_spec: merge_spec["analysis"] = []
			merge_spec["analysis"].insert(0, self.path_handler.out_path_pattern)
		else:
			merge_spec = {"analysis": [self.path_handler.out_path_pattern]}
		return merge_spec

	def _split_template(self, template_text):
		before, after = re.split(self.insert_pattern, template_text, maxsplit=1, flags=re.MULTILINE)
		return (before, after)

	def _DE_id_dict_from_path(self, path):
		"""
		return a dict of contrast title to ID,
		constructed from the config file corresponding to a specific analysis
		which is described by the respective path pattern
		"""
		if path not in self._id_cache:
			config_dict = self.load_config_from_path(path, self.path_handler)
			id_dict = dict(contrast=self.path_handler.get_contrast_id_dict(config_dict["contrasts"]["contrast_list"]))
			self._id_cache[path] = id_dict
		return self._id_cache[path]

	def _mapping_id_dict_from_path(self, path):
		"""
		return a dict of sample name to ID (currently identical)
		that is constructed from the config file (corresponding to a specific analysis)
		which is described by the respective path pattern
		"""
		if path not in self._id_cache:
			config_file = self.path_handler.file_path("pipeline_report", "yaml", fix="all", path_pattern=path)
			with open(config_file, "r") as stream:
				try:
					config_dict = yaml.safe_load(stream)
				except yaml.YAMLError as exc:
					print(exc)
			id_dict = dict(sample={s: s for s in list(config_dict["sample_info"])})
			self._id_cache[path] = id_dict
		return self._id_cache[path]

	def _make_id(self, name, results_path):
		entry_name, templ_name = name
		id_dicts = self.id_dict_for_analysis(results_path)
		if templ_name in id_dicts:
			return id_dicts[templ_name][entry_name]
		else:
			return self._make_name(entry_name)

	def _insert_entry_name(self, text, name, results_path):
		entry_name = name[0]
		return text.replace(
			self.entry_name_wildcard[0], entry_name
		).replace(
			self.entry_name_wildcard[1], self._make_id(name, results_path)
		)

	def _get_entry_heading_code(self, template_text):
		return re.search(self.entry_heading_pattern, template_text).group(1)

	def _get_entry_heading(self, heading_code, name):
		return heading_code.replace(self.entry_name_wildcard[0], name)

	def _rem_entry_heading_code(self, template_text):
		return re.sub(self.entry_heading_pattern, "", template_text)

	def _insert_file_paths(self, text, path_pattern):
		wildcards_orig = re.findall("{{([^{}]+)}}", text)
		wildcards_prep = [rwo.replace(".", "_") for rwo in wildcards_orig]
		for rwo, rwp in zip(wildcards_orig, wildcards_prep):
			text = text.replace(rwo, rwp)
			parts = rwo.split("-")
			if len(parts) == 2:
				text = text.replace("{{" + rwp + "}}", self.path_handler.file_path(parts[0], parts[1], contrast="all",
																				   path_pattern=path_pattern))
			elif len(parts) == 3:
				text = text.replace("{{" + rwp + "}}",
									self.path_handler.file_path(parts[0], parts[1], contrast=parts[2],
																path_pattern=path_pattern))
		return text

	def _get_entry_list_from_str(self, entries, results_path):
		""" expand predefined placeholders to entry lists """
		if entries == "__contrasts__":
			return list(self._DE_id_dict_from_path(results_path)["contrast"])
		elif entries == "__timeseries__":
			return list(self.config["time_series"])
		elif entries == "__timeseries_comp__":
			return list(self.config["dodr"]["comparisons"])
		elif entries == "__samples__":
			return list(self._mapping_id_dict_from_path(results_path)["sample"])
		else:
			return [entries]

	def get_id_suffix(self, tag, num):
		if len(self.use_results[tag]) == 1 and tag != "analysis": num = ""
		ana_id_var = "_{}{}".format(tag, num) if self.merge_mode else ""
		ana_id_tit = " -- {}{}".format(tag, num) if self.merge_mode else ""
		return (ana_id_var, ana_id_tit)

	def _edit_template(self, text, path_pattern, tag, num):
		ana_id_var, ana_id_tit = self.get_id_suffix(tag, num)
		text = text.replace("file_tab", "file_tab{}".format(ana_id_var))
		text = text.replace("config", "config{}".format(ana_id_var))
		text = re.sub("#+ [^{\n]+", r"\g<0>{} ".format(ana_id_tit), text, count=1)

		def subs(m):
			return f"{m.group(0)}" if "=" in m.group(0) else f"{m.group(0)}_{ana_id_var}"
		text = re.sub("{ ?r +[^,\n ]+(?! *\n)", subs, text)
		return self._insert_file_paths(text, path_pattern)

	# -------------------------------------------------- recursive assembly --------------------------------------------------#

	def _assemble_entries(self, entries, path, snippet_name, entry_heading_code, results_key):
		""" assemble entry list (e.g. list of contrasts) """

		entry_text = []

		# add subsection entries
		results_path = self.use_results[results_key[1]][results_key[0]]
		if type(entries) is str: entries = self._get_entry_list_from_str(entries, results_path)
		for entry in entries:
			if type(entry) is str:

				entry_name = entry
				snippet_list_defaults = snippet_name + "_list"
				if snippet_list_defaults not in self.report_snippet_defaults:
					raise KeyError("Error compiling report snippets for {snip} {entr}! (no snippets provided and "
								   "key {snip} not found in config defaults)".format(snip=snippet_list_defaults,
																					 entr=entry))
				sub_snippet_list = self.report_snippet_defaults[snippet_list_defaults]

			elif isinstance(entry, dict):

				assert len(entry)==1
				entry_name = list(entry.keys())[0]

				sub_snippet_list  = entry[entry_name]

			else:
				raise TypeError(
					"Error in report snippet building plan! (expected str or dict, got {})".format(type(entry))
				)

			entry_text += [
				self._get_entry_heading(entry_heading_code, entry_name),
				"\n\n" ,
				self._assemble_template(
					sub_snippet_list, path, snippet_name, entry_heading_code, (entry_name, snippet_name), results_key
				)
			]

		return "".join(entry_text)


	def _search_snippet_path(self, snippet, searchpath):
		""" given a search path, find where the snippet is, return the full path """
		for path in searchpath:
			pp = path / snippet
			if pp.exists():
				return pp
		raise KeyError("Snippet '{}' not found in path '{}'".format(snippet, ':'.join([ str(p) for p in searchpath ])))
		

	def _assemble_template(self, snippet_list, path, snippet_name, entry_heading_code, entry=("", ""),
						   results_key=(-1, "analysis")):
		""" assemble snippet list """

		snippet_text = []

		if type(snippet_list) is str:
			snippet_list = (
				self.report_snippet_defaults[snippet_name]
				if snippet_list=="__defaults__"
				else [snippet_list]
			)
		for snippet in snippet_list:
			# --- add snippets
			if type(snippet) is str:

				#snippet_file = path / snippet
				#snippet_cont = snippet_file.read_text()
				snippet_cont = self._search_snippet_path(snippet, path).read_text()

				for i, results_path in (tup for tup in enumerate(self.use_results[results_key[1]]) if tup[0]==results_key[0] or results_key[0]<0):
					snippet_prep = self._insert_entry_name(snippet_cont, entry, results_path)
					requirements = re.findall("(?<=#REQUIRE)\s+{{(\S+)}}", snippet_prep)
					snippet_prep = re.sub(    "#REQUIRE\s+{{\S+}}\n+", "", snippet_prep)

					if all(Path(self.path_handler.file_path( **dict(zip( self.req_fields, req.split("-") ), path_pattern=results_path) )).exists() for req in requirements):
						snippet_text.append(self._edit_template(snippet_prep, results_path, results_key[1], i))

			elif isinstance(snippet, dict):
				assert len(snippet) == 1
				snippet_key = list(snippet.keys())[0]
				snippet_value = list(snippet.values())[0]

				# --- add list of entries
				if snippet_key == "__list__":
					for i, results_path in (tup for tup in enumerate(self.use_results[results_key[1]]) if
											tup[0] == results_key[0] or results_key[0] < 0):
						snippet_text.append(
							self._assemble_entries(snippet_value, path, snippet_name, entry_heading_code,
												   results_key=(i, results_key[1])))
				# --- change results folder
				elif snippet_key in self.use_results:
					if not isinstance(snippet_value, list): snippet_value = [snippet_value]
					add_txt = self._assemble_template(snippet_value, path, snippet_name, entry_heading_code, entry,
													  results_key=(-1, snippet_key))
					snippet_text.append(add_txt)
				# --- open sub-template folder
				else:
					sub_template_path = self._search_snippet_path(Path(snippet_key) / (snippet_key + "_main_template.Rmd"), path)
					# sub_template_path  = path / snippet_key / (snippet_key + "_main_template.Rmd")
					sub_template_text  = sub_template_path.read_text()

					entry_heading_code = self._get_entry_heading_code(sub_template_text)
					sub_template_text  = self._rem_entry_heading_code(sub_template_text)

					# add sub-section
					for i, results_path in (
						tup
						for tup in enumerate(self.use_results[results_key[1]])
						if tup[0]==results_key[0] or results_key[0]<0
					):
						temp_begin, temp_end = self._split_template(
							self._edit_template(sub_template_text, results_path, results_key[1], i)
						)
						sub_path = [ self._search_snippet_path(snippet_key, path) ]
						sub_section_text = self._assemble_template(
							snippet_value, sub_path, snippet_key, entry_heading_code, entry, results_key=(i,results_key[1])
						)
						all_text = temp_begin + sub_section_text + temp_end

						all_text_prep = self._insert_entry_name(all_text, entry, results_path)
						requirements  = re.findall("(?<=#REQUIRE)\s+{{(\S+)}}", all_text_prep)
						all_text_prep = re.sub(    "#REQUIRE\s+{{\S+}}\n+", "", all_text_prep)

						if all(
							Path(
								self.path_handler.file_path(
									**dict(zip(self.req_fields, req.split("-")), path_pattern=results_path)
								)
							).exists()
							for req in requirements
						):
							snippet_text.append(all_text_prep)
			else:
				raise TypeError(
					"Error in report snippet building plan! (expected str or dict, got {})".format(type(snippet))
				)

		return "\n".join(snippet_text)

	# ---------------------------------------------------- access methods ----------------------------------------------------#

	def generate_report(self):
		"""
		Generate a report by assembling Rmd snippets as defined in the config.
		
		Starts in report directory and recursively alternates between concatenating lists of snippets 
		and creating sub-lists of entries (e.g. contrasts).
		"""
		template_path = self._search_snippet_path(self.start_template, self.snippet_path)
		template_text = template_path.read_text()

		# generate report
		temp_begin, temp_end = self._split_template(template_text)
		report_text = self._assemble_template(
			self.report_snippet_building_plan,
			path=self.snippet_path,
			snippet_name="other",
			entry_heading_code="# {{ENTRY_NAME}}"
		)

		return temp_begin + report_text + temp_end
