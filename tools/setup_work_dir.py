if __name__ == '__main__':
	import sys, time, shutil
	from pathlib import Path

	## run: setup_work_dir.py [<dir_name> <config_file1> <config_file2> ...]

	############# PATHS OF CONFIG FILES
	configs = dict(DE      = "DE_config.yaml",
		       mapping = "mapping_config.yaml")

	choice = sys.argv[2:] if len(sys.argv)>2 else configs.keys()

	print("available config files:")
	for key, val in configs.items(): print("\t{}\t{}".format(key,val))
	print("copy: {}".format(", ".join(choice)))
	print("-----------------------------------")
	print("starting...")


	config_files = [Path(sys.path[0]) / ".." / val for key,val in configs.items() if key in choice]


	############# CREATE WORKING DIRECTORY
	if len(sys.argv)>1:
		fstr = sys.argv[1]
	else:
		fstr = "results_%Y_%m_%d_%X/"


	working_dir = Path(time.strftime(fstr))

	try:
		working_dir.mkdir(parents=True)
		print("working directory {} created...".format(str(working_dir)))
	except FileExistsError:
		print("...directory {} already exists".format(str(working_dir)))
		print("...stopping.")
		raise 


	############# COPY CONFIG FILES
	for configf in config_files:
		shutil.copy(str(configf), str(working_dir / configf.name))
		print("config file {} copied...".format(str(configf.name)))
	print("...finished.")


	############# SYMLINK TO WRAPPER
	Path("sea-snap").symlink_to(Path(sys.path[0]) / "sea-snap.py")
