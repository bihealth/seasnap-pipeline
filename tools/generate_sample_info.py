if __name__ == '__main__':
	import sys
	from pipeline_tools import SampleInfoTool

	config_files = sys.argv[1:]


	sit = SampleInfoTool(*config_files)

	# update info
	sit.update_sample_info(library_default="unstranded") # <-- EDIT file after generation to replace default

	# write to file
	sit.write_yaml("sample_info.yaml")


	# # additional read/write options:
	# sit.read_yaml("sample_info.yaml")
	# sit.write_table("sample_info.tsv")
	# sit.read_table("sample_info.tsv")
