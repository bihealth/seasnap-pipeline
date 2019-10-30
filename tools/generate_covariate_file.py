if __name__ == '__main__':
	import sys
	from pipeline_tools import CovariateFileTool

	# e.g. step=salmon extension=sf
	step, extension = sys.argv[1:3]
	config_files    = sys.argv[3:]


	cft = CovariateFileTool(*config_files)

	# fill 5 mandatory columns
	cft.update_covariate_data(step, extension)

	# add custom columns <-- EDIT HERE
	cft.add_column("condition", {"classical": ["SRR4053824","SRR4053795"], "nonclassical": ["SRR4053802","SRR4053812"]})

	# write to file
	cft.write_covariate_file("covariate_file.txt")
