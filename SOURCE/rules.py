wildcard_constraints:
	protein_name = config["protein_name"]


rule parse_input:
	input:
		sample_list = "INPUT/sample_list.csv"
	output:
		numQueries = "OUTPUT/tmp/numQueries.RData",
		rawNames = "OUTPUT/tmp/rawNames.RData"
	params:
		protein_name = config["protein_name"]
	conda:
		"dependencies.yaml"
	script:
		"1_parseInput.R"


rule split_quantitation:
	input:
		sample_list = "INPUT/sample_list.csv",
		rawNames = "OUTPUT/tmp/rawNames.RData"
	output:
		split_quantitation = expand("OUTPUT/tmp/{protein_name}/split_quantitation.txt",
			protein_name=config["protein_name"])
	params:
		protein_name = config["protein_name"]
	conda:
		"dependencies.yaml"
	script:
		"2_splitReplicates_quantitation.R"


rule extract_peptides:
	input:
		sample_list = "INPUT/sample_list.csv",
		rawNames = "OUTPUT/tmp/rawNames.RData"
	output:
		extract_peptides = expand("OUTPUT/tmp/{protein_name}/extract_peptides.txt",
			protein_name=config["protein_name"])
	params:
		protein_name = config["protein_name"],
		KK = config["KK"],
		IONscore = config["IONscore"],
		thrDif = config["thrDif"],
		thrDifpcp = config["thrDifpcp"]
	conda:
		"dependencies.yaml"
	script:
		"3_extractPeptides.R"


rule reassign_queries:
	input:
		sample_list = "INPUT/sample_list.csv",
		numQueries = "OUTPUT/tmp/numQueries.RData",
		extract_peptides = expand("OUTPUT/tmp/{protein_name}/extract_peptides.txt",
			protein_name=config["protein_name"]),
		split_quantitation = expand("OUTPUT/tmp/{protein_name}/split_quantitation.txt",
			protein_name=config["protein_name"])
	output:
		reassign_queries = "OUTPUT/tmp/reassign_queries.txt"
	params:
		protein_name = config["protein_name"]
	conda:
		"dependencies.yaml"
	script:
		"4_reassignQuantitationQuerries.R"


rule merge_replicates:
	input:
		sample_list = "INPUT/sample_list.csv",
		reassign_queries = "OUTPUT/tmp/reassign_queries.txt"
	output:
		mapping = expand("OUTPUT/{protein_name}/mapping.RData",
			protein_name=config["protein_name"]),
		merge_replicates = "OUTPUT/tmp/merge_replicates.txt"
	params:
		protein_name = config["protein_name"],
		PCPthresh = config["PCPthresh"],
		PSPthresh = config["PSPthresh"]
	conda:
		"dependencies.yaml"
	script:
		"5_mergeReassignedReplicates.R"


rule parse_intensities:
	input:
		sample_list = "INPUT/sample_list.csv",
		mapping = expand("OUTPUT/{protein_name}/mapping.RData",
			protein_name=config["protein_name"])
	output:
		quantity = expand("OUTPUT/{protein_name}/quantity.RData",
			protein_name=config["protein_name"]),
		quantity_df = expand("OUTPUT/{protein_name}/quantity.csv",
			protein_name=config["protein_name"])
	params:
		protein_name = config["protein_name"]
	conda:
		"dependencies.yaml"
	script:
		"6_parseIntensitiesForallReplicates.R"


rule collect_replicates:
	input:
		sample_list = "INPUT/sample_list.csv",
		quantity = expand("OUTPUT/{protein_name}/quantity.RData",
			protein_name=config["protein_name"])
	output:
		quantity_rep = expand("OUTPUT/{protein_name}/quantity_rep.RData",
			protein_name=config["protein_name"])
	params:
		protein_name = config["protein_name"]
	conda:
		"dependencies.yaml"
	script:
		"7_collectReplicates.R"


rule process_kinetics:
	input:
		sample_list = "INPUT/sample_list.csv",
		quantity_rep = expand("OUTPUT/{protein_name}/quantity_rep.RData",
			protein_name=config["protein_name"])
	output:
		results = expand("OUTPUT/{protein_name}/results.RData",
			protein_name=config["protein_name"])
	params:
		protein_name = config["protein_name"]
	conda:
		"dependencies.yaml"
	script:
		"8_postProcessKinetics.R"


rule filter_kinetics:
	input:
		sample_list = "INPUT/sample_list.csv",
		results = expand("OUTPUT/{protein_name}/results.RData",
			protein_name=config["protein_name"])
	output:
		filteredResults = expand("OUTPUT/{protein_name}/filteredResults.RData",
			protein_name=config["protein_name"]),
		filteredMeans = expand("OUTPUT/{protein_name}/filteredMeans.RData",
			protein_name=config["protein_name"])
	params:
		protein_name = config["protein_name"]
	conda:
		"dependencies.yaml"
	script:
		"9_filterKinetics.R"


rule complete_assignment:
	input:
		sample_list = "INPUT/sample_list.csv",
		filteredResults = expand("OUTPUT/{protein_name}/filteredResults.RData",
			protein_name=config["protein_name"]),
		filteredMeans = expand("OUTPUT/{protein_name}/filteredMeans.RData",
			protein_name=config["protein_name"])
	output:
		finalfilteredResults = expand("OUTPUT/{protein_name}/filteredResults_final.RData",
			protein_name=config["protein_name"]),
		finalfilteredMeans = expand("OUTPUT/{protein_name}/filteredMeans_final.RData",
			protein_name=config["protein_name"]),
		KineticsDB = "OUTPUT/KineticsDB.csv"
	params:
		protein_name = config["protein_name"],
		rm_tmp = config["rm_tmp"]
	conda:
		"dependencies.yaml"
	script:
		"10_getCompleteTypeAssignment.R"


