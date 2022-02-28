import numpy as np
import pandas as pd
import igraph as ig
import pyreadr
from genetic_algorithm.coverage_evaluation import coverage_sum

class EvaluationFunction(object):
	_instance = None

	def __new__(
			cls,
			smiles_distance_path,
			moa_distance_path,
			graph_distance_path,
			ppi_path,
			graph_rank_path,
			drug_targets_path
		):
		# Singleton design pattern
		if cls._instance is None:
			cls._instance = super(EvaluationFunction, cls).__new__(cls)
			cls._instance.smiles_distances = pyreadr.read_r(smiles_distance_path)[None]
			cls._instance.moa_distances = pyreadr.read_r(moa_distance_path)[None]
			cls._instance.paths_distances = pyreadr.read_r(graph_distance_path)[None]
			cls._instance.ppi_network = ig.read(ppi_path)
			cls._instance.graph_rank = pyreadr.read_r(graph_rank_path)[None]
			cls._instance.L1000_drug_targets = pd.read_csv(drug_targets_path, sep="\t")
			cls._instance.drug_names = cls._instance.smiles_distances.index.to_list()
		return cls._instance

	def __call__(self, individual):
		if np.sum(individual) <= 1:
			return (0, 0, 0, 1, len(self.drug_names))
		candidate_drugs = [self.drug_names[i] for i, bit in enumerate(individual) if bit == 1]
		smiles_submat = self.smiles_distances.loc[candidate_drugs, candidate_drugs]
		moa_submat = self.moa_distances.loc[candidate_drugs, candidate_drugs]
		paths_submat = self.paths_distances.loc[candidate_drugs, candidate_drugs]
		_, coverage_pval = coverage_sum(
			candidate_drugs,
			self.L1000_drug_targets,
			self.ppi_network,
			self.graph_rank
		)
		n_drugs = np.sum(individual)
		#if (n_drugs > 100) or (n_drugs <= 0):
		#	n_drugs = 316
		linear_index = np.triu_indices(n_drugs, k=1)
		return (
			np.mean(smiles_submat.values[linear_index]),
			np.mean(moa_submat.values[linear_index]),
			np.mean(paths_submat.values[linear_index]),
			coverage_pval,
			n_drugs
		)
