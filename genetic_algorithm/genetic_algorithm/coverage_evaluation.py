import pandas as pd
import numpy as np
import itertools as it
from scipy.stats import norm

def compute_cov_score(targets, ppi_network, graph_rank, order=1):
	target_neighbor_list = dict()
	for t in targets:
		target_neighbor_list[t] = [t] + [ppi_network.vs["name"][n] for n in ppi_network.neighborhood(t, order=order)]
	# ogni vicino di ogni target viene considerato una sola volta
	covered_nodes_overall = set(it.chain(*[v for k,v in target_neighbor_list.items()]))
	# calcolo il rank mediano ti tutti i vicini considerati according to il ranking nel file graph_rank
	inform_rank = graph_rank.loc[covered_nodes_overall].median()[0]
	# calcola la proporzione di coverage della network di prior knowledge
	coverage_overall = float(len(covered_nodes_overall)) / len(ppi_network.vs)
	#chiedi a Giovanni se va bene, altrimenti usiamo solo il coverage
	return coverage_overall/inform_rank

def generate_random_targets(targets, deg_dist):
	drugNodes_bins = deg_dist.loc[targets].groupby(by="bin").count()
	random_targets = []
	for row in drugNodes_bins.itertuples():
		bin_id, freq = row
		if freq == 0:
			continue
		candidates = deg_dist[deg_dist["bin"] == bin_id]
		random_targets.extend(candidates.sample(n=freq, replace=False).index.to_list())
	return random_targets

def get_drug_input_coverage(candidate_drugs, ppi_network, drug_target_df, graph_rank, deg_dist, order=1):
	# estrai i target delle drug selezionate dall'individuo corrente del GA
	targets = set(drug_target_df.loc[drug_target_df["dat.drug.molecule_name"].isin(candidate_drugs), "dat.target.gene_info.symbol"])
	# cerca solo i target che sono annotati nel grafo di prior knowledge
	targets = list(set(ppi_network.vs["name"]) & targets)
	observed = compute_cov_score(targets, ppi_network, graph_rank, order=1)
	simulated = []
	for i in range(100):
		random_targets = generate_random_targets(targets, deg_dist)
		simulated.append(compute_cov_score(random_targets, ppi_network, graph_rank, order))
	z_score = (observed - np.mean(simulated)) / np.std(simulated)
	return z_score, norm.cdf(-abs(z_score))

def coverage_sum(candidate_drugs, L1000_drug_targets, ppi_network, graph_rank):
	# prendi da L1000 i target del set di 316 drug che stiamo ottimizzando
	drug_target_df = L1000_drug_targets[L1000_drug_targets["dat.drug.molecule_name"].isin(candidate_drugs)]
	# costruisco un DF con tutti i nodi e i loro degree in ordine crescente
	deg_dist = pd.DataFrame(data = ppi_network.degree(), index=ppi_network.vs["name"], columns=["dist"])
	deg_dist.sort_values(by="dist", inplace=True)
	deg_dist["bin"] = pd.qcut(deg_dist.dist, q=20, labels=range(20))
	return get_drug_input_coverage(
		candidate_drugs=candidate_drugs, 
		ppi_network=ppi_network, 
		drug_target_df=drug_target_df, 
		graph_rank=graph_rank, 
		deg_dist=deg_dist
	)
