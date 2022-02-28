import random
import pandas as pd
import numpy as np
import igraph as ig
from deap import base
from deap import creator
from deap import tools
from deap import algorithms
import matplotlib.pyplot as plt
import sys
import click
from genetic_algorithm import GA
from genetic_algorithm import EvaluationFunction

def load_drugs(drug_names_path):
	with open(drug_names_path, "r") as f:
		drug_names = [line.strip() for line in f]
	return drug_names

def setup_ga(
	dataset,
	population_size=100,
	attribute_init_prob=0.3,
	attribute_mutation_prob=0.1,
	n_hall_of_fame=10
	):

	datadir = "./data/TCGA-{}/".format(dataset)
	# setup evaluation function
	evaluate = EvaluationFunction(
		datadir + "SMILES_{}_distance_matrix.rds".format(dataset),
		datadir + "SMILES_{}_distance_matrix.rds".format(dataset),
		datadir + "shortestpath_{}_distance_matrix.rds".format(dataset),
		datadir + "prior_network_all_3_layers_alisa_{}.gml".format(dataset),
		datadir + "graph_rank_{}_mat.rds".format(dataset),
		datadir + "filter_by_L1000_drugs_on_prior_{}.txt".format(dataset)
	)

	creator.create("Fitness", base.Fitness, weights=(1.0, 1.0, 1.0, -1.0, -1.0))
	creator.create("Individual", list, fitness=creator.Fitness)

	IND_SIZE = len(evaluate.drug_names)

	toolbox = base.Toolbox()
	toolbox.register("init_attribute", lambda: np.random.random() < attribute_init_prob)
	toolbox.register("individual", tools.initRepeat, creator.Individual,
					toolbox.init_attribute, n=IND_SIZE)
	toolbox.register("population", tools.initRepeat, list, toolbox.individual, n=population_size)

	toolbox.register("evaluate", evaluate)
	toolbox.register("mutate", tools.mutShuffleIndexes, indpb=attribute_mutation_prob)
	toolbox.register("mate", tools.cxOnePoint)
	toolbox.register("select", tools.selNSGA2)

	stats = tools.Statistics(lambda ind: ind.fitness.values)
	stats.register("avg_smile", lambda x: np.mean(x, axis=0)[0])
	stats.register("avg_moa", lambda x: np.mean(x, axis=0)[1])
	stats.register("avg_paths", lambda x: np.mean(x, axis=0)[2])
	stats.register("avg_coverage", lambda x: np.mean(x, axis=0)[3])
	stats.register("n_drugs", lambda x: np.mean(x, axis=0)[4])

	hof = tools.HallOfFame(n_hall_of_fame)

	return toolbox, stats, hof

def print_individual(individual, drug_names):
	drugs = [drug_names[i] for i, attr in enumerate(individual) if attr]
	print(drugs)
	smiles_dist, moa_dist, graph_dist, coverage, n_drugs = individual.fitness.values
	print("\tNr. of drugs", n_drugs)
	print("\tSmiles Distance", smiles_dist)
	print("\tMoa Distance", moa_dist)
	print("\tGraph Distance", graph_dist)
	print("\tTarget coverage", coverage)


def plot_trace(dataset, log, output_dir):
	gen, n_drugs, avg_smile, avg_moa, avg_paths, avg_coverage = log.select(
		"gen", "n_drugs", "avg_smile", "avg_moa", "avg_paths", "avg_coverage"
	)

	fig, drugs_ax = plt.subplots()
	drugs_ax.set_xlabel("Generation")
	drugs_ax.set_ylabel("Nr. Drugs")
	drugs_plot, = drugs_ax.plot(gen, n_drugs, label="n drugs", color="C0")
	drugs_ax.yaxis.label.set_color(drugs_plot.get_color())

	smiles_ax = drugs_ax.twinx()
	smiles_ax.set_ylabel("SMILES distances")
	smiles_plot, = smiles_ax.plot(gen, avg_smile, label="average smiles", color="C1")
	smiles_ax.yaxis.label.set_color(smiles_plot.get_color())

	moa_ax = drugs_ax.twinx()
	moa_ax.set_ylabel("MOA distances")
	moa_plot, = moa_ax.plot(gen, avg_moa, label="average moa", color="C2")
	moa_ax.yaxis.label.set_color(moa_plot.get_color())
	moa_ax.spines['right'].set_position(('outward', 60))
	#moa_ax.xaxis.set_ticks([])

	paths_ax = drugs_ax.twinx()
	paths_ax.set_ylabel("Shortest Path distances")
	paths_plot, = paths_ax.plot(gen, avg_paths, label="average paths", color="C3")
	paths_ax.yaxis.label.set_color(paths_plot.get_color())
	paths_ax.spines['right'].set_position(('outward', 120))
	#paths_ax.xaxis.set_ticks([])

	coverage_ax = drugs_ax.twinx()
	coverage_ax.set_ylabel("Target coverage")
	coverage_plot, = coverage_ax.plot(gen, avg_coverage, label="average coverage", color="C4")
	coverage_ax.yaxis.label.set_color(coverage_plot.get_color())
	coverage_ax.spines['right'].set_position(('outward', 180))

	fig.legend(loc="lower right")
	fig.tight_layout()
	fig.savefig("{}/trace_{}.pdf".format(output_dir, dataset))
	return fig

@click.command()
@click.option('--dataset', default="BRCA", help="Dataset to use")
@click.option('--output_dir', default=".", help="Output directory")
def main(dataset, output_dir):
	toolbox, stats, hof = setup_ga(dataset)
	drug_names = toolbox.evaluate.func.drug_names
	pop = toolbox.population()
	pop = toolbox.select(pop, len(pop))
	pop, log = GA(pop, toolbox, cxpb=0.7, mutpb=0.3, ngen=2500, stats=stats, verbose=True, halloffame=hof)

	print("Hall of fame")
	for ind in hof:
		print_individual(ind, drug_names)

	fig = plot_trace(dataset, log, output_dir)
	plt.show()

if __name__ == "__main__":
	main()
