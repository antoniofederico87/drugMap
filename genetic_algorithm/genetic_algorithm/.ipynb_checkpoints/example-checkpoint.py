import random
import numpy as np

from deap import base
from deap import creator
from deap import tools
from deap import algorithms

def asymmFlipBit(individual, indpb_on, indpb_off):
    for i in range(len(individual)):
        if individual[i] == 1 and random.random() < indpb_off:
            individual[i] = type(individual[i])(not individual[i])
        if individual[i] == 0 and random.random() < indpb_on:
            individual[i] = type(individual[i])(not individual[i])

    return individual,

def BasicGA(population, toolbox, cxpb, mutpb, ngen, stats=None, halloffame=None, verbose=__debug__):
	logbook = tools.Logbook()
	logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

	# Evaluate the individuals with an invalid fitness
	invalid_ind = [ind for ind in population if not ind.fitness.valid]
	fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
	for ind, fit in zip(invalid_ind, fitnesses):
		ind.fitness.values = fit

	if halloffame is not None:
		halloffame.update(population)

	record = stats.compile(population) if stats else {}
	logbook.record(gen=0, nevals=len(invalid_ind), **record)
	if verbose:
		print(logbook.stream)

	# Begin the generational process
	for gen in range(1, ngen + 1):
		# Generate new offspring
		offspring = algorithms.varOr(population, toolbox, 20, cxpb, mutpb)

		# evaluate the offspring
		fitnesses = toolbox.map(toolbox.evaluate, offspring)
		for ind, fit in zip(offspring, fitnesses):
			ind.fitness.values = fit

		# Select the next generation of fittest individuals
		population = toolbox.select(population+offspring, len(population))

		# Update the hall of fame with the generated individuals
		if halloffame is not None:
			halloffame.update(offspring)

		# Append the current generation statistics to the logbook
		record = stats.compile(population) if stats else {}
		logbook.record(gen=gen, nevals=len(invalid_ind), **record)
		if verbose:
			print(logbook.stream)

	return population, logbook

smiles_distances = np.loadtxt("data/smiles_distances.csv", delimiter=",")
moa_distances = np.loadtxt("data/moa_distances.csv", delimiter=",")
paths_distances = np.loadtxt("data/shortest_paths_distances.csv", delimiter=",")

def evaluate_smiles(individual):
	if np.sum(individual) <= 1:
		return (0, 0, 0, 316)
	submat_index = np.ix_(individual, individual)
	smiles_submat = smiles_distances[submat_index]
	moa_submat = moa_distances[submat_index]
	paths_submat = paths_distances[submat_index]
	n_drugs = np.sum(individual)
	#if (n_drugs > 100) or (n_drugs <= 0):
	#	n_drugs = 316
	linear_index = np.triu_indices(n_drugs, k=1)
	return (
		np.mean(smiles_submat[linear_index]),
		np.mean(moa_submat[linear_index]),
		np.mean(paths_submat[linear_index]),
		n_drugs
	)
	#return ( 0 , np.sum(individual))
	#return np.sum(individual), 

creator.create("FitnessMax", base.Fitness, weights=(1.0, 1.0, 1.0, -1.0))
creator.create("Individual", list, fitness=creator.FitnessMax)

with open("data/drugs.txt", "r") as f:
	drugs = [line.strip() for line in f]

IND_SIZE = len(drugs)

toolbox = base.Toolbox()
#toolbox.register("attr_bin", random.randint,0,1)
toolbox.register("attr_bin", lambda: np.random.random() < 0.3)
toolbox.register("individual", tools.initRepeat, creator.Individual,
                 toolbox.attr_bin, n=IND_SIZE)
toolbox.register("population", tools.initRepeat, list, toolbox.individual, n=100)

toolbox.register("evaluate", evaluate_smiles)
toolbox.register("mutate", tools.mutShuffleIndexes, indpb=0.1)
toolbox.register("mate", tools.cxOnePoint)
toolbox.register("select", tools.selNSGA2)

stats = tools.Statistics(lambda ind: ind.fitness.values)
stats.register("avg_smile", lambda x: np.mean(x, axis=0)[0])
stats.register("avg_moa", lambda x: np.mean(x, axis=0)[1])
stats.register("avg_paths", lambda x: np.mean(x, axis=0)[2])
stats.register("n_drugs", lambda x: np.mean(x, axis=0)[3])
hof = tools.HallOfFame(10)
pop = toolbox.population()
pop = toolbox.select(pop, len(pop))
pop, log = BasicGA(pop, toolbox, cxpb=0.7, mutpb=0.3, ngen=2000, stats=stats, verbose=True, halloffame=hof)


print("Hall of fame")
for ind in hof:
	print(ind.fitness)

print("Population")
for ind in pop:
	print(ind.fitness)

import matplotlib.pyplot as plt
gen, n_drugs, avg_smile, avg_moa, avg_paths = log.select("gen", "n_drugs", "avg_smile", "avg_moa", "avg_paths")

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

fig.legend(loc="lower right")
fig.tight_layout()
plt.show()


