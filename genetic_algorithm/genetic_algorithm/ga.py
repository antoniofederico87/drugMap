from deap import tools
from deap import algorithms

#def asymmFlipBit(individual, indpb_on, indpb_off):
#    for i in range(len(individual)):
#        if individual[i] == 1 and random.random() < indpb_off:
#            individual[i] = type(individual[i])(not individual[i])
#        if individual[i] == 0 and random.random() < indpb_on:
#            individual[i] = type(individual[i])(not individual[i])
#    return individual,

def GA(population, toolbox, cxpb, mutpb, ngen, stats=None, halloffame=None, verbose=__debug__):
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
