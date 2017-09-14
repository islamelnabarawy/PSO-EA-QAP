#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <random>
#include <unordered_map>
#include <tuple>

#include "Problem.h"

class Optimizer
{
public:
	enum NeighborhoodType {
		HYPERCUBE, HIERARCHY
	};
	enum NeighborSelectionType {
		LOCAL, FULLY_INFORMED
	};
	enum CrossoverNeighborhood {
		FULL_POPULATION, LOCAL_NEIGHBORS
	};
	struct Settings {
		double probability_crossover;
		double probability_mutation;
		double probability_directed_inital;
		double probability_directed_end;
		bool adaptive_directed;
		bool mutually_exclusive;
		NeighborhoodType neighborhood_type;
		NeighborSelectionType neighbor_selection;
		CrossoverNeighborhood crossover_neighborhood;
		bool use_fitness_map;
		int hierarchy_degree;

		Settings();
	};
	struct Counters {
		int crossover;
		int mutation;
		int directed_mutation;

		Counters();
	};

	Optimizer(const Problem& problem, int populationSize, Settings settings = Settings{});
	~Optimizer();

	Permutation run(int numIterations);

	void dumpPopulation(std::ostream& output) const;
	void dumpFitnessMap(std::ostream& output) const;
	void dumpBestFitnessInMap(std::ostream& output) const;
	
	void dumpHierarchyIX(std::ostream& output) const;
	void dumpHierarchyBF(std::ostream& output) const;
	void dumpHierarchyDF(std::ostream& output) const;

	void dumpCounters(std::ostream& output) const;
	Counters getCounters() const;

private:
	std::vector<Permutation> m_population;
	Problem m_problem;
	Settings m_settings;
	Counters m_counters;

	using fitness_map = std::unordered_map<Permutation::hash_type, Permutation::fitness_type>;
	fitness_map m_fitnessMap;

	// hierarchy node (parent, sibling, first child)
	using hierarchy_node = std::tuple<int, int, int>;
	std::vector<hierarchy_node> m_hierarchy;

	Permutation runExclusive(int numIterations);
	Permutation runInclusive(int numIterations);

	int getFittest(const std::vector<Permutation>& list) const;
	int rouletteWheelSelection(const std::vector<Permutation>& list) const;

	void crossoverMethod(const Permutation& parent1, const Permutation& parent2, Permutation& child1, Permutation& child2);

	Permutation swapMutate(const Permutation& child);

	Permutation directedMutationlb(const Permutation& child, int index);
	Permutation directedMutationfi(const Permutation& child, int index);

	void updatePopulation();

	void updateHierarchy();

	double directedProbability(double progress);

	void updateFitness(Permutation& p);

	std::vector<int> getCrossoverNeighbors(size_t index) const;
	std::vector<int> getMutationNeighbors(size_t index) const;
	// get the hypercube neighbors of "index"
	std::vector<int> getNeighborsHypercube(size_t index) const;
	// get the crossover hierarchy neighbors of "index"
	std::vector<int> getCrossoverNeighborsHierarchy(size_t index) const;
	// get the mutation hierarchy neighbors of "index"
	std::vector<int> getMutationNeighborsHierarchy(size_t index) const;
	
	void initHierarchy();
	void dumpHierarchyNode(std::ostream& output, int index, int tabs = 0) const;

	// random number generator
	static std::default_random_engine s_generator;
};

#endif // OPTIMIZER_H
