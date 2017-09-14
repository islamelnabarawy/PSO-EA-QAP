#include "Optimizer.h"

#include <iomanip>
#include <numeric>
#include <map>
#include <queue>
#include <algorithm>

#include <time.h>

Optimizer::Settings::Settings()
: probability_crossover{ 0.5 }
, probability_mutation{ 0.5 }
, probability_directed_inital{ 0.5 }
, probability_directed_end{ 0.5 }
, adaptive_directed{ true }
, mutually_exclusive{ false }
, neighborhood_type{ NeighborhoodType::HIERARCHY }
, neighbor_selection{ NeighborSelectionType::FULLY_INFORMED }
, crossover_neighborhood{ CrossoverNeighborhood::LOCAL_NEIGHBORS }
, use_fitness_map{ true }
, hierarchy_degree{ 3 }
{ }

Optimizer::Counters::Counters()
: crossover{ 0 }
, mutation{ 0 }
, directed_mutation{ 0 }
{ }

Optimizer::Optimizer(const Problem& problem, int populationSize, Optimizer::Settings settings)
: m_problem{ problem }
, m_population{}
, m_settings{ settings }
, m_counters{}
{
	// initialize the population
	for (auto index = 0; index < populationSize; ++index) {
		Permutation perm = Permutation{ m_problem.size() };
		updateFitness(perm);
		m_population.push_back(perm);
	}
	if (m_settings.neighborhood_type == NeighborhoodType::HIERARCHY) {
		// initialize the hierarchy
		initHierarchy();
		
	}
	updatePopulation();
}

Optimizer::~Optimizer()
{
}

Permutation Optimizer::run(int numIterations)
{
	if (m_settings.mutually_exclusive)
		return runExclusive(numIterations);
	else
		return runInclusive(numIterations);
}

Permutation Optimizer::runInclusive(int numIterations)
{
	auto uniform = std::uniform_real_distribution<double>{0.0, 1.0};

	auto bestSolution = Permutation{ m_problem.size(), -1 };
	auto bestFitness = std::numeric_limits<Permutation::fitness_type>::max();

	for (auto i = 0; i < numIterations; ++i) {
		for (size_t index = 0; index < m_population.size(); ++index) {
			auto list = std::vector<Permutation>{};
			// add the two parents to the list
			list.push_back(m_population[index]);

			if (m_settings.crossover_neighborhood == FULL_POPULATION) {
				list.push_back(m_population[rouletteWheelSelection(m_population)]);
			} else {
				auto neighbors = getCrossoverNeighbors(index);
				auto neighborhood = std::vector<Permutation>{};
				for (auto n : neighbors) { 
					neighborhood.push_back(m_population[n]); 
				}
				list.push_back(neighborhood[rouletteWheelSelection(neighborhood)]);
			}

			auto first_index = 0;

			auto prob = uniform(s_generator);
			if (prob < m_settings.probability_crossover) {
				// do crossover
				++m_counters.crossover;
				auto child1 = Permutation{ m_problem.size(), -1 };
				auto child2 = Permutation{ m_problem.size(), -1 };
				crossoverMethod(list[first_index], list[first_index + 1], child1, child2);
				// add the results to the list
				list.push_back(child1);
				list.push_back(child2);
				first_index += 2;
			}

			prob = uniform(s_generator);
			if (prob < m_settings.probability_mutation) {
				//do mutate
				++m_counters.mutation;
				auto child1 = swapMutate(list[first_index]);
				auto child2 = swapMutate(list[first_index + 1]);
				// add the results to the list
				list.push_back(child1);
				list.push_back(child2);
				first_index += 2;
			}

			prob = uniform(s_generator);
			auto directedprobabiliy = directedProbability((double)((double)i / (double)numIterations));
			if (prob < directedprobabiliy) {
				//do directed mutate
				++m_counters.directed_mutation;
				auto child1 = Permutation{};
				auto child2 = Permutation{};
				if (m_settings.neighbor_selection == NeighborSelectionType::LOCAL) {
					child1 = directedMutationlb(list[first_index], index);
					child2 = directedMutationlb(list[first_index+1], index);
				} else {
					child1 = directedMutationfi(list[first_index], index);
					child2 = directedMutationfi(list[first_index + 1], index);
				}
				// add the results to the list
				list.push_back(child1);
				list.push_back(child2);
				first_index += 2;
			}

			for (auto item : list) {
				if (item.getFitness() < bestFitness) {
					bestSolution = item;
					bestFitness = item.getFitness();

					if (bestFitness <= m_problem.getBestFitness()) {
						return bestSolution;
					}
				}
			}

			m_population[index] = list[rouletteWheelSelection(list)];
		}

		updatePopulation();

	}
	return bestSolution;
}

double Optimizer::directedProbability(double progress)
{
	double probability = m_settings.probability_directed_inital;
	if (m_settings.adaptive_directed)
	{
		probability += ((m_settings.probability_directed_end - m_settings.probability_directed_inital) * progress);
	}
	return probability;
}

Permutation Optimizer::runExclusive(int numIterations)
{
	auto uniform = std::uniform_real_distribution<double>{0.0, 1.0};

	auto bestSolution = Permutation{ m_problem.size(), -1 };
	auto bestFitness = std::numeric_limits<Permutation::fitness_type>::max();

	for (auto i = 0; i < numIterations; ++i) {
		for (size_t index = 0; index < m_population.size(); ++index) {
			auto list = std::vector<Permutation>{};
			// add the two parents to the list
			list.push_back(m_population[index]);

			if (m_settings.crossover_neighborhood == FULL_POPULATION) {
				list.push_back(m_population[rouletteWheelSelection(m_population)]);
			} else {
				auto neighbors = getCrossoverNeighbors(index);
				auto neighborhood = std::vector<Permutation>{};
				for (auto n : neighbors) {
					neighborhood.push_back(m_population[n]);
				}
				list.push_back(neighborhood[rouletteWheelSelection(neighborhood)]);
			}

			auto first_index = 0;

			auto probabilty = uniform(s_generator);
			double sum = m_settings.probability_crossover;

			if (sum < probabilty) {
				// do crossover
				++m_counters.crossover;
				auto child1 = Permutation{ m_problem.size(), -1 };
				auto child2 = Permutation{ m_problem.size(), -1 };
				crossoverMethod(list[first_index], list[first_index + 1], child1, child2);
				// add the results to the list
				list.push_back(child1);
				list.push_back(child2);
				first_index += 2;
			} else {
				sum += m_settings.probability_mutation;
				if (sum < probabilty) {
					//do mutate
					++m_counters.mutation;
					auto child1 = swapMutate(list[first_index]);
					auto child2 = swapMutate(list[first_index + 1]);
					// add the results to the list
					list.push_back(child1);
					list.push_back(child2);
					first_index += 2;
				} else {
					auto directedprobabiliy = directedProbability((double)((double)i / (double)numIterations));
					sum += directedprobabiliy;
					if (sum < probabilty) {
						//do directed mutate
						++m_counters.directed_mutation;
						auto child1 = Permutation{};
						auto child2 = Permutation{};
						if (m_settings.neighbor_selection == NeighborSelectionType::LOCAL) {
							child1 = directedMutationlb(list[first_index], index);
							child2 = directedMutationlb(list[first_index + 1], index);
						} else {
							child1 = directedMutationfi(list[first_index], index);
							child2 = directedMutationfi(list[first_index + 1], index);
						}
						// add the results to the list
						list.push_back(child1);
						list.push_back(child2);
						first_index += 2;
					}
				}
			}

			for (auto item : list) {
				if (item.getFitness() < bestFitness) {
					bestSolution = item;
					bestFitness = item.getFitness();

					if (bestFitness <= m_problem.getBestFitness()) {
						return bestSolution;
					}
				}
			}

			m_population[index] = list[rouletteWheelSelection(list)];
		}
	}
	return bestSolution;
}

int Optimizer::getFittest(const std::vector<Permutation>& list) const
{
	auto index = -1;
	auto max = std::numeric_limits<Permutation::fitness_type>::max();

	for (size_t i = 0; i < list.size(); ++i) {
		if (max > list[i].getFitness()) {
			index = i;
			max = list[i].getFitness();
		}
	}

	return index;
}

int Optimizer::rouletteWheelSelection(const std::vector<Permutation>& list) const
{
	// sum the fitness using std::accumulate and a lambda expression
	auto sumFitness = std::accumulate(list.begin(), list.end(), Permutation::fitness_type{ 0 },
		[](const Permutation::fitness_type& sum, const Permutation& item) { return sum + item.getFitness(); });

	if (sumFitness == 0) {
		return 0;
	}

	auto randProbability = std::uniform_int_distribution<Permutation::fitness_type>{0, sumFitness - 1}(s_generator);

	auto sumProbability = Permutation::fitness_type{ 0 };
	auto index = 0;

	while (sumProbability <= randProbability) {
		sumProbability += list[index].getFitness();
		++index;
	}

	return index - 1;
}

void Optimizer::crossoverMethod(const Permutation& parent1, const Permutation& parent2, Permutation& child1, Permutation& child2)
{
	auto bitIndex = std::uniform_int_distribution<int>{0, m_problem.size()}(s_generator);

	for (auto i = 0; i < bitIndex; i++) {
		child1[i] = parent1[i];
		child2[i] = parent2[i];
	}

	auto count1 = bitIndex;
	auto count2 = bitIndex;
	for (auto i = 0; i < m_problem.size(); i++) {
		if (find(child1.begin(), child1.end(), parent2[i]) == child1.end()) {
			child1[count1] = parent2[i];
			count1++;
		}
		if (find(child2.begin(), child2.end(), parent1[i]) == child2.end()) {
			child2[count2] = parent1[i];
			count2++;
		}
	}
	//Recalculate the fittness and deltafittness
	updateFitness(child1);
	updateFitness(child2);

}

Permutation Optimizer::swapMutate(const Permutation& child)
{
	auto dist = std::uniform_int_distribution<int>{0, child.size() - 1};

	int bitIndex1 = dist(s_generator);
	int bitIndex2 = dist(s_generator);

	auto copy = Permutation{ child };
	copy.swap(bitIndex1, bitIndex2);
	updateFitness(copy);
	return copy;
}

//directed mutation
Permutation Optimizer::directedMutationlb(const Permutation& child, int index)
{
	auto neighbors = getMutationNeighbors(index);
	std::vector<Permutation> list;
	for (auto n : neighbors) {
		list.push_back(m_population[n]);
	}
	auto bindex = getFittest(list);

	int bitIndex1 = std::uniform_int_distribution<int>{0, m_problem.size() - 1}(s_generator);
	int bitIndex2 = m_population[index].findIndex(m_population[bindex][bitIndex1]);

	auto copy = Permutation{ child };
	copy.swap(bitIndex1, bitIndex2);
	updateFitness(copy);
	return copy;
}

//directed mutation
Permutation Optimizer::directedMutationfi(const Permutation& child, int index)
{
	auto neighbors = getMutationNeighbors(index);
	std::vector<Permutation> list;
	auto copy = Permutation{ child };

	for (auto n : neighbors) {
		int bitIndex1 = std::uniform_int_distribution<int>{0, m_problem.size() - 1}(s_generator);
		int bitIndex2 = m_population[index].findIndex(m_population[n][bitIndex1]);
		copy.swap(bitIndex1, bitIndex2);
	}

	updateFitness(copy);
	return copy;
}

// get the neighbors of "index"
std::vector<int> Optimizer::getCrossoverNeighbors(size_t index) const
{
	if (m_settings.neighborhood_type == HIERARCHY) {
		return getCrossoverNeighborsHierarchy(index);
	}
	return getNeighborsHypercube(index);
}

// get the neighbors of "index"
std::vector<int> Optimizer::getMutationNeighbors(size_t index) const
{
	if (m_settings.neighborhood_type == HIERARCHY) {
		return getMutationNeighborsHierarchy(index);
	}
	return getNeighborsHypercube(index);
}

// get the hypercube neighbors of "index"
std::vector<int> Optimizer::getNeighborsHypercube(size_t index) const
{
	auto list = std::vector<int>{};

	const size_t sz = sizeof(index)* 8;

	for (size_t i = 0; i < sz; ++i) {
		// flip the ith bit
		auto x = index ^ (1 << i);
		// add it to the list if it's within bounds
		if (x < m_population.size()) {
			list.push_back(x);
		}
	}

	return list;
}

// get the crossover hierarchy neighbors of "index"
std::vector<int> Optimizer::getMutationNeighborsHierarchy(size_t index) const
{
	auto list = std::vector<int>{};

	if (index == 0)
		list.push_back(0);
		

	// get the current element's parent
	auto c = std::get<0>(m_hierarchy[index]);

	// add the ancestors to the list
	while (c != -1) {
		list.push_back(c);
		// get the next parent
		c = std::get<0>(m_hierarchy[c]);
	}

	return list;
}

// get the mutation hierarchy neighbors of "index"
std::vector<int> Optimizer::getCrossoverNeighborsHierarchy(size_t index) const
{
	auto list = std::vector<int>{};

	// get the current element's parent
	auto c = std::get<0>(m_hierarchy[index]);

	if (c != -1)
		list.push_back(c);
	else{
		list.push_back(rouletteWheelSelection(m_population));
	}
		

	return list;
}

void Optimizer::updatePopulation()
{
	if (m_settings.neighborhood_type == HIERARCHY) {
		updateHierarchy();
	}
}

void Optimizer::updateHierarchy()
{
	std::queue<int> q;
	q.push(0);
	
	while (!q.empty()) {
		auto node = q.front();
		q.pop();

		// temp holds the population-based indices of the children
		std::vector<int> temp;
		// list holds the actual permutation of the children
		auto list = std::vector<Permutation> {};

		temp.push_back(node);
		list.push_back(m_population[node]);

		// get the current element's first child
		auto c = std::get<2>(m_hierarchy[node]);

		// push the children to the queue
		while (c != -1) {
			temp.push_back(c);
			list.push_back(m_population[c]);
			
			q.push(c);
			// get the next sibling
			c = std::get<1>(m_hierarchy[c]);
		}

		// get the list-index of the fittest permutation
		auto index = getFittest(list);
		
		// temp[index] has the corresponding population-based index
		if (node != temp[index]) {
			// if a child is better than its parent, swap them
			std::swap(m_population[temp[index]], m_population[node]);
		}
	}
}

void Optimizer::initHierarchy()
{
	// add the root node statically
	m_hierarchy.push_back(hierarchy_node{ -1, -1, 1 });
	// handy values
	const auto degree = m_settings.hierarchy_degree;
	const auto maxNodes = static_cast<int>(m_population.size());
	// start with level 1
	auto level = 1;
	// nodes in the first level = degree of tree
	int nodes = degree;
	// start with node 1
	auto currentNode = 1;
	// keep adding nodes until tree is full
	while (currentNode < maxNodes) {
		// add the nodes in the current level
		for (int i = 0; i < nodes && currentNode < maxNodes; ++i) {
			// calculate the index of the current node's parent dynamically
			auto parent = (currentNode - 1) / degree;
			// calculate the index of the next node
			auto next = (currentNode % degree > 0) ? currentNode + 1 : -1;
			// calculate the index of the current node's first child dynamically
			auto child = nodes + i * degree + (currentNode - i);
			// add the node and increment the currentNode index 
			m_hierarchy.push_back(hierarchy_node{
				parent,
				(next < maxNodes) ? next : -1,
				(child < maxNodes) ? child : -1
			});
			++currentNode;
		}
		// move to the next level
		++level;
		// there are "degree" nodes in the next level for each node in this level
		// so the number of nodes increases exponentially each level
		nodes *= degree;
	}
}

void Optimizer::updateFitness(Permutation& p)
{
	Permutation::fitness_type fitness = std::numeric_limits<Permutation::fitness_type>::max();
	if (m_settings.use_fitness_map) {
		auto hash = p.hash();

		// look for that permutation in the list
		auto elem = m_fitnessMap.find(hash);
		if (elem != m_fitnessMap.end()) {
			// get it from the list instead of re-computing it
			fitness = elem->second;
		} else {
			// it's not in the list, compute it and put it there
			fitness = m_problem.getFitness(p);
			m_fitnessMap[hash] = fitness;
		}
	} else {
		fitness = m_problem.getFitness(p);
	}

	p.setFitness(fitness);
}

void Optimizer::dumpPopulation(std::ostream& output) const
{
	for (auto elem : m_population) {
		output << elem << " -> " << elem.getFitness() << std::endl;
	}
}

void Optimizer::dumpFitnessMap(std::ostream& output) const
{
	std::map<Permutation::fitness_type, Permutation::hash_type> reverseMap{};
	for (auto elem : m_fitnessMap) {
		reverseMap[elem.second] = elem.first;
	}
	for (auto elem : reverseMap) {
		output << "( ";
		for (auto c : elem.second) { output << static_cast<int>(c - 1) << " "; }
		output << ") -> " << elem.first << std::endl;
	}
}

void Optimizer::dumpBestFitnessInMap(std::ostream& output) const
{
	auto bestFitness = std::numeric_limits<Permutation::fitness_type>::max();
	auto bestPermutation = Permutation::hash_type{};

	for (auto elem : m_fitnessMap) {
		if (bestFitness > elem.second) {
			bestPermutation = elem.first;
			bestFitness = elem.second;
		}
	}

	output << "( ";
	for (auto c : bestPermutation) { output << static_cast<int>(c - 1) << " "; }
	output << ") -> " << bestFitness << std::endl;
}

void Optimizer::dumpHierarchyNode(std::ostream& output, int index, int tabs) const
{
	auto h = m_hierarchy[index];
	auto p = m_population[index];
	output << std::setw(3) << index << ":";
	for (auto t = 0; t < tabs; ++t) output << '\t';
	output << std::setw(3) << std::get<0>(h)
		<< " " << std::setw(3) << std::get<1>(h)
		<< " " << std::setw(3) << std::get<2>(h)
		<< " ==> " << p << " -> " << p.getFitness()
		<< std::endl;
}

void Optimizer::dumpHierarchyIX(std::ostream& output) const
{
	auto nodes = 1;
	auto tabs = 0;
	auto i = 0;
	for (size_t x = 0; x < m_hierarchy.size(); ++x) {
		dumpHierarchyNode(output, x, tabs);
		if (++i == nodes) {
			nodes *= m_settings.hierarchy_degree;
			++tabs;
			i = 0;
		}
	}
}

void Optimizer::dumpHierarchyBF(std::ostream& output) const
{
	using hier_pair = std::pair<int, int>; // index, tabs
	std::queue<hier_pair> q;
	q.push(hier_pair{ 0, 0 });
	while (!q.empty()) {
		auto node = q.front();
		q.pop();
		// dump the current node
		dumpHierarchyNode(output, node.first, node.second);

		// get the current element's first child
		auto c = std::get<2>(m_hierarchy[node.first]);

		// push the children to the queue
		while (c != -1) {
			q.push(hier_pair{ c, node.second + 1 });
			// get the next sibling
			c = std::get<1>(m_hierarchy[c]);
		}
	}
}

void Optimizer::dumpHierarchyDF(std::ostream& output) const
{
	using hier_pair = std::pair<int, int>; // index, tabs
	std::vector<hier_pair> stack;
	stack.push_back(hier_pair{ 0, 0 });
	while (!stack.empty()) {
		auto node = stack.back();
		stack.pop_back();
		// dump the current node
		dumpHierarchyNode(output, node.first, node.second);
		
		// get the current element's first child
		auto c = std::get<2>(m_hierarchy[node.first]);
		
		// push the children to the stack
		std::vector<hier_pair> temp;
		while (c != -1) {
			temp.push_back(hier_pair{ c, node.second + 1 });
			// get the next sibling
			c = std::get<1>(m_hierarchy[c]);
		}
		while (!temp.empty()) {
			stack.push_back(temp.back());
			temp.pop_back();
		}
	}
}

void Optimizer::dumpCounters(std::ostream& output) const
{
	output << "- Crossover: " << m_counters.crossover << std::endl
		<< "- Mutation: " << m_counters.mutation << std::endl
		<< "- Directed Mutation: " << m_counters.directed_mutation << std::endl;
}

Optimizer::Counters Optimizer::getCounters() const
{
	return m_counters;
}

std::default_random_engine Optimizer::s_generator{ static_cast<unsigned int>(time(nullptr)) };
