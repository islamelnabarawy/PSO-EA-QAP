// main.cpp
#include <iostream>
#include <iomanip>
#include <fstream> 

#include "Permutation.h"
#include "Problem.h"
#include "Optimizer.h"

#include "utils.h"


// Command line options
#include <gflags/gflags.h>
DEFINE_string(problem_input, "test.dat", "QAP problem input file (QAPLIB format)");
DEFINE_string(problem_solution, "test.sln", "QAP problem solution file (QAPLIB format)");
DEFINE_int32(population_size, 50, "population size");
DEFINE_int32(num_iterations, 10000, "number of iterations");
DEFINE_double(mutation, 0, "mutation probability");
DEFINE_double(directed_inital, 0.1, "directed mutation probability");
DEFINE_double(directed_end, 0.1, "directed mutation probability");
DEFINE_double(crossover, 0.7, "crossover probability");
DEFINE_bool(local_best, false , "local best neighbor selection");
DEFINE_bool(adaptive, false, "local best neighbor selection");
DEFINE_bool(mutually_exclusive, false, "mutual exclusion between mutation, crossover, and directed mutation");
DEFINE_bool(local_crossover, false, "select crossover parents from local neighborhood");
DEFINE_string(out_dir, "output/EA", "output directory for files");
DEFINE_string(out_file, "test.txt", "output filename");
DEFINE_string(dump_file, "dump.txt", "dump output filename");
DEFINE_bool(dump_fitness, false, "dump fitness map (takes a long time and a lot of memory!)");
DEFINE_bool(use_fitness_map, false, "track all the generated permutations in a fitness map");
DEFINE_bool(hierarchical_topology, false, "use a hierarchical topology instead of the hypercube");
DEFINE_int32(hierarchy_degree, 4, "the degree of the tree in the hierarchical topology (only relevant if hierarchical_topology is true)");

void print_flag_values(std::ostream& output) {
	std::vector<google::CommandLineFlagInfo> flags;
	google::GetAllFlags(&flags);

	auto max = std::max_element(flags.begin(), flags.end(), 
		[](const google::CommandLineFlagInfo& f1, const google::CommandLineFlagInfo& f2) {
			return f1.name.length() < f2.name.length();
		});
	auto max_width = max->name.length();

	for (auto f : flags) {
		if (f.filename == "main.cpp") {
			output << "\t" << std::setw(max_width) << f.name << " = " << f.current_value << std::endl;
		}
	}
}

// main function:

int main(int argc, char** argv) {
	// using Gflags for command line parameter parsing
	google::SetVersionString("1.0.0. Written by Ayah Helal <ayah.helal@aucegypt.edu> and Islam Elnabarawy <islam.o@aucegypt.edu>.");
	google::SetUsageMessage("Run experiment.");
	google::ParseCommandLineFlags(&argc, &argv, true);
	
	// create the output directory
	if (0 != utils::mkdir_recursive(FLAGS_out_dir.c_str())) {
		std::cerr << "Failed to create output folder: " << FLAGS_out_dir.c_str() << std::endl;
		return -1;
	}
		
	std::ofstream ofs(FLAGS_out_dir + utils::SEPARATOR  + FLAGS_out_file, std::ofstream::out);

	std::ofstream dump_output(FLAGS_out_dir + utils::SEPARATOR + FLAGS_dump_file, std::ofstream::out);

	ofs << "Reading problem instance from: " << FLAGS_problem_input << " and " << FLAGS_problem_solution << std::endl;
	
	ofs << "Command line flags:" << std::endl;
	print_flag_values(ofs);

	// read problem instance
	auto instance = Problem{ FLAGS_problem_input, FLAGS_problem_solution };

	// prepare the settings for the optimizer
	auto settings = Optimizer::Settings{};
	settings.mutually_exclusive = FLAGS_mutually_exclusive;
	settings.probability_crossover = FLAGS_crossover;
	settings.probability_mutation = FLAGS_mutation;
	settings.probability_directed_inital = FLAGS_directed_inital;
	settings.probability_directed_end = FLAGS_directed_end;
	settings.adaptive_directed = FLAGS_adaptive;
	settings.neighbor_selection = FLAGS_local_best ? Optimizer::LOCAL : Optimizer::FULLY_INFORMED;
	settings.crossover_neighborhood = FLAGS_local_crossover ? Optimizer::LOCAL_NEIGHBORS : Optimizer::FULL_POPULATION;
	settings.use_fitness_map = FLAGS_use_fitness_map;
	settings.neighborhood_type = FLAGS_hierarchical_topology ? Optimizer::HIERARCHY : Optimizer::HYPERCUBE;
	settings.hierarchy_degree = FLAGS_hierarchy_degree;

	// initialize the optimizer
	auto optimizer = Optimizer{ instance, FLAGS_population_size, settings };

	dump_output << "Initial population:" << std::endl;
	optimizer.dumpPopulation(dump_output);
	dump_output << std::endl;

	// run the optimizer
	auto bestSolution = optimizer.run(FLAGS_num_iterations);

	ofs << "Found fitness: " << bestSolution.getFitness() << std::endl
		<< "Found solution: " << bestSolution << std::endl
		<< " -- compare to --" << std::endl
		<< "Best Known Fitness: " << instance.getBestFitness() << std::endl
		<< "Best Known Solution: " << instance.getBestSolution() << std::endl;

	ofs.close();

	dump_output << "Final population:" << std::endl;
	optimizer.dumpPopulation(dump_output);
	dump_output << std::endl;

	dump_output << "Counters:" << std::endl;
	optimizer.dumpCounters(dump_output);
	dump_output << std::endl;

	if (FLAGS_use_fitness_map) {
		if (FLAGS_dump_fitness) {
			dump_output << "Fitness map:" << std::endl;
			optimizer.dumpFitnessMap(dump_output);
			dump_output << std::endl;
		} else {
			dump_output << "Best in fitness map:" << std::endl;
			optimizer.dumpBestFitnessInMap(dump_output);
			dump_output << std::endl;
		}
	}

	if (FLAGS_hierarchical_topology) {
		dump_output << "Hierarchy breadth-first:" << std::endl;
		optimizer.dumpHierarchyBF(dump_output);
		dump_output << std::endl;
		dump_output << "Hierarchy depth-first:" << std::endl;
		optimizer.dumpHierarchyDF(dump_output);
		dump_output << std::endl;
	}

	dump_output.close();

	return 0;
}
