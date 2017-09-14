#include "Problem.h"

#include <cassert>
#include <fstream>

Problem::Problem()
: m_size{ 0 }
, m_flow{}
, m_distance{}
, m_bestFitness{ 0 }
, m_bestSolution{}
{ }

Problem::Problem(std::string problemFilename)
: Problem{} // call the default constructor first
{
	auto input = std::ifstream{};

	input.open(problemFilename);

	if (input) {
		// read problem size n
		input >> m_size;

		// read flow matrix A
		m_flow = Eigen::MatrixXi{ m_size, m_size };
		for (auto row = 0; row < m_size; ++row) {
			for (auto col = 0; col < m_size; ++col) {
				input >> m_flow(row, col);
			}
		}

		// read distance matrix B
		m_distance = Eigen::MatrixXi{ m_size, m_size };
		for (auto row = 0; row < m_size; ++row) {
			for (auto col = 0; col < m_size; ++col) {
				input >> m_distance(row, col);
			}
		}
	}

	input.close();
}

Problem::Problem(std::string problemFilename, std::string solutionFilename)
: Problem{ problemFilename } // call the problem file constructor first
{
	auto input = std::ifstream{};

	input.open(solutionFilename);

	if (input) {
		auto temp = 0;
		// read problem size n
		input >> temp;
		// make sure problem and solution have the same size
		assert(temp == m_size);
		// make sure our hash_char_type can hold the permutations
		assert(m_size <= std::numeric_limits<Permutation::hash_char_type>::max());
		// read solution fitness
		input >> m_bestFitness;
		// read solution permutation
		m_bestSolution = Permutation{ m_size, -1 };
		for (auto i = 0; i < m_size; ++i) {
			input >> m_bestSolution[i];
		}
		m_bestSolution.setFitness(m_bestFitness);
	}

	input.close();
}

Problem::~Problem()
{ }

int Problem::size() const
{
	return m_size;
}

Permutation::fitness_type Problem::getFitness(const Permutation& p) const
{
	long long fitness{ 0 };
	for (auto i = 0; i < m_size; ++i) {
		for (auto j = 0; j < m_size; ++j) {
			fitness = fitness + m_flow(i, j) * m_distance(p[i], p[j]);
		}
	}
	return fitness;
}

Permutation::fitness_type Problem::getBestFitness() const
{
	return m_bestFitness;
}

const Permutation& Problem::getBestSolution() const
{
	return m_bestSolution;
}
