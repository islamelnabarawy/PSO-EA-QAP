#ifndef PROBLEM_H
#define PROBLEM_H

#include <Dense>

#include "Permutation.h"

class Problem
{
public:
	Problem();
	explicit Problem(std::string problemFilename);
	Problem(std::string problemFilename, std::string solutionFilename);
	~Problem();

	int size() const;

	Permutation::fitness_type getFitness(const Permutation& p) const;
	Permutation::fitness_type getBestFitness() const;
	const Permutation& getBestSolution() const;

private:
	int m_size;
	Eigen::MatrixXi m_flow, m_distance;

	Permutation::fitness_type m_bestFitness;
	Permutation m_bestSolution;
};

#endif // PROBLEM_H
