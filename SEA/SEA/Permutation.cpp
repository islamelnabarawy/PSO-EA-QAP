#include "Permutation.h"

#include <algorithm>
#include <sstream>
#include <cassert>

#include <time.h>

Permutation::Permutation()
: m_state{}
, m_fitness{ std::numeric_limits<Permutation::fitness_type>::max() }
{ }

Permutation::Permutation(int size) 
: m_state{ size, -1 }
, m_fitness{ std::numeric_limits<Permutation::fitness_type>::max() }
{
	auto dist = std::uniform_int_distribution<int>{0, size - 1};

	auto checklist = std::vector<bool>(size, false);

	for (auto i = 0; i < size;) {
		auto num = dist(s_generator);
		if (!checklist[num]) {
			checklist[num] = true;
			m_state[i] = num;
			++i;
		}
	}
}

Permutation::Permutation(int size, int value) 
: m_state{ size, value }
, m_fitness{ std::numeric_limits<long>::max() }
{ }

Permutation::~Permutation()
{ }

int Permutation::operator[](int index) const
{
	return m_state[index];
}

int& Permutation::operator[](int index)
{
	return m_state[index];
}

Permutation::iterator Permutation::begin()
{
	return m_state.begin();
}

Permutation::iterator Permutation::end()
{
	return m_state.end();
}

Permutation::const_iterator Permutation::begin() const
{
	return m_state.cbegin();
}

Permutation::const_iterator Permutation::end() const
{
	return m_state.cend();
}

int Permutation::size() const
{
	return m_state.size();
}

Permutation::hash_type Permutation::hash() const
{
	hash_stream_type result;
	for (auto i : m_state) {
		result << static_cast<hash_char_type>(i + 1);
	}
	return result.str();
}

Permutation::fitness_type Permutation::getFitness() const
{
	return m_fitness;
}

void Permutation::setFitness(Permutation::fitness_type fitness)
{
	m_fitness = fitness;
}

void Permutation::swap(int ix1, int ix2)
{
	std::swap(m_state[ix1], m_state[ix2]);
}

int Permutation::findIndex(int facility)
{
	int index = -1;

	auto it = find(m_state.begin(), m_state.end(), facility);
	if (it != m_state.end()) {
		index = it - m_state.begin();
	}

	return index;
}


std::ostream& operator<<(std::ostream& out, const Permutation& p) {
	out << "( ";
	for (auto x : p.m_state) {
		out << x << " ";
	}
	out << ")";
	return out;
}

std::default_random_engine Permutation::s_generator{ static_cast<unsigned int>(time(nullptr)) };