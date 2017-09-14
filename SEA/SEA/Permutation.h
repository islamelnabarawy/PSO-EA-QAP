#ifndef PERMUTATION_H
#define PERMUTATION_H

#include <random>
#include <vector>

class Permutation
{
public:
	Permutation();
	explicit Permutation(int size);
	Permutation(int size, int value);
	~Permutation();

	int operator[](int index) const;
	int& operator[](int index);

	using fitness_type = unsigned long long int;
	using list_type = std::vector<int>;
	using iterator = list_type::iterator;
	using const_iterator = list_type::const_iterator;

	iterator begin();
	iterator end();

	const_iterator begin() const;
	const_iterator end() const;

	using hash_char_type = wchar_t;
	using hash_type = std::wstring;
	using hash_stream_type = std::wstringstream;

	int size() const;
	hash_type hash() const;

	fitness_type getFitness() const;
	void setFitness(fitness_type fitness);

	void swap(int ix1, int ix2);

	// find the index of the facility in permutation
	int findIndex(int facility);

	friend std::ostream& operator<<(std::ostream& out, const Permutation& p);

private:
	list_type m_state;
	fitness_type m_fitness;

	static std::default_random_engine s_generator;
};

#endif // PERMUTATION_H
