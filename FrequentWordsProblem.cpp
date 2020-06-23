#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

int get_code(char letter);

char get_letter(int code);

int count_mers(std::vector<int>& mers_counter, int k, const std::string& database);

int count_frequence(int num_of_mers, const std::vector<int>& mers_counter, std::vector<int>& most_frequent_index);

int convert_to_letters(std::vector<std::string>& most_freq_mers, const std::vector<int>& most_frequent_index, int mer_length);

int short_term_processing(int mer_length, const std::string& database, std::vector<std::string>& most_freq_mers);

int long_term_processing(int mer_length, const std::string& database, std::vector<std::string>& most_freq_mers);

int count_mer(const std::string& database, int mer_start, int mer_length);

const int Large_term_length = 12;

int main()
{
	std::string database;

	std::ifstream input("cholera_genome.txt");

	input >> database;

	int base_len = database.length();

	int mer_length = 0;

	std::cin >> mer_length;

	std::vector<std::string> most_freq_mers;

	if (mer_length < Large_term_length)
		short_term_processing(mer_length, database, most_freq_mers);
	else
		long_term_processing(mer_length, database, most_freq_mers);
	
	for (const std::string& mer : most_freq_mers)
		std::cout << mer << " ";

	return 0;
}

int short_term_processing(int mer_length, const std::string& database, std::vector<std::string>& most_freq_mers)
{
	int base_len = database.length();

	int num_of_mers = 1;

	for (int i = 0; i < mer_length; ++i)
		num_of_mers *= 4;

	std::vector<int> mers_counter(num_of_mers);

	count_mers(mers_counter, mer_length, database);

	std::vector<int> most_frequent_index;

	count_frequence(num_of_mers, mers_counter, most_frequent_index);

	convert_to_letters(most_freq_mers, most_frequent_index, mer_length);

	std::sort(most_freq_mers.begin(), most_freq_mers.end());

	return 0;
}

int long_term_processing(int mer_length, const std::string& database, std::vector<std::string>& most_freq_mers)
{
	int base_len = database.length();

	std::vector<int> mers_counter(1);

	int max_freq = 0;

	for (int cur_mer_start = 0; cur_mer_start < database.length() - mer_length + 1; ++cur_mer_start)
	{
		int cur_mer_freq = count_mer(database, cur_mer_start, mer_length);

		if (cur_mer_freq > max_freq)
		{
			max_freq = cur_mer_freq;
			mers_counter.resize(0);
			mers_counter.push_back(cur_mer_start);
		}
		else if (cur_mer_freq == max_freq)
			mers_counter.push_back(cur_mer_start);
	}

	std::string most_freq_mer;

	for (int cur_mer_start : mers_counter)
	{
		for (int i = 0; i < mer_length; ++i)
			most_freq_mer += database[cur_mer_start + i];

		most_freq_mers.push_back(most_freq_mer);

		most_freq_mer.resize(0);
	}

	std::sort(most_freq_mers.begin(), most_freq_mers.end());

	most_freq_mers.erase(std::unique(most_freq_mers.begin(), most_freq_mers.end()), most_freq_mers.end());

	return 0;
}

int count_mer(const std::string& database, int mer_start, int mer_length)
{
	int mer_freq = 0;

	for (int cur_mer_start = 0; cur_mer_start < database.length() - mer_length + 1; ++cur_mer_start)
	{
		int match_checker = 0;

		for (int i = 0; i < mer_length; ++i)
		{
			if (database[cur_mer_start + i] == database[mer_start + i])
				++match_checker;
			else
				continue;
		}

		if (match_checker == mer_length)
			++mer_freq;
	}

	return mer_freq;
}

int get_code(char letter)
{
	switch (letter) {
	case 'A':
		return 0;
	case 'C':
		return 1;
	case 'G':
		return 2;
	case 'T':
		return 3;
	}

	return 10;
}

char get_letter(int code)
{
	switch (code) {
	case 0:
		return 'A';
	case 1:
		return 'C';
	case 2:
		return 'G';
	case 3:
		return 'T';
	}

	return 10;
}

int count_mers(std::vector<int>& mers_counter, int k, const std::string& database)
{
	int base_len = database.length();

	for (int mer_start = 0; mer_start < base_len - k + 1; ++mer_start)
	{
		int cur_mer = 0;

		for (int cur_letter = 0; cur_letter < k; ++cur_letter)
		{
			int letters_code = get_code(database[mer_start + cur_letter]);

			for (int i = 0; i < k - cur_letter - 1; ++i)
				letters_code *= 4;

			cur_mer += letters_code;
		}

		if (k < 12)
			++mers_counter[cur_mer];
		else
			mers_counter.push_back(cur_mer);
	}

	return 0;
}

int count_frequence(int num_of_mers, const std::vector<int>& mers_counter, std::vector<int>& most_frequent_index)
{
	int max_freq = 0;

	for (int i = 0; i < num_of_mers; ++i)
		if (mers_counter[i] > max_freq)
			max_freq = mers_counter[i];

	for (int i = 0; i < num_of_mers; ++i)
		if (mers_counter[i] == max_freq)
			most_frequent_index.push_back(i);

	return 0;
}

int convert_to_letters(std::vector<std::string>& most_freq_mers, const std::vector<int>& most_frequent, int mer_length)
{
	std::string one_mer = "";

	for (int cur_mer : most_frequent)
	{
		while (cur_mer)
		{
			one_mer += get_letter(cur_mer % 4);
			cur_mer -= (cur_mer % 4);
			cur_mer /= 4;
		}

		int len_without_first_A = one_mer.length();

		for (int i = 0; i < mer_length - len_without_first_A; ++i)
			one_mer += "A";

		std::reverse(one_mer.begin(), one_mer.end());

		most_freq_mers.push_back(one_mer);

		one_mer.resize(0);
	}

	return 0;
}
