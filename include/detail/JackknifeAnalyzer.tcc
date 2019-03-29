/*
 * JackknifeAnalyzer.tcc
 *
 *  Created on: 22.02.2019
 *      Author: reisinger
 */

#include <map>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <type_traits>
#include <functional>

#include <helper_functions.hh>
#include <JackknifeAnalyzer.hh>

template<typename K, typename T>
JackknifeAnalyzer<K, T>::JackknifeAnalyzer() {
	init();
}

template<typename K, typename T>
JackknifeAnalyzer<K, T>::JackknifeAnalyzer(const K& Xkey, const std::vector<T>& Xsamples) {
	init();
	resample(Xkey, Xsamples);
}

template<typename K, typename T>
void JackknifeAnalyzer<K, T>::add_resampled(const K& Xkey, const std::vector<T>& Xjackknife_samples, const T& mu_X) {
	if (verify_X(Xkey, Xjackknife_samples)) {
		Xs_reduced_samples[Xkey] = Xjackknife_samples;
		Xs_mu[Xkey] = mu_X;
	}
}

template<typename K, typename T>
void JackknifeAnalyzer<K, T>::resample(const K& Xkey, const std::vector<T>& Xsamples) {
	if (verify_X(Xkey, Xsamples)) {
		T sum_samples = 0;
		for (T d : Xsamples)
			sum_samples += d;
		Xs_mu[Xkey] = sum_samples / (T) N_samples;

		std::vector<T> red_samples;
		for (T d : Xsamples)
			red_samples.push_back((sum_samples - d) / (T) (N_samples - 1));

		Xs_reduced_samples[Xkey] = red_samples;
	}
}

template<typename K, typename T>
template<typename Function>
void JackknifeAnalyzer<K, T>::add_function(const K& Fkey, Function F, const std::vector<K>& F_arg_keys) {
	static_assert(std::is_convertible<Function, std::function<T(std::vector<T>)> >::value,
			"JackknifeAnalyzer::add_function invalid function");

	if (Xs_reduced_samples.count(Fkey) == 0 && Xs_mu.count(Fkey) == 0) {
		std::vector<T> args_mu;
		for (const K& key : F_arg_keys)
			args_mu.push_back(Xs_mu.at(key));
		Xs_mu[Fkey] = F(args_mu);

		std::vector<T> F_jackknife_samples;
		for (int i = 0; i < N_samples; ++i) {
			std::vector<T> args_red_samples;
			for (const K& key : F_arg_keys)
				args_red_samples.push_back(Xs_reduced_samples.at(key).at(i));
			F_jackknife_samples.push_back(F(args_red_samples));
		}
		Xs_reduced_samples[Fkey] = F_jackknife_samples;
	}
}

template<typename K, typename T>
template<typename Function, typename ... Ks>
void JackknifeAnalyzer<K, T>::add_function(const K& Fkey, Function F, const Ks& ... F_arg_keys) {
	static_assert(and_type<std::is_convertible<Ks, K>::value ...>::value,
			"JackknifeAnalyzer::add_function invalid key type");
	static_assert(std::is_convertible<Function, std::function<T(decltype(Xs_mu[F_arg_keys])...)> >::value,
			"JackknifeAnalyzer::add_function invalid function");

	if (Xs_reduced_samples.count(Fkey) == 0 && Xs_mu.count(Fkey) == 0) {
		Xs_mu[Fkey] = F(Xs_mu.at(F_arg_keys)...);

		std::vector<T> F_jackknife_samples;
		for (int i = 0; i < N_samples; ++i)
			F_jackknife_samples.push_back(F(Xs_reduced_samples.at(F_arg_keys).at(i)...));
		Xs_reduced_samples[Fkey] = F_jackknife_samples;
	}
}

template<typename K, typename T>
void JackknifeAnalyzer<K, T>::remove(const K& Xkey) {
	Xs_mu.erase(Xkey);
	Xs_reduced_samples.erase(Xkey);
}

template<typename K, typename T>
std::vector<K> JackknifeAnalyzer<K, T>::keys() {
	std::vector<K> ks;
	for (auto& key_mu : Xs_mu)
		ks.push_back(key_mu.first);
	return ks;
}

template<typename K, typename T>
double JackknifeAnalyzer<K, T>::mu(const K& Xkey) {
	return Xs_mu.at(Xkey);
}

template<typename K, typename T>
double JackknifeAnalyzer<K, T>::sigma(const K& Xkey) {
	double sigma = 0.0;
	for (T d : Xs_reduced_samples.at(Xkey))
		sigma += pow(d - Xs_mu.at(Xkey), (T) 2);
	return sqrt((((T) (N_samples - 1)) / ((T) N_samples)) * sigma);
}

template<typename K, typename T>
bool JackknifeAnalyzer<K, T>::jackknife(const K& Xkey, T& mu_X, T& sigma_X) {
	if (Xs_reduced_samples.count(Xkey) && Xs_mu.count(Xkey)) {
		mu_X = Xs_mu[Xkey];

		sigma_X = 0;
		for (T d : Xs_reduced_samples[Xkey])
			sigma_X += pow(d - mu_X, (T) 2);
		sigma_X = sqrt((((T) (N_samples - 1)) / ((T) N_samples)) * sigma_X);
		return true;
	} else
		return false;
}

template<typename K, typename T>
std::vector<T> JackknifeAnalyzer<K, T>::samples(const K& Xkey) {
	return Xs_reduced_samples.at(Xkey);
}

// ************************************** private **************************************

/**
 * initialize
 */
template<typename K, typename T>
void JackknifeAnalyzer<K, T>::init() {
	static_assert(std::is_arithmetic<T>::value, "JackknifeAnalyzer data type is not arithmetic");
	N_samples = 0;
}

/**
 * Ensure equal and > 1 number of samples for all variables, otherwise an exception is thrown.
 * Return true if number of samples is equal to existing data and data with the key Xkey does not exist, false otherwise.
 */
template<typename K, typename T>
bool JackknifeAnalyzer<K, T>::verify_X(const K& Xkey, const std::vector<T>& Xsamples) {
	if (Xs_reduced_samples.count(Xkey) || Xs_mu.count(Xkey))
		return false;

	if (N_samples == 0) {
		if (Xsamples.size() > 1)
			N_samples = Xsamples.size();
		else
			throw std::runtime_error("trying to add dataset with less than 2 samples.");
	} else if (Xsamples.size() != N_samples)
		throw std::runtime_error("trying to add dataset with different number of samples than already existing ones.");

	return true;
}
