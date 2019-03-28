/*
 * JackknifeAnalyzer.hh
 *
 *  Created on: 22.02.2019
 *      Author: reisinger
 */

#ifndef INCLUDE_JACKKNIFEANALYZER_HH_
#define INCLUDE_JACKKNIFEANALYZER_HH_

#include <map>
#include <vector>

template<typename K, typename T>
class JackknifeAnalyzer {
public:

	/**
	 * Create an empty JackknifeAnalyzer with no datasets. Data is stored with arithmetic type T.
	 * Jackknife-resampled datasets can be added directly or computed from non-resampled datasets or a function.
	 * Jackknife samples and means are stored internally using keys of type K.
	 * Mean and jackknife error of a stored variable can be computed by calling jackknife(...), mu(...), sigma(...) with according key.
	 */
	JackknifeAnalyzer();

	/**
	 * Same as calling
	 * JackknifeAnalyzer();
	 * resample(Xkey, Xsamples);
	 */
	JackknifeAnalyzer(const K& Xkey, const std::vector<T>& Xsamples);

	/**
	 * Store the given jackknife samples of X under the key Xkey for use in further computations.
	 * Must also provide a mean mu_X of X to store, since Xmu cannot be computed from reduced samples if X depends on
	 * other quantities in a non-linear way.
	 * Does nothing if key Xkey already exists. Throws if the number of samples does not match already existing datasets.
	 */
	void add_resampled(const K& Xkey, const std::vector<T>& Xjackknife_samples, const T& mu_X);

	/**
	 * Resamples the given samples of X and stores these and the mean of X under the key Xkey for use in further computations.
	 * If Xkey already exists, does nothing. Throws if the number of samples does not match already existing datasets.
	 */
	void resample(const K& Xkey, const std::vector<T>& Xsamples);

	/**
	 * Computes and stores jackknife samples and mean of a variable F which is a function taking a vector of data type T
	 * of variables with keys in F_arg_keys.
	 * Does nothing if key Fkey already exists.
	 * Throws if Fkey and one or more keys in F_arg_keys do not exist.
	 *
	 * This version of add_function can be used if the number of arguments of F is unknown at compile time, but
	 * requires vector input for the keys and function arguments.
	 */
	template<typename Function>
	void add_function(const K& Fkey, Function F, const std::vector<K>& F_arg_keys);

	/**
	 * Computes and stores jackknife samples and mean of a variable F, which is a function of variables of data type T
	 * with keys F_arg_keys, under the key Fkey.
	 * Does nothing if key Fkey already exists.
	 * Throws if Fkey and one or more keys F_arg_keys do not exist.
	 *
	 * Accepts a function F taking any number of arguments of data type T, which requires knowing the number of arguments at compile time.
	 */
	template<typename Function, typename ... Ks>
	void add_function(const K& Fkey, Function F, const Ks& ... F_arg_keys);

	/**
	 * Removes the variable with key Xkey from the JackknifeAnalyzer.
	 * Does nothing if Xkey does not exist.
	 */
	void remove(const K& Xkey);

	/**
	 * Returns a vector of keys of all variables in the JackknifeAnalyzer.
	 */
	std::vector<K> keys();

	/**
	 * Returns the mean of the variable with key Xkey.
	 * Throw if Xkey does not exist.
	 */
	double mu(const K& Xkey);

	/**
	 * Returns the jackknife error of the variable with key Xkey.
	 * Throws if Xkey does not exist.
	 */
	double sigma(const K& Xkey);

	/**
	 * If a key Xkey does not exist, does nothing and return false, otherwise
	 * assigns the mean / jackknife error of the variable with key Xkey to mu_X / sigma_X and returns true.
	 */
	bool jackknife(const K& Xkey, T& mu_X, T& sigma_X);

	/**
	 * Returns a copy of the jackknife samples of the variable with key Xkey.
	 */
	std::vector<T> samples(const K& Xkey);

private:

	int N_samples;
	void init();
	bool verify_X(const K& Xkey, const std::vector<T>& Xsamples);

	std::map<K, std::vector<T> > Xs_reduced_samples;
	std::map<K, T> Xs_mu;

};

#include <detail/JackknifeAnalyzer.tcc>

#endif /* INCLUDE_JACKKNIFEANALYZER_HH_ */
