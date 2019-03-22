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
	 * Mean and jackknife error of a stored variable can be computed by calling jackknife(...) with according key.
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
	 * Must also provide a mean mu_X of X to store, since Xmu cannot be computed from reduced samples if X depends on other quantities in a non-linear way.
	 * Does nothing if key XKey already exists. Throws if the number of samples does not match already existing datasets.
	 */
	void add_resampled(const K& Xkey, const std::vector<T>& Xjackknife_samples, const T& mu_X);

	/**
	 * Resample the given samples of X and store these and the mean of X under the key Xkey for use in further computations.
	 * If Xkey already exists, does nothing. Throws if the number of samples does not match already existing datasets.
	 */
	void resample(const K& Xkey, const std::vector<T>& Xsamples);

	/**
	 * Computes and stores jackknife samples and mean of a variable F which is a function of variables X with keys Xkeys.
	 * Does nothing if key FKey already exists.
	 * Throws if Fkey and one or more keys Xkeys do not exist.
	 */
	template<typename ... Ks, typename ... Ts>
	void add_function(const K& Fkey, T (*F)(Ts ...), const Ks& ... Xkeys);

	/**
	 * If a key Xkey does not exist, do nothing and return false, otherwise
	 * get mean / compute jackknife error of the variable with key Xkey and return true.
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
