#include <map>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cmath>
#include <type_traits>
#include <functional>

#include <helper_functions.hh>
#include <JackknifeAnalyzer.hh>

namespace de_uni_frankfurt_itp {
namespace reisinger {
namespace jackknife_analyzer_0219 {

template<typename K, typename T>
JackknifeAnalyzer<K, T>::JackknifeAnalyzer(std::size_t bin_size) :
		N_bins { 0 }, bin_size { bin_size } {

	static_assert(std::is_arithmetic<T>::value, "JackknifeAnalyzer data type is not arithmetic");
}

template<typename K, typename T>
JackknifeAnalyzer<K, T>::JackknifeAnalyzer(const K& Xkey, const std::vector<T>& Xsamples, std::size_t bin_size) :
		JackknifeAnalyzer<K, T> { bin_size } {
	resample(Xkey, Xsamples);
}

template<typename K, typename T>
void JackknifeAnalyzer<K, T>::add_resampled(const K& Xkey, const std::vector<T>& Xjackknife_samples, const T& mu_X) {
	if (Xs_reduced_samples.count(Xkey) == 0 && Xs_mu.count(Xkey) == 0) {
		init_or_verify_N(Xjackknife_samples, true);

		Xs_reduced_samples[Xkey] = Xjackknife_samples;
		Xs_mu[Xkey] = mu_X;
	}
}

template<typename K, typename T>
void JackknifeAnalyzer<K, T>::resample(const K& Xkey, const std::vector<T>& Xsamples) {
	if (Xs_reduced_samples.count(Xkey) == 0 && Xs_mu.count(Xkey) == 0) {
		init_or_verify_N(Xsamples, false);

		T sum_samples = 0;
		for (T d : Xsamples)
			sum_samples += d;
		Xs_mu[Xkey] = sum_samples / static_cast<T>(Xsamples.size());

		std::vector<T> red_samples;
		red_samples.reserve(N_bins);

		for (std::size_t b = 0; b < N_bins; ++b) {
			T red_sample = sum_samples;
			const auto next_bin_first_sample = (b + 1) * bin_size;
			for (std::size_t i = b * bin_size; i < next_bin_first_sample; ++i)
				red_sample -= Xsamples[i];

			red_samples.push_back(red_sample / static_cast<T>(Xsamples.size() - bin_size));
		}

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
		for (std::size_t i = 0; i < N_bins; ++i) {
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
	static_assert(tools::helper::and_type<std::is_convertible<Ks, K>::value ...>::value,
			"JackknifeAnalyzer::add_function invalid key type");
	static_assert(std::is_convertible<Function, std::function<T(decltype(Xs_mu[F_arg_keys])...)> >::value,
			"JackknifeAnalyzer::add_function invalid function");

	if (Xs_reduced_samples.count(Fkey) == 0 && Xs_mu.count(Fkey) == 0) {
		const T F_mu = F(Xs_mu.at(F_arg_keys)...); // before [] to be exception safe
		Xs_mu[Fkey] = F_mu;

		std::vector<T> F_jackknife_samples;
		for (std::size_t i = 0; i < N_bins; ++i)
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
std::vector<K> JackknifeAnalyzer<K, T>::keys() const {
	std::vector<K> ks;
	for (const auto& key_mu : Xs_mu)
		ks.push_back(key_mu.first);
	return ks;
}

template<typename K, typename T>
T JackknifeAnalyzer<K, T>::mu(const K& Xkey) const {
	return Xs_mu.at(Xkey);
}

template<typename K, typename T>
T JackknifeAnalyzer<K, T>::sigma(const K& Xkey) const {
	double sigma = 0.0;
	for (const T& d : Xs_reduced_samples.at(Xkey))
		sigma += pow(d - Xs_mu.at(Xkey), (T) 2);
	return sqrt((((T) (N_bins - 1)) / ((T) N_bins)) * sigma);
}

template<typename K, typename T>
bool JackknifeAnalyzer<K, T>::jackknife(const K& Xkey, T& mu_X, T& sigma_X) const {
	if (Xs_reduced_samples.count(Xkey) && Xs_mu.count(Xkey)) {
		mu_X = Xs_mu.at(Xkey);

		sigma_X = 0;
		for (const T& d : Xs_reduced_samples.at(Xkey))
			sigma_X += pow(d - mu_X, (T) 2);
		sigma_X = sqrt((((T) (N_bins - 1)) / ((T) N_bins)) * sigma_X);
		return true;
	} else
		return false;
}

template<typename K, typename T>
std::vector<T> JackknifeAnalyzer<K, T>::samples(const K& Xkey) const {
	return Xs_reduced_samples.at(Xkey);
}

// ************************************** private **************************************

template<typename K, typename T>
bool JackknifeAnalyzer<K, T>::init_or_verify_N(const std::vector<T>& Xsamples, bool binned) {
	const auto num_bins = Xsamples.size() / (binned ? 1 : bin_size);

	if (N_bins == 0) {
		if (num_bins > 1)
			N_bins = num_bins;
		else
			throw std::runtime_error("trying to add dataset with less than 2 bins.");
	} else if (num_bins != N_bins)
		throw std::runtime_error("trying to add dataset with different number of bins than already existing ones.");

	return true;
}

}
}
}
