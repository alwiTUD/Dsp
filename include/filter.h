#pragma once

#include <vector>

#include "Signal.h"

/// @brief Functions to apply and design digital filters
namespace dsp::filter
{
	/// @brief Filters the input data x using a rational transfer function defined by the numerator and denominator coefficients b and a.
	/// @tparam T Data type of the samples
	/// @param b Numerator coefficients
	/// @param a Denominator coefficients
	/// @param x Input signal
	/// @return Filtered signal
	template<class T>
	std::vector<T> filter(std::vector<T> b, std::vector<T> a, 
		const std::vector<T>& x);
	
	/// @brief Returns linear prediction filter coefficients
	/// @tparam T Data type of the samples
	/// @param x Input vector
	/// @param N Order of the linear predictor
	/// @return Linear predictor coefficients
	template<class T>
	std::vector<T> lpc(const std::vector<T>& x, unsigned N);

	/// @brief Perform a median filter on a vector.
	/// @tparam T Data type of the samples
	/// @param x Input vector
	/// @param kernel_size Size of the median filter window. Should be odd! Default: 3
	/// @return Vector containing the median filtered samples (same size as input).
	template<class T>
	std::vector<T> medianfilter(const std::vector<T>& x, size_t kernel_size = 3);

	/// @brief Perform a median filter on a Signal.
	/// @tparam T Data type of the samples
	/// @param x Input Signal
	/// @param kernel_size Size of the median filter window. Should be odd! Default: 3
	/// @return Signal containing the median filtered samples (same size as input).
	template<class T>
	dsp::Signal<T> medianfilter(const Signal<T>& x, typename Signal<T>::size_type kernel_size = 3);

	/// @brief Applies a mel filterbank to a magnitude spectrogram.
	/// @tparam T Data type of the spectrogram
	/// @param spectrogram Input spectrogram
	/// @param numFilters Number of triangular filters
	/// @param cutoffFreqs_Hz Lower and upper border frequencies of the mel filterbank [Hz]
	/// @param samplingRate Sampling rate [Hz]
	/// @return MelSpectrogram
	template<class T>
	std::vector<std::vector<T>> melFilterbank(const std::vector<std::vector<T>>& spectrogram, int numFilters, std::pair<double, double> cutoffFreqs_Hz, double samplingRate);

	/// @brief Find border frequencies for given number of filters and cutoff frequencies [Hz]
	/// Border frequencies will be equally spaced in mel domain
	/// @param numFilters Number of triangular filters
	/// @param cutoffFreqs_Hz Lower and upper border frequencies of the mel filterbank [Hz]
	/// @return Border frequencies of the triangular filters [Hz]
	std::vector<double> findBorderFreqs_Hz(int numFilters, std::pair<double, double> cutoffFreqs_Hz);

	/// @brief Find heights of the triangular filters so that all filters will have an area equal to 1
	/// @param borderFreqs_Hz Border frequencies [Hz]
	/// @return Filter heights
	std::vector<double> findFilterHeights(std::vector<double> borderFreqs_Hz);

	/// @brief Find the weights of each frequency bin for each filter of the Mel filter bank
	/// @param borderFreqs_Hz Border frequencies of the triangular filters [Hz]
	/// @param filterHeights Filter heights of the mel filter bank
	/// @param samplingRate Sampling rate [Hz]
	/// @param numPositiveFrequencyBins equal to nFFT/2 + 1
	/// @return A vector of vectors. The "matrix" contains the filter weights of the mel filter bank [numFilters x numPositiveFrequencyBins]
	std::vector<std::vector<double>> findFilterWeights(std::vector<double>& borderFreqs_Hz, std::vector<double>& filterHeights, double samplingRate, int numPositiveFrequencyBins);

	/// @brief Applies the mel filterbank to the spectrogram
	/// @param wFilt The filter weights of the mel filter bank
	/// @param spectrogram The spectrogram
	/// @return The values of the energy within each mel filter for each time frame
	std::vector<std::vector<double>> filterAndSum(std::vector<std::vector<double>>& wFilt, std::vector<std::vector<double>>& spectrogram);
}