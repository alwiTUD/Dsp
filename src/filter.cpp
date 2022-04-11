#include "filter.h"

#include "dsp.h"

template <class T>
std::vector<T> dsp::filter::filter(std::vector<T> b,
	std::vector<T> a, const std::vector<T>& x)
{
	// Normalize filter coefficients
	auto a0 = a[0];
	for (auto& bi : b)
	{
		bi /= a0;
	}
	for (auto& ai : a)
	{
		ai /= a0;
	}

	if (a.size() == 1)
	{
		// This is the simple case without a denominator in the transfer function
		return dsp::convolve(x, b);
	}

	// This is the more complicated case with both numerator and denominator

	std::vector<T> y(x.size());  // The filtered output
	auto n = std::max(b.size(), a.size());
	a.resize(n);
	b.resize(n);
	std::vector<T> w_m(n + 1, 0);  // w(m)
	std::vector<T> w_mm1(n + 1, 0);  // w(m-1)
	for (size_t m = 0; m < x.size(); ++m)
	{
		y[m] = b[0] * x[m] + w_mm1[1];
		for (size_t k = n-1; k > 0; --k)
		{
			w_m[k] = b[k] * x[m] + w_mm1[k + 1] /* zero for k = n-1 */ - a[k] * y[m];
		}		
		w_mm1 = w_m;
	}

	return y;
}

template<class T>
std::vector<T> dsp::filter::lpc(const std::vector<T>& x, unsigned N)
{
	std::vector<T> r;
	auto numCoefficients = N + 1;
	r.resize(numCoefficients);

	for (size_t i = 0; i <= N; i++)
	{
		r[i] = 0.0;
		for (size_t j = 0; j < x.size() - i; ++j)
		{
			r[i] += x[j] * x[j + i];
		}
	}

	// Levinson-Durbin
	T E = r[0];
	std::vector<T> alpha(numCoefficients, 0);
	alpha[0] = 1.0;
	std::vector<T> beta(numCoefficients, 0);
	std::vector<T> z(numCoefficients, 0);

	for (unsigned p = 1; p <= N; ++p)
	{
		T q = 0.0;
		for (unsigned i = 0; i < p; ++i) { q += alpha[i] * r[p - i]; }
		if (E == 0.0) { E = static_cast<T>(0.0001); }
		z[p] = -q / E;
		alpha[p] = 0.0;
		for (unsigned i = 0; i <= p; ++i) { beta[i] = alpha[i] + z[p] * alpha[p - i]; }
		for (unsigned i = 0; i <= p; ++i) { alpha[i] = beta[i]; }

		E = E * (static_cast<T>(1.0) - z[p] * z[p]);
	}

	std::vector<T> coeff;
	coeff.resize(N + 1);
	coeff[0] = 1;
	for (unsigned i = 1; i <= N; ++i) { coeff[i] = -alpha[i]; }

	return coeff;
}

template <class T>
std::vector<T> dsp::filter::medianfilter(const std::vector<T>& x, size_t kernel_size)
{
	if (kernel_size % 2 != 1) { throw std::runtime_error("Kernel size should be odd!"); }

	// Pad vector with zeros at the edges
	std::vector<T> padded_x = pad(x, { kernel_size / 2, kernel_size / 2 });
	std::vector<T> out;

	// Find median in each kernel
	for (size_t k = 0; k < x.size(); ++k)
	{
		out.push_back(median<T>(padded_x.begin() + k, padded_x.begin() + k + kernel_size));
	}
	
	return out;
}

template <class T>
dsp::Signal<T> dsp::filter::medianfilter(const Signal<T>& x, typename Signal<T>::size_type kernel_size)
{
	return Signal(x.getSamplingRate_Hz(), medianfilter<T>(x.getSamples(), kernel_size));
}

template<class T>
std::vector<std::vector<T>> dsp::filter::melFilterbank(const std::vector<std::vector<T>>& spectrogram, int numFilters, std::pair<double, double> cutoffFreqs_Hz, double samplingRate)
{
	// Find the border frequencies of the mel filter bank
	std::vector<double> borderFreqs_Hz;
	borderFreqs_Hz.reserve(static_cast<size_t>(numFilters + 2.0));
	borderFreqs_Hz = findBorderFreqs_Hz(numFilters, cutoffFreqs_Hz);

	// Find the heights for each of the filters of the mel filter bank
	std::vector<double> filterHeights; 
	filterHeights.reserve(numFilters);
	filterHeights = findFilterHeights(borderFreqs_Hz);

	// Find the weights of each filter of the mel filter bank in relation to the frequency bins of the input spectrogram
	const int numPositiveFrequencyBins = static_cast<int>(spectrogram.size());
	std::vector<std::vector<double>> filterWeights(numFilters, std::vector<double>(numPositiveFrequencyBins));
	filterWeights = findFilterWeights(borderFreqs_Hz, filterHeights, samplingRate, numPositiveFrequencyBins);

	// Applies the mel filter bank to the spectrogram, resulting in a mel spectrogram
	std::vector<std::vector<double>> melSpectrogram;
	melSpectrogram.resize(numFilters, std::vector<double>(static_cast<int>(spectrogram[0].size())));
	melSpectrogram = filterAndSum(filterWeights, spectrogram);

	// Applies log to the mel spectrogram, resulting in a mel log spectrogram
	std::vector<std::vector<double>> melLogSpectrogram;
	melLogSpectrogram.resize(numFilters, std::vector<double>(static_cast<int>(spectrogram[0].size())));
	melLogSpectrogram = linearToLog(melSpectrogram);
	
	return melLogSpectrogram;
}

std::vector<double> dsp::filter::findBorderFreqs_Hz(int numFilters, std::pair<double, double> cutoffFreqs_Hz)
{
	// Convert border frequencies from Hz to Mel
	double fmin_Mel = dsp::convert::hz2mel(cutoffFreqs_Hz.first);
	double fmax_Mel = dsp::convert::hz2mel(cutoffFreqs_Hz.second);
	double freqDiff_Mel = fmax_Mel - fmin_Mel;
	// Calculate width of individual filters (all same width in Mel domain)
	double freqDelta_Mel = freqDiff_Mel / static_cast<double>(numFilters + 1.0);

	std::vector<double> borderFreqs;
	borderFreqs.reserve(static_cast<size_t>(numFilters + 2));
	// Fill vector borderFreqs with border frequencies of individual filters in Hz
	for (int filterIdx = 0; filterIdx < numFilters + 2; filterIdx++)
	{
		borderFreqs.push_back(dsp::convert::mel2hz(fmin_Mel + filterIdx * freqDelta_Mel));
	}
	return borderFreqs;
}

std::vector<double> dsp::filter::findFilterHeights(std::vector<double> borderFreqs_Hz)
{
	std::vector<double> filterHeights;
	filterHeights.reserve(borderFreqs_Hz.size() - 2);

	// Calculates the height of each filter. The area of each triangular filter should be equal to 1.
	// (f_upper - f_lower)*h/2 = 1, i.e., h = 2/(f_upper - f_lower)
	for (int borderFreqIdx = 0; borderFreqIdx < borderFreqs_Hz.size() - 2; ++borderFreqIdx)
	{
		filterHeights.push_back(2.0 / (borderFreqs_Hz[borderFreqIdx + 2] - borderFreqs_Hz[borderFreqIdx]));
	}
	
	return filterHeights;
}

std::vector<std::vector<double>> dsp::filter::findFilterWeights(std::vector<double>& borderFreqs_Hz, std::vector<double>& filterHeights, double samplingRate, int numPositiveFrequencyBins)
{
	const int numFilters = static_cast<int>(filterHeights.size());
	std::vector<std::vector<double>> wFilt(numFilters, std::vector<double>(numPositiveFrequencyBins));

	for (int melFilt = 0; melFilt < numFilters; ++melFilt)
	{
		// Border frequencies of each filter of the mel filter bank
		const double fLower = borderFreqs_Hz[melFilt];
		const double fUpper = borderFreqs_Hz[melFilt + 2];
		const double fCenter = borderFreqs_Hz[melFilt + 1];
		const double filterHeight = filterHeights[melFilt];
		for (int frequencyBinIdx = 0; frequencyBinIdx < numPositiveFrequencyBins; ++frequencyBinIdx)
		{
			// Calculates the current frequency in Hz (between 0 and the Nyquist frequency)
			const double freq = samplingRate * (static_cast<double>(frequencyBinIdx) / (2*(numPositiveFrequencyBins - 1)));
			if (freq >= fLower && freq < fCenter)
			{
				// Equation for the filter weight: left side
				wFilt[melFilt][frequencyBinIdx] = (freq - fLower) * (filterHeight / (fCenter - fLower));
			}
			else if (freq <= fUpper && freq >= fCenter)
			{
				// Equation for the filter weight: right side
				wFilt[melFilt][frequencyBinIdx] = filterHeight - (freq - fCenter) * (filterHeight / (fUpper - fCenter));
			}
			else
			{
				// Zero (outside of the filter's scope)
				wFilt[melFilt][frequencyBinIdx] = 0.0;
			}
		}
	}	
	
	return wFilt;
}

std::vector<std::vector<double>> dsp::filter::filterAndSum(const std::vector<std::vector<double>>& wFilt, const std::vector<std::vector<double>>& spectrogram)
{
	int numFilt = static_cast<int>(wFilt.size());
	int numberOfFrames = static_cast<int>(spectrogram[0].size());
	int numberOfPositiveFrequencyBins = static_cast<int>(wFilt[0].size());

	// Initializes a vector filled with zeros to store the values of the Mel-spectrogram
	std::vector<std::vector<double>> melSpectrogram;
	melSpectrogram.resize(numFilt, std::vector<double>(numberOfFrames, 0.0));
	
	// Applies mel filter bank to the original spectrogram
	for (int melFilt = 0; melFilt < numFilt; ++melFilt) {
		for (int frame = 0; frame < numberOfFrames; ++frame) {
			for (int frequencyBin = 0; frequencyBin < numberOfPositiveFrequencyBins; ++frequencyBin) {
				melSpectrogram[melFilt][frame] += wFilt[melFilt][frequencyBin] * spectrogram[frequencyBin][frame];
			}
		}
	}

	return melSpectrogram;
}

std::vector<std::vector<double>> dsp::filter::linearToLog(const std::vector<std::vector<double>>& linMelSpectrogram) {
	std::vector<std::vector<double>> logMelSpectrogram;
	logMelSpectrogram.reserve(linMelSpectrogram.size());

	std::vector<double> logMelSpectrogramFrame;
	logMelSpectrogramFrame.resize(linMelSpectrogram[0].size());

	// Dimensions of variables
	// linMelSpectrogram = [numFilter x numFrames]
	// filterOutput = [1 x numFrames]
	// logMelSpectrogramFrame = [1 x numFrames]
	// logMelSpectrogram = [numFilter x numFrames]

	// Applies log to each value of the mel spectrogram (a very small epsilon value is added to prevent log(0))
	for (auto& filterOutput : linMelSpectrogram) {
		std::transform(filterOutput.begin(), filterOutput.end(), logMelSpectrogramFrame.begin(), [](double val) -> double {return std::log10(val + std::numeric_limits<double>::epsilon()); }); // EPSILON equals the one from MATLAB
		logMelSpectrogram.push_back(logMelSpectrogramFrame);
	}

	return logMelSpectrogram;
}


// Explicit template instantiation
template std::vector<float> dsp::filter::filter(std::vector<float> b,
	std::vector<float> a, const std::vector<float>& x);
template std::vector<double> dsp::filter::filter(std::vector<double> b,
	std::vector<double> a, const std::vector<double>& x);
template std::vector<long double> dsp::filter::filter(std::vector<long double> b,
	std::vector<long double> a, const std::vector<long double>& x);

template std::vector<float> dsp::filter::lpc(const std::vector<float>& x, unsigned N);
template std::vector<double> dsp::filter::lpc(const std::vector<double>& x, unsigned N);
template std::vector<long double> dsp::filter::lpc(const std::vector<long double>& x, unsigned N);

template std::vector<float> dsp::filter::medianfilter(const std::vector<float>& x, size_t kernel_size);
template std::vector<double> dsp::filter::medianfilter(const std::vector<double>& x, size_t kernel_size);
template std::vector<long double> dsp::filter::medianfilter(const std::vector<long double>& x, size_t kernel_size);

template dsp::Signal<float> dsp::filter::medianfilter(const dsp::Signal<float>& x, Signal<float>::size_type kernel_size);
template dsp::Signal<double> dsp::filter::medianfilter(const dsp::Signal<double>& x, Signal<double>::size_type kernel_size);
template dsp::Signal<long double> dsp::filter::medianfilter(const dsp::Signal<long double>& x, Signal<long double>::size_type kernel_size);

template std::vector<std::vector<double>> dsp::filter::melFilterbank(const std::vector<std::vector<double>>& spectrogram, int numFilters, std::pair<double, double> cutoffFreqs_Hz, double samplingRate);