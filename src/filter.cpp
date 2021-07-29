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
	std::vector<double> borderFreqs_Hz;
	borderFreqs_Hz.reserve(numFilters + 2);
	
	borderFreqs_Hz = findBorderFreqs_Hz(numFilters, cutoffFreqs_Hz);

	std::vector<double> filterHeights; 
	filterHeights.reserve(numFilters);
	filterHeights = findFilterHeights(borderFreqs_Hz);

	// TODO: Allocation of FFT bins to corresponding filter inside of the filter bank
	// TODO: Write function findFilterWeights()
	// TODO: Test filter weights with matlab benchmark: designAuditoryFilterBank(16000, "FFTLength", 1024, "NumBands", 15, "FrequencyRange", [200 3700]); 
	// TODO: Weighting and summation for each filter ---> filterAndSum();
	// TODO: Logarithmierung

	const int numPositiveFrequencyBins = static_cast<int>(spectrogram.size()[0]);
	std::vector<std::vector<double>> wFilt(numFilters, std::vector<double>(numPositiveFrequencyBins));
	wFilt = findFilterWeights(borderFreqs_Hz, filterHeights, samplingRate, numPositiveFrequencyBins);

	std::vector<std::vector<double>> melSpectrogram;
	melSpectrogram.resize(numFilters, std::vector<double>(static_cast<int>(spectrogram[0].size())));
	melSpectrogram = filterAndSum(wFilt, spectrogram);

	//melSpectrogram = filterAndSum(x, y, z)
	//melLogSpectrogram = melLog(x, y, z)
	
	return std::vector<std::vector<T>>();
}

std::vector<double> dsp::filter::findBorderFreqs_Hz(int numFilters, std::pair<double, double> cutoffFreqs_Hz)
{
	double fmin_Mel = dsp::convert::hz2mel(cutoffFreqs_Hz.first);
	double fmax_Mel = dsp::convert::hz2mel(cutoffFreqs_Hz.second);
	double freqDiff_Mel = fmax_Mel - fmin_Mel;
	double freqDelta_Mel = freqDiff_Mel / (numFilters + 1);

	std::vector<double> borderFreqs;
	borderFreqs.reserve(numFilters + 2);
	for (int ii = 0; ii < numFilters + 2; ii++)
	{
		if (dsp::convert::mel2hz(fmin_Mel + ii * freqDelta_Mel) < cutoffFreqs_Hz.first)
		{
			borderFreqs.push_back(cutoffFreqs_Hz.first);
		}
		else if (dsp::convert::mel2hz(fmin_Mel + ii * freqDelta_Mel) > cutoffFreqs_Hz.second)
		{
			borderFreqs.push_back(cutoffFreqs_Hz.second);
		}
		else
		{
			borderFreqs.push_back(dsp::convert::mel2hz(fmin_Mel + ii * freqDelta_Mel));
		}
	}
	return borderFreqs;
}

std::vector<double> dsp::filter::findFilterHeights(std::vector<double> borderFreqs_Hz)
{
	std::vector<double> filterHeights;
	filterHeights.reserve(borderFreqs_Hz.size() - 2);

	for (int ii = 0; ii < borderFreqs_Hz.size() - 2; ++ii)
	{
		filterHeights.push_back(2.0 / (borderFreqs_Hz[ii + 2] - borderFreqs_Hz[ii]));
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

	std::vector<std::vector<double>> melSpectrogram;
	melSpectrogram.resize(numFilt, std::vector<double>(numberOfFrames, 0.0));

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

	for (auto& frame : linMelSpectrogram) {
		std::transform(frame.begin(), frame.end(), logMelSpectrogramFrame.begin(), [](double val) -> double {return std::log10(val + std::numeric_limits<double>::epsilon()); }); // EPSILON equals the one from MATLAB
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