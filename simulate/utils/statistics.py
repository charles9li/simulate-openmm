from scipy import stats
import numpy as np


class StatisticalInformation(object):

    def __init__(self, data, warmup=True):
        self.data = np.array(data)
        self.warmup_data = np.empty([0])
        if warmup:
            self._detect_warmup_mser_5()
        self._compute_statistics()

    def _detect_warmup_mser_5(self):
        self._detect_warmup_mser_m(5)

    def _detect_warmup_mser_m(self, m):
        # Pre-block the data in chunks of 5 steps, and pad data if necessary
        N = len(self.data)
        if N % m != 0:
            data_padded = np.pad(self.data, (0, m - N % m), mode='mean', stat_length=N % m)
            data_reshaped = data_padded.reshape(int(N/m) + 1, m)
        else:
            data_reshaped = self.data.reshape(int(N/m), m)
        data_block_averaged = np.mean(data_reshaped, axis=1)

        # Compute standard errors of mean
        std_err_mean_list = []
        full_length = len(data_block_averaged)
        while True:
            std_err_mean = stats.sem(data_block_averaged)
            std_err_mean_list.append(std_err_mean)
            N_remaining = len(data_block_averaged)
            if N_remaining <= 0.25*full_length or N_remaining < 5:
                break
            data_block_averaged = data_block_averaged[1:]

        # Find index that minimizes the variance
        num_blocks_remove = np.argmin(std_err_mean_list)
        index_remove = num_blocks_remove*m

        if index_remove > 0:
            warmup_data = self.data[:index_remove - 1]
        else:
            warmup_data = np.empty([0])
        production_data = self.data[index_remove:]

        self.data = production_data
        self.warmup_data = warmup_data

    def _compute_statistics(self):

        # Calculate statistics
        nobs, minmax, mean, variance, skewness, kurtosis = stats.describe(self.data)

        # Standard error of mean
        sem = stats.sem(self.data)

        # Compute the autocorrelation function
        data_shifted = self.data - mean
        correlation = np.correlate(data_shifted, data_shifted, mode='same')/variance
        autocorrelation = correlation[int(correlation.size/2):]
        autocorrelation = autocorrelation/np.arange(nobs - 1, nobs - 1 - autocorrelation.size, -1)

        # Choose where to cut off the autocorrelation time sum
        cutoff = autocorrelation.size
        j = 0
        while j < cutoff:
            if autocorrelation[j] < np.sqrt(2./(nobs - j)):
                cutoff = np.minimum(cutoff, 5*j)
            j += 1

        # Compute correlation time
        kappa = 1.0 + 2.0*np.sum(autocorrelation[1:int(2.0*cutoff/5.0)])

        # Update the standard error of the mean for a correlation correction
        sem_cc = sem*np.sqrt(kappa)

        self.nobs = nobs
        self.min, self.max = minmax
        self.mean = mean
        self.sem = sem
        self.sem_cc = sem_cc
        self.kappa = kappa
        self.variance = variance
        self.autocorrelation = autocorrelation
