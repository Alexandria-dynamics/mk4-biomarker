"""
MK4 Engine - Core frequency-domain biomarker algorithm.

MK4 (Molecular Kinetics 4th-order) uses frequency decomposition
of RNA expression data to compute chaos and coherence scores
that distinguish disease from healthy states.
"""

import numpy as np
from scipy import fft
from scipy import stats


def mk4_transform(expression_vector):
    """
    Apply MK4 frequency-domain transform to a single sample.

    Computes the power spectrum of gene expression values
    and extracts frequency-domain features.

    Parameters
    ----------
    expression_vector : np.ndarray
        1D array of gene expression values for one sample.

    Returns
    -------
    dict
        'power_spectrum': power spectral density,
        'frequencies': corresponding frequencies,
        'spectral_entropy': entropy of normalized spectrum,
        'dominant_freq': frequency with highest power.
    """
    expression_vector = np.asarray(expression_vector, dtype=float)

    # Remove mean (DC component)
    centered = expression_vector - np.mean(expression_vector)

    # FFT
    n = len(centered)
    fft_vals = fft.rfft(centered)
    power = np.abs(fft_vals) ** 2
    freqs = fft.rfftfreq(n)

    # Normalize power spectrum
    total_power = np.sum(power)
    if total_power > 0:
        norm_power = power / total_power
    else:
        norm_power = np.ones_like(power) / len(power)

    # Spectral entropy
    nonzero = norm_power[norm_power > 0]
    spectral_entropy = -np.sum(nonzero * np.log2(nonzero))

    # Dominant frequency
    dominant_idx = np.argmax(power[1:]) + 1  # skip DC
    dominant_freq = freqs[dominant_idx]

    return {
        'power_spectrum': power,
        'frequencies': freqs,
        'spectral_entropy': spectral_entropy,
        'dominant_freq': dominant_freq,
    }


def calculate_chaos(expression_vector):
    """
    Calculate chaos score for a single sample.

    The chaos score quantifies the disorder/irregularity
    in gene expression patterns using spectral entropy.
    Higher chaos = more disordered expression = potential disease.

    Parameters
    ----------
    expression_vector : np.ndarray
        1D array of gene expression values.

    Returns
    -------
    float
        Chaos score (spectral entropy, typically 5-15 range).
    """
    result = mk4_transform(expression_vector)
    return result['spectral_entropy']


def calculate_coherence(expression_vector):
    """
    Calculate coherence score for a single sample.

    Coherence measures how organized the expression pattern is.
    It is the inverse complement of chaos: coherence = max_entropy - chaos.
    Higher coherence = more organized = typically healthy.

    Parameters
    ----------
    expression_vector : np.ndarray
        1D array of gene expression values.

    Returns
    -------
    float
        Coherence score.
    """
    expression_vector = np.asarray(expression_vector, dtype=float)
    n_freqs = len(fft.rfft(expression_vector))
    max_entropy = np.log2(n_freqs) if n_freqs > 0 else 1.0
    chaos = calculate_chaos(expression_vector)
    return max_entropy - chaos


class MK4Analyzer:
    """
    Main MK4 analysis class.

    Provides a high-level interface for analyzing expression datasets
    using the MK4 frequency-domain method.

    Parameters
    ----------
    normalize : bool
        Whether to z-score normalize each sample before analysis.

    Examples
    --------
    >>> analyzer = MK4Analyzer()
    >>> results = analyzer.analyze_dataset(expression_matrix, labels)
    >>> print(results['chaos'].mean())
    """

    def __init__(self, normalize=True):
        self.normalize = normalize
        self._fitted = False

    def analyze_sample(self, expression_vector):
        """
        Analyze a single sample.

        Parameters
        ----------
        expression_vector : np.ndarray
            1D array of gene expression values.

        Returns
        -------
        dict
            'chaos': chaos score,
            'coherence': coherence score,
            'transform': full MK4 transform result.
        """
        vec = np.asarray(expression_vector, dtype=float)
        if self.normalize:
            std = np.std(vec)
            if std > 0:
                vec = (vec - np.mean(vec)) / std

        transform = mk4_transform(vec)
        chaos = transform['spectral_entropy']
        n_freqs = len(transform['power_spectrum'])
        max_entropy = np.log2(n_freqs) if n_freqs > 0 else 1.0
        coherence = max_entropy - chaos

        return {
            'chaos': chaos,
            'coherence': coherence,
            'transform': transform,
        }

    def analyze_dataset(self, expression_matrix, labels=None):
        """
        Analyze an entire expression matrix.

        Parameters
        ----------
        expression_matrix : np.ndarray
            2D array (n_samples, n_genes).
        labels : np.ndarray, optional
            Sample labels (0=control, 1=disease).

        Returns
        -------
        dict
            'chaos': array of chaos scores,
            'coherence': array of coherence scores,
            'labels': labels (if provided).
        """
        expression_matrix = np.asarray(expression_matrix, dtype=float)
        n_samples = expression_matrix.shape[0]

        chaos_scores = np.zeros(n_samples)
        coherence_scores = np.zeros(n_samples)

        for i in range(n_samples):
            result = self.analyze_sample(expression_matrix[i])
            chaos_scores[i] = result['chaos']
            coherence_scores[i] = result['coherence']

        self._fitted = True

        output = {
            'chaos': chaos_scores,
            'coherence': coherence_scores,
        }
        if labels is not None:
            output['labels'] = np.asarray(labels)

        return output
