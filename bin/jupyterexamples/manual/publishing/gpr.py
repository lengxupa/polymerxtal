import warnings
warnings.simplefilter(action="ignore", category=UserWarning)
import GPy
import numpy as np


class GPR:

    @staticmethod
    def _check_and_regularize_1d_array(x):
        """
        Checks if `x` is a `np.ndarray`. If `x.ndim` is 1,
        then, it turns it into a column vector.
        Otherwise, it leaves it as is.
        """
        x = np.array(x)
        if x.ndim == 1:
            x = x[:, None]
        assert x.ndim == 2
        return x

    def __init__(self, x, y):
        """
        :param x:               The observed inputs (1D `np.ndarray`)
        :param y:               The observed outputs (1D `np.ndarray`).
        """
        self.x = GPR._check_and_regularize_1d_array(x)
        self.y = GPR._check_and_regularize_1d_array(y)

    def run(self, variance=1.,
            length_scale=1.,
            noise_variance=1e-3,
            optimize=False, plot=True):
        """
        Perform 1D regression.
        :param variance:        The signal strength of the square exponential
                            covariance function (positive float).
        :param length_scale:    The length scale of the square
                            exponential covariance function (positive float).
        :param noise_variance:  The noise of the model (non-negative float).
        :param optimize:        If `True` then the model is optimized
                                by maximizing the marginal likelihood
                                with respect to the hyper-parameters (bool).
        :returns:               A dictionary containing the following elements:
                                    + x_eval:       points on which the predictive
                                                    distribution is actually evaluated
                                                    (1D `np.ndarray`)
                                    + y_mu:         the mean of the predictive distribution
                                                    (1D `np.ndarray` of size x_eval.shape[0])
                                    + y_var:        the predictive variance
                                                    (1D `np.ndarray` of size x_eval.shape[0])
                                    + y_025:        the 2.5% lower quantile of the predictive
                                                    distribution of the GP
                                                    (1D `np.ndarray` of size x_eval.shape[0])
                                    + y_975:        the 97.5% lower quantile of the predictive
                                                    distribution of the GP
                                                    (1D `np.ndarray` of size x_eval.shape[0])
                                    + y_s:          samples from the predictive distribution
                                                    of the Gaussian process
                                                    (2D `np.ndarray` of size
                                                     num_samples x x_eval.shape[0])
                                    + model:        the trained gaussian process model
        """
        variance = float(variance)
        length_scale = float(length_scale)
        noise_variance = float(noise_variance)
        optimize = bool(optimize)

        k = GPy.kern.RBF(input_dim=1, lengthscale=length_scale, variance=variance)

        model = GPy.models.GPRegression(self.x, self.y, k)
        model.Gaussian_noise.variance = noise_variance
        if optimize:
            model.optimize(messages=False)
        if plot:
            model.plot()
        else:
            return model

