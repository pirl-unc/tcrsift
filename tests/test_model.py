# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests for the ridge-regularized logistic frequency model."""

import math

import numpy as np
import pytest

from tcrsift.model import (
    LogisticFrequencyModel,
    analytical_logistic_threshold,
    fit_frequency_logistic_model,
)


class TestLogisticFrequencyModel:
    """Tests for LogisticFrequencyModel and fit_frequency_logistic_model."""

    def test_fit_positive_slope_and_predict(self):
        x = np.linspace(0.01, 0.9, 60)
        y = (x > 0.3).astype(float)
        model = fit_frequency_logistic_model(x, y)
        assert isinstance(model, LogisticFrequencyModel)
        assert model.converged
        assert float(model.params[0]) > 0.0
        probs = model.predict([0.0, 0.5, 1.0])
        assert probs[0] == pytest.approx(0.5)
        assert probs[1] < probs[2]

    def test_fit_validation_errors(self):
        with pytest.raises(ValueError, match="empty"):
            fit_frequency_logistic_model([], [])
        with pytest.raises(ValueError, match="same shape"):
            fit_frequency_logistic_model([0.1, 0.2], [1.0])
        with pytest.raises(ValueError, match="finite"):
            fit_frequency_logistic_model([0.1, np.nan], [0.0, 1.0])
        with pytest.raises(ValueError, match="positive and negative"):
            fit_frequency_logistic_model([0.1, 0.2], [1.0, 1.0])


class TestAnalyticalLogisticThreshold:
    """Tests for analytical_logistic_threshold."""

    def test_exact_inversion_of_sigmoid(self):
        weight = 10.0
        for fdr in (0.15, 0.1, 0.01):
            t = analytical_logistic_threshold(weight, fdr)
            expected = math.log((1.0 - fdr) / fdr) / weight
            assert t == pytest.approx(expected)
            pred = 1.0 / (1.0 + math.exp(-weight * t))
            assert pred == pytest.approx(1.0 - fdr)

    def test_stricter_fdr_requires_higher_threshold(self):
        fdrs = [0.15, 0.1, 0.05, 0.01]
        thresholds = [analytical_logistic_threshold(12.0, f) for f in fdrs]
        assert thresholds == sorted(thresholds)

    def test_unreachable_or_nonpositive_returns_inf(self):
        assert math.isinf(analytical_logistic_threshold(0.0, 0.1))
        assert math.isinf(analytical_logistic_threshold(-2.0, 0.1))
        assert math.isinf(analytical_logistic_threshold(10.0, 0.0))
        # Exact threshold for weight=2, fdr=0.001 is ln(999)/2 ~= 3.45 > max_frequency=1.0
        assert math.isinf(analytical_logistic_threshold(2.0, 0.001, max_frequency=1.0))

    def test_min_value_floor_respected(self):
        assert analytical_logistic_threshold(50.0, 0.2, min_value=0.25) == 0.25
        assert analytical_logistic_threshold(10.0, 1.0, min_value=0.05) == 0.05
