import numpy as np

import genulens


def test_python_binding():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=10, seed=1234)
    result = genulens.simulate(cfg)
    arr = result.to_numpy()
    assert isinstance(arr, np.ndarray)
    assert arr.shape[0] == 10
    assert "tE" in result.columns


def test_python_callable_likelihood():
    cfg = genulens.Config(l=1.0, b=-3.9, n_simu=10, seed=1234)

    def my_like(event):
        return 1.0 if event.tE > 10 else 0.0

    result = genulens.simulate(cfg, likelihood=my_like)
    assert result.to_numpy().shape[0] == 10

