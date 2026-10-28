"""Headless UI smoke test: inputs -> run -> all tabs render without exceptions."""
import pathlib
from streamlit.testing.v1 import AppTest

APP = str(pathlib.Path(__file__).resolve().parents[1] / 'app.py')


def test_app_runs_and_designs():
    at = AppTest.from_file(APP, default_timeout=120).run()
    assert not at.exception
    next(b for b in at.sidebar.button if b.label == 'Run design').click()
    at.run()
    assert not at.exception, at.exception
    labels = [m.label for m in at.metric]
    assert 'Methane (Nm³/d)' in labels and 'pH' in labels
    assert 'Impeller Ø (m)' in labels and 'Coil area (m²)' in labels
    assert any(b.label == 'Run scenarios' for b in at.button)
    next(b for b in at.button if b.label == 'Run scenarios').click()
    at.run()
    assert not at.exception, at.exception
    next(b for b in at.button if b.label == 'Prepare design report').click()
    at.run()
    assert not at.exception, at.exception
