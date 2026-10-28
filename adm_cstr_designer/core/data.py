"""Library loading (kinetics and feedstock databases)."""
import json
import pathlib
import copy

DATA = pathlib.Path(__file__).resolve().parents[1] / 'data'


def load_kinetics() -> dict:
    return json.loads((DATA / 'kinetics_db.json').read_text())['sets']


def load_feedstocks() -> dict:
    return json.loads((DATA / 'feedstocks.json').read_text())['feedstocks']


def feedstock_meta() -> str:
    return json.loads((DATA / 'feedstocks.json').read_text())['_meta']['status']


def clone(d):
    return copy.deepcopy(d)
