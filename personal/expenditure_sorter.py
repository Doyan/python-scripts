"""
Created on Sat Mar 06 21:52 2021

Sorter for expenditures, to help me get a grip on my personal economy...
With Gui in streamlit.


@author: Gabriel Gustafsson
"""
#%%
import numpy as np
import pandas as pd

from ruyaml import YAML

from pathlib import Path

path = Path(__file__)

config_filename = "tracker_config.yml"
config_path = path.parent / config_filename

yaml = YAML(typ="safe")

with open(config_path, "r", encoding="UTF-8") as config_file:
    config = yaml.load(config_file)

# ==============================================================================


def read_file(path_to_file):
    df = pd.read_csv(path_to_file, delimiter=";")
    df["Bokföringsdag"] = pd.to_datetime(df["Bokföringsdag"], errors="coerce")
    df = df.dropna(subset=["Bokföringsdag"])
    df = df.sort_values(by="Bokföringsdag", ascending=True)
    return df


def read_records(records_path: Path) -> pd.DataFrame:

    records = [read_file(file) for file in records_path.iterdir()]
    if not all([all(r.keys() == records[0].keys()) for r in records]):
        raise KeyError("Supplied csvs dont have matching headers")

    records.reverse()
    df = records.pop()
    while len(records) > 0:
        df = pd.concat([df, records.pop()]).drop_duplicates().reset_index(drop=True)

    return df


def categorise_based_on_Rubrik(df, contains, category):
    df.category = np.where(df.Rubrik.str.contains(contains), category, df.category)
    return df


categorisation_funcs = {"Rubrik": categorise_based_on_Rubrik}


#%%
# main

df = read_records(Path(config["records_folder"]))

df = df.assign(category="uncategorised")

df["Belopp"] = df["Belopp"].str.replace(",", ".")
df["Belopp"] = pd.to_numeric(df["Belopp"])

income = df[df["Belopp"] > 0].copy()
expenses = df[df["Belopp"] < 0].copy()

#%%
for filter_type, income_filters in config["income_filters"].items():
    for category, entries in income_filters.items():
        for entry in entries:
            income = categorisation_funcs[filter_type](income, entry, category)

income[income.category == "uncategorised"]
#%%

for filter_type, expense_filters in config["expense_filters"].items():
    for category, entries in expense_filters.items():
        for entry in entries:
            expenses = categorisation_funcs[filter_type](expenses, entry, category)

expenses[expenses.category == "uncategorised"]
