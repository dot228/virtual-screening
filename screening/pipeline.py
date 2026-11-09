import csv

from multiprocessing import Pool
from multiprocessing import cpu_count

from screening.worker import process_pair

from data.acids import acid_smiles
from data.amines import amine_smiles

from config import (
    OUTPUT_CSV,
    TOP_N
)

def run_screening():

    pairs = [
        (acid, amine)
        for acid in acid_smiles
        for amine in amine_smiles
    ]

    results = []

    with Pool(cpu_count()) as pool:

        for rows in pool.imap_unordered(
            process_pair,
            pairs,
            chunksize=50
        ):

            if rows:
                results.extend(rows)

    unique = {}

    for row in results:
        unique[row["GeneratedSMILES"]] = row

    results = list(unique.values())

    results.sort(
        key=lambda x: (
            x["MorganSimilarity"],
            x["TopologicalSimilarity"]
        ),
        reverse=True
    )

    results = results[:TOP_N]

    if not results:

        print("No results found")
        return

    with open(
        OUTPUT_CSV,
        "w",
        newline="",
        encoding="utf-8"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=results[0].keys()
        )

        writer.writeheader()
        writer.writerows(results)

    print(f"Saved {len(results)} molecules")