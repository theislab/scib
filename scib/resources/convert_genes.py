import requests
from tqdm import tqdm


def get_gene_name_from_ensembl(gene_ids: list, species: str):
    base_url = "https://rest.ensembl.org"
    gene_names = []

    for gene_id in tqdm(gene_ids):
        response = requests.get(
            f"{base_url}/lookup/id/{gene_id}?expand=1;species={species}",
            headers={"Content-Type": "application/json"},
        )

        if response.status_code == 200:
            data = response.json()
            gene_name = data.get("display_name", gene_id)
            gene_names.append(gene_name)
        else:
            print(f"Error: {response.status_code}, skipping gene {gene_id}...")
            gene_names.append(gene_id)

    return gene_names


def get_gene_id_from_ensembl(gene_names: list, species: str):
    base_url = "https://rest.ensembl.org"
    gene_ids = []

    for gene_name in tqdm(gene_names):
        response = requests.get(
            f"{base_url}/xrefs/symbol/{species}/{gene_name}?expand=1",
            headers={"Content-Type": "application/json"},
        )

        if response.status_code == 200:
            data = response.json()
            if data:
                gene_id = data[0].get(
                    "id", gene_name
                )  # Get the first result's Ensembl ID
                gene_ids.append(gene_id)
            else:
                gene_ids.append("Not Found")
        else:
            print(f"Error: {response.status_code}, skipping gene {gene_name}...")
            gene_ids.append(gene_name)

    return gene_ids


if __name__ == "__main__":
    from pathlib import Path

    import pandas as pd

    root = Path(__file__).parent

    cc_files = {
        "mus_musculus": "https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt",
        "homo_sapiens": "https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt",
        "caenorhabditis_elegans": "https://raw.githubusercontent.com/hbc/tinyatlas/refs/heads/master/cell_cycle/Caenorhabditis_elegans.csv",
        "danio_rerio": "https://raw.githubusercontent.com/hbc/tinyatlas/refs/heads/master/cell_cycle/Danio_rerio.csv",
    }

    # Tirosh mouse and human
    # processed according to https://github.com/scverse/scanpy_usage/blob/master/180209_cell_cycle/cell_cycle.ipynb
    for organism in ["mus_musculus", "homo_sapiens"]:
        print(f"Organism: {organism}")

        # read file
        gene_names = pd.read_csv(cc_files[organism], header=None)[0]

        if organism == "mus_musculus":
            gene_names = gene_names.str.capitalize()

        # convert gene names
        gene_ids = get_gene_id_from_ensembl(gene_names, species=organism)

        # create gene map
        gene_map = pd.DataFrame(dict(gene_name=gene_names, gene_id=gene_ids))

        # set cell cycle phase
        gene_map.loc[:43, "phase"] = "S"
        gene_map.loc[43:, "phase"] = "G2/M"

        # write to file
        print(gene_map)
        gene_map.to_csv(
            root / f"cell_cycle_genes_{organism}.tsv", sep="\t", index=False
        )

    # Tinyatlas gene sets
    # https://github.com/hbc/tinyatlas/tree/master/cell_cycle
    for organism in ["caenorhabditis_elegans", "danio_rerio"]:
        print(f"Organism: {organism}")

        # read file
        gene_map = pd.read_csv(cc_files[organism])
        gene_map["gene_id"] = gene_map["geneID"]
        del gene_map["geneID"]

        # get gene names
        gene_map["gene_name"] = get_gene_name_from_ensembl(
            gene_map["gene_id"], species=organism
        )

        # write to file
        print(gene_map)
        gene_map.to_csv(
            root / f"cell_cycle_genes_{organism}.tsv", sep="\t", index=False
        )
