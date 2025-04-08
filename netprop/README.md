# Network Propagation

This folder (`netprop`) contains Python scripts for performing network propagation analysis. It includes two main scripts:

* `netprop_ctrp.py`: This script applies network propagation to identify genes relevant to drugs from the Cancer Therapeutics Response Portal (CTRP) dataset.
* `netprop_prism.py`: This script applies network propagation to identify genes relevant to drugs from the PRISM Repurposing dataset.

## Running the Scripts

Both scripts follow a similar structure and depend on the availability of specific data files. Here's a general guide on how to run them:

**1. Prerequisites:**

* **Python 3** is required.
* The following Python libraries need to be installed. You can install them using pip:
    ```bash
    pip install pandas numpy networkx tqdm
    ```

**2. Data Setup:**

The scripts expect the following data files to be located in the `../data/` directory with the following subdirectory structure:

* `../data/string/`:
    * `9606.protein.info.v12.0.txt`: Protein information from the STRING database.
        * **Source URL:** [https://string-db.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz](https://string-db.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz) (You'll need to unzip this file)
    * `9606.protein.links.full.v12.0.txt`: Protein-protein interaction (PPI) data from the STRING database.
        * **Source URL:** [https://string-db.org/download/protein.links.full.v12.0/9606.protein.links.full.v12.0.txt.gz](https://string-db.org/download/protein.links.full.v12.0/9606.protein.links.full.v12.0.txt.gz) (You'll need to unzip this file)
* `../data/ctrp/`: (Required for `netprop_ctrp.py`)
    * `v20.meta.per_compound.txt`: Metadata for compounds in the CTRP dataset, containing drug-target information.
        * **Source URL:** You will need to obtain this file from the CTRP project. A direct download link might not be publicly available and may require registration or access through specific channels. Please refer to the CTRP documentation or website for data access instructions.
* `../data/prism/`: (Required for `netprop_prism.py`)
    * `secondary-screen-replicate-collapsed-treatment-info.csv`: Treatment information from the PRISM Repurposing dataset, containing drug-target information.
        * **Source URL:** You will need to obtain this file from the PRISM Repurposing project. Similar to the CTRP data, a direct public link might not be available. Please refer to the Broad Institute's website or PRISM documentation for data access details.

**3. Running the Scripts:**

Navigate to the `netprop` directory within your repository in your terminal.

* **To run `netprop_ctrp.py`:**
    ```bash
    python netprop_ctrp.py
    ```
    This script will:
    * Read the STRING PPI data and the CTRP drug-target information.
    * Construct a protein-protein interaction network.
    * Perform network propagation for each drug based on its known targets.
    * Identify the top 20 propagated genes for each drug.
    * Save the results to `../data/wrangled/selected_genes_ctrp.csv`.

* **To run `netprop_prism.py`:**
    ```bash
    python netprop_prism.py
    ```
    This script will:
    * Read the STRING PPI data and the PRISM drug-target information.
    * Construct a protein-protein interaction network.
    * Perform network propagation for each drug based on its known targets.
    * Identify the top 20 propagated genes for each drug.
    * Save the results to `../data/wrangled/selected_genes_prism.csv`.

**4. Output:**

Both scripts will generate a CSV file in the `../data/wrangled/` directory:

* `selected_genes_ctrp.csv`: Contains two columns, `genes` and `drugs`, listing the top 20 propagated genes for each drug in the CTRP dataset.
* `selected_genes_prism.csv`: Contains two columns, `genes` and `drugs`, listing the top 20 propagated genes for each drug in the PRISM dataset.

**Data Directories:**

It is crucial to maintain the specified directory structure (`../data/string/`, `../data/ctrp/` and `../data/prism/`) and place the corresponding data files in their respective locations for the scripts to run correctly.