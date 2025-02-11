# Repo for SI Retrosynthesis Challenge 2024 - DFS

## Setup Instructions
1. Clone the repository
2. Create a conda environment using the provided environment.yml file using the below command:
```bash
conda env create -f environment.yml
```
3. Activate the environment using the below command:
```bash
conda activate dfs_si_challenge
```
4. download the Aizynthfinder models using the below command:
```bash
mkdir aizynthfinder/models
download_public_data aizynthfinder/models/
```
5. Download the reaction classification model from the [link](https://drive.google.com/file/d/1sQVTGDkKdhf1KiTe4xqZuIMh5WR0ZKQZ/view?usp=sharing) and place it in the `reaction_prediction` folder.
6. Make sure to setup the .env file with the correct paths to the data and models. A template is provided in the .env.template file.
7. Run the jupyter notebook to generate the retrosynthesis predictions.

## ToDos
Refer to github issues for the same.