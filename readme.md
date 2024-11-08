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
5. Make sure to setup the .env file with the correct paths to the data and models. A template is provided in the .env.template file.
6. Run the jupyter notebook to generate the retrosynthesis predictions.

## ToDos
- [ ] Add Flask API for retrosynthesis predictions
- [ ] Add Dockerfile for easy deployment
- [ ] Add tests for the retrosynthesis predictions
- [ ] Add CI/CD pipeline for the project
- [ ] Add documentation for the project
- [ ] Add clean logging for the project