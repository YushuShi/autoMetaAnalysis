# Meta-Analysis Tool

## Overview
This tool searches PubMed for scientific articles related to a specific disease and exposure, extracts effect sizes (OR/RR/HR) and confidence intervals, and performs a random effects meta-analysis. It also generates a forest plot.

## Setup
1.  Ensure you have Python installed.
2.  Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```
3.  (Optional) Add your PubMed API key/Email in `mykey.env`.

## Usage
Run the script:
```bash
python meta_analysis.py
```
You will be prompted to enter the **Disease** and **Exposure**.
Default values:
- Disease: `Lung Cancer`
- Exposure: `Smoking`

## Output
- `meta_analysis_results.csv`: Extracted data from PubMed articles.
- `forest_plot.png`: Forest plot of the meta-analysis.
- Console output: Summary of the random effects model.

## Note on Extraction
This tool uses heuristics (Regex) to extract effect sizes from abstracts. It is not perfect and may miss data or pick up incorrect numbers. Always verify the results against the original papers.
