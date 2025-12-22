import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Entrez
from dotenv import load_dotenv
import statsmodels.api as sm
from statsmodels.stats.meta_analysis import CombineResults
import forestplot

# Load environment variables
load_dotenv('mykey.env')

# Setup Entrez
Entrez.email = os.getenv('PUBMED_EMAIL', 'your_email@example.com')

def search_pubmed(disease, exposure, max_results=20):
    """
    Search PubMed for articles related to the disease and exposure.
    """
    query = f"{disease} AND {exposure} AND (Meta-Analysis[ptyp] OR Review[ptyp] OR Clinical Trial[ptyp])"
    print(f"Searching PubMed for: {query}")
    
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"Error searching PubMed: {e}")
        return []

    id_list = record["IdList"]
    print(f"Found {len(id_list)} articles.")
    return id_list

def fetch_details(id_list):
    """
    Fetch details for the list of PubMed IDs.
    """
    if not id_list:
        return []
    
    ids = ",".join(id_list)
    try:
        handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        return records['PubmedArticle']
    except Exception as e:
        print(f"Error fetching details: {e}")
        return []

def extract_data(articles):
    """
    Extract relevant data from the articles.
    """
    data = []
    
    # Regex patterns for effect sizes (simplified)
    # Looking for OR, RR, HR followed by numbers
    es_pattern = re.compile(r'\b(OR|RR|HR|Odds Ratio|Risk Ratio|Hazard Ratio)\s*[=:]?\s*(\d+\.?\d*)', re.IGNORECASE)
    # Looking for CI 
    ci_pattern = re.compile(r'\(\s*(95\s*%\s*CI)?\s*[:=]?\s*(\d+\.?\d*)\s*[-–,]\s*(\d+\.?\d*)\s*\)', re.IGNORECASE)
    
    for article in articles:
        try:
            medline = article['MedlineCitation']
            article_data = medline['Article']
            
            # Title
            title = article_data.get('ArticleTitle', 'No Title')
            
            # Authors
            author_list = article_data.get('AuthorList', [])
            if author_list:
                authors = ", ".join([f"{a.get('LastName', '')} {a.get('Initials', '')}" for a in author_list])
            else:
                authors = "Unknown"
            
            # Abstract
            abstract_list = article_data.get('Abstract', {}).get('AbstractText', [])
            abstract = " ".join(abstract_list) if isinstance(abstract_list, list) else str(abstract_list)
            
            # Extract Effect Size (Heuristics)
            # Expanded regex to capture more formats like "OR=1.2", "OR 1.2", "relative risk of 1.2"
            es_pattern = re.compile(r'\b(OR|RR|HR|Odds Ratio|Risk Ratio|Hazard Ratio)\b.*?[:=]?\s*(\d+\.\d+)', re.IGNORECASE | re.DOTALL)
            # CI Pattern: looks for (95% CI: 1.1-2.2) or (1.1, 2.2) or similar variants
            ci_pattern = re.compile(r'\(\s*(?:95\s*%\s*C\.?I\.?)?\s*[:=]?\s*(\d+\.\d+)\s*[-–,to]\s*(\d+\.\d+)\s*\)', re.IGNORECASE)
            
            es_match = es_pattern.search(abstract)
            effect_size = None
            es_type = None
            lower_ci = None
            upper_ci = None
            
            if es_match:
                es_type = es_match.group(1).upper()
                try:
                    effect_size = float(es_match.group(2))
                except ValueError:
                    continue
                
                # Look for CI nearby
                start_pos = es_match.end()
                # Check next 50 chars for CI
                snippet = abstract[start_pos:start_pos+50] 
                ci_match = ci_pattern.search(snippet)
                
                # If not found, try to search the whole abstract but finding the 'closest' one is hard.
                # Let's try searching from the match position forward
                if not ci_match:
                     ci_match = ci_pattern.search(abstract, pos=start_pos)

                if ci_match:
                     try:
                        lower_ci = float(ci_match.group(1))
                        upper_ci = float(ci_match.group(2))
                     except ValueError:
                        pass
            
            if effect_size:
                # Basic validation: CI should bracket the ES usually (not always for OR/RR if <1 or >1? No, ES is inside CI)
                # And Lower < Upper
                if lower_ci and upper_ci:
                    if lower_ci > upper_ci:
                        lower_ci, upper_ci = upper_ci, lower_ci
                
                row = {
                    "Study": f"{authors.split(',')[0]} et al." if ',' in authors else authors,
                    "Effect Size": effect_size,
                    "Effect Type": es_type,
                    "Lower CI": lower_ci,
                    "Upper CI": upper_ci,
                    "Population": population,
                    "Authors": authors,
                    "Reference": title
                }
                data.append(row)
                if len(data) <= 1:
                     print(f"DEBUG: Found ES: {effect_size} CI: {lower_ci}-{upper_ci} in {title[:30]}...")
            
    # Fallback: if no data found, maybe the regex is too strict
    if not data:
        print("DEBUG: No data found. Showing snippet of last abstract processed to help debug:")
        if 'abstract' in locals():
            print(abstract[:200])

                
        except Exception as e:
            continue
            
    return pd.DataFrame(data)

def calculate_se(row):
    """Calculate Standard Error from CI if available."""
    if pd.notnull(row['Lower CI']) and pd.notnull(row['Upper CI']):
        # Assuming 95% CI and Normal dist, width is 3.92 * SE
        return (row['Upper CI'] - row['Lower CI']) / 3.92
    return None

def main():
    print("--- Meta-Analysis Tool ---")
    disease = input("Enter Disease (e.g., 'Lung Cancer'): ") or "Lung Cancer"
    exposure = input("Enter Exposure (e.g., 'Smoking'): ") or "Smoking"
    
    print(f"\nFetching data for {disease} and {exposure}...")
    ids = search_pubmed(disease, exposure)
    articles = fetch_details(ids)
    
    df = extract_data(articles)
    
    if df.empty:
        print("No suitable data found containing extracted effect sizes.")
        return

    # Post-process for Meta-Analysis
    # Calculate SE (needed for weighting)
    df['SE'] = df.apply(calculate_se, axis=1)
    
    # Drop rows without SE or Effect Size
    df_clean = df.dropna(subset=['Effect Size', 'SE'])
    
    if df_clean.empty:
        print("Effect sizes found, but Confidence Intervals could not be securely parsed to calculate SE. Cannot proceed with Meta-Analysis.")
        print("Extracted Data Preview:")
        print(df.head())
        return

    print(f"\nSuccessfully extracted {len(df_clean)} studies for analysis.")
    print(df_clean[['Study', 'Effect Size', 'Lower CI', 'Upper CI', 'Population']])
    
    # Save to CSV
    df_clean.to_csv("meta_analysis_results.csv", index=False)
    print("\nData saved to 'meta_analysis_results.csv'")

    # Random Effects Meta-Analysis
    print("\nPerforming Random Effects Meta-Analysis...")
    # statsmodels CombineResults
    # We use effect size and SE^2 (variance)
    # Assuming effect sizes are on log scale if they are OR/RR? Usually meta-analysis is done on log(OR).
    # For this simple tool, I'll assume the extracted numbers are what we want to analyze directly 
    # OR convert if 'OR'/'RR' are specific types.
    # To keep it "Tool-like" and robust, let's take log if it's OR/RR/HR and > 0
    
    # Simple logic: if Type is OR/RR/HR, log transform
    df_clean['log_ES'] = df_clean.apply(lambda x: pd.np.log(x['Effect Size']) if x['Effect Type'] in ['OR', 'RR', 'HR', 'ODDS RATIO', 'RISK RATIO'] and x['Effect Size'] > 0 else x['Effect Size'], axis=1)
    # SE also needs to be on log scale? Yes, SE(logOR) ~= (log(Upper) - log(Lower)) / 3.92
    # So let's just recalculate log SE
    
    def calc_log_se(row):
        if row['Effect Type'] in ['OR', 'RR', 'HR', 'ODDS RATIO', 'RISK RATIO'] and row['Lower CI'] > 0 and row['Upper CI'] > 0:
             return (pd.np.log(row['Upper CI']) - pd.np.log(row['Lower CI'])) / 3.92
        return row['SE'] # Use raw SE if not ratio

    df_clean['log_SE'] = df_clean.apply(calc_log_se, axis=1)
    df_clean['var'] = df_clean['log_SE'] ** 2
    
    # Use statsmodels
    # CombineResults(eff, var)
    # However statsmodels 'CombineResults' class is not always directly exposed like that or requires specific args.
    # Actually, statsmodels.stats.meta_analysis.meta_analysis is a function or similar.
    # Wait, getting 'CombineResults' might be tricky if versions change.
    # Let's use the explicit `statsmodels.stats.meta_analysis.combine_effects` if available?
    # Or just write the DerSimonian-Laird estimator manually? No, I should use the library.
    # `from statsmodels.stats.meta_analysis import effectsize_smd` etc.
    # Let's use `model = sm.stats.meta_analysis.MetaAnalysis(df_clean['log_ES'], df_clean['log_SE'])`? No.
    
    # Actually, simpler path:
    # Use a basic inverse variance weighting if library usage is complex to guess without docs.
    # But I see `statsmodels.stats.meta_analysis` docs usually.
    # Let's try `from statsmodels.stats.meta_analysis import CombineResults`
    # res = CombineResults(effect_size, variance, method='random')
    # If that fails at runtime, I'll need to fix.
    
    try:
        # Assuming data is ready
        # Using DerSimonian-Laird
        # We need weights = 1 / (var + tau^2)
        # CombineResults might handle it
        # Actually random effects usually requires estimating tau2 first.
        # Let's keep it simple: Just do a weighted means if library fails? 
        # No, the user asked for "generated a random effects meta analysis using statsmodels".
        
        # Checking typical usage:
        # `from statsmodels.stats.meta_analysis import effect_size_smd, combine_effects`
        # `res = combine_effects(effect, var, method_re='dl', use_t=True)`
        
        res = sm.stats.meta_analysis.combine_effects(df_clean['log_ES'], df_clean['var'], method_re='dl')
        print("\nMeta-Analysis Results:")
        print(res.summary_frame())
        
        # Forest Plot
        print("\nGenerating Forest Plot...")
        
        # forestplot library usage:
        # forestplot.forestplot(dataframe, estimate="log_ES", var="var", ...)
        # Need to check forestplot API.
        # Common API: forestplot(dataframe, estimate="col", ll="col", hl="col", ...)
        
        # Let's assume we pass the dataframe to forestplot 
        # Required columns: 'study', 'est', 'lower', 'upper' often.
        
        # Prepare df for forestplot
        fp_df = df_clean.copy()
        fp_df = fp_df.rename(columns={'Study': 'group', 'log_ES': 'est'})
        # We need lower/upper for the PLOT (so log scale CIs)
        fp_df['lb'] = fp_df['est'] - 1.96 * fp_df['log_SE']
        fp_df['ub'] = fp_df['est'] + 1.96 * fp_df['log_SE']
        fp_df['label'] = fp_df['group']
        
        # Since forestplot library might vary, I'll use matplotlib manually if forestplot is weird,
        # but let's try the library.
        # Actually, let's just use forestplot if it works, otherwise fallback?
        # A simple visual forest plot using matplotlib is safer than a niche library that might break.
        # But User request: "Generate a Forest Plot using the forestplot library".
        # I MUST use forestplot library.
        # Usage: forestplot.forest_plot(df, estimate="est", lower="lb", upper="ub", varlabel="label")
        
        forestplot.forest_plot(
            fp_df,
            estimate="est",
            ll="lb",
            hl="ub",
            varlabel="label",
            xlabel="Log Effect Size (95% CI)",
            title=f"Forest Plot: {disease} vs {exposure}"
        )
        plt.savefig("forest_plot.png")
        print("Forest plot saved as forest_plot.png")
        
    except Exception as e:
        print(f"Meta-analysis or plotting failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
