import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg') # Set non-interactive backend for Flask
import matplotlib.pyplot as plt
from Bio import Entrez
from dotenv import load_dotenv
import statsmodels.api as sm
from statsmodels.stats.meta_analysis import CombineResults, combine_effects
import forestplot

# ... (rest of imports)

def get_analysis_data(disease, exposure, outcome="Incidence", exclude_meta=False):
    """
    Main entry point for web app. Returns a dict with results.
    """
    print(f"Analyzing: {disease} vs {exposure} (Outcome: {outcome}, Exclude Meta/Reviews: {exclude_meta})")
    ids = search_pubmed(disease, exposure, outcome=outcome, exclude_meta=exclude_meta, max_results=100)
    articles = fetch_details(ids)
    
    df = extract_data(articles)
    
    if df.empty:
        return {"error": "No suitable data found extraction effect sizes."}

    df['SE'] = df.apply(calculate_se, axis=1)
    df_clean = df.dropna(subset=['Effect Size', 'SE'])
    
    if df_clean.empty:
        return {"error": "Effect sizes found but no Confidence Intervals to calculate SE."}

    return perform_meta_analysis(df_clean, disease, exposure)

def perform_meta_analysis(df_clean, disease, exposure):
    """
    Performs random-effects meta-analysis on the provided DataFrame.
    """
    # Sort by Effect Size for consistent display (Table and Plot)
    df_clean = df_clean.sort_values(by='Effect Size', ascending=True)

    # Meta-Analysis
    # Log transformation logic
    df_clean['log_ES'] = df_clean.apply(lambda x: np.log(x['Effect Size']) if x['Effect Type'] in ['OR', 'RR', 'HR', 'ODDS RATIO', 'RISK RATIO'] and x['Effect Size'] > 0 else x['Effect Size'], axis=1)
    
    def calc_log_se(row):
        if row['Effect Type'] in ['OR', 'RR', 'HR', 'ODDS RATIO', 'RISK RATIO'] and row['Lower CI'] > 0 and row['Upper CI'] > 0:
             return (np.log(row['Upper CI']) - np.log(row['Lower CI'])) / 3.92
        return row['SE']

    df_clean['log_SE'] = df_clean.apply(calc_log_se, axis=1)
    df_clean['var'] = df_clean['log_SE'] ** 2
    
    # Set index to Study for better summary labels
    # Set index to Study for better summary labels
    analysis_df = df_clean.set_index('Study')
    
    try:
        if len(analysis_df) < 2:
            # Handling Single Study: Use original raw values to match exactly
            row = analysis_df.iloc[0]
            pooled_es = row['Effect Size']
            pooled_lower = row['Lower CI']
            pooled_upper = row['Upper CI']
            
            # Create a dummy summary_df just to satisfy variable existence for summary_html (even if unused)
            summary_df = pd.DataFrame({
                'Effect': [pooled_es],
                '95% CI lower': [pooled_lower],
                '95% CI upper': [pooled_upper]
            }, index=['Pooled Result (Single Study)'])
            summary = summary_df.to_html(classes='table table-striped', header=True)

            # Interpretation logic
            is_significant = (pooled_lower > 1) or (pooled_upper < 1) 
            # Note: The above significance check assumes Ratio (null=1). 
            # If not ratio (e.g. 0), it should be diff > 0. 
            # Use confidence interval crossing null hypothesis check based on CI signs?
            # Actually, standard way: if lower and upper are on same side of Null.
            # Ratios are always > 0.
            
            # Let's improve significance check based on type
            if row['Effect Type'] in ['OR', 'RR', 'HR', 'ODDS RATIO', 'RISK RATIO']:
                 is_significant = (pooled_lower > 1) or (pooled_upper < 1)
                 log_eff = np.log(pooled_es) # For direction check
            else:
                 # Linear scale, null is 0
                 is_significant = (pooled_lower > 0 and pooled_upper > 0) or (pooled_lower < 0 and pooled_upper < 0)
                 log_eff = pooled_es # Just for direction
            
            interpretation = "Statistically Significant" if is_significant else "Not Statistically Significant"
            if is_significant:
                 direction = "Increased Risk/Odds" if log_eff > 0 else "Decreased Risk/Odds"
                 # Correction: if ratio, log_eff > 0 means ES > 1. Correct.
                 # If linear, log_eff > 0 means ES > 0. Correct.
                 interpretation += f" ({direction})"
            
            headline = {
                 "pooled_es": float(round(pooled_es, 2)),
                 "ci_low": float(round(pooled_lower, 2)),
                 "ci_upp": float(round(pooled_upper, 2)),
                 "interpretation": interpretation
            }
            
        else:
            res = combine_effects(analysis_df['log_ES'], analysis_df['var'], method_re='dl')
            summary_df = res.summary_frame()
            
            # Explicitly update the index
            n_studies = len(analysis_df)
            new_index = list(analysis_df.index) + list(summary_df.index[n_studies:])
            summary_df.index = new_index
            
            # Cleanup table columns
            cols_to_drop = ['w_fe', 'w_re']
            summary_df = summary_df.drop(columns=cols_to_drop, errors='ignore')
            
            # Rename columns
            summary_df = summary_df.rename(columns={
                'eff': 'Effect',
                'sd_eff': 'SD Effect',
                'ci_low': '95% CI lower',
                'ci_upp': '95% CI upper'
            })
            
            # Rename rows
            summary_df = summary_df.rename(index={
                'random effect wls': 'Random-effects meta-analysis (WLS)',
                'fixed effect wls': 'Fixed-effect meta-analysis (WLS)'
            })

            summary_df = summary_df.round(4)
            
            # Display DF
            display_df = summary_df.copy()
            rows_to_drop = ['Fixed-effect meta-analysis (WLS)', 'Random-effects meta-analysis (WLS)']
            rows_to_drop.extend(['fixed effect wls', 'random effect wls', 'fixed effect', 'random effect'])
            display_df = display_df.drop(index=[r for r in rows_to_drop if r in display_df.index], errors='ignore')
            
            summary = display_df.to_html(classes='table table-striped', header=True)
            
            # Extract Headline
            try:
                # Handle potential row naming variations
                if 'Random-effects meta-analysis (WLS)' in summary_df.index:
                    re_row = summary_df.loc['Random-effects meta-analysis (WLS)']
                elif 'random effect wls' in summary_df.index:
                    re_row = summary_df.loc['random effect wls']
                elif 'random effect' in summary_df.index:
                    re_row = summary_df.loc['random effect']
                else:
                    re_row = summary_df.iloc[-1]
            
                log_eff = re_row['Effect']
                log_ci_low = re_row['95% CI lower']
                log_ci_upp = re_row['95% CI upper']
                
                pooled_es = np.exp(log_eff)
                pooled_lower = np.exp(log_ci_low)
                pooled_upper = np.exp(log_ci_upp)
                
                is_significant = (log_ci_low > 0) or (log_ci_upp < 0)
                
                interpretation = "Statistically Significant" if is_significant else "Not Statistically Significant"
                
                if is_significant:
                    direction = "Increased Risk/Odds" if log_eff > 0 else "Decreased Risk/Odds"
                    interpretation += f" ({direction})"
                
                headline = {
                    "pooled_es": float(round(pooled_es, 2)),
                    "ci_low": float(round(pooled_lower, 2)),
                    "ci_upp": float(round(pooled_upper, 2)),
                    "interpretation": interpretation
                }

            except Exception as e:
                print(f"Error parsing summary stats: {e}")
                headline = None
        
        # Plot
        plt.figure(figsize=(10, 6))

        fp_df = df_clean.copy()
        fp_df = fp_df.rename(columns={'Study': 'group', 'log_ES': 'est'})
        fp_df['lb'] = fp_df['est'] - 1.96 * fp_df['log_SE']
        fp_df['ub'] = fp_df['est'] + 1.96 * fp_df['log_SE']
        fp_df['label'] = fp_df['group']
        
        forestplot.forestplot(
            fp_df,
            estimate="est",
            ll="lb",
            hl="ub",
            varlabel="label",
            xlabel="Log Effect Size (95% CI)",
            title=f"Forest Plot: {disease} vs {exposure}"
        )
        
        plot_path = os.path.join("static", "forest_plot.png")
        if not os.path.exists("static"):
            os.makedirs("static")
        plt.savefig(plot_path, bbox_inches='tight')
        plt.close() 
        
        # Convert df to records
        studies_data = df_clean[['Study', 'Effect Size', 'Lower CI', 'Upper CI', 'Population', 'Reference', 'Authors', 'Journal', 'Year', 'Link', 'Effect Type', 'SE', 'Sample Size']].to_dict(orient='records')
        
        return {
            "success": True,
            "studies": studies_data,
            "summary_html": summary,
            "headline": headline,
            "plot_url": "static/forest_plot.png?t=" + str(np.random.randint(0,10000)) 
        }
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        return {"error": f"Meta-analysis failed: {str(e)}"}

# Keep main for CLI usage but renamed/refactored if needed, or just let the new function handle it.
# We will modify the existing main to use this new function if we wanted to keep CLI, 
# but for now I'm just injecting the function to be used by Flask.


# Load environment variables
load_dotenv('mykey.env')

# Setup Entrez
Entrez.email = os.getenv('PUBMED_EMAIL', 'your_email@example.com')

def search_pubmed(disease, exposure, outcome="Incidence", exclude_meta=False, max_results=20):
    """
    Search PubMed for articles related to the disease and exposure.
    """
    # Define outcome terms
    if outcome == "Survival":
        outcome_terms = '(survival OR mortality OR prognosis OR "overall survival" OR "OS" OR "hazard ratio" OR "death")'
    elif outcome == "Disease-Free Survival":
        outcome_terms = '("disease-free survival" OR "DFS" OR "recurrence-free survival" OR "RFS" OR "relapse-free survival")'
    else:
        # Default to Incidence
        outcome_terms = '(incidence OR risk OR development OR "associated with" OR "odds ratio")'

    if exclude_meta:
        # Exclude Meta-Analysis and Reviews, prioritize primary studies
        query = f"{disease} AND {exposure} AND {outcome_terms} AND (Clinical Trial[ptyp] OR Randomized Controlled Trial[ptyp] OR Journal Article[ptyp]) NOT (Meta-Analysis[ptyp] OR Review[ptyp] OR Systematic Review[ptyp])"
    else:
        # Original inclusive query
        query = f"{disease} AND {exposure} AND {outcome_terms} AND (Meta-Analysis[ptyp] OR Review[ptyp] OR Clinical Trial[ptyp] OR Journal Article[ptyp])"
        
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
    # Expanded regex to capture more formats like "OR=1.2", "OR 1.2", "relative risk of 1.2"
    es_pattern = re.compile(r'\b(OR|RR|HR|Odds Ratio|Risk Ratio|Hazard Ratio)\b.*?[:=]?\s*(\d+\.\d+)', re.IGNORECASE | re.DOTALL)
    # CI Pattern: looks for (95% CI: 1.1-2.2) or (1.1, 2.2) or similar variants
    ci_pattern = re.compile(r'\(\s*(?:95\s*%\s*C\.?I\.?)?\s*[:=]?\s*(\d+\.\d+)\s*[-–,to]\s*(\d+\.\d+)\s*\)', re.IGNORECASE)
    
    last_abstract_debug = ""

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
            last_abstract_debug = abstract

            # Extract Effect Size (Heuristics)
            es_match = es_pattern.search(abstract)
            effect_size = None
            es_type = None
            lower_ci = None
            upper_ci = None
            
            if es_match:
                raw_type = es_match.group(1).upper()
                if "ODDS" in raw_type: es_type = "OR"
                elif "RISK" in raw_type or "RR" in raw_type: es_type = "RR"
                elif "HAZARD" in raw_type or "HR" in raw_type: es_type = "HR"
                else: es_type = raw_type
                try:
                    effect_size = float(es_match.group(2))
                except ValueError:
                    continue
                
                # Look for CI nearby - Extended window
                start_pos = es_match.end()
                snippet = abstract[start_pos:start_pos+150] # Extended to 150
                
                # Improved CI regex: allows for square brackets, various separators, optional '95% CI' prefix nearby
                # (?: ... ) is non-capturing group
                # We interpret CI as two numbers separated by -, to, or ,
                # We assume they appear within parens or brackets OR just after "95% CI"
                # This is tricky regex.
                # Let's try matching "number [sep] number" that are close to "CI" or in brackets
                
                # Regex for "1.23-4.56" or "1.23, 4.56" or "1.23 to 4.56"
                # enclosed in parens/brackets OR preceded by CI
                
                # Case 1: (...) or [...] containing two numbers
                ci_pattern_1 = re.compile(r'[(\[]\s*(?:95\s*%\s*C\.?I\.?[:\s]*)?(\d+\.\d+)\s*[-–,;to]+\s*(\d+\.\d+)\s*[)\]]', re.IGNORECASE)
                
                # Case 2: "95% CI 1.23-4.56" (no parens around numbers)
                ci_pattern_2 = re.compile(r'95\s*%\s*C\.?I\.?[:\s]*(\d+\.\d+)\s*[-–,;to]+\s*(\d+\.\d+)', re.IGNORECASE)
                
                ci_match = ci_pattern_1.search(snippet)
                if not ci_match:
                    ci_match = ci_pattern_2.search(snippet)
                    
                if not ci_match:
                     # Fallback: search wide in snippet for just two numbers if "CI" is mentioned
                     if "CI" in snippet or "confidence interval" in snippet.lower():
                         # Just find two floats
                         nums = re.findall(r'(\d+\.\d+)', snippet)
                         if len(nums) >= 2:
                             # Assume first two are the CI if they bracket the ES? 
                             # Or just take them.
                             try:
                                 v1, v2 = float(nums[0]), float(nums[1])
                                 ci_match = type('Match', (object,), {'group': lambda s, i: v1 if i==1 else v2})()
                             except:
                                 pass

                if ci_match:
                     try:
                        lower_ci = float(ci_match.group(1))
                        upper_ci = float(ci_match.group(2))
                     except ValueError:
                        pass
                else:
                    # Debug print for missed CI
                    if len(data) < 5:
                       print(f"DEBUG: Found ES {effect_size} in '{title[:20]}...' but NO CI in snippet: '{snippet}'")

            # Journal and Year
            journal_info = article_data.get('Journal', {})
            journal_title = journal_info.get('Title', 'Unknown Journal')
            # Year can be tricky in Medline (sometimes in MedlineDate)
            pub_date = journal_info.get('JournalIssue', {}).get('PubDate', {})
            year = pub_date.get('Year', '')
            if not year:
                # Try MedlineDate
                medline_date = pub_date.get('MedlineDate', '')
                year_match = re.search(r'\d{4}', medline_date)
                if year_match:
                    year = year_match.group(0)
                else:
                    year = "Unknown"

            if effect_size:
                # Basic validation:
                if lower_ci and upper_ci:
                    if lower_ci > upper_ci:
                        lower_ci, upper_ci = upper_ci, lower_ci
                    
                    # Validate that Effect Size is within the CI (allowing small epsilon for rounding)
                    # If ES is outside CI, it's likely a parsing error of unrelated numbers
                    if not (lower_ci <= effect_size <= upper_ci):
                         # print(f"DEBUG: Discarding {title[:30]}... ES {effect_size} not in CI {lower_ci}-{upper_ci}")
                         continue
                
                # Format Study with Year
                short_author = f"{authors.split(',')[0]} et al." if ',' in authors else authors
                study_label = f"{short_author} ({year})"
                
                # Construct PubMed Link
                pmid = medline.get('PMID', '')
                pmid_link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "#"
                
                # Attempt to extract Sample Size
                sample_size = "N/A"
                # Patterns to look for sample size:
                # 1. "n = 123" or "N = 123"
                # 2. "123 patients" or "123 participants" or "123 cases"
                # 3. "total of 123"
                
                # We prioritize specific patterns
                n_match = re.search(r'\b[nN]\s*=\s*(\d+(?:,\d{3})*)', abstract)
                if n_match:
                    sample_size = n_match.group(1)
                else:
                    # Look for larger numbers followed by participants/patients
                    # Avoid years like 2020
                    part_match = re.search(r'\b(\d+(?:,\d{3})*)\s+(participants|patients|subjects|cases|women|men|individuals)', abstract, re.IGNORECASE)
                    if part_match:
                        # Simple check to avoid years (e.g. 1990-2010 participants... wait, 2010 participants is valid)
                        # Let's assume if it is > 2025 or < 1900 it's likely a number, OR if formatted with comma
                        val_str = part_match.group(1).replace(',', '')
                        val = int(val_str)
                        if val > 10 and (val < 1900 or val > 2030 or ',' in part_match.group(1)):
                             sample_size = part_match.group(1)
                
                row = {
                    "Study": study_label,
                    "Effect Size": effect_size,
                    "Effect Type": es_type,
                    "Sample Size": sample_size,
                    "Lower CI": lower_ci,
                    "Upper CI": upper_ci,
                    "Population": "General",
                    "Authors": authors,
                    "Reference": title,
                    "Journal": journal_title,
                    "Year": year,
                    "Link": pmid_link
                }
                
                # Attempt Population extraction
                if "children" in abstract.lower(): row["Population"] = "Children"
                elif "adults" in abstract.lower(): row["Population"] = "Adults"
                elif "patients" in abstract.lower(): row["Population"] = "Patients"

                data.append(row)

            
        except Exception as e:
            continue

    # Fallback: if no data found
    if not data:
        print("DEBUG: No data found. Showing snippet of last abstract processed to help debug:")
        if last_abstract_debug:
            print(last_abstract_debug[:200])
            
    return pd.DataFrame(data)

def calculate_se(row):
    """Calculate Standard Error from CI if available."""
    if pd.notnull(row['Lower CI']) and pd.notnull(row['Upper CI']):
        # Assuming 95% CI and Normal dist, width is 3.92 * SE
        return (row['Upper CI'] - row['Lower CI']) / 3.92
    return None

def main():
    print("--- Meta-Analysis Tool ---")
    disease = input("Enter Disease (e.g., 'Breast Cancer'): ") or "Breast Cancer"
    exposure = input("Enter Exposure (e.g., 'Coffee'): ") or "Coffee"
    
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
    df_clean['log_ES'] = df_clean.apply(lambda x: np.log(x['Effect Size']) if x['Effect Type'] in ['OR', 'RR', 'HR', 'ODDS RATIO', 'RISK RATIO'] and x['Effect Size'] > 0 else x['Effect Size'], axis=1)
    # SE also needs to be on log scale? Yes, SE(logOR) ~= (log(Upper) - log(Lower)) / 3.92
    # So let's just recalculate log SE
    
    def calc_log_se(row):
        if row['Effect Type'] in ['OR', 'RR', 'HR', 'ODDS RATIO', 'RISK RATIO'] and row['Lower CI'] > 0 and row['Upper CI'] > 0:
             return (np.log(row['Upper CI']) - np.log(row['Lower CI'])) / 3.92
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
    # Let's try `from statsmodels.stats.meta_analysis import combine_effects`
    
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
        
        res = combine_effects(df_clean['log_ES'], df_clean['var'], method_re='dl')
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
        
        forestplot.forestplot(
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
