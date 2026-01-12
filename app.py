from flask import Flask, render_template, request, jsonify, send_from_directory
import meta_analysis
import os
import pandas as pd

# Define absolute paths for templates and static files based on current directory
try:
    TEMPLATE_DIR = os.path.abspath('templates')
    STATIC_DIR = os.path.abspath('static')
    if not os.path.exists(TEMPLATE_DIR):
        os.makedirs(TEMPLATE_DIR)
    if not os.path.exists(STATIC_DIR):
        os.makedirs(STATIC_DIR)
        
    app = Flask(__name__, template_folder=TEMPLATE_DIR, static_folder=STATIC_DIR)
except Exception as e:
    # Fallback if abspath fails (rare)
    app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/analyze', methods=['POST'])
def analyze():
    data = request.json
    disease = data.get('disease', 'Mesothelioma')
    exposure = data.get('exposure', 'Smoking')
    outcome = data.get('outcome', 'Incidence')
    exclude_meta = data.get('exclude_meta', False)
    
    result = meta_analysis.get_analysis_data(disease, exposure, outcome=outcome, exclude_meta=exclude_meta)
    return jsonify(result)

@app.route('/reanalyze', methods=['POST'])
def reanalyze():
    data = request.json
    studies = data.get('studies', [])
    disease = data.get('disease', 'Custom Analysis')
    exposure = data.get('exposure', 'Custom Exposure')
    
    if not studies:
        return jsonify({"error": "No studies provided for analysis."})
    
    # Convert list of dicts to DataFrame
    df = pd.DataFrame(studies)
    
    # Ensure numeric columns are floats
    numeric_cols = ['Effect Size', 'Lower CI', 'Upper CI', 'SE']
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
            
    # Calculate SE if missing (though the frontend should pass it back ideally, 
    # but we can re-calculate if needed, provided we have CIs)
    # However, get_analysis_data calculated SE effectively. 
    # The frontend will send back the full study objects which include 'SE'.
    
    result = meta_analysis.perform_meta_analysis(df, disease, exposure)
    return jsonify(result)

@app.route('/static/<path:filename>')
def serve_static(filename):
    return send_from_directory('static', filename)

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port)
