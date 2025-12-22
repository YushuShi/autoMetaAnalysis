from flask import Flask, render_template, request, jsonify, send_from_directory
import meta_analysis
import os

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

@app.route('/analyze', methods=['POST'])
def analyze():
    data = request.json
    disease = data.get('disease', 'Lung Cancer')
    exposure = data.get('exposure', 'Smoking')
    
    result = meta_analysis.get_analysis_data(disease, exposure)
    return jsonify(result)

@app.route('/static/<path:filename>')
def serve_static(filename):
    return send_from_directory('static', filename)

if __name__ == '__main__':
    app.run(debug=True, port=5000)
