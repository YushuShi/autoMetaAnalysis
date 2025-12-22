from meta_analysis import get_analysis_data

print("Running test...")
result = get_analysis_data("Lung Cancer", "Smoking")
print("Headline:", result.get('headline'))
if result.get('headline') is None:
    print("Error in headline generation.")
else:
    print("Success!")
