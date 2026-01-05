let currentStudies = [];

// Exposure Dropdown Logic
document.addEventListener('DOMContentLoaded', async () => {
    const exposureInput = document.getElementById('exposure');
    const exposureOptions = document.getElementById('exposure-options');
    let exposures = [];

    // Fetch exposures
    try {
        const res = await fetch('/static/exposures.json');
        if (res.ok) {
            exposures = await res.json();
            console.log(`Loaded ${exposures.length} exposures`);
        } else {
            console.error("Failed to load exposures.json");
        }
    } catch (e) {
        console.error("Error fetching exposures:", e);
    }

    // Input event to filter options
    exposureInput.addEventListener('input', () => {
        const val = exposureInput.value.toLowerCase();
        if (val.length === 0) {
            exposureOptions.style.display = 'none';
            return;
        }

        const matches = exposures.filter(e => e.toLowerCase().includes(val));

        // Populate options
        exposureOptions.innerHTML = '';
        if (matches.length > 0) {
            matches.forEach(match => {
                const div = document.createElement('div');
                div.textContent = match;
                div.addEventListener('click', () => {
                    exposureInput.value = match;
                    exposureOptions.style.display = 'none';
                });
                exposureOptions.appendChild(div);
            });
            exposureOptions.style.display = 'block';
        } else {
            exposureOptions.style.display = 'none';
        }
    });

    // Hide dropdown when clicking outside
    document.addEventListener('click', (e) => {
        if (!exposureInput.contains(e.target) && !exposureOptions.contains(e.target)) {
            exposureOptions.style.display = 'none';
        }
    });

    // Show options on focus if there's input or show all?
    // Let's show search results on focus if existing value?
    exposureInput.addEventListener('focus', () => {
        if (exposureInput.value.trim().length > 0) {
            exposureInput.dispatchEvent(new Event('input'));
        }
    });
});

// Sorting State
let currentSort = { field: null, direction: 'asc' };

document.getElementById('analyze-btn').addEventListener('click', async () => {
    const disease = document.getElementById('disease').value;
    const exposure = document.getElementById('exposure').value;
    const outcome = document.getElementById('outcome').value;
    const excludeMeta = document.getElementById('exclude-meta').checked;
    const loading = document.getElementById('loading');
    const results = document.getElementById('results');
    const errorMsg = document.getElementById('error-message');

    if (!disease || !exposure) {
        alert("Please enter both disease and exposure.");
        return;
    }

    // UI Reset
    loading.classList.remove('hidden');
    results.classList.add('hidden');
    errorMsg.classList.add('hidden');

    try {
        const response = await fetch('/analyze', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ disease, exposure, outcome, exclude_meta: excludeMeta })
        });

        const data = await response.json();

        if (data.error) {
            errorMsg.textContent = data.error;
            errorMsg.classList.remove('hidden');
        } else {
            currentStudies = data.studies; // Store for re-analysis
            updateResultsUI(data);

            // Initial Render
            renderStudiesTable(); // Use the new render function

            results.classList.remove('hidden');
        }

    } catch (e) {
        errorMsg.textContent = "An internal error occurred. Please check the server logs.";
        errorMsg.classList.remove('hidden');
        console.error(e);
    } finally {
        loading.classList.add('hidden');
    }
});

// Re-analyze click handler
document.getElementById('update-btn').addEventListener('click', async () => {
    const disease = document.getElementById('disease').value;
    const exposure = document.getElementById('exposure').value;
    const updateBtn = document.getElementById('update-btn');

    // Get selected indices based on CURRENT rendering order
    // But since we sort, indices might change visually.
    // We should rely on unique IDs or map back to object.
    // For simplicity, let's just grab the displayed objects again?
    // Or simpler: The checkboxes can store the ORIGINAL index? 
    // NO, if we sort, we want to select the study itself.

    // Better approach:
    // Store 'selected' state in currentStudies objects?
    // Or iterate over the checkboxes and find the study in currentStudies that matches?
    // Since we don't have unique IDs, we rely on array index?
    // If we sort currentStudies array directly, then indices match the new order.

    // STRATEGY: We sort `currentStudies` in place (or copy). 
    // `renderStudiesTable` renders `currentStudies`.
    // Checkboxes correspond to `currentStudies[i]`.

    const checkboxes = document.querySelectorAll('.study-checkbox');
    const selectedStudies = [];

    checkboxes.forEach((cb, index) => {
        if (cb.checked) {
            selectedStudies.push(currentStudies[index]);
        }
    });

    if (selectedStudies.length === 0) {
        alert("Please select at least one study to analyze.");
        return;
    }

    updateBtn.textContent = "Updating...";
    updateBtn.disabled = true;

    try {
        const response = await fetch('/reanalyze', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                studies: selectedStudies,
                disease: disease,
                exposure: exposure
            })
        });

        const data = await response.json();

        if (data.error) {
            alert(data.error);
        } else {
            updateResultsUI(data);
        }

    } catch (e) {
        console.error(e);
        alert("Error updating analysis.");
    } finally {
        updateBtn.textContent = "Update Analysis";
        updateBtn.disabled = false;
    }
});

// Render Function
function renderStudiesTable() {
    const tbody = document.querySelector('#studies-table tbody');
    tbody.innerHTML = '';

    currentStudies.forEach((study, index) => {
        const tr = document.createElement('tr');
        tr.innerHTML = `
            <td><input type="checkbox" class="study-checkbox" data-index="${index}" checked></td>
            <td>${index + 1}</td>
            <td>${study.Study}</td>
            <td>${study['Effect Type'] ? study['Effect Type'] + ': ' : ''}${study['Effect Size']}</td>
            <td>${study['Lower CI']}, ${study['Upper CI']}</td>
            <td>${study.Population}</td>
            <td>${study['Sample Size'] || 'N/A'}</td>
            <td style="font-size: 0.85em; opacity: 0.8;">${study.Reference}</td>
            <td>${study.Journal || '-'} (${study.Year || '-'})</td>
            <td><a href="${study.Link}" target="_blank" style="color: var(--accent); text-decoration: none; font-weight: 600;">Link</a></td>
        `;
        tbody.appendChild(tr);
    });

    // Update header icons
    updateSortIcons();
}

// Sorting Function
window.handleSort = (field) => {
    // Toggle direction
    if (currentSort.field === field) {
        currentSort.direction = currentSort.direction === 'asc' ? 'desc' : 'asc';
    } else {
        currentSort.field = field;
        currentSort.direction = 'asc'; // Default new sort to asc
    }

    currentStudies.sort((a, b) => {
        let valA, valB;

        if (field === 'Effect Size') {
            valA = parseFloat(a['Effect Size']);
            valB = parseFloat(b['Effect Size']);
        } else if (field === 'Year') {
            // Handle "Unknown" or strings
            valA = parseInt(a['Year']) || 0;
            valB = parseInt(b['Year']) || 0;
        }

        if (valA < valB) return currentSort.direction === 'asc' ? -1 : 1;
        if (valA > valB) return currentSort.direction === 'asc' ? 1 : -1;
        return 0;
    });

    renderStudiesTable();
};

function updateSortIcons() {
    // Reset classes
    document.querySelectorAll('th.sortable').forEach(th => {
        th.classList.remove('sort-asc', 'sort-desc');
    });

    // Add class to active
    let id = '';
    if (currentSort.field === 'Effect Size') id = 'th-es';
    if (currentSort.field === 'Year') id = 'th-year';

    if (id) {
        const th = document.getElementById(id);
        if (th) th.classList.add(currentSort.direction === 'asc' ? 'sort-asc' : 'sort-desc');
    }
}


function updateResultsUI(data) {
    // Update Plot
    // Force cache bust update
    const imgInfo = data.plot_url.split('?');
    const finalUrl = `/${imgInfo[0]}?t=${new Date().getTime()}`;
    document.getElementById('forest-plot').src = finalUrl;

    // Update Headline
    if (data.headline) {
        const hl = document.getElementById('headline-result');
        hl.classList.remove('hidden');
        document.getElementById('pooled-es').textContent = data.headline.pooled_es;
        document.getElementById('pooled-ci').textContent = `${data.headline.ci_low}, ${data.headline.ci_upp}`;

        const interpEl = document.getElementById('interpretation');
        interpEl.textContent = data.headline.interpretation;

        // Color coding
        if (data.headline.interpretation.includes("Not")) {
            interpEl.style.color = "#8b949e"; // Grey for null
        } else if (data.headline.interpretation.includes("Increased")) {
            interpEl.style.color = "#ff7b72"; // Red for risk
        } else {
            interpEl.style.color = "#238636"; // Green for protective (assuming disease risk)
        }
    }
}
