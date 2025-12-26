let currentStudies = [];

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

            // Generate Table with Checkboxes
            const tbody = document.querySelector('#studies-table tbody');
            tbody.innerHTML = '';

            data.studies.forEach((study, index) => {
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
    const errorMsg = document.getElementById('error-message');
    const updateBtn = document.getElementById('update-btn');

    // Get selected indices
    const checkboxes = document.querySelectorAll('.study-checkbox'); // Fixed selector class
    const selectedStudies = [];

    checkboxes.forEach(cb => {
        if (cb.checked) {
            const index = parseInt(cb.getAttribute('data-index'));
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

function updateResultsUI(data) {
    // Update Plot
    // Force cache bust update
    const imgInfo = data.plot_url.split('?');
    const finalUrl = `/${imgInfo[0]}?t=${new Date().getTime()}`;
    document.getElementById('forest-plot').src = finalUrl;

    // Update Summary
    // Since summary_html is full HTML table, we can just dump it
    // Wait, we removed 'summary-stats' div previously?
    // User asked to "exclude the 'Statistical Details' table".
    // So 'summary-stats' div might be gone or empty?
    // Let's check index.html from previous edit.
    // The div with id 'summary-stats' was REMOVED.
    // So we don't need to update it!

    // We only update the Headline results.

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
