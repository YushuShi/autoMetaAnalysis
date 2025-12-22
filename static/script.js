document.getElementById('analyze-btn').addEventListener('click', async () => {
    const disease = document.getElementById('disease').value;
    const exposure = document.getElementById('exposure').value;
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
            body: JSON.stringify({ disease, exposure })
        });

        const data = await response.json();

        if (data.error) {
            errorMsg.textContent = data.error;
            errorMsg.classList.remove('hidden');
        } else {
            // Update Plot
            document.getElementById('forest-plot').src = data.plot_url;

            // Update Summary
            document.getElementById('summary-stats').innerHTML = data.summary_html;

            // Highlight 'random effect wls' row - UPDATED to look for 'Random-effects meta-analysis (WLS)'
            const summaryTable = document.querySelector('#summary-stats table');
            if (summaryTable) {
                const rows = summaryTable.querySelectorAll('tr');
                rows.forEach(row => {
                    const firstCell = row.cells[0];
                    // Case insensitive check for the new label
                    if (firstCell && firstCell.textContent.trim().toLowerCase().includes('random-effects meta-analysis (wls)')) {
                        row.classList.add('highlight-row');
                    }
                });
            }

            // Update Headline
            if (data.headline) {
                const hl = document.getElementById('headline-result');
                hl.classList.remove('hidden');
                document.getElementById('pooled-es').textContent = data.headline.pooled_es;
                document.getElementById('pooled-ci').textContent = `${data.headline.ci_low} - ${data.headline.ci_upp}`;

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

            // Update Table
            const tbody = document.querySelector('#studies-table tbody');
            tbody.innerHTML = '';

            data.studies.forEach((study, index) => {
                const tr = document.createElement('tr');
                tr.innerHTML = `
                    <td>${index + 1}</td>
                    <td>${study.Study}</td>
                    <td>${study['Effect Size']}</td>
                    <td>${study['Lower CI']} - ${study['Upper CI']}</td>
                    <td>${study.Population}</td>
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
