# Full Circle website

Static site (no build step): `index.html`, `assets/style.css`, `assets/app.js`, and the JSON in
`data/` written by `python/scripts/export_website.py` from `python/results/final/`.

    cd python && uv run python scripts/export_website.py   # rebuild data/ (runs the export checks)
    python3 -m http.server 8000 -d website                 # view at http://localhost:8000

`data/csv/` holds the per-House result files offered for download (CC BY 4.0).
Featured votes per House: edit `python/scripts/featured_votes.json`, then rerun the export.
