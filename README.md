# Stalling Profile Explorer

An interactive website for exploring replisome coverage (stalling) profiles across a
library of ~20,000 sequences. Each profile counts, at every position of the reference,
the reads whose last incorporated nucleotide is at that position, so peaks mark where the
polymerase stalled.

It is a static site (HTML, CSS and JavaScript, plotting with Plotly), so it runs on GitHub
Pages with no server. The public version shows a **simulated demo dataset**; real data is
loaded from your own computer and never leaves the browser.

## Views

| View | What it shows |
|---|---|
| **Overview** | Mean normalised profile ± SEM for the filtered sequences, one line per sample. Split into panels by any metadata column (e.g. `motif.type`), or flip it to one panel per sample with one line per group. Below: the fraction of each profile's reads that end in the stall window, by group and sample. |
| **Heatmap** | One row per sequence, one column per position, ordered by window fraction, peak position, total reads, the change between two samples, or grouped by a metadata column. With a second sample selected, colour shows B − A. Click a row to open that sequence. |
| **Sequences** | Sortable table of window fraction, peak and reads per sample, and Δ between two samples. Click a row for its profile in every sample, its metadata, and the sequence with the STR, stall window and peak highlighted. ← / → step through the list. |

The sidebar holds sample selection, filters (free-text search plus one filter per metadata
column, built automatically from whatever columns you export), a minimum-reads cut-off,
the normalisation, the STR start position and the stall window.

Plots download as SVG from the camera icon in each plot's toolbar; the overview summary and the
sequence table download as CSV.

## Using your own data

1. Build `cov_long` in R as usual (read the `.depth` files, join `STR.info.light.df`).
2. Run `R/export_for_web.R` in the same session. It writes
   `sequences.tsv.gz`, `profiles.tsv.gz` and `manifest.json` into `data/local/`.
3. Then either:
   - **Drop the two `.tsv.gz` files onto the page** (local copy or the GitHub Pages site).
     They are read in the browser; nothing is uploaded.
   - **Or serve the folder locally** so the data loads automatically:
     ```bash
     cd ~/Projects/stalling-profile-explorer
     python3 -m http.server 8000
     ```
     and open <http://localhost:8000>.

`data/local/` is in `.gitignore`, so the real data is never committed.

### File formats

- `sequences` — `seq.id` plus any metadata columns. A column called `seq` is shown as the
  sequence. Columns with up to 40 distinct values become checkbox filters; other numeric
  columns become range filters; everything else is searchable.
- `profiles` (wide) — `seq.id`, `sample`, `p1` … `p150`.
- Also accepted: a long table with `seq.id`, `pos`, `depth`, `sample`, or one `*.depth`
  column per sample. Tab- or comma-separated, plain or gzip-compressed.

### Definitions

- **Distance from STR start**: `dist = STR start − pos` (default start 100), as in the R script.
- **÷ max** normalisation divides each sequence × sample profile by its highest position
  (the `norm` column in R). **÷ total** gives the fraction of the profile's reads at each
  position. Profiles with zero reads are left out of averages.
- **Window fraction**: reads ending inside the stall window ÷ all reads of that profile,
  computed from raw counts. **Δ** = window fraction in B − in A.

## Publishing on GitHub Pages

```bash
cd ~/Projects/stalling-profile-explorer
git init -b main
git add .
git commit -m "Stalling Profile Explorer"
git remote add origin https://github.com/federicofassetta/stalling-profile-explorer.git
git push -u origin main
```

Then on GitHub: **Settings → Pages → Build and deployment → Deploy from a branch →
`main` / `(root)`**. The site appears at
`https://federicofassetta.github.io/stalling-profile-explorer/`.

## Structure

```
index.html            page layout
css/style.css         styles (palette shared with the RCARU site)
js/data.js            file parsing and dataset assembly (no DOM; testable in Node)
js/demo.js            seeded simulated demo dataset
js/app.js             views, plots and interaction
R/export_for_web.R    exports cov_long for the website
data/local/           your real data (git-ignored)
```
