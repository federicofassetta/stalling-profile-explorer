Put your exported data here (`R/export_for_web.R` writes `sequences.tsv.gz`,
`profiles.tsv.gz` and `manifest.json` into this folder).

Everything in this folder except this note is ignored by git, so unpublished
data is never pushed to GitHub. When the site is opened through a local web
server (`python3 -m http.server 8000` from the repository root) it loads the
data from here automatically.
