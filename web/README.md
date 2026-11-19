# qqman web app

A static, client-side web app for Manhattan and Q-Q plots, served at [qqman.stephenturner.us](https://qqman.stephenturner.us). It has no build step and no dependencies. Data are parsed and plotted in the browser and never sent anywhere.

## Run locally

ES modules and the Web Worker need to be served over HTTP, so opening `index.html` from disk won't work. From this directory:

```sh
python3 -m http.server 8000
```

Then open <http://localhost:8000>. Add `#example` to the URL to load the gwasResults example on page load.

## Deploy on Netlify

In the site settings, set the base directory to `web`. There is no build command, and `netlify.toml` sets the publish directory and response headers.

## Layout

- `index.html`, `css/style.css`: page and styles (light and dark themes).
- `js/app.js`: UI state, options, interaction, export.
- `js/worker.js`: streaming parsing of files (including `.gz`), URLs and pasted text, plus the simulated example data. Runs off the main thread.
- `js/parse.js`: delimiter and header detection, column storage. No DOM, so it also runs in Node.
- `js/prepare.js`: column auto-detection, data checks, chromosome normalization, sorting, λGC.
- `js/plot.js`: Manhattan and Q-Q drawing against a small painter interface with canvas and SVG back ends, so the screen, PNG and SVG output match.
- `js/stats.js`: normal and beta quantiles, ticks, number formatting.
- `js/rcode.js`: builds the equivalent qqman R call.
- `examples/`: `gwasResults` and `snpsOfInterest` exported from the R package.
