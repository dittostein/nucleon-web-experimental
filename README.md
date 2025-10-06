# Nucleon Web Experimental

This project is a browser-based "Universe Sandbox" for nucleons and electrons. The backend is a tiny Flask server that hosts the static client and offers optional JSON APIs for saving and loading scenarios, while the frontend performs all rendering and physics in the browser.

## Getting started

1. Create a virtual environment and install dependencies:

   ```bash
   python -m venv .venv
   source .venv/bin/activate
   pip install -r requirements.txt
   ```

   For local experimentation you can also simply `pip install flask`.

2. Launch the Flask development server:

   ```bash
   flask --app app run
   ```

   The simulator will be available at <http://127.0.0.1:5000>.

## Frontend overview

* Fixed-timestep velocity-Verlet integrator with electron sub-stepping.
* Softened Coulomb and toy nuclear forces with per-pair clamping for stability.
* Toolbar-driven user interface with drag-to-fling placement, delete tool, presets, pause/step/reset, and optional overlays.
* Orbit assist feature and dedicated orbiting-electron placement tool to help build quasi-stable atoms.
* Boundary flashes and kill band for wall deletions, plus optional debug panel for live constant tuning.

## Backend API

The backend exposes two optional JSON endpoints:

* `POST /api/save_scenario` → `{ "name", "description", "constants", "particles" }` → `{ "id" }`
* `GET  /api/load_scenario?id=...` → saved scenario payload or 404

Scenarios are kept in-memory for the development server to keep the MVP simple.

## Development notes

* Physics buffers are stored in typed arrays for efficiency and determinism.
* Orbit clustering uses single-link union-find grouping among nucleons with a configurable link distance.
* Electron drag and velocity blending are applied on each substep to avoid runaway orbits.
* Constants can be tweaked from the debug panel and reset to defaults.
