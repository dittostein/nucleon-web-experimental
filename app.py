from __future__ import annotations

import os
import threading
import uuid
from dataclasses import dataclass, field
from typing import Dict, List

from flask import Flask, jsonify, request, send_from_directory


@dataclass
class Scenario:
    name: str
    description: str = ""
    constants: Dict[str, object] = field(default_factory=dict)
    particles: List[Dict[str, object]] = field(default_factory=list)


class ScenarioStore:
    """Thread-safe in-memory scenario storage."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._data: Dict[str, Scenario] = {}

    def save(self, scenario: Scenario) -> str:
        scenario_id = uuid.uuid4().hex
        with self._lock:
            self._data[scenario_id] = scenario
        return scenario_id

    def load(self, scenario_id: str) -> Scenario | None:
        with self._lock:
            return self._data.get(scenario_id)


app = Flask(__name__, static_folder="static", static_url_path="")
store = ScenarioStore()


@app.route("/")
def index() -> object:
    return send_from_directory(app.static_folder, "index.html")


@app.route("/api/save_scenario", methods=["POST"])
def save_scenario() -> object:
    data = request.get_json(force=True, silent=True)
    if not isinstance(data, dict):
        return jsonify({"error": "Invalid JSON payload."}), 400

    name = str(data.get("name", ""))
    description = str(data.get("description", ""))
    constants = data.get("constants", {})
    particles = data.get("particles", [])

    if not name:
        return jsonify({"error": "Scenario name is required."}), 400

    if not isinstance(constants, dict) or not isinstance(particles, list):
        return jsonify({"error": "Invalid scenario structure."}), 400

    scenario = Scenario(name=name, description=description, constants=constants, particles=particles)
    scenario_id = store.save(scenario)
    return jsonify({"id": scenario_id})


@app.route("/api/load_scenario")
def load_scenario() -> object:
    scenario_id = request.args.get("id")
    if not scenario_id:
        return jsonify({"error": "Missing id parameter."}), 400

    scenario = store.load(scenario_id)
    if scenario is None:
        return jsonify({"error": "Scenario not found."}), 404

    return jsonify({
        "name": scenario.name,
        "description": scenario.description,
        "constants": scenario.constants,
        "particles": scenario.particles,
    })


@app.route("/api/ping")
def ping() -> object:
    return jsonify({"status": "ok"})


if __name__ == "__main__":
    port = int(os.environ.get("PORT", "5000"))
    app.run(debug=True, port=port)
