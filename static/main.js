const Tool = Object.freeze({
  PROTON: "proton",
  NEUTRON: "neutron",
  ELECTRON: "electron",
  ORBITING_ELECTRON: "orbiting-electron",
  DELETE: "delete",
});

const ParticleType = Object.freeze({
  PROTON: 0,
  NEUTRON: 1,
  ELECTRON: 2,
});

const DEFAULT_CONSTANTS = Object.freeze({
  masses: { proton: 1.0, neutron: 1.0, electron: 0.05 },
  charges: { proton: 1.0, neutron: 0.0, electron: -1.0 },
  coulomb: { kC: 1.01e5, epsC: 4.0 },
  nuclear: {
    A: 10000.0,
    lambda_attract: 15.0,
    B: 8000.0,
    lambda_repulse: 5.0,
  },
  integrator: { dt: 1 / 240, electron_substeps: 4, velverlet: true },
  drag: { electron_gamma: 0.003 },
  caps: { Fmax: 3000.0, v_spawn_cap: 300.0 },
  spawn: { k_drag: 2.0 },
  orbit: {
    assist_eta_per_dt: 0.02,
    cluster_link_distance: 20.0,
    nucleus_pick_radius: 60.0,
  },
});

const PRESETS = [
  {
    name: "Hydrogen-like orbit",
    description: "Single proton with one orbiting electron",
    setup: (system, helpers) => {
      system.clear();
      const centerX = system.bounds.width / 2;
      const centerY = system.bounds.height / 2;
      system.addParticle(ParticleType.PROTON, centerX, centerY, 0, 0);
      const electronX = centerX + 50;
      const electronY = centerY;
      const r = 50;
      const { kC, epsC } = system.constants.coulomb;
      const mass = system.constants.masses.electron;
      const charge = Math.abs(
        system.constants.charges.electron * system.constants.charges.proton
      );
      const denom = Math.pow(r * r + epsC * epsC, 1.5);
      const speed = Math.sqrt((kC * charge * r * r) / (mass * denom));
      system.addParticle(ParticleType.ELECTRON, electronX, electronY, 0, speed);
    },
  },
  {
    name: "p-n dimer",
    description: "Proton and neutron in a nuclear bound state",
    setup: (system) => {
      system.clear();
      const cx = system.bounds.width / 2;
      const cy = system.bounds.height / 2;
      system.addParticle(ParticleType.PROTON, cx - 5, cy, 0, 0);
      system.addParticle(ParticleType.NEUTRON, cx + 5, cy, 0, 0);
    },
  },
  {
    name: "Two protons repelling",
    description: "Repulsion test",
    setup: (system) => {
      system.clear();
      const cx = system.bounds.width / 2;
      const cy = system.bounds.height / 2;
      system.addParticle(ParticleType.PROTON, cx - 20, cy, 0, 0);
      system.addParticle(ParticleType.PROTON, cx + 20, cy, 0, 0);
    },
  },
  {
    name: "Electron slingshot",
    description: "Fast electron hitting boundary",
    setup: (system) => {
      system.clear();
      const y = system.bounds.height / 2;
      system.addParticle(ParticleType.ELECTRON, 100, y, 250, 0);
    },
  },
  {
    name: "Three-body chaos",
    description: "Chaotic nucleon cluster",
    setup: (system) => {
      system.clear();
      const cx = system.bounds.width / 2;
      const cy = system.bounds.height / 2;
      system.addParticle(ParticleType.PROTON, cx - 25, cy - 10, 20, 0);
      system.addParticle(ParticleType.NEUTRON, cx + 15, cy + 10, -15, 10);
      system.addParticle(ParticleType.PROTON, cx + 20, cy - 20, 5, -10);
      system.addParticle(ParticleType.ELECTRON, cx, cy + 50, -40, 0);
    },
  },
];

const COLORS = {
  [ParticleType.PROTON]: "#ff4d4f",
  [ParticleType.NEUTRON]: "#b5bac6",
  [ParticleType.ELECTRON]: "#3fa9f5",
};

const RADII = {
  [ParticleType.PROTON]: 5,
  [ParticleType.NEUTRON]: 5,
  [ParticleType.ELECTRON]: 3.5,
};

const OUTLINE_COLOR = "rgba(255,255,255,0.25)";

class ParticleSystem {
  constructor(constants, bounds) {
    this.constants = JSON.parse(JSON.stringify(constants));
    this.bounds = bounds;
    this.capacity = 128;
    this.count = 0;
    this.killBand = 2;
    this._ensureBuffers();
    this.tempAx = new Float32Array(this.capacity);
    this.tempAy = new Float32Array(this.capacity);
    this.lastClusterFrame = -1;
    this.clusters = [];
    this.clusterLookup = new Map();
    this.frameIndex = 0;
    this.orbitAssistEnabled = true;
    this.orbitChirality = 1;
  }

  _ensureBuffers() {
    if (this.x) return;
    const cap = this.capacity;
    this.x = new Float32Array(cap);
    this.y = new Float32Array(cap);
    this.vx = new Float32Array(cap);
    this.vy = new Float32Array(cap);
    this.ax = new Float32Array(cap);
    this.ay = new Float32Array(cap);
    this.mass = new Float32Array(cap);
    this.charge = new Float32Array(cap);
    this.type = new Int8Array(cap);
    this.alive = new Uint8Array(cap);
  }

  _grow() {
    this.capacity *= 2;
    this.x = ParticleSystem._growArray(this.x, this.capacity);
    this.y = ParticleSystem._growArray(this.y, this.capacity);
    this.vx = ParticleSystem._growArray(this.vx, this.capacity);
    this.vy = ParticleSystem._growArray(this.vy, this.capacity);
    this.ax = ParticleSystem._growArray(this.ax, this.capacity);
    this.ay = ParticleSystem._growArray(this.ay, this.capacity);
    this.mass = ParticleSystem._growArray(this.mass, this.capacity);
    this.charge = ParticleSystem._growArray(this.charge, this.capacity);
    this.type = ParticleSystem._growArray(this.type, this.capacity, Int8Array);
    this.alive = ParticleSystem._growArray(this.alive, this.capacity, Uint8Array);
    this.tempAx = ParticleSystem._growArray(this.tempAx, this.capacity);
    this.tempAy = ParticleSystem._growArray(this.tempAy, this.capacity);
  }

  static _growArray(arr, capacity, ctor) {
    const Ctor = ctor ?? arr.constructor;
    const next = new Ctor(capacity);
    next.set(arr);
    return next;
  }

  clear() {
    this.count = 0;
    this.alive.fill(0);
    this.frameIndex = 0;
  }

  addParticle(type, x, y, vx, vy) {
    if (this.count >= this.capacity) {
      this._grow();
    }
    const idx = this.count++;
    this.x[idx] = x;
    this.y[idx] = y;
    this.vx[idx] = vx;
    this.vy[idx] = vy;
    this.ax[idx] = 0;
    this.ay[idx] = 0;
    this.type[idx] = type;
    const masses = this.constants.masses;
    const charges = this.constants.charges;
    switch (type) {
      case ParticleType.PROTON:
        this.mass[idx] = masses.proton;
        this.charge[idx] = charges.proton;
        break;
      case ParticleType.NEUTRON:
        this.mass[idx] = masses.neutron;
        this.charge[idx] = charges.neutron;
        break;
      case ParticleType.ELECTRON:
        this.mass[idx] = masses.electron;
        this.charge[idx] = charges.electron;
        break;
      default:
        throw new Error("Unknown particle type");
    }
    this.alive[idx] = 1;
    return idx;
  }

  removeParticle(index) {
    if (index < 0 || index >= this.count || !this.alive[index]) return;
    const last = this.count - 1;
    if (index !== last) {
      this._swap(index, last);
    }
    this.alive[last] = 0;
    this.count = Math.max(0, this.count - 1);
  }

  _swap(i, j) {
    if (i === j) return;
    const keys = ["x", "y", "vx", "vy", "ax", "ay", "mass", "charge", "type", "alive"];
    for (const key of keys) {
      const arr = this[key];
      const tmp = arr[i];
      arr[i] = arr[j];
      arr[j] = tmp;
    }
  }

  computeAccelerations() {
    const n = this.count;
    const ax = this.ax;
    const ay = this.ay;
    for (let i = 0; i < n; i++) {
      ax[i] = 0;
      ay[i] = 0;
    }

    const { kC, epsC } = this.constants.coulomb;
    const { A, lambda_attract: la, B, lambda_repulse: lr } = this.constants.nuclear;
    const Fmax = this.constants.caps.Fmax;

    for (let i = 0; i < n - 1; i++) {
      if (!this.alive[i]) continue;
      const xi = this.x[i];
      const yi = this.y[i];
      for (let j = i + 1; j < n; j++) {
        if (!this.alive[j]) continue;
        let dx = this.x[j] - xi;
        let dy = this.y[j] - yi;
        let r2 = dx * dx + dy * dy;
        const eps2 = epsC * epsC;
        const rSoft = Math.sqrt(r2 + eps2);
        const chargeProduct = this.charge[i] * this.charge[j];
        let fx = 0;
        let fy = 0;
        if (rSoft > 0) {
          const fc = (kC * chargeProduct) / (rSoft * rSoft * rSoft);
          fx += fc * dx;
          fy += fc * dy;
        }

        const typeI = this.type[i];
        const typeJ = this.type[j];
        const nucleonI = typeI === ParticleType.PROTON || typeI === ParticleType.NEUTRON;
        const nucleonJ = typeJ === ParticleType.PROTON || typeJ === ParticleType.NEUTRON;
        if (nucleonI && nucleonJ && r2 > 0) {
          const r = Math.sqrt(r2);
          const fr = -(A / la) * Math.exp(-r / la) + (B / lr) * Math.exp(-r / lr);
          if (r > 0) {
            const invR = 1 / r;
            fx += fr * dx * invR;
            fy += fr * dy * invR;
          }
        }

        const mag = Math.hypot(fx, fy);
        if (mag > Fmax && mag > 0) {
          const scale = Fmax / mag;
          fx *= scale;
          fy *= scale;
        }

        const invMi = 1 / this.mass[i];
        const invMj = 1 / this.mass[j];
        ax[i] -= fx * invMi;
        ay[i] -= fy * invMi;
        ax[j] += fx * invMj;
        ay[j] += fy * invMj;
      }
    }
  }

  integrate(dt, orbitManager, flashes) {
    const { electron_substeps } = this.constants.integrator;
    const dtSub = dt / electron_substeps;
    const dragGamma = this.constants.drag.electron_gamma;
    const dragFactor = Math.exp(-dragGamma * dtSub);
    const etaSub = this.orbitAssistEnabled
      ? this.constants.orbit.assist_eta_per_dt / electron_substeps
      : 0;

    this.computeAccelerations();
    const n = this.count;
    const oldAx = this.tempAx;
    const oldAy = this.tempAy;
    for (let i = 0; i < n; i++) {
      oldAx[i] = this.ax[i];
      oldAy[i] = this.ay[i];
    }

    // Update nucleon positions first
    for (let i = 0; i < n; i++) {
      if (!this.alive[i]) continue;
      const type = this.type[i];
      if (type === ParticleType.ELECTRON) continue;
      const axi = oldAx[i];
      const ayi = oldAy[i];
      this.x[i] += this.vx[i] * dt + 0.5 * axi * dt * dt;
      this.y[i] += this.vy[i] * dt + 0.5 * ayi * dt * dt;
    }

    // Electron substeps
    for (let step = 0; step < electron_substeps; step++) {
      for (let i = 0; i < n; i++) {
        if (!this.alive[i]) continue;
        if (this.type[i] !== ParticleType.ELECTRON) continue;
        const axi = oldAx[i];
        const ayi = oldAy[i];
        this.x[i] += this.vx[i] * dtSub + 0.5 * axi * dtSub * dtSub;
        this.y[i] += this.vy[i] * dtSub + 0.5 * ayi * dtSub * dtSub;
      }

      this.computeAccelerations();

      for (let i = 0; i < n; i++) {
        if (!this.alive[i]) continue;
        const type = this.type[i];
        if (type === ParticleType.ELECTRON) {
          const axiPrev = oldAx[i];
          const ayiPrev = oldAy[i];
          const axi = this.ax[i];
          const ayi = this.ay[i];
          this.vx[i] += 0.5 * (axiPrev + axi) * dtSub;
          this.vy[i] += 0.5 * (ayiPrev + ayi) * dtSub;
          this.vx[i] *= dragFactor;
          this.vy[i] *= dragFactor;

          if (etaSub > 0) {
            const assist = orbitManager.computeAssist(i);
            if (assist) {
              const targetVx = assist.targetVx * this.orbitChirality;
              const targetVy = assist.targetVy * this.orbitChirality;
              this.vx[i] = (1 - etaSub) * this.vx[i] + etaSub * targetVx;
              this.vy[i] = (1 - etaSub) * this.vy[i] + etaSub * targetVy;
            }
          }

          oldAx[i] = axi;
          oldAy[i] = ayi;
        }
      }
    }

    // After final electron substep, velocities for nucleons
    for (let i = 0; i < n; i++) {
      if (!this.alive[i]) continue;
      const type = this.type[i];
      if (type === ParticleType.ELECTRON) continue;
      this.vx[i] += 0.5 * (oldAx[i] + this.ax[i]) * dt;
      this.vy[i] += 0.5 * (oldAy[i] + this.ay[i]) * dt;
    }

    this.frameIndex++;
    this.handleBounds(flashes);
  }

  handleBounds(flashes) {
    const minX = 0 - this.killBand;
    const minY = 0 - this.killBand;
    const maxX = this.bounds.width + this.killBand;
    const maxY = this.bounds.height + this.killBand;

    for (let i = this.count - 1; i >= 0; i--) {
      if (!this.alive[i]) continue;
      const x = this.x[i];
      const y = this.y[i];
      if (x < minX || x > maxX || y < minY || y > maxY) {
        flashes.spawn(this.x[i], this.y[i]);
        this.removeParticle(i);
      }
    }
  }

  findNearest(x, y, radius) {
    let bestIndex = -1;
    let bestDist2 = radius * radius;
    for (let i = 0; i < this.count; i++) {
      if (!this.alive[i]) continue;
      const dx = this.x[i] - x;
      const dy = this.y[i] - y;
      const dist2 = dx * dx + dy * dy;
      if (dist2 <= bestDist2) {
        bestIndex = i;
        bestDist2 = dist2;
      }
    }
    return bestIndex;
  }
}

class OrbitManager {
  constructor(system) {
    this.system = system;
    this.clusters = [];
  }

  updateClusters() {
    const system = this.system;
    this.clusters = [];

    const linkDistance = system.constants.orbit.cluster_link_distance;
    const link2 = linkDistance * linkDistance;

    const nucleonIndices = [];
    for (let i = 0; i < system.count; i++) {
      if (!system.alive[i]) continue;
      const type = system.type[i];
      if (type === ParticleType.PROTON || type === ParticleType.NEUTRON) {
        nucleonIndices.push(i);
      }
    }

    const n = nucleonIndices.length;
    if (n === 0) {
      return;
    }

    const parent = new Int32Array(n);
    for (let i = 0; i < n; i++) parent[i] = i;

    const find = (i) => {
      while (parent[i] !== i) {
        parent[i] = parent[parent[i]];
        i = parent[i];
      }
      return i;
    };
    const union = (a, b) => {
      const ra = find(a);
      const rb = find(b);
      if (ra !== rb) parent[rb] = ra;
    };

    for (let i = 0; i < n - 1; i++) {
      const idxI = nucleonIndices[i];
      const xi = system.x[idxI];
      const yi = system.y[idxI];
      for (let j = i + 1; j < n; j++) {
        const idxJ = nucleonIndices[j];
        const dx = system.x[idxJ] - xi;
        const dy = system.y[idxJ] - yi;
        const dist2 = dx * dx + dy * dy;
        if (dist2 <= link2) {
          union(i, j);
        }
      }
    }

    const clusters = new Map();
    for (let i = 0; i < n; i++) {
      const root = find(i);
      if (!clusters.has(root)) {
        clusters.set(root, {
          indices: [],
          massSum: 0,
          charge: 0,
          comX: 0,
          comY: 0,
        });
      }
      const data = clusters.get(root);
      const idx = nucleonIndices[i];
      data.indices.push(idx);
      data.massSum += system.mass[idx];
      if (system.type[idx] === ParticleType.PROTON) {
        data.charge += 1;
      }
      data.comX += system.mass[idx] * system.x[idx];
      data.comY += system.mass[idx] * system.y[idx];
    }

    this.clusters = Array.from(clusters.values()).map((cluster) => {
      const invMass = cluster.massSum > 0 ? 1 / cluster.massSum : 0;
      cluster.comX *= invMass;
      cluster.comY *= invMass;
      return cluster;
    });
  }

  computeAssist(electronIndex) {
    this.updateClusters();
    const system = this.system;
    const pickRadius = system.constants.orbit.nucleus_pick_radius;
    const pick2 = pickRadius * pickRadius;
    const ex = system.x[electronIndex];
    const ey = system.y[electronIndex];

    let bestCluster = null;
    let bestDist2 = pick2;
    for (const cluster of this.clusters) {
      if (cluster.charge <= 0) continue;
      const dx = ex - cluster.comX;
      const dy = ey - cluster.comY;
      const dist2 = dx * dx + dy * dy;
      if (dist2 < bestDist2) {
        bestDist2 = dist2;
        bestCluster = cluster;
      }
    }

    if (!bestCluster) {
      let bestProton = -1;
      for (let i = 0; i < system.count; i++) {
        if (!system.alive[i]) continue;
        if (system.type[i] !== ParticleType.PROTON) continue;
        const dx = ex - system.x[i];
        const dy = ey - system.y[i];
        const dist2 = dx * dx + dy * dy;
        if (dist2 < bestDist2) {
          bestDist2 = dist2;
          bestProton = i;
        }
      }
      if (bestProton >= 0) {
        bestCluster = {
          comX: system.x[bestProton],
          comY: system.y[bestProton],
          charge: 1,
        };
      }
    }

    if (!bestCluster) {
      return null;
    }

    const dx = ex - bestCluster.comX;
    const dy = ey - bestCluster.comY;
    const r2 = dx * dx + dy * dy;
    if (r2 === 0) {
      return null;
    }
    const r = Math.sqrt(r2);
    const { kC, epsC } = system.constants.coulomb;
    const m = system.mass[electronIndex];
    const q = Math.abs(bestCluster.charge * system.charge[electronIndex]);
    const denom = Math.pow(r2 + epsC * epsC, 1.5);
    const speed = Math.sqrt((kC * q * r2) / (m * denom));
    const invR = 1 / r;
    const perpX = -dy * invR;
    const perpY = dx * invR;

    const result = {
      targetVx: perpX * speed,
      targetVy: perpY * speed,
      speed,
    };
    return result;
  }
}

class FlashManager {
  constructor() {
    this.entries = [];
    this.duration = 0.1; // seconds
  }

  spawn(x, y) {
    const now = performance.now() / 1000;
    this.entries.push({ x, y, born: now });
  }

  prune() {
    const now = performance.now() / 1000;
    this.entries = this.entries.filter((entry) => now - entry.born < this.duration);
  }
}

class Renderer {
  constructor(canvas, system, flashes, orbitManager) {
    this.canvas = canvas;
    this.ctx = canvas.getContext("2d");
    this.system = system;
    this.flashes = flashes;
    this.orbitManager = orbitManager;
    this.overlayVelocity = false;
    this.overlayAccel = false;
    this.overlayCom = false;
    this.dragState = null;
  }

  setDragState(state) {
    this.dragState = state;
  }

  render(dt, accumulator) {
    const ctx = this.ctx;
    const { width, height } = this.canvas;
    ctx.clearRect(0, 0, width, height);

    this.drawKillBand(ctx);
    this.drawParticles(ctx);
    if (this.overlayVelocity) this.drawVelocityVectors(ctx);
    if (this.overlayAccel) this.drawAccelerationVectors(ctx);
    if (this.overlayCom) this.drawClusterCenters(ctx);
    this.drawFlashes(ctx);
    this.drawDragArrow(ctx);
  }

  drawKillBand(ctx) {
    const band = this.system.killBand;
    ctx.save();
    ctx.strokeStyle = "rgba(255,255,255,0.05)";
    ctx.lineWidth = band * 2;
    ctx.strokeRect(band, band, this.canvas.width - band * 2, this.canvas.height - band * 2);
    ctx.restore();
  }

  drawParticles(ctx) {
    const system = this.system;
    for (let i = 0; i < system.count; i++) {
      if (!system.alive[i]) continue;
      const type = system.type[i];
      const radius = RADII[type];
      ctx.beginPath();
      ctx.fillStyle = COLORS[type];
      ctx.arc(system.x[i], system.y[i], radius, 0, Math.PI * 2);
      ctx.fill();
      ctx.lineWidth = 1;
      ctx.strokeStyle = OUTLINE_COLOR;
      ctx.stroke();
    }
  }

  drawVelocityVectors(ctx) {
    const system = this.system;
    ctx.save();
    ctx.strokeStyle = "rgba(63,169,245,0.6)";
    ctx.lineWidth = 1;
    for (let i = 0; i < system.count; i++) {
      if (!system.alive[i]) continue;
      ctx.beginPath();
      ctx.moveTo(system.x[i], system.y[i]);
      ctx.lineTo(system.x[i] + system.vx[i] * 0.1, system.y[i] + system.vy[i] * 0.1);
      ctx.stroke();
    }
    ctx.restore();
  }

  drawAccelerationVectors(ctx) {
    const system = this.system;
    ctx.save();
    ctx.strokeStyle = "rgba(255,146,0,0.6)";
    ctx.lineWidth = 1;
    for (let i = 0; i < system.count; i++) {
      if (!system.alive[i]) continue;
      ctx.beginPath();
      ctx.moveTo(system.x[i], system.y[i]);
      ctx.lineTo(system.x[i] + system.ax[i] * 0.01, system.y[i] + system.ay[i] * 0.01);
      ctx.stroke();
    }
    ctx.restore();
  }

  drawClusterCenters(ctx) {
    this.orbitManager.updateClusters();
    ctx.save();
    ctx.fillStyle = "rgba(255,255,255,0.8)";
    for (const cluster of this.orbitManager.clusters) {
      ctx.beginPath();
      ctx.arc(cluster.comX, cluster.comY, 3, 0, Math.PI * 2);
      ctx.fill();
      ctx.fillText(`+${cluster.charge}`, cluster.comX + 4, cluster.comY - 4);
    }
    ctx.restore();
  }

  drawFlashes(ctx) {
    const now = performance.now() / 1000;
    const duration = this.flashes.duration;
    for (const flash of this.flashes.entries) {
      const age = now - flash.born;
      if (age >= duration) continue;
      const t = age / duration;
      const radius = 6 + 20 * t;
      const alpha = 1 - t;
      ctx.beginPath();
      ctx.strokeStyle = `rgba(255,255,255,${alpha.toFixed(3)})`;
      ctx.lineWidth = 2;
      ctx.arc(flash.x, flash.y, radius, 0, Math.PI * 2);
      ctx.stroke();
    }
  }

  drawDragArrow(ctx) {
    if (!this.dragState || !this.dragState.active) return;
    const { startX, startY, currentX, currentY, kDrag, cap } = this.dragState;
    const dx = currentX - startX;
    const dy = currentY - startY;
    const speed = Math.min(Math.hypot(dx, dy) * kDrag, cap);
    const dir = Math.atan2(dy, dx) + Math.PI;
    const length = (speed / cap) * 80;
    const endX = startX + Math.cos(dir) * length;
    const endY = startY + Math.sin(dir) * length;

    ctx.save();
    ctx.strokeStyle = "rgba(255,255,255,0.8)";
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(startX, startY);
    ctx.lineTo(endX, endY);
    ctx.stroke();

    ctx.beginPath();
    ctx.arc(startX, startY, 4, 0, Math.PI * 2);
    ctx.stroke();

    ctx.fillStyle = "rgba(8,10,15,0.85)";
    ctx.font = "12px sans-serif";
    ctx.fillText(`${speed.toFixed(1)} px/s`, endX + 6, endY + 4);
    ctx.restore();
  }
}

class UIController {
  constructor(system, renderer, orbitManager, flashes) {
    this.system = system;
    this.renderer = renderer;
    this.orbitManager = orbitManager;
    this.flashes = flashes;

    this.canvas = document.getElementById("sim-canvas");
    this.fpsLabel = document.getElementById("hud-fps");
    this.countLabel = document.getElementById("hud-particles");
    this.orbitToggle = document.getElementById("orbit-assist-toggle");
    this.debugPanel = document.getElementById("debug-panel");
    this.debugConstants = document.getElementById("debug-constants");
    this.resetConstantsButton = document.getElementById("reset-constants");

    this.toolButtons = Array.from(document.querySelectorAll("#tool-buttons button"));
    this.currentTool = Tool.PROTON;
    this.dragState = { active: false };
    this.renderer.setDragState(this.dragState);

    this.paused = false;
    this.accumulator = 0;
    this.lastTimestamp = performance.now();
    this.deleteInterval = null;

    this.hudArrow = null;

    this.overlayVelocityCheckbox = document.getElementById("overlay-velocity");
    this.overlayAccelCheckbox = document.getElementById("overlay-accel");
    this.overlayComCheckbox = document.getElementById("overlay-com");

    this.initUI();
    this.initDebug();
    this.initKeyboard();
    this.populatePresets();
  }

  initUI() {
    this.toolButtons.forEach((btn) => {
      btn.addEventListener("click", () => {
        this.toolButtons.forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        this.currentTool = btn.dataset.tool;
      });
    });

    document.getElementById("pause-button").addEventListener("click", () => {
      this.paused = !this.paused;
      document.getElementById("pause-button").textContent = this.paused ? "Resume" : "Pause";
    });

    document.getElementById("step-button").addEventListener("click", () => {
      if (!this.paused) return;
      this.stepSimulation();
    });

    document.getElementById("reset-button").addEventListener("click", () => {
      this.system.clear();
    });

    this.orbitToggle.addEventListener("change", () => {
      this.system.orbitAssistEnabled = this.orbitToggle.checked;
    });

    this.overlayVelocityCheckbox.addEventListener("change", () => {
      this.renderer.overlayVelocity = this.overlayVelocityCheckbox.checked;
    });
    this.overlayAccelCheckbox.addEventListener("change", () => {
      this.renderer.overlayAccel = this.overlayAccelCheckbox.checked;
    });
    this.overlayComCheckbox.addEventListener("change", () => {
      this.renderer.overlayCom = this.overlayComCheckbox.checked;
    });

    this.canvas.addEventListener("mousedown", (ev) => this.handlePointerDown(ev));
    window.addEventListener("mousemove", (ev) => this.handlePointerMove(ev));
    window.addEventListener("mouseup", (ev) => this.handlePointerUp(ev));
    this.canvas.addEventListener("mouseleave", () => this.handlePointerUp());

    this.resetConstantsButton.addEventListener("click", () => {
      this.resetConstants();
    });
  }

  initKeyboard() {
    window.addEventListener("keydown", (ev) => {
      if (ev.key === " ") {
        ev.preventDefault();
        this.paused = !this.paused;
        document.getElementById("pause-button").textContent = this.paused ? "Resume" : "Pause";
      } else if (ev.key === "s" || ev.key === "S") {
        ev.preventDefault();
        if (this.paused) this.stepSimulation();
      } else if (ev.key === "d" || ev.key === "D") {
        this.debugPanel.classList.toggle("hidden");
      } else if (ev.key === "o" || ev.key === "O") {
        this.orbitToggle.checked = !this.orbitToggle.checked;
        this.system.orbitAssistEnabled = this.orbitToggle.checked;
      } else if (ev.key === "c" || ev.key === "C") {
        this.overlayVelocityCheckbox.checked = !this.overlayVelocityCheckbox.checked;
        this.renderer.overlayVelocity = this.overlayVelocityCheckbox.checked;
      } else if (ev.key === "r" || ev.key === "R") {
        this.system.clear();
      } else if (ev.key === "f" || ev.key === "F") {
        this.system.orbitChirality *= -1;
      }
    });
  }

  initDebug() {
    const constContainer = document.createElement("div");
    constContainer.className = "constants-container";
    const buildInputs = (obj, path = []) => {
      const entries = Object.entries(obj);
      return entries
        .map(([key, value]) => {
          const currentPath = [...path, key];
          if (typeof value === "object" && value !== null) {
            const wrapper = document.createElement("div");
            wrapper.innerHTML = `<strong>${currentPath.join(".")}</strong>`;
            const child = buildInputs(value, currentPath);
            child.forEach((el) => wrapper.appendChild(el));
            return wrapper;
          }
          const label = document.createElement("label");
          label.textContent = `${currentPath.join(".")}`;
          const input = document.createElement("input");
          input.type = "number";
          input.value = value;
          input.step = "0.001";
          input.addEventListener("change", () => {
            this.setConstant(currentPath, parseFloat(input.value));
          });
          label.appendChild(input);
          return label;
        })
        .flat();
    };
    const inputs = buildInputs(this.system.constants);
    inputs.forEach((el) => constContainer.appendChild(el));
    this.debugConstants.appendChild(constContainer);
  }

  setConstant(path, value) {
    let target = this.system.constants;
    for (let i = 0; i < path.length - 1; i++) {
      target = target[path[i]];
    }
    target[path[path.length - 1]] = value;
  }

  resetConstants() {
    this.system.constants = JSON.parse(JSON.stringify(DEFAULT_CONSTANTS));
    this.debugConstants.innerHTML = "";
    this.initDebug();
  }

  populatePresets() {
    const select = document.getElementById("preset-select");
    PRESETS.forEach((preset, index) => {
      const option = document.createElement("option");
      option.value = index;
      option.textContent = preset.name;
      select.appendChild(option);
    });
    document.getElementById("load-preset").addEventListener("click", () => {
      const index = parseInt(select.value, 10);
      const preset = PRESETS[index];
      if (preset) {
        preset.setup(this.system, {
          computeOrbitTarget: (x, y) => this.computeOrbitTarget(x, y),
        });
      }
    });
  }

  computeOrbitTarget(x, y) {
    const helper = new OrbitManager(this.system);
    helper.updateClusters();
    const electronIndex = this.system.findNearest(x, y, 1);
    if (electronIndex === -1) return null;
    return helper.computeAssist(electronIndex);
  }

  handlePointerDown(ev) {
    const rect = this.canvas.getBoundingClientRect();
    if (rect.width === 0 || rect.height === 0) return;
    const scaleX = this.canvas.width / rect.width;
    const scaleY = this.canvas.height / rect.height;
    const x = (ev.clientX - rect.left) * scaleX;
    const y = (ev.clientY - rect.top) * scaleY;
    if (this.currentTool === Tool.DELETE) {
      this.deleteAt(x, y);
      this.deleteInterval = setInterval(() => this.deleteAt(x, y), 120);
      return;
    }
    this.dragState.active = true;
    this.dragState.startX = x;
    this.dragState.startY = y;
    this.dragState.currentX = x;
    this.dragState.currentY = y;
    this.dragState.kDrag = this.system.constants.spawn.k_drag;
    this.dragState.cap = this.system.constants.caps.v_spawn_cap;
  }

  handlePointerMove(ev) {
    if (!this.dragState.active) return;
    const rect = this.canvas.getBoundingClientRect();
    if (rect.width === 0 || rect.height === 0) return;
    const scaleX = this.canvas.width / rect.width;
    const scaleY = this.canvas.height / rect.height;
    this.dragState.currentX = (ev.clientX - rect.left) * scaleX;
    this.dragState.currentY = (ev.clientY - rect.top) * scaleY;
  }

  handlePointerUp(ev) {
    if (this.deleteInterval) {
      clearInterval(this.deleteInterval);
      this.deleteInterval = null;
    }
    if (!this.dragState.active) return;
    const rect = this.canvas.getBoundingClientRect();
    const scaleX = rect.width === 0 ? 1 : this.canvas.width / rect.width;
    const scaleY = rect.height === 0 ? 1 : this.canvas.height / rect.height;
    const releaseX = ev
      ? (ev.clientX - rect.left) * scaleX
      : this.dragState.startX;
    const releaseY = ev
      ? (ev.clientY - rect.top) * scaleY
      : this.dragState.startY;
    const dx = releaseX - this.dragState.startX;
    const dy = releaseY - this.dragState.startY;
    const speedScale = this.system.constants.spawn.k_drag;
    const cap = this.system.constants.caps.v_spawn_cap;
    const speed = Math.min(Math.hypot(dx, dy) * speedScale, cap);
    const angle = Math.atan2(dy, dx) + Math.PI;
    const vx = Math.cos(angle) * speed;
    const vy = Math.sin(angle) * speed;

    this.spawnParticle(this.dragState.startX, this.dragState.startY, vx, vy);
    this.dragState.active = false;
  }

  spawnParticle(x, y, vx, vy) {
    const bounds = this.system.bounds;
    if (x < 0 || y < 0 || x > bounds.width || y > bounds.height) return;

    let type;
    switch (this.currentTool) {
      case Tool.PROTON:
        type = ParticleType.PROTON;
        break;
      case Tool.NEUTRON:
        type = ParticleType.NEUTRON;
        break;
      case Tool.ELECTRON:
      case Tool.ORBITING_ELECTRON:
        type = ParticleType.ELECTRON;
        break;
      default:
        return;
    }

    const index = this.system.addParticle(type, x, y, vx, vy);
    if (this.currentTool === Tool.ORBITING_ELECTRON) {
      const assist = this.orbitManager.computeAssist(index);
      if (assist) {
        this.system.vx[index] = assist.targetVx * this.system.orbitChirality;
        this.system.vy[index] = assist.targetVy * this.system.orbitChirality;
        const cap = this.system.constants.caps.v_spawn_cap;
        const mag = Math.hypot(this.system.vx[index], this.system.vy[index]);
        if (mag > cap) {
          const s = cap / mag;
          this.system.vx[index] *= s;
          this.system.vy[index] *= s;
        }
      }
    }
  }

  deleteAt(x, y) {
    const idx = this.system.findNearest(x, y, 10);
    if (idx >= 0) {
      this.flashes.spawn(this.system.x[idx], this.system.y[idx]);
      this.system.removeParticle(idx);
    }
  }

  stepSimulation() {
    const dt = this.system.constants.integrator.dt;
    this.orbitManager.updateClusters();
    this.system.integrate(dt, this.orbitManager, this.flashes);
    this.flashes.prune();
  }

  updateHUD(dt, fps) {
    this.fpsLabel.textContent = `FPS: ${fps.toFixed(1)}`;
    this.countLabel.textContent = `Particles: ${this.system.count}`;
  }

  tick() {
    const now = performance.now();
    const delta = (now - this.lastTimestamp) / 1000;
    this.lastTimestamp = now;

    if (!this.paused) {
      this.accumulator += delta;
      const dt = this.system.constants.integrator.dt;
      let steps = 0;
      const maxSteps = 6;
      while (this.accumulator >= dt && steps < maxSteps) {
        this.stepSimulation();
        this.accumulator -= dt;
        steps++;
      }
      if (steps === maxSteps) {
        this.accumulator = 0;
      }
    }

    this.renderer.render();
    const fps = 1 / Math.max(delta, 1 / 120);
    this.updateHUD(delta, fps);
    requestAnimationFrame(() => this.tick());
  }
}

function main() {
  const canvas = document.getElementById("sim-canvas");
  const bounds = { width: canvas.width, height: canvas.height };
  const system = new ParticleSystem(DEFAULT_CONSTANTS, bounds);
  const syncCanvasSize = () => {
    const displayWidth = Math.floor(canvas.clientWidth || canvas.width);
    const displayHeight = Math.floor(canvas.clientHeight || canvas.height);
    if (!displayWidth || !displayHeight) return;
    if (canvas.width !== displayWidth || canvas.height !== displayHeight) {
      canvas.width = displayWidth;
      canvas.height = displayHeight;
    }
    system.bounds.width = canvas.width;
    system.bounds.height = canvas.height;
  };
  syncCanvasSize();
  const flashes = new FlashManager();
  const orbitManager = new OrbitManager(system);
  const renderer = new Renderer(canvas, system, flashes, orbitManager);
  const ui = new UIController(system, renderer, orbitManager, flashes);
  window.addEventListener("resize", syncCanvasSize);
  requestAnimationFrame(() => ui.tick());
}

window.addEventListener("load", main);
