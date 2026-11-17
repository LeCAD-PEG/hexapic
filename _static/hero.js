// Copyright 2026 David Lajevec — HEXAPIC
//
// The home page's background: a real 1D electrostatic PIC run of the
// two-stream instability, drawn in phase space (x across, v up).  Two
// counter-streaming electron beams over a fixed ion background, cloud-in-cell
// deposit, the field from Gauss's law on a periodic grid, leapfrog push; the
// units are the plasma frequency and the Debye length.  The beams' sinusoidal
// perturbation grows until the phase space rolls up into vortices, and then
// the run starts over.  (The previous website's hero, in this site's colours:
// BEAM_A, BEAM_B.)
//
// ?hero=<t> starts at time t, without the entrance (a still for a screenshot).
(function () {
  "use strict";

  var L = 50, NX = 256, DX = L / NX;   // the box and its grid
  var N = 14000;                        // electrons, half in each beam
  var VB = 3, VTH = 1, AMP = 0.1;       // the beam speed, thermal speed, perturbation
  var DT = 0.03, T_END = 95;            // per frame (a quarter of the first 0.12: a slow roll-up); one run's length
  var VMAX = 8.5;                       // the plot's velocity range
  var BEAM_A = "#f2f2f2", BEAM_B = "#9085e9";   // the beams: the ink and the accent

  function gauss(rand) {               // Box-Muller
    var u = 1 - rand(), v = rand();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  function mulberry32(a) {             // a small seeded generator: the same run each visit
    return function () {
      a |= 0; a = (a + 0x6d2b79f5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  function Run(seed) {
    var rand = mulberry32(seed);
    this.x = new Float32Array(N);
    this.v = new Float32Array(N);
    this.n = new Float32Array(NX);
    this.E = new Float32Array(NX);
    for (var i = 0; i < N; i++) {
      var x = rand() * L;
      var v = VTH * gauss(rand) + (i < N / 2 ? VB : -VB);
      this.x[i] = x;
      this.v[i] = v * (1 + AMP * Math.sin((2 * Math.PI * x) / L));
    }
    this.t = 0;
    this.field();
    for (i = 0; i < N; i++) this.v[i] -= 0.5 * DT * -this.at(this.x[i]);   // v at t - dt/2
  }

  Run.prototype.field = function () {
    var n = this.n, E = this.E, i, j, f, w = NX / N;   // the density in units of the mean
    n.fill(0);
    for (i = 0; i < N; i++) {
      f = this.x[i] / DX; j = Math.floor(f); f -= j;
      n[j % NX] += (1 - f) * w; n[(j + 1) % NX] += f * w;
    }
    // dE/dx = 1 - n (ions minus electrons), E with zero mean on the periodic box
    var acc = 0, mean = 0;
    for (j = 0; j < NX; j++) { acc += (1 - n[j]) * DX; E[j] = acc; mean += acc; }
    mean /= NX;
    for (j = 0; j < NX; j++) E[j] -= mean;
  };

  Run.prototype.at = function (x) {
    var f = x / DX, j = Math.floor(f); f -= j;
    return (1 - f) * this.E[j % NX] + f * this.E[(j + 1) % NX];
  };

  Run.prototype.step = function () {
    for (var i = 0; i < N; i++) {
      this.v[i] += DT * -this.at(this.x[i]);          // electrons: a = -E
      var x = this.x[i] + DT * this.v[i];
      x %= L; if (x < 0) x += L;
      this.x[i] = x;
    }
    this.field();
    this.t += DT;
  };

  function start() {
    var hero = document.querySelector(".hx-hero");
    var canvas = hero && hero.querySelector(".hx-hero__canvas");
    if (!canvas) return;
    var ctx = canvas.getContext("2d");
    var reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    var skip = parseFloat(new URLSearchParams(window.location.search).get("hero")) || 0;
    var run = new Run(7), seed = 7, raf = 0, visible = true, fade = 1;
    var w = 0, h = 0, dpr = 1;

    function resize() {
      dpr = Math.min(window.devicePixelRatio || 1, 2);
      w = canvas.clientWidth; h = canvas.clientHeight;
      canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }

    function draw() {
      ctx.clearRect(0, 0, w, h);
      var sx = w / L, sy = h / (2 * VMAX), half = N / 2, s = 1.6;
      ctx.globalAlpha = 0.85 * fade;
      ctx.fillStyle = BEAM_A;
      for (var i = 0; i < half; i++)
        ctx.fillRect(run.x[i] * sx, h / 2 - run.v[i] * sy, s, s);
      ctx.fillStyle = BEAM_B;
      for (i = half; i < N; i++)
        ctx.fillRect(run.x[i] * sx, h / 2 - run.v[i] * sy, s, s);
      ctx.globalAlpha = 1;
    }

    function frame() {
      raf = 0;
      if (!visible || document.hidden) return;
      run.step();
      if (run.t > T_END) {                             // fade out, start a new run
        fade -= 0.02;
        if (fade <= 0) { seed += 1; run = new Run(seed); fade = 0; }
      } else if (fade < 1) fade = Math.min(1, fade + 0.02);
      draw();
      raf = requestAnimationFrame(frame);
    }

    function resume() { if (!raf && !reduced) raf = requestAnimationFrame(frame); }

    resize();
    if (skip > 0) { hero.classList.add("hx-hero--still"); while (run.t < skip) run.step(); }
    if (reduced) {                                      // one still frame, the vortices formed
      while (run.t < 55) run.step();
      draw();
    } else {
      draw();
      resume();
    }

    window.addEventListener("resize", function () { resize(); draw(); });
    document.addEventListener("visibilitychange", function () { if (!document.hidden) resume(); });
    new IntersectionObserver(function (es) {
      visible = es[0].isIntersecting; if (visible) resume();
    }).observe(canvas);
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start);
  else start();
})();
