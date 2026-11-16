// Black Box recorder (Group 5). Records runs of BlackBox_N.php from inside the page itself.
//
// How to use:
//   1. Open https://casci.binghamton.edu/academics/ssie501/blackbox/BlackBox_N.php in a browser.
//   2. Open the browser's developer console on that page and paste this whole file.
//   3. Run, for example:   await REC.batch([{presses: 2000, n: 1, label: "single_01"}])
//      Each finished run is saved by the browser as blackbox_<label>_<date>.json (normally in Downloads).
//   4. Move the files to data/group5_runs/ and run analysis/12_fresh_runs.py.
//
// What it does: it requests the same URLs as the page's own buttons, one after another, on the page's own
// session: "?reset=1&cycles_input=1" (Reset) and "?cycles=<n>" (Next n Step). After every request it reads
// the 400 cells and the "Current step:" counter from the returned page. Every request is checked: the counter
// must rise by exactly n, and the page must show 400 cells; anything else is written to run.errors.
// A run file holds the Reset grid, then for every press only the cells that changed ([cell index, new value],
// cells numbered 0-399 row by row), the step counter after every press, and the final grid as a check.

window.REC = { runs: [], status: "idle", saved: new Set() };

REC.page = async (query) => {
  const r = await fetch("BlackBox_N.php" + query, { cache: "no-store", credentials: "same-origin" });
  const d = new DOMParser().parseFromString(await r.text(), "text/html");
  const v = [...d.querySelectorAll("td")].map((t) => {
    const s = t.textContent.trim();
    return s === "" ? -1 : Number(s);
  });
  const m = d.body.innerText.match(/Current step:\s*(\d+)/);
  return { v, step: m ? Number(m[1]) : null };
};

REC.run = async (presses, n, label) => {
  const g0 = await REC.page("?reset=1&cycles_input=1");
  if (g0.v.length !== 400) throw new Error("reset grid has " + g0.v.length + " cells");
  const run = { label, cycles: n, started: new Date().toISOString(), reset_step: g0.step,
                reset: g0.v, changes: [], steps: [], errors: [] };
  let prev = g0;
  for (let k = 0; k < presses; k++) {
    const g = await REC.page("?cycles=" + n);
    if (g.v.length !== 400) run.errors.push([k, "cells", g.v.length]);
    if (prev.step !== null && g.step !== prev.step + n) run.errors.push([k, "step", prev.step, g.step]);
    const ch = [];
    for (let i = 0; i < 400; i++) if (g.v[i] !== prev.v[i]) ch.push([i, g.v[i]]);
    run.changes.push(ch);
    run.steps.push(g.step);
    prev = g;
    REC.progress = label + ": " + (k + 1) + "/" + presses;
  }
  run.final = prev.v;
  run.ended = new Date().toISOString();
  REC.runs.push(run);
  REC.save(run);
  return run;
};

REC.save = (run) => {
  const blob = new Blob([JSON.stringify(run)], { type: "application/json" });
  const a = document.createElement("a");
  const day = run.started.slice(0, 10).replaceAll("-", "");
  a.href = URL.createObjectURL(blob);
  a.download = "blackbox_" + run.label + "_" + day + ".json";
  document.body.appendChild(a);
  a.click();
  a.remove();
  REC.saved.add(run.label);
  return { file: a.download, bytes: blob.size };
};

REC.batch = async (specs) => {
  REC.status = "running";
  try {
    for (const s of specs) await REC.run(s.presses, s.n, s.label);
    REC.status = "done";
  } catch (e) {
    REC.status = "error: " + e.message;
  }
};
