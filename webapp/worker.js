/* Runs the engine off the page's main thread, so the sliders stay responsive while a case solves. */
importScripts("engine.js", "presets.js", "checks.js", "jobs.js");
const JOBS = self.EJC_JOBS.make(self.EJC, self.EJC_PRESETS, self.EJC_CHECKS.make(self.EJC, self.EJC_PRESETS));
self.onmessage = (ev) => {
  const job = ev.data;
  if (job.type === "ping") { self.postMessage({ type: "pong" }); return; }
  let out;
  try {
    out = JOBS.handle(job, (p) => self.postMessage(Object.assign({ type: "progress", id: job.id, job: job.type }, p)));
  } catch (e) {
    out = { type: "error", id: job.id, job: job.type, error: String((e && e.message) || e) };
  }
  self.postMessage(out);
};
