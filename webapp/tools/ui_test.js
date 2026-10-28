// Drive the page in Chromium: no console errors, the default case solves to the thesis values, every tab draws,
// inputs and presets update the results, and the layout holds at phone width and in dark mode. Screenshots go to
// tools/shots/.
// Usage (from webapp/): python3 tools/make_preview.py && python3 -m http.server 8765 --bind 127.0.0.1 &
//                        node tools/ui_test.js
// Offline (no access to the CDN): PLOTLY_JS=/path/to/plotly-2.35.2.min.js node tools/ui_test.js serves that file
// for the CDN URL and stubs the Google Fonts requests. Needs the playwright package and its Chromium.
"use strict";
const fs = require("fs");
const path = require("path");
const { chromium } = require("playwright");

const URL = process.env.URL || "http://127.0.0.1:8765/preview.html";
const SHOTS = path.join(__dirname, "shots");
fs.mkdirSync(SHOTS, { recursive: true });
const PLOTLY = process.env.PLOTLY_JS ? fs.readFileSync(process.env.PLOTLY_JS) : null;

async function newPage(browser, opts) {
  const ctx = await browser.newContext(opts);
  if (PLOTLY) {
    await ctx.route("https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.35.2/plotly.min.js",
      (r) => r.fulfill({ status: 200, contentType: "application/javascript", body: PLOTLY }));
    await ctx.route("https://fonts.googleapis.com/**", (r) => r.fulfill({ status: 200, contentType: "text/css", body: "" }));
    await ctx.route("https://fonts.gstatic.com/**", (r) => r.abort());
  }
  const page = await ctx.newPage();
  const problems = [];
  page.on("console", (m) => { if (m.type() === "error" || m.type() === "warning") problems.push(`${m.type()}: ${m.text()}`); });
  page.on("pageerror", (e) => problems.push(`pageerror: ${e.message}`));
  return { ctx, page, problems };
}
const text = (page, sel) => page.$eval(sel, (el) => el.textContent.trim());
async function waitSolved(page) {
  await page.waitForFunction(() => /Up to date/.test(document.querySelector("#runstate").textContent), null, { timeout: 120000 });
}
async function check(cond, label) { console.log(`${cond ? "ok  " : "FAIL"} ${label}`); if (!cond) process.exitCode = 1; }

(async () => {
  const browser = await chromium.launch();
  // ---------------- desktop, light
  let { ctx, page, problems } = await newPage(browser, { viewport: { width: 1440, height: 1000 }, colorScheme: "light" });
  await page.goto(URL);
  await waitSolved(page);
  await page.waitForSelector("#chart-3d .main-svg", { timeout: 60000 });
  await check(await text(page, "#fos3") === "2.706024", `default 3D FoS ${await text(page, "#fos3")}`);
  await check(await text(page, "#fos2") === "2.447911", `default 2D FoS ${await text(page, "#fos2")}`);
  await check(/Accepted/i.test(await text(page, "#chip3")), "3D accepted chip");
  await check(/6 of 76 \(8%\) \/ 1 of 30 \(3%\)/.test(await text(page, "#facts")), "summary: tensile bases 6 of 76 (3D) and 1 of 30 (2D)");
  await check(/Tensile bases: 3D 6 of 76 \(8%\), 2D 1 of 30 \(3%\)/.test(await text(page, "#notes")), "summary: tensile-base note (info)");
  await check(!/Vertical/.test(await text(page, "#notes")), "summary: no vertical-face note for the main sphere");
  console.log("     runstate:", await text(page, "#runstate"), "| worker:", await page.evaluate(() => document.querySelector("#runstate").textContent.includes("on this page") ? "fallback" : "worker"));
  await page.waitForTimeout(600);
  await page.screenshot({ path: path.join(SHOTS, "w1-3d-light.png"), fullPage: true });

  for (const tab of ["cross-section", "free-body", "checks"]) {
    await page.click(`#tab-${tab}`);
    await page.waitForTimeout(900);
    await page.screenshot({ path: path.join(SHOTS, `w1-${tab}-light.png`), fullPage: true });
  }
  await check((await page.$$eval("#panel-checks td .res.fail", (e) => e.length)) === 0, "no failed final-state checks for the default case");
  await page.click("#tab-free-body");
  await page.waitForTimeout(400);
  console.log("     column:", await text(page, "#col-title"), "|", (await text(page, "#col-balance")).slice(0, 120));
  console.log("     slice:", await text(page, "#slice-title"), "|", (await text(page, "#slice-balance")).slice(0, 120));
  await check(/PASS/.test(await text(page, "#col-balance")) && /PASS/.test(await text(page, "#slice-balance")), "vertical balance PASS for column and slice");
  // every force on the sliding mass: 2D section, then the 3D full mass
  await page.waitForSelector("#chart-mass .main-svg", { timeout: 60000 });
  let mb = await text(page, "#mass-balance");
  console.log("     mass 2D:", mb.replace(/\s+/g, " ").slice(0, 260));
  await check((mb.match(/PASS/g) || []).length === 2 && /1\.89% of the weight/.test(mb), "whole mass 2D: vertical and moment balance PASS, 1.89% horizontal");
  await check((await page.$$eval("#mass-table tbody tr", (e) => e.length)) === 5, "whole mass 2D: table of W, N′, U, S and their sum");
  await page.check("#massview-3");
  await page.waitForTimeout(1500);
  mb = await text(page, "#mass-balance");
  await check(/FH balances ΣS = ΣWκ instead/.test(mb) && /1\.50% of the weight/.test(mb) && (mb.match(/PASS/g) || []).length === 1, "whole mass 3D: F_H equation, 1.50% horizontal");
  await page.screenshot({ path: path.join(SHOTS, "w1-free-body-3d-light.png"), fullPage: false });
  await page.check("#massview-2");
  await page.waitForTimeout(600);
  // change the column with its slider
  await page.$eval("#sl-k3", (el) => { el.value = 5; el.dispatchEvent(new Event("input", { bubbles: true })); });
  await page.waitForTimeout(500);
  await check(/Column 5 /.test(await text(page, "#col-title")), `column slider selects column 5: ${await text(page, "#col-title")}`);

  // verification tab runs by itself
  await page.click("#tab-verification");
  await page.waitForFunction(() => /Finished/.test(document.querySelector("#ver-status").textContent), null, { timeout: 120000 });
  console.log("     verification:", await text(page, "#ver-summary"), "|", await text(page, "#ver-status"));
  await check(/82 pass/i.test(await text(page, "#ver-summary")) && /0 fail/i.test(await text(page, "#ver-summary")), "verification 82 pass / 0 fail");
  await page.screenshot({ path: path.join(SHOTS, "w1-verification-light.png"), fullPage: false });

  // studies: a short nu sweep and mesh study
  await page.click("#tab-studies");
  await page.fill("#nu0", "0.2"); await page.fill("#nu1", "0.3"); await page.fill("#dnu", "0.05");
  await page.click("#nu-run");
  await page.waitForSelector("#chart-nu .main-svg", { timeout: 60000 });
  for (const n of [40, 80]) await page.uncheck(`#m3-${n}`);
  for (const n of [80, 160, 320]) await page.uncheck(`#m2-${n}`);
  await page.click("#mesh-run");
  await page.waitForSelector("#chart-mesh3 .main-svg", { timeout: 60000 });
  await page.waitForTimeout(500);
  console.log("     mesh study:", (await text(page, "#mesh-out")).slice(0, 220).replace(/\s+/g, " "));
  // critical-surface search (2D), then make the critical circle the prescribed surface
  await page.click("#srch-run");
  await page.waitForSelector("#chart-search .main-svg", { timeout: 180000 });
  const srch = await text(page, "#srch-out");
  console.log("     search:", (await text(page, "#srch-status")), "|", srch.slice(0, 260).replace(/\s+/g, " "));
  await check(/FoS on the critical surface\s*2\.118587/.test(srch) && /−13\.5 %/.test(srch) && /160 slices:\s*2\.125346/.test(srch),
    "search: 2D critical FoS 2.118587 (−13.5 %), mesh check 2.125346 on 160 slices");
  await page.waitForTimeout(500);
  await page.screenshot({ path: path.join(SHOTS, "w1-studies-light.png"), fullPage: true });
  await page.click('[data-use="2"]');
  await page.waitForFunction(() => /Up to date/.test(document.querySelector("#runstate").textContent) && document.querySelector("#fos2").textContent === "2.118587", null, { timeout: 60000 });
  await check(await text(page, "#fos2") === "2.118587" && await page.$eval("#nb-R", (el) => el.value) === "592.030235019", "Use this circle: the analysis reruns on the critical circle");
  await page.click("#reset");
  await waitSolved(page);

  // number box: R = 760 m -> new FoS; then back to the preset
  await page.click("#tab-view-3d");
  await page.fill("#nb-R", "760");
  await page.dispatchEvent("#nb-R", "change");
  await page.waitForFunction(() => /Up to date/.test(document.querySelector("#runstate").textContent) && document.querySelector("#fos3").textContent !== "2.706024", null, { timeout: 60000 });
  console.log("     R = 760:", await text(page, "#fos3"), await text(page, "#fos2"), "|", (await text(page, "#notes")).slice(0, 160));
  // slider: beta
  await page.$eval("#sl-beta", (el) => { el.value = 45; el.dispatchEvent(new Event("input", { bubbles: true })); });
  await waitSolved(page);
  console.log("     beta = 45:", await text(page, "#fos3"), await text(page, "#fos2"), "| preset state:", await text(page, "#preset-state"));
  await page.click("#reset");
  await waitSolved(page);
  await check(await text(page, "#fos3") === "2.706024", "reset to preset restores 2.706024");

  // presets: native v15 wet (63 iterations), D'Esterre with floor columns, original closure, a slow preset
  await page.selectOption("#preset", "thesis_native_wet");
  await waitSolved(page);
  const facts = await text(page, "#facts");
  await check(await text(page, "#fos3") === "1.272331" && /63 \//.test(facts), `native wet: ${await text(page, "#fos3")} | ${facts.slice(0, 60)}`);
  await page.selectOption("#preset", "desterre_DE-1.5-c3-ru0.5");
  await waitSolved(page);
  await check((await text(page, "#fos3")).startsWith("4.706") && (await text(page, "#fos2")).startsWith("3.84"), `D'Esterre 1.5-c3-ru0.5: ${await text(page, "#fos3")} / ${await text(page, "#fos2")}`);
  await page.waitForTimeout(600);
  await page.screenshot({ path: path.join(SHOTS, "w1-desterre-3d-light.png"), fullPage: false });
  await page.click("#tab-cross-section");
  await page.waitForTimeout(700);
  await page.screenshot({ path: path.join(SHOTS, "w1-desterre-section-light.png"), fullPage: false });
  await page.selectOption("#preset", "thesis_orig_0271");
  await waitSolved(page);
  console.log("     original closure:", (await text(page, "#orig")).slice(0, 200).replace(/\s+/g, " "));
  await page.screenshot({ path: path.join(SHOTS, "w1-original-light.png"), fullPage: false });
  await page.selectOption("#preset", "thesis_finest_dry");
  await page.waitForTimeout(300);
  console.log("     slow preset:", await text(page, "#runstate"));
  // an input error: a circle behind the toe with the default grid
  await page.selectOption("#preset", "thesis_dry_05");
  await waitSolved(page);
  await page.fill("#nb-R", "700"); await page.dispatchEvent("#nb-R", "change");
  await waitSolved(page);
  console.log("     R = 700:", (await text(page, "#notes")).slice(0, 220));
  await page.fill("#nb-R", "600"); await page.dispatchEvent("#nb-R", "change");
  await waitSolved(page);
  await check(/require R > z_c .*use the base-circle option/.test(await text(page, "#notes")), "R < z_c: the error names the base-circle option");
  // a circle whose arc ends below the ground at its back extreme: the summary warns about the vertical face
  for (const [k, v] of [["xc", "25"], ["zc", "250"], ["R", "251"]]) { await page.fill(`#nb-${k}`, v); await page.dispatchEvent(`#nb-${k}`, "change"); }
  await waitSolved(page);
  await check(/Vertical back face\..*19\.34 m high.*\(269\.34 m\)/.test(await text(page, "#notes")) && await text(page, "#fos3") === "2.758886",
    "vertical back face warning (19.34 m), FoS unchanged");
  await page.screenshot({ path: path.join(SHOTS, "w1-vertical-face-light.png"), fullPage: false });
  await page.selectOption("#preset", "thesis_wet_05");
  await waitSolved(page);
  await check(/Tensile bases: 3D 23 of 76 \(30%\).*More than a quarter/.test(await text(page, "#notes")), "wet baseline: tensile-base warning (23 of 76)");
  await page.selectOption("#preset", "thesis_dry_05");
  await waitSolved(page);
  await page.fill("#nb-R", "700"); await page.dispatchEvent("#nb-R", "change");
  await waitSolved(page);
  await page.click("#tab-cross-section");
  await page.waitForTimeout(600);
  await page.screenshot({ path: path.join(SHOTS, "w1-error-section-light.png"), fullPage: false });
  // an out-of-range value
  await page.fill("#nb-gamma", "55"); await page.dispatchEvent("#nb-gamma", "change");
  await page.waitForTimeout(300);
  console.log("     gamma = 55:", await text(page, "#runstate"), "|", await text(page, "#err-gamma"));
  await page.fill("#nb-gamma", "27"); await page.dispatchEvent("#nb-gamma", "change");
  console.log("     problems (desktop light):", problems.length ? problems : "none");
  await ctx.close();

  // ---------------- desktop, dark (a fresh browser context, so no stored inputs)
  ({ ctx, page, problems } = await newPage(browser, { viewport: { width: 1440, height: 1000 }, colorScheme: "dark" }));
  await page.goto(URL);
  await waitSolved(page);
  await page.waitForSelector("#chart-3d .main-svg", { timeout: 60000 });
  await page.waitForTimeout(800);
  await page.screenshot({ path: path.join(SHOTS, "w2-3d-dark.png"), fullPage: true });
  for (const tab of ["cross-section", "free-body"]) {
    await page.click(`#tab-${tab}`);
    await page.waitForTimeout(900);
    await page.screenshot({ path: path.join(SHOTS, `w2-${tab}-dark.png`), fullPage: true });
  }
  console.log("     problems (desktop dark):", problems.length ? problems : "none");
  await ctx.close();

  // ---------------- phone
  ({ ctx, page, problems } = await newPage(browser, { viewport: { width: 390, height: 844 }, colorScheme: "light", deviceScaleFactor: 2, isMobile: true, hasTouch: true }));
  await page.goto(URL);
  await waitSolved(page);
  await page.waitForTimeout(800);
  const overflow = await page.evaluate(() => [document.documentElement.scrollWidth, document.documentElement.clientWidth]);
  await check(overflow[0] <= overflow[1], `no horizontal scroll at 390 px (${overflow})`);
  await page.screenshot({ path: path.join(SHOTS, "w3-phone-top.png"), fullPage: false });
  await page.screenshot({ path: path.join(SHOTS, "w3-phone-full.png"), fullPage: true });
  await page.evaluate(() => window.scrollTo(0, 1400));
  await page.waitForTimeout(500);
  await check(await page.$eval("#minibar", (el) => el.classList.contains("show")), "phone: FoS bar shows after scrolling past the summary");
  await page.screenshot({ path: path.join(SHOTS, "w3-phone-scrolled.png"), fullPage: false });
  await page.evaluate(() => window.scrollTo(0, 0));
  await page.waitForTimeout(400);
  await check(!(await page.$eval("#minibar", (el) => el.classList.contains("show"))), "phone: FoS bar hidden at the top");
  await page.click("#tab-free-body");
  await page.waitForTimeout(800);
  const overflow2 = await page.evaluate(() => [document.documentElement.scrollWidth, document.documentElement.clientWidth]);
  await check(overflow2[0] <= overflow2[1], `no horizontal scroll on the FBD tab at 390 px (${overflow2})`);
  await page.click("#tab-studies");
  await page.click("#srch-run");
  await page.waitForSelector("#chart-search .main-svg", { timeout: 180000 });
  await page.waitForTimeout(500);
  const overflow3 = await page.evaluate(() => [document.documentElement.scrollWidth, document.documentElement.clientWidth]);
  await check(overflow3[0] <= overflow3[1], `no horizontal scroll on the Studies tab with search results at 390 px (${overflow3})`);
  await page.screenshot({ path: path.join(SHOTS, "w3-phone-search.png"), fullPage: true });
  console.log("     problems (phone):", problems.length ? problems : "none");
  await ctx.close();
  await browser.close();
})().catch((e) => { console.error(e); process.exit(1); });
