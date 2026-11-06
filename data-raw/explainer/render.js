// Save each panel of explainer.html as vignettes/articles/figures/how-<id>.png (2x). Needs Playwright with Chromium.
const { chromium } = require('playwright');
const path = require('path');
(async () => {
  const out = path.join(__dirname, '..', '..', 'vignettes', 'articles', 'figures');
  const b = await chromium.launch();
  const p = await b.newPage({ viewport: { width: 900, height: 600 }, deviceScaleFactor: 2 });
  await p.emulateMedia({ colorScheme: 'light' });
  await p.goto('file://' + path.join(__dirname, 'explainer.html'));
  await p.evaluate(() => document.fonts.ready); await p.waitForTimeout(300);
  for (const id of await p.$$eval('figure.panel', fs => fs.map(f => f.id))) {
    await p.locator('#' + id).screenshot({ path: path.join(out, 'how-' + id + '.png') });
  }
  await b.close();
})();
