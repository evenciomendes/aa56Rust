/**
 * Servidor Node.js para aa56-rust (DE404, domoon, conjunct, moonrise, star.cat, orbit.cat).
 */

const http = require('http');
const { execFile } = require('child_process');
const fs = require('fs');
const path = require('path');
const url = require('url');

const PORT = process.env.PORT || 3070;

function getRustBinaryPath() {
  const candidates = [
    path.join(__dirname, '..', 'aa56-rust', 'target', 'release', 'aa56-cli.exe'),
    path.join(__dirname, '..', 'aa56-rust', 'target', 'release', 'aa56-cli'),
    path.join(__dirname, '..', 'aa56-rust', 'target', 'debug', 'aa56-cli.exe'),
    path.join(__dirname, '..', 'aa56-rust', 'target', 'debug', 'aa56-cli'),
  ];
  for (const bin of candidates) {
    if (fs.existsSync(bin)) return bin;
  }
  return null;
}

const { computeEphemerisFallback } = require('./aa56_fallback');

function callEngine(params, modeOverride) {
  return new Promise((resolve, reject) => {
    const mode = modeOverride || params.mode || 'ephemeris';
    const binPath = getRustBinaryPath();

    if (!binPath) {
      try {
        const data = computeEphemerisFallback(params);
        data.engine_mode = 'Espelho JS de Homologação (Código Rust pronto em aa56-rust)';
        return resolve(data);
      } catch (e) {
        return reject(e);
      }
    }

    const args = ['--mode', mode];
    if (params.year) args.push('--year', String(params.year));
    if (params.month) args.push('--month', String(params.month));
    if (params.day) args.push('--day', String(params.day));
    if (params.hour) args.push('--hour', String(params.hour));
    if (params.minute) args.push('--minute', String(params.minute));
    if (params.second) args.push('--second', String(params.second));
    if (params.dt) args.push('--dt', String(params.dt));
    if (params.lon) args.push('--lon', String(params.lon));
    if (params.lat) args.push('--lat', String(params.lat));
    if (params.days) args.push('--days', String(params.days));
    if (params.cycles) args.push('--cycles', String(params.cycles));

    const t0 = process.hrtime.bigint();
    execFile(binPath, args, { timeout: 15000, maxBuffer: 15 * 1024 * 1024 }, (err, stdout, stderr) => {
      const elapsedMs = Number(process.hrtime.bigint() - t0) / 1e6;
      if (err) {
        return reject(new Error(`Erro ao executar aa56-cli: ${err.message} ${stderr || ''}`));
      }
      try {
        const sanitized = stdout.replace(/(\d)"",/g, '$1\\"",');
        const data = JSON.parse(sanitized);
        data.node_roundtrip_ms = Number(elapsedMs.toFixed(3));
        data.engine_mode = 'Binário Nativo Rust DE404 (aa56-cli)';
        resolve(data);
      } catch (parseErr) {
        reject(new Error(`Falha ao ler JSON do Rust: ${parseErr.message}\nSaída: ${stdout}`));
      }
    });
  });
}

const server = http.createServer(async (req, res) => {
  const parsedUrl = url.parse(req.url, true);
  const pathname = parsedUrl.pathname;

  if (pathname === '/healthz') {
    res.setHeader('Content-Type', 'application/json; charset=utf-8');
    res.writeHead(200);
    res.end(JSON.stringify({ status: 'ok', service: 'aa56-rust' }));
    return;
  }

  const apiModes = {
    '/api/ephemeris': 'ephemeris',
    '/api/domoon': 'ephemeris',
    '/api/conjunct': 'conjunct',
    '/api/moonrise': 'moonrise',
    '/api/stars': 'stars',
    '/api/stars1950': 'stars1950',
    '/api/messier': 'messier',
    '/api/orbits': 'orbits',
  };

  if (apiModes[pathname]) {
    res.setHeader('Content-Type', 'application/json; charset=utf-8');
    try {
      const result = await callEngine(parsedUrl.query, parsedUrl.query.mode || apiModes[pathname]);
      res.writeHead(200);
      res.end(JSON.stringify(result, null, 2));
    } catch (err) {
      res.writeHead(500);
      res.end(JSON.stringify({ error: err.message }));
    }
    return;
  }

  // Frontend SPA
  const htmlPath = path.join(__dirname, 'public', 'index.html');
  fs.readFile(htmlPath, 'utf8', (err, content) => {
    if (err) {
      res.writeHead(500);
      res.end('Erro ao carregar public/index.html');
      return;
    }
    res.setHeader('Content-Type', 'text/html; charset=utf-8');
    res.writeHead(200);
    res.end(content);
  });
});

server.listen(PORT, '0.0.0.0', () => {
  console.log('==============================================================');
  console.log(`🚀 Painel de Testes AA-56 (DE404 + domoon + Node.js) rodando em:`);
  console.log(`👉 http://localhost:${PORT}`);
  console.log('==============================================================');
});
