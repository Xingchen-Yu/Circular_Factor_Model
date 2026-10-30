// Full Circle worker: serves the static site (website/) and one tiny API, a counter of
// visitors asking for a Senate version. Stores only a number: no IPs, no cookies.
//   GET  /api/senate -> {"count": n}
//   POST /api/senate -> adds one, returns {"count": n}
export default {
  async fetch(request, env) {
    const url = new URL(request.url);
    if (url.pathname !== "/api/senate") return env.ASSETS.fetch(request);
    await env.DB.prepare("CREATE TABLE IF NOT EXISTS counters (name TEXT PRIMARY KEY, n INTEGER NOT NULL)").run();
    if (request.method === "POST") {
      await env.DB.prepare("INSERT INTO counters (name, n) VALUES ('senate', 1) ON CONFLICT(name) DO UPDATE SET n = n + 1").run();
    } else if (request.method !== "GET") {
      return new Response("Method not allowed", { status: 405 });
    }
    const row = await env.DB.prepare("SELECT n FROM counters WHERE name = 'senate'").first();
    return Response.json({ count: row ? row.n : 0 }, { headers: { "cache-control": "no-store" } });
  },
};
