"""
systemd entrypoint for antibiotic-api on orbcloud.dev.

Wraps the `handler` class from state.py (shared with the original Vercel
function) in a plain stdlib HTTPServer bound to 127.0.0.1:8091. nginx on
the VPS proxies /use-cases/antibiotic-scientist/api/* here after stripping
the prefix, so handler sees paths like /api/state exactly as it did on
Vercel. Keeping state.py byte-identical to the upstream makes future
sync trivial.

Env:
  ORB_API_KEY    required — bearer token for api.orbcloud.dev
  ORB_BASE_URL   optional — default https://api.orbcloud.dev
  BIND_HOST      optional — default 127.0.0.1 (never expose publicly)
  BIND_PORT      optional — default 8091
"""

from __future__ import annotations

import os
import sys
from http.server import HTTPServer, ThreadingHTTPServer

from state import handler

HOST = os.environ.get("BIND_HOST", "127.0.0.1")
PORT = int(os.environ.get("BIND_PORT", "8091"))


def main() -> int:
    if not os.environ.get("ORB_API_KEY", "").strip():
        print("antibiotic-api: ORB_API_KEY not set — refusing to start", file=sys.stderr)
        return 1
    server = ThreadingHTTPServer((HOST, PORT), handler)
    print(f"antibiotic-api listening on http://{HOST}:{PORT}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
