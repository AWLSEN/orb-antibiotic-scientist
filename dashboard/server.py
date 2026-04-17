#!/usr/bin/env python3
"""
Tiny local HTTP server for the dashboard.

Serves dashboard/index.html at `/` and exposes the repository root as
additional static mount so the page can `fetch()` files under
`findings/` and `findings/loop-health/`. On Orb Cloud, you can point a
reverse proxy at port 8000 instead.

Usage:
  python dashboard/server.py [--port 8000] [--host 0.0.0.0]
"""

from __future__ import annotations

import argparse
import http.server
import os
import socketserver
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
INDEX = Path(__file__).parent / "index.html"


class Handler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(ROOT), **kwargs)

    def end_headers(self):
        # No-cache for live JSON; the page fetches with ?t=timestamp but
        # belt-and-braces.
        self.send_header("Cache-Control", "no-store")
        self.send_header("Access-Control-Allow-Origin", "*")
        super().end_headers()

    def do_GET(self):
        if self.path == "/" or self.path.startswith("/?"):
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            body = INDEX.read_bytes()
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        super().do_GET()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=int(os.environ.get("PORT", 8000)))
    ap.add_argument("--host", default=os.environ.get("HOST", "0.0.0.0"))
    args = ap.parse_args()

    with socketserver.TCPServer((args.host, args.port), Handler) as srv:
        srv.allow_reuse_address = True
        print(f"serving dashboard on http://{args.host}:{args.port}")
        try:
            srv.serve_forever()
        except KeyboardInterrupt:
            return 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
