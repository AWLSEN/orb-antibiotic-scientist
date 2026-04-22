# deploy/vps — Antibiotic Scientist API on orbcloud.dev

Deploys the same `state.py` that Vercel runs (`dashboard/vercel/api/state.py`)
to the orbcloud.dev VPS as a systemd service on `127.0.0.1:8091`. Nginx
proxies `orbcloud.dev/use-cases/antibiotic-scientist/api/*` here.

One canonical state.py lives at `dashboard/vercel/api/state.py`. The
`serve.py` wrapper in this directory imports it and binds a plain stdlib
`ThreadingHTTPServer`. The deploy workflow ships both files to
`/opt/antibiotic-api/` on the VPS.

## Files

| File | Goes to | What |
|---|---|---|
| `../../dashboard/vercel/api/state.py` | `/opt/antibiotic-api/state.py` | Request handler — shared with Vercel |
| `serve.py` | `/opt/antibiotic-api/serve.py` | stdlib `HTTPServer` wrapper that runs `handler` on :8091 |
| `antibiotic-api.service` | `/etc/systemd/system/antibiotic-api.service` | Hardened systemd unit |
| `nginx.conf.snippet` | paste into `/etc/nginx/sites-enabled/orbcloud` | Replaces the Vercel `proxy_pass` with a static + `/api/` split |

## One-time VPS setup

All commands run as root on `178.104.44.143`.

```bash
# 1. Dedicated user + dir
useradd --system --home /opt/antibiotic-api --shell /usr/sbin/nologin antibiotic-api
install -d -o antibiotic-api -g antibiotic-api -m 0755 /opt/antibiotic-api

# 2. Env file (chmod 600 — holds the org API key)
#    Grab ORB_API_KEY from either the current Vercel project env or mint
#    a new read-only key scoped to the antibiotic org.
cat > /etc/antibiotic-api.env <<'EOF'
ORB_API_KEY=<paste-key-here>
EOF
chmod 600 /etc/antibiotic-api.env
chown root:antibiotic-api /etc/antibiotic-api.env

# 3. First-time source install (subsequent pushes handled by CI below)
install -o antibiotic-api -g antibiotic-api -m 0644 \
  dashboard/vercel/api/state.py /opt/antibiotic-api/state.py
install -o antibiotic-api -g antibiotic-api -m 0755 \
  deploy/vps/serve.py /opt/antibiotic-api/serve.py

# 4. Install + start the service
install -m 0644 deploy/vps/antibiotic-api.service \
  /etc/systemd/system/antibiotic-api.service
systemctl daemon-reload
systemctl enable --now antibiotic-api
systemctl status antibiotic-api    # must show "active (running)"

# 5. Smoke test locally — expect JSON with a "timestamp" field
curl -s http://127.0.0.1:8091/api/state | head -c 200

# 6. Wire nginx. Edit /etc/nginx/sites-enabled/orbcloud — replace the
#    existing `location ^~ /use-cases/antibiotic-scientist { proxy_pass
#    <vercel>; }` block with the TWO blocks from nginx.conf.snippet.
vim /etc/nginx/sites-enabled/orbcloud     # or $EDITOR of choice
nginx -t
systemctl reload nginx

# 7. Smoke test the public path — expect the same JSON as step 5
curl -s https://orbcloud.dev/use-cases/antibiotic-scientist/api/state | head -c 200

# 8. Archive Vercel
#    Browser: https://vercel.com/<team>/orb-antibiotic-scientist → Settings →
#    pause/archive.
```

## Updates

After one-time setup, push-to-main on `state.py` or `deploy/vps/serve.py`
auto-deploys via `.github/workflows/deploy-vps.yml` — scp + systemd
restart, no manual step.

## Acceptance (cross-check with the migration spec)

- [ ] `systemctl status antibiotic-api` → `active (running)`, zero restart loop
- [ ] `curl -s http://127.0.0.1:8091/api/state` returns JSON locally
- [ ] `curl -sI https://orbcloud.dev/use-cases/antibiotic-scientist/` shows `Server: nginx`, not Vercel
- [ ] `curl -s https://orbcloud.dev/use-cases/antibiotic-scientist/api/state | jq '.computers[0].agent_state'` returns a real state
- [ ] Browser load shows orb·cloud nav + candidate cards within 3s
- [ ] Vercel project paused/archived
