#!/usr/bin/env bash
# ============================================================================
# orb_deploy.sh — reproducible end-to-end deploy of orb-antibiotic-scientist
# to Orb Cloud via the public REST API at https://api.orbcloud.dev/.
#
# The script is idempotent: every step checks whether the resource already
# exists and skips when it does, so you can safely re-run it.
#
# Required env vars:
#   ORB_API_KEY          — Bearer token from https://app.orbcloud.dev/
#                          OR leave empty to auto-register a new key.
#   ANTHROPIC_AUTH_TOKEN — Z.AI GLM Coding Plan key (default provider)
#                          OR ANTHROPIC_API_KEY if --provider anthropic.
#                          Forwarded into the agent via org_secrets.
#
# Optional env vars:
#   ORB_LLM_PROVIDER     — "zai" (default) or "anthropic". Controls which
#                          secret name is forwarded to the agent.
#   ORB_COMPUTER_NAME    — default "orb-antibiotic-scientist"
#   ORB_RUNTIME_MB       — default 4096 (match orb.toml [resources])
#   ORB_DISK_MB          — default 10240
#   ORB_ORB_TOML         — default ./orb.toml
#
# The script writes the computer-id and agent-id to .orb-state/ so later
# runs can `promote`, `demote`, or `logs` without guessing.
#
# Usage:
#   ORB_API_KEY=orb_... ANTHROPIC_API_KEY=sk-... ./scripts/orb_deploy.sh deploy
#   ./scripts/orb_deploy.sh status
#   ./scripts/orb_deploy.sh promote
#   ./scripts/orb_deploy.sh demote
#   ./scripts/orb_deploy.sh logs
# ============================================================================

set -euo pipefail

BASE_URL="${ORB_BASE_URL:-https://api.orbcloud.dev}"
COMPUTER_NAME="${ORB_COMPUTER_NAME:-orb-antibiotic-scientist}"
RUNTIME_MB="${ORB_RUNTIME_MB:-4096}"
DISK_MB="${ORB_DISK_MB:-10240}"
ORB_TOML_PATH="${ORB_ORB_TOML:-./orb.toml}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STATE_DIR="${REPO_ROOT}/.orb-state"
mkdir -p "$STATE_DIR"

COMPUTER_FILE="$STATE_DIR/computer-id"
AGENT_FILE="$STATE_DIR/agent-id"

log() { printf "[orb_deploy] %s\n" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

require() {
    command -v "$1" >/dev/null 2>&1 || die "missing required binary: $1"
}

require curl
require jq

# --------------------------------------------------------------------------
# Auth
# --------------------------------------------------------------------------

ensure_api_key() {
    if [[ -n "${ORB_API_KEY:-}" ]]; then
        return 0
    fi
    log "ORB_API_KEY unset — registering a new key at $BASE_URL/api/v1/auth/register"
    : "${ORB_REGISTER_EMAIL:?ORB_REGISTER_EMAIL required to auto-register}"
    local resp
    resp=$(curl -fsS -X POST "$BASE_URL/api/v1/auth/register" \
        -H 'Content-Type: application/json' \
        -d "{\"email\":\"$ORB_REGISTER_EMAIL\"}")
    ORB_API_KEY=$(printf '%s' "$resp" | jq -r .api_key)
    [[ -n "$ORB_API_KEY" && "$ORB_API_KEY" != "null" ]] || die "register did not return an api_key"
    printf '%s' "$ORB_API_KEY" > "$STATE_DIR/api-key"
    chmod 600 "$STATE_DIR/api-key"
    log "saved API key to $STATE_DIR/api-key"
    export ORB_API_KEY
}

auth_hdr=( )
auth_hdr_init() {
    ensure_api_key
    auth_hdr=( -H "Authorization: Bearer ${ORB_API_KEY}" )
}

# --------------------------------------------------------------------------
# Resource helpers
# --------------------------------------------------------------------------

find_or_create_computer() {
    auth_hdr_init
    local listing cid
    listing=$(curl -fsS "${auth_hdr[@]}" "$BASE_URL/v1/computers")
    cid=$(printf '%s' "$listing" \
        | jq -r --arg n "$COMPUTER_NAME" '.computers[]? | select(.name==$n) | .id' \
        | head -1)
    if [[ -n "$cid" && "$cid" != "null" ]]; then
        log "reusing existing computer $cid ($COMPUTER_NAME)"
        printf '%s' "$cid" > "$COMPUTER_FILE"
        return 0
    fi
    log "creating computer $COMPUTER_NAME (runtime=${RUNTIME_MB}MB, disk=${DISK_MB}MB)"
    local resp
    resp=$(curl -fsS "${auth_hdr[@]}" -X POST "$BASE_URL/v1/computers" \
        -H 'Content-Type: application/json' \
        -d "{\"name\":\"$COMPUTER_NAME\",\"runtime_mb\":$RUNTIME_MB,\"disk_mb\":$DISK_MB}")
    cid=$(printf '%s' "$resp" | jq -r .id)
    [[ -n "$cid" && "$cid" != "null" ]] || die "create computer returned no id: $resp"
    printf '%s' "$cid" > "$COMPUTER_FILE"
    log "created computer $cid"
}

upload_config() {
    local cid; cid=$(cat "$COMPUTER_FILE")
    [[ -f "$ORB_TOML_PATH" ]] || die "orb.toml not found at $ORB_TOML_PATH"
    log "uploading orb.toml to computer $cid"
    curl -fsS "${auth_hdr[@]}" -X POST "$BASE_URL/v1/computers/$cid/config" \
        -H 'Content-Type: application/toml' \
        --data-binary "@$ORB_TOML_PATH" \
        | jq -c '.'
}

trigger_build() {
    local cid; cid=$(cat "$COMPUTER_FILE")
    log "triggering build on $cid (this blocks ~30-120s)"
    local provider="${ORB_LLM_PROVIDER:-zai}"

    # Resolve Z.AI key: env var first, cached file as fallback.
    local zai_key="${ANTHROPIC_AUTH_TOKEN:-}"
    if [[ -z "$zai_key" && -f "$STATE_DIR/zai-key" ]]; then
        zai_key=$(cat "$STATE_DIR/zai-key")
    fi

    # Build org_secrets map. Orb resolves "${NAME}" placeholders in
    # orb.toml [agent.env] against these values at build time.
    local secrets='{}'
    if [[ -n "${GITHUB_TOKEN:-}" ]]; then
        secrets=$(jq -n --arg t "$GITHUB_TOKEN" '{GITHUB_TOKEN: $t}')
    fi
    case "$provider" in
        zai)
            [[ -n "$zai_key" ]] || die "ANTHROPIC_AUTH_TOKEN (Z.AI key) required"
            secrets=$(jq -n --argjson base "$secrets" --arg t "$zai_key" \
                '$base + {ANTHROPIC_AUTH_TOKEN: $t}')
            ;;
        anthropic)
            : "${ANTHROPIC_API_KEY:?ANTHROPIC_API_KEY required for provider=anthropic}"
            secrets=$(jq -n --argjson base "$secrets" --arg k "$ANTHROPIC_API_KEY" \
                '$base + {ANTHROPIC_API_KEY: $k}')
            ;;
    esac

    local payload
    payload=$(jq -n --argjson s "$secrets" '{org_secrets: $s}')
    curl -fsS "${auth_hdr[@]}" -X POST "$BASE_URL/v1/computers/$cid/build" \
        -H 'Content-Type: application/json' \
        -d "$payload" \
        | jq -c '.success as $s | {success: $s, steps_n: (.steps | length)}'
}

start_agent() {
    local cid; cid=$(cat "$COMPUTER_FILE")
    local provider="${ORB_LLM_PROVIDER:-zai}"
    local payload

    # Env is plumbed via orb.toml [agent.env] + org_secrets on BUILD.
    # /agents only needs to ask Orb to start the process.
    log "starting agent on $cid (provider=$provider; env applied from orb.toml [agent.env])"
    payload='{"task":"start"}'

    local resp
    resp=$(curl -fsS "${auth_hdr[@]}" -X POST "$BASE_URL/v1/computers/$cid/agents" \
        -H 'Content-Type: application/json' \
        -d "$payload")
    printf '%s\n' "$resp" | jq -c '.'
    local aid; aid=$(printf '%s' "$resp" | jq -r .id 2>/dev/null || true)
    if [[ -n "$aid" && "$aid" != "null" ]]; then
        printf '%s' "$aid" > "$AGENT_FILE"
    fi
    local short_id="${cid:0:8}"
    log "live URL will be https://${short_id}.orbcloud.dev"
}

status() {
    auth_hdr_init
    [[ -f "$COMPUTER_FILE" ]] || die "no computer registered yet; run deploy first"
    local cid; cid=$(cat "$COMPUTER_FILE")
    log "computer $cid"
    curl -fsS "${auth_hdr[@]}" "$BASE_URL/v1/computers/$cid" | jq '.'
    log "agents:"
    curl -fsS "${auth_hdr[@]}" "$BASE_URL/v1/computers/$cid/agents" | jq '.'
}

promote() {
    auth_hdr_init
    [[ -f "$COMPUTER_FILE" ]] || die "no computer registered"
    local cid; cid=$(cat "$COMPUTER_FILE")
    log "promote (wake) agent on $cid"
    curl -fsS "${auth_hdr[@]}" -X POST "$BASE_URL/v1/computers/$cid/agents/promote" \
        -H 'Content-Type: application/json' -d '{"port":8000}' | jq '.'
}

demote() {
    auth_hdr_init
    [[ -f "$COMPUTER_FILE" ]] || die "no computer registered"
    local cid; cid=$(cat "$COMPUTER_FILE")
    log "demote (sleep) agent on $cid"
    curl -fsS "${auth_hdr[@]}" -X POST "$BASE_URL/v1/computers/$cid/agents/demote" \
        -H 'Content-Type: application/json' -d '{"port":8000}' | jq '.'
}

usage_info() {
    auth_hdr_init
    curl -fsS "${auth_hdr[@]}" "$BASE_URL/v1/usage" | jq '.'
}

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

cmd="${1:-deploy}"
case "$cmd" in
    deploy)
        auth_hdr_init
        find_or_create_computer
        upload_config
        trigger_build
        start_agent
        log "✔ deploy complete"
        ;;
    status)  status   ;;
    promote) promote  ;;
    demote)  demote   ;;
    usage)   usage_info ;;
    *)
        cat <<EOF >&2
Usage: $(basename "$0") {deploy|status|promote|demote|usage}
  deploy   create/reuse computer, upload orb.toml, build, start agent
  status   print computer + agents JSON
  promote  wake the sleeping agent (manual)
  demote   sleep the running agent (manual)
  usage    print Orb usage/cost
EOF
        exit 2
        ;;
esac
