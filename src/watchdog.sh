#!/bin/bash
# ============================================================
# orb-antibiotic-scientist -- Watchdog
#
# Runs the Agent SDK-based agent (auto-compacts, runs forever).
# If the process dies for any reason, restart it.
# Stall detection: if no new files in 5 hours, kill and restart.
# ============================================================

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
LOG_DIR="$PROJECT_DIR/logs"
PID_DIR="$LOG_DIR/pids"
LOCK_FILE="$LOG_DIR/watchdog.lock"

# Stall detection: kill if no new files in this many seconds.
# Default 30 min — a healthy run writes candidate/docking files every few
# minutes; anything longer almost certainly means the SDK stream wedged
# (e.g. Z.AI timeout mid-stream) and agent.py's per-message timeout
# already tried to recover. Override via STALL_TIMEOUT env.
STALL_TIMEOUT="${STALL_TIMEOUT:-1800}"

mkdir -p "$PID_DIR" "$LOG_DIR"

# Prevent multiple watchdogs. PID 1 is treated as stale (old container init).
if [ -f "$LOCK_FILE" ]; then
    EXISTING_PID=$(cat "$LOCK_FILE")
    if [ "$EXISTING_PID" != "1" ] && [ "$EXISTING_PID" != "$$" ] && kill -0 "$EXISTING_PID" 2>/dev/null; then
        if ps -p "$EXISTING_PID" -o comm= 2>/dev/null | grep -q bash; then
            echo "[watchdog] Already running (PID: $EXISTING_PID). Exiting."
            exit 1
        fi
    fi
    rm -f "$LOCK_FILE"
fi

echo $$ > "$LOCK_FILE"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [watchdog] $1" | tee -a "$LOG_DIR/watchdog.log"
}

cleanup() {
    log "Shutting down..."
    if [ -f "$PID_DIR/agent-sdk.pid" ]; then
        AGENT_PID=$(cat "$PID_DIR/agent-sdk.pid")
        kill "$AGENT_PID" 2>/dev/null
        sleep 2
        kill -9 "$AGENT_PID" 2>/dev/null
        rm -f "$PID_DIR/agent-sdk.pid"
    fi
    rm -f "$LOCK_FILE"
    log "Stopped."
    exit 0
}

trap cleanup SIGINT SIGTERM

log "============================================"
log "orb-antibiotic-scientist -- Agent SDK Mode"
log "Auto-compaction enabled. Runs indefinitely."
log "Stall detection: restart if idle >${STALL_TIMEOUT}s."
log "============================================"

while true; do
    log "Starting agent (Agent SDK with auto-compaction)..."

    cd "$PROJECT_DIR"

    python3 "$PROJECT_DIR/src/agent.py" &

    AGENT_PID=$!
    log "Agent PID: $AGENT_PID"

    LAST_ACTIVITY=$(date +%s)

    while kill -0 "$AGENT_PID" 2>/dev/null; do
        sleep 60

        RECENT=$(find "$PROJECT_DIR/findings" -type f -mmin -5 2>/dev/null | head -1)

        if [ -n "$RECENT" ]; then
            LAST_ACTIVITY=$(date +%s)
        fi

        NOW=$(date +%s)
        IDLE_TIME=$(( NOW - LAST_ACTIVITY ))

        if [ "$IDLE_TIME" -ge "$STALL_TIMEOUT" ]; then
            log "STALL DETECTED: No new files in $(( IDLE_TIME / 3600 ))h. Killing and restarting..."
            kill "$AGENT_PID" 2>/dev/null
            sleep 2
            kill -9 "$AGENT_PID" 2>/dev/null
            break
        fi
    done

    wait "$AGENT_PID" 2>/dev/null
    EXIT_CODE=$?

    CANDIDATES=$(find "$PROJECT_DIR/findings/candidates" -maxdepth 1 -type d 2>/dev/null | tail -n +2 | wc -l | tr -d ' ')
    DOCKING=$(find "$PROJECT_DIR/findings/docking" -type f 2>/dev/null | wc -l | tr -d ' ')
    ADMET=$(find "$PROJECT_DIR/findings/admet" -type f 2>/dev/null | wc -l | tr -d ' ')
    NOVELTY=$(find "$PROJECT_DIR/findings/novelty" -type f 2>/dev/null | wc -l | tr -d ' ')

    log "Agent exited (code: ${EXIT_CODE}). Candidates: ${CANDIDATES} | Docking: ${DOCKING} | ADMET: ${ADMET} | Novelty: ${NOVELTY}"

    # Auto-commit + push findings (optional, only if GITHUB_TOKEN + remote configured)
    if [ -n "$GITHUB_TOKEN" ] && [ -d "$PROJECT_DIR/.git" ]; then
        cd "$PROJECT_DIR"
        git config user.email "agent@orbcloud.dev" 2>/dev/null
        git config user.name "orb-antibiotic-scientist" 2>/dev/null
        REMOTE_URL=$(git remote get-url origin 2>/dev/null)
        if [ -n "$REMOTE_URL" ]; then
            AUTH_URL=$(echo "$REMOTE_URL" | sed "s|https://github.com|https://${GITHUB_TOKEN}@github.com|")
            git remote set-url origin "$AUTH_URL" 2>/dev/null
        fi
        git add findings/ 2>/dev/null
        CHANGED=$(git diff --cached --numstat 2>/dev/null | wc -l | tr -d ' ')
        if [ "$CHANGED" -gt 0 ]; then
            git commit -m "agent run: ${CANDIDATES} candidates, ${DOCKING} docking, ${ADMET} ADMET, ${NOVELTY} novelty" 2>&1 | tail -3 | while read line; do log "  git: $line"; done
            git push origin HEAD 2>&1 | tail -3 | while read line; do log "  push: $line"; done
        else
            log "No changes to commit."
        fi
    fi

    log "Restarting in 30s..."
    sleep 30
done
