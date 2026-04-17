#!/bin/bash
# ============================================================
# orb-antibiotic-scientist -- Control Script
# Usage: ./control.sh {start|stop|status|logs|findings}
# ============================================================

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
PID_DIR="$PROJECT_DIR/logs/pids"
LOG_DIR="$PROJECT_DIR/logs"
LOCK_FILE="$LOG_DIR/watchdog.lock"

case "$1" in
    start)
        if [ -f "$LOCK_FILE" ] && kill -0 "$(cat "$LOCK_FILE")" 2>/dev/null; then
            echo "Agent already running (PID: $(cat "$LOCK_FILE"))"
            exit 1
        fi
        echo "Starting orb-antibiotic-scientist..."
        nohup bash "$PROJECT_DIR/src/watchdog.sh" >> "$LOG_DIR/watchdog.log" 2>&1 &
        sleep 1
        if [ -f "$LOCK_FILE" ]; then
            echo "Agent started (PID: $(cat "$LOCK_FILE"))"
            echo "Logs: tail -f $LOG_DIR/watchdog.log"
        else
            echo "Failed to start. Check $LOG_DIR/watchdog.log"
        fi
        ;;

    stop)
        if [ -f "$LOCK_FILE" ]; then
            PID=$(cat "$LOCK_FILE")
            if kill -0 "$PID" 2>/dev/null; then
                echo "Stopping agent (PID: $PID)..."
                kill "$PID"
                sleep 2
                if kill -0 "$PID" 2>/dev/null; then
                    kill -9 "$PID"
                fi
                rm -f "$LOCK_FILE"
                echo "Agent stopped."
            else
                echo "Agent not running (stale PID file). Cleaning up."
                rm -f "$LOCK_FILE"
            fi
        else
            echo "Agent not running."
        fi
        for pid_file in "$PID_DIR"/*.pid; do
            if [ -f "$pid_file" ]; then
                PID=$(cat "$pid_file")
                kill "$PID" 2>/dev/null
                rm -f "$pid_file"
            fi
        done
        ;;

    status)
        echo "=== orb-antibiotic-scientist Status ==="
        echo ""
        if [ -f "$LOCK_FILE" ] && kill -0 "$(cat "$LOCK_FILE")" 2>/dev/null; then
            echo "Watchdog:     RUNNING (PID: $(cat "$LOCK_FILE"))"
        else
            echo "Watchdog:     STOPPED"
        fi
        for pid_file in "$PID_DIR"/*.pid; do
            if [ -f "$pid_file" ]; then
                NAME=$(basename "$pid_file" .pid)
                PID=$(cat "$pid_file")
                if kill -0 "$PID" 2>/dev/null; then
                    echo "$NAME:  RUNNING (PID: $PID)"
                else
                    echo "$NAME:  STOPPED (stale PID)"
                fi
            fi
        done
        echo ""
        CAND=$(find "$PROJECT_DIR/findings/candidates" -maxdepth 1 -type d 2>/dev/null | tail -n +2 | wc -l | tr -d ' ')
        DOCK=$(find "$PROJECT_DIR/findings/docking" -type f 2>/dev/null | wc -l | tr -d ' ')
        ADMET=$(find "$PROJECT_DIR/findings/admet" -type f 2>/dev/null | wc -l | tr -d ' ')
        NOV=$(find "$PROJECT_DIR/findings/novelty" -type f 2>/dev/null | wc -l | tr -d ' ')
        MECH=$(find "$PROJECT_DIR/findings/mechanism" -type f 2>/dev/null | wc -l | tr -d ' ')
        RED=$(find "$PROJECT_DIR/findings/red-team" -type f 2>/dev/null | wc -l | tr -d ' ')
        RETRO=$(find "$PROJECT_DIR/findings/retrosynthesis" -type f 2>/dev/null | wc -l | tr -d ' ')
        WEEK=$(find "$PROJECT_DIR/findings/weekly-reports" -type f 2>/dev/null | wc -l | tr -d ' ')
        echo "--- Output Stats ---"
        echo "Candidates:           $CAND"
        echo "Docking runs:         $DOCK"
        echo "ADMET reports:        $ADMET"
        echo "Novelty reports:      $NOV"
        echo "Mechanism reports:    $MECH"
        echo "Red-team critiques:   $RED"
        echo "Retrosynthesis plans: $RETRO"
        echo "Weekly reports:       $WEEK"
        echo ""
        if [ -f "$LOG_DIR/watchdog.log" ]; then
            echo "--- Last 5 watchdog entries ---"
            tail -5 "$LOG_DIR/watchdog.log"
        fi
        ;;

    logs)
        WHICH="${2:-watchdog}"
        LOG_FILE="$LOG_DIR/$WHICH.log"
        if [ -f "$LOG_FILE" ]; then
            tail -f "$LOG_FILE"
        else
            echo "No log file: $LOG_FILE"
            echo "Available: watchdog, agent-sdk"
        fi
        ;;

    findings)
        echo "=== Recent Findings ==="
        echo ""
        for dir in candidates docking admet novelty mechanism red-team retrosynthesis weekly-reports; do
            COUNT=$(find "$PROJECT_DIR/findings/$dir" -type f 2>/dev/null | wc -l | tr -d ' ')
            echo "--- $dir ($COUNT files) ---"
            for f in "$PROJECT_DIR/findings/$dir/"*.md; do
                if [ -f "$f" ]; then
                    echo ""
                    head -6 "$f"
                fi
            done
            echo ""
        done
        if [ -f "$PROJECT_DIR/findings/leaderboard.json" ]; then
            echo "=== Leaderboard (top 10) ==="
            python3 -c "
import json, sys
try:
    d = json.load(open('$PROJECT_DIR/findings/leaderboard.json'))
    entries = d.get('candidates', d) if isinstance(d, dict) else d
    for i, e in enumerate(entries[:10], 1):
        name = e.get('id') or e.get('name') or '?'
        score = e.get('rigor_score') or e.get('composite_score') or '-'
        print(f'  {i:>2}. {name:<20}  score={score}')
except Exception as exc:
    print(f'leaderboard.json not readable: {exc}')
"
        fi
        ;;

    *)
        echo "orb-antibiotic-scientist -- Control"
        echo ""
        echo "Usage: $0 {start|stop|status|logs [component]|findings}"
        echo ""
        echo "  start              Start the antibiotic research agent"
        echo "  stop               Stop the antibiotic research agent"
        echo "  status             Show agent status and stats"
        echo "  logs [component]   Tail logs (watchdog|agent-sdk)"
        echo "  findings           Show recent findings and leaderboard"
        ;;
esac
