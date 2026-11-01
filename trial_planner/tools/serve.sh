#!/usr/bin/env bash
# Start / stop the planner locally:  tools/serve.sh start | stop   (http://localhost:8787)
cd "$(dirname "$0")/.." || exit 1
PID=logs/app.pid
case "$1" in
  start) mkdir -p logs
         if curl -s -o /dev/null http://127.0.0.1:8787; then echo "port 8787 already in use; run: $0 stop"; exit 1; fi
         nohup ~/miniconda/envs/trialplanner-r/bin/Rscript -e 'shiny::runApp("app", port = 8787, host = "0.0.0.0", launch.browser = FALSE)' > logs/app.log 2>&1 &
         echo $! > "$PID"
         for i in $(seq 1 30); do curl -s -o /dev/null http://127.0.0.1:8787 && { echo "up: http://localhost:8787"; exit 0; }; sleep 1; done
         echo "did not start; see logs/app.log"; exit 1 ;;
  stop)  [ -f "$PID" ] && kill "$(cat "$PID")" 2>/dev/null && rm -f "$PID" && echo stopped || echo "not running" ;;
  *)     echo "usage: $0 start|stop"; exit 2 ;;
esac
