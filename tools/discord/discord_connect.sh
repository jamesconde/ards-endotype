#!/bin/bash
# Bring up the ARDS endotypes Discord bridge. Idempotent.
#
# Launches the listener DETACHED so it survives the Claude session closing.
# After running this, the Claude session Monitors secrets/discord_inbox.txt and
# replies with tools/discord/discord_send.py.
cd "$(dirname "$0")/../.." || exit 1

PY=/home/james/miniconda/envs/ards/bin/python
SEC=secrets

for f in discord_bot_token.txt discord_webhook.txt discord_channel.txt; do
  if [ ! -s "$SEC/$f" ]; then
    echo "ERROR: $SEC/$f missing or empty — see tools/discord/README.md"
    exit 1
  fi
done

if ! "$PY" -c "import discord" 2>/dev/null; then
  echo "installing discord.py into the ards env..."
  "$PY" -m pip install -q discord.py || { echo "pip install failed"; exit 1; }
fi

# Match on this repo's path so we never mistake the MolMAGE listener for ours.
if pgrep -f "ards-endotypes/tools/discord/discord_b[o]t\.py" >/dev/null; then
  echo "ARDS listener already running (pid $(pgrep -f 'ards-endotypes/tools/discord/discord_b[o]t\.py' | head -1))"
else
  touch "$SEC/discord_inbox.txt"
  DISCORD_BOT_TOKEN=$(cat "$SEC/discord_bot_token.txt") \
  DISCORD_CHANNEL=$(cat "$SEC/discord_channel.txt") \
  nohup setsid "$PY" tools/discord/discord_bot.py \
      > "$SEC/discord_bot.log" 2>&1 < /dev/null &
  disown 2>/dev/null || true
  sleep 4
  echo "launched. Log:"
  tail -4 "$SEC/discord_bot.log"
fi

echo ""
echo "SEND:  $PY tools/discord/discord_send.py \"message\""
echo "READ:  (in the Claude session) Monitor secrets/discord_inbox.txt"
echo ""
if pgrep -f "molmage/scripts/discord_b[o]t\.py" >/dev/null; then
  echo "NOTE: the MolMAGE listener is also running and is NOT channel-restricted,"
  echo "      so it will also log this channel's messages into its own inbox."
  echo "      Harmless, but to stop the cross-talk add to molmage/scripts/discord_connect.sh:"
  echo "        DISCORD_CHANNEL=<molmage channel id>"
fi
