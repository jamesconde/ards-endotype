#!/usr/bin/env python
"""Discord listener for the ARDS endotypes session — CHANNEL-RESTRICTED.

Appends incoming messages from ONE channel to secrets/discord_inbox.txt, which a
Claude Code session Monitors. Replies go out through the webhook
(tools/discord/discord_send.py).

Why the channel restriction is not optional here
------------------------------------------------
The MolMAGE bridge runs the same bot with no channel filter, so it captures every
channel the bot can see. Without a filter on this side too, the two projects'
conversations interleave in both inboxes and neither session can tell which
project a message was meant for. `DISCORD_CHANNEL` is therefore REQUIRED, and
this script refuses to start without it.

Run:
    DISCORD_BOT_TOKEN=xxx DISCORD_CHANNEL=<id> python tools/discord/discord_bot.py
Normally you just use:
    bash tools/discord/discord_connect.sh
"""
import datetime
import os
import sys
from pathlib import Path

import discord

REPO = Path(__file__).resolve().parent.parent.parent
INBOX = REPO / "secrets" / "discord_inbox.txt"

ONLY_CHANNEL = os.environ.get("DISCORD_CHANNEL", "").strip()

intents = discord.Intents.default()
intents.message_content = True
client = discord.Client(intents=intents)


@client.event
async def on_ready():
    print(f"bot online as {client.user}", flush=True)
    print(f"  channel filter : {ONLY_CHANNEL}", flush=True)
    print(f"  inbox          : {INBOX}", flush=True)


@client.event
async def on_message(msg):
    if msg.author == client.user:
        return
    # Ignore our own webhook posts, or replies echo back and loop.
    if msg.webhook_id is not None:
        return
    if str(msg.channel.id) != ONLY_CHANNEL:
        return
    line = (f"{datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')} "
            f"#{msg.channel} {msg.author}: {msg.content}\n")
    INBOX.parent.mkdir(parents=True, exist_ok=True)
    with INBOX.open("a") as f:
        f.write(line)
    print("inbox <- " + line.strip(), flush=True)


def main():
    token = os.environ.get("DISCORD_BOT_TOKEN")
    if not token:
        sys.exit("set DISCORD_BOT_TOKEN")
    if not ONLY_CHANNEL:
        sys.exit(
            "DISCORD_CHANNEL is required for this bridge.\n"
            "Running unfiltered would mix this project's messages into the MolMAGE "
            "inbox and vice versa. Put the channel id in secrets/discord_channel.txt "
            "(Discord: User Settings > Advanced > Developer Mode, then right-click the "
            "channel > Copy Channel ID)."
        )
    client.run(token)


if __name__ == "__main__":
    main()
