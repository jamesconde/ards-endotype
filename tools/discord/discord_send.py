#!/usr/bin/env python
"""Post a message to the ARDS endotypes Discord channel via its webhook.

Usage:
    python tools/discord/discord_send.py "message"
    echo "message" | python tools/discord/discord_send.py

Discord 403s urllib's default User-Agent, so one is set explicitly.
Messages are capped at 2000 characters by Discord; long text is split rather
than truncated, because a silently cut-off status report is worse than two.
"""
import json
import sys
import time
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
WEBHOOK = REPO / "secrets" / "discord_webhook.txt"
LIMIT = 1900


def chunks(text: str, n: int = LIMIT):
    """Split on line boundaries where possible, so code blocks stay readable."""
    while text:
        if len(text) <= n:
            yield text
            return
        cut = text.rfind("\n", 0, n)
        if cut <= 0:
            cut = n
        yield text[:cut]
        text = text[cut:].lstrip("\n")


def send(msg: str) -> int:
    if not WEBHOOK.exists():
        sys.exit(f"{WEBHOOK} missing — see tools/discord/README.md")
    url = WEBHOOK.read_text().strip()
    sent = 0
    for part in chunks(msg):
        req = urllib.request.Request(
            url,
            data=json.dumps({"content": part}).encode(),
            headers={"Content-Type": "application/json",
                     "User-Agent": "ards-endotypes/1.0 (+local)"},
        )
        try:
            urllib.request.urlopen(req)
            sent += 1
        except Exception as e:  # noqa: BLE001
            print(f"send failed: {e}", file=sys.stderr)
            return 1
        time.sleep(0.4)   # stay under the webhook rate limit
    print(f"sent ({sent} message{'s' if sent != 1 else ''})")
    return 0


if __name__ == "__main__":
    text = " ".join(sys.argv[1:]) or sys.stdin.read()
    text = text.strip()
    if not text:
        sys.exit("nothing to send")
    sys.exit(send(text))
