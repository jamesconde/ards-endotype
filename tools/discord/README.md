# Discord bridge — ARDS endotypes session

Two-way bridge so this Claude Code session can be driven from Discord (phone-friendly).
Mirrors the MolMAGE bridge, with one deliberate difference: **this listener is
channel-restricted and refuses to start without a channel id.**

## Why a separate channel is necessary, and not sufficient

The MolMAGE listener runs the *same bot* with **no channel filter**, so it logs every
channel the bot can see into its own inbox. Creating a new channel therefore separates
where *you* type, but not where messages *land* — MolMAGE's inbox would still receive
this project's traffic.

Full separation needs three things, and this bridge provides the second and third:

1. A separate channel (you create it).
2. A channel-restricted listener writing to **this** repo's inbox (`discord_bot.py`).
3. A webhook belonging to **that** channel, so replies go back to the right place.

Optionally, add `DISCORD_CHANNEL=<molmage channel id>` to
`molmage/scripts/discord_connect.sh` so MolMAGE stops logging this channel too. That
change is on the MolMAGE side and has deliberately not been made from here.

## One-time setup

1. **Create the channel** in your server, e.g. `#ards-endotypes`.

2. **Get the channel id.** Discord → User Settings → Advanced → enable **Developer
   Mode**. Then right-click the channel → **Copy Channel ID**.
   ```bash
   echo '<channel id>' > secrets/discord_channel.txt
   ```

3. **Create a webhook for that channel.** Channel → Edit Channel → Integrations →
   Webhooks → New Webhook → Copy Webhook URL.
   ```bash
   echo '<webhook url>' > secrets/discord_webhook.txt
   ```

4. **Bot token** — reuse the existing bot, already in the server:
   ```bash
   cp /home/james/molmage/outputs/discord_bot_token.txt secrets/discord_bot_token.txt
   ```
   Make sure the bot has **View Channel** and **Read Message History** on the new
   channel. It already has the Message Content Intent enabled.

5. **Start it.**
   ```bash
   bash tools/discord/discord_connect.sh
   ```
   Expect `bot online as ...` in `secrets/discord_bot.log`.

## Everyday use

| Action | Command |
|---|---|
| Send to Discord | `python tools/discord/discord_send.py "message"` |
| Read from Discord | the Claude session Monitors `secrets/discord_inbox.txt` |
| Restart listener | `bash tools/discord/discord_connect.sh` (idempotent) |
| Stop listener | `pkill -f 'ards-endotypes/tools/discord/discord_bot.py'` |

## Limitations, stated plainly

- **The listener does not wake Claude up.** It only appends to the inbox. Nothing
  responds unless a live Claude session is Monitoring that file. Messages sent while no
  session is running are queued in the file and read when a session next watches it.
- Long replies are split across multiple Discord messages rather than truncated, since a
  silently cut-off status report is worse than two messages.

## Secrets

`secrets/` is gitignored, and `tools/repo_guard.py` additionally scans file *contents*
for Discord webhook URLs and bot tokens on every commit — a gitignore is one `git add -f`
away from being bypassed. Both files are credentials: the webhook lets anyone post to the
channel, and the token controls the bot.
