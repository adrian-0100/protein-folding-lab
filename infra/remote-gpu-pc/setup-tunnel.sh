#!/usr/bin/env bash
# Cloudflare Tunnel for the guest account. Runs ON the remote GPU PC as the administrator's own
# user (no sudo). You normally call it through `infra/remote-gpu-pc/remote tunnel-install`.
#
#   bash setup-tunnel.sh install quick                    free, no account, no domain; random address
#   bash setup-tunnel.sh install ssh-remote.example.com   permanent address; needs a domain in Cloudflare
#   bash setup-tunnel.sh status
#   bash setup-tunnel.sh address                          print just the current address
#   bash setup-tunnel.sh uninstall
#
# How it works: `cloudflared` on this Mac opens an outgoing connection to Cloudflare and keeps it
# open. Someone who runs `cloudflared access ssh --hostname <address>` on their computer is
# carried through that connection to this Mac's SSH. No port is opened on the router and
# the Mac's own address is never published.
#
# Two kinds:
#   quick   Cloudflare hands out a random name under trycloudflare.com. Nothing to sign up for.
#           The name changes whenever cloudflared restarts (a reboot, a crash), and Cloudflare
#           offers quick tunnels for testing, with no promise of availability.
#   named   A name you choose under a domain whose DNS is in Cloudflare. Permanent.
#
# Safe only together with setup-guest.sh: its SSH rules make sure that anything arriving this way
# can log in as the guest account and nothing else. This script refuses to run without them.

set -euo pipefail

ACCOUNT=foldlab
TUNNEL=foldlab-remote
LABEL=com.foldlab.cloudflared
SSHD_DROPIN="/etc/ssh/sshd_config.d/200-$ACCOUNT.conf"
DIR="$HOME/.cloudflared"
CONFIG="$DIR/$TUNNEL.yml"
PLIST="$HOME/Library/LaunchAgents/$LABEL.plist"
LOG="$HOME/Library/Logs/$TUNNEL.log"
ORIGIN="ssh://127.0.0.1:22"
BREW=/opt/homebrew/bin/brew
CF=/opt/homebrew/bin/cloudflared

step() { printf '\n\033[1m▸ %s\033[0m\n' "$*"; }
ok()   { printf '  \033[32m✓\033[0m %s\n' "$*"; }
die()  { printf '\033[31m✗ %s\033[0m\n' "$*" >&2; exit 1; }

[[ $EUID -ne 0 ]] || die "Run this as your own user, not with sudo."

# Prints the named tunnel's id, or nothing if it does not exist (or we never logged in).
tunnel_id() {
  [[ -f "$DIR/cert.pem" ]] || return 0
  "$CF" tunnel list --name "$TUNNEL" --output json 2>/dev/null | python3 -c '
import json, sys
try:
    tunnels = json.load(sys.stdin) or []
    print(tunnels[0]["id"] if tunnels else "")
except Exception:
    print("")' || true
}

# The address guests use right now: the name in the config (named) or the newest one in the log (quick).
address() {
  if [[ -f "$CONFIG" ]]; then
    awk '/hostname:/ {print $3; exit}' "$CONFIG"
  else
    grep -oE '[a-z0-9-]+\.trycloudflare\.com' "$LOG" 2>/dev/null | tail -1 || true
  fi
}

running() { launchctl print "gui/$(id -u)/$LABEL" >/dev/null 2>&1 || launchctl print "user/$(id -u)/$LABEL" >/dev/null 2>&1; }

stop_agent() {
  for domain in "gui/$(id -u)" "user/$(id -u)"; do
    launchctl bootout "$domain/$LABEL" >/dev/null 2>&1 || true
  done
}

# start_agent <arguments to cloudflared...>
start_agent() {
  local arg domain loaded=""
  mkdir -p "$(dirname "$PLIST")" "$(dirname "$LOG")"
  {
    cat <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>Label</key><string>$LABEL</string>
  <key>ProgramArguments</key>
  <array>
    <string>$CF</string>
EOF
    for arg in "$@"; do echo "    <string>$arg</string>"; done
    cat <<EOF
  </array>
  <key>RunAtLoad</key><true/>
  <key>KeepAlive</key><true/>
  <key>ThrottleInterval</key><integer>30</integer>
  <key>StandardOutPath</key><string>$LOG</string>
  <key>StandardErrorPath</key><string>$LOG</string>
</dict>
</plist>
EOF
  } >"$PLIST"
  stop_agent
  sleep 1
  : >"$LOG"
  # With the graphical session open the agent lives in gui/, otherwise in user/.
  for domain in "gui/$(id -u)" "user/$(id -u)"; do
    if launchctl bootstrap "$domain" "$PLIST" >/dev/null 2>&1; then loaded="$domain"; break; fi
  done
  [[ -n "$loaded" ]] || die "launchd did not accept the service. On the Mac: launchctl bootstrap gui/\$(id -u) $PLIST"
  ok "launchd service $LABEL in $loaded (starts at login, restarts if it stops)"
}

wait_connected() {
  step "Waiting for the tunnel to connect"
  for _ in $(seq 1 30); do
    if grep -q "Registered tunnel connection" "$LOG" 2>/dev/null && [[ -n "$(address)" ]]; then
      ok "connected"
      echo
      echo "Guests reach this Mac at: $(address)"
      return 0
    fi
    sleep 2
  done
  die "The tunnel did not connect in 60 s. See $LOG"
}

prepare() {
  step "SSH rules"
  grep -q "^Match Address 127.0.0.0/8,::1" "$SSHD_DROPIN" 2>/dev/null && grep -q "AllowUsers $ACCOUNT" "$SSHD_DROPIN" \
    || die "The guest account's SSH rules are missing. Run 'remote admin-install' first: without them a tunnel would expose every account on this Mac."
  ok "connections through the tunnel can only be $ACCOUNT, with a key"

  step "cloudflared"
  [[ -x "$CF" ]] || "$BREW" install cloudflared
  ok "$("$CF" --version | head -1)"
}

install_quick() {
  prepare
  # cloudflared refuses to open a quick tunnel when a default config file exists.
  [[ ! -f "$DIR/config.yml" && ! -f "$DIR/config.yaml" ]] || die "$DIR/config.yml exists; a quick tunnel cannot run next to it."
  rm -f "$CONFIG"
  step "Quick tunnel (no account, random address)"
  start_agent tunnel --no-autoupdate --url "$ORIGIN"
  wait_connected
  echo "This address changes whenever cloudflared restarts. After a reboot, run: remote tunnel-status"
}

install_named() {
  local hostname="$1" id
  [[ "$hostname" =~ ^[a-z0-9]([a-z0-9-]*[a-z0-9])?(\.[a-z0-9]([a-z0-9-]*[a-z0-9])?)+$ ]] || die "'$hostname' is not a hostname. Use 'quick' for a free random address."
  prepare

  step "Cloudflare login"
  if [[ ! -f "$DIR/cert.pem" ]]; then
    echo "  Open the address below in a browser, sign in to Cloudflare and pick the domain that"
    echo "  $hostname belongs to. This window continues by itself afterwards."
    "$CF" tunnel login
  fi
  [[ -f "$DIR/cert.pem" ]] || die "Cloudflare login did not finish."
  ok "authorised for one domain ($DIR/cert.pem; keep that file private)"

  step "Tunnel"
  id="$(tunnel_id)"
  if [[ -z "$id" ]]; then
    "$CF" tunnel create "$TUNNEL" >/dev/null
    id="$(tunnel_id)"
  fi
  [[ -n "$id" && -f "$DIR/$id.json" ]] || die "The tunnel's credentials file is not on this Mac ($DIR/$id.json). Delete the tunnel '$TUNNEL' in Cloudflare and run this again."
  chmod 600 "$DIR/$id.json" "$DIR/cert.pem"
  ok "$TUNNEL ($id)"

  step "DNS"
  # Fails, on purpose, if the name already points somewhere else.
  "$CF" tunnel route dns "$TUNNEL" "$hostname" 2>&1 | sed 's/^/  /'

  cat >"$CONFIG" <<EOF
# Managed by protein-folding-lab/infra/remote-gpu-pc/setup-tunnel.sh
tunnel: $id
credentials-file: $DIR/$id.json
ingress:
  - hostname: $hostname
    service: $ORIGIN
  - service: http_status:404
EOF
  "$CF" tunnel --config "$CONFIG" ingress validate >/dev/null || die "cloudflared rejected $CONFIG"

  step "Keeping it running"
  start_agent tunnel --no-autoupdate --config "$CONFIG" run
  wait_connected
}

status() {
  if running; then ok "service $LABEL is loaded"; else die "The tunnel service is not running. Install it with: remote tunnel-install"; fi
  if [[ -f "$CONFIG" ]]; then echo "  kind: named (permanent address)"; else echo "  kind: quick (address changes on restart)"; fi
  echo "  address: $(address)"
  echo "  last lines of $LOG:"
  tail -4 "$LOG" 2>/dev/null | cut -c1-160 | sed 's/^/    /' || true
}

uninstall() {
  step "Stopping the tunnel"
  stop_agent
  local named=""
  [[ -f "$CONFIG" ]] && named=yes
  rm -f "$PLIST" "$CONFIG"
  if [[ -n "$(tunnel_id)" ]]; then
    "$CF" tunnel cleanup "$TUNNEL" >/dev/null 2>&1 || true
    "$CF" tunnel delete "$TUNNEL" 2>&1 | sed 's/^/  /'
  fi
  : >"$LOG" 2>/dev/null || true
  ok "tunnel removed; the Mac is no longer reachable from outside the local network"
  if [[ -n "$named" ]]; then
    echo "  The DNS record stays in Cloudflare (it now points at nothing). Delete it in the dashboard:"
    echo "  your domain → DNS → Records."
  fi
}

case "${1:-}" in
  install)
    case "${2:-}" in
      "")    die "Usage: setup-tunnel.sh install quick | <hostname>" ;;
      quick) install_quick ;;
      *)     install_named "$2" ;;
    esac ;;
  status)    status ;;
  address)   address ;;
  uninstall) uninstall ;;
  *) sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
