#!/usr/bin/env bash
# Guest compute account on the remote GPU PC (macOS). Runs ON that machine, as root.
# You normally do not call this directly: `infra/remote-gpu-pc/remote admin-install` copies it over and runs it.
#
#   sudo bash setup-guest.sh install <keys-file>   create the account, lock it down, install the keys
#   sudo bash setup-guest.sh keys <keys-file>      replace the list of people who may log in
#   sudo bash setup-guest.sh check                 prove the lock-down works
#   sudo bash setup-guest.sh uninstall             remove the SSH and firewall rules (keeps the account's files)
#   sudo bash setup-guest.sh uninstall --delete-account   ...and delete the account and everything in its home
#
# What "locked down" means here:
#   - a standard account (not admin, not in `staff`), so it cannot read other users' homes
#   - no password that anyone knows; SSH keys only, and the key list is owned by root
#   - no port forwarding, so the account cannot be used as a tunnel into this machine or the local network
#   - SSH connections that arrive from this machine itself (that is how a Cloudflare Tunnel delivers
#     them) may only log in as this account, with a key: the tunnel can never reach another account
#   - a firewall rule for this one account: it can reach the internet, but not this machine's own
#     services (127.0.0.1), the local network, or a tailnet
#
# What it does not do: keep guests apart from each other. Everyone you grant shares this one account.

set -euo pipefail

ACCOUNT=foldlab
FULL_NAME="Folding Lab (guest compute)"
HOME_DIR="/Users/$ACCOUNT"
KEYS="/etc/ssh/${ACCOUNT}_authorized_keys"
SSHD_DROPIN="/etc/ssh/sshd_config.d/200-$ACCOUNT.conf"
PF_ANCHOR="/etc/pf.anchors/$ACCOUNT"
PF_CONF="/etc/pf.$ACCOUNT.conf"
PF_DAEMON="/Library/LaunchDaemons/com.$ACCOUNT.pf.plist"

step() { printf '\n\033[1m▸ %s\033[0m\n' "$*"; }
ok()   { printf '  \033[32m✓\033[0m %s\n' "$*"; }
bad()  { printf '  \033[31m✗\033[0m %s\n' "$*"; FAILED=1; }
die()  { printf '\033[31m✗ %s\033[0m\n' "$*" >&2; exit 1; }

[[ "$(uname -s)" == Darwin ]] || die "This script is for macOS."
[[ $EUID -eq 0 ]] || die "Run with sudo."

uid_of() { dscl . -read "/Users/$ACCOUNT" UniqueID 2>/dev/null | awk '{print $2}'; }

# --- keys ----------------------------------------------------------------------------------------
install_keys() {
  local source="${1:?keys file missing}" tmp count=0
  [[ -r "$source" ]] || die "Cannot read $source"
  tmp="$(mktemp)"
  while IFS= read -r line || [[ -n "$line" ]]; do
    [[ -z "$line" || "$line" == \#* ]] && continue
    # Each line must be exactly one valid public key: "type base64 comment".
    printf '%s\n' "$line" | ssh-keygen -l -f /dev/stdin >/dev/null 2>&1 || die "Not a public key: ${line:0:60}..."
    [[ "$line" == ssh-* || "$line" == ecdsa-* || "$line" == sk-* ]] || die "Key options are not allowed: ${line:0:40}..."
    # `restrict` turns off every forwarding feature for this key; `pty` gives back an interactive shell.
    printf 'restrict,pty %s\n' "$line" >>"$tmp"
    count=$((count + 1))
  done <"$source"
  [[ $count -gt 0 ]] || die "$source has no keys; refusing to lock everyone out."
  install -o root -g wheel -m 644 "$tmp" "$KEYS"
  rm -f "$tmp"
  ok "$count key(s) may log in as $ACCOUNT:"
  ssh-keygen -l -f "$KEYS" | sed 's/^/      /'
}

# --- account -------------------------------------------------------------------------------------
install_account() {
  local gid uid
  if ! dscl . -read "/Groups/$ACCOUNT" >/dev/null 2>&1; then
    gid=600
    while dscl . -list /Groups PrimaryGroupID | awk '{print $2}' | grep -qx "$gid"; do gid=$((gid + 1)); done
    dseditgroup -o create -i "$gid" -r "$FULL_NAME" "$ACCOUNT"
  fi
  gid="$(dscl . -read "/Groups/$ACCOUNT" PrimaryGroupID | awk '{print $2}')"

  if ! id "$ACCOUNT" >/dev/null 2>&1; then
    # A random password that is never shown or stored: nobody can log in with it.
    sysadminctl -addUser "$ACCOUNT" -fullName "$FULL_NAME" -home "$HOME_DIR" -shell /bin/zsh \
      -password "$(openssl rand -base64 30)" 2>&1 | sed 's/^/      /'
  fi
  uid="$(uid_of)"
  [[ -n "$uid" ]] || die "Could not create the account."

  # Its own primary group instead of `staff`: other users' homes are group-readable by staff.
  dscl . -create "/Users/$ACCOUNT" PrimaryGroupID "$gid"
  dseditgroup -o edit -d "$ACCOUNT" -t user staff >/dev/null 2>&1 || true
  dseditgroup -o edit -d "$ACCOUNT" -t user admin >/dev/null 2>&1 || true
  dscl . -create "/Users/$ACCOUNT" IsHidden 1    # not shown on the login window

  mkdir -p "$HOME_DIR"
  chown -R "$uid:$gid" "$HOME_DIR"
  chmod 700 "$HOME_DIR"

  # macOS only lets members of this group use Remote Login when it is restricted to "only these users".
  if dscl . -read /Groups/com.apple.access_ssh >/dev/null 2>&1; then
    dseditgroup -o edit -a "$ACCOUNT" -t user com.apple.access_ssh
  fi
  ok "account $ACCOUNT (uid $uid, group $gid), home $HOME_DIR"
}

# --- sshd ----------------------------------------------------------------------------------------
install_sshd() {
  cat >"$SSHD_DROPIN" <<EOF
# Managed by protein-folding-lab/infra/remote-gpu-pc/setup-guest.sh.
Match User $ACCOUNT
    AuthorizedKeysFile $KEYS
    AuthenticationMethods publickey
    PasswordAuthentication no
    KbdInteractiveAuthentication no
    AllowTcpForwarding no
    AllowStreamLocalForwarding no
    AllowAgentForwarding no
    PermitTunnel no
    X11Forwarding no
    PermitUserRC no

# A Cloudflare Tunnel hands its connections to sshd from 127.0.0.1, so "from this machine" means
# "possibly from the internet". Those connections may only be $ACCOUNT, and never with a password.
Match Address 127.0.0.0/8,::1
    AllowUsers $ACCOUNT
    PasswordAuthentication no
    KbdInteractiveAuthentication no
Match all
EOF
  chmod 644 "$SSHD_DROPIN"
  if ! /usr/sbin/sshd -t 2>/tmp/sshd-check.$$; then
    rm -f "$SSHD_DROPIN"
    cat /tmp/sshd-check.$$ >&2
    rm -f /tmp/sshd-check.$$
    die "sshd rejected the new rules; they were removed and nothing changed."
  fi
  rm -f /tmp/sshd-check.$$
  ok "sshd: $ACCOUNT is key-only with no forwarding; connections from this machine itself are $ACCOUNT-only"
  echo "      (macOS starts sshd per connection, so this is already live)"
}

# --- firewall ------------------------------------------------------------------------------------
install_firewall() {
  local uid; uid="$(uid_of)"
  cat >"$PF_ANCHOR" <<EOF
# Managed by protein-folding-lab/infra/remote-gpu-pc/setup-guest.sh.
# Connections opened by the $ACCOUNT account (uid $uid) may go to the internet and nowhere else:
# not to this machine's own services, not to the local network, not to a tailnet.
table <${ACCOUNT}_private> const { 127.0.0.0/8, 10.0.0.0/8, 172.16.0.0/12, 192.168.0.0/16, 169.254.0.0/16, 100.64.0.0/10, ::1, fe80::/10, fc00::/7, fd7a:115c:a1e0::/48 }
block return out quick proto { tcp, udp } from any to <${ACCOUNT}_private> user $uid
# Traffic to any of this machine's own addresses travels over lo0, whatever the address family.
block return out quick on lo0 proto { tcp, udp } user $uid
EOF
  # The stock /etc/pf.conf, plus our anchor. Kept as a separate file so macOS updates do not erase it.
  cat >"$PF_CONF" <<EOF
scrub-anchor "com.apple/*"
nat-anchor "com.apple/*"
rdr-anchor "com.apple/*"
dummynet-anchor "com.apple/*"
anchor "$ACCOUNT"
anchor "com.apple/*"
load anchor "com.apple" from "/etc/pf.anchors/com.apple"
load anchor "$ACCOUNT" from "$PF_ANCHOR"
EOF
  /sbin/pfctl -n -f "$PF_CONF" 2>/dev/null || { rm -f "$PF_ANCHOR" "$PF_CONF"; die "pf rejected the firewall rules; they were removed."; }

  cat >"$PF_DAEMON" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>Label</key><string>com.$ACCOUNT.pf</string>
  <key>ProgramArguments</key>
  <array>
    <string>/bin/sh</string>
    <string>-c</string>
    <string>/sbin/pfctl -f $PF_CONF; /sbin/pfctl -e; true</string>
  </array>
  <key>RunAtLoad</key><true/>
</dict>
</plist>
EOF
  chown root:wheel "$PF_DAEMON"; chmod 644 "$PF_DAEMON"
  launchctl bootout system "$PF_DAEMON" >/dev/null 2>&1 || true
  launchctl bootstrap system "$PF_DAEMON"
  /sbin/pfctl -f "$PF_CONF" 2>/dev/null || die "Could not load the firewall rules."
  /sbin/pfctl -e 2>/dev/null || true    # "already enabled" is fine
  ok "firewall: $ACCOUNT can reach the internet only (reloaded at every boot)"
}

# --- check ---------------------------------------------------------------------------------------
as_guest() { sudo -n -u "$ACCOUNT" -H "$@"; }

check() {
  FAILED=0
  local uid groups effective home port
  uid="$(uid_of)"; [[ -n "$uid" ]] || die "The $ACCOUNT account does not exist."

  step "Account"
  groups="$(id -Gn "$ACCOUNT")"
  if [[ " $groups " == *" admin "* || " $groups " == *" staff "* ]]; then bad "is in admin or staff: $groups"; else ok "not in admin or staff"; fi
  if as_guest sudo -n true 2>/dev/null; then bad "can use sudo"; else ok "cannot use sudo"; fi

  step "Other people's files"
  for home in /Users/*; do
    [[ "$home" == "$HOME_DIR" || "$home" == /Users/Shared || ! -d "$home" ]] && continue
    if as_guest ls "$home" >/dev/null 2>&1; then bad "can list $home"; else ok "cannot list $home"; fi
  done

  step "SSH rules for $ACCOUNT"
  effective="$(/usr/sbin/sshd -T -C "user=$ACCOUNT,host=x,addr=203.0.113.1" 2>/dev/null)"
  for want in "passwordauthentication no" "kbdinteractiveauthentication no" "allowtcpforwarding no" \
              "allowagentforwarding no" "authenticationmethods publickey"; do
    if grep -qix "$want" <<<"$effective"; then ok "$want"; else bad "expected: $want"; fi
  done
  [[ -s "$KEYS" ]] && ok "$(grep -c . "$KEYS") key(s) in $KEYS (root-owned)" || bad "no keys installed"

  step "Network reach of $ACCOUNT"
  for port in $(lsof -nP -iTCP -sTCP:LISTEN 2>/dev/null | awk 'NR>1 {n=split($9,a,":"); print a[n]}' | sort -un); do
    [[ "$port" == 22 ]] && continue
    if as_guest nc -z -G 2 127.0.0.1 "$port" >/dev/null 2>&1; then bad "can connect to local port $port"; else ok "local port $port is blocked"; fi
  done
  if as_guest curl -fsS -m 15 -o /dev/null https://rest.uniprot.org/uniprotkb/P0CG48.fasta; then ok "internet works (UniProt)"; else bad "cannot reach the internet"; fi

  step "Other accounts"
  local other="${SUDO_USER:-root}"
  effective="$(/usr/sbin/sshd -T -C "user=$other,host=x,addr=127.0.0.1" 2>/dev/null)"
  if grep -qix "allowusers $ACCOUNT" <<<"$effective" && grep -qix "passwordauthentication no" <<<"$effective"; then
    ok "through a tunnel (from 127.0.0.1), $other is refused: only $ACCOUNT is allowed, keys only"
  else
    bad "through a tunnel (from 127.0.0.1), accounts other than $ACCOUNT could log in"
  fi
  effective="$(/usr/sbin/sshd -T -C "user=$other,host=x,addr=192.168.1.20" 2>/dev/null | grep -iE '^(passwordauthentication|allowtcpforwarding) ' | tr '\n' ' ')"
  if /usr/sbin/sshd -T -C "user=$other,host=x,addr=192.168.1.20" 2>/dev/null | grep -qi '^allowusers '; then
    bad "$other is restricted on the local network too; that was not intended"
  else
    ok "from the local network, $other is unchanged: $effective"
  fi

  echo
  if [[ $FAILED -eq 0 ]]; then printf '\033[32mAll checks passed.\033[0m\n'; else printf '\033[31mSome checks failed. Do not share access until they pass.\033[0m\n'; exit 1; fi
}

# --- uninstall -----------------------------------------------------------------------------------
uninstall() {
  step "Removing SSH and firewall rules"
  launchctl bootout system "$PF_DAEMON" >/dev/null 2>&1 || true
  rm -f "$PF_DAEMON" "$PF_CONF" "$PF_ANCHOR" "$SSHD_DROPIN" "$KEYS"
  /sbin/pfctl -f /etc/pf.conf 2>/dev/null || true
  ok "rules removed; nobody can log in as $ACCOUNT any more"
  if [[ "${1:-}" == --delete-account ]]; then
    step "Deleting the account and $HOME_DIR"
    sysadminctl -deleteUser "$ACCOUNT" 2>&1 | sed 's/^/      /'
    dseditgroup -o delete "$ACCOUNT" >/dev/null 2>&1 || true
  else
    echo "  The account and its files in $HOME_DIR are still there. Add --delete-account to remove them."
  fi
}

case "${1:-}" in
  install)
    # The firewall goes in before any key works, so the account is never reachable unboxed.
    step "Account";  install_account
    step "Firewall"; install_firewall
    step "Keys";     install_keys "${2:-}"
    step "SSH";      install_sshd
    check
    ;;
  keys)      step "Keys"; install_keys "${2:-}" ;;
  check)     check ;;
  uninstall) uninstall "${2:-}" ;;
  *) sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
