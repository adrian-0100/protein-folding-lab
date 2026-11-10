# Running on the remote GPU PC

The laptop (M3, 8 GB) needs 38 minutes for ten predictions of ubiquitin. The remote
GPU PC (an Apple-silicon Mac: M1 Ultra, 20 cores, 128 GB) does twenty in 14, and
leaves the laptop free. How that number was reached is in
[Using the whole machine](#using-the-whole-machine). This folder sets up a
**guest account** there that you can hand to other people, and one command to send
work to it.

The scripts are for a remote machine running macOS (they use `dscl`, `pf` and
`launchd`). Nothing in them is tied to one particular machine: its address and
administrator login live in `remote.env`, which is not in git.

## The design in one paragraph

Nobody gets the remote machine's real login. There is one extra account, `foldlab`,
that exists only to run this project. Each person you let in has their own SSH key
for that account, so you can remove one person without touching the others. The
account is boxed in: it cannot read anyone else's files, cannot use `sudo`, cannot
be used as a tunnel, and a firewall rule stops it from talking to anything on the
machine itself or on the local network. It can reach the internet, which the
pipeline needs (UniProt, the PDB, the MSA server, PyPI).

That last rule is the important one. A shared machine usually runs other services,
and services often trust anything that arrives from the machine itself
(`127.0.0.1`). Without the firewall rule, anyone with a shell there could talk to
them.

## One-time setup

**0. Tell the command where the machine is.** Copy
[remote.env.example](remote.env.example) to `remote.env` in this folder and fill in
the address and the administrator login. Then record the machine's identity and your
own public key (both files are machine-specific and gitignored):

```bash
ssh-keyscan -t ed25519 <address> | sed 's/^[^ ]*/foldlab-remote/' > infra/remote-gpu-pc/known_hosts
```

```bash
awk '{print $1, $2, "me"}' ~/.ssh/id_ed25519.pub > infra/remote-gpu-pc/authorized_keys
```

Check the fingerprint in `known_hosts` against the one the machine itself shows
(`ssh-keygen -lf /etc/ssh/ssh_host_ed25519_key.pub`, run there): that is what stops
anyone else from posing as it later.

**1. Create the account** (you, in your own terminal, on the local network):

```bash
infra/remote-gpu-pc/remote admin-install
```

It copies [setup-guest.sh](setup-guest.sh) to the remote and runs it with `sudo`, so
it asks for the administrator's password. Read the script first; it is the only thing
here that changes the machine's system settings. It ends by testing its own work and
prints `All checks passed` or tells you what failed.

What it changes on the remote:

| Change | Where | Undo |
|---|---|---|
| Account and group `foldlab` (standard, hidden, random password nobody has) | Directory Services, `/Users/foldlab` | `uninstall --delete-account` |
| Who may log in | `/etc/ssh/foldlab_authorized_keys` (owned by root) | `uninstall` |
| SSH rules for `foldlab`: keys only, no forwarding | `/etc/ssh/sshd_config.d/200-foldlab.conf` | `uninstall` |
| SSH rule for connections that come from the machine itself (`127.0.0.1`): only `foldlab`, keys only | same file | `uninstall` |
| Firewall rule for `foldlab` only, loaded at boot | `/etc/pf.anchors/foldlab`, `/etc/pf.foldlab.conf`, `/Library/LaunchDaemons/com.foldlab.pf.plist` | `uninstall` |

Nothing else on the machine is modified, and the administrator logs in from the
local network exactly as before. One thing does change for the administrator:
`ssh <admin>@localhost` typed on the machine itself stops working. That is
deliberate. A Cloudflare Tunnel delivers its connections to SSH from `127.0.0.1`,
so that rule is what guarantees the tunnel can never reach the administrator's
account, whatever key or password someone holds.

The firewall step turns on macOS's packet filter (`pf`), which is off by default;
the only rule added applies to the `foldlab` user.

**2. Install the project there** (uses your key; no password):

```bash
infra/remote-gpu-pc/remote bootstrap
```

**3. Confirm the box holds, from the inside:**

```bash
infra/remote-gpu-pc/remote check
```

Set `REMOTE_PROBE_PORTS` in `remote.env` to the ports of the services that run on
the machine, so the check tries exactly those.

## Daily use

```bash
infra/remote-gpu-pc/remote run experiments/001-ubiquitin-L67A
```

Copies the experiment to the remote, runs the pipeline there, streams the log, and
copies predictions and results back into the same folder on the laptop. Add a step
name (`predict`, `analyze`...) to run just that step. Ctrl-C detaches without
stopping the run; the same command reattaches. Runs there record the machine as
`remote-gpu-pc` in `results/provenance.json` (`REMOTE_MACHINE` changes it), not its
real hostname.

To log in yourself: `infra/remote-gpu-pc/remote ssh`. To use plain
`ssh foldlab-remote`, paste the output of `infra/remote-gpu-pc/remote ssh-config`
into `~/.ssh/config`.

## Using the whole machine

AlphaFold2 runs on the CPU on any Mac: neither it nor ColabFold has a backend for
Apple's GPU, and ColabFold's own speed-ups are for NVIDIA cards only. Measured on
the remote with ubiquitin (76 residues):

| How it is run | 20 predictions | Cores busy (of 20) | Memory |
|---|---|---|---|
| One process at a time (the first run) | 41 min | about 6 | 2.8 GB |
| Ten processes side by side (now the default) | 14 min | 19.8 | about 25 GB |

One prediction alone leaves two thirds of the machine idle, so `remote run` starts
one process per variant and network and runs ten at once (`REMOTE_JOBS`). The
predictions are bit-for-bit the same as in the one-at-a-time run. That is 2.9 times
faster, and close to the ceiling: the same work would take about 12 minutes if all
20 cores were used with no overhead at all.

Things that were tried and made no difference: running without `nice` (same
speed), and limiting each process to one thread (the setting is ignored by this
version of JAX).

What would go further is the GPU, and that means a different predictor, not a
setting. There are structure predictors with an Apple-GPU (MLX) backend, such as
Apple's SimpleFold and community ports of OpenFold3 and MiniFold. None of them is
AlphaFold2, so their output is a different experiment, not a faster copy of this
one. Not tried here.

While a run is going, the machine's CPU is fully used. The runs are started with
`nice`, so other work gets the cores first, but tell whoever shares the machine.
Set `REMOTE_JOBS=3` to leave most of it free.

## Letting someone else in

Best: they run `ssh-keygen -t ed25519` on their own computer and send you the
`.pub` file. Their private key never leaves their machine.

```bash
infra/remote-gpu-pc/remote grant ~/Downloads/maria.pub maria
```

Easiest, for someone who has never used SSH: you make the key pair and send them a
folder with the key, a config file and instructions (in Spanish).

```bash
infra/remote-gpu-pc/remote invite maria
```

The private key in that folder has no passphrase and travels over whatever you send
it with, so use this for people and channels you trust, and delete your copy after.

```bash
infra/remote-gpu-pc/remote who
```

```bash
infra/remote-gpu-pc/remote revoke maria
```

`grant`, `invite` and `revoke` each ask for the administrator's password: changing
who can get in is an administrator's decision, and the key list on the remote is
owned by root so that a guest cannot add their friends.

The local list of public keys is `authorized_keys` in this folder: the record of who
has access. It is kept out of git because it describes one particular machine.

## Reaching it from outside the local network: Cloudflare Tunnel

On the local network the guest account is reached directly at `REMOTE_HOST`. If the
machine is also on a tailnet and Tailscale's own SSH answers there, set
`REMOTE_TAILNET`: when you are away and there is no tunnel, `remote` hops through
your administrator login on the tailnet and from there opens the guest account with
your key, automatically. Guests cannot use that route.

From anywhere else the connection travels through a Cloudflare Tunnel. A small
program on the remote (`cloudflared`) dials out to Cloudflare and holds the line
open; a guest runs the same program on their computer and is carried through it to
the remote's SSH. No port is opened on the router and the machine's own address is
never published. SSH still runs end to end inside the tunnel: Cloudflare carries it
but cannot read it, and the guest's computer checks the machine's identity against
`known_hosts`.

### Free, no account: a quick tunnel

```bash
infra/remote-gpu-pc/remote tunnel-install
```

Run it from the local network. It runs [setup-tunnel.sh](setup-tunnel.sh) on the
remote as the administrator, without `sudo`: installs `cloudflared` with Homebrew,
starts it as a launchd service, and saves the address Cloudflare assigned (a random
name under `trycloudflare.com`) in `infra/remote-gpu-pc/remote.env`. It refuses to
run if the SSH rules from `admin-install` are not in place. There is nothing to sign
up for.

What you give up for free:

- **The address is random and changes** whenever `cloudflared` restarts: a reboot,
  a crash, a macOS update. While it stays up, it keeps its address. After a restart,
  `infra/remote-gpu-pc/remote tunnel-status` reads the new one, saves it, and prints
  the one line your guests have to change in their `~/.ssh/config`.
- **No promise it stays up.** Cloudflare offers quick tunnels for testing and can
  limit or drop them. Good enough for a study group; not something to depend on.

Test it without leaving the local network, by forcing the outside route:

```bash
REMOTE_ROUTE=tunnel infra/remote-gpu-pc/remote check
```

Besides the usual checks, this tries to log in as the administrator through the
tunnel with your own key, the one that works on the local network, and confirms it
is refused.

After that every command picks its route by itself: direct when the remote answers
on the local network, the tunnel otherwise. Invitations created with
`remote invite` after the tunnel exists are written for the tunnel and work from
anywhere, the local network included. Guests need `cloudflared` installed; the
instructions in their folder say how.

### If you later want a permanent address

A domain with its DNS in Cloudflare (about 10 USD a year for a `.com`; Cloudflare's
plan itself is free) turns the random address into one you choose and that never
changes:

```bash
infra/remote-gpu-pc/remote tunnel-install ssh-remote.yourdomain.com
```

This variant prints a Cloudflare address to open in your browser, where you sign in
and pick the domain; it then creates the tunnel and the DNS record. Everything else
works the same. With a domain you can also add Cloudflare Access, a sign-in in
front of the hostname so that only email addresses you list ever reach SSH at all.

### What the tunnel adds on the remote

All inside the administrator's home; nothing needs `sudo`.

| Added | Where | Undo |
|---|---|---|
| `cloudflared` | Homebrew | `brew uninstall cloudflared` |
| Service that keeps it running | `~/Library/LaunchAgents/com.foldlab.cloudflared.plist` | `remote tunnel-uninstall` |
| Its log, which holds the current address | `~/Library/Logs/foldlab-remote.log` | `remote tunnel-uninstall` |
| Named tunnels only: Cloudflare authorisation, tunnel secret, settings, DNS record | `~/.cloudflared/`, your Cloudflare domain | `remote tunnel-uninstall`, then delete the DNS record in the dashboard |

### Other free ways in, if the quick tunnel is not enough

- **Port forward on the router.** A permanent address and nothing for guests to
  install, but you need the router's admin page and a public IP address, and the
  machine's SSH becomes visible to the whole internet. The SSH rules here would need
  one more line first, so that outside addresses can only ever be `foldlab`.
- **Tailscale.** Sharing just this machine with someone is free, if the machine is
  already on a tailnet and its SSH policy admits them.

## What this does not protect against

- **Guests from each other.** Everyone shares one account: they can read, change and
  delete each other's files there, and see each other's processes. Fine for a study
  group; not for people who should not see each other's work.
- **A guest using the machine for something else.** They have a shell, 20 cores and
  internet. The box stops them reaching your data, not mining coins. Grant access to
  people, not to the public.
- **Persistence after revoking.** Removing a key stops new logins. It does not end a
  session that is open, or a job they left running. To be sure after removing
  someone you no longer trust, restart the machine or run `sudo pkill -u foldlab` there.
- **Files that are readable by everyone on the machine.** The check covers home
  folders. External drives and anything under `/Users/Shared` are not covered.
- **Programs that do their own DNS.** Name lookups by normal tools go through
  macOS and work. A tool that queries the router's DNS directly would be blocked by
  the firewall rule; if a download fails with a DNS error, that is why.

## Removing everything

```bash
infra/remote-gpu-pc/remote admin-uninstall --delete-account
```

Without `--delete-account` it removes the access and the rules but keeps the files
in `/Users/foldlab`.
