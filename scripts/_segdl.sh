#!/usr/bin/bash
# Robust segmented parallel downloader over HTTP range requests + gzip integrity check.
# usage: _segdl.sh <url> <outfile> [nparts]
# Each segment is re-requested (appending) until it reaches its exact expected byte span,
# because the EBI FTP occasionally returns a short range (truncated body).
set -u
url="$1"; out="$2"; nparts="${3:-10}"
if [ -s "$out" ] && gzip -t "$out" 2>/dev/null; then echo "SKIP_CACHED $(basename "$out")"; exit 0; fi
size=$(curl -sI "$url" | grep -i '^content-length' | tr -dc '0-9')
if [ -z "$size" ] || [ "$size" = "0" ]; then echo "NO_SIZE $url"; exit 3; fi
rm -f "${out}".part*
step=$(( (size + nparts - 1) / nparts ))
echo "START $(basename "$out") size=$size parts=$nparts step=$step $(date +%T)"
worker(){ # s e pp
  local s=$1 e=$2 pp=$3 want=$(( e - s + 1 )) f="${out}.part${pp}" tries=0
  while :; do
    local cur=0; [ -f "$f" ] && cur=$(stat -c%s "$f")
    [ "$cur" -ge "$want" ] && break
    tries=$((tries+1)); [ $tries -gt 8 ] && { echo "PARTFAIL $pp cur=$cur want=$want"; return 1; }
    local rs=$(( s + cur ))
    local got_before=$cur
    curl -s --max-time 1200 -r "${rs}-${e}" "$url" >> "$f" 2>/dev/null
    local after=$(stat -c%s "$f")
    [ "$after" -le "$got_before" ] && sleep 3   # no progress, brief backoff
  done
  # enforce exact length (trim any overrun, should not happen)
  if [ "$(stat -c%s "$f")" -gt "$want" ]; then head -c "$want" "$f" > "${f}.trim"; mv "${f}.trim" "$f"; fi
}
pids=""
for i in $(seq 0 $((nparts-1))); do
  s=$(( i * step )); e=$(( s + step - 1 ))
  [ "$e" -ge "$size" ] && e=$(( size - 1 ))
  [ "$s" -ge "$size" ] && break
  pp=$(printf '%03d' "$i")
  worker "$s" "$e" "$pp" &
  pids="$pids $!"
done
wait $pids
files=$(ls "${out}".part* 2>/dev/null | sort)
[ -z "$files" ] && { echo "NO_PARTS $out"; exit 4; }
cat $files > "$out"
if gzip -t "$out" 2>/dev/null; then
  rm -f "${out}".part*
  echo "OK $(basename "$out") got=$(stat -c%s "$out") want=$size $(date +%T)"
else
  echo "GZIPFAIL $(basename "$out") got=$(stat -c%s "$out") want=$size $(date +%T)"
fi
