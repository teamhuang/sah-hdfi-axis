import os, gzip, tarfile, time, requests
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scrna")
os.makedirs(DATA, exist_ok=True)
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311852/suppl/GSE311852_RAW.tar"
HDR = {"User-Agent": "Mozilla/5.0"}
path = os.path.join(DATA, "GSE311852_RAW.tar")
exp = int(requests.head(URL, headers=HDR, timeout=60).headers.get("Content-Length", 0))
print("expected", exp, flush=True)
for attempt in range(1, 8):
    try:
        if os.path.exists(path) and os.path.getsize(path) == exp:
            break
        resume = os.path.getsize(path) if os.path.exists(path) else 0
        h = dict(HDR); h["Range"] = f"bytes={resume}-" if resume else ""
        with requests.get(URL, headers=h, stream=True, timeout=180) as r:
            if r.status_code not in (200, 206): raise IOError(r.status_code)
            mode = "ab" if resume and r.status_code == 206 else "wb"
            with open(path, mode) as f:
                for c in r.iter_content(1 << 20):
                    if c: f.write(c)
        got = os.path.getsize(path)
        print(f"try{attempt}: {got}/{exp}", flush=True)
        if got == exp: break
    except Exception as e:
        print("err", repr(e)[:100], flush=True); time.sleep(5)
print("download complete" if os.path.getsize(path) == exp else "INCOMPLETE", flush=True)
if os.path.getsize(path) == exp:
    with tarfile.open(path) as t:
        t.extractall(DATA)
    print("extracted:", sorted(os.listdir(DATA))[:20], flush=True)
