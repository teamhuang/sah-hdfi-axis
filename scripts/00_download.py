import os, gzip, time, requests

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
os.makedirs(DATA, exist_ok=True)

BASE = "https://ftp.ncbi.nlm.nih.gov/geo"
FILES = {
    # discovery: acute SAH (ruptured aneurysm) blood, GPL10558
    "GSE36791_series_matrix.txt.gz": f"{BASE}/series/GSE36nnn/GSE36791/matrix/GSE36791_series_matrix.txt.gz",
    # validation: past aSAH blood, same platform
    "GSE73378-GPL10558_series_matrix.txt.gz": f"{BASE}/series/GSE73nnn/GSE73378/matrix/GSE73378-GPL10558_series_matrix.txt.gz",
    # probe annotation (GPL10558 -> gene symbol)
    "GPL10558.annot.gz": f"{BASE}/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz",
}

HDR = {"User-Agent": "Mozilla/5.0 (bioinformatics pipeline)"}

def expected_size(url):
    try:
        r = requests.head(url, headers=HDR, timeout=30, allow_redirects=True)
        return int(r.headers.get("Content-Length", 0))
    except Exception:
        return 0

def download(name, url, max_tries=5):
    path = os.path.join(DATA, name)
    if os.path.exists(path):
        try:
            with gzip.open(path, "rb") as f:  # validate gzip integrity
                while f.peek(1 << 20):
                    f.read(1 << 20)
            print(f"[skip] {name} already complete ({os.path.getsize(path)} bytes)")
            return True
        except Exception:
            os.remove(path)
    exp = expected_size(url)
    for attempt in range(1, max_tries + 1):
        try:
            with requests.get(url, headers=HDR, stream=True, timeout=120) as r:
                r.raise_for_status()
                tmp = path + ".part"
                done = 0
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1 << 16):
                        if chunk:
                            f.write(chunk)
                            done += len(chunk)
            if exp and done < exp:
                raise IOError(f"incomplete {done}/{exp}")
            os.replace(tmp, path)
            print(f"[ok] {name} {done} bytes (expected {exp})")
            return True
        except Exception as e:
            print(f"[try {attempt}] {name}: {repr(e)[:120]}")
            time.sleep(min(3 * attempt, 15))
    return False

if __name__ == "__main__":
    ok = all(download(n, u) for n, u in FILES.items())
    print("ALL_OK" if ok else "SOME_FAILED")
