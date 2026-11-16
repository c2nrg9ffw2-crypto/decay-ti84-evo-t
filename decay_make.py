# PC only: builds the data block in decay.py from the IAEA nuclide chart.
# 1. Download the table (ground states of all nuclei):
#    curl -o gs.csv "https://nds.iaea.org/relnsd/v1/data?fields=ground_states&nuclides=all"
# 2. python3 decay_make.py gs.csv
import csv
import math
import sys

YEAR = 31556926                  # seconds, same year as the IAEA table
STABLE_Y = 1e16                  # longer than this = "practically stable"
# IAEA decay name -> one letter (must match MOVE in decay.py)
CODE = {"A": "a", "B-": "b", "EC+B+": "c", "B+": "c", "EC": "e",
        "P": "p", "2P": "q", "N": "n", "2N": "m", "B-N": "d", "B-2N": "t",
        "ECP": "g", "B+P": "g", "B-A": "h", "SF": "f", "ECSF": "f",
        "2B-": "l", "2B+": "l", "2EC": "l", "B+A": "i", "ECA": "i",
        "SF+EC+B+": "f"}
SIDE_MIN = 0.1                   # smaller side branches (in %) are left out
# Same order and k ranges as UNITS in decay.py
UNITS = [("s", 1, -24, -10), ("ns", 1e-9, 0, 2), ("us", 1e-6, 0, 2),
         ("ms", 1e-3, 0, 2), ("s", 1, 0, 1), ("min", 60, 0, 1),
         ("h", 3600, 0, 1), ("d", 86400, 0, 2), ("y", YEAR, 0, 32)]
ABC = ("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
       "abcdefghijklmnopqrstuvwxyz!#%()*+,-./:;<=>?[]_{|}")
assert len(ABC) == 85 and len(set(ABC)) == 85


def pick(t):                     # which unit fits t seconds
    for i, top in enumerate([1e-9, 1e-6, 1e-3, 1, 60, 3600, 86400, YEAR]):
        if t < top:
            return i
    return 8


def round2(v):                   # 2 digits: v = m * 10^(k-1)
    k = math.floor(math.log10(v))
    m = round(v / 10 ** (k - 1))
    if m == 100:
        m, k = 10, k + 1
    return m, k


def code_hl(t):                  # seconds -> number (0 = unknown)
    if not t or t <= 0:
        return 0
    for _ in range(3):           # rounding can push it to the next unit
        u = pick(t)
        m, k = round2(t / UNITS[u][1])
        t2 = m * 10 ** (k - 1) * UNITS[u][1]
        if pick(t2) == u:
            break
        t = t2
    k = max(UNITS[u][2], min(UNITS[u][3], k))
    v = 1
    for i in range(u):
        v += (UNITS[i][3] - UNITS[i][2] + 1) * 90
    return v + (k - UNITS[u][2]) * 90 + (m - 10)


def code_pct(p):                 # side branch % -> one letter of ABC
    if p < 1.5:
        return ABC[max(0, min(9, round(p * 10) - 1))]   # 0.1 .. 1.0 %
    return ABC[min(84, round(p) + 8)]                   # 2 .. 76 %


def num(s):
    try:
        return float(s)
    except ValueError:
        return None


def build(path):
    rows = {}
    sym = {}
    for r in csv.DictReader(open(path)):
        z, n = int(r["z"]), int(r["n"])
        sym[z] = r["symbol"]
        t = num(r["half_life_sec"])
        if r["half_life"] == "STABLE":
            m = "s"
        else:
            m = CODE.get(r["decay_1"], "u")
            if t and t > STABLE_Y * YEAR:
                m = "l"
        side = []
        if m in "abcepqnmdtghi":
            for k in ("2", "3"):
                c, p = CODE.get(r["decay_" + k], ""), num(r["decay_" + k + "_%"])
                if c and c in "abcepqnmdtghif" and p and p >= SIDE_MIN:
                    side.append((c, p))
        rows[(z, n)] = (m, code_hl(t) if m != "s" else 0, side)
    zmax = max(z for z, n in rows)
    nmin, cnt, modes, hl, br = [], [], [], [], []
    for z in range(zmax + 1):
        ns = [n for (zz, n) in rows if zz == z]
        lo, hi = min(ns), max(ns)
        nmin.append(lo)
        cnt.append(hi - lo + 1)
        for n in range(lo, hi + 1):
            m, h, side = rows.get((z, n), ("u", 0, []))
            for c, p in side:
                i = len(modes)
                br.append(ABC[i // 85] + ABC[i % 85] + c + code_pct(p))
            modes.append(m)
            hl.append(ABC[h // 85] + ABC[h % 85])
    return (nmin, cnt, "".join(modes), "".join(hl), "".join(br),
            [sym[z] for z in range(zmax + 1)])


def chunks(s, w=60):
    return "\n".join('    "' + s[i:i + w] + '"' for i in range(0, len(s), w))


def numbers(name, v):
    out, line = [], name + " = ["
    for x in v:
        p = str(x) + ", "
        if len(line) + len(p) > 70:
            out.append(line.rstrip())
            line = "    "
        line += p
    out.append(line.rstrip(", ") + "]")
    return "\n".join(out)


def main():
    nmin, cnt, modes, hl, br, sym = build(sys.argv[1] if len(sys.argv) > 1 else "gs.csv")
    block = "\n".join([
        numbers("NMIN", nmin),
        numbers("CNT", cnt),
        'SYM = ("' + " ".join(sym[:60]) + '"\n    " ' + " ".join(sym[60:]) + '")',
        "MODES = (\n" + chunks(modes) + ")",
        "HL = (\n" + chunks(hl) + ")",
        "# Side branches: 2 letters place, 1 letter decay, 1 letter %",
        "BR = (\n" + chunks(br) + ")",
    ])
    src = open("decay.py").read()
    a = src.index("# --- DATA START")
    a = src.index("\n", a) + 1
    b = src.index("# --- DATA END")
    open("decay.py", "w").write(src[:a] + block + "\n" + src[b:])
    print(len(modes), "nuclei,", len(br) // 4, "side branches,", len(src[:a] + block + src[b:]), "bytes in decay.py")


main()
