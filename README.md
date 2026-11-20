# Decay Chain Solver for the TI-84 Evo-T

Type a nucleus like **U 238**, or its number of **protons** and **neutrons**.
The calculator follows its **radioactive decay** step by step on the nuclide chart,
until it reaches a stable nucleus.

It is written in Python for the **TI-84 Evo-T** graphing calculator.

> **Made with AI:** this program and this README were created with the help of AI (Claude by Anthropic).
> Tested on a real TI-84 Evo-T.
>
> The file is big (about 27 KB). If the calculator runs out of memory, please open an issue.

## What it shows

For U-238 (92 protons, 146 neutrons):

```
U-238  alpha 4.5E9 y
Th-234 beta- 24 d
Pa-234 beta- 6.7 h
U-234  alpha 250000 y
...
Po-210 alpha 140 d
Pb-206 stable

Result: Pb-206
Main path: 14 steps
8x alpha
6x beta minus
```

- **Every step:** the nucleus, how it decays and its **half-life**.
- **The result:** the final nucleus, with Z, N and the mass number A.
- **Step count:** how many alpha and beta steps it took.

## More than one way to decay

Some nuclei can decay in two or more ways.
The chain follows the **most likely** way, and shows the other ways below it:

```
Bi-212 beta- 1 h
   or alpha 36%
```

After the result, each other way gets its own page. It follows that way step by step,
until it reaches its end or joins the main path again:

```
Other way 1 of 1:
Bi-212 alpha 36%
Tl-208 beta- 3.1 min
Pb-208 stable
```

When the paths end in **different** nuclei, a last page lists all possible ends and their chance:

```
All possible ends:
Ca-40  89%
Ar-40  11%
```

Side paths under 0.1% are left out.

## Decay types

| On screen | Meaning | Protons | Neutrons |
|---|---|---|---|
| alpha | alpha decay | −2 | −2 |
| beta- | beta minus | +1 | −1 |
| beta+ | beta plus | −1 | +1 |
| EC | electron capture | −1 | +1 |
| p / 2p | 1 or 2 protons fly out | −1 / −2 | 0 |
| n / 2n | 1 or 2 neutrons fly out | 0 | −1 / −2 |
| b-,n | beta minus, then a neutron | +1 | −2 |
| fission | the nucleus splits in two | – | – |

**~stable** means the half-life is longer than 10^16 years (for example Bi-209).
The program stops there, like most school nuclide charts.

## What you need

- A **TI-84 Evo-T** calculator
- A USB cable to connect it to your computer
- **Google Chrome** or **Microsoft Edge** (Safari and Firefox do not work with the TI tool)
- The file [`decay.py`](decay.py) from this repo

## Setup

1. Download [`decay.py`](decay.py) to your computer.
2. Connect the calculator to your computer with the USB cable.
   Leave the calculator on the **home screen**.
3. Open **https://connectevo.ti.com** in **Chrome** or **Edge**.
4. Connect to your calculator in the page and allow access when the browser asks.
5. Drag `decay.py` onto the page to send it to the calculator.
6. **Choose RAM.** The Python app on the Evo-T only shows files that are in RAM.

> **Already sent it to Archive?** On the **home screen** (not in the Python app) press **2nd → mem → Mem Management**,
> find **DECAY** and press **enter**. The `*` in front of the name disappears, which means it is now in RAM.

## How to use it

1. Open the **Python** app and run **DECAY**.
2. `Nucleus or Z?` – type the element symbol and the mass number, like `U 238`.
   `U-238`, `U238` and `238U` also work.
3. Or type only the protons, like `92` (or only `U`). Then `Neutrons N?` asks for the neutrons, like `146`.
4. Press **any key** to see the next page.
5. Press **enter** on an empty line to quit.

## Good to know

- The data covers **3,466 nuclei**: all ground states on the IAEA nuclide chart.
- Only **ground states** are used. Example: Pa-234 shows 6.7 h, not the 1.2 min excited state.
- Half-lives are rounded to 2 digits.
- The Evo-T **deletes RAM** when it stays turned off for a while.
  Keep the program in **Archive**, and move it back to RAM (see above) before you use it.

## Data source

Decay modes, branch percentages and half-lives come from the
[IAEA LiveChart of Nuclides](https://www-nds.iaea.org/livechart/) (Nuclear Data Section, International Atomic Energy Agency).

To rebuild the data on a computer (normal Python 3):

```
curl -o gs.csv "https://nds.iaea.org/relnsd/v1/data?fields=ground_states&nuclides=all"
python3 decay_make.py gs.csv
```

[`decay_make.py`](decay_make.py) packs the table into short text strings and writes them into `decay.py`.
