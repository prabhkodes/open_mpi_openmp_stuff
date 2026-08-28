# Interview cheat-sheet — DKRZ RSE, HPC Code Performance Optimisation

> Built on Day 6 from notes/day1-6.md. Keep it to two pages. If you cannot
> say a line out loud in under 30 seconds, cut it.

## The 90-second version

<!-- What you built this week and what it demonstrates. Lead with the
     coupling work -- that is the job. -->

## Numbers I can quote

| Claim | Number | Where it came from |
|---|---|---|
| Halo exchange, blocking vs non-blocking crossover |  | Ex 01 |
| False-sharing penalty |  | Ex 02 |
| Comm/compute overlap efficiency |  | Ex 03 |
| Edge cut, linear vs Hilbert at 64 parts |  | Ex 07 |
| Best nproma, and why |  | Ex 08 |
| Conservation error, conservative vs linear remap |  | Ex 10 |
| Load imbalance, wasted core-hours |  | Ex 16 |
| GPU data-residency speedup |  | Ex 17 |
| Collective vs independent parallel I/O |  | Ex 20 |
| Ice-sheet coupling: sync vs concurrent throughput |  | Ex 23 |

## The three diagrams I can draw in 30 seconds

1. **Coupled timeline, synchronous vs concurrent** (Ex 23) — the most
   important one. Three bars, mark the stall.
2. **Roofline with the stencil on it** (Ex 14) — ridge point, where the
   kernel sits, which way blocking moves it.
3. **Halo schedule inversion** (Ex 06) — "I know what I need; the sender
   does not" and the two collectives that fix it.

## Questions I have for them

<!-- Have four. Good ones for this role:
     - Which model(s) would I be working on, and what is the current
       biggest performance pain point?
     - Is the ice-sheet coupling work YAC-based, and is it new development
       or extending existing infrastructure?
     - How much of the role is GPU porting vs CPU optimisation vs
       coupling infrastructure?
     - Which European systems is the code being adapted for right now? -->

## Things to be honest about

<!-- Where your experience is thin, and what you did about it. Interviewers
     respect "I had not worked with unstructured grids, so I built X and
     learned Y" far more than a bluff they can see through. -->
