---
layout: default
title: Outreach Collision Player
---

# Outreach Collision Player

ICATS also includes a small outreach player for introducing molecular
collisions visually. Pick two molecules below and the page loads a short
pre-rendered trajectory-style video for that pair.

[Launch the outreach collision player on the published manual site](https://robertopherson.github.io/ICATS/outreach/)

If that link gives a 404, GitHub Pages still needs to be enabled for the
repository: use `Settings -> Pages`, choose the `main` branch, and select the
`/docs` folder. Clicking `docs/outreach/index.html` inside the normal GitHub
repository view will show the HTML source code, not the running player.

<iframe src="outreach/index.html"
        title="ICATS molecular collision outreach player"
        loading="lazy"
        width="100%"
        height="790"
        style="border:1px solid #d0d7de; border-radius:8px; background:#0b0f14;">
  <p><a href="outreach/index.html">Open the outreach collision player.</a></p>
</iframe>

## What It Shows

The player is deliberately simple. It is a visual companion to the manual,
rather than a scientific diagnostic. The molecule buttons select one of the
available pair videos, including self-collisions such as water-water or
nitrogen-nitrogen.

The videos are included in this repository. The bundled set contains compressed
MP4 files for all 15 combinations of ammonia, carbon dioxide, water, methane,
and nitrogen. The raw frame-generation material is not included because it is
much larger and is not needed for the web player.
