---
layout: default
title: Outreach Collision Player
---

# Outreach Collision Player

ICATS also includes a small outreach player for introducing molecular
collisions visually. Pick two molecules below and the page loads a short
pre-rendered trajectory-style video for that pair.

[Launch the outreach collision player](outreach/index.html)

If this manual is being read through the GitHub repository file viewer, the
embedded window below may not be displayed because GitHub blocks iframes in
Markdown previews. The launch link above opens the same HTML player directly.

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

The current manual bundle includes compressed MP4 files for all 15 combinations
of ammonia, carbon dioxide, water, methane, and nitrogen. These are suitable for
quick web viewing; the raw frame-generation material should stay outside the
main ICATS repository.

## File Layout

The embedded player is stored under:

```text
docs/outreach/index.html
docs/outreach/videos/vid_*.mp4
```

Video names follow the same convention as the player code:

```text
vid_<molecule1>_<molecule2>.mp4
```

The molecule ids are sorted alphabetically in the filename. For example,
selecting water and nitrogen loads `vid_h2o_n2.mp4`.

## Updating The Videos

When replacing a video, keep the manual copy small enough for a normal GitHub
clone. The current files were made with a 960-pixel-wide H.264 encode:

```bash
ffmpeg -i source.mp4 \
  -vf scale=960:-2 \
  -an -c:v libx264 -preset medium -crf 32 -movflags +faststart \
  docs/outreach/videos/vid_h2o_n2.mp4
```

If the outreach material grows into raw videos, frame dumps, scripts, or
multiple versions, it should move to a separate repository and this page can
embed that published page instead.
