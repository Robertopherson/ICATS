---
layout: default
title: School and Outreach Videos
---

# School and Outreach Videos

ICATS is intended to be useful both as a research tool and as a teaching aid.
The outreach video demo will be maintained separately from the code repository
so that the main ICATS clone stays small.

The planned demo will provide a simple web interface where a student can choose
two molecules and play a short trajectory-style visualization of the collision.
It is not meant to replace the scientific tutorials; it is a visual introduction
to the idea of molecular scattering.

The video repository should contain only the files needed for the web demo:

```text
scatervid.html
vid_*.mp4
README.md
```

The raw video set is roughly 100 MB. That is technically possible to store in a
GitHub repository, but it will make the repository heavier for every clone.

Recommended publishing plan:

1. Create a separate `ICATS-videos` repository.
2. Keep only the HTML interface and compressed MP4 files.
3. Publish that repository with GitHub Pages.
4. Link to the video page from this manual.

Before uploading, compress one representative video and compare quality. If the
whole set can be reduced substantially, direct GitHub Pages hosting becomes more
reasonable.
