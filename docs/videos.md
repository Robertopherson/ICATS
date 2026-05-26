# School and Outreach Videos

The outreach videos are best kept separate from the main code repository. The
current candidate material is in:

```text
/home/chris/work/initc/new_videos
```

The useful minimal set is:

```text
scatervid.html
vid_*.mp4
README.md
```

The raw video set is roughly 100 MB. That is technically possible to store in a
GitHub repository, but it will make the repository heavier for every clone.

Recommended plan:

1. Create a separate `ICATS-videos` repository.
2. Keep only the HTML interface and compressed MP4 files.
3. Publish that repository with GitHub Pages.
4. Link to the video page from this manual.

Before uploading, compress one representative video and compare quality. If the
whole set can be reduced substantially, direct GitHub Pages hosting becomes more
reasonable.
