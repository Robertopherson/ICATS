# ICATS Working Notes

## Next Session

1. Draft project/funding/authorship statements:
   - who Roberto/Chris is in relation to ICATS and the manuscript,
   - what the human scientific contribution was,
   - what AI-assisted coding/documentation contributed,
   - funding/support acknowledgement language suitable for the repository,
     manual, and possibly the manuscript.

2. Build the outreach/video side:
   - decide whether the video material belongs in a separate GitHub repository,
   - identify the minimal HTML/video assets needed from
     `/home/chris/work/initc/new_videos`,
   - make the manual link to or embed the outreach pages cleanly,
   - avoid copying unnecessarily large raw video/project files into ICATS.

3. Publication-quality Wang-Landau figure:
   - regenerate NH3 + H2O `J`, `L`, and umbrella panels from final scripts/data,
   - archive the input file, `wang.pkl` metadata, and plotting script,
   - use the screenshots currently in the manual only as diagnostic/manual
     figures unless replaced by regenerated publication panels.

4. Revisit Wang-Landau automation:
   - automatically estimate the sampled `Jab` distribution and useful `Jab`
     range before building the WL umbrella,
   - use that estimate to choose the WL range without asking the user to guess
     it manually,
   - extend the WL-corrected region somewhat beyond the important `Jab` mixing
     range,
   - for `J > 1.5 * max(Jab)` or a similar automatically chosen threshold,
     consider switching to a flat/asymptotic correction because `J` mostly
     follows `L`,
   - document and test this carefully so the automatic choice is conservative
     and visible in `wang.pkl` metadata/log output.

5. Finish leading-Wigner naming cleanup:
   - the sampled vibrational radial law is now documented as the
     energy-matched leading-Wigner sampler, not the textbook Husimi function,
   - a focused sampler moment test exists in
     `tests/test_leading_wigner_sampling.py`, covering
     `<0.5*(Q^2+P^2)> = v + 0.5`,
   - later, consider renaming compatibility functions such as `HarmHusimi` and
     `HusimiFuncICDF` behind aliases, but only with a careful deprecation path
     because old generated metadata may reference those names.
