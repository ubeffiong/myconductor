# Regulatory-region input

The JSON object has schema `myconductor.regulatory.v1`, a nonempty `version`,
boolean `illustrative`, and `regions`. Each region follows
`RegulatoryRegion` in ../regulatory.py. Coordinates are either null or
`{chrom,start,end,strand}`, with **1-based inclusive** positions in the
explicit reference assembly. Unknown coordinates remain null.

Matching uses an explicit variant-label list first, coordinates second, and
otherwise an annotated promoter/intergenic gene association. The latter is
region membership only. No transcription-start site or expression effect is
inferred. CATALOGUED describes this imported region record, not a WHO
resistance grade. It cannot establish resistance. Candidate records are
report-only. Retractions use a new version and status `retracted`; retain
both versions and use the report manifest/change-impact audit.

The distributed table is empty. Existing graded promoter variants remain in
the main catalogue; this file creates no additional resistance assertions.
