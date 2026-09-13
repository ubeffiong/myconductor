"""Build the self-contained HTML report for a benchmark run.

Report rules, all of them about not overstating
-----------------------------------------------
* **Say what the run establishes, first.** The top of the page states whether
  accuracy was measured at all. A concordance-only run says so in its headline,
  because "96% agreement with TB-Profiler" reads like accuracy to anyone
  skimming, and it is not.
* **A missing value is never a zero.** A drug with no evaluable phenotype
  renders as "not measured", not as 0.000.
* **Underpowered is not failure.** A drug below its pre-registered minimum
  isolate count gets its own state, distinct from a drug that was measured and
  missed its bound.
* **Abstention is shown beside accuracy.** Every accuracy row carries the call
  rate, because a predictor that declines every hard isolate posts excellent
  numbers on the rest.
* **Failures stay visible.** Stage tables show failed and skipped isolates with
  their reasons rather than quietly reporting only what worked.
* Offline and single-file: no CDN, no external fonts, no network at view time.
"""
from __future__ import annotations

import csv
import html
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

from . import __version__

#: One definition per term, used twice: the glossary section renders from it,
#: and every first use in the page is annotated from it. "what" is the
#: measurement; "why" is the half a label can never carry.
GLOSSARY: dict[str, dict[str, str]] = {
    "Concordance": {
        "what": "Fraction of isolates where two engines reached the same "
                "verdict, computed only over isolates both engines called.",
        "why": "It measures agreement, not correctness. Two tools reading the "
               "same catalogue will agree with each other whether or not "
               "either is right, so a high figure here is not evidence of "
               "accuracy.",
    },
    "VME": {
        "what": "Very major error: predicted susceptible, phenotypically "
                "resistant, as a fraction of phenotypically resistant "
                "isolates.",
        "why": "The error that reaches the patient. A drug called susceptible "
               "when the organism is resistant is a drug that will be given "
               "and will not work, so this is the tightest bound in the "
               "pre-registration.",
    },
    "ME": {
        "what": "Major error: predicted resistant, phenotypically susceptible, "
                "as a fraction of phenotypically susceptible isolates.",
        "why": "A usable drug withheld. Harmful but recoverable, so its bound "
               "is looser than VME's.",
    },
    "Call rate": {
        "what": "Fraction of evaluable isolates for which a verdict was "
                "reached at all.",
        "why": "Accuracy here is conditional on a call being made. A predictor "
               "that abstains on every difficult isolate would post near "
               "perfect accuracy on the easy remainder, so this number and the "
               "error rates must be read together.",
    },
    "Abstention": {
        "what": "An isolate/drug where the verdict was NOT_ASSESSED, "
                "INDETERMINATE or NO_CALL rather than resistant or "
                "susceptible.",
        "why": "Folding these into 'susceptible' is the specific defect the "
               "interpretation engine was rebuilt to remove, so the scoring "
               "here excludes them from error rates and counts them "
               "separately.",
    },
    "Evaluable": {
        "what": "An isolate/drug with a laboratory phenotype at accepted "
                "quality.",
        "why": "Measuring a predictor against a low-quality phenotype measures "
               "the phenotype. Only phenotypes at the accepted grade enter the "
               "denominators.",
    },
    "Effective n": {
        "what": "Inverse Simpson index over BioProject sizes.",
        "why": "213 isolates drawn from 9 studies, one of which supplies half, "
               "are not 213 independent observations. This is a descriptive "
               "summary of that clustering, not a correction that licenses "
               "treating them as independent.",
    },
    "Underpowered": {
        "what": "Fewer evaluable isolates than the pre-registered minimum for "
                "that drug.",
        "why": "Not having measured something is different from having "
               "measured it and found it wanting, so it gets its own state "
               "rather than being reported as a failure.",
    },
    "Species control": {
        "what": "A non-tuberculous mycobacterium the pipeline is required to "
                "refuse.",
        "why": "A drug call against the MTBC profile for an organism that is "
               "not MTBC is a pipeline failure. These isolates exist to be "
               "declined, and a verdict on one is the result being wrong.",
    },
}


def esc(value) -> str:
    return html.escape("" if value is None else str(value))


def read_tsv(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with open(path, newline="", encoding="utf-8") as handle:
        stream = (line for line in handle
                  if line.strip() and not line.lstrip().startswith("#"))
        return list(csv.DictReader(stream, delimiter="\t"))


def read_json(path: Path) -> dict:
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def pct(value: Optional[str | float], digits: int = 1) -> str:
    if value in (None, ""):
        return '<span class="na">not measured</span>'
    try:
        return f"{float(value) * 100:.{digits}f}%"
    except (TypeError, ValueError):
        return esc(value)


def num(value, digits: int = 3) -> str:
    if value in (None, ""):
        return '<span class="na">not measured</span>'
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return esc(value)


def band(value: Optional[str], good_high: bool = True) -> str:
    """Colour band for a rate. Returns a CSS class, never a silent zero."""
    if value in (None, ""):
        return "b-na"
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "b-na"
    if not good_high:
        v = 1 - v
    if v >= 0.95:
        return "b-good"
    if v >= 0.85:
        return "b-ok"
    if v >= 0.70:
        return "b-warn"
    return "b-bad"


VERDICT_CLASS = {"pass": "v-pass", "fail": "v-fail",
                 "underpowered": "v-under", "not-targeted": "v-na"}


@dataclass
class ReportData:
    results_dir: Path
    cohort: list[dict]
    manifest: dict
    concordance: list[dict]
    accuracy: list[dict]
    controls: list[dict]
    engine_calls: list[dict]
    my_calls: list[dict]
    stage_status: dict[str, list[dict]]

    @property
    def accuracy_measured(self) -> bool:
        return bool(self.accuracy)

    @property
    def n_isolates(self) -> int:
        return len({r["sample_id"] for r in self.my_calls}) or len(self.cohort)


def load(results_dir: Path, samples: Optional[Path]) -> ReportData:
    stage_status = {}
    for stage in ("download", "profile", "interpret"):
        rows = read_tsv(results_dir / f"{stage}_status.tsv")
        if rows:
            stage_status[stage] = rows
    cohort = []
    if samples and Path(samples).is_file():
        cohort = read_tsv(Path(samples))
    return ReportData(
        results_dir=results_dir, cohort=cohort,
        manifest=read_json(results_dir / "run_manifest.json"),
        concordance=read_tsv(results_dir / "concordance.tsv"),
        accuracy=read_tsv(results_dir / "accuracy.tsv"),
        controls=read_tsv(results_dir / "species_control.tsv"),
        engine_calls=read_tsv(results_dir / "engine_calls.tsv"),
        my_calls=read_tsv(results_dir / "myconductor_calls.tsv"),
        stage_status=stage_status)


# -- sections -------------------------------------------------------------
def headline(data: ReportData) -> str:
    if data.accuracy_measured:
        verdicts = [r.get("verdict", "") for r in data.accuracy]
        failed = verdicts.count("fail")
        passed = verdicts.count("pass")
        under = verdicts.count("underpowered")
        klass = "hl-mixed" if failed else "hl-ok"
        body = (f"Accuracy was measured against laboratory phenotypes for "
                f"<strong>{len(data.accuracy)}</strong> drug(s): "
                f"<strong>{passed} met</strong> their pre-registered bounds, "
                f"<strong>{failed} missed</strong> them, and "
                f"<strong>{under}</strong> were underpowered.")
    else:
        klass = "hl-warn"
        body = ("<strong>No accuracy figure was computed in this run.</strong> "
                "No paired laboratory phenotypes were supplied, so the numbers "
                "below measure <em>agreement with TB-Profiler</em> and nothing "
                "more. Agreement is not correctness: two engines reading the "
                "same catalogue agree with each other whether or not either is "
                "right. To measure accuracy, build a phenotyped cohort with "
                "<code>mycobench fetch-phenotypes</code>.")
    return f'<div class="headline {klass}"><p>{body}</p></div>'

def cohort_section(data: ReportData) -> str:
    if not data.cohort:
        return ""
    projects: dict[str, int] = {}
    geos: dict[str, int] = {}
    tiers: dict[str, int] = {}
    controls = 0
    phenotyped = 0
    for row in data.cohort:
        projects[row.get("bioproject", "")] = projects.get(row.get("bioproject", ""), 0) + 1
        geos[row.get("geo_loc_name", "") or "(absent)"] = \
            geos.get(row.get("geo_loc_name", "") or "(absent)", 0) + 1
        tiers[row.get("cohort_tier", "")] = tiers.get(row.get("cohort_tier", ""), 0) + 1
        if (row.get("expected_outcome") or "").strip():
            controls += 1
        if (row.get("phenotype_source") or "").strip():
            phenotyped += 1

    total = len(data.cohort)
    effective = (1.0 / sum((c / total) ** 2 for c in projects.values())
                 if total and projects else 0.0)
    largest = max(projects.items(), key=lambda kv: kv[1]) if projects else ("", 0)
    share = largest[1] / total if total else 0

    warning = ""
    if share > 0.30:
        warning = (f'<p class="warn">One BioProject ({esc(largest[0])}) '
                   f'supplies {share:.0%} of this panel. Its '
                   f'<span class="term" data-term="Effective n">effective n'
                   f'</span> is {effective:.1f}, not {total}. Treat intervals '
                   f'and significance claims accordingly.</p>')

    rows = "".join(
        f"<tr><td>{esc(p)}</td><td class='n'>{c}</td>"
        f"<td class='n'>{c / total:.0%}</td></tr>"
        for p, c in sorted(projects.items(), key=lambda kv: -kv[1]))
    geo_rows = "".join(
        f"<tr><td><code>{esc(g)}</code></td><td class='n'>{c}</td></tr>"
        for g, c in sorted(geos.items(), key=lambda kv: -kv[1]))

    return f"""
<section id="cohort"><h2>Cohort</h2>
<div class="tiles">
  <div class="tile"><span class="k">isolates</span><span class="v">{total}</span></div>
  <div class="tile"><span class="k">BioProjects</span><span class="v">{len(projects)}</span></div>
  <div class="tile"><span class="k"><span class="term" data-term="Effective n">effective n</span></span><span class="v">{effective:.1f}</span></div>
  <div class="tile"><span class="k">with phenotypes</span><span class="v">{phenotyped}</span></div>
  <div class="tile"><span class="k"><span class="term" data-term="Species control">controls</span></span><span class="v">{controls}</span></div>
  <div class="tile"><span class="k">tiers</span><span class="v">{esc(', '.join(f'{t}:{n}' for t, n in sorted(tiers.items())))}</span></div>
</div>
{warning}
<div class="two">
<div><h3>Study clustering</h3><div class="tw"><table>
<thead><tr><th>BioProject</th><th class="n">isolates</th><th class="n">share</th></tr></thead>
<tbody>{rows}</tbody></table></div></div>
<div><h3>Deposited origin</h3><div class="tw"><table>
<thead><tr><th>geo_loc_name</th><th class="n">isolates</th></tr></thead>
<tbody>{geo_rows}</tbody></table></div>
<p class="note">Origin is the submitter's deposited value, confirmed against
the controlled <code>geo_loc_name</code> field rather than a free-text search.
It is not curator-verified provenance.</p></div>
</div></section>"""


def stages_section(data: ReportData) -> str:
    if not data.stage_status:
        return ""
    blocks = []
    for stage, rows in data.stage_status.items():
        counts = {"ok": 0, "failed": 0, "skipped": 0}
        for row in rows:
            counts[row.get("status", "ok")] = counts.get(row.get("status", "ok"), 0) + 1
        problems = [r for r in rows if r.get("status") in ("failed", "skipped")]
        detail = ""
        if problems:
            items = "".join(
                f"<tr><td><code>{esc(r['sample_id'])}</code></td>"
                f"<td><span class='pill s-{esc(r['status'])}'>{esc(r['status'])}</span></td>"
                f"<td>{esc(r.get('reason', ''))}</td></tr>" for r in problems)
            detail = (f"<details><summary>{len(problems)} isolate(s) not "
                      f"completed</summary><div class='tw'><table><thead><tr>"
                      f"<th>isolate</th><th>status</th><th>reason</th></tr>"
                      f"</thead><tbody>{items}</tbody></table></div></details>")
        blocks.append(
            f"<div class='stage'><h3>{esc(stage)}</h3>"
            f"<p><span class='pill s-ok'>{counts['ok']} ok</span> "
            f"<span class='pill s-failed'>{counts['failed']} failed</span> "
            f"<span class='pill s-skipped'>{counts['skipped']} skipped</span></p>"
            f"{detail}</div>")
    return (f"<section id='stages'><h2>Stage completion</h2>"
            f"<p class='note'>A failing isolate is recorded and skipped, never "
            f"silently dropped and never scored as a zero.</p>"
            f"<div class='stages'>{''.join(blocks)}</div></section>")


def concordance_section(data: ReportData) -> str:
    if not data.concordance:
        return ("<section id='concordance'><h2>Concordance</h2>"
                "<p class='na'>No concordance table was produced.</p></section>")
    rows = []
    for row in sorted(data.concordance, key=lambda r: r["drug"]):
        rate = row.get("agreement_rate", "")
        both = int(row.get("agree", 0) or 0) + int(row.get("disagree", 0) or 0)
        bar = ""
        if rate:
            width = max(1, round(float(rate) * 100))
            bar = (f"<div class='bar'><span class='{band(rate)}' "
                   f"style='width:{width}%'></span></div>")
        rows.append(
            f"<tr><td><strong>{esc(row['drug'])}</strong></td>"
            f"<td class='n'>{esc(row.get('compared'))}</td>"
            f"<td class='n'>{both}</td>"
            f"<td class='n'>{esc(row.get('agree'))}</td>"
            f"<td class='n'>{esc(row.get('disagree'))}</td>"
            f"<td class='n'>{pct(rate)}{bar}</td>"
            f"<td class='n'>{esc(row.get('only_myconductor'))}</td>"
            f"<td class='n'>{esc(row.get('only_engine'))}</td>"
            f"<td class='n'>{esc(row.get('neither'))}</td></tr>")
    return f"""
<section id="concordance"><h2><span class="term" data-term="Concordance">Concordance</span> with TB-Profiler</h2>
<p>Agreement is computed only over isolates <em>both</em> engines called. An
engine abstaining is recorded separately and never counted as agreement: two
tools both declining to answer have not agreed about anything.</p>
<div class="tw"><table>
<thead><tr><th>drug</th><th class="n">isolates</th><th class="n">both called</th>
<th class="n">agree</th><th class="n">disagree</th><th class="n">agreement</th>
<th class="n">only Myconductor</th><th class="n">only TB-Profiler</th>
<th class="n"><span class="term" data-term="Abstention">neither</span></th></tr></thead>
<tbody>{''.join(rows)}</tbody></table></div>
<p class="note">Myconductor declines a drug whose loci were not shown callable.
A high "only TB-Profiler" count therefore reflects the coverage gate doing its
job, not a missing feature.</p></section>"""


def accuracy_section(data: ReportData) -> str:
    if not data.accuracy:
        return f"""
<section id="accuracy"><h2>Accuracy</h2>
<div class="headline hl-warn"><p><strong>Not measured in this run.</strong>
Accuracy requires paired laboratory phenotypes, and none were supplied. The
public Nigerian BioSamples carry no DST or MIC metadata of any kind, so no
accuracy figure can be derived from them.</p></div>
<p>To measure accuracy, build the phenotyped cohort:</p>
<pre><code>mycobench fetch-phenotypes --out-dir data/cryptic \\
  --cohort-out cohorts/cryptic-phenotyped.tsv
mycobench run --cohort cryptic-phenotyped \\
  --phenotypes data/cryptic/cryptic_phenotypes.tsv</code></pre>
<p class="note">Nigerian isolates with paired DST would be the single most
valuable addition to this benchmark: they would turn the concordance track into
an accuracy track for African lineages, where published tools are least
evaluated.</p></section>"""

    rows = []
    for row in sorted(data.accuracy, key=lambda r: r["drug"]):
        verdict = row.get("verdict", "")
        klass = VERDICT_CLASS.get(verdict, "v-na")
        rows.append(
            f"<tr><td><strong>{esc(row['drug'])}</strong></td>"
            f"<td><span class='pill {klass}'>{esc(verdict)}</span></td>"
            f"<td class='n'>{esc(row.get('evaluable'))}</td>"
            f"<td class='n'>{pct(row.get('call_rate'))}</td>"
            f"<td class='n'>{esc(row.get('abstained'))}</td>"
            f"<td class='n mono'>{esc(row.get('tp'))}/{esc(row.get('fp'))}/"
            f"{esc(row.get('tn'))}/{esc(row.get('fn'))}</td>"
            f"<td class='n'>{num(row.get('sensitivity'))}</td>"
            f"<td class='n'>{num(row.get('specificity'))}</td>"
            f"<td class='n {band(row.get('vme_rate'), good_high=False)}'>"
            f"{num(row.get('vme_rate'))}</td>"
            f"<td class='n {band(row.get('me_rate'), good_high=False)}'>"
            f"{num(row.get('me_rate'))}</td></tr>")
        reasons = row.get("reasons", "")
        if reasons:
            rows.append(f"<tr class='sub'><td></td><td colspan='9'>"
                        f"{esc(reasons)}</td></tr>")

    registration = (data.manifest.get("pre_registration") or {})
    return f"""
<section id="accuracy"><h2>Accuracy against laboratory phenotype</h2>
<p>Every figure is conditional on a call being made, and is reported beside the
<span class="term" data-term="Call rate">call rate</span>. Neither means
anything alone.</p>
<p class="note">Bounds were fixed before the data were examined:
<code>{esc(registration.get('description', 'pre-registration not recorded'))}</code></p>
<div class="tw"><table>
<thead><tr><th>drug</th><th>verdict</th>
<th class="n"><span class="term" data-term="Evaluable">evaluable</span></th>
<th class="n">call rate</th>
<th class="n"><span class="term" data-term="Abstention">abstained</span></th>
<th class="n">TP/FP/TN/FN</th><th class="n">sens</th><th class="n">spec</th>
<th class="n"><span class="term" data-term="VME">VME</span></th>
<th class="n"><span class="term" data-term="ME">ME</span></th></tr></thead>
<tbody>{''.join(rows)}</tbody></table></div>
<p class="note"><span class="term" data-term="Underpowered">Underpowered</span>
is deliberately distinct from failure. A drug below its minimum evaluable count
makes no claim in either direction.</p></section>"""


def controls_section(data: ReportData) -> str:
    if not data.controls:
        return ""
    failed = [r for r in data.controls if r.get("passed") != "yes"]
    klass = "hl-warn" if failed else "hl-ok"
    verdict = (f"<strong>{len(failed)} of {len(data.controls)} species "
               f"control(s) FAILED</strong>: a drug verdict was produced for a "
               f"non-MTBC genome against the MTBC profile."
               if failed else
               f"All <strong>{len(data.controls)}</strong> species controls were "
               f"correctly refused: no drug verdict was produced for any "
               f"non-MTBC genome.")
    rows = "".join(
        f"<tr><td><code>{esc(r['sample_id'])}</code></td>"
        f"<td><em>{esc(r.get('organism'))}</em></td>"
        f"<td><span class='pill {'v-pass' if r.get('passed') == 'yes' else 'v-fail'}'>"
        f"{'refused' if r.get('passed') == 'yes' else 'FAILED'}</span></td>"
        f"<td>{esc(r.get('detail'))}</td></tr>"
        for r in data.controls)
    return f"""
<section id="controls"><h2><span class="term" data-term="Species control">Species controls</span></h2>
<div class="headline {klass}"><p>{verdict}</p></div>
<div class="tw"><table><thead><tr><th>isolate</th><th>organism</th>
<th>outcome</th><th>detail</th></tr></thead><tbody>{rows}</tbody></table></div>
</section>"""


def glossary_section() -> str:
    items = "".join(
        f"<div class='gl'><h4>{esc(term)}</h4>"
        f"<p><strong>What.</strong> {esc(body['what'])}</p>"
        f"<p><strong>Why it matters.</strong> {esc(body['why'])}</p></div>"
        for term, body in GLOSSARY.items())
    return (f"<section id='keys'><h2>What these terms mean</h2>"
            f"<div class='glossary'>{items}</div></section>")


def boundaries_section(data: ReportData) -> str:
    extra = ""
    if not data.accuracy_measured:
        extra = ("<li>Nothing in this report is an accuracy measurement. No "
                 "laboratory phenotype was supplied.</li>")
    return f"""
<section id="boundaries"><h2>What this report does not establish</h2>
<ul>
{extra}
<li>No component of Myconductor has been clinically evaluated. These numbers
are research output, not evidence of fitness for patient care.</li>
<li>Concordance with TB-Profiler is agreement between two implementations. Both
can be wrong together, and sharing a catalogue makes that correlated rather
than independent.</li>
<li>Isolates from one BioProject share a sampling frame, a laboratory and often
a transmission chain. Interval widths computed as though they were independent
are too narrow.</li>
<li>Deposited metadata is the submitter's claim. Geographic origin, collection
date and isolation source were not curator-verified.</li>
<li>Lineage balance is not controlled. Lineage confounds both resistance
prediction and catalogue transfer, so performance on this panel need not
transfer to another population.</li>
<li>The TB-Profiler results-JSON parser is written from a documented schema and
has not been validated against real output; the collate parser has. Which
parser produced a call is recorded in <code>engine_calls.tsv</code>.</li>
</ul></section>"""


def provenance_section(data: ReportData) -> str:
    manifest = data.manifest
    if not manifest:
        return ""
    images = manifest.get("container_images") or {}
    catalogue = manifest.get("catalogue") or {}
    phenotypes = manifest.get("phenotypes") or {}
    rows = [
        ("mycobench", manifest.get("mycobench_version", "")),
        ("myconductor", manifest.get("myconductor_version", "")),
        ("generated", manifest.get("generated_utc", "")),
        ("cohort", (manifest.get("cohort") or {}).get("path", "")),
        ("catalogue", catalogue.get("catalogue_version")
         or catalogue.get("source", "")),
        ("catalogue commit", catalogue.get("commit", "")),
        ("phenotypes", phenotypes.get("release") or phenotypes.get("note", "")),
        ("pre-registration", (manifest.get("pre_registration") or {}).get("sha256", "")),
    ]
    rows += [(f"image: {name}", digest) for name, digest in images.items()]
    body = "".join(
        f"<tr><td>{esc(k)}</td><td class='mono'>{esc(v) or '<span class=na>not recorded</span>'}</td></tr>"
        for k, v in rows if k)
    warning = ""
    if catalogue.get("warning"):
        warning = f"<p class='warn'>{esc(catalogue['warning'])}</p>"
    return (f"<section id='provenance'><h2>Provenance</h2>{warning}"
            f"<div class='tw'><table><tbody>{body}</tbody></table></div>"
            f"</section>")


def artifacts_section(data: ReportData) -> str:
    names = ["run_manifest.json", "concordance.tsv", "accuracy.tsv",
             "species_control.tsv", "engine_calls.tsv",
             "myconductor_calls.tsv", "download_status.tsv",
             "profile_status.tsv", "interpret_status.tsv"]
    rows = []
    for name in names:
        path = data.results_dir / name
        if path.is_file():
            size = path.stat().st_size
            rows.append(f"<tr><td><a href='{esc(name)}'>{esc(name)}</a></td>"
                        f"<td class='n'>{size:,} bytes</td></tr>")
    if not rows:
        return ""
    return (f"<section id='artifacts'><h2>Underlying data</h2>"
            f"<p>Every number above comes from these files. They sit beside "
            f"this report so a reader can recompute rather than trust it.</p>"
            f"<div class='tw'><table><thead><tr><th>file</th>"
            f"<th class='n'>size</th></tr></thead><tbody>"
            f"{''.join(rows)}</tbody></table></div></section>")


STYLE = """
:root{color-scheme:light;--bg:#f6f7f9;--surface:#fff;--sunk:#eef0f4;
--ink:#171b24;--ink2:#414a5a;--ink3:#6d7889;--rule:#d8dde5;--rule2:#e7eaf0;
--accent:#2c5f92;--accent-soft:#e4edf6;--good:#2b6b52;--good-soft:#e3f0ea;
--warn:#8a6510;--warn-soft:#faf1dc;--bad:#a8235f;--bad-soft:#fae8f0;
--mono:"IBM Plex Mono",ui-monospace,Consolas,monospace;
--sans:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){color-scheme:dark;
--bg:#10131a;--surface:#171b24;--sunk:#1e232e;--ink:#e8ebf1;--ink2:#b3bcca;
--ink3:#808b9c;--rule:#2c333f;--rule2:#222834;--accent:#7fb2e0;
--accent-soft:#17293b;--good:#72c3a0;--good-soft:#142a21;--warn:#d8ab4e;
--warn-soft:#2e2413;--bad:#ef7fae;--bad-soft:#341623}}
:root[data-theme=dark]{color-scheme:dark;--bg:#10131a;--surface:#171b24;
--sunk:#1e232e;--ink:#e8ebf1;--ink2:#b3bcca;--ink3:#808b9c;--rule:#2c333f;
--rule2:#222834;--accent:#7fb2e0;--accent-soft:#17293b;--good:#72c3a0;
--good-soft:#142a21;--warn:#d8ab4e;--warn-soft:#2e2413;--bad:#ef7fae;
--bad-soft:#341623}
*{box-sizing:border-box}
body{background:var(--bg);color:var(--ink);font-family:var(--sans);
font-size:15px;line-height:1.6;margin:0;-webkit-font-smoothing:antialiased}
.shell{max-width:1180px;margin:0 auto;padding:0 22px 80px}
header.mast{padding:44px 0 22px;border-bottom:2px solid var(--ink);margin-bottom:28px}
.eyebrow{font-family:var(--mono);font-size:11px;letter-spacing:.15em;
text-transform:uppercase;color:var(--ink3)}
h1{font-size:clamp(1.8rem,4.5vw,2.7rem);line-height:1.08;margin:14px 0 0;
letter-spacing:-.02em}
h2{font-size:1rem;letter-spacing:.08em;text-transform:uppercase;color:var(--ink3);
margin:0 0 16px;padding-bottom:8px;border-bottom:1px solid var(--rule)}
h3{font-size:1.1rem;margin:0 0 10px}
h4{font-size:.98rem;margin:0 0 6px}
section{margin:0 0 44px}
p{max-width:74ch}
code{font-family:var(--mono);font-size:.86em;background:var(--sunk);
padding:1px 5px;border-radius:2px;word-break:break-word}
pre{background:var(--sunk);padding:12px 14px;overflow-x:auto;border-radius:3px}
pre code{background:none;padding:0}
a{color:var(--accent)}
.note{color:var(--ink3);font-size:.92rem}
.warn{background:var(--warn-soft);border-left:3px solid var(--warn);
padding:10px 14px;color:var(--ink)}
.na{color:var(--ink3);font-style:italic}
.headline{padding:16px 20px;margin:0 0 22px;border-left:3px solid var(--accent);
background:var(--accent-soft)}
.headline p{margin:0;max-width:80ch}
.hl-ok{background:var(--good-soft);border-left-color:var(--good)}
.hl-warn{background:var(--warn-soft);border-left-color:var(--warn)}
.hl-mixed{background:var(--bad-soft);border-left-color:var(--bad)}
.tiles{display:grid;gap:1px;background:var(--rule);border:1px solid var(--rule);
grid-template-columns:repeat(auto-fit,minmax(140px,1fr));margin:0 0 18px}
.tile{background:var(--surface);padding:12px 14px;display:flex;
flex-direction:column;gap:4px}
.tile .k{font-family:var(--mono);font-size:10px;letter-spacing:.12em;
text-transform:uppercase;color:var(--ink3)}
.tile .v{font-size:1.4rem;font-variant-numeric:tabular-nums}
.two{display:grid;gap:26px}
@media(min-width:840px){.two{grid-template-columns:1fr 1fr}}
.tw{overflow-x:auto;border:1px solid var(--rule);background:var(--surface)}
table{border-collapse:collapse;width:100%;font-size:.9rem}
th,td{text-align:left;padding:8px 12px;border-bottom:1px solid var(--rule2);
vertical-align:top}
th{font-family:var(--mono);font-size:10px;letter-spacing:.09em;
text-transform:uppercase;color:var(--ink3);background:var(--sunk);white-space:nowrap}
td.n{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
td.mono,.mono{font-family:var(--mono);font-size:.84em}
tr.sub td{color:var(--ink3);font-size:.86rem;padding-top:0;border-bottom:1px solid var(--rule2)}
tr:last-child td{border-bottom:0}
.pill{display:inline-block;font-family:var(--mono);font-size:10px;font-weight:600;
letter-spacing:.05em;padding:2px 7px;white-space:nowrap}
.v-pass,.s-ok{background:var(--good-soft);color:var(--good)}
.v-fail,.s-failed{background:var(--bad-soft);color:var(--bad)}
.v-under,.s-skipped{background:var(--warn-soft);color:var(--warn)}
.v-na{background:var(--sunk);color:var(--ink3)}
.bar{height:4px;background:var(--sunk);margin-top:4px;min-width:56px}
.bar span{display:block;height:100%}
.b-good{background:var(--good)}.b-ok{background:var(--accent)}
.b-warn{background:var(--warn)}.b-bad{background:var(--bad)}
.b-na{background:var(--ink3)}
td.b-good{color:var(--good)}td.b-ok{color:var(--accent)}
td.b-warn{color:var(--warn)}td.b-bad{color:var(--bad);font-weight:600}
.stages{display:grid;gap:16px}
@media(min-width:760px){.stages{grid-template-columns:repeat(3,1fr)}}
.stage{background:var(--surface);border:1px solid var(--rule);padding:14px 16px}
.stage p{margin:0 0 8px;display:flex;gap:6px;flex-wrap:wrap}
details{margin-top:8px}
summary{cursor:pointer;font-size:.88rem;color:var(--accent)}
.glossary{display:grid;gap:14px}
@media(min-width:760px){.glossary{grid-template-columns:1fr 1fr}}
.gl{background:var(--surface);border:1px solid var(--rule);padding:14px 16px}
.gl p{margin:0 0 6px;font-size:.92rem;max-width:none}
.term{border-bottom:1px dotted var(--ink3);cursor:help}
ul{max-width:78ch}li{margin-bottom:6px}
footer{margin-top:52px;padding-top:22px;border-top:2px solid var(--ink);
font-family:var(--mono);font-size:11px;line-height:1.7;color:var(--ink3)}
"""

SCRIPT = """
(function(){
  var g = JSON.parse(document.getElementById('gl-data').textContent);
  document.querySelectorAll('.term').forEach(function(el){
    var t = el.dataset.term, d = g[t];
    if (d) { el.title = t + ' \\u2014 ' + d.what + '  Why: ' + d.why; }
  });
})();
"""


def build_report(results_dir: str | Path, samples: Optional[str | Path] = None,
                 out: str | Path = "mycobench_report.html") -> Path:
    results_dir = Path(results_dir)
    data = load(results_dir, Path(samples) if samples else None)

    cohort_name = "unknown cohort"
    if samples:
        cohort_name = Path(samples).stem
    manifest = data.manifest

    body = "".join([
        headline(data),
        cohort_section(data),
        stages_section(data),
        concordance_section(data),
        accuracy_section(data),
        controls_section(data),
        glossary_section(),
        boundaries_section(data),
        provenance_section(data),
        artifacts_section(data),
    ])

    page = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>mycobench report — {esc(cohort_name)}</title>
<style>{STYLE}</style></head><body><div class="shell">
<header class="mast">
  <div class="eyebrow">mycobench {esc(__version__)} &middot;
    myconductor {esc(manifest.get('myconductor_version', '?'))} &middot;
    {esc(manifest.get('generated_utc', ''))}</div>
  <h1>{esc(cohort_name)}</h1>
  <p>Validation run over {data.n_isolates} isolate(s). Research output from a
  scaffold that has not been clinically evaluated.</p>
</header>
{body}
<footer>
Generated by mycobench {esc(__version__)}. Research and decision-support
scaffold — NOT a clinical device. No component has been clinically evaluated,
and no number here should inform the treatment of a patient.<br>
Underlying tables sit beside this file; every figure can be recomputed from
them.
</footer>
</div>
<script id="gl-data" type="application/json">{json.dumps(GLOSSARY)}</script>
<script>{SCRIPT}</script>
</body></html>
"""
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page, encoding="utf-8")
    return out
