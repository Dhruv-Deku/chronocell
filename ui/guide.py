"""
Guide workspace: the whole app in plain language, with analogies, a guided tour, a glossary,
data formats, how to use your own AI key, and what is real versus simulated.
"""

from __future__ import annotations

import streamlit as st

from chronocell import accuracy as ACC
from ui import predict_view
from ui.common import html

ss = st.session_state


def _go(page: str, **extra) -> None:
    ss.workspace = page
    for k, v in extra.items():
        ss[k] = v


def render(version: str) -> None:
    html('<p class="cc-eyebrow" style="margin-top:14px">Guide · ChronoCell-5D ' + version + '</p>'
         '<h1 class="cc-title">The genome, folded — in plain words</h1>'
         '<p class="cc-lead">Every human cell holds about two metres of DNA, packed into a nucleus a hundred times '
         'thinner than a hair. <b>How</b> it is folded decides which genes can be read. In cancer, ageing and some brain '
         'diseases the folding goes wrong. ChronoCell-5D rebuilds that fold in 3D, watches it change over time (the 4th '
         'dimension), and explains what the changes mean (the 5th: interpretation).</p>')
    html('<div class="cc-analogy"><b>The big analogy.</b> Think of a chromosome as a very long necklace stuffed into a '
         'small box. Beads that sit near the lid are easy to reach (genes that are <b>switched on</b>); beads crushed at '
         'the bottom are hard to reach (genes that are <b>silenced</b>). Disease can repack the box.</div>')

    html('<h2 class="cc-h2">A 2-minute tour</h2>')
    html('<ol class="cc-steps">'
         '<li><b>Turn on demo patients</b> in the sidebar (Biological state → <i>Load demo patients</i>). You get a '
         'synthetic healthy, cancer and senescent chr22, clearly labelled as demo data.</li>'
         '<li><b>01 3D structure</b>: rotate the chromosome, read its size and signal in the cards above it.</li>'
         '<li>Pick <b>Disease State / Cancer</b> in the sidebar and watch the numbers change against Healthy.</li>'
         '<li><b>03 Compare</b>: healthy on the left, cancer on the right. Rotate one; the other follows. Dark red = '
         'what moved.</li>'
         '<li><b>04 Drug lab</b>: pick a drug class and drag the dose slider under the fold. See how far it moves the '
         'fold back toward healthy.</li>'
         '<li><b>05 Genes</b>: type a gene (BCR, NF2, CHEK2). Is it open (likely active) or buried (likely silenced)? '
         'Which genes does it touch in 3D?</li>'
         '<li>Scroll to <b>🤖 ChronoAgent</b> under any page, press <i>Analyse</i>, and download the Markdown report or '
         'the PDF dossier.</li></ol>')
    cols = st.columns(5)
    for col, (label, page) in zip(cols, [("Open 3D structure", "3D structure"), ("Open Compare", "Compare"),
                                         ("Open Drug lab", "Drug lab"), ("Open Genes", "Genes"),
                                         ("Open 4D dynamics", "4D dynamics")]):
        col.button(label, key=f"guide_go_{page}", on_click=_go, args=(page,), width="stretch")

    html('<h2 class="cc-h2">What each page does</h2>')
    left, right = st.columns(2, gap="large")
    with left:
        html('<div class="cc-callout"><h4>01 · 3D structure — the map</h4>The chromosome as a tube you can rotate. '
             'Each bead is a short stretch of DNA (10,000 letters on chr22). Colour it by position, activity marks '
             '(<i>Epigenomic Signal Heatmap</i>: blue quiet, magenta busy), neighbourhoods (TADs) or the classic '
             'rainbow. Panels on the right measure it like an engineer measures a coil.</div>')
        html('<div class="cc-callout"><h4>02 · 4D dynamics — the time-lapse</h4>Plays a time course, morphs one '
             'condition into another, or simulates a DNA accident (a deletion, or the swap that makes the '
             '"Philadelphia chromosome" in leukaemia) and lets the fold settle. Like a time-lapse of a room being '
             're-arranged.</div>')
        html('<div class="cc-callout"><h4>03 · Compare — spot the difference</h4>Two states side by side with '
             'linked cameras, the structures overlaid as well as possible, and every bead coloured by how far it '
             'moved. The table says what changed and the list says which genes sit where it changed most.</div>')
    with right:
        html('<div class="cc-callout"><h4>04 · Drug lab — a flight simulator for drugs</h4>Epigenetic drugs work by '
             'loosening or tightening DNA packing in specific places. Pick one; the lab applies <i>where</i> it acts and '
             '<i>which way</i> it pushes, over a dose range, and measures how far the fold returns to healthy. It is a '
             'simulator for ideas, not a prediction for patients.</div>')
        html('<div class="cc-callout"><h4>05 · Genes — who is awake, who is asleep</h4>All 19,386 human genes are '
             'placed on the fold. A gene whose start sits in open, well-marked chromatin is flagged <i>predicted '
             'active</i>; one buried in a crowded clump <i>predicted silenced</i>. Add RNA-seq to check the '
             'prediction against measured activity.</div>')
        html('<div class="cc-callout"><h4>🤖 ChronoAgent — the expert on call</h4>Reads every number on screen and '
             'writes a short report: what the shape means, which drug mechanisms could be tested, which genes are '
             'likely affected. Offline rules work instantly; with your own AI key it answers free-form questions. '
             'Exports a Markdown report and a PDF dossier.</div>')

    html('<h2 class="cc-h2">Reading the numbers</h2>')
    html('<dl class="cc-gloss">'
         '<dt>R<sub>g</sub> (radius of gyration)</dt><dd>How big the ball of DNA is: the average distance of every bead '
         'from the centre. Bigger = more open. <i>Like measuring a ball of wool.</i></dd>'
         '<dt>Max 3D span</dt><dd>The widest distance across the fold, bead to bead.</dd>'
         '<dt>ν (compaction exponent)</dt><dd>How the size grows with length. ≈ 0.33: tightly crumpled like paper '
         'squeezed in a fist (normal for chromatin). ≈ 0.5: a loose coil like cooked spaghetti. Higher: stretched out.</dd>'
         '<dt>Packing fraction</dt><dd>How much of its own space the DNA fills (0.19 for our reference fold).</dd>'
         '<dt>Crowding</dt><dd>How many other beads touch a bead. Crowded spots behave like silenced '
         'heterochromatin.</dd>'
         '<dt>Signal (H3K27ac)</dt><dd>A chemical tag that marks active switches (enhancers) — like sticky notes saying '
         '"read me". Your own track can replace it.</dd>'
         '<dt>TAD</dt><dd>A self-contained neighbourhood of DNA that mostly touches itself: rooms in a house.</dd>'
         '<dt>A / B compartment</dt><dd>The busy city centre (A, active) versus the quiet suburbs (B, inactive).</dd>'
         '<dt>Contact decay γ, P(s) slope</dt><dd>How fast the chance that two pieces of DNA touch drops as they get '
         'further apart along the chain. Slope ≈ −1 for normal crumpled chromatin; steeper = looser.</dd>'
         '<dt>Restoration (%)</dt><dd>In the drug lab: how much closer to the healthy shape the treated fold is.</dd>'
         '<dt>Reference model</dt><dd>A synthetic stand-in shown until real data arrive — a mannequin in the shop '
         'window until the real person walks in.</dd></dl>')

    html('<h2 class="cc-h2">How accurate is it? Two scores, never mixed</h2>')
    html('<div class="cc-analogy"><b>The exam analogy.</b> <i>Contact-map fit</i> is like checking a student against '
         'the homework they copied from: a high mark only shows they copied carefully. <i>Microscopy accuracy</i> is the '
         'real exam: questions they never saw. We report both, side by side, and never add them together.</div>')
    models = (ACC.load_benchmark() or {}).get("models", {})

    def pct(name: str) -> str:
        v = models.get(name, {}).get("overall_percent_of_ceiling")
        return "not measured on this install" if v is None else f"{v:.1f} %"

    pred = models.get("predicted_sequence_ctcf")
    caveat = predict_view.control_caveat(pred) if pred else None
    ev = ACC.interval_evidence() or {}
    rng = lambda r: f"{r[0]:.0f}–{r[1]:.0f} %"                                  # noqa: E731
    hic_txt = (f"stated 90 % ranges held {rng(ev['imaging_recalibrated_90'])} after recalibration with imaging-derived "
               "contacts" if ev.get("imaging_recalibrated_90") else "not yet measured for imaging-derived contacts")
    if ev.get("hic_verdict"):
        hic_txt += (f"; with sequencing Hi-C, {rng(ev['hic_recalibrated_90'])} after recalibration (pre-registered test "
                    f"{'passed' if ev['hic_verdict'] == 'pass' else 'not passed: read the range as model spread only'})")
    elif ev.get("hic_raw_90"):
        hic_txt += f"; with sequencing Hi-C only {rng(ev['hic_raw_90'])} ({ev['hic_source']})"
    html('<dl class="cc-gloss">'
         '<dt>Contact-map fit</dt><dd>How well the 3D model reproduces the contact data it was built from (Spearman '
         '\u03c1, measured live on your window). It shows the fit converged, not that the shape is right.</dd>'
         '<dt>Microscopy accuracy</dt><dd>How well the <b>method</b> predicts distances measured under a microscope in '
         'cells it never saw (chromatin tracing, Bintu et al. 2018; Su et al. 2020). On the held-out test sets the '
         f'windowed <b>population model</b> recovers {pct("ensemble_v3_3")} of the folding pattern the experiment can '
         f'reproduce, the older single-structure model {pct("single_structure_v3_2")}, and the whole-chromosome model '
         f'built from sequencing Hi-C {pct("population_v4")}. It is a property of the method, not a measurement on '
         'your data. Every number is read from the result files; details are in <code>validation/RESULTS.md</code>.</dd>'
         '<dt>Population model</dt><dd>One chromosome folds differently in every cell, like a crowd of people each '
         'standing a little differently. The population model builds the whole crowd (100 trajectories), not one '
         'average person. Find it under 01 3D structure \u2192 03 Model &amp; convergence: the windowed model up to '
         '400 beads, the whole-chromosome model beyond.</dd>'
         '<dt>Predicted from sequence + CTCF (no contact data)</dt><dd>For a window with no Hi-C, the population model '
         'can be built from a <i>prediction</i>: CTCF ChIP-seq peaks of the cell type and the DNA sequence. Like '
         'guessing a room layout from where the doors are. It is labelled "predicted" everywhere. Held-out test: '
         f'{pct("predicted_sequence_ctcf")} of the reproducible pattern, against the numbers above for models built '
         'from contacts: a prior, not a measurement.' + (f' {caveat}' if caveat else '') + '</dd>'
         '<dt>Predicted next to built from contacts</dt><dd>When one window has both models, 03 Model &amp; convergence '
         'shows their distance maps side by side with how well they agree. That is agreement between two models, not '
         'accuracy, and the two are never blended: mixing a prediction with contact data was not tested.</dd>'
         '<dt>Measure and Slicing plane</dt><dd><i>Measure</i> (above the view) gives the distance between any two '
         'beads, and for the population model the typical range across cells. Under it, the app says how often such a '
         f'range held real single-cell distances in held-out tests <i>for your kind of input</i>: {hic_txt}, and untested '
         'for a prediction. <i>Display \u2192 Slicing plane</i> cuts the fold open like slicing a cake, to see '
         'inside.</dd></dl>')

    html('<h2 class="cc-h2">Bring your own data</h2>')
    st.dataframe([
        {"You have": "3D coordinates", "Formats": ".npy (N×3), .pdb, .npz bundle, .xyz, .csv",
         "Where": "sidebar → Add files to a state, Data → Coordinates, or coordinates/<chrom>/"},
        {"You have": "An activity / ChIP / ATAC track", "Formats": ".npy (one value per bead), .bedGraph, .bed, .bigWig*",
         "Where": "sidebar → Add files to a state"},
        {"You have": "Micro-C / Hi-C contacts", "Formats": ".cool, .mcool, .hic*, text table (bin bin count)",
         "Where": "sidebar (per state) or Data → Graph or contact map"},
        {"You have": "RNA-seq expression", "Formats": ".csv / .tsv (gene, value)", "Where": "05 Genes → Add measured expression"},
        {"You have": "Raw sequencing (FASTA + bigWig + mcool)", "Formats": "Colab notebook on a T4 GPU",
         "Where": "colab/ChronoCell5D_Colab.ipynb → unzip output into coordinates/"},
    ], hide_index=True, width="stretch")
    html('<p class="cc-note">* bigWig is read with pyBigWig when installed, otherwise with the built-in reader; .hic is '
         'read with hic-straw when installed, otherwise with the built-in reader. '
         'Files are matched to Healthy / Disease / Senescent by words in their name or folder (healthy, control, tumour, '
         'cancer, senescent…), or by uploading them into a state.</p>')

    html('<h2 class="cc-h2">Using your own AI key</h2>')
    html('<ol class="cc-steps">'
         '<li>Get a free key: <b>Google AI Studio</b> (aistudio.google.com → Get API key; starts with <code>AIza</code>) '
         'or <b>OpenRouter</b> (openrouter.ai → Keys; starts with <code>sk-or-</code>).</li>'
         '<li><b>Quick way:</b> paste it into the sidebar field <i>AI API Key</i>. It stays in this browser tab only.</li>'
         '<li><b>Permanent way:</b> copy <code>.streamlit/secrets.toml.example</code> to <code>.streamlit/secrets.toml</code> '
         'and put your key there (<code>GEMINI_API_KEY = "…"</code>). The app picks it up automatically; the file is '
         'git-ignored.</li>'
         '<li>Press <b>Analyse</b> in the ChronoAgent panel. One press = one request. Without a key, the offline '
         'engine answers instantly.</li></ol>')
    html('<div class="cc-callout"><h4>Keep the key private</h4>The key is sent only to the provider you choose, inside a '
         'request header, and never appears in reports or on screen. If you deploy the app publicly with a key in '
         'secrets, every visitor spends <i>your</i> quota: keep public deployments keyless or protected.</div>')

    html('<h2 class="cc-h2">What is real, and what is simulated</h2>')
    html('<div class="cc-callout"><b>Real:</b> the human genome map (chromosome sizes, bands, gaps, 19,386 genes from '
         'RefSeq), the physics formulas, the reconstruction algorithm, and anything computed from files you provide.'
         '<br><b>Synthetic / simulated:</b> the reference model and the demo patients (labelled everywhere), the 4D '
         'disease scenarios, and the drug lab (a mechanism simulator). Gene "active / silenced" labels are predictions '
         'from structure. ChronoAgent\'s therapy section lists research ideas, never medical advice.</div>')

    _phase_b()
    _quantum()


def _phase_b() -> None:
    """New analyses (Phase B), each with its measured standing read from the validation result files."""
    from chronocell import report_html as RH
    std = {g: RH.gate_standing(g)[0] for g in ("4c", "4d", "6", "7")}
    mouse = ACC.mouse_predictor_evidence()
    html('<h2 class="cc-h2">Analyses for your own maps</h2>')
    left, right = st.columns(2, gap="large")
    with left:
        html('<div class="cc-callout"><h4>Analysis suite · 01 → 06</h4>Finds <b>loops</b> (two distant pieces of DNA held '
             'together, like a bead pinned to another bead), <b>neighbourhoods</b> (TADs) and the active / quiet '
             f'<b>compartments</b> in the contacts of the window, and exports them for IGV and Juicebox. {std["6"]}.</div>')
        html('<div class="cc-callout"><h4>Differential analysis · 03 → third tab</h4>Two conditions (say untreated and '
             'treated), each with its replicate maps: which contacts, loops, boundaries and compartments change. With two '
             f'or more replicates per condition each change gets a false-discovery-controlled test. {std["7"]}.</div>')
    with right:
        html('<div class="cc-callout"><h4>Variant impact engine v2 · 02 → 05</h4>Describe a rearrangement by how its '
             'broken ends are re-joined (or by copy number), even across two chromosomes, and see contacts, '
             'neighbourhoods, genes and enhancer–promoter pairs change, ranked by how likely they are to disrupt gene '
             f'regulation. A <b>mechanism simulator, not validated</b>: {std["4d"]}; the cohesin-loss model it shares '
             f'the ensemble with: {std["4c"]}.</div>')
        mtxt = ("passed on two mouse ES-cell loci (" + "; ".join(f"{v['percent_of_ceiling']:.0f} % of the pattern"
                                                               for v in mouse["loci"].values()) + ")"
                if mouse and mouse["verdict"] == "pass" else "not offered (the mouse test has not passed)")
        html('<div class="cc-callout"><h4>Platform</h4><b>Research mode</b> (sidebar) hides the mechanism simulators and '
             'the rule-based state labels when off. <b>Projects</b> save a session; <b>Jobs</b> run heavy work in the '
             'background, one GPU job at a time. <b>Mouse</b> prediction without contacts: ' + mtxt + '.</div>')



def _quantum() -> None:
    """07 Quantum lab in plain words, with Gate Q's standing (Research mode only)."""
    if not st.session_state.get("research_mode", True):
        return
    from ui import quantum_lab as QL
    html('<h2 class="cc-h2">Quantum lab, in plain words</h2>')
    html('<div class="cc-callout">A quantum computer stores information in <b>qubits</b>, which can be 0 and 1 at the '
         'same time, and is good at searching many possible answers at once. ChronoCell turns some of its questions '
         'into the forms a quantum computer takes and solves them on a <b>simulator</b> running on this computer (not '
         'real quantum hardware), always next to the normal (classical) answer. No speed-up is claimed: at the sizes a '
         'simulator can handle, the classical methods are faster.</div>')
    left, right = st.columns(2, gap="large")
    with left:
        html('<div class="cc-callout"><h4>Words you will see</h4><b>QUBO</b>: a puzzle of yes/no choices with a cost for '
             'each choice and for each pair of choices; every choice becomes one qubit. <b>QAOA</b>: a quantum algorithm '
             'that tunes a circuit until measuring it gives low-cost answers. <b>VQE</b>: the same idea for the energy '
             'of a molecule. <b>Shots</b>: how many times the circuit is run and measured. <b>Noise</b>: today’s '
             'machines make errors on every gate; the noise switch shows what that does.</div>')
        html('<div class="cc-callout"><h4>Where to find it</h4><b>07 Quantum lab</b> has everything in one place, plus '
             'a sizes chart and how to run a circuit on a real quantum computer (download it as OpenQASM). Each '
             'workspace also has a Quantum section: 01 → 07 domain walls and lattice fold; 02 → 06 quantum walk and '
             'variant set; 03 → Quantum similarity; 04 → drug combination and molecule energy; 05 → gene classifier '
             'and gene group.</div>')
    with right:
        r = QL.gate_q()
        if r and "overall" in r:
            items = "".join(f"<li>{QL.PART_LABEL[k]}: <b>{'pass' if v else 'fail'}</b></li>" for k, v in r["overall"].items())
            html(f'<div class="cc-callout"><h4>Measured standing (Gate Q, held-out data)</h4><ul>{items}</ul>The numbers '
                 'are in validation/RESULTS.md.</div>')
        else:
            html('<div class="cc-callout"><h4>Measured standing</h4>Gate Q has not been run.</div>')
        html('<div class="cc-callout"><h4>What it is good for</h4>Learning and showing how quantum algorithms would '
             'tackle chromatin questions, with honest comparisons. It is not a faster or more accurate way to analyse '
             'your data today: use the classical tools in the other workspaces for that.</div>')
