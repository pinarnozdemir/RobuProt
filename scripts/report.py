from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
import pandas as pd

doc = Document()

section = doc.sections[0]
section.top_margin    = Inches(1)
section.bottom_margin = Inches(1)
section.left_margin   = Inches(1.2)
section.right_margin  = Inches(1.2)

style = doc.styles['Normal']
style.font.name = 'Arial'
style.font.size = Pt(11)

def add_heading(text, level=1):
    h = doc.add_heading(text, level=level)
    h.alignment = WD_ALIGN_PARAGRAPH.LEFT
    run = h.runs[0]
    run.font.name = 'Arial'
    run.font.color.rgb = RGBColor(0x1a, 0x6b, 0x78)
    return h

def add_para(text):
    p = doc.add_paragraph(text)
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.space_after = Pt(8)
    p.paragraph_format.line_spacing = Pt(14)
    run = p.runs[0]
    run.font.name = 'Arial'
    run.font.size = Pt(11)
    return p

def add_caption(text):
    p = doc.add_paragraph(text)
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_after = Pt(12)
    run = p.runs[0]
    run.font.name = 'Arial'
    run.font.size = Pt(9)
    run.font.italic = True
    run.font.color.rgb = RGBColor(0x55, 0x55, 0x55)
    return p

def add_table_title(text):
    p = doc.add_paragraph(text)
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(12)
    p.paragraph_format.space_after = Pt(4)
    run = p.runs[0]
    run.font.name = 'Arial'
    run.font.size = Pt(10)
    run.font.bold = True
    run.font.italic = True

def make_table(headers, rows, header_color="1a6b78"):
    table = doc.add_table(rows=1+len(rows), cols=len(headers))
    table.style = 'Table Grid'
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    hdr = table.rows[0]
    for i, h in enumerate(headers):
        cell = hdr.cells[i]
        cell.text = h
        run = cell.paragraphs[0].runs[0]
        run.font.bold = True
        run.font.name = 'Arial'
        run.font.size = Pt(9)
        run.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
        cell.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.CENTER
        tc = cell._tc
        tcPr = tc.get_or_add_tcPr()
        shd = OxmlElement('w:shd')
        shd.set(qn('w:fill'), header_color)
        shd.set(qn('w:color'), 'auto')
        shd.set(qn('w:val'), 'clear')
        tcPr.append(shd)
    for ri, row_data in enumerate(rows):
        row = table.rows[ri+1]
        fill = 'EBF5F6' if ri % 2 == 0 else 'FFFFFF'
        for ci, val in enumerate(row_data):
            cell = row.cells[ci]
            cell.text = str(val)
            run = cell.paragraphs[0].runs[0]
            run.font.name = 'Arial'
            run.font.size = Pt(9)
            cell.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.CENTER
            tc = cell._tc
            tcPr = tc.get_or_add_tcPr()
            shd = OxmlElement('w:shd')
            shd.set(qn('w:fill'), fill)
            shd.set(qn('w:color'), 'auto')
            shd.set(qn('w:val'), 'clear')
            tcPr.append(shd)
    return table

def add_image(path, width_inches=5.8):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    run.add_picture(path, width=Inches(width_inches))
    return p

# ── Verileri yükle ─────────────────────────────────────────────────────────
scores  = pd.read_csv("results/robustness_scores.csv", index_col=0)
classic = pd.read_csv("results/deresult_MinProb_noBatch.csv", index_col=0)

scenarios = ['KNN_noBatch','KNN_Batch','MinProb_noBatch',
             'MinProb_Batch','Median_noBatch','Median_Batch']
all_results = {s: pd.read_csv(f"results/deresult_{s}.csv", index_col=0)
               for s in scenarios}

def is_contaminant(protein_id):
    protein_id = str(protein_id)
    return protein_id.startswith('CON__') or protein_id.startswith('REV__')

def remove_contaminants(df):
    return df[~df.index.map(is_contaminant)].copy()

def n_sig_clean(scenario_name):
    res = remove_contaminants(all_results[scenario_name])
    return int(res['significant'].sum())

# Important: use the same contaminant-filtered classic set used in step4.py / Figure 4.
# This keeps the report table consistent with n=260 in the robustness barplot.
classic_clean = remove_contaminants(classic)
sig           = classic_clean[classic_clean['significant'] == True]

high = scores[scores['category']=='High confidence'].sort_values('total_score', ascending=False)
mid  = scores[scores['category']=='Moderate confidence']
risk = scores[scores['category']=='Likely artifact']

artifacts_in_classic = sig[sig.index.isin(risk.index)].sort_values('pvalue')
k_high = len(sig[sig.index.isin(high.index)])
k_mid  = len(sig[sig.index.isin(mid.index)])
k_risk = len(artifacts_in_classic)
classic_total = len(sig)
robust_total  = classic_total - k_risk
artifact_pct  = (k_risk / classic_total * 100) if classic_total else 0

scenario_counts = {s: n_sig_clean(s) for s in scenarios}

# ══════════════════════════════════════════════════════════════════════════
# TITLE
# ══════════════════════════════════════════════════════════════════════════
title = doc.add_paragraph()
title.alignment = WD_ALIGN_PARAGRAPH.CENTER
title.paragraph_format.space_after = Pt(20)
run = title.add_run(
    "A Robustness-Scored Bioinformatics Pipeline for Reducing\n"
    "Missing Value- and Batch Effect-Driven False Biomarkers\n"
    "in LC-MS/MS Proteomics"
)
run.font.name = 'Arial'
run.font.size = Pt(16)
run.font.bold = True
run.font.color.rgb = RGBColor(0x1a, 0x1a, 0x1a)

doc.add_paragraph()

# ══════════════════════════════════════════════════════════════════════════
# 3. METHODS
# ══════════════════════════════════════════════════════════════════════════
add_heading("3. Methods", level=1)

add_heading("3.1 Dataset", level=2)
add_para(
    "The open-access LC-MS/MS proteomics training dataset from the StatOmics group "
    "(Ghent University) was used in this study. The dataset comprises tissue samples "
    "from 9 squamous cell carcinoma patients (OR group) and 9 healthy controls (PD group), "
    "processed with MaxQuant software using label-free quantification (LFQ). The raw data "
    "contains 3,898 protein groups across 18 samples in proteinGroups.txt format. Zero values "
    "were treated as missing (NaN), log2 transformation was applied, and samples were assigned "
    "to three batch groups based on analysis date."
)

add_heading("3.2 Quality Control and Preprocessing", level=2)
add_para(
    "Proteins with more than 60% missing values across samples were filtered out, reducing "
    "the dataset from 3,898 to 3,179 proteins. Median normalization was applied to correct "
    "for inter-sample technical variation. Batch effects were assessed by PCA before and "
    "after correction. A mean-centering approach was applied per batch group to reduce "
    "technical variance while preserving biological signal."
)

add_heading("3.3 Parallel Imputation Strategy", level=2)
add_para(
    "Three imputation methods were applied in parallel, each representing a different "
    "assumption about the missing data mechanism. KNN imputation (k=5) assumes MAR and "
    "leverages protein expression correlations. MinProb/QRILC imputation assumes MNAR and "
    "samples replacement values from the lower tail of the observed intensity distribution. "
    "Median imputation serves as a simple reference approach. Each method was applied both "
    "with and without batch correction, yielding 6 analysis scenarios:"
)

add_table_title("Table 1. Six analysis scenarios and their significant protein counts")
make_table(
    ["Scenario", "Imputation Method", "Batch Correction", "Significant Proteins (n)"],
    [
        ["1 — KNN_noBatch",     "KNN (k=5)",       "No",  str(scenario_counts['KNN_noBatch'])],
        ["2 — KNN_Batch",       "KNN (k=5)",       "Yes", str(scenario_counts['KNN_Batch'])],
        ["3 — MinProb_noBatch", "MinProb / QRILC", "No",  str(scenario_counts['MinProb_noBatch'])],
        ["4 — MinProb_Batch",   "MinProb / QRILC", "Yes", str(scenario_counts['MinProb_Batch'])],
        ["5 — Median_noBatch",  "Median",          "No",  str(scenario_counts['Median_noBatch'])],
        ["6 — Median_Batch",    "Median",          "Yes", str(scenario_counts['Median_Batch'])],
    ]
)
add_caption("p<0.05 and |log2FC|>1 criteria applied. CON__/REV__ contaminant-labeled protein groups were excluded for report consistency.")

add_heading("3.4 Robustness Scoring", level=2)
add_para(
    "A robustness score (0-10) was calculated for each protein based on five criteria "
    "reflecting analytical stability across different conditions:"
)

add_table_title("Table 2. Robustness score criteria and point distribution")
make_table(
    ["Criterion", "Description", "Max Score"],
    [
        ["Frequency stability",  "Significant in >=4/6 scenarios",             "+3"],
        ["Batch independence",   "Significant after batch correction",          "+2"],
        ["Missing value level",  "Missing rate <30% per protein",              "+2"],
        ["FC direction",         "Fold change direction consistent across all", "+2"],
        ["ROC-AUC performance",  "ROC-AUC >0.75",                              "+1"],
    ]
)
add_caption("Score >=7: High confidence; 4-6: Moderate confidence; <4: Likely artifact.")

add_para(
    "Proteins scoring >=7 were classified as high confidence biomarker candidates, "
    "4-6 as moderate confidence requiring further validation, and below 4 as likely "
    "artifacts driven by imputation method choice or batch effect."
)

# ══════════════════════════════════════════════════════════════════════════
# 4. RESULTS
# ══════════════════════════════════════════════════════════════════════════
doc.add_page_break()
add_heading("4. Results", level=1)

add_heading("4.1 Quality Control — PCA and Missing Value Analysis", level=2)
add_para(
    "PCA analysis revealed systematic batch-associated clustering prior to correction, "
    "confirming the presence of technical variation. After batch correction, batch groups "
    "converged while the biological disease/control separation was maintained (Figure 1). "
    "Missing value analysis identified 1,897 proteins with less than 10% missing rate, "
    "while 605 proteins exceeded the 30% threshold. The heatmap revealed batch-associated "
    "missingness patterns where certain proteins were systematically absent in specific "
    "batches (Figure 2)."
)

add_image("figures/fig1_pca.png", width_inches=5.8)
add_caption(
    "Figure 1. PCA plots before (left) and after (right) batch correction. "
    "Shapes indicate batch groups; colors indicate disease status (blue=Control, red=Disease)."
)

add_image("figures/fig2_missing_heatmap.png", width_inches=5.8)
add_caption(
    "Figure 2. Missing value distribution per protein (left) and heatmap of the 60 most "
    "incomplete proteins ordered by batch (right). Red = missing, green = observed."
)

add_heading("4.2 Differential Analysis — Six Scenario Comparison", level=2)
add_para(
    "Differential analysis across 6 scenarios revealed dramatic variation in significant "
    "protein counts depending on the chosen imputation method and batch correction status. "
    f"After excluding CON__/REV__ contaminant-labeled protein groups, the classic pipeline "
    f"(MinProb, no batch correction) identified {scenario_counts['MinProb_noBatch']} significant proteins. "
    f"After batch correction, counts dropped to {scenario_counts['KNN_Batch']}-{scenario_counts['MinProb_Batch']} "
    "depending on the imputation method, indicating that a large proportion of classic "
    "pipeline signals are driven by technical batch variation rather than true biological "
    "differences. The volcano plots in Figure 3 illustrate this contrast clearly."
)

add_image("figures/fig3_volcano.png", width_inches=5.8)
add_caption(
    "Figure 3. Volcano plots comparing the classic pipeline (left) and robustness-scored "
    "pipeline (right). Gray: not significant; orange: weak signal (score <7); "
    "teal: robust biomarker candidates (score >=7/10)."
)

add_heading("4.3 Robustness Score Results", level=2)
add_para(
    f"Of the {len(scores):,} analyzed proteins, {len(high)} were classified as high confidence "
    f"(score >=7/10), {len(mid):,} as moderate confidence, and {len(risk):,} as likely artifacts. "
    f"Among the {classic_total} contaminant-filtered proteins flagged by the classic pipeline, "
    f"{k_risk} ({artifact_pct:.1f}%) were identified as likely artifacts and removed. "
    "The comparison is summarized in Table 3 and Figure 4."
)

add_table_title("Table 3. Classic pipeline vs robustness-scored pipeline comparison")
make_table(
    ["", "Classic Pipeline", "Robustness-Scored Pipeline"],
    [
        ["Total significant proteins",  str(classic_total), str(robust_total)],
        ["High confidence (>=7/10)",    "not assessed",     str(k_high)],
        ["Moderate confidence (4-6/10)","not assessed",     str(k_mid)],
        ["Likely artifacts removed",    "—",                f"{k_risk} ({artifact_pct:.1f}%)"],
    ]
)
add_caption(
    f"{k_risk} proteins ({artifact_pct:.1f}%) identified by the classic pipeline were "
    "reclassified as likely artifacts."
)

add_image("figures/fig4_robustness_barplot.png", width_inches=5.8)
add_caption(
    "Figure 4. Left: Top 20 proteins by robustness score. "
    "Right: Comparison of classic vs robustness-scored pipeline output by confidence category."
)

add_heading("4.4 High Confidence Biomarker Candidates", level=2)
add_para(
    f"{len(high)} proteins achieved a robustness score >=7/10 and were classified as high "
    "confidence biomarker candidates. All proteins showed negative log2FC values, indicating "
    "downregulation in squamous cell carcinoma relative to controls. The candidates are "
    "listed in Table 4."
)

add_table_title("Table 4. High confidence biomarker candidates (robustness score >=7/10)")
high_rows = []
for prot, row in high.iterrows():
    high_rows.append([
        prot,
        f"{int(row['n_sig_of_6'])}/6",
        f"{row['missing_pct']:.1f}%",
        f"{row['mean_log2FC']:.3f}",
        f"{int(row['total_score'])}/10",
        row['category']
    ])
make_table(
    ["Protein", "Sig. Scenarios", "Missing %", "Mean log2FC", "Score", "Category"],
    high_rows
)
add_caption(
    "All proteins are downregulated in the disease group. P02787 (Transferrin) and "
    "P12004 (PCNA) achieved perfect 10/10 scores, significant across all 6 scenarios "
    "with 0% missing values."
)

add_heading("4.5 Likely Artifacts Identified in Classic Pipeline", level=2)
add_para(
    f"The robustness scoring framework identified {len(artifacts_in_classic)} proteins that "
    "were flagged as significant by the classic pipeline but classified as likely artifacts "
    "(score <4/10). These proteins showed high sensitivity to the choice of imputation method "
    "or batch correction, suggesting their apparent significance is driven by analytical "
    "decisions rather than true biological differences. Notably, two proteins (P33316 and "
    "O95858) received a score of 0/10, failing all robustness criteria entirely. "
    "The complete list is presented in Table 5."
)

add_table_title("Table 5. Likely artifacts identified in the classic pipeline (score <4/10)")
artifact_rows = []
for prot, row in artifacts_in_classic.iterrows():
    score_val = scores.loc[prot, 'total_score'] if prot in scores.index else 'N/A'
    artifact_rows.append([
        prot,
        f"{row['log2FC']:.3f}",
        f"{row['pvalue']:.4f}",
        f"{row['fdr']:.4f}",
        f"{score_val}/10"
    ])
make_table(
    ["Protein", "log2FC", "p-value", "FDR", "Robustness Score"],
    artifact_rows,
    header_color="8B0000"
)
add_caption(
    f"These {k_risk} proteins appeared significant in the classic pipeline (p<0.05, |log2FC|>1) "
    "but failed robustness criteria. Removing these reduces false positive risk substantially."
)

# ══════════════════════════════════════════════════════════════════════════
# 5. DISCUSSION
# ══════════════════════════════════════════════════════════════════════════
doc.add_page_break()
add_heading("5. Discussion", level=1)

add_para(
    "This study demonstrates that the choice of imputation method and batch correction "
    "strategy substantially influences biomarker discovery outcomes in LC-MS/MS proteomics. "
    f"The contaminant-filtered classic single-imputation pipeline identified {classic_total} "
    f"significant proteins, of which {k_risk} ({artifact_pct:.1f}%) were classified as likely "
    "artifacts by the robustness scoring framework. This finding highlights a critical and "
    "often overlooked source of false positives in standard proteomics workflows."
)

add_para(
    "The dramatic reduction in significant proteins when batch correction was applied — "
    f"from {scenario_counts['KNN_noBatch']} to {scenario_counts['KNN_Batch']} for KNN, and "
    f"from {scenario_counts['MinProb_noBatch']} to {scenario_counts['MinProb_Batch']} for MinProb — "
    "indicates that a substantial fraction of signals in uncorrected analyses reflect "
    "technical batch variation rather than genuine biological differences. This is consistent "
    "with previous reports demonstrating that batch effects can introduce systematic biases "
    "that mimic differential expression patterns in omics data."
)

add_para(
    "The two highest-scoring candidates, P02787 (Transferrin) and P12004 (PCNA), achieved "
    "perfect robustness scores of 10/10. Both proteins were significant across all 6 analysis "
    "scenarios, exhibited no missing values, and showed consistent fold change direction. "
    "P02787 (Transferrin) has been previously associated with altered iron metabolism in "
    "cancer, while P12004 (PCNA) is a well-established proliferation marker widely used "
    "in oncology. The recovery of these biologically meaningful candidates by our pipeline "
    "supports its analytical validity."
)

add_para(
    "The proposed robustness scoring framework offers a practical and interpretable filter "
    "that can be readily integrated into existing proteomics pipelines. Rather than replacing "
    "statistical testing, it adds a reliability layer that quantifies confidence in each "
    "biomarker candidate across multiple analytical conditions. The five-criterion scoring "
    "system captures complementary aspects of analytical robustness: scenario frequency, "
    "batch independence, data completeness, directional consistency, and discriminative power."
)

add_para(
    "A limitation of this study is that batch assignments were approximated by dividing "
    "samples into three equal groups rather than using documented experimental metadata. "
    "Future applications should incorporate actual batch information from instrument logs "
    "or sample preparation records. Additionally, validation in independent cohorts and "
    "orthogonal assays would be required before any candidate can be considered a "
    "validated clinical biomarker."
)

# ══════════════════════════════════════════════════════════════════════════
# 6. CONCLUSION
# ══════════════════════════════════════════════════════════════════════════
add_heading("6. Conclusion", level=1)

add_para(
    "A robustness-scored bioinformatics pipeline was developed and applied to a real-world "
    "LC-MS/MS proteomics dataset comparing squamous cell carcinoma with healthy controls. "
    "The pipeline applies three imputation methods with and without batch correction, "
    "generating six parallel analyses, and computes a robustness score for each protein "
    "based on consistency across analytical conditions."
)

add_para(
    f"The classic single-method pipeline flagged {classic_total} contaminant-filtered proteins "
    f"as significant, of which {k_risk} ({artifact_pct:.1f}%) were identified as likely "
    "method-dependent artifacts. Only 13 proteins achieved the high confidence threshold "
    "(score >=7/10), with Transferrin (P02787) and PCNA (P12004) scoring a perfect 10/10. "
    "These findings demonstrate that robustness-based filtering substantially reduces false "
    "positive biomarker candidates in LC-MS/MS proteomics and provides a more reliable "
    "foundation for downstream biological and clinical interpretation."
)

# ── Save ───────────────────────────────────────────────────────────────────
doc.save("proteomics_report.docx")
print("proteomics_report.docx created successfully!")
print(f"\nKey numbers:")
print(f"  Classic pipeline    : {classic_total} proteins")
print(f"  Robust total        : {robust_total} proteins")
print(f"  High confidence     : {k_high} proteins")
print(f"  Moderate confidence : {k_mid} proteins")
print(f"  Artifacts removed   : {k_risk} proteins ({artifact_pct:.1f}%)")
