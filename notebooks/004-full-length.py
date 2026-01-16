#!/usr/bin/env python
# coding: utf-8

# In[290]:


from glob import glob
import pandas as pd


# In[291]:


def num_sort_key(s):
    s_original = s
    if "/" in s:
        s = s.split("/")[-1]
    if len(s) > 1 and s[0] == "P":
        s = s[1:]
    if s.isdigit():
        return int(s)
    else:
        return s_original
    

    
for peptide_num_path in sorted(glob("../cellranger-vdj-outputs/*"), key=num_sort_key):
    peptide_num = peptide_num_path.split("/")[-1]
    print(peptide_num, peptide_num_path)
    


# In[292]:


ls ../cellranger-vdj-outputs/P1


# In[293]:


get_ipython().system('head ../cellranger-vdj-outputs/P1/clonotypes.csv')


# In[294]:


get_ipython().system('head ../cellranger-vdj-outputs/P1/consensus.fasta')


# In[295]:


get_ipython().system('head ../cellranger-vdj-outputs/P1/clonotypes.csv')


# In[296]:


get_ipython().system('head ../cellranger-vdj-outputs/P1/filtered_contig_annotations.csv')


# In[297]:


def num_sort_key(s):
    s_original = s
    if "/" in s:
        s = s.split("/")[-1]
    if len(s) > 1 and s[0] == "P":
        s = s[1:]
    if s.isdigit():
        return int(s)
    else:
        return s_original

peptide_num_to_filtered_contig_annots = {}
for peptide_num_path in sorted(glob("../cellranger-vdj-outputs/*"), key=num_sort_key):
    peptide_num = peptide_num_path.split("/")[-1]
    df_filtered_contig_annots = pd.read_csv("../cellranger-vdj-outputs/%s/filtered_contig_annotations.csv" % peptide_num)
    df_filtered_contig_annots["sample_peptide"] = peptide_num
    df_filtered_contig_annots["sample_path"] = peptide_num_path
    peptide_num_to_filtered_contig_annots[peptide_num] = df_filtered_contig_annots
    print(peptide_num, len(df_filtered_contig_annots))

df = pd.concat(list(peptide_num_to_filtered_contig_annots.values()))
df


# In[298]:


def translate_dna(dna_seq, codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}):
    seq_len = len(dna_seq) 
    seq_len_trimmed = (seq_len // 3) * 3
    
    if seq_len != seq_len_trimmed:
        ragged_nt = dna_seq[seq_len_trimmed:]
        dna_seq = dna_seq[:seq_len_trimmed]
    else:
        ragged_nt = ""
        
    aa_seq = "".join([
            codon_table.get(dna_seq[i:i+3], 'X')
            for i in range(0, len(dna_seq), 3)
        ]
    )
    assert "X" not in aa_seq
    if "*" in aa_seq:
        ragged_nt = ""
        aa_seq = aa_seq[:aa_seq.index("*")]
        
    return aa_seq, ragged_nt


# In[300]:


# combine TCR sequence parts into full sequence
tcr_segment_aa_cols = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
tcr_segment_nt_cols = [c + "_nt" for c in tcr_segment_aa_cols]

df["combined_vdj_aa_sequence"] = df[tcr_segment_aa_cols].agg(''.join, axis=1)
df["combined_vdj_nt_sequence"] = df[tcr_segment_nt_cols].agg(''.join, axis=1)

df


# In[301]:


from tqdm import tqdm
aa_seqs = []
ragged_nts = []
for nt_seq in tqdm(df["combined_vdj_nt_sequence"]):
    aa_seq, ragged_nt = translate_dna(nt_seq)
    aa_seqs.append(aa_seq)
    ragged_nts.append(ragged_nt)
df["combined_vdj_translated_sequence"] = aa_seqs
df["combined_vdj_ragged_3p_nt"] = ragged_nt


# In[302]:


assert 1.0 == (df["combined_vdj_aa_sequence"] == df["combined_vdj_translated_sequence"]).mean()


# In[303]:


for a, b, in zip(df["combined_vdj_aa_sequence"], df['combined_vdj_translated_sequence']):
    print("---")
    print(a)
    print(b)
    break 


# In[304]:


# Are there doublets? (yes!)
index_cols = ['sample_path', 'sample_peptide', 'barcode']
index_cols_with_chain = index_cols + ['chain']
df[index_cols].value_counts()


# In[305]:


# sort the dataframe to get chains with more UMIs on top


# Assuming your original dataframe is named 'df'
# and the columns include 'barcode', 'chain_type', 'umis', 'reads', and other TCR-related columns

# Step 1: Sort entries within each barcode and chain_type by UMI count

sort_keys = index_cols_with_chain + ["umis", "reads"]
sort_ascending = [True] * len(index_cols_with_chain) + [False, False]
df = df.sort_values(sort_keys, ascending=sort_ascending)
df[sort_keys]


# In[306]:


# Create a unique identifier for each TRA/TRB entry, now sorted by UMI count
df['entry_id'] = df.groupby(index_cols_with_chain).cumcount() + 1
df


# In[307]:


df[index_cols].value_counts()


# In[308]:


index_cols_with_chain_and_entry = index_cols + ['chain', 'entry_id']
df_pivoted = df.pivot_table(
    index=index_cols, 
    columns=['chain', 'entry_id'],
    values=[col for col in df.columns if col not in index_cols_with_chain_and_entry],
    aggfunc='first'  # In case of any duplicates
)
df_pivoted


# In[309]:


# Flatten and clean up the column names
df_pivoted.columns = [f'{col[1]}_{col[2]}_{col[0]}'.strip() for col in df_pivoted.columns.values]
df_pivoted 


# In[310]:


# Add boolean flags for TRA and TRB presence

# Count the number of TRA and TRB entries
df_pivoted['TRA_count'] = df_pivoted[["TRA_1_umis", "TRA_2_umis"]].notna().sum(axis=1)
df_pivoted['TRB_count'] = df_pivoted[["TRB_1_umis", "TRB_2_umis"]].notna().sum(axis=1)

# Note the presence of any TRA and TRB chains as well as a separate boolean flag for doublets
df_pivoted['has_TRA'] = df_pivoted['TRA_count'] > 0
df_pivoted['has_TRB'] = df_pivoted['TRB_count'] > 0

df_pivoted['multi_TRA'] = df_pivoted['TRA_count'] > 1
df_pivoted['multi_TRB'] = df_pivoted['TRB_count'] > 1
df_pivoted['multi_either_TRA_or_TRB'] = df_pivoted['multi_TRA']  | df_pivoted['multi_TRB']


# Reset the index to make 'barcode' a regular column
df_final = df_pivoted.reset_index()
df_final


# In[311]:


list(df_final.columns)


# In[312]:


df_final.TRA_1_combined_vdj_nt_sequence


# In[313]:


cols_to_save = [c for c in df_final.columns if not c.endswith("_nt")]
print(cols_to_save)
df_final.to_csv("PFO004.vdj-tcr-sequences.aa-only.no-leader.no-constant.csv")

for (tra_count, trb_count), sub_df in df_final.groupby(["TRA_count", "TRB_count"]):
    sub_df.to_csv("PFO004.vdj-tcr-sequences.%d-tra.%d-trb.no-leader.no-constant.csv" % (tra_count, trb_count))


# In[314]:


ls ../../tcr-selection-cdr3-ab/v4/*.txt


# In[315]:


selected_tcr_cdr3_abs = []
cols = {"alpha": [], "beta": []}
with open("../../tcr-selection-cdr3-ab/v4/cd8-selected-tcrs-union-of-all-counting-methods.txt") as f:
    for l in f:
        l = l.strip()
        if l:
            a, b = l.split("/")
            selected_tcr_cdr3_abs.append((a,b))
            cols["alpha"].append(a)
            cols["beta"].append(b)
df_selected_cd8_cdr3_ab = pd.DataFrame(cols)


# In[316]:


"""
pd.merge(
    left: 'DataFrame | Series',
    right: 'DataFrame | Series',
    how: 'MergeHow' = 'inner',
    on: 'IndexLabel | AnyArrayLike | None' = None,
    left_on: 'IndexLabel | AnyArrayLike | None' = None,
    right_on: 'IndexLabel | AnyArrayLike | None' = None,
    left_index: 'bool' = False,
    right_index: 'bool' = False,
    sort: 'bool' = False,
    suffixes: 'Suffixes' = ('_x', '_y'),
    copy: 'bool | None' = None,
    indicator: 'str | bool' = False,
    validate: 'str | None' = None,
) -> 'DataFrame'
"""

df_join = pd.merge(
    left=df_final, 
    right=df_selected_cd8_cdr3_ab, 
    how='inner', 
    left_on=['TRA_1_cdr3', 'TRB_1_cdr3'], 
    right_on=['alpha', 'beta'])
df_join


# In[317]:


df_join[[c for c in df_join.columns if "gene" in c]]


# In[318]:


from collections import defaultdict
selected_cols = defaultdict(list)

for i, ((a, b), sub_df) in enumerate(df_join.sort_values("sample_peptide").groupby(['TRA_1_cdr3', 'TRB_1_cdr3'])):
    selected_cols["CDR3_alpha"] += [a]
    selected_cols["CDR3_beta"] += [b]
    peptide = sub_df["sample_peptide"].mode().iloc[0]
    selected_cols["Peptide"] += [peptide]
    
    sub_df_peptide = sub_df[sub_df.sample_peptide == peptide]
    
    
    full_alphas = sub_df_peptide["TRA_1_combined_vdj_aa_sequence"]
    full_betas = sub_df_peptide["TRB_1_combined_vdj_aa_sequence"]

    full_combined = full_alphas + "/" + full_betas
    most_common_combined = full_combined.mode().iloc[0]
    
    selected_cols["VDJ_alpha_beta_aa"] += [most_common_combined]
    alpha_aa = most_common_combined.split("/")[0]
    selected_cols["VDJ_alpha_aa"] += [alpha_aa]
    beta_aa = most_common_combined.split("/")[1]
    selected_cols["VDJ_beta_aa"] += [beta_aa]
    
    #### 
    # get the subset of cells that have the most common alpha/beta pairing
    ###
    mask = (full_combined == most_common_combined)
    sub_df_most_common = sub_df_peptide[mask]


    # counts for cells that support this alpha/beta CDR3
    selected_cols["Cell_Count"] += [len(sub_df_most_common)]
    selected_cols["Cell_Barcodes"] += [list(sub_df_most_common.barcode)]
    selected_cols["Alpha_Contig_IDs"] += [list(sub_df_most_common.TRA_1_contig_id)] 
    selected_cols["Beta_Contig_IDs"] += [list(sub_df_most_common.TRB_1_contig_id)] 

    selected_cols["Alpha_Reads_Mean"] += [sub_df_most_common.TRA_1_reads.mean()]
    selected_cols["Beta_Reads_Mean"] += [sub_df_most_common.TRB_1_reads.mean()]
    selected_cols["Alpha_UMIs_Mean"] += [sub_df_most_common.TRA_1_umis.mean()]
    selected_cols["Beta_UMIs_Mean"] += [sub_df_most_common.TRB_1_umis.mean()]
    
    # make sure alpha chain DNA sequence gets us back the same amino acid seq
    def modal_value(column):
        """Returns most common value for a column in the subset of cells that match
        the alpha/beta CDR3 we have grouped on"""
        col_values = sub_df_most_common[column]
        if col_values.isnull().all():
            return None
        return col_values.mode().iloc[0]
        
    alpha_dna = modal_value("TRA_1_combined_vdj_nt_sequence")
    selected_cols["VDJ_alpha_dna"] += [alpha_dna]
    alpha_dna_translated, alpha_dna_ragged_nt = translate_dna(alpha_dna)
    assert alpha_dna_translated == alpha_aa
    selected_cols["VDJ_alpha_dna_ragged_3p_nt"] += [alpha_dna_ragged_nt]

    # make sure beta chain DNA sequence gets us back the same amino acid seq
    beta_dna = modal_value("TRB_1_combined_vdj_nt_sequence")
    
    selected_cols["VDJ_beta_dna"] += [beta_dna]
    beta_dna_translated, beta_dna_ragged_nt = translate_dna(beta_dna)
    assert beta_dna_translated == beta_aa
    selected_cols["VDJ_beta_dna_ragged_3p_nt"] += [beta_dna_ragged_nt]


    
    for chain, prefix in [("alpha", "TRA_1"), ("beta", "TRB_1")]:
        # CDR1 and CDR2 for alpha/beta chains
        for other_cdr in ["cdr1", "cdr2"]:
            selected_cols["%s_%s" % (other_cdr.upper(), chain)] += [modal_value("%s_%s" % (prefix, other_cdr))]
        # include gene usage
        for gene in "vdjc":
            if gene == "d" and chain == "alpha":
                # no D gene / diversity segments in the alpha chain of the TCR
                continue 
            col_name = "%s_%s_gene" % (prefix, gene)
         
            selected_cols["%s_%s_gene" % (chain, gene)] += [modal_value(col_name)]
df_selected_full_seqs = pd.DataFrame(selected_cols);
df_selected_full_seqs


# In[319]:


df_selected_full_seqs.Cell_Barcodes


# In[320]:


sub_df_peptide


# In[321]:


df_selected_full_seqs["Cell_Count"].describe()


# In[322]:


df_selected_full_seqs.to_csv("PFO004.selected.CD8.vdj-tcr-sequences.no-leader.no-constant.csv")


# In[323]:


df_final[["TRA_count", "TRB_count"]].value_counts()


# In[324]:


from collections import Counter
counter = Counter()

for _, row in df_final[["TRA_1_combined_vdj_aa_sequence", "TRA_2_combined_vdj_aa_sequence", "TRB_1_combined_vdj_aa_sequence", "TRB_2_combined_vdj_aa_sequence"]].iterrows():
    xs = row.values
    if all([type(x) is str for x in xs]):
        key = tuple(sorted(xs))
        counter[key] += 1
        


# In[325]:


counter.most_common(10)


# In[326]:


df_final[
    df_final.TRA_1_combined_vdj_aa_sequence == "AQTVTQSQPEMSVQEAETVTLSCTYDTSESDYYLFWYKQPPSRQMILVIRQEAYKQQNATENRFSVNFQKAAKSFSLKISDSQLGDAAMYFCAYRSARGGSEKLVFGKGTKLTVNP"
]


# In[327]:


df_final[
    df_final.TRB_1_combined_vdj_aa_sequence == "SAVISQKPSRDICQRGTSLTIQCQVDSQVTMMFWYRQQPGQSLTLIATANQGSEATYESGFVIDKFPISRPNLTFSTLTVSNMSPEDSSIYLCSVDYPRQGLSPLHFGNGTRLTVT"
]


# In[328]:


df_selected_full_seqs


# In[329]:


def num_sort_key(s):
    s_original = s
    if "/" in s:
        s = s.split("/")[-1]
    if len(s) > 1 and s[0] == "P":
        s = s[1:]
    if s.isdigit():
        return int(s)
    else:
        return s_original

peptide_num_to_clonotypes = {}
for peptide_num_path in sorted(glob("../cellranger-vdj-outputs/*"), key=num_sort_key):
    peptide_num = peptide_num_path.split("/")[-1]
    df_peptide_clonotypes = pd.read_csv("../cellranger-vdj-outputs/%s/clonotypes.csv" % peptide_num)
    df_peptide_clonotypes["sample_peptide"] = peptide_num
    df_peptide_clonotypes["sample_path"] = peptide_num_path
    peptide_num_to_clonotypes[peptide_num] = df_peptide_clonotypes
    print(peptide_num, len(df_peptide_clonotypes))

df_clonotypes = pd.concat(list(peptide_num_to_clonotypes.values()))
df_clonotypes 


# In[330]:


df_clonotypes.inkt_evidence.isnull().mean()


# In[331]:


df_clonotypes.mait_evidence.isnull().mean()


# In[332]:


df_clonotypes[~df_clonotypes.mait_evidence.isnull()]


# In[333]:


get_ipython().system('ls ../cellranger-vdj-outputs/P1/')


# In[334]:


get_ipython().system('head ../cellranger-vdj-outputs/P1/filtered_contig.fasta')


# In[335]:


from tqdm import tqdm 

def num_sort_key(s):
    s_original = s
    if "/" in s:
        s = s.split("/")[-1]
    if len(s) > 1 and s[0] == "P":
        s = s[1:]
    if s.isdigit():
        return int(s)
    else:
        return s_original

def parse_fasta(path):
    results = {}
    curr_id = None
    lines = []
    def add_entry():
      
        if lines and curr_id:
            assert curr_id not in results
            results[curr_id] = "".join(lines)
            lines.clear()
        
    with open(path) as f:
        for l in tqdm(f.read().split("\n")):
            l = l.strip()
            if l.startswith(">"):
                if lines and curr_id:
                    add_entry()
                curr_id = l[1:]
            else:
                lines.append(l)
    add_entry()
    return results

peptide_num_to_contig_id_to_contig = {}
for peptide_num_path in sorted(glob("../cellranger-vdj-outputs/*"), key=num_sort_key):
    peptide_num = peptide_num_path.split("/")[-1]
    print("Reading %s" % peptide_num)
    peptide_num_to_contig_id_to_contig[peptide_num] = parse_fasta(
        "../cellranger-vdj-outputs/%s/all_contig.fasta" % peptide_num)
    


# In[336]:


df_selected_full_seqs


# In[337]:


from tqdm import tqdm
for _, row in tqdm(df_selected_full_seqs.iterrows()):
    
    for (barcode, alpha_contig_id, beta_contig_id) in zip(row.Cell_Barcodes, row.Alpha_Contig_IDs, row.Beta_Contig_IDs):
        alpha_contig = peptide_num_to_contig_id_to_contig[row.Peptide].get(alpha_contig_id)
        beta_contig = peptide_num_to_contig_id_to_contig[row.Peptide].get(beta_contig_id)
        
        
        assert barcode == alpha_contig_id.split("_")[0] == beta_contig_id.split("_")[0]
        if row.VDJ_alpha_dna not in alpha_contig:
            print("WARNING: alpha chain DNA sequence for cell %s in culture %s not contained in original assembly" % (
                barcode,
                row.Peptide
            ))
            print("Alpha AA:", row.VDJ_alpha_aa)
            print("> Alpha chain DNA:")
            print(row.VDJ_alpha_dna)
            print("> Assembly contig %s" % alpha_contig_id)
            print(alpha_contig)
        if row.VDJ_beta_dna not in beta_contig:
            print("WARNING: beta chain DNA sequence for cell %s in culture %s not contained in original assembly" % (
                barcode,
                row.Peptide
            ))
            print("Beta AA:", row.VDJ_beta_aa)
            print("> Beta chain DNA:")
            print(row.VDJ_beta_dna)
            print("> Assembly contig %s" % beta_contig_id)
            print(beta_contig)

        


# In[338]:


def find_all_start_offsets(dna_seq):
    return [i for i in range(len(dna_seq)) if dna_seq[i:i+3] == "ATG"]
    
def translate_longest_orf(dna_seq):
    longest_aa = None
    longest_offset = None
    longest_ragged_3p = None
    
    for start_offset in find_all_start_offsets(dna_seq):
        subseq = dna_seq[start_offset:]
        aa, ragged_nt = translate_dna(subseq)

        if longest_aa is None or len(aa) > len(longest_aa):
            longest_aa = aa
            longest_offset = start_offset
            longest_ragged_3p = ragged_nt
    return (longest_aa, longest_offset, longest_ragged_3p)



# In[339]:


translate_longest_orf("AGATCAGAAGAGGAGGCTTCTCACCCTGCAGCAGGGACCTGTGAGCATGGCATGCCCTGGCTTCCTGTGGGCACTTGTGATCTCCACCTGTCTTGAATTTAGCATGGCTCAGACAGTCACTCAGTCTCAACCAGAGATGTCTGTGCAGGAGGCAGAGACCGTGACCCTGAGCTGCACATATGACACCAGTGAGAGTGATTATTATTTATTCTGGTACAAGCAGCCTCCCAGCAGGCAGATGATTCTCGTTATTCGCCAAGAAGCTTATAAGCAACAGAATGCAACAGAGAATCGTTTCTCTGTGAACTTCCAGAAAGCAGCCAAATCCTTCAGTCTCAAGATCTCAGACTCACAGCTGGGGGATGCCGCGATGTATTTCTGTGCTTATCGGGATAGCAACTATCAGTTAATCTGGGGCGCTGGGACCAAGCTAATTATAAAGCCAGATATCCAGAACCCTGACCCTGCCGTGTACCAGCTGAGAGACT")


# In[340]:


df_selected_full_seqs.Cell_Barcodes


# In[341]:


from collections import Counter, defaultdict

extra_cols = defaultdict(list)
for _, row in df_selected_full_seqs.iterrows():
    alpha_leader_peptides = Counter()
    beta_leader_peptides = Counter()
    alpha_constant_peptides = Counter()
    beta_constant_peptides = Counter()
    
    ragged_alpha_3p_counts = Counter()
    ragged_beta_3p_counts = Counter()
    alpha_vdj_dna_counts = Counter()
    beta_vdj_dna_counts = Counter()
    
    aa_to_dna_counts = defaultdict(Counter)
    
    cell_barcodes = row.Cell_Barcodes
    if type(cell_barcodes) is str:
        cell_barcodes = cell_barcodes.split(";")

    alpha_contig_ids = row.Alpha_Contig_IDs
    if type(alpha_contig_ids) is str:
        alpha_contig_ids = alpha_contig_ids.split(";")

    beta_contig_ids = row.Beta_Contig_IDs
    if type(beta_contig_ids) is str:
        beta_contig_ids = beta_contig_ids.split(";")
    
    for (barcode, alpha_contig_id, beta_contig_id) in zip(cell_barcodes, alpha_contig_ids, beta_contig_ids):
        alpha_contig = peptide_num_to_contig_id_to_contig[row.Peptide].get(alpha_contig_id)
        beta_contig = peptide_num_to_contig_id_to_contig[row.Peptide].get(beta_contig_id)
        assert barcode == alpha_contig_id.split("_")[0] == beta_contig_id.split("_")[0]
        if not alpha_contig:
            print("No alpha sequence for %s!" % (barcode,))
        if not beta_contig:
            print("No beta sequence for %s!" % (barcode,))
            
        alpha_trans, alpha_trans_offset, alpha_trans_ragged_3p = translate_longest_orf(alpha_contig)
        beta_trans, beta_trans_offset, beta_trans_ragged_3p = translate_longest_orf(beta_contig)

        
        if row.VDJ_alpha_aa not in alpha_trans:
            print("WARNING: alpha chain for cell %s in culture %s not contained in translated assembly" % (
                barcode,
                row.Peptide
            ))
            print("Alpha AA:", row.VDJ_alpha_aa)
            print("Translation:", alpha_trans)
            
            print("> Alpha chain DNA:")
            print(row.VDJ_alpha_dna)
            print("> Assembly contig %s" % alpha_contig_id)
            print(alpha_contig)
        if row.VDJ_beta_aa not in beta_trans:
            print("WARNING: beta chain DNA sequence for cell %s in culture %s not contained in original assembly" % (
                barcode,
                row.Peptide
            ))
            print("Beta AA:", row.VDJ_beta_aa)
            print("Translation:", beta_trans)
            print("> Beta chain DNA:")
            print(row.VDJ_beta_dna)
            print("> Assembly contig %s" % beta_contig_id)
            print(beta_contig)
        
        alpha_contig_subseq = alpha_contig[alpha_trans_offset:]
        if alpha_dna_ragged_nt:
            alpha_contig_subseq = alpha_contig_subseq[:-len(alpha_dna_ragged_nt)]
        beta_contig_subseq = beta_contig[beta_trans_offset:]
        if beta_dna_ragged_nt:
            beta_contig_subseq = beta_contig_subseq[:-len(beta_dna_ragged_nt)]
        
        alpha_parts = alpha_trans.split(row.VDJ_alpha_aa)
        alpha_leader = alpha_parts[0]
        if len(alpha_parts) > 1:
            alpha_constant_fragment = alpha_parts[1]
        else:
            alpha_constant_fragment = ""
            
        
        beta_parts = beta_trans.split(row.VDJ_beta_aa)
        beta_leader = beta_parts[0]
        if len(beta_parts) > 1:
            beta_constant_fragment = beta_parts[1]
        else:
            beta_constant_fragment = ""
        
        
        alpha_leader_peptides[alpha_leader] += 1
        beta_leader_peptides[beta_leader] += 1
        alpha_constant_peptides[alpha_constant_fragment] += 1
        beta_constant_peptides[beta_constant_fragment] += 1
        
        alpha_leader_num_nt = 3 * len(alpha_leader)
        beta_leader_num_nt = 3 * len(beta_leader)
        alpha_vdj_num_nt = 3 * len(alpha_trans)
        assert alpha_vdj_num_nt % 3 == 0
        beta_vdj_num_nt = 3 * len(beta_trans)
        assert beta_vdj_num_nt % 3 == 0
        alpha_leader_dna = alpha_contig_subseq[:alpha_leader_num_nt]
        beta_leader_dna = beta_contig_subseq[:beta_leader_num_nt]
        
        # make sure leader sequence is in-frame and starts with start codon
        assert len(alpha_leader_dna) % 3 == 0, \
            "Leader peptide coding sequence in %s is length %d, not a multiple of 3" % (
                alpha_contig_id,
                len(alpha_leader_dna)    
            )
        assert len(beta_leader_dna) % 3 == 0, \
            "Leader peptide coding sequence in %s is length %d, not a multiple of 3" % (
                beta_contig_id,
                len(beta_leader_dna)    
            )
        
        
        assert alpha_leader_dna[:3] == "ATG", "Expected start codon"
        assert beta_leader_dna[:3] == "ATG", "Expected start codon"
        
        aa_to_dna_counts[alpha_leader][alpha_leader_dna] += 1
        aa_to_dna_counts[beta_leader][beta_leader_dna] += 1
        
        alpha_vdj_dna = alpha_contig_subseq[alpha_leader_num_nt:alpha_vdj_num_nt]
        
        beta_vdj_dna = beta_contig_subseq[beta_leader_num_nt:beta_vdj_num_nt]
        
        assert len(alpha_vdj_dna) % 3 == 0, \
            "VDJ alpha coding sequence in %s is length %d, not a multiple of 3" % (
                alpha_contig_id,
                len(alpha_vdj_dna)    
            )
        assert len(beta_vdj_dna) % 3 == 0, \
            "VDJ beta coding sequence in %s is length %d, not a multiple of 3" % (
                beta_contig_id,
                len(beta_vdj_dna)    
            )
        
        
        alpha_vdj_dna_counts[alpha_vdj_dna] += 1
        beta_vdj_dna_counts[beta_vdj_dna] += 1
        
        aa_to_dna_counts[row.VDJ_alpha_aa][alpha_vdj_dna] += 1
        aa_to_dna_counts[row.VDJ_beta_aa][beta_vdj_dna] += 1
        
        ragged_alpha_3p_counts[alpha_trans_ragged_3p] += 1
        ragged_beta_3p_counts[beta_trans_ragged_3p] += 1
        
    most_common_alpha_leader = alpha_leader_peptides.most_common(1)[0][0]
    extra_cols["assembly_alpha_leader_aa"].append(most_common_alpha_leader)
    
    most_common_beta_leader = beta_leader_peptides.most_common(1)[0][0]
    extra_cols["assembly_beta_leader_aa"].append(most_common_beta_leader)
    
    extra_cols["assembly_alpha_leader_dna"].append(
        aa_to_dna_counts[most_common_alpha_leader].most_common(1)[0][0])
    
    extra_cols["assembly_beta_leader_dna"].append(
        aa_to_dna_counts[most_common_beta_leader].most_common(1)[0][0])

    assembly_alpha_vdj_dna = alpha_vdj_dna_counts.most_common(1)[0][0]
    extra_cols["assembly_alpha_vdj_dna"].append(assembly_alpha_vdj_dna)
    assembly_beta_vdj_dna = beta_vdj_dna_counts.most_common(1)[0][0]
    extra_cols["assembly_beta_vdj_dna"].append(assembly_beta_vdj_dna)
    
    extra_cols["assembly_alpha_vdj_ragged_3p_dna"].append(ragged_alpha_3p_counts.most_common(1)[0][0])
    extra_cols["assembly_beta_vdj_ragged_3p_dna"].append(ragged_beta_3p_counts.most_common(1)[0][0])

    assembly_alpha_constant_fragment = alpha_constant_peptides.most_common(1)[0][0]
    extra_cols["assembly_alpha_constant_fragment"].append(assembly_alpha_constant_fragment)
    assembly_beta_constant_fragment = beta_constant_peptides.most_common(1)[0][0]
    extra_cols["assembly_beta_constant_fragment"].append(assembly_beta_constant_fragment)

    num_aa_assembly_alpha_constant_fragment = len(assembly_alpha_constant_fragment) 
    num_nt_assembly_alpha_constant_fragment = 3 * num_aa_assembly_alpha_constant_fragment

    num_aa_assembly_beta_constant_fragment = len(assembly_beta_constant_fragment) 
    num_nt_assembly_beta_constant_fragment = 3 * num_aa_assembly_beta_constant_fragment

    # the VDJ sequence has a few amino acids of the constant gene for both alpha and beta chains
    # so we need to split those off into their own dna sequence
    assembly_alpha_vdj_only_constant_fragment_dna = assembly_alpha_vdj_dna[-num_nt_assembly_alpha_constant_fragment:]
    assembly_alpha_vdj_without_constant_fragment_dna = assembly_alpha_vdj_dna[:-num_nt_assembly_alpha_constant_fragment]

    assembly_beta_vdj_without_constant_fragment_dna = assembly_beta_vdj_dna[:-num_nt_assembly_beta_constant_fragment]
    assembly_beta_vdj_only_constant_fragment_dna = assembly_beta_vdj_dna[-num_nt_assembly_beta_constant_fragment:]
    
    
    extra_cols["assembly_alpha_vdj_only_constant_fragment_dna"].append(assembly_alpha_vdj_only_constant_fragment_dna)
    extra_cols["assembly_alpha_vdj_without_constant_fragment_dna"].append(assembly_alpha_vdj_without_constant_fragment_dna)
    extra_cols["assembly_beta_vdj_only_constant_fragment_dna"].append(assembly_beta_vdj_only_constant_fragment_dna)
    extra_cols["assembly_beta_vdj_without_constant_fragment_dna"].append(assembly_beta_vdj_without_constant_fragment_dna)
    
    
for (col_name, col_values) in extra_cols.items():
    df_selected_full_seqs[col_name] = col_values



# In[342]:


df_selected_full_seqs.assembly_alpha_vdj_ragged_3p_dna.value_counts()


# In[343]:


from pyensembl import ensembl_grch38


# In[344]:


def find_first_stop_codon(seq, offset=0):
    seq_slice = seq[offset:] if offset else seq
    for i in range(0, len(seq_slice), 3):
        codon = seq_slice[i:i+3]
        if codon in {"TAA", "TAG", "TGA"}:
            return offset + i
    return None

print(":: alpha constant ::")
trac = ensembl_grch38.genes_by_name("TRAC")[0]
trac_transcript = trac.transcripts[0]
trac_seq = trac_transcript.sequence 
trac_stop_codon_index = find_first_stop_codon(trac_seq, offset=2)
assert trac_stop_codon_index is not None 

trac_coding_seq = trac_seq[:trac_stop_codon_index + 3]
trbc1_3p_utr = trac_seq[trac_stop_codon_index + 3:]
print("--alpha seq-->")
print(trac_seq)
print("--alpha coding seq-->")
print(trac_coding_seq)

print(":: beta constant ::")
trbc_names = ["TRBC1", "TRBC2"]

trbc_name_to_gene = {name: ensembl_grch38.genes_by_name(name)[0] for name in trbc_names}
trbc_name_to_transcript = {name: gene.transcripts[0] for (name, gene) in trbc_name_to_gene.items()}
trbc_name_to_seq = {name: t.sequence for (name, t) in trbc_name_to_transcript.items()}
trbc_name_to_stop_codon_index = {name: find_first_stop_codon(seq, offset=2) for (name, seq) in trbc_name_to_seq.items()}
for (name, stop_codon_index) in trbc_name_to_stop_codon_index.items():
    print(name, stop_codon_index)
    assert stop_codon_index is not None 
trbc_name_to_coding_seq = {name: seq[:trbc_name_to_stop_codon_index[name] + 3] for (name, seq) in trbc_name_to_seq.items()}
trbc_name_to_3p_utr = {name:  seq[trbc_name_to_stop_codon_index[name] + 3:] for (name, seq) in trbc_name_to_seq.items()}

for name in trbc_names:
    
    print("--beta seq: %s-->" % name)
    print(trbc_name_to_seq[name])
    print("--beta coding seq: %s-->" % name)
    print(trbc_name_to_coding_seq[name])


# In[345]:


translate_dna("G" + trbc_name_to_coding_seq["TRBC1"])


# In[346]:


df_selected_full_seqs["alpha_constant_dna_from_reference"] = [
    nt + trac_coding_seq
    for nt in df_selected_full_seqs.assembly_alpha_vdj_ragged_3p_dna.str.slice(0, 1)
]

df_selected_full_seqs["alpha_constant_dna"] = [
    constant_from_vdj_assembly + constant_from_reference[len(constant_from_vdj_assembly):]
    for (constant_from_vdj_assembly, constant_from_reference) in zip(
        df_selected_full_seqs.assembly_alpha_vdj_only_constant_fragment_dna,
        df_selected_full_seqs.alpha_constant_dna_from_reference)
    
]


df_selected_full_seqs["alpha_constant_translated"] = [
    translate_dna(dna)[0] for dna in df_selected_full_seqs["alpha_constant_dna"]
]
df_selected_full_seqs["beta_constant_dna_from_reference"] = [
    nt + trbc_name_to_coding_seq[trbc_gene_name] 
    for (nt, trbc_gene_name) in zip(
        df_selected_full_seqs.assembly_beta_vdj_ragged_3p_dna.str.slice(0, 1),
        df_selected_full_seqs.beta_c_gene)
        
]

df_selected_full_seqs["beta_constant_dna"] = [
    constant_from_vdj_assembly + constant_from_reference[len(constant_from_vdj_assembly):]
    for (constant_from_vdj_assembly, constant_from_reference) in zip(
        df_selected_full_seqs.assembly_beta_vdj_only_constant_fragment_dna,
        df_selected_full_seqs.beta_constant_dna_from_reference)
    
]


df_selected_full_seqs["beta_constant_translated"] = [
    translate_dna(dna)[0] for dna in df_selected_full_seqs["beta_constant_dna"]
]


df_selected_full_seqs["full_alpha_dna"] = df_selected_full_seqs[
    ["assembly_alpha_leader_dna", "assembly_alpha_vdj_without_constant_fragment_dna", "alpha_constant_dna"]].agg("".join, axis=1)

df_selected_full_seqs["full_alpha_translated"] = [
    translate_dna(dna)[0] for dna in df_selected_full_seqs["full_alpha_dna"]
]


df_selected_full_seqs["full_beta_dna"] = df_selected_full_seqs[
    ["assembly_beta_leader_dna", "assembly_beta_vdj_without_constant_fragment_dna", "beta_constant_dna"]].agg("".join, axis=1)


df_selected_full_seqs["full_beta_translated"] = [
    translate_dna(dna)[0] for dna in df_selected_full_seqs["full_beta_dna"]
]

# create versions of the coding sequence with stop codons so they can be concatenated

# make sure full sequences start with ATG, end with a stop codon, and are a multiple of 3 length

df_selected_full_seqs["full_alpha_dna_has_start_codon"] = \
    df_selected_full_seqs.full_alpha_dna.str.startswith("ATG")
df_selected_full_seqs["full_beta_dna_has_start_codon"] = \
    df_selected_full_seqs.full_beta_dna.str.startswith("ATG")

stop_codons = ("TGA", "TAG", "TAA")

df_selected_full_seqs["full_alpha_dna_ends_in_stop_codon"] = \
    df_selected_full_seqs.full_alpha_dna.str.endswith(stop_codons)
df_selected_full_seqs["full_beta_dna_ends_in_stop_codon"] = \
    df_selected_full_seqs.full_beta_dna.str.endswith(stop_codons)

df_selected_full_seqs["full_alpha_dna_frame"] = \
    df_selected_full_seqs.full_alpha_dna.str.len() % 3
df_selected_full_seqs["full_beta_dna_frame"] = \
    df_selected_full_seqs.full_beta_dna.str.len() % 3 

df_selected_full_seqs["full_alpha_dna_without_last_codon"] = \
    df_selected_full_seqs.full_alpha_dna.str.slice(0, -3)
df_selected_full_seqs["full_alpha_without_last_codon_translated"] = [translate_dna(s)[0] for s in df_selected_full_seqs.full_alpha_dna_without_last_codon]

df_selected_full_seqs["full_beta_dna_without_last_codon"] = df_selected_full_seqs.full_beta_dna.str.slice(0, -3)
df_selected_full_seqs["full_beta_without_last_codon_translated"] = [translate_dna(s)[0] for s in df_selected_full_seqs.full_beta_dna_without_last_codon]
# concatenate beta : linker : alpha chains
T2A_linker = "AGAGCCGAGGGCAGGGGAAGTCTTCTAACATGCGGGGACGTGGAGGAAAATCCCGGGCCG"
df_selected_full_seqs["T2A_linker"] = T2A_linker
df_selected_full_seqs["T2A_linker_translated"] = translate_dna(T2A_linker)[0]



df_selected_full_seqs["single_chain_full_alpha_and_beta_with_linker"] = df_selected_full_seqs[
    ["full_beta_dna_without_last_codon", "T2A_linker", "full_alpha_dna"]].agg("".join, axis=1)

df_selected_full_seqs["single_chain_full_alpha_and_beta_with_linker_translated"] = [
    translate_dna(dna)[0] for dna in df_selected_full_seqs.single_chain_full_alpha_and_beta_with_linker
]
df_selected_full_seqs


# In[347]:


df_selected_full_seqs.beta_constant_translated


# In[348]:


df_selected_full_seqs.full_alpha_translated.str.len()


# In[349]:


df_selected_full_seqs.full_alpha_translated


# In[350]:


df_selected_full_seqs.assembly_beta_leader_aa


# In[351]:


df_selected_full_seqs.full_beta_dna


# In[352]:


df_selected_full_seqs.full_beta_translated


# In[353]:


df_selected_full_seqs.full_beta_translated.str.len()


# In[354]:


df_selected_full_seqs.full_beta_dna.str.slice(-3, None)


# In[355]:


df_selected_full_seqs[["full_alpha_translated", "full_beta_translated", "single_chain_full_alpha_and_beta_with_linker_translated"]].map(len)


# In[356]:


df_selected_full_seqs.full_beta_dna.str.len().map(lambda x: x%3)


# In[357]:


for c in df_selected_full_seqs.columns:
    print(c)


# In[358]:


for _, row in df_selected_full_seqs.iterrows():
    assert row.VDJ_alpha_aa in row.single_chain_full_alpha_and_beta_with_linker_translated
    assert row.VDJ_beta_aa in row.single_chain_full_alpha_and_beta_with_linker_translated


# In[359]:


import pandas as pd
import numpy as np

df_selected_full_seqs = df_selected_full_seqs.map(
    lambda x: ";".join(map(str, x)) if isinstance(x, list) else x)

# sort by peptide number and cell count
df_selected_full_seqs["Peptide_Num"] = [
    num_sort_key(p) for p in list(df_selected_full_seqs.Peptide)
]

df_selected_full_seqs = df_selected_full_seqs.sort_values(by=["Peptide_Num", "Cell_Count"])
df_selected_full_seqs


# In[360]:


df_selected_full_seqs.columns


# In[361]:


SEQUENCE_COLUMNS = ["full_beta_without_last_codon_translated", "T2A_linker_translated", "full_alpha_translated"]
KEY_COLUMNS = ["Peptide_Num", "CDR3_alpha", "CDR3_beta"]
QC_KMER_SIZE = 10
MIN_REPEAT = 2
MIN_CDR3_LENGTH = 5
MAX_CDR3_LENGTH = 40
MIN_CHAIN_LENGTH = 200
MAX_CHAIN_LENGTH = 450
MIN_READS = 20

def extract_full_seq(row, sequence_columns=SEQUENCE_COLUMNS):
    return "".join([row[col] for col in sequence_columns])


def repeated_kmers_in_seq(seq, k=QC_KMER_SIZE, min_repeat=MIN_REPEAT):
    kmer_counts = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        kmer_counts[kmer] += 1
    return [(kmer, c) for (kmer, c) in kmer_counts.most_common() if c >= min_repeat]

def repeated_kmers_in_row(row, sequence_columns=SEQUENCE_COLUMNS, k=QC_KMER_SIZE, min_repeat=MIN_REPEAT):
    concat_seq = extract_full_seq(row, sequence_columns=sequence_columns)
    return repeated_kmers_in_seq(concat_seq, k=k, min_repeat=min_repeat)

def repeated_kmers_in_df(df, key_columns=KEY_COLUMNS, sequence_columns=SEQUENCE_COLUMNS, k=QC_KMER_SIZE, min_repeat=MIN_REPEAT):
    result = {}
    for _, row in df.iterrows():
        key = tuple(row[key_columns].values)
        kmers = repeated_kmers_in_row(row, sequence_columns=sequence_columns, k=k, min_repeat=min_repeat)
        if kmers:
            result[key] = kmers
    return result



def qc_check(df, key_columns=KEY_COLUMNS, sequence_columns=SEQUENCE_COLUMNS, k=QC_KMER_SIZE, min_repeat=MIN_REPEAT):

   
    qc_cols = ["full_alpha_dna_has_start_codon", "full_beta_dna_has_start_codon", "full_alpha_dna_ends_in_stop_codon", "full_beta_dna_ends_in_stop_codon"]
    for col in qc_cols:
        print("Checking %s..." % col)
        values = df[col]
        if not values.all():
            raise ValueError("Found False value in %s" % (col,))

    nonzero_cols = ["Cell_Count", "Alpha_UMIs_Mean", "Beta_UMIs_Mean"]
    for col in nonzero_cols:
        
        print("Checking % s > 0..." % col)
        values = df[col]
        if not (values > 0).all():
            raise ValueError("Found 0 value in %s" % (col,))

    read_count_cols = ["Alpha_Reads_Mean", "Beta_Reads_Mean"]
    for col in read_count_cols:
        
        print("Checking % s > %d..." % (col, MIN_READS))
        values = df[col]
        if not (values > MIN_READS).all():
            raise ValueError("Found value < %d in %s" % (MIN_READS, col,))

    for chain in ["alpha", "beta"]:
        for cdr in ["CDR1", "CDR2", "CDR3"]:
            col_name = "%s_%s" % (cdr, chain)
            print("Checking %s in full construct..." % (col_name))
            for _, row in df.iterrows():
                seq = extract_full_seq(row, sequence_columns=sequence_columns)
                if row[col_name] not in seq:
                    raise ValueError("%s not in full sequence!" % col_name)
                if seq.count(row[col_name]) > 1:
                    raise ValueError("%s occurs more than once!" % col_name)
            
    print("Checking CDR3_alpha length...")
    for _, row in df.iterrows():
        cdr3a = row["CDR3_alpha"]
        if len(cdr3a) < MIN_CDR3_LENGTH:
            raise ValueError("CDR3_alpha too short!")
        if len(cdr3a) > MAX_CDR3_LENGTH:
            raise ValueError("CDR3_alpha too long!")
    
    print("Checking CDR3_beta length...")
    for _, row in df.iterrows():
        cdr3b = row["CDR3_beta"]
        if len(cdr3b) < MIN_CDR3_LENGTH:
            raise ValueError("CDR3_beta too short!")
        if len(cdr3b) > MAX_CDR3_LENGTH:
            raise ValueError("CDR3_beta too long!")
    
    print("Checking for repeated kmers (k=%d)..." % (k,))
    if repeated_kmers_in_df(df, key_columns=key_columns, sequence_columns=sequence_columns, k=k, min_repeat=min_repeat):
        raise ValueError("Repeated kmers in amino acid sequences!")
    
    print("Checking full chain length...")
    for chain in ["alpha", "beta"]:
        col_name = "full_%s_translated" % chain
        for seq in df[col_name]:
            if len(seq) < MIN_CHAIN_LENGTH:
                raise ValueError("%s chain too short!" % chain)
            if len(seq) > MAX_CHAIN_LENGTH:
                raise ValueError("%s chain too long!" % chain)

    for chain in ["alpha", "beta"]:
        print("Checking %s constant chain sequence..." % chain)
        gene_col_name = "%s_c_gene" % chain
        constant_seq_col_name = "%s_constant_translated" % chain
        full_seq_col_name = "full_%s_translated" % chain
        for gene_name, constant_seq, chain_seq, full_seq  in zip(df[gene_col_name], df[constant_seq_col_name], df[full_seq_col_name], df['single_chain_full_alpha_and_beta_with_linker_translated']):
            expected = None

            expected = {"TRAC": "LLMTLRLWSS", "TRBC1": "VKRKDF", "TRBC2": "VKRKDSRG"}.get(gene_name)
            if not constant_seq.endswith(expected):
                raise ValueError("Wrong constant sequence for %s (%s) in translated constant sequence" % (chain, gene_name))
            if not chain_seq.endswith(expected):
                raise ValueError("Wrong constant sequence for %s (%s) in full chain sequence" % (chain, gene_name))
            if not expected:
                raise ValueError("Unrecognized constant gene: %s" % gene_name)
            if expected not in full_seq:
                raise ValueError("Constant sequence for %s not in single_chain_full_alpha_and_beta_with_linker_translated" % gene_name) 
                
            
    print("---")
    print("QC check: OK")

    


# In[362]:


qc_check(df_selected_full_seqs)


# In[366]:


df_selected_full_seqs.Beta_UMIs_Mean


# In[372]:


# order the 
SORT_COLUMNS = ["Peptide_Num", "Cell_Count", "CDR3_alpha", "CDR3_beta", "Beta_UMIs_Mean", "Alpha_UMIs_Mean"]
df_selected_full_seqs = df_selected_full_seqs.sort_values(by=SORT_COLUMNS).reset_index().drop(columns=["index"])
df_selected_full_seqs["TCR_Num"] = np.arange(len(df_selected_full_seqs)) + 1
print(df_selected_full_seqs[["Peptide_Num", "CDR3_alpha", "CDR3_beta", "Cell_Count"]])
df_selected_full_seqs.to_csv("PFO004.selected.full-tcr-sequences.v5.csv", index=False)
df_selected_full_seqs


# In[371]:





# In[373]:


from collections import defaultdict
from reportlab.pdfgen import canvas
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from itertools import cycle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from collections import Counter


def create_tcr_pdf(df, output_path, sequence_font_size = 18, label_font_size = 12, title_font_size = 14, chars_per_line = 45,
                   color_list = [colors.darkblue, colors.blue, colors.lightblue, 
                                 colors.orange, 
                                 colors.darkgreen, colors.green, colors.lightgreen, 
                                 colors.purple, colors.brown, colors.pink, colors.darkgray, colors.red],
                   sequence_columns={
                        "assembly_beta_leader_aa": "beta_leader",
                        "VDJ_beta_aa": "beta_VDJ",
                        "beta_constant_translated": "beta_constant",
                        "T2A_linker_translated": "T2A_linker",
                        "assembly_alpha_leader_aa": "alpha_leader",
                        "VDJ_alpha_aa": "alpha_VDJ",
                        "alpha_constant_translated": "alpha_constant"}):
    """
    Creates a PDF with color-coded TCR sequences using monospaced font.
    """
    # Register Courier as our monospaced font
    pdfmetrics.getFont('Courier')
    
   
    color_cycle = cycle(color_list)
    color_map = defaultdict(lambda: next(color_cycle))
    
    c = canvas.Canvas(output_path, pagesize=letter)
    width, height = letter
    

    char_width = sequence_font_size * 0.6  # Fixed width for monospaced font

        
    for _, row in df.iterrows():
        y_position = height - 50
        def write(s, x=20, y_offset=20):
            nonlocal y_position
            c.drawString(x, y_position, s)
            y_position -= y_offset
        def blank():
            nonlocal y_position
            y_position -= 30
            
        # TCR identifier
        c.setFont("Helvetica-Bold", title_font_size)
        c.setFillColor("red")
        write(f"TCR #{row.TCR_Num}")
        c.setFillColor("black")
        write(f"Peptide {row.Peptide}")
        write(f"CDR3 alpha = {row.CDR3_alpha}")
        write(f"CDR3 beta = {row.CDR3_beta}")

        # normal font weight for smaller bits of information
        c.setFont("Helvetica", title_font_size)
        write(f"Length = {sum([len(row[col]) for col in sequence_columns])}", x=40)
    
        for col in df.columns:
            if col.endswith("_gene"):
                write(f"{col} = {row[col]}", x=40)
        
        blank()
        
        # Color legend
        legend_x = width - 200
        legend_y = height - 50
        c.setFont("Helvetica", label_font_size)
        for col, alias in sequence_columns.items():
            color = color_map[col]
            c.setFillColor(color)
            c.drawString(legend_x, legend_y, alias)
            legend_y -= label_font_size * 1.5
        
        # Write wrapped sequence with color coding
        x_position = 50
        current_line_width = 0

        combined = ""
    
        for col in sequence_columns.keys():
            
            sequence = row[col]
            color = color_map[col]
            
            for char in sequence:
                combined += char
                if current_line_width >= chars_per_line:
                    y_position -= sequence_font_size * 1.2
                    current_line_width = 0
                    x_position = 50
                
                c.setFont("Courier", sequence_font_size)
                c.setFillColor(color)
                c.drawString(x_position, y_position, char)
                x_position += char_width
                current_line_width += 1
        if combined != row['single_chain_full_alpha_and_beta_with_linker_translated']:
            raise ValueError("Expected '%s' but got '%s'" % (row['single_chain_full_alpha_and_beta_with_linker'], combined))
        
        repeated_kmers_with_counts = repeated_kmers_in_row(row, sequence_columns)
        if repeated_kmers_with_counts and repeated_kmers_with_counts[0][1] > 1:
            y_position -= 50
            for kmer, count in repeated_kmers_with_counts[:10]:
                c.setFont("Helvetica", title_font_size)
                c.setFillColor("red")
                message = "Warning! Repeated k-mer (n=%d): %s" % (count, kmer)
                c.drawString(20, y_position, message)
                print(message)
                y_position -= 20
        c.showPage()
    
    c.save()


# In[374]:


create_tcr_pdf(
    df_selected_full_seqs, 
    output_path="PFO004.selected.report.v5.pdf", 
   )


# In[375]:


df_selected_full_seqs.columns


# In[376]:


[c for c in df_selected_full_seqs.columns if "constant" in c]


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




