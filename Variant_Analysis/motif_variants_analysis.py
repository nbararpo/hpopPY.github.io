class motif_matrix:
    # class variables
    base_position = {'A':0,'C':1,'G':2,'T':3}
    
    def __init__(self,motif_file = None):
        self.name = ''
        self.pfm = np.array([])
        self.str_rep = ''
        self.length = 0
        if motif_file is not None:
            self.read_motif(motif_file)
    
    def display(self):
        print(f'Motif name: {self.name}')
        print(f'Motif string: {self.str_rep}')
        print(f'Motif length: {self.length}')
        print(f'Motif Position Frequency Matrix (PFM): \n{self.pfm}')
        print(f'Maximum score possible with the PFM: {self.max_possible_score}')
    
    def read_motif(self,file):
        with open(file,'r') as f:
            header = f.readline()
            matrix = f.readlines()
        header = header.strip().split('\t')
        self.str_rep = header[0].strip('>')
        self.name = header[1]
        matrix = np.array([[float(p) if float(p) > 0 else float(p) + 0.001 for p in row.strip().split('\t')] for row in matrix])
        matrix = np.array([r / sum(r) for r in matrix])
        self.pfm = matrix
        self.length = len(self.pfm)
        self.best_match_scores = [np.log2(max(x)/0.25) for x in self.pfm]
        self.max_possible_score = sum(self.best_match_scores)
    
    @staticmethod
    def score_base(base,p_list,base_position = base_position):
        base_p = p_list[base_position[base]]
        return(np.log2(base_p/0.25))
    
    def score_sequence(self, seq, base_position = base_position):
        if len(seq) != self.length:
            raise BaseException(f'Cannot calculate sequence score: sequence need to be {self.length} base long.')
        seq = seq.upper()
        base_prob_pairs = zip(seq,self.pfm)
        score = sum([self.score_base(x[0],x[1]) for x in base_prob_pairs])
        return(score)

def seq_query(chrom,start,end, genome = 'hg38'):
    url = f'https://api.genome.ucsc.edu/getData/sequence?genome={genome};chrom={chrom};start={start};end={end}'
    result = urllib.request.urlopen(url).read().decode().strip()
    result = json.loads(result)
    seq = result['dna']
    return(seq)

def rev_comp(seq):
    comp = {'A':'T','C':'G','G':'C','T':'A'}
    rc_seq = [comp[n] for n in reversed(seq.upper())]
    return(''.join(rc_seq))

def construct_alt(df_row, annotation_base = '1-based'):
    chrom, motif_start, motif_end, motif_strand = df_row['chr'], df_row['motif_start'], df_row['motif_end'], df_row['motif_strand']
    if annotation_base == '0-based': # convert bed annotation convention for data query
        motif_start += 1
    expected_len = motif_end - motif_start
    variant_start, variant_end = df_row['variant_start'], df_row['variant_end']
    variants = df_row['ref'].split(',') + df_row['alt'].split(',')
    variants = ['' if v == '*' else v for v in variants]
    # prepare motif and flanking region sequences
    max_indel_len = max([len(v) for v in variants]) - min([len(v) for v in variants])
    query_start, query_end = (motif_start - max_indel_len), (motif_end + max_indel_len)
    motif_seq = seq_query(chrom,query_start,query_end)
    motif_seq_parts = [motif_seq[:min((variant_start - query_start),(motif_start - query_start))], #0 left flanking
                       motif_seq[min((variant_start - query_start),(motif_start - query_start)):(variant_start - query_start)], #1 left of variant
                       motif_seq[(variant_start - query_start):(variant_end - query_start)], #2 variant
                       motif_seq[(variant_end - query_start):max((variant_end - query_start),(motif_end - query_start))], #3 right of variant
                       motif_seq[max((variant_end - query_start),(motif_end - query_start)):]] #4 right_flanking
    seqs = {'0':[motif_seq[(motif_start - query_start):(motif_end - query_start)]]}
    for i in range(1,len(variants)):
        alt_seq_parts = motif_seq_parts
        alt_seq_parts[2] = variants[i]
        if len(''.join(alt_seq_parts[1:4])) == expected_len: # snps
            seqs.update({str(i):[''.join(alt_seq_parts[1:4])]})
        else: # indel
            left_preserve = ''.join(alt_seq_parts[1:5])[:expected_len]
            right_preserve = ''.join(alt_seq_parts[0:4])[-expected_len:]
            seqs.update({str(i) : [left_preserve,right_preserve]})
    if motif_strand == '-':
        for i in seqs:
            seqs.update({str(i):[rev_comp(s) for s in seqs[i]]})
    return(seqs)
    
def calculate_scores(df_row,annotation_base = '0-based',motif_file_dir = None):
    if motif_file_dir is None:
        return(None)
    motif = motif_matrix(Path(motif_file_dir) / df_row['motif_file'])
    seqs = construct_alt(df_row, annotation_base = annotation_base)
    scores = dict()
    for i in seqs:
        seqs_scores = {s:motif.score_sequence(s) for s in seqs[i]}
        max_score = max(seqs_scores.values())
        scores.update({i:[(seq,score) for seq,score in seqs_scores.items() if score == max_score][0]})
    return(scores)

def read_data(data_tables,n_info_cols = 21):
    out_df = pd.DataFrame()
    for f in data_tables:
        df = pd.read_table(f)
        out_df = pd.concat([out_df,df], ignore_index = True)
    genotype_df = out_df.iloc[:,n_info_cols:].copy()
    info_df = out_df[['chr','motif_start', 'motif_end', 'motif_name', 'motif_strand', 'motif_file',
                      'enh_start', 'enh_end','gene_id', 'gene_name', 'gene_strand', 'cell_lines',
                      'start','end', 'ref', 'alt','features']].copy()
    info_df.columns = ['chr','motif_start', 'motif_end', 'motif_name', 'motif_strand', 'motif_file',
                       'enhancer_start', 'enhancer_end','gene_id', 'gene_name', 'gene_strand', 'cell_lines',
                       'variant_start','variant_end', 'ref', 'alt','features']
    out_df = pd.concat([info_df,genotype_df],axis = 1)
    out_df.sort_values(by = ['chr','variant_start'], inplace = True)
    out_df.reset_index(inplace = True, drop = True)
    return(out_df)

def read_ethnic_dict(file,header = True):
    eth_dict = dict()
    with open(file,'r') as f:
        if header == True:
            l = f.readline() # Skipping header
        while True:
            l = f.readline()
            if len(l) == 0:
                break
            l = l.strip().split(',')
            eth_dict.update({l[0]:l[1]})
    return(eth_dict)

def split_ethnic(df,eth_dict,n_info_cols = 17):
    info_cols = df.iloc[:,:n_info_cols].copy()
    sample_cols = df.iloc[:,n_info_cols:].copy()
    all_samples = df.columns[n_info_cols:]
    eth_list = np.array([eth_dict[x] for x in all_samples])
    eth_set = set(eth_list)
    output_dict = dict()
    for e in eth_set:
        include_columns = all_samples[eth_list == e]
        e_sample_df = sample_cols.loc[:,include_columns].copy()
        e_df = pd.concat([info_cols,e_sample_df],axis = 1)
        output_dict.update({e:e_df})
    return(output_dict)

def count_genotypes(df_row,n_info_cols = 6):
    n_alleles = len(df_row['alt'].split(',')) + 1
    allele_counts = dict(zip([str(x) for x in range(n_alleles)],[0] * n_alleles))
    sample_genotypes = list(df_row)[n_info_cols:]
    genotypes, counts = np.unique(sample_genotypes,return_counts=True)
    genotype_counts = dict(zip(genotypes,counts))
    if './.' in genotype_counts:
        del genotype_counts['./.']
    for g,n in genotype_counts.items():
        alleles = g.split('/')
        for a in alleles:
            allele_counts[a] += n
    return([genotype_counts,allele_counts])

def ethnic_genotype_count(df,eth_dict,n_info_cols = 17):
    info_cols = df.iloc[:,:n_info_cols].copy()
    sample_cols = df.iloc[:,n_info_cols:].copy()
    all_samples = df.columns[n_info_cols:]
    eth_list = np.array([eth_dict[x] for x in all_samples])
    eth_set = set(eth_list)
    out_df = info_cols.copy()
    for e in eth_set:
        include_columns = all_samples[eth_list == e]
        e_sample_df = sample_cols.loc[:,include_columns].copy()
        e_sample_df = pd.concat([info_cols,e_sample_df],axis = 1)
        counts = e_sample_df.apply(lambda x: count_genotypes(x,n_info_cols = n_info_cols), axis = 1)
        out_df[f'{e}_genotype_counts'] = [c[0] for c in counts]
        out_df[f'{e}_allele_counts'] = [c[1] for c in counts]
    return(out_df)

def compare_ethnic(df,eth1,eth2):
    eth_1_allele_counts = df[f'{eth1}_allele_counts'].values
    eth_2_allele_counts = df[f'{eth2}_allele_counts'].values
    allele_counts = list(zip(eth_1_allele_counts, eth_2_allele_counts))
    cont_tables = [[[x[1]+0.01 for x in sorted(c1.items())],[x[1]+0.01 for x in sorted(c2.items())]] for c1,c2 in allele_counts]
    p_list = [scipy.stats.chi2_contingency(x).pvalue for x in cont_tables]
    return(p_list)

def fdr(ps):
    ranked_ps = scipy.stats.rankdata(ps)
    fdr = [p * len(ps) / r for p,r in zip(ps, ranked_ps)]
    fdr = [1 if x > 1 else x for x in fdr]
    return(fdr)

help_msg = '''
Requires python3

Usage: python motif_variants_analysis.py -f [data_file_list] -o [output_file] --ethnic_table [ethnic_table]
Required arguments:
    -f: data files, join multiple files by comma-separated list (no spaces) 
        Same as "--files"
        
    -o: output file prefix, all output files will share the same prefix
        Same as "--output"
        
    --ethnic_table: two-column table defining sample ethnicity. 
    
    --contrast: specify ethnicities to compare. e.g. "Asian,Caucasian". Default is "all" which compares all ethnicity pair-wise
    
    --config: Use external file to set up all inputs, one line per argument
        e.g. motif_bed_dir = /path/to/directory
    
Options:
    --motif_matrix_dir: Path to directory with motif matrix files. Required to calculate binding scores.
    
    -h: Display this message.
        Same as "--help"
'''

def main():
    arg_list = sys.argv[1:]
    short_opts = 'f:o:h'
    long_opts = ['help','files=','output=','motif_matrix_dir=','contrast=','ethnic_table=','config=']
    try:
        opt_list = getopt.getopt(arg_list, short_opts, long_opts)[0]
    except getopt.error as error:
        sys.exit(error)
    
    if (('--help','') in opt_list) or (('-h','') in opt_list) or len(arg_list) == 0:
        print(help_msg)
        sys.exit(0)

    data_files = None
    ethnic_table, motif_matrix_dir = None, None # These can have default values
    contrast = 'all'
    save_name = None
    config_file = None
    
    for current_arg, current_val in opt_list:
        if current_arg in ['-f','--files']:
            data_files = current_val.split(',')
        elif current_arg in ['-o','--output']:
            save_name = current_val
        elif current_arg == '--motif_matrix_dir':
            motif_matrix_dir = current_val
        elif current_arg == '--ethnic_table':
            ethnic_table = current_val
        elif current_arg == '--contrast':
            contrast = current_val
        elif current_arg == '--config':
            config_file = current_val
            
    # Use config file to set up inputs without having to define them in the command line
    if config_file is not None: 
        with open(config_file) as f:
            configs = f.readlines()
        configs = [l.split('=') for l in configs]
        configs = [[k.strip(),v.strip()] for k,v in configs]
        for k,v in configs:
            if k in ['f','files']:
                data_files = v.split(',')
            elif k in ['o','output']:
                save_name = v
            elif k == 'motif_matrix_dir':
                motif_matrix_dir = v
            elif k == 'ethnic_table':
                ethnic_table = v
            elif k == 'contrast':
                contrast = v
    
    data_df = read_data(data_files)
    eth = read_ethnic_dict(ethnic_table)
    e_df = ethnic_genotype_count(data_df,eth,n_info_cols = 17)
    
    if motif_matrix_dir is not None:
        print('   Calculating scores for ref and alt alleles')
        scores_dict = dict(e_df.apply(lambda x: calculate_scores(x,annotation_base = '0-based',motif_file_dir = motif_matrix_dir), axis = 1))
        e_df['ref_seq'] = [v['0'][0] for v in scores_dict.values()]
        e_df['ref_score'] = [v['0'][1] for v in scores_dict.values()]
        e_df['alt_seqs'] = [','.join([seq[0] for a,seq in v.items() if a != '0']) for v in scores_dict.values()]
        e_df['alt_scores'] = [','.join([str(seq[1]) for a,seq in v.items() if a != '0']) for v in scores_dict.values()]
    
    
    if contrast == 'all':
        contrast_list = [['Asian','Caucasian'],
                         ['Asian','Indian'],
                         ['Asian','Other'],
                         ['Caucasian','Indian'],
                         ['Caucasian','Other'],
                         ['Indian','Other']]
    else:
        contrast_list = [contrast.split(',')]
    for a,b in contrast_list:
        out_df = e_df.copy()
        print(f'Processing contrast {a}-{b}')
        output_file_name = save_name + '_' + a + '-' + b + '.txt'
        out_df['contrast'] = f'{a}-{b}'
        out_df['p_value'] = compare_ethnic(e_df,a,b)
        out_df['fdr'] = fdr(out_df['p_value'].values)
        out_df.to_csv(output_file_name,sep = '\t',index = False)
    return()



if __name__ == '__main__':
    import os,sys,getopt,glob, json
    import pandas as pd
    from pathlib import Path
    import numpy as np
    import scipy.stats
    import urllib.request
    main()