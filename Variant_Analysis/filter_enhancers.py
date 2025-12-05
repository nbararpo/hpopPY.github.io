help_msg = '''
Requires python3

Usage: python filter_enhancers.py -g [gene_list] -e [enhancer_files] -o [output_file]
Required arguments:
    -g: Gene list. Can be comma separated list or a file containing gene names (one gene per row)
        Same as "--gene_list"
    -e: Folder containing enhancer annotations. This script will iterate through all files within the folder with pandas.read_table. Any file that cannot be read by the function will break the script.
        Same as "--enhancers"
    -o: output file prefix, all output files will share the same prefix
        Same as "--output"
Options:
    -h: Display this message.
        Same as "--help"
    --id_type: if gene names is official gene symbols (gene_name) or Ensembl ID (gene_id). Default is 'gene_name'.
'''

def read_gene_list(file_name):
    if not Path(file_name).exists():
        sys.exit(f'Cannot find file: {file_name}')
    out_list = []
    with open(file_name) as f:
        while True:
            l = f.readline()
            if len(l) == 0:
                break
            g = l.strip().split(',')
            out_list = out_list + g
    return(out_list)

def read_enhancers(enhancer_file,cell_line = 'auto_detect'):
    if cell_line == 'auto_detect':
        file_name = Path(enhancer_file).stem
        cell_line = file_name.split('_')[0]
    df = pd.read_table(enhancer_file,header = None)
    df = df.iloc[:,0:6]
    df.columns = ['chr','start','end','id','gene_name','strand']
    df['cell_line'] = cell_line
    return(df)

def filter_enhancers(enh, genes, id_type = 'gene_name'):
    if id_type == 'gene_name':
        col = 'gene_name'
    else:
        col = 'id'
    filt_list = enh[col].apply(lambda x: x in genes)
    out_df = enh[filt_list]
    out_df.reset_index(inplace = True, drop = True)
    return(out_df)

def merge_overlaps(df):
    df.sort_values(by = ['chr','start'], inplace = True,ignore_index = True)
    out_df = pd.DataFrame()
    for (c,i,g,s), c_df in df.groupby(['chr','id','gene_name','strand']):
        start_array = [(x[0], 1, x[1]) for x in c_df[['start','cell_line']].values]
        end_array = [(x[0], -1, x[1]) for x in c_df[['end','cell_line']].values]
        pos_array = start_array + end_array
        pos_array = sorted(pos_array, key = lambda x: x[0])
        n = 0
        cell_lines = []
        for x in pos_array:
            if n == 0: # this is a starting position
                start = x[0]
            if x[1] == 1:
                cell_lines.append(x[2])
            n += x[1]
            if n == 0: # this is an ending position
                end = x[0]
                new_line = pd.DataFrame({'chr':[c],
                                         'start':[start],
                                         'end':[end],
                                         'id':[i],
                                         'gene_name':[g],
                                         'strand':[s],
                                         'cell_line':[','.join(set(cell_lines))]})
                out_df = pd.concat([out_df,new_line],ignore_index = True)
                cell_lines = []
    return(out_df)

def main():
    arg_list = sys.argv[1:]
    short_opts = 'g:e:o:h'
    long_opts = ['gene_list=','enhancers=','output=','id_type=','help']
    try:
        opt_list = getopt.getopt(arg_list, short_opts, long_opts)[0]
    except getopt.error as error:
        sys.exit(error)
    
    if (('--help','') in opt_list) or (('-h','') in opt_list) or len(arg_list) == 0:
        print(help_msg)
        sys.exit(0)
    
    gene_list, enhancer_dir, save_name = None, None, None
    id_type = 'gene_name'

    for current_arg, current_val in opt_list:
        if current_arg in ['-g', '--gene_list']:
            if Path(current_val).exists():
                genes = read_gene_list(current_val)
            else:
                genes = current_val.strip().split(',')
                genes = [g.strip() for g in genes]
        elif current_arg in ['-o','--output']:
            save_name = current_val
        elif current_arg in ['-e','--enhancers']:
            enhancer_dir = current_val
        elif current_arg == '--id_type':
            id_type = current_val
    
    if id_type == 'gene_id':
        genes = [g.split('.')[0] for g in genes]
    
    enh_files = glob.glob(str(Path(enhancer_dir) / '*'))
    all_enhancers = pd.DataFrame()
    for e in enh_files:
        enh_pd = read_enhancers(e,cell_line = 'auto_detect')
        filtered_enh = filter_enhancers(enh_pd, genes, id_type = id_type)
        all_enhancers = pd.concat([all_enhancers,filtered_enh], ignore_index = True)
    
    all_enhancers.sort_values(by = ['chr','start'], inplace = True,ignore_index = True)
    out_df = merge_overlaps(all_enhancers)
    out_df.to_csv(save_name + '.csv',index = None)
    
    
    
if __name__ == '__main__':
    import os,sys,getopt,glob
    import pandas as pd
    from pathlib import Path
    import numpy as np
    main()