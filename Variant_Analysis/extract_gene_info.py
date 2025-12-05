#!/usr/bin/python3

help_msg = '''
Requires python3

Usage: python extract_gene_info.py -f [gtf_file] -g [gene_list] -o [output_file]
Required arguments:
    -f: gtf file to search in. 
        Same as "--gtf-file"
    -g: Gene list. Can be comma separated list or a file containing gene names (one gene per row)
        Same as "--genes"
    -o: output file prefix, all output files will share the same prefix
        Same as "--output"
Options:
    
    --flanking: Define the direction of flanking regions to include. Valid options are 'none', 'up' (default), 'down', or 'both'.
    --flanking-size: Define the size of flanking regions in base pairs. Default 5000
    --id-type: gene_id (ENSG ids) or gene_name (official gene symbols), default gene_name
    -h: Display this message.
        Same as "--help"
'''



def make_info_dict(l):
    dict_items = [i.split(' ') for i in l.split('; ')]
    dict_keys = [i[0] for i in dict_items]
    dict_values = [i[1].strip('"') for i in dict_items]
    out_dict = dict(zip(dict_keys,dict_values))
    return(out_dict)

def is_gene(l,gene_list,id_type):
    l = l.strip().split('\t')
    l[8] = make_info_dict(l[8])
    if (l[8][id_type] in gene_list):
        return([l[0],l[2],l[3],l[4],l[6],l[8]])
    else:
        return(None)

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

# directions is one of 'up', 'down', or 'both'
def flank(lst,direction = 'up', size = 5000):
    c, in_feat, start, end, strand, gene_name, gene_id, transcript_name, exon_num = lst
    if in_feat != 'gene':
        return()
    left_flank = [str(int(start) - size),start]
    right_flank = [end,str(int(end) + size)]
    if strand == '+':
        upstream = (c, 'upstream', left_flank[0], left_flank[1], strand, gene_name, gene_id,transcript_name, exon_num)
        downstream = (c, 'downstream', right_flank[0], right_flank[1], strand, gene_name, gene_id, transcript_name, exon_num)
    if strand == '-':
        upstream = (c, 'upstream', right_flank[0], right_flank[1], strand, gene_name, gene_id, transcript_name, exon_num)
        downstream = (c, 'downstream', left_flank[0], left_flank[1], strand, gene_name, gene_id, transcript_name, exon_num)
    if strand == '.':
        upstream = (c, 'left_flanking', left_flank[0], left_flank[1], strand, gene_name, gene_id, transcript_name, exon_num)
        downstream = (c, 'right_flanking', right_flank[0], right_flank[1], strand, gene_name, gene_id, transcript_name, exon_num)
    
    if direction == 'both':
        return([upstream,downstream])
    if direction == 'up':
        return([upstream])
    if direction == 'down':
        return([downstream])
    

def main():
    arg_list = sys.argv[1:]
    short_opts = 'g:f:o:h'
    long_opts = ['help','genes=','gtf-file=','output=','flanking=','flanking-size=','id-type=']

    try:
        opt_list = getopt.getopt(arg_list, short_opts, long_opts)[0]
    except getopt.error as error:
        sys.exit(error)
    
    if (('--help','') in opt_list) or (('-h','') in opt_list) or len(arg_list) == 0:
        print(help_msg)
        sys.exit(0)

    file,genes,outname = None, None, None
    flanking, flanking_size = 'up', 5000
    id_type = 'gene_name'
    for current_arg, current_val in opt_list:
        if current_arg in ['-f','--gtf-file']:
            file = current_val
            if not Path(file).exists():
                sys.exit(f'Cannot find file: {file}')
            print(f'gtf file: {os.path.abspath(file)}')
        elif current_arg in ['-o','--output']:
            outname = current_val
        elif current_arg in ['-g','--genes']:
            if Path(current_val).exists():
                genes = read_gene_list(current_val)
            else:
                genes = current_val.strip().split(',')
                genes = [g.strip() for g in genes]
        elif current_arg == '--flanking':
            flanking = current_val
        elif current_arg == '--flanking-size':
            flanking_size = int(current_val)
        elif current_arg == '--id-type':
            id_type = current_val

    if any([file is None, genes is None, outname is None]):
        print('Missing at least one required argument.')
    
    annotation = outname + '_annotations.txt'
    bed = outname + '.bed'
    
    anno = open(annotation, 'w')
    b = open(bed, 'w')
    header = ['chr','feature','start','end','strand','gene_name','gene_id','transcript_name','exon_number']
    anno.write('\t'.join(header) + '\n')
    b.write('#chr\tstart\tend\tgene\n')
    current_chr = ''
    i = 0
    gene_name_list = [id_type + ' "' + g + '"' for g in genes]
    with open(file) as f:
        l = f.readline()
        while l.startswith('#'):
            l = f.readline()
            i += 1
        while True:
            i += 1
            if i % 100000 == 0:
                print(f'{i} rows traversed')
            if len(l) == 0:
                print('All rows traversed.')
                break
            if current_chr != l.split('\t')[0]:
                current_chr = l.split('\t')[0]
                print(f'Traversing {current_chr}')
            if any([(g in l) for g in gene_name_list]):
                parsed_line = is_gene(l,genes,id_type)
            else:
                l = f.readline()
                continue
            if parsed_line is None:
                l = f.readline()
                continue
            out_list = parsed_line[:5]
            info_dict = parsed_line[5]
            gene_name = info_dict['gene_name']
            gene_id = info_dict['gene_id']
            try:
                transcript_name = info_dict['transcript_name']
            except:
                transcript_name = 'NA'
            try:
                exon_number = info_dict['exon_number']
            except:
                exon_number = 'NA'
            out_list = out_list + [gene_name,gene_id,transcript_name,exon_number]
            anno.write('\t'.join(out_list) + '\n')
            if out_list[1] == 'gene':
                flanks = flank(out_list,direction = flanking, size = flanking_size)
                b.write('\t'.join([out_list[0],out_list[2],out_list[3],out_list[5]]) + '\n')
                for r in flanks:
                    anno.write('\t'.join(r) + '\n')
                    b.write('\t'.join([r[0],r[2],r[3],r[5] + '_' + r[1]]) + '\n')
            l = f.readline()
    anno.close()
    b.close()

if __name__ == '__main__':
    import os,sys,getopt
    from pathlib import Path
    main()
