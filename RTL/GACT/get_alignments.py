import sys

args = sys.argv[1:]

ref_file = args[0]
query_file = args[1]
num_tb_files = int(args[2])

f = open(ref_file, 'r')
text = f.read()
ref_seqs  = text.splitlines()
f.close()

f = open(query_file, 'r')
text = f.read()
query_seqs  = text.splitlines()
f.close()

for i in range(1,num_tb_files+1):
    tb_file = './results/test.' + str(i) + '.out'
    f = open(tb_file, 'r')
    text = f.read()
    tb_content  = text.splitlines()
    f.close()
    
    l = len(tb_content)
    tb_pointers = tb_content[4:l-2]
    
    max_score = int(tb_content[l-2].replace(',','').split()[8])
    max_ref_pos = int(tb_content[l-2].replace(',','').split()[2])
    max_query_pos = int(tb_content[l-2].replace(',','').split()[5])

    ref_pos = max_ref_pos     
    query_pos = max_query_pos 
    
    ref = ref_seqs[i-1]
    query = query_seqs[i-1]
    
    ref_aligned = ''
    query_aligned = ''
    
    for tb in tb_pointers:
        if tb == '3':
            ref_aligned = ref[ref_pos] + ref_aligned
            query_aligned = query[query_pos] + query_aligned
            ref_pos -= 1
            query_pos -= 1
        elif tb == '1':
            ref_aligned = '-' + ref_aligned
            query_aligned = query[query_pos] + query_aligned
            query_pos -= 1
        elif tb == '2':
            ref_aligned = ref[ref_pos] + ref_aligned
            query_aligned = '-' + query_aligned
            ref_pos -= 1
    
    print(ref_aligned)
    print(query_aligned)
    print('')

