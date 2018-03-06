
seq_len = 0
seq_count = 0
with open("data/1519_alignment.fasta") as f_align:
    for line in f_align:
        if line.startswith('>'):
            seq_count += 1

        else:
            if seq_count == 2:
                break
            else:
                seq_len += len(line.strip())

print(seq_len)
# 17213
