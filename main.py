# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
 

    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    fa_list=[gg_seq, mm_seq, br_seq, tt_seq]
    fa_name=['Gallus_gallus', 'Mus_musculus', 'Balaeniceps_rex', 'tursiops_truncatus']
    fa_score={}

    for i in range(0,4):
        NW=NeedlemanWunsch(sub_matrix_file='./substitution_matrices/BLOSUM62.mat', gap_open=-10, gap_extend=-1)
        score, _, _=NW.align(hs_seq, fa_list[i])
        fa_score[fa_name[i]]=score
        #print all alignment scores
        print("{} Alignment Score: {}".format(fa_name[i], score))
    
    

    print("Species similarity order to humans (most to least):")
    for key, value in sorted(fa_score.items(), reverse=True):
        print(key)

    

if __name__ == "__main__":
    main()
