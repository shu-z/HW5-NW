# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    NW=NeedlemanWunsch(sub_matrix_file='./substitution_matrices/BLOSUM62.mat', gap_open=-10, gap_extend=-1)
    seq1seq2_score, seq1_aligned, seq2_aligned=NW.align(seq1,seq2)

    expected_align=np.array([[0., -np.inf, -np.inf, -np.inf ],
                            [-np.inf,   5, -11, -13],
                            [-np.inf, -12,   4,  -8],
                            [-np.inf, -12,  -1,   5], 
                            [-np.inf, -14,  -6,   4]])

    
    expected_A=np.array([[-10, -np.inf, -np.inf, -np.inf],
                            [-11, -12,  -6,  -7,],
                            [-12, -13, -14,  -7,],
                            [-13, -14, -15, -12,], 
                            [-14, -15, -16, -17]])
    

    expected_B=np.array([[-10, -11, -12, -13],
                            [-np.inf,  -12, -13, -14,],
                            [-np.inf, -6, -14, -15,],
                            [-np.inf, -7,  -7, -16,], 
                            [-np.inf, -8,  -8,  -6]])

    assert np.allclose(NW._align_matrix, expected_align), 'alignment matrix not as expected'
    assert np.allclose(NW._gapA_matrix, expected_A), 'gap a matrix not as expected'
    assert np.allclose(NW._gapB_matrix, expected_B), 'gap b matrix not as expected'
   
    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    NW=NeedlemanWunsch(sub_matrix_file='./substitution_matrices/BLOSUM62.mat', gap_open=-10, gap_extend=-1)
    seq3seq4_score, seq3_aligned, seq4_aligned=NW.align(seq3,seq4)

    assert seq3seq4_score==17, "alignment score incorrect"
    assert seq3_aligned=="MAVHQLIRRP", "seq3 is aligned incorrectly"
    assert seq4_aligned=="M---QLIRHP", "seq4 is aligned incorrectly"




