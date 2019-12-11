import pytest
import os
from netZooPy.panda.panda import Panda
import pandas as pd
import numpy as np

def test_panda():
    #print(os.getcwd())
    print('Start Panda run ...')
    ppi            ='tests/puma/ToyData/ToyPPIData.txt'
    motif          ='tests/puma/ToyData/ToyMotifData.txt'
    expression_data='tests/puma/ToyData/ToyExpressionData.txt'
    lioness_file   =''
    rm_missing     = False
    output_file    ='travis_test_panda.txt'
    gt_file        ='tests/panda/union_test_panda.txt'

    #1. Union
    panda_obj      = Panda(expression_data, motif, ppi, save_tmp=False, remove_missing=rm_missing,
                      keep_expression_matrix=bool(lioness_file), modeProcess='union')
    panda_obj.save_panda_results(output_file)
    res=pd.read_csv(output_file, sep=' ', header=None)
    gt =pd.read_csv(gt_file, sep=' ', header=None)
    pd.testing.assert_frame_equal(res,gt,check_less_precise=False,check_exact=False)

    #2. Legacy
    panda_obj = Panda(expression_data, motif, ppi, save_tmp=True, remove_missing=rm_missing,
                      keep_expression_matrix=True, save_memory=True, modeProcess='legacy')
    panda_obj.save_panda_results(output_file)
    gt_file         ='tests/panda/legacy_test_panda.txt'
    res = pd.read_csv(output_file, sep=' ', header=None)
    gt = pd.read_csv(gt_file, sep=' ', header=None)
    pd.testing.assert_frame_equal(res,gt,check_less_precise=False,check_exact=False)
    print('Test panda passed was successful!')
