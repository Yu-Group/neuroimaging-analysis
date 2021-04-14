

% figure 4 permutation test
rt = '/Users/Gokul/Dropbox/Manuscript_ExLLSM/GokulWorkspace/Figure4_PermutationTest/';
cd(rt)
load('Fig4_PermutationTest.mat')
opts = {'CmpFunction', @mean};


[~, pValue_S16_pmasDiameter] = permTest(PMAS_avgDia_S16_L5, PMAS_avgDia_S16_L6, opts{:});
[~, pValue_4J] = permTest(somaVol_4J_L5, somaVol_4J_L6, opts{:});
[~, pValue_4K] = permTest(somaVol_4K_L6_cont, somaVol_4K_L6_inter, opts{:});
[~, pValue_4L_L5_L6_all] = permTest(PMASLength_4L_L5, [PMASLength_4L_L6_cont;PMASLength_4L_L6_inter], opts{:});
[~, pValue_4L_L5_L6_cont] = permTest(PMASLength_4L_L5, PMASLength_4L_L6_cont, opts{:});
[~, pValue_4L_L5_L6_inter] = permTest(PMASLength_4L_L5, PMASLength_4L_L6_inter, opts{:});
[~, pValue_4L_L6_cont_L6_inter] = permTest(PMASLength_4L_L6_cont, PMASLength_4L_L6_inter, opts{:});
